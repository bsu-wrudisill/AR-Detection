from skimage.graph import route_through_array
from scipy.spatial import distance
import numpy as np


class Centerline():
    '''
    Centerline object
      1. Finds the two most distance points within object
      2. Creates cost Array based on IVT mag
      3. Finds least cost path between points
      4. Finds landfalling location along path
    '''        
    def __init__(self, ivt, lats_mesh):

        #These things should not change for each timestep!! Then we call for each blob

        #-------------------------------------------------------------------------#
        # Global Values 
        #-------------------------------------------------------------------------#
        land       = np.load('../data/land.npy')
        coastline  = np.load('../data/west_coast.npy')
        #-------------------------------------------------------------------------#
        self.ivt       = ivt
        self.lats_mesh = lats_mesh

    def indices_to_lat_lon(index): 
        #input a TUPLE () of grid values 
        x = self.lats_mesh[index]   #
        y = self.lons_mesh[index]   # 
        return x,y

    def distance_matrix(self,label_indices_alt):
        # Find # TODO: he two points with the greatest distance from a shape
        # The entry array is a binary array where 1 == object
        # returns the start and end points of shape 
        try: 
            dmat   = distance.cdist(label_indices_alt, label_indices_alt)   # Scipy euclidian distance matrix
            points = np.argwhere(dmat == dmat.max())                        # Coords where there is max distance  
                
            ######### WHY WOULD THIS FAIL????? ########
            if len(points) > 0:
                points = points[0]
                start  = label_indices_alt[points[0]]
                end    = label_indices_alt[points[1]]
                passing = True
        
            else:
                start = None
                end   = None
                passing = False

        # Set start and end as the same ....
        except Exception as e:
            start = None
            end   = None
            passing = False
            raise e 

        # return things 
        return {'start':start, 'end':end, 'passing':passing }

    def least_cost(self, label_indices, start, end):
   
        # Cost array creation. Set all values to 9999.999
        cost_arr   = np.ones_like(self.ivt)*9999.999 

        # Replace 9999 values of AR structure w/ cost
        cost_arr[label_indices] = 5000.0 - self.ivt[label_indices]

        #start, end are tuples of indices (row,col)
        nrow,ncol      = self.ivt.shape
        indices, cost  = route_through_array(cost_arr, start, end)
        indices        = np.array(indices).T

        # When we index np arrays, we want a tuple of 1-d arrays
        center_indices_tuple              = (indices[0], indices[1])
        
        # global grid. 1 == center line, everyhwere else == 0 
        center_grid                       = np.zeros_like(cost_arr)
        center_grid[center_indices_tuple] = 1
        
        # return a dictionary                    
        return {'center_indices_tuple':center_indices_tuple, 'center_grid':center_grid}
    

    def landfall_locator(self, center_indices_tuple):
        # a/c are the list of array grid coords A = [np.ndarray[x,y], .... ] for 
        ###### Could also do something like this...
        # the center path and the land mask ...
        land         = map(tuple, self.coastline)
        intersection = list(set(land).intersection(center_indices_tuple))
        print intersection

        #map(np.array, list(set(a).intersection(c)))
        # p       = self.land_mask[self.ind]
        # p_diff  = np.diff(p)
        # p_where = np.argwhere(abs(p_diff) == 1)
        # self.p_where = p_where
        # self.landfall_location = [] 
        # try:
        #     for i in p_where:
        #         self.landfall_location.append((self.ind[0][i[0]],self.ind[1][i[0]]))
        #     # returns a tuple of grid indices of landfall location
        # except:
        #     self.landfall_location = None
        

    def path_length(self, center_indices_tuple):
        # Calculate the dimensions of each grid cell. Assumes spherical earth 
        earth_rad = 6371.0
        h       = np.cos(np.abs(self.lats_mesh)*np.pi/180.)*earth_rad  # the radius of the great circle by latitude
        gc      = np.pi*2.*h                                      # the circumference of the greate circle
        grid_dx = gc/720.0                                        # horizontal grid cell distance 
        # calculate vertical grid cell distance; this is the same for all lats 
        grid_dy = earth_rad*2.0*np.pi/720.0


        # Calculate distance along the AR track 
        x,y = center_indices_tuple
        distance = []
        for i in xrange(0, len(x)-1, 1):              # loop through each grid cell in centerline
            x_ij    = grid_dx[x[i],y[i]] * 0.5        # x_ij
            x_ijP   = grid_dx[x[i+1],y[i+1]] * 0.5    # x_ij + 1 
            scale_x = np.abs(x[i] - x[i+1])           # difference between x indice
            scale_y = np.abs(y[i] - y[i+1])           # difference between y indice
            DX = (x_ij + x_ijP)**2 * scale_x          
            DY = grid_dy**2  * scale_y 
            distance.append(np.sqrt(DX + DY))         # distance formula
        path_distance = sum(distance)
        return path_distance


    def __call__(self, label_indices, label_indices_alt):
        # Distance matrix to get start/end pts
        # --> Least-cost(start, end) to get centerline
        # --> Path_length --> landfall_locator
        dm = self.distance_matrix(label_indices_alt)
        lc = self.least_cost(label_indices, dm['start'], dm['end'])
        pl = self.path_length(lc['center_indices_tuple'])
        ll = self.landfall_locator(lc['center_indices_tuple'])
        output = {'path_length':pl, 'll':ll, 'start':dm['start']}


        # find landfall point