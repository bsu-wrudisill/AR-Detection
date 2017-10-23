from skimage.graph import route_through_array
from scipy.spatial import distance
import numpy as np
#from atmDetect import Logger

############ CHANGE ME ##########
############ CHANGE ME ##########


'''
I should use a call method rather than what I'm doing here; create class object using static data,
call method to apply to each AR object... 

'''

############ CHANGE ME ##########
############ CHANGE ME ##########


class Centerline():
    '''
    Centerline object
      1. Finds the two most distance points within object
      2. Creates cost Array based on IVT mag
      3. Finds least cost path between points
      4. Finds landfalling location along path
    '''        

    # this seems dumb   
    def Indices_to_Lat_Lon(self, index): 
        #input a TUPLE () of grid values 
        x = self.lats_mesh[index]   #
        y = self.lons_mesh[index]   # 
        return x,y

    def Distance_Matrix(self, label_indices_alt):
        # Find # TODO: he two points with the greatest distance from a shape
        # The entry array is a binary array where 1 == object
        # returns the start and end points of shape 
        try: 
            dmat   = distance.cdist(label_indices_alt, label_indices_alt)       # Scipy euclidian distance matrix
            points = np.argwhere(dmat == dmat.max())                                   # Coords where there is max distance  
                
            if len(points) > 0:
                points = points[0]
                start  = label_indices_alt[points[0]]
                end    = label_indices_alt[points[1]]
                self.start = start 
                self.end   = end 
                self.passing = True
        
            else:
                self.start = None
                self.end   = None
                self.passing = False

        # Set start and end as the same ....
        except Exception as e:
            self.start = None
            self.end   = None
            self.passing = False

    def Least_Cost(self, label_indices, ivt):
        # Cost array creation. Set all values to 9999.999

        # The squre of the max ivt value, times 1e3. This guarenteed 
        # the background is at least 1e3 times greater than the lowest value that
        # can possibly be met w/in the AR region
        cost_arr   = np.ones_like(ivt)*ivt.max()**2*1e3

        # Create cost map w/in AR structure. The difference between the cell and the max cell
        # within the object, squared. 
        cost_arr[label_indices] = (ivt[label_indices].max() - ivt[label_indices])**2


        # Apply route thru array method
        # Takes a starting, end point and a cost grid 
        nrow,ncol      = ivt.shape
        indices, cost  = route_through_array(cost_arr, self.start, self.end)

        # The grid indices of the least-cost path
        indices        = np.array(indices).T

        # When we index np arrays, we want a tuple of 1-d arrays
        ind       = (indices[0], indices[1])
        path      = np.zeros_like(cost_arr)
        path[ind] = 1
        path_len  = len(indices[0])
                
        self.path      = path
        self.cost      = cost
        self.path_len  = path_len
        self.ind       = ind 


    def Landfall_Locator(self):
        # a/c are the list of array grid coords A = [np.ndarray[x,y], .... ] for 
        ###### Could also do something like this...
        # the center path and the land mask ...
        # land = map(tuple, a)
        # c = map(tuple, c)
        # list(set(a).intersection(c))
        # map(np.array, list(set(a).intersection(c)))

        p       = self.coastline[self.ind]
        p_diff  = np.diff(p)
        p_where = np.argwhere(abs(p_diff) == 1)
        self.landfall_location = [] 
        try:
            for i in p_where:
                self.landfall_location.append((self.ind[0][i[0]],self.ind[1][i[0]]))
            # returns a tuple of grid indices of landfall location
        except:
            self.landfall_location = None

        

    def Measr_pLength(self):

        # Calculate the dimensions of each grid cell. Assumes spherical earth 
        earth_rad = 6371.0
        h       = np.cos(np.abs(self.lats_mesh)*np.pi/180.)*earth_rad  # the radius of the great circle by latitude
        gc      = np.pi*2.*h                                      # the circumference of the greate circle
        grid_dx = gc/720.0                                        # horizontal grid cell distance 
        # calculate vertical grid cell distance; this is the same for all lats 
        grid_dy = earth_rad*2.0*np.pi/720.0


        # Calculate distance along the AR track 
        x,y = self.ind
        distance = []
        for i in xrange(0, len(x)-1, 1):              # loop through each grid cell in centerline
            x_ij    = grid_dx[x[i],y[i]] * 0.5        # x_ij
            x_ijP   = grid_dx[x[i+1],y[i+1]] * 0.5    # x_ij + 1 

            scale_x = np.abs(x[i] - x[i+1])           # difference between x indice
            scale_y = np.abs(y[i] - y[i+1])           # difference between y indice
            DX = (x_ij + x_ijP)**2 * scale_x          
            DY = grid_dy**2  * scale_y 
            distance.append(np.sqrt(DX + DY))         # distance formula
        
        self.path_length = sum(distance)

        # find landfall point



        # NOTE: Double Underscores
    def __init__(self, label_indices_alt, label_indices, ivt, lats_mesh):

        #-------------------------------------------------------------------------#
        # Global Values 
        #-------------------------------------------------------------------------#
        land       = np.load('../data/land.npy')
        self.coastline  = np.load('../data/west_coast.npy')
        #-------------------------------------------------------------------------#

        self.Distance_Matrix(label_indices_alt)
        self.lats_mesh = lats_mesh

        if self.passing == True:
            self.Least_Cost(label_indices, ivt)
            self.Landfall_Locator()
            self.label_indices = label_indices
            self.Landfall_Locator()
            self.Measr_pLength()
        else:
            self.landfall_location = 'Failure'
            self.path              = 'Failure'
            self.path_len          = 'Failure'
            self.ind               = 'Failure'








