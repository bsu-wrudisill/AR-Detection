from skimage.graph import route_through_array
from scipy.spatial import distance
import numpy as np
from atmDetect import Logger

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
    
    def indices_to_lat_lon(index): 
        #input a TUPLE () of grid values 
        x = lats_mesh[index]   #
        y = lons_mesh[index]   # 
        return x,y

    def distance_matrix(self, label_indices_alt):
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
            Logger().failure(e)

    def least_cost(self, label_indices, ivt):
   
        # Cost array creation. Set all values to 9999.999
        cost_arr   = np.ones_like(ivt)*9999.999 

        # Replace 9999 values of AR structure w/ cost
        cost_arr[label_indices] = 5000.0 - ivt[label_indices]

        #start, end are tuples of indices (row,col)
        nrow,ncol      = ivt.shape
        indices, cost  = route_through_array(cost_arr, self.start, self.end)
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


    def landfall_locator(self, land_mask):
        # a/c are the list of array grid coords A = [np.ndarray[x,y], .... ] for 
        ###### Could also do something like this...
        # the center path and the land mask ...
        # land = map(tuple, a)
        # c = map(tuple, c)
        # list(set(a).intersection(c))
        # map(np.array, list(set(a).intersection(c)))


        p       = land_mask[self.ind]
        p_diff  = np.diff(p)
        p_where = np.argwhere(abs(p_diff) == 1)
        self.p_where = p_where
        self.landfall_location = [] 
        try:
            for i in p_where:
                self.landfall_location.append((self.ind[0][i[0]],self.ind[1][i[0]]))
            # returns a tuple of grid indices of landfall location
        except:
            self.landfall_location = None
        



        # NOTE: Double Underscores
    def __init__(self, label_indices_alt, label_indices, ivt, land_mask):
        self.distance_matrix(label_indices_alt)
        if self.passing == True:
            self.least_cost(label_indices, ivt)
            self.landfall_locator(land_mask)
            self.label_indices = label_indices
        else:
            self.landfall_location = None
            self.p_where           = None
            self.path              = None
            self.path_len          = None
            self.ind               = None








