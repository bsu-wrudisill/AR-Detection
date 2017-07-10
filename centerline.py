from skimage.graph import route_through_array
from scipy.spatial import distance
import numpy as np

class Centerline():
                
        '''
        Centerline object
        '''        
        
        def distance_matrix(self, label_indices_alt):
                # Find the two points with the greatest distance from a shape
                # The entry array is a binary array where 1 == object
                # returns the start and end points of shape 

            dmat   = distance.cdist(label_indices_alt, label_indices_alt)    # Scipy euclidian distance matrix
            points = np.argwhere(dmat == dmat.max())                                   # Coords where there is max distance  
                
            if len(points) > 0:
                points = points[0]
                start  = label_indices_alt[points[0]]
                end    = label_indices_alt[points[1]]
                self.start = start 
                self.end   = end 
                        
            else:
                self.start = 0
                self.end   = 0

        def least_cost(self, label_indices, ivt):
   
                # Cost array creation. Set all values to -9999
                cost_arr   = np.ones_like(ivt)*9999.999 

                # Replace -9999 values of AR structure w/ cost
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

        # NOTE: Double Underscores
        def __init__(self, label_indices_alt, label_indices, ivt):
                self.distance_matrix(label_indices_alt)
                self.least_cost(label_indices, ivt)
