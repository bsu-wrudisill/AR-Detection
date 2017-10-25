import numpy as np
from skimage.graph import route_through_array
from scipy.spatial import distance
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy.ma as ma




fn = '../data/outfiles/1996-05-16_12_OBJECT_1.0.npy'
ivt = np.load(fn)[0]


label_indices = np.where(ivt > 0)
label_indices_alt = np.argwhere(ivt > 0)

def Distance_Matrix(label_indices_alt):
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
            
        else:
            start = None
            end   = None
            passing = False

    # Set start and end as the same ....
    except Exception as e:
        start = None
        end   = None
        passing = False
    return start, end

start, end = Distance_Matrix(label_indices_alt)

def Least_Cost(label_indices, ivt, start, end):
    # Cost array creation. Set all values to 9999.999

    # The squre of the max ivt value, times 1e3. This guarenteed 
    # the background is at least 1e3 times greater than the lowest value that
    # can possibly be met w/in the AR region

    cost_arr   = np.ones_like(ivt)*ivt.max()**2*1e3

    # Create cost map w/in AR structure. The difference between the cell and the max cell
    # within the object, squared. 

    cost_arr[label_indices] = (ivt[label_indices].max() - ivt[label_indices])**2

    #start, end are tuples of indices (row,col)
    nrow,ncol      = ivt.shape
    indices, cost  = route_through_array(cost_arr, start, end)
    indices        = np.array(indices).T

    # When we index np arrays, we want a tuple of 1-d arrays
    ind       = (indices[0], indices[1])
    path      = np.zeros_like(cost_arr)
    path[ind] = 1
    path_len  = len(indices[0])
    return path


# cost_arr   = np.zeros_like(ivt)
# cost_arr[label_indices] = round((ivt[label_indices].max() - ivt[label_indices])**2.0)


path = Least_Cost(label_indices, ivt, start, end)
path = ma.masked_equal(path, 0)

fig, ax = plt.subplots(1,1)
land  = np.load('../data/land.npy') # Check me 
land = ma.masked_equal(land, 0)
levels = np.linspace(100., 2000., 20)#39

ax.imshow(land, alpha = 1.0, cmap='Greys',zorder=1) 
ivtclr = ax.contourf(ivt, cmap='magma_r',alpha = .8, levels=levels, zorder=2)
ax.imshow(path, alpha = 1.0, interpolation='None', cmap='gray', zorder=4)













    