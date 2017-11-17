import numpy as np
from scipy import ndimage
import json
import matplotlib.pyplot as plt

####################
# Script creates W coast 
####################

class collect_points():
    omega = []
    def __init__(self,array):
        self.array = array
        #self.omega = array??
    def onclick(self,event):
        # print 'xdata=%f, ydata=%f'%(
        #     int(round(event.xdata)),   int(round(event.ydata)))
        self.omega.append((int(round(event.ydata)), int(round(event.xdata))))
    def indices(self):
        plot = plt.imshow(self.array, cmap = plt.cm.hot, interpolation = 'nearest', origin= 'upper')
        fig = plt.gcf()
        ax = plt.gca()
        zeta = fig.canvas.mpl_connect('button_press_event', self.onclick)
        plt.colorbar()
        plt.show()
        return self.omega


land 	  = np.load('../data/land.npy')
coastline = np.zeros_like(land) 

clicked_points = collect_points(land).indices()


# # ------------ Create a grid -------------#
# lats = np.linspace(-90.,90.,361.)[::-1]
# _lons = np.linspace(-180.,179.5,720.)
# lons = np.where (_lons < 0.0, _lons + 180. , _lons - 180.)
# lons_mesh,lats_mesh = np.meshgrid(lons,lats)



# with open('/Users/will/Desktop/mygeodata/coast.json') as json_d:
# 	coast = json.load(json_d)


# points = coast['features'][0]['geometry']['coordinates']

# land 	  = np.load('../data/land.npy')
# coastline = np.zeros_like(land) 

# for i in points:
# 	lon = np.round(i[0]*2)/2.  # round to the nearest .5 
# 	lat = np.round(i[1]*2)/2.  # round to the nearest .5 	
# 	idx = np.where((lats_mesh == lat) & (lons_mesh == lon))
# 	coastline[idx] = 10.





# #get edges of land forms
# sx = ndimage.sobel(land, axis=0, mode='constant')
# sy = ndimage.sobel(land, axis=1, mode='constant')
# sob = np.where(np.hypot(sx, sy) > 0, 1, 0)

# land[np.where(sob == 1)] = 2 # assign border cells 2
# idx, n_idx= roi()
# coastline[np.where(sob ==1)] = 1
# coastline[n_idx]             = 0    # indices of just the california coast

np.save('../west_coast.npy', coastline)


# create distance map 





# def roi():
# 	# please fix me to take lat/lon arguments
# 	# creates a box of values on a grid defined by ur, ul, ll, lr 
# 	idx = np.where((lons_mesh > -125.) & (lons_mesh < -112.) & (lats_mesh > 30.) & (lats_mesh < 50.))
# 	foo = np.zeros_like(lons_mesh) 
# 	foo[idx] = 1
# 	n_idx = np.where(foo != 1)
# 	return idx, n_idx

