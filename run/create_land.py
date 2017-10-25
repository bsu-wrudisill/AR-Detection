import numpy as np
from scipy import ndimage


####################
# Script creates W coast 
####################


def roi():
	# please fix me to take lat/lon arguments
	# creates a box of values on a grid defined by ur, ul, ll, lr 
	idx = np.where((lons_mesh > -125.) & (lons_mesh < -112.) & (lats_mesh > 30.) & (lats_mesh < 50.))
	foo = np.zeros_like(lons_mesh) 
	foo[idx] = 1
	n_idx = np.where(foo != 1)
	return idx, n_idx


from reorg import calc_ivt
ivt_timeslice = calc_ivt('../data/pgbhnl.gdas.19960516-19960520.nc', 0)



ivt = ivt_timeslice['ivt']
vgrd = ivt_timeslice['vgrd']
ugrd = ivt_timeslice['ugrd']	
lons_mesh   = ivt_timeslice['lons_mesh']
lats_mesh   = ivt_timeslice['lats_mesh']


land 	  = np.load('../data/land.npy')
coastline = np.zeros_like(land) 

#get edges of land forms
sx = ndimage.sobel(land, axis=0, mode='constant')
sy = ndimage.sobel(land, axis=1, mode='constant')
sob = np.where(np.hypot(sx, sy) > 0, 1, 0)

land[np.where(sob == 1)] = 2 # assign border cells 2
idx, n_idx= roi()
coastline[np.where(sob ==1)] = 1
coastline[n_idx]             = 0    # indices of just the california coast

np.save('../west_coast.npy', coastline)


# create distance map 





