import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import numpy as np                
from netCDF4 import Dataset
import numpy as np
import glob
import json
from scipy.ndimage import filters, morphology, measurements
from scipy import ndimage
from datetime import datetime, timedelta
import logging
import gc 
import time 
import traceback
from skimage.measure import regionprops
from centerline import Centerline


if 'fig' in locals():  # Delete fig object if it exists 
    del fig


def lat_lon_to_indices(lat,lon):
	# lat and lon are geographic coordinates 

    lat = round(lat * 2)/2   # rounds to nearest .5
    lon = round(lon * 2)/2   # rounds to nearest .5   
    x = np.where(lats_mesh == lat)[0]
    y = np.where(lons_mesh == lon)[1]
    return x,y

def indices_to_lat_lon(index): 
	#input a TUPLE () of grid values 
    x = lats_mesh[index]   #
    y = lons_mesh[index]   # 
    return x,y

def roi():
	idx = np.where((lons_mesh > -125.) & (lons_mesh < -112.) & (lats_mesh > 30.) & (lats_mesh < 50.))
	foo = np.zeros_like(lons_mesh) 
	foo[idx] = 1
	n_idx = np.where(foo != 1)
	return idx, n_idx


# x,y = lat_lon_to_indices(-49.029864, -69.162492)
def calc_ivt(filename, time):
	# calculate IVT and other things...
	ds                      = Dataset(filename, format='NETCDF4_CLASSIC')        
	lons                    = ds.variables['lon_0'][:]
	lats                    = ds.variables['lat_0'][:]
	_lons_mesh, lats_mesh   = np.meshgrid(lons, lats)
	lons_mesh               = np.where (_lons_mesh > 180.0, _lons_mesh - 360, _lons_mesh) 

	ivt =  np.zeros((20, 361, 720))
	g = 9.81 
	dp = 25.00 #pa

	for i in range(20):
	    spfh = ds.variables['SPFH_P0_L100_GLL0'][time,i,:,:]
	    vgrd = ds.variables['VGRD_P0_L100_GLL0'][time,i,:,:]
	    ugrd = ds.variables['UGRD_P0_L100_GLL0'][time,i,:,:]
	    phi =np.arctan2(vgrd,ugrd)
	    tV = (ugrd**2)*(vgrd**2)**.5
	    ivt[i,:,:] = spfh*tV/g*dp
	ivt_integral = np.ndarray.sum(ivt, axis=0)
	vgrd =  ds.variables['VGRD_P0_L100_GLL0'][time,:,:,:]
	ugrd =  ds.variables['UGRD_P0_L100_GLL0'][time,:,:,:]
	ds.close()
	return ivt_integral, vgrd, ugrd, lons_mesh, lats_mesh, lons, lats

#-------------------------------------------------------------------------#
# Global Values 
#-------------------------------------------------------------------------#

ivt_min    = 120.0                     # Minimum IVT value in kg/ms to be retained
size_mask  = 100000.0                  # Min Grid cell size of object
cell_to_km = 50.0                   # km
land 	   = np.load('../data/land.npy')
coastline  = np.load('../data/west_coast.npy')

#-------------------------------------------------------------------------#


# Calculate IVT
ivt, vgrd, ugrd, lons_mesh, lats_mesh, lons, lats = calc_ivt('../data/pgbhnl.gdas.19960516-19960520.nc', 10)


# Convert to binary. 1 is above min value
threshold_array = np.where(ivt > ivt_min, 1, 0)

# Labeled components. given values of 1:n
# Each 'blob' identified is assigned a label 
label_array, num_labels = ndimage.measurements.label(threshold_array)

# List of the label values 
labels = np.array(range(num_labels + 1))

# List size, in pixels, of components 
sizes = ndimage.sum(label_array, label_array, range(num_labels + 1))
   
# Labels to keep (corresponding w/ objects gt. than size threshold)
keep_label = labels[sizes > size_mask]


# '''
# Go through AR identifier criteria
# '''

# AR_EXISTS = False

#-------------------------------------------------------------------------#
#  Test Flags; set to true if passing 
#-------------------------------------------------------------------------#
TC_a = False      # Mean IVT 
TC_b = False      # Coherency in IVT Direction  (variance)
TC_c = False      # Object Mean Meridonal IVT (flowing towards a pole)
TC_d = False      # Object/IVT direction Consistency
TC_e = False      # Length/Width Ratio
TC_f = False      # Interior USA Landfalling 
TC_g = False      # NOT crossing equator (we don't want AR crossing eq)
#-------------------------------------------------------------------------#

output = np.zeros_like(land)

for label in keep_label:

	# Indicies of blob of interest
	label_indices     = np.where(label_array == label)
	label_indices_alt = np.argwhere(label_array == label)  # indices in [(a,b), ....] fmt

	# Remove all elements except blob of interest
	sub_array     = np.where(label_array == label, 1, 0)
	blob          = regionprops(sub_array)[0]                    # set to 0; only 1 region        
	blob_length   = blob.major_axis_length*cell_to_km
	blob_width    = blob.minor_axis_length*cell_to_km
	blob_size     = blob.filled_area*cell_to_km


	mean_ivt  	  = np.mean(ivt[label_indices])              # Mean IVT of blob
	# mean_u_i        = np.mean(u_i[label_indices])
	# mean_v_i        = np.mean(v_i[label_indices])
	# mean_wind_dir_  = np.arctan2(mean_v_i,mean_u_i)/np.pi * 180
	# mean_wind_dir   = np.where(mean_wind_dir_< 0, 360 + mean_wind_dir_ , mean_wind_dir_)

	#get center of blob
	cntr = tuple(map(int, map(round, blob.centroid)))
	sub_array[cntr] = 10
	center   = Centerline(label_indices_alt, label_indices, ivt, land)
	output = output + center.path
	print center.landfall_location
# fig, ax1 = plt.subplots(1)
# ax1.imshow(sub_array + center.path)

#ax2.imshow(ivt, vmin=100.0, vmax = 1000.00)



