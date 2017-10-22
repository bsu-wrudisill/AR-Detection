from netCDF4 import Dataset, num2date, date2num
import numpy as np
from datetime import datetime, timedelta
from glob import glob
from matplotlib import pyplot as plt

# calculate IVT and other things...

#for time in range(20):
time=0 

filename = '../data/pgbhnl.gdas.19960516-19960520.nc'

ds                      = Dataset(filename, format='NETCDF4_CLASSIC')        
lons                    = ds.variables['lon_0'][:]
lats                    = ds.variables['lat_0'][:]
_lons_mesh, lats_mesh   = np.meshgrid(lons, lats)
lons_mesh               = np.where (_lons_mesh > 180.0, _lons_mesh - 360, _lons_mesh) 
ivt =  np.zeros((20, 361, 720))
g = 9.81 
dp = 25.00 #pa

# loop through pressure levels, calculate IVT 
# for i in range(20):
#     spfh = ds.variables['SPFH_P0_L100_GLL0'][time,i,:,:]
#     vgrd = ds.variables['VGRD_P0_L100_GLL0'][time,i,:,:]
#     ugrd = ds.variables['UGRD_P0_L100_GLL0'][time,i,:,:]
#     phi =np.arctan2(vgrd,ugrd)
#     tV = (ugrd**2)*(vgrd**2)**.5
#     ivt[i,:,:] = spfh*tV/g*dp

# # sum through pressure levels 
# ivt_integral = np.ndarray.sum(ivt, axis=0)

# U and V winds respectively

spfh  = ds.variables['SPFH_P0_L100_GLL0'][time,:,:,:]
vgrd  =  ds.variables['VGRD_P0_L100_GLL0'][time,:,:,:]
ugrd  =  ds.variables['UGRD_P0_L100_GLL0'][time,:,:,:]
speed = np.sqrt(vgrd**2 + ugrd**2) 
ds.close()


def lat_lon_to_indices(lat,lon):
    # lat and lon are geographic coordinates 
    lat = round(lat * 2)/2   # rounds to nearest .5
    lon = round(lon * 2)/2   # rounds to nearest .5   
    x = np.where(lats_mesh == lat)[0][0]
    y = np.where(lons_mesh == lon)[1][0]
    return x,y


#     A = [34., -123.]
#     B = [47., -123.]
#     ax, ay = lat_lon_to_indices(*A)
#     bx, by = lat_lon_to_indices(*B)

#     fig, ax = plt.subplots(2,1)
#     ax[0].imshow(speed[:,86:112, 474])
#     ax[1].imshow(spfh[:,86:112, 474])
#     plt.show()

# plt.quiver(vgrd[0, 80:120, 440:480], ugrd[0, 80:120, 440:480])


def uv2deg(v,u): # v (vertical), u (horizontal)
    # by convention, positve v flows from the south to the north
    # positive u flows from the west to the east 
    deg = np.arctan2(u,v) * 180./np.pi 
    deg = np.where(deg < 0, 360.0 + deg, deg)
    return deg 


def ScaleBYq(v,u,q):
    v_mean = np.sum(v*q, axis=0)/np.sum(q, axis=0)
    u_mean = np.sum(u*q, axis=0)/np.sum(q, axis=0)
    return uv2deg(v_mean,u_mean)











#np.where(ivt_integral > 250.0, ivt_integral, 0)

