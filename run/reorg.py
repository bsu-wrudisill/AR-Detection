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


#fname                   = '../data/IVT_19960516-19960520.nc'
fname                   = '../data/global_veg.nc'
ds                      = Dataset(fname, format='NETCDF4_CLASSIC')        
lons                    = ds.variables['lon_0'][:]
lats                    = ds.variables['lat_0'][:]
_lons_mesh, lats_mesh   = np.meshgrid(lons, lats)
lons_mesh               =  np.where (_lons_mesh > 180.0, _lons_mesh - 360, _lons_mesh) 



veg                     = ds.variables['VEG_P0_L1_GLL0']
# ivt                     = ds.variables['ivt'][0,0,:,:]
# wnd                     = ds.variables['w_dir'][0,0,:,:] # In Units of Radians
# wnd_                    = wnd/np.pi * 180.0
# wnd_360                 = np.where(wnd_ < 0, 360 + wnd_, wnd_)

# # u and v components (coordinates of unit vector)
# u_i = np.cos(wnd)
# v_i = np.sin(wnd)

# #Components of IVT 
# u_ivt = ivt*u_i/1000.0
# v_ivt = ivt*v_i/1000.0  

fig, [ax1, ax2] = plt.subplots(2,1)

vv = veg[0, :,:] 


def lat_lon_to_indices(lat,lon):
    lat = round(lat * 2)/2   # rounds to nearest .5
    lon = round(lon * 2)/2   # rounds to nearest .5   
    x = np.where(lats_mesh == lat)[0]
    y = np.where(lons_mesh == lon)[1]
    return x,y

x,y = lat_lon_to_indices(-49.029864, -69.162492)

vv[:,y] = 1000.0
vv[x,:]  = 1000.0


ax1.imshow(vv)
levels = np.linspace(-180, 180, 10)
im = ax2.contourf(lons_mesh, interpolation='None', levels=levels)

plt.colorbar(im)
plt.show()



# m = Basemap(projection='cyl', 
#             resolution='c', 
#             llcrnrlat= -90,
#             urcrnrlat= 90,
#             llcrnrlon= 0,
#             urcrnrlon= 360)


# #shifted_lons = m.shiftdata(lons, lon_0=180)
# lon, lat = np.meshgrid(lons, lats)
# xi, yi = m(lon, lat)

# #-------------------------------------------------------------------------#
# # 1.a Draw State/country Borders
# #-------------------------------------------------------------------------#

# m.drawparallels(np.arange(-80., 81., 20.), labels=[1,0,0,0], fontsize=5)
# m.drawmeridians(np.arange(-180., 181., 20.), labels=[0,0,0,1], fontsize=5)
# m.drawcoastlines()
# m.drawstates()
# m.drawcountries()


#-------------------------------------------------------------------------#
# 1.b Plot 2d Field 
#-------------------------------------------------------------------------#

# Mask out 0 Values using numpy mask

#varmask = np.ma.masked_less(ivt, 0)
#cs = m.pcolor(xi,yi,varmask,latlon=True, vmin =250, vmax=1000.0)

#-------------------------------------------------------------------------#
# 1.c Plot Wind Barbs
#-------------------------------------------------------------------------#
        
# Create a sparse array for plotting barbs to prevent clutter
#yy = np.arange(0, yi.shape[0], 3)
#xx = np.arange(0, xi.shape[1], 3)
#points  = np.meshgrid(yy, xx)                                 
#m.quiver(xi[points], yi[points], u_ivt[points], v_ivt[points], scale = 100, pivot='mid', width =0.001, color='grey') 

# class FxToArray(np.ndarray):

#     def __new__(cls, input_array, info=None):
#         # Input array is an already formed ndarray instance
#         # We first cast to be our class type
#         obj = np.asarray(input_array).view(cls)
#         # add the new attribute to the created instance
#         obj.info = info
#         # Finally, we must return the newly created object:
#         return obj



#     def __array_finalize__(self, obj):
#         # see InfoArray.__array_finalize__ for comments
#         if obj is None: return
#         self.info = getattr(obj, 'info', None)





