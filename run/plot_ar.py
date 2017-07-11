import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import numpy as np                

#-------------------------------------------------------------------------#
# Feature: Write out to a netcdf file.
#ds = Dataset(fname, 'r+') # not working; HDF5 error
#AR_nc = ds.createVariable('Atmospheric_River', np.float32, ('time', 'level', 'lat', 'lon'))
#AR_nc.grid_type = 'latitude/longitude'
#AR_nc.units = 'none'
#AR_nc.long_name = 'Atmospheric River Structure'
#AR_nc[0,:,:,:] = new_arr
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
# 1.0 Plot with Basemap
#-------------------------------------------------------------------------#

def make_plot(lons, lats, data_array, u_ivt, v_ivt, hr_time_str, save_me):

    #------------------------------------------------------------------------#-
    # FUNCTION PARAMETERS
    #------------------------------------------------------------------------#-
    # data_array:  2d field to plot (probable ivt)
    # u_ivt     :  IVT in u direction; used to calc wind vectors
    # v_ivt     :  IVT in v direction; used to calc wind vectors
    # lons
    # lats 
    #------------------------------------------------------------------------#-

    # Basemap object, cylindrical projection
    m = Basemap(projection='cyl', 
            resolution='c', 
            llcrnrlat= -90,
            urcrnrlat= 90,
            llcrnrlon= 0,
            urcrnrlon= 360)

    # shift lons
    lons = m.shiftdata(lons, lon_0=180)
    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)

    #-------------------------------------------------------------------------#
    # 1.a Draw State/country Borders
    #-------------------------------------------------------------------------#
    
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1,0,0,0], fontsize=5)
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0,0,0,1], fontsize=5)
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()


    #-------------------------------------------------------------------------#
    # 1.b Plot 2d Field 
    #-------------------------------------------------------------------------#
    
    # Mask out 0 Values using numpy mask

    varmask = np.ma.masked_less(data_array, 1)
    cs = m.pcolor(xi,yi,varmask,latlon=True, vmin = 200.0, vmax=1000.0)

    #-------------------------------------------------------------------------#
    # 1.c Plot Wind Barbs
    #-------------------------------------------------------------------------#
                
    # Create a sparse array for plotting barbs to prevent clutter
#    yy = np.arange(0, yi.shape[0], 3)
#    xx = np.arange(0, xi.shape[1], 3)
#    points  = np.meshgrid(yy, xx)                                 
#    m.quiver(xi[points], yi[points], u_ivt[points], v_ivt[points], scale = 100, pivot='mid', width =0.001, color='grey') 
                


    #-------------------------------------------------------------------------#
    # 1.d Save Figure
    #-------------------------------------------------------------------------#
    if save_me == True:
    # Set Colorbar
        cbar = m.colorbar(cs,location='bottom',pad="5%")
        cbar.set_label('mm')
               
        # Save Plot
            
        plt.title(hr_time_str)
        plt.savefig('AR_'+hr_time_str ,format='png', dpi = 700)
        plt.close()
                
