"""
===============
Title: Plot IVT Field
Author: Will Rudisill
Description: Plots global IVT field from Netcdf onto Basemap, saves to PNG
I: Netcdf Files of Integrated Vapor Transport 
O: PNG file
===============
"""


from netCDF4 import Dataset, num2date
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import rcParams

ncfile = 'ivt_files/pgbhnl.gdas.20100601-20100605'

def saveplots(fname):
	fil = fname
	rcParams['xtick.labelsize']=1
	rcParams['ytick.labelsize']=1
	ds = Dataset(fil)
	lons = ds.variables['longitude'][:]
	lats = ds.variables['latitude'][:]
	ivt = ds.variables['ivt'] #precipitable water 
	time = ds.variables['time'][:]
	lon_0 = -140
	lat_0 = 20
	lon, lat = np.meshgrid(lons, lats)
	# Loops through each time dimension
	for i in range(20): 
		#ivt at current timestep
		ivt_intime = ivt[i][0]
		#convert the time var to a readable format to use as a title
		header = num2date(time[i], 'hours since 1801-01-01 00:00:00.0', 'gregorian').strftime('%Y%m%d_%H')
		#start a basemap
		m = Basemap(projection='cyl', resolution='c')
#, lat_0 = lat_0, lon_0=lon_0, urcrnrlat=60, urcrnrlon=-80, llcrnrlat=-20, llcrnrlon=-280)
		xi, yi = m(lon, lat)
		m.drawparallels(np.arange(-80., 81., 20.), labels=[1,0,0,0], fontsize=5)
		# m.drawparallels(np.arange(-80., 81., 20.), labels=[0,0,0,0], fontsize=5)
		m.drawmeridians(np.arange(-180., 181., 20.), labels=[0,0,0,1], fontsize=5)
		m.drawmeridians(np.arange(-180., 181., 20.), labels=[0,0,0,0], fontsize=5)
		# # Add Coastlines, States, and Country Boundaries
		m.drawcoastlines()
		m.drawstates()
		m.drawcountries()
		cs = m.pcolor(xi,yi,ivt_intime,latlon=True,vmin=0, vmax=6000)
		cb = m.colorbar(cs, "bottom", size='5%', pad='2%')
                cb.set_label("IVT kgm-1s-1")
                #plt.savefig(header+'IVT.png', format='png', dpi=300, transparent=True)
		plt.show()
        plt.close()
	print 'done'
	ds.close()




saveplots(ncfile)




