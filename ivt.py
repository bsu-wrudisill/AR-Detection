from netCDF4 import Dataset, num2date, date2num
import numpy as np
from datetime import datetime, timedelta
import pandas as pd 


def IVT(infile, outfile):
 	# Infile is the CFS Reanalysis data to convert, out is the name for 
	# the new file. CRSR downloaded from NCAR RDA
	name = outfile
	dataset = Dataset(name,'w' ,format='NETCDF4_CLASSIC')
	level = dataset.createDimension('level', 1)
	lat = dataset.createDimension('lat', 361)
	lon = dataset.createDimension('lon', 720)
	time = dataset.createDimension('time', None)

	# create NetCDF variables
	times = dataset.createVariable('time', np.float64, ('time',))
	times.units = "hours since 1801-01-01 00:00:00.0" 
	times.calendar = "gregorian"
	levels = dataset.createVariable('level', np.int32, ('level',))
	latitudes = dataset.createVariable('latitude', np.float32, ('lat',))
	longitudes = dataset.createVariable('longitude', np.float32, ('lon',))

	#set lat attrbutes
	latitudes.La1 = '90.f'
	latitudes.Lo1 = '0.f'
	latitudes.La2 = '-90.f'
	latitudes.Lo1 = '359.5f'
	latitudes.Di = '0.5f'
	latitudes.Dj = '0.5f'
	latitudes.units = "degrees north"
	latitudes.grid_type = "Latitude/Longitude"
	latitudes.long_name = "latitude" 

	#set lon attrbutes
	longitudes.La1 = '90.f'
	longitudes.Lo1 = '0.f'
	longitudes.La2 = '-90.f'
	longitudes.Lo1 = '359.5f'
	longitudes.Di = '0.5f'
	longitudes.Dj = '0.5f'
	longitudes.units = "degrees east"
	longitudes.grid_type = "Latitude/Longitude"
	longitudes.long_name = "longitude" 
        
        # Create IVT var
	# To assign ivt a var, ivt[i,:,:,:] = var, i is time dim 
	ivt = dataset.createVariable('ivt', np.float32,('time','level','lat','lon'))
	ivt.grid_type = "Latitude/longitude"
	ivt.units = 'kg m-1 s-1'
	ivt.long_name = 'Integrated Vapor Transport'
        
        # Create Wind Direction Var 
        wnd = dataset.createVariable('w_dir', np.float32,('time','level','lat','lon'))
	wnd.grid_type = "Latitude/longitude"
	wnd.units = 'Degrees'
	wnd.long_name = 'Mean Wind Direction'


	#create geo dimensions.we know a-priori its a .5 degree grid 
	latitudes[:] =  np.arange(-90,90.5,.5)
	longitudes[:] =  np.arange(-180,180,.5)


	# Dataset of interest
	subset = Dataset(infile)

	# Take initial time variable from orig data, transform to a useable format
	dt = reduce(lambda x,y: x+y, subset.variables['initial_time0'][0])

	# Create the right time steps for data
	dates = [datetime(int(dt[6:10]),int(dt[0:2]),int(dt[3:5]))+ n*timedelta(hours = 6) for n in range(20)]
	times[:] = date2num(dates,units=times.units,calendar=times.calendar)


	#Loop through pressure and Time levels, calculate IVT
	for i in range(20):
		Q = subset['SPFH_P0_L100_GLL0'][i,:]
		u = subset['UGRD_P0_L100_GLL0'][i,:]
                v = subset['VGRD_P0_L100_GLL0'][i,:]
                tV = (u**2)*(v**2)**.5
                phi =np.arctan2(v,u)/np.pi*180 #first arg is the y dir, second x dir (meridonal, zonal)                
		g = 9.81 
		dp = 25.00 #pa
		ivtval = np.ndarray.sum(Q*tV, axis=0)*1/g*dp
                wndval = np.ndarray.mean(phi, axis=0)
                ivt[i,:,:,:] = ivtval
                wnd[i,:,:,:] = wndval

        subset.close()
	dataset.close()


#IVT('sample/pgbhnl.gdas.20101226-20101231.nc', 'sample/20101226-20101231_IVT.nc')

IVT('input_files/pgbhnl.gdas.20100601-20100605.nc', 'sample/20100601-20100605_IVT.nc')

print 'Im done'
