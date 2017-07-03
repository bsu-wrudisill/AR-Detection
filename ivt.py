from netCDF4 import Dataset, num2date, date2num
import numpy as np
from datetime import datetime, timedelta
from glob import glob
import os
import logging

def IVT(infile,outfile):
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
        
        # Create World Landcover Var
        land = dataset.createVariable('landcover', np.float32, ('time', 'level', 'lat', 'lon'))
	land.grid_type = "Latitude/longitude"
	land.units = 'Degrees'
	land.long_name = 'Landcover Type; from Modis'

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


        # read land data so we don't read it in a loop;
        land_input = Dataset('../global_veg.nc')
        landval = land_input['VEG_P0_L1_GLL0'][0, ::-1, :]


	#Loop through pressure and Time levels, calculate IVT
	for i in range(20):
		Q = subset['SPFH_P0_L100_GLL0'][i,:]
		u = subset['UGRD_P0_L100_GLL0'][i,:]
                v = subset['VGRD_P0_L100_GLL0'][i,:]
                v_mean = np.ndarray.mean(v, axis=0)
                u_mean = np.ndarray.mean(u, axis=0)
                
                # IVT Calc
                tV = (u**2)*(v**2)**.5
		g = 9.81 
		dp = 25.00 #pa
		ivtval = np.ndarray.sum(Q*tV, axis=0)*1/g*dp

                # Wind Calc 
                phi =np.arctan2(v_mean,u_mean) # Keep in pi units. first arg is the y dir, second x dir (meridonal, zonal)  
                wndval = phi
                
                # Assign to variables (time, pressure (1), lat, lon)
                ivt[i,:,:,:]   = ivtval
                wnd[i,:,:,:]   = wndval
                land[i,:,:,:]  = landval

        subset.close()
	dataset.close()
        # End Funcion 

def try_IVT(infile,outfile):
        try: 
                return IVT(infile,outfile)
        except Exception as e:
                logging.exception('message')
                

os.chdir('cfsr_nc/')

file_list     = glob('*')
out_file_list = map(lambda X: '../ivt_files/'+X, file_list)
map(try_IVT, file_list, out_file_list)

os.chdir('../.')


