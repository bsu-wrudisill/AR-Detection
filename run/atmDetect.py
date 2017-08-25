
"""
===============
Title: Atmospheric River Detection 
Author: Will Rudisill
Description: Finds Atmospheric Rivers from CFS reanalysis (or other) dataset
I: Netcdf Files of integrated vapor transport 
O: TBD
===============
"""

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


#-------------------------------------------------------------------------#
# TO DO List
# 1. Land mask
# 2. Subtract seasonal means
# 3. Properly calculate grid cell size
# 4. Fix poleward component test; northern needs to be positive, southern negative
# 5. Check on len/width ratio math 

#-------------------------------------------------------------------------#
np.seterr(all='print')


#-------------------------------------------------------------------------#
#  Functions 
#-------------------------------------------------------------------------#

class Logger():        
        '''
        Logger Class. Log write to log based on failure/success

        '''
        def __init__(self):
                pass 

        def success(self, obj):
                log_name = datetime.now().strftime('%Y-%m-%d_%H:%M:%S') + '_AR_Detect.log'
                logger = logging.getLogger('TK')
                hdlr = logging.FileHandler('log_name')
                formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
                hdlr.setFormatter(formatter)
                logger.addHandler(hdlr) 
                logger.setLevel(logging.INFO)
                logger.info('completed_%s', obj)

        def failure(self,obj):
                log_name = datetime.now().strftime('%Y-%m-%d_%H:%M:%S') + '_AR_Detect.log'
                logger = logging.getLogger(log_name)
                hdlr = logging.FileHandler(log_name)
                formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
                hdlr.setFormatter(formatter)
                logger.addHandler(hdlr) 
                logger.setLevel(logging.DEBUG)
                logger.exception('failure_%s', obj )
        

def smallest_angle(a,b):
        # Returns smalles angle in degrees for [-180, 180]
        tup = (360 - abs(a - b), a-b)
        return min(map(abs, tup))

def format_date(hours):
        t0 = datetime(1801, 01, 01, 00)
        dt = t0 + timedelta(hours=hours)
        return dt


 def timeit(method):
        def timed(*args, **kw):
                ts = time.time()
                result = method(*args, **kw)
                te = time.time()
                print  method.__name__,  str(te-ts) + 's'
                
                return result
        return timed


#@timeit
def FindAR_Wrapper(fname):
        #Xnp.seterr(all='print')
        
        # Open nc dataset
        ds  = Dataset(fname, format='NETCDF4_CLASSIC')        
        
        # Read in the static fields, pass them into FX

        # Get lat/lon for plotting purposes
        lons  = ds.variables['lon_0'][:]
        lats  = ds.variables['lat_0'][:]
        
        # meshgrid
        lons_mesh, lats_mesh  = np.meshgrid(lons, lats)
        
        # Subset of IVT
        ivt   = ds.variables['ivt']
        
        # Subset of wind_direction; 
        wnd   = ds.variables['w_dir'] # In Units of Radians

        # Landcover 
        global_cover_class  = ds.variables['landcover'][0, 0, :, :]
        land_mask           = np.where(global_cover_class < 1000, 1, 0)

        # Number of timesteps in file
        tlen = range(len(ds.variables['time'][:]))
        
        # Map function to timelist 
        # Note that FindAR exists w/in function scope... 
        try:
                map(lambda X: FindAR(ds, fname, land_mask, lons_mesh, lats_mesh, wnd, ivt, lons, lats,  X), tlen)
                Logger().success(fname)
                
        ## THIS PART IS IMPORTANT; Mutliprocessing will not give error traceback
        except Exception as e:
               Logger().failure(fname)
               raise e

        #Close dataset
        ds.close()
        # Log Success 


def FindAR(dataset, fname, land_mask, lons_mesh, lats_mesh, wnd, ivt, lons, lats, time):

        '''
        FUNCTION INPUTS
        -------------------------------------------------------------------------#
        dataset -- opened netcdf dataset
        fname -- netcdf file of IVT (string) (this is just used as as ref)
        time  -- arraay index correspongding to timestamp (integer)
        -------------------------------------------------------------------------#
        ''' 

        # Enable Garbage Collection
        gc.enable()




        #-------------------------------------------------------------------------#
        # Global Values 
        #-------------------------------------------------------------------------#
        ivt_min = 250                     # Minimum IVT value in kg/ms to be retained
        size_mask = 1000                  # Min Grid cell size of object
        cell_to_km = 50                   # km
        #-------------------------------------------------------------------------#

        # AR logic flag; false initially 
        AR_EXISTS = False
        
        #Let's pass in the opened dataset rather than read each loop 
        ds = dataset 
                
        # Get the date time and convert it to Human Readable
        hr_time      = format_date(ds.variables['time'][time])
        hr_time_str  = hr_time.strftime("%Y-%m-%d_%H")
        
        # Subset of IVT
        ivt   = ivt[time, 0, :, :] 
        # Subset of wind_direction; 
        wnd   = wnd[time, 0, :, :]

        #-------------------------------------------------------------------------#
        # Wind Direction Calculations. 
        # Nested in Try loop since there may be overflow.... 
        #-------------------------------------------------------------------------#
        try:
                # Convert Wind Dir to [0, 359.9]
                wnd_           = wnd/np.pi * 180.0
                wnd_360        = np.where(wnd_ < 0, 360 + wnd_, wnd_)
                
                # u and v components (coordinates of unit vector)
                u_i = np.cos(wnd)
                v_i = np.sin(wnd)

                #Components of IVT 
                u_ivt = ivt*u_i/1000.0
                v_ivt = ivt*v_i/1000.0

        except Exception as e:
                Logger().failure(fname)
                raise e 
                return 


        #-------------------------------------------------------------------------#
        # Label Regions (Blobs)
        # Threshold IVT 
        #-------------------------------------------------------------------------#     
        
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

        # If there is nothing to keep, write to dic and continue
        if not keep_label.any():
		info = 'No ATM Rivers'
		#does no execute the rest of the loop, starts back at top
                return

	#total_blob_size = sum(sizes[sizes>size_mask]) ## why do i need to know this?

        #labeled components to keep 
        blob_num = int(0)

        #-------------------------------------------------------------------------#
        # BEGIN Blob-Loop: Loop through each blob in the list of canditate blobs
        #-------------------------------------------------------------------------#
        from scratch import Blob_Tester
        
        try:
                map(lambda X: Blob_Tester(blob_num,X,label_array,u_i, v_i, wnd_360, ivt,land_mask,lats, lons,lats_mesh, lons_mesh, hr_time_str, fname), keep_label)
                
                
        except Exception as e:
                Logger().failure(fname)
                raise e
        
        
        # Now Make a Plot


#----------------------------------------------------------------------------------#
# FIGURE OUT HOW TO PLOT THINGS HERE
#----------------------------------------------------------------------------------#
#        from plot_ar import make_plot
#        make_plot(lons, lats, zero_arr_1, u_ivt, v_ivt, hr_time_str, save_me=True)

        gc.collect()
        



#----------------------------------------------------------------------------------#
# Runs Script if called directly, i.e. python atmDetect.py
#----------------------------------------------------------------------------------#

if __name__ == '__main__':
        path = '../foo.nc'
        ivt_min = 250                     # Minimum IVT value in kg/ms to be retained
        size_mask  = 1000                  # Min Grid cell size of object
        cell_to_km = 50                   # km
        FindAR_Wrapper(path)



        
        


        
