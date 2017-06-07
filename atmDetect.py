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
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from skimage.measure import regionprops
import glob
import json
#from scipy.ndimage.interpolation import rotate
from skimage.morphology import skeletonize
from scipy.ndimage import filters, morphology, measurements
#from skimage.graph import route_through_array
from scipy import ndimage


#-------------------------------------------------------------------------#
# Global Values, Filepaths, etc.
#-------------------------------------------------------------------------#

path = "/home/wrudisill/scratch/Find_ARs/sample/20101226-20101231_IVT.nc"
ivt_min = 250                     # Minimum IVT value in kg/ms to be retained
size_mask = 1000                  # Min Grid cell size of object
cell_to_km = 50                   # km
Results_dictionary = {}           # output dictionary


#-------------------------------------------------------------------------#
# Helper Functions 
#-------------------------------------------------------------------------#

def subtract_angle(a,b):
        # Returns the smallest angle between two angles between [-180, 180]
        tup = (360 - (a - b), a-b)
        return abs(min(tup))



#-------------------------------------------------------------------------#
# BEGIN File-Loop:Loop through nc files within path
#-------------------------------------------------------------------------#

for fname in glob.glob(path):

        # Name of stuff
        Results_dictionary[fname[40:]] = {}
        AR_EXISTS = False

        # Open netcdf dataset
        ds = Dataset(fname, format='NETCDF4_CLASSIC')

	# subset of IVT
        ivt =ds.variables['ivt']
        ivt = ivt[0, 0, :, :] # Right now we're only looking at one time!!!
        
        # subset of wind_direction; 
        wnd = ds.variables['w_dir']
        wnd = wnd[0, 0, :, :]
                
	#subtract mean. TODO: subtract seasonal variation
	#m_out = out - np.mean(out)

	#convert to binary. 1 is above min value
	threshold_array = np.where(ivt > ivt_min, 1, 0)
        
	#labeled components. given values of 1:n
	label_array, num_labels = ndimage.measurements.label(threshold_array)
        
        #each 'blob' identified is assigned a label 

	#list of the label values 
	labels = np.array(range(num_labels + 1))

	#list size, in pixels, of components 
	sizes = ndimage.sum(label_array, label_array, range(num_labels + 1))
       
	#labels to keep (corresponding w/ objects gt. than size threshold)
	keep_label = labels[sizes > size_mask]

        #if there is nothing to keep, write to dic and continue
	
        if not keep_label.any():
		info = 'No ATM Rivers'
		Results_dictionary[fname[41:]] = info 
		#does no execute the rest of the loop, starts back at top
                continue

	#total_blob_size = sum(sizes[sizes>size_mask]) ## why do i need to know this?

        #labeled components to keep 
        blob_num = int(0)
        
        # Array of Zeros; fill up w/ AR that get ID'd
        new_arr = np.zeros_like(label_array)        

        #-------------------------------------------------------------------------#
        # BEGIN Blob-Loop: Loop through each blob in the list of canditate blobs
        #-------------------------------------------------------------------------#
        for label in keep_label:               
 
                #-------------------------------------------------------------------------#
                #  Test Flags; set to true if passing 
                #-------------------------------------------------------------------------#
                TC_a = False      # Mean IVT 
                TC_b = False      # Coherency in IVT Direction (i.e. variance)
                TC_c = False      # Object Mean Meridonal IVT 
                TC_d = False      # Object/IVT direction Consistency
                TC_e = False      # Length/Width Ratio
                TC_f = True       # Lanfalling 
                #-------------------------------------------------------------------------#
                

                # Indicies of blob of interest
                label_indices = np.where(label_array == label)

                # Remove all elements except blob of interest
                sub_array = np.where(label_array == label, 1, 0)

                #-------------------------------------------------------------------------#
                # Region Properties of blob of interest
                # From scipy.morphology.regionprops
                #-------------------------------------------------------------------------#
                
                blob = regionprops(sub_array)[0]                    # set to 0, since there is only one region

                #-------------------------------------------------------------------------#
                # Lenght/Width Criteria
                # TC_e: Length must be gt. 2000km, width gt. 
                #-------------------------------------------------------------------------#                
                blob_length = blob.major_axis_length*cell_to_km
                blob_width = blob.minor_axis_length*cell_to_km

                #-------------------------------------------------------------------------#
                # Mean IVT Calculation
                # TC_a: Mean IVT must be greater than 250 kg/m/s
                #-------------------------------------------------------------------------#
                mean_ivt = np.mean(ivt[label_indices])              # Mean IVT of blob
                
                #-------------------------------------------------------------------------#
                # Wind Direction calculations
                # TC_b:  No more than 1/2 of gridcells in blob should deviate more than 45* from objects mean IVT direction
                # TC_c:  Mean wind dir and Object Orientation (blob_dir) Should not deviate by more than 45*
                # TC_d:  Mean wind dir must have 'significant' poleward (meridonal) component;    
                #-------------------------------------------------------------------------#
                
                # Blob Orientation
                blob_dir = blob.orientation/np.pi*180               
                mean_wind_dir = np.mean(wnd[label_indices])         

                # Angular difference between object and mean wind dir

                angle_diff = map(lambda X: subtract_angle(X, mean_wind_dir), wnd[label_indices])
                angle_gt_mean = map(lambda x: x>45, angle_diff).count(True) # Counts angles gt 45 deg from mean dir
        

                #----------------------------------------------------------------------------------#
                # Object Testing 
                #----------------------------------------------------------------------------------#
                if mean_ivt > 250.0:
                        TC_a = True

                if angle_gt_mean < len(angle_diff)/2:
                        TC_b = True
                
                if subtract_angle(mean_wind_dir, blob_dir) < 45.0: 
                        TC_c = True
                
                if (mean_wind_dir > 0) & (mean_wind_dir < 90.0):
                        TC_d = True
                
                if blob_length > 2000.0:           # Later add a width component...maybe this does not matter
                        TC_e = True
           
                
                # Add landfalling later
                
                #---------------------------------------------------------------------------------#
                # Write output to dictionary or do Nothing
                #---------------------------------------------------------------------------------#
                                
                if TC_a + TC_b + TC_c + TC_d + TC_e == 5:  # In Python, True == 1 
                        
                        AR_EXISTS = True
                        
                        #-------------------------------------------------------------------------#
                        # Write output to dictionary 
                        # optional: save output plot
                        # 
                        #-------------------------------------------------------------------------#
                                
                        # Name of AR 
                        AR_Name = 'AR_' + str(blob_num)
                        blob_num = blob_num + 1 

                        # Dictionary Entry 
                        info = {'object_orientation_direction': blob_dir, 
                                'object_length': blob_length, 
                                'object_width':blob_width, 
                                'mean_IVT': mean_ivt, 
                                'mean_IVT_dir': mean_wind_dir}

                        Results_dictionary[fname[40:]][AR_Name] = info
                        
                        # Add AR to output Array
                        new_arr[label_indices] = blob_num * 100
                        
                else:
                        continue

        #----------------------------------------------------------------------------------------#
        # END Blob-Loop
        #---------------------------------------------------------------------------------#
                
        #---------------------------------------------------------------------------------#
        #  If there are AR:
        #    a. Print out summary statistics
        #    b. Plot the ARs w/ basemps
        #    c. In progress --- write AR array out to a NCfile
        #---------------------------------------------------------------------------------#

        if AR_EXISTS == True:
                print Results_dictionary
                
                #-------------------------------------------------------------------------#
                # Get lat/lon for plotting pu
                lons = ds.variables['longitude'][:]
                lats = ds.variables['latitude'][:]
                ds.close() # close Netcdf File 
                #-------------------------------------------------------------------------#


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
                # Plot with Basemap
                #-------------------------------------------------------------------------#
                m = Basemap(projection='cyl', resolution='c')
                lon, lat = np.meshgrid(lons, lats)
                xi, yi = m(lon, lat)
                m.drawparallels(np.arange(-80., 81., 20.), labels=[1,0,0,0], fontsize=5)
                m.drawmeridians(np.arange(-180., 181., 20.), labels=[0,0,0,1], fontsize=5)
                m.drawcoastlines()
                m.drawstates()
                m.drawcountries()
                cs = m.pcolor(xi,yi,new_arr,latlon=True)
                plt.savefig('testmap',format='png', dpi = 1000)
                plt.close()
                
#-----------------------------------------------------------------------------------------#
# END File-Loop: Done looping through files in path
#-----------------------------------------------------------------------------------------#

print "Finished"



#with open('results.json', 'w') as outfile:
#        json.dump(Results_dictionary, outfile)
#            blob_skeleton = skeletonize(sub_array) # Perhaps use this to calc the lengths...


