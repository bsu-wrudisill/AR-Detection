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


### Global Values 

path = "/home/wrudisill/scratch/Find_ARs/sample/20101226-20101231_IVT.nc"
ivt_min = 250                     # Minimum IVT value in kg/ms to be retained
size_mask = 1000                  # Min Grid cell size of object
cell_to_km = 50                   # km

### Helper Functions

def subtract_angle(a,b):
        # Returns the smallest angle between two angles between [-180, 180]
        tup = (360 - (a - b), a-b)
        return abs(min(tup))


# Outile 
Results_dictionary = {}

for fname in glob.glob(path):

	Results_dictionary[fname[41:]] = {}
        
        # open netcdf dataset
        ds = Dataset(fname)

	# subset of IVT
        ivt =ds.variables['ivt']
        ivt = ivt[0, 0, :, :]
        ivt[180:360, :] = 0
        
        # subset of wind_direction; 
        wnd = ds.variables['w_dir']
        wnd = wnd[0, 0, :, :]
        wnd[180:360, :] = 0
        
        
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
        new_arr = np.zeros_like(label_array)        

        #-------------------------------------------------------------------------#
        # 1.0 Loop through each blob in the list of canditate blobs
        #-------------------------------------------------------------------------#

        for label in keep_label:               
 
                #-------------------------------------------------------------------------#
                #  Test Flags; set to true if passing 
                #-------------------------------------------------------------------------#
                TC_a = False      # Mean IVT 
                TC_b = False      # Coherency in IVT Direction
                TC_c = False      # Object Mean Meridonal IVT 
                TC_d = False      # Object/IVT direction Consistency
                TC_e = False      # Length/Width Ratio
                TC_f = True      # Lanfalling 
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
                # TC_d:  Mean wind dir must have 'significant' poleward (meridonal) component;    
                # TC_c:  Mean wind dir and Object Orientation (blob_dir) Should not deviate by more than 45*
                #-------------------------------------------------------------------------#
                
                # Blob Orientation
                blob_dir = blob.orientation/np.pi*180               
                mean_wind_dir = np.mean(wnd[label_indices])         

                # Angular difference between object and mean wind dir

                angle_diff = map(lambda X: subtract_angle(X, mean_wind_dir), wnd[label_indices])
                angle_gt_mean = map(lambda x: x>45, angle_diff).count(True) # Counts angles gt 45 deg from mean dir
        

                #-------------------------------------------------------------------------#
                # Object Testing 
                #-------------------------------------------------------------------------#
                if mean_ivt > 250.0:
                        TC_a = True

                if angles_gt_mean < len(angle_diff)/2:
                        TC_b = True
                
                if subtract_angle(mean_wind_dir, blob_dir) < 45.0: 
                        TC_c = True
                
                if (mean_wind_dir > 0) & (mean_wind_dir < 90.0):
                        TC_d = True
                
                if blob_length > 2000.0:           # Later add a width component...maybe this does not matter
                        TC_e = True
           
                
                #-------------------------------------------------------------------------#
                # Write output to dictionary 
                # 
                # 
                #-------------------------------------------------------------------------#
                
                blob_lab = 'AR_Blob'+str(blob_num)
                info = {'orientation': blob.orientation/np.pi*180 , 'length': blob.major_axis_length*cell_to_km, 'width':blob.minor_axis_length*cell_to_km, 'mean_ivt' : mean_ivt, 'mean_wind_dir' : mean_wind_dir}
                Results_dictionary[fname[41:]][blob_lab] = info
                
                #new_arr[label_indices] = blob.orientation
                blob_num = blob_num + 1 
        
                
#       plt.imshow(new_arr)
#       plt.colorbar()
#       plt.show()

        
print Results_dictionary



#with open('results.json', 'w') as outfile:
#        json.dump(Results_dictionary, outfile)

#            blob_skeleton = skeletonize(sub_array) # Perhaps use this to calc the lengths...



