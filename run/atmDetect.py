
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
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap, cm
from skimage.measure import regionprops
import glob
import json
#from scipy.ndimage.interpolation import rotate
from skimage.morphology import skeletonize
from scipy.ndimage import filters, morphology, measurements
from scipy import ndimage
from datetime import datetime, timedelta
from centerline import Centerline

#-------------------------------------------------------------------------#
# TO DO List
# 1. Land mask
# 2. Subtract seasonal means
# 3. Properly calculate grid cell size
# 4. Fix poleward component test; northern needs to be positive, southern negative
# 5. Check on len/width ratio math 

#-------------------------------------------------------------------------#



#-------------------------------------------------------------------------#
# Helper Functions 
#-------------------------------------------------------------------------#


def smallest_angle(a,b):
        # Returns smalles angle in degrees for [-180, 180]
        tup = (360 - abs(a - b), a-b)
        return min(map(abs, tup))

def format_date(hours):
        t0 = datetime(1801, 01, 01, 00)
        dt = t0 + timedelta(hours=hours)
        return dt


#-------------------------------------------------------------------------#
# Main 
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
# BEGIN Main_Function
#-------------------------------------------------------------------------#

def FindAR(fname, time):
        #-------------------------------------------------------------------------#
        # FUNCTION INPUTS
        #-------------------------------------------------------------------------#
        #  Fname -- netcdf file of IVT (string)
        #  Time  -- arraay index correspongding to timestamp (integer)
        #-------------------------------------------------------------------------#


        #-------------------------------------------------------------------------#
        # Global Values 
        #-------------------------------------------------------------------------#
        ivt_min = 250                     # Minimum IVT value in kg/ms to be retained
        size_mask = 1000                  # Min Grid cell size of object
        cell_to_km = 50                   # km
        #-------------------------------------------------------------------------#

        # AR logic flag; false initially 
        AR_EXISTS = False

        #-------------------------------------------------------------------------#
        # Open netcdf dataset
        # Wind, IVT Field Calcs
        #-------------------------------------------------------------------------#
        ds  = Dataset(fname, format='NETCDF4_CLASSIC')        
        
        # Get the date time and convert it to Human Readable
        hr_time      = format_date(ds.variables['time'][time])
        hr_time_str  = hr_time.strftime("%Y-%m-%d_%H")

        #-------------------------------------------------------------------------#
        # Create Output Dictionary; one for each file 
        #-------------------------------------------------------------------------#

#        Results_dictionary[hr_time_str] = {}

        # Get lat/lon for plotting pu
        lons  = ds.variables['lon_0'][:]
        lats  = ds.variables['lat_0'][:]
        
        # meshgrid
        lons_mesh, lats_mesh  = np.meshgrid(lons, lats)
        
        # Subset of IVT
        ivt   = ds.variables['ivt']
        ivt   = ivt[time, 0, :, :] # Right now we're only looking at one time!!!
        
        # Subset of wind_direction; 
        wnd   = ds.variables['w_dir'] # In Units of Radians
        wnd   = wnd[time, 0, :, :]

        # Convert Wind Dir to [0, 359.9]
        wnd_           = wnd/np.pi * 180.0
        wnd_360        = np.where(wnd_ < 0, 360 + wnd_, wnd_)

        # Landcover 
        global_cover_class  = ds.variables['landcover'][0, 0, :, :]
        land_mask           = np.where(global_cover_class < 1000, 1, 0)
        
        # Close Dataset
        ds.close()

        # u and v components (coordinates of unit vector)
        u_i = np.cos(wnd)
        v_i = np.sin(wnd)

        #Components of IVT 
        u_ivt = ivt*u_i/1000.0
        v_ivt = ivt*v_i/1000.0
        
        #-------------------------------------------------------------------------#
        # Label Regions (Blobs)
        # Threshold IVT 
        #-------------------------------------------------------------------------#     
        
	#subtract mean. TODO: subtract seasonal variation
	#m_out = out - np.mean(out)

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
		#Results_dictionary[fname[41:]] = info 
		#does no execute the rest of the loop, starts back at top
                return

	#total_blob_size = sum(sizes[sizes>size_mask]) ## why do i need to know this?

        #labeled components to keep 
        blob_num = int(0)

        #-------------------------------------------------------------------------#
        # Placeholder Arrays; Fill up with useful things 
        # 
        #-------------------------------------------------------------------------#
        
        zero_arr_0 = np.zeros_like(label_array)        
        zero_arr_1 = np.zeros_like(label_array)        

        #-------------------------------------------------------------------------#
        # BEGIN Blob-Loop: Loop through each blob in the list of canditate blobs
        #-------------------------------------------------------------------------#
        for label in keep_label:               
 
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
                
                
                # Indicies of blob of interest
                label_indices     = np.where(label_array == label)
                label_indices_alt = np.argwhere(label_array == label) 


                # Remove all elements except blob of interest
                sub_array = np.where(label_array == label, 1, 0)

                #-------------------------------------------------------------------------#
                # Region Properties of blob of interest 
                #    
                #    From scipy.morphology.regionprops
                #    Lenght/Width Criteria                
                #    TC_e: Length must be gt. 2000km and l:w gt 2. 
                #-------------------------------------------------------------------------#

                blob          = regionprops(sub_array)[0]                    # set to 0; only 1 region        
                blob_length   = blob.major_axis_length*cell_to_km
                blob_width    = blob.minor_axis_length*cell_to_km
                blob_size     = blob.filled_area*cell_to_km
                blob_bbox     = map(int, blob.bbox)

                #-------------------------------------------------------------------------#
                # Blob Direction Calulations;
                # We must flip blob dir (*-1) since the array is flipped w.r.t lat
                #-------------------------------------------------------------------------#
                blob_dir_raw         = blob.orientation / np.pi * 180
                blob_dir_corrected   = blob_dir_raw * -1 

                # Blob Orientation; we must deal with ambiguity of blob orientation w.r.t to wind dir
                # So we give it two possible directions
                # Recalls: orientation is between -90 and 90. So we convert that to [0,359.9]

                if  blob_dir_corrected >= 0:                                        # Positive blob dir
                        blob_dir = (blob_dir_corrected, 180+blob_dir_corrected)
                else:                                                               # Negative blob dir
                        blob_dir = (360 - abs(blob_dir_corrected), 360 - abs(blob_dir_corrected) - 180)
                        
                
                #-------------------------------------------------------------------------#
                # Blob Width/Length 
                #-------------------------------------------------------------------------#

                # Break loop if there is no width...

                if blob_width == 0: 
                        break
                blob_length_width_ratio = blob_length/blob_width

        
                #-------------------------------------------------------------------------#
                # AR Perimeter Boundary Grid Cells
                #-------------------------------------------------------------------------#
                # scipy.ndimage.morphology.distance_transform_edt
                # outputs euclidian distance to the outer object
                # Use to find AR perimeter


                #-------------------------------------------------------------------------#
                #    Mean IVT Calculation
                #    
                #   TC_a: Mean IVT must be greater than 250 kg/m/s
                #-------------------------------------------------------------------------#

                mean_ivt = np.mean(ivt[label_indices])              # Mean IVT of blob
                

                #-------------------------------------------------------------------------#
                #    Wind Direction calculations
                #
                #    TC_b:  No more than 1/2 of gridcells in blob should deviate more than 45* 
                #           from objects mean IVT direction
                #    TC_c:  Mean wind dir and Object Orientation (blob_dir) Should not deviate by more than 45*
                #    TC_d:  Mean wind dir must have 'significant' poleward (meridonal) component;    
                #    
                #    ----------Mean Wind direction-----------
                #    Calculated by finding the mean of the u and v components respectively, then 
                #    Finding the angle between them and converting to degree coordinates
                #    And finally converting to the range [0, 359.0]
                #--------------------------------------------------------------------

                mean_u_i        = np.mean(u_i[label_indices])
                mean_v_i        = np.mean(v_i[label_indices])
                mean_wind_dir_  = np.arctan2(mean_v_i,mean_u_i)/np.pi * 180
                mean_wind_dir   = np.where(mean_wind_dir_< 0, 360 + mean_wind_dir_ , mean_wind_dir_)

                #    ----------Angular difference between object and mean wind dir-------
                angle_diff = map(lambda X: smallest_angle(X, mean_wind_dir), wnd_360[label_indices])

                # Counts angles gt 45 deg from mean dir
                angle_gt_mean = map(lambda x: x>45, angle_diff).count(True) 

                
                #    ----------Poleward IVT----------
                poleward_IVT = mean_ivt*np.sin(mean_wind_dir/180*np.pi)


        
                #--------------------------------------------------------------------
                # Hemispere
                #--------------------------------------------------------------------
                if any(n < 0 for n in lats[label_indices[0]]) == True:
                        
                        Hemisphere = 'Southern'
                        
                else:
                        Hemisphere = 'Northern'

                #----------------------------------------------------------------------------------#
                # Landfalling 
                #----------------------------------------------------------------------------------#
                
                if any( n == 1 for n in land_mask[label_indices]) == True:
                        Landfalling = True 
                else:
                        Landfalling = False


                
                #------------------------------------------------------------------
                # Does AR make it into interior West? Between 120 W and 100 W 
                #--------------------------------------------------------------------

                if any( (n > 240) and (n < 260 ) for n in lons_mesh[label_indices]) == True:
                        Interior_Landfalling = True

                else:
                        Interior_Landfalling = False



                        
                #-------------------------------------------------------------------------#
                # Blog Length Calc -- route through array method
                # And Landfalling Location 
                #-------------------------------------------------------------------------#
                
                center   = Centerline(label_indices_alt, label_indices, ivt, land_mask)
                
                if center.landfall_location:
                        landfall_point = []
                        for i in center.landfall_location:
                                landfall_point.append((lats_mesh[i],lons_mesh[i]))
                
                else:
                        landfall_point = None




                #----------------------------------------------------------------------------------#
                # Object Testing 
                #----------------------------------------------------------------------------------#

                if mean_ivt > 350.0:
                        TC_a = True

                if angle_gt_mean < len(angle_diff)/2:
                        TC_b = True
                
                if smallest_angle(mean_wind_dir, blob_dir[0]) < 50.0 or smallest_angle(mean_wind_dir, blob_dir[1]) < 50.0:
                        TC_c = True
                
                if Hemisphere == 'Northern':
                        if poleward_IVT > 50.0:     
                                TC_d = True

                elif abs(poleward_IVT) > 50:
                        TC_d = True 
                

                if (blob_length > 2000.0) and (blob_length_width_ratio > 2.0): # Later add a width component...maybe this does not matter
                        TC_e = True
           
                # if Interior_Landfalling == True:
                #        TC_f = True
                
                TC_f = True         
                        
                # Add landfalling later

                
                #---------------------------------------------------------------------------------#
                # Write output to dictionary or do Nothing
                #---------------------------------------------------------------------------------#
                                
                if TC_a + TC_b + TC_c + TC_d + TC_e + TC_f == 6:  # In Python, True == 1 
                        
                        AR_EXISTS = True

                        #-------------------------------------------------------------------------#
                        # Write output to dictionary 
                        # optional: save output plot
                        #-------------------------------------------------------------------------#

                        # ID of AR blob feature
                        blob_num = blob_num + 1                                
                
                        # Name of AR 
                        AR_Name = 'AR_' + str(blob_num)
                        

                        # Dictionary Entry 
                       

                        info =  {
                                'hr_time_str'                      : hr_time_str,
                                'AR_Name'                          : AR_Name,
                                'fname'                            : fname,
                                'object_length'                    : float(blob_length),
                                'object_width'                     : float(blob_width),
                                'mean_IVT'                         : float(mean_ivt),
                                'mean_wind_dir'                    : float(mean_wind_dir),
                                'poleward_IVT'                     : float(poleward_IVT),
                                'object_orientation_direction_a'   : float(blob_dir[0]),
                                'object_orientation_direction_b'   : float(blob_dir[1]),
                                'length_to_width'                  : float(blob_length_width_ratio),
                                'Hemispere'                        : str(Hemisphere),       
                                'Path_len'                         : float(center.path_len),
                                'ID_Landfalling'                   : str(Interior_Landfalling),
                                'Landfalling'                      : str(Landfalling),
                                'Landfall_point'                   : str(landfall_point)
                        }
                        


                        #-----------------------------------------------------------------------# 
                        # import make_dbase function;
                        # log AR info to sqli
                        #-----------------------------------------------------------------------# 
                        import make_dbase
                        make_dbase.make_db(**info)



                        #-----------------------------------------------------------------------# 
                        # Create Output Arrays for plotting; fill zero Arrays
                        #-----------------------------------------------------------------------# 
                        zero_arr_0                = zero_arr_0 + center.path           # Least-cost path
                        zero_arr_1[label_indices] = ivt[label_indices]
                        
                else:
                        continue

        #----------------------------------------------------------------------------------------#
        # END Blob-Loop
        #---------------------------------------------------------------------------------#
#        lon, lat = np.meshgrid(lons, lats)

        
        #---------------------------------------------------------------------------------#
        #  If there are AR:
        #    a. Print out summary statistics
        #    b. Plot the ARs w/ basemps
        #    c. In progress --- write AR array out to a NCfile
        #---------------------------------------------------------------------------------#

#        if AR_EXISTS == True:
#                from plot_ar import make_plot
#                make_plot(lons, lats, zero_arr_1, u_ivt, v_ivt, hr_time_str, save_me=True)
#                print "Finished"
                

#-----------------------------------------------------------------------------------------#
# END Main_Function
#-----------------------------------------------------------------------------------------#



# Runs Script if called directly, i.e. python atmDetect.py
if __name__ == '__main__':
#        path = '/home/wrudisill/scratch/Find_ARs/ivt_files/pgbhnl.gdas.20000201-20000205.nc'
        path = 'foo.nc'
        ivt_min = 250                     # Minimum IVT value in kg/ms to be retained
        size_mask  = 1000                  # Min Grid cell size of object
        cell_to_km = 50                   # km
        for i in range(20):
                FindAR(path, i)



        
        


        
