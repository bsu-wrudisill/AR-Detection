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
from mpl_toolkits.basemap import Basemap, cm
from skimage.measure import regionprops
import glob
import json
#from scipy.ndimage.interpolation import rotate
from skimage.morphology import skeletonize
from scipy.ndimage import filters, morphology, measurements
#from skimage.graph import route_through_array
from scipy import ndimage

#-------------------------------------------------------------------------#
# TO DO List
# 1. Land mask
# 2. Subtract seasonal means
# 3. Properly calculate grid cell size
# 4. Fix poleward component test; northern needs to be positive, southern negative
#-------------------------------------------------------------------------#




#-------------------------------------------------------------------------#
# Global Values, Filepaths, etc.
#-------------------------------------------------------------------------#

path = '/home/wrudisill/scratch/Find_ARs/sample/20101226-20101231_IVT.nc'

ivt_min = 250                     # Minimum IVT value in kg/ms to be retained
size_mask = 1000                  # Min Grid cell size of object
cell_to_km = 50                   # km
Results_dictionary = {}           # output dictionary


#-------------------------------------------------------------------------#
# Helper Functions 
#-------------------------------------------------------------------------#


def smallest_angle(a,b):
        # Returns smalles angle in degrees for [-180, 180]
        tup = (360 - abs(a - b), a-b)
        return min(map(abs, tup))


#-------------------------------------------------------------------------#
# Main 
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
# BEGIN Main_Function
#-------------------------------------------------------------------------#

def FindAR(fname, time):
        #  Fname -- netcdf file of IVT (string)
        #  Time  -- arraay index correspongding to timestamp (integer)
     
        # Output Dictionary
        Results_dictionary[fname[40:]] = {}
        AR_EXISTS = False

        #-------------------------------------------------------------------------#
        # Open netcdf dataset
        # Wind, IVT Field Calcs
        #-------------------------------------------------------------------------#
        ds    = Dataset(fname, format='NETCDF4_CLASSIC')
        
        # Get lat/lon for plotting pu
        lons  = ds.variables['longitude'][:]
        lats  = ds.variables['latitude'][:]
        
        # Subset of IVT
        ivt   = ds.variables['ivt']
        ivt   = ivt[time, 0, :, :] # Right now we're only looking at one time!!!
        
        # Subset of wind_direction; 
        wnd   = ds.variables['w_dir'] # In Units of Radians
        wnd   = wnd[0, 0, :, :]

        # Convert Wind Dir to [0, 359.9]
        wnd_           = wnd/np.pi * 180.0
        wnd_360        = np.where(wnd_ < 0, 360 + wnd_, wnd_)

        # Landcover 
        global_cover_class  = ds.variables['landcover'][0, 0, :, :]
        land_mask           = np.where(global_cover_class > 0)
        water_mask          = np.where(global_cover_class == 0)

        # Close Dataset
        ds.close()

        # u and v components (coordinates of unit vector)
        u_i = np.cos(wnd)
        v_i = np.sin(wnd)

        # IVT magnitude
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
		Results_dictionary[fname[41:]] = info 
		#does no execute the rest of the loop, starts back at top
                return

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
                TC_b = False      # Coherency in IVT Direction  (variance)
                TC_c = False      # Object Mean Meridonal IVT (flowing towards a pole)
                TC_d = False       # Object/IVT direction Consistency
                TC_e = False      # Length/Width Ratio
                TC_f = True       # Lanfalling 
                #-------------------------------------------------------------------------#
                

                # Indicies of blob of interest
                label_indices = np.where(label_array == label)

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
                        blob_dir = (360 - abs(blob_dir_corrected), 90 + abs(blob_dir_corrected))
        
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
                #-------------------------------------------------------------------------#
                
                # Mean Wind direction
                # Calculated by finding the mean of the u and v components respectively, then 
                # Finding the angle between them and converting to degree coordinates
                # And finally converting to the range [0, 359.0]
                
                mean_u_i        = np.mean(u_i[label_indices])
                mean_v_i        = np.mean(v_i[label_indices])
                mean_wind_dir_  = np.arctan2(mean_v_i,mean_u_i)/np.pi * 180
                mean_wind_dir   = np.where(mean_wind_dir_< 0, 360 + mean_wind_dir_ , mean_wind_dir_)

                # Angular difference between object and mean wind dir
                angle_diff = map(lambda X: smallest_angle(X, mean_wind_dir), wnd_360[label_indices])
                angle_gt_mean = map(lambda x: x>45, angle_diff).count(True) # Counts angles gt 45 deg from mean dir
                
                # Poleward IVT 
                poleward_IVT = mean_ivt*np.sin(mean_wind_dir/180*np.pi)

                #----------------------------------------------------------------------------------#
                # Landfalling 
                #----------------------------------------------------------------------------------#
        

                #----------------------------------------------------------------------------------#
                # Object Testing 
                #----------------------------------------------------------------------------------#

                if mean_ivt > 250.0:
                        TC_a = True

                if angle_gt_mean < len(angle_diff)/2:
                        TC_b = True
                
                if smallest_angle(mean_wind_dir, blob_dir[0]) < 50.0 or smallest_angle(mean_wind_dir, blob_dir[1]) < 50.0:
                        TC_c = True
                
                if abs(poleward_IVT) > 50.0:      # Should southern/northern hemispher require south/north directed flow?
                        TC_d = True
                
                if (blob_length > 2000.0) and (blob_length_width_ratio > 2.0): # Later add a width component...maybe this does not matter
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
                        #-------------------------------------------------------------------------#

                        # ID of AR blob feature
                        blob_num = blob_num + 1                                
                        # Name of AR 
                        AR_Name = 'AR_' + str(blob_num)
 

                        # Dictionary Entry 
                        info = {'object_orientation_direction': blob_dir, 
                                'object_length': blob_length, 
                                'object_width':blob_width, 
                                'mean_IVT': mean_ivt, 
                                'mean_wind_dir': mean_wind_dir,
                                'poleward_IVT': poleward_IVT,
                                'lenght_to_width':blob_length_width_ratio,
                        } 

                        Results_dictionary[fname[40:]][AR_Name] = info
                        
                        #Add AR to output Array
                        #new_arr[label_indices] = ivt[label_indices]
                        new_arr[label_indices] = ivt[label_indices]
                        
#                else:
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

                m = Basemap(projection='cyl', 
                            resolution='c', 
                            llcrnrlat= -90,
                            urcrnrlat= 90,
                            llcrnrlon= 0,
                            urcrnrlon= 360)

                lons = m.shiftdata(lons, lon_0=180)
                lon, lat = np.meshgrid(lons, lats)
                xi, yi = m(lon, lat)

                #-------------------------------------------------------------------------#
                # 1.a Draw State/country Borders
                #-------------------------------------------------------------------------#

#                m.drawparallels(np.arange(-80., 81., 20.), labels=[1,0,0,0], fontsize=5)
#                m.drawmeridians(np.arange(-180., 181., 20.), labels=[0,0,0,1], fontsize=5)
                m.drawcoastlines()
                m.drawstates()
                m.drawcountries()


                #-------------------------------------------------------------------------#
                # 1.b Plot Wind Barbs
                #-------------------------------------------------------------------------#
                
                # Mask out 0 Values using numpy mask
                varmask = np.ma.masked_less(ivt, 250)
                cs = m.pcolor(xi,yi,varmask,latlon=True)


                # Option: plot with contours
                clevs = range(0,100,1)                
                cs = m.contourf(xi,yi,global_cover_class,clevs,cmap='plasma')

                #-------------------------------------------------------------------------#
                # 1.c Plot Wind Barbs
                #-------------------------------------------------------------------------#
                
                # Create a sparse array for plotting barbs to prevent clutter
                yy = np.arange(0, yi.shape[0], 3)
                xx = np.arange(0, xi.shape[1], 3)
                points = np.meshgrid(yy, xx) 
                m.quiver(xi[points], yi[points], u_ivt[points], v_ivt[points], scale = 100, pivot='mid', width =0.001, color='grey') 



                #-------------------------------------------------------------------------#
                # 1.d Save Figure
                #-------------------------------------------------------------------------#

                # Set Colorbar
                cbar = m.colorbar(cs,location='bottom',pad="5%")
                cbar.set_label('mm')
                
                # Format Name string; pad with Zeros
                if len(str(time)) == 1:
                        ts = '0'+str(time)
                else:
                        ts = str(time)

                # Save Plots
                plt.savefig('testmap' + ts ,format='png', dpi = 700)
                plt.close()
                
                return varmask
        print "Finished"

#-----------------------------------------------------------------------------------------#
# END Main_Function
#-----------------------------------------------------------------------------------------#


for i in range(0,1):
        var = FindAR(path, i)
        print 'done with %s' %i





#with open('results.json', 'w') as outfile:
#        json.dump(Results_dictionary, outfile)
#            blob_skeleton = skeletonize(sub_array) # Perhaps use this to calc the lengths...


