
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
from skimage.graph import route_through_array
from scipy import ndimage
from datetime import datetime, timedelta
from scipy.spatial import distance


#-------------------------------------------------------------------------#
# TO DO List
# 1. Land mask
# 2. Subtract seasonal means
# 3. Properly calculate grid cell size
# 4. Fix poleward component test; northern needs to be positive, southern negative
# 5. Check on len/width ratio math 
#-------------------------------------------------------------------------#




#-------------------------------------------------------------------------#
# Global Values, Filepaths, etc.
#-------------------------------------------------------------------------#
#path = '/home/wrudisill/scratch/Find_ARs/ivt_files/pgbhnl.gdas.20000201-20000205.nc'
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

def format_date(hours):
        t0 = datetime(1801, 01, 01, 00)
        dt = t0 + timedelta(hours=hours)
        return dt

def least_cost(array, start, end):
	#start, end are tuples of indices (row,col)
	nrow,ncol = array.shape
	indices, cost = route_through_array(array, start, end)
	indices = np.array(indices).T
	path = np.zeros_like(array)
	path[indices[0], indices[1]] = 1
        path_len = len(indices[0])
	return path, cost, path_len 

def distance_matrix(array):
        # Find the two points with the greatest distance from a shape
        # The entry array is a binary array where 1 == object
        # returns the start and end points of shape 
        ind = np.argwhere(array == 1)     # Argwhere returns coords as array ([[j,k], [.,.]...)
        dmat = distance.cdist(ind, ind)   # Scipy euclidian distance matrix
        points = np.argwhere(dmat == dmat.max())     # Coords where there is max distance  
        if len(points) > 0:
                points = points[0]
        start = ind[points[0]]
        end =   ind[points[1]]
        return start, end 

#-------------------------------------------------------------------------#
# Main 
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
# BEGIN Main_Function
#-------------------------------------------------------------------------#

def FindAR(fname, time):
        #  Fname -- netcdf file of IVT (string)
        #  Time  -- arraay index correspongding to timestamp (integer)
     
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
                TC_f = True       # Lanfalling 
                TC_g = False      # NOT crossing equator (we don't want AR crossing eq)
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
                # Blog Lenght Calc -- route throgh array method
                #-------------------------------------------------------------------------#
                try:
                        
                        cost_arr   = np.ones_like(sub_array)*9999.999
                        cost_arr[label_indices] = max(ivt[label_indices]) - ivt[label_indices]
                                
#                        start      = np.array([blob_bbox[0], blob_bbox[1]])
#                        end        = np.array([blob_bbox[2], blob_bbox[3]])
                        start, end = distance_matrix(sub_array)
                        path, cost, path_len = least_cost(cost_arr, start, end)
                        path_len = path_len * cell_to_km
                        
                except Exception:
                        print hr_time_str
                        continue


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
                # Hemisphere -- North or South 
                # Throw out if the object crosses the equator
                #----------------------------------------------------------------------------------#
                
                if any(n < 0 for n in lats[label_indices[0]]) == True:
                        
                        Hemisphere = 'Southern'
                        
                        #if any(n > 0 for n in lats[label_indices[0]]) == True:
                        #        break
                        #else: 
                        #        continue
                else:
                        Hemisphere = 'Northern'
                        TC_f       = True


                #----------------------------------------------------------------------------------#
                # Landfalling 
                #----------------------------------------------------------------------------------#
        

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
                        info = {'AR_Name': AR_Name,
                                'hr_time_str': hr_time_str,
                                'object_orientation_direction': blob_dir, 
                                'object_length': blob_length, 
                                'object_width':  blob_width, 
                                'mean_IVT':  mean_ivt, 
                                'mean_wind_dir': mean_wind_dir,
                                'poleward_IVT':  poleward_IVT,
                                'length_to_width': blob_length_width_ratio,
                                'Hemispere':  Hemisphere,
                                'Path_len' :  path_len,
                                'fname':  fname
                        } 

#                        Results_dictionary[hr_time_str][AR_Name] = info
                        

                        import make_dbase
                        make_dbase.log_AR(**info)


                        
                        #-----------------------------------------------------------------------# 
                        # Create Output Arrays for plotting; fill zero Arrays
                        #-----------------------------------------------------------------------# 
                        zero_arr_0                = zero_arr_0 + path           # Least-cost path
                        zero_arr_1[label_indices] = ivt[label_indices]

                        #Add AR to output Array
                        #new_arr[label_indices] = ivt[label_indices]
                        new_arr[label_indices] = ivt[label_indices]
                        ar_index_list.append(label_indices)
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
                        print 'bar'

        
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

                m.drawparallels(np.arange(-80., 81., 20.), labels=[1,0,0,0], fontsize=5)
#                m.drawmeridians(np.arange(-180., 181., 20.), labels=[0,0,0,1], fontsize=5)
                m.drawcoastlines()
                m.drawstates()
                m.drawcountries()


                #-------------------------------------------------------------------------#
                # 1.b Plot Wind Barbs
                #-------------------------------------------------------------------------#
                
                # Mask out 0 Values using numpy mask

                varmask = np.ma.masked_less(zero_arr_1, 1)
                cs = m.pcolor(xi,yi,varmask,latlon=True, vmin = 200.0, vmax=1000.0)

                varmask0 = np.ma.masked_less(zero_arr_0, 1)
                cs0 = m.pcolor(xi,yi,varmask0,latlon=True)


                # Option: plot with contours
#                clevs = range(0,100,1)                
#                cs = m.contourf(xi,yi,global_cover_class,clevs,cmap='plasma')

                #-------------------------------------------------------------------------#
                # 1.c Plot Wind Barbs
                #-------------------------------------------------------------------------#
                
                # Create a sparse array for plotting barbs to prevent clutter
                yy = np.arange(0, yi.shape[0], 3)
                xx = np.arange(0, xi.shape[1], 3)
                points  = np.meshgrid(yy, xx)                                 
                
#                m.quiver(xi[points], yi[points], u_ivt[points], v_ivt[points], scale = 100, pivot='mid', width =0.001, color='grey') 
                
#                label_ = np.where(new_arr > 0)
#                m.quiver(xi[label_], yi[label_], u_ivt[label_], v_ivt[label_], scale = 100, pivot='mid', width =0.001, color='grey')           


                #-------------------------------------------------------------------------#
                # 1.d Save Figure
                #-------------------------------------------------------------------------#

                # Set Colorbar
                #cbar = m.colorbar(cs,location='bottom',pad="5%")
                #cbar.set_label('mm')
                
                # Format Name string; pad with Zeros
                #if len(str(time)) == 1:
                #        ts = '0'+str(time)
                #else:
#                        ts = str(time)

                # Save Plot
                #plt.title(hr_time)
                #plt.savefig('AR_'+hr_time_str ,format='png', dpi = 700)
                #plt.close()
                

        print "Finished"

#-----------------------------------------------------------------------------------------#
# END Main_Function
#-----------------------------------------------------------------------------------------#







