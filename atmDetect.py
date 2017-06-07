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


# fname = "/Users/will/Desktop/pgbf2010053100.01.2010053100.nc"
path = "/home/wrudisill/scratch/Find_ARs/sample/20101226-20101231_IVT.nc"

###globals
ivt_min = 250 
size_mask = 1000
cell_to_km = 50 #km

#out 
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
        
        #-------------------------------------------------------------------------#
        # 1.0 Loop through each blob in the list of canditate blobs
        #   a. Get shape statistics (length, width, orientation)
        #   b. Find mean IVT w/in the blob region
        #   c. Find mean wind direction w/in blob region
        #   d. Find blob skeleton shape (maybe this isn't important...)
        #   e. More to come.... 
        #-------------------------------------------------------------------------#

        new_arr = np.zeros_like(label_array)

        for label in keep_label:                
                label_indices = np.where(label_array == label)
                sub_array = np.where(label_array == label, 1, 0)
                
                #properties of AR components
                blob_props = regionprops(sub_array)
                blob_num = blob_num + 1
                blob_skeleton = skeletonize(sub_array)
                
                # 
                mean_ivt = np.mean(ivt[label_indices])
                mean_wind_dir = np.mean(wnd[label_indices])
                
                # distance matrix                                
                blob = blob_props[0]
                blob_dir = blob.orientation
                blob_lab = 'AR_Blob'+str(blob_num)
                info = {'orientation': blob.orientation , 'length': blob.major_axis_length*cell_to_km, 'width':blob.minor_axis_length*cell_to_km, 'mean_ivt' : mean_ivt, 'mean_wind_dir' : mean_wind_dir}
                Results_dictionary[fname[41:]][blob_lab] = info
                
                new_arr[label_indices] = blob.orientation
                

                
        plt.imshow(new_arr)
        plt.colorbar()
        plt.show()

        
print Results_dictionary

#with open('results.json', 'w') as outfile:
#        json.dump(Results_dictionary, outfile)






# write_out = orientation + ',' + length + ',' + width

# f = open('results.csv','w')
# f.write('\n') #Give your csv text here.
# ## Python will convert \n to os.linesep
# f.close()



# sk = skeletonize(np.where(m_out>ivt_min, 1, 0))
# ir,ic = np.nonzero(sk)

# rl = max(ir) - min(ir)
# cl = max(ic) - min(ic) 

# if rl > cl:
# 	ar_flag = "vertical"
# elif rl < cl:
# 	ar_flag = "horizontal"


# duplicates = []
# newr,newc = [],[]

# for ind,row in enumerate(ir):
#   if np.count_nonzero(np.where(ir == ir[ind])) > 1 and row not in duplicates:
#     duplicates.append(row)
#     newr.append(ir[ind])
#     newc.append(np.min(ic[np.where(ir == ir[ind])]))
#   elif row in duplicates:
#     continue
#   else:
#     newr.append(ir[ind])
#     newc.append(ic[ind])

# snewr, snewc = (list(el) for el in zip(*sorted(zip(newr, newc))))
# z = np.polyfit(snewr,snewc,3)
# p = np.poly1d(z)
# y = np.arange(np.min(snewr),np.max(snewr)+1,1)
# # y = np.where(y >= bounds[3], bounds[3]-0.1, y)
# x= p(y)

# bounds = (67., 5., 252., 118.) # As: N, S, E, W
# lats = np.arange(-90, 90.5, 0.5)
# lons = np.arange(30, 390, 2./3)
# ur_lat_ind,ur_lat = min(enumerate(lats), key=lambda i: abs(i[1]-bounds[0]))
# ur_lon_ind,ur_lon = min(enumerate(lons), key=lambda i: abs(i[1]-bounds[2]))
# ll_lat_ind,ll_lat = min(enumerate(lats), key=lambda i: abs(i[1]-bounds[1]))
# ll_lon_ind,ll_lon = min(enumerate(lons), key=lambda i: abs(i[1]-bounds[3]))
# sub_lats = lats[ll_lat_ind:ur_lat_ind]
# sub_lons = lons[ll_lon_ind:ur_lon_ind]
# nrows = len(sub_lats)
# ncols = len(sub_lons)


# y_vals = np.copy(y)
# y_pixel_km = (1./2)*(111.2) # Approx. pixel extent in km
# piecewise_length = []

# for i in xrange(0,len(x)-1,1):
#     # Visiting Pythagoras to estimate length of skeleton in km
#     x_pixel_km = (2./3)*(np.pi/180.0)*6378.137*abs(np.cos((np.pi/180.0)*sub_lats[y_vals[i]]))
#     len_km = np.sqrt((abs(y_vals[i+1]-y_vals[i])*x_pixel_km)**2+(abs(x[i+1]-x[i])*y_pixel_km)**2)
#     piecewise_length.append(len_km)





# cost = abs(out-np.max(out))


# ends0 = setEndCol0(cost)
# l_cost = least_cost(ends0)
# # # ishow(l_cost[0]*50 + cost_out)
# out = spatial_filter(l_cost[0], cost, 4, 1)
# # foo = unsub(np.where(out>-1, 1, 0),  pwat, height, width)

# ishow(np.where(foo==1, 9, pwat))









####I/O

####Step1: Read prob File
# percentile_file = '/Users/will/Desktop/Atm-Riv/ProbabilityPWAT.nc'	
# ds = Dataset(percentile_file).variables
# ks = ds.keys()[2:21]
# dic = {}
# #subsets prob array, stores in dic
# map(lambda V: subhelper(ds, V, h, w), ks)

# ####Read in files
# #atm precip water 
# fname = ('/Users/will/Desktop/Atm-Riv/figures/pgbf2010060500.01.2010053100.nc', 'PWAT_P0_2L108_GLL0')
# pwat, lat ,lon = ReadNc(fname[0], fname[1])
# foo = iter(dic, pwat, h, w)
# ####
# cost_out =  iter(fi, pwat)
