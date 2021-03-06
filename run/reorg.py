import numpy as np                
from netCDF4 import Dataset
import numpy as np
import glob
from scipy.ndimage import filters, morphology, measurements
from scipy import ndimage
from datetime import datetime, timedelta
import logging
import gc 
import time 
import traceback
from skimage.measure import regionprops
from centerline import Centerline
from make_dbase import make_db
from AR_Object import AR_Object



if 'fig' in locals():  # Delete fig object if it exists 
    del fig

def format_date(hours):
    t0 = datetime(1800, 01, 01, 00)
    dt = t0 + timedelta(hours=hours)
    return dt


def GrepTimeDim(filename):
    ds = Dataset(filename)
    return len(ds.variables['initial_time0_hours'])

def lat_lon_to_indices(lat,lon):
	# lat and lon are geographic coordinates 
    lat = round(lat * 2)/2   # rounds to nearest .5
    lon = round(lon * 2)/2   # rounds to nearest .5   
    x = np.where(lats_mesh == lat)[0]
    y = np.where(lons_mesh == lon)[1]
    return x,y

def write_array(array, name, filepath):
	np.save(filepath+name, array)


def calc_85thp(timestring):
	# accepts a string of time format 'YYYY-mm-dd_hh'
	# weighted averages two monthly files;
	# assumes 15th of the month represents the monthly statistic used
	dt       = datetime.strptime(timestring, '%Y-%m-%d_%H')
	month_a  = dt.month
	day      = dt.day
	# we make the simplifying assumption that months are 30 days long....
	diff = 15. - day  
	#weights 
	wgt1    = (30. - abs(diff))/15.    #wgt for current month
	wgt2    =  2. - wgt1			   #wgt for next(previous) month
        
        #decide which month to use to calc wghted avg
	if diff <= 0:
            if month_a == 1:
                month_b = 12
            else:
                month_b = month_a - 1.
                
	else:
            if month_a == 12:
                month_b = 1
            else:
                month_b = month_a + 1

	# convert month to string to read file
	def zero_pad(mo):
            if mo < 10:
                mostr = '0'+str(int(mo))
            else:
                mostr = str(int(mo))
            return mostr

	path = '/home/wrudisill/scratch/AR-Detection/data/' # path to files
	filea = np.load(path+'Month_'+zero_pad(month_a)+'.npy')
	fileb = np.load(path+'Month_'+zero_pad(month_b)+'.npy')

	p85 = (wgt1*filea + wgt2*fileb)/2.  # weighted average
	return p85   # returns grid


# def blob_dir_correct(blob_orientation):
# 	blob_dir_raw         = blob_orientation / np.pi * 180
# 	blob_dir_corrected   = blob_dir_raw * -1 
# 	# Blob Orientation; we must deal with ambiguity of blob orientation w.r.t to wind dir
# 	# So we give it two possible directions
# 	# Recalls orientation is between -90 and 90. So we convert that to [0,359.9]
# 	if  blob_dir_corrected >= 0:                                        # Positive blob dir
# 		blob_dir = (blob_dir_corrected, 180+blob_dir_corrected)
# 	else:                                                               # Negative blob dir
# 		blob_dir = (360 - abs(blob_dir_corrected), 360 - abs(blob_dir_corrected) - 180)
	# return blob_dir # note: we get two values here 




def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    All args must be of equal length.    
    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2.
    c = 2. * np.arcsin(np.sqrt(a))
    km = 6367. * c
    return km


# x,y = lat_lon_to_indices(-49.029864, -69.162492)
def calc_ivt(filename, time):
	# calculate IVT and other things...
	ds                      = Dataset(filename, format='NETCDF4_CLASSIC')        
	lons                    = ds.variables['lon_0'][:]
	lats                    = ds.variables['lat_0'][:]
	_lons_mesh, lats_mesh   = np.meshgrid(lons, lats)
	lons_mesh               = np.where (_lons_mesh > 180.0, _lons_mesh - 360, _lons_mesh) 
	hr_time                 = format_date(ds.variables['initial_time0_hours'][time])
	hr_time_str             = hr_time.strftime("%Y-%m-%d_%H")
	ivt =  np.zeros((20, 361, 720))
	g = 9.81 
	dp = 25.00 #pa

	# loop through pressure levels, calculate IVT 
	for i in range(20):
	    spfh = ds.variables['SPFH_P0_L100_GLL0'][time,i,:,:]
	    vgrd = ds.variables['VGRD_P0_L100_GLL0'][time,i,:,:]
	    ugrd = ds.variables['UGRD_P0_L100_GLL0'][time,i,:,:]
	    phi =np.arctan2(vgrd,ugrd)
	    tV = (ugrd**2)*(vgrd**2)**.5
	    ivt[i,:,:] = spfh*tV/g*dp

	# sum through pressure levels 
	ivt_integral = np.ndarray.sum(ivt, axis=0)

	# U and V winds respectively
	vgrd =  ds.variables['VGRD_P0_L100_GLL0'][time,:,:,:]
	ugrd =  ds.variables['UGRD_P0_L100_GLL0'][time,:,:,:]
	spfh =  ds.variables['SPFH_P0_L100_GLL0'][time,:,:,:]			
	ds.close()

	output_dictionary = {'ivt':ivt_integral,
						 'vgrd':vgrd,
						 'ugrd':ugrd,
						 'spfh':spfh,
						 'lons_mesh':lons_mesh,
						 'lats_mesh':lats_mesh,
						 'lons':lons,
						 'lats':lats,
						 'hr_time_str':hr_time_str,
						 'filename': filename}
	# return things 
	return output_dictionary


def blob_tester(ivt_timeslice, **kwargs):
	'''
	Test potential AR features and log results.
	'''
	def pacific_region():
		# needs lons_mesh and lats_mesh to exist w/in scope
		# idx = np.where((lons_mesh > -125.) & (lons_mesh < -112.) & (lats_mesh > 30.) & (lats_mesh < 50.))
		idxa = np.where((lons_mesh >= -180.) & (lons_mesh < -120.) & (lats_mesh > 0.) & (lats_mesh < 60.))	
		idxb = np.where((lons_mesh > 120.) & (lons_mesh <= 180.) & (lats_mesh > 0.) & (lats_mesh < 60.))
		foo  = np.zeros_like(lons_mesh) 
		foo[idxa] = 1
		foo[idxb] = 1 
		idx   = np.argwhere(foo == 1 )
		return idx

	def indices_to_lat_lon(index): 
		#input a TUPLE () of grid values 
	    x = lats_mesh[index]   #
	    y = lons_mesh[index]   # 
	    return x,y


	def uv2deg(v,u): # v (vertical), u (horizontal)
	    # by convention, positve v flows from the south to the north
	    # positive u flows from the west to the east 
	    deg = np.arctan2(u,v) * 180./np.pi 
	    deg = np.where(deg < 0, 360.0 + deg, deg)
	    return deg 

	def CircVar(v,u):
		# Circular variance
		c  = np.hypot(v,u)
		vv = np.sum(v/c)/v.size   # equivalent to (SUM cos(*theta*) ) / n 
		uu = np.sum(u/c)/u.size   # equivalent to (SUM sin(*theta*) ) / n 
		r  = np.hypot(uu,vv)	
		dispersion = 1. - r  # should be between zero and one 
		return dispersion

	def ScaleBYq(v,u,q):
		# weighted mean of u,v. flattens along pressure dimension (axis 0)
		v_wgt_mn = np.sum(v*q, axis=0)/np.sum(q, axis=0)
		u_wgt_mn = np.sum(u*q, axis=0)/np.sum(q, axis=0)
		#	    return uv2deg(v_mean,u_mean)
		return v_wgt_mn, u_wgt_mn

	# kwarg options defaults are set 
	ivt_min         = kwargs.get('ivt_min' , 250.0)                 
	size_mask       = kwargs.get('size_mask', 150.0)     		      # Minimum size to retain, in px 
#	min_orientation = kwargs.get('min_orientation',  0.95)     # Minimum orientation (taken as np.abs())
	min_length      = kwargs.get('min_length', 25.) 			      # Shortest feature length (in pixels)
	min_eccentricity= kwargs.get('min_eccentricity', .87)     # Minimum eccentricity of a retained feature


	# unpack from input
	ivt = ivt_timeslice['ivt']
	vgrd = ivt_timeslice['vgrd']
	ugrd = ivt_timeslice['ugrd']
	spfh = ivt_timeslice['spfh']
	lons_mesh   = ivt_timeslice['lons_mesh']
	lats_mesh   = ivt_timeslice['lats_mesh']
  	hr_time_str = ivt_timeslice['hr_time_str'] 

	# -----grid size calcs -----
	earth_rad = 6371.0
	h       = np.cos(np.abs(lats_mesh)*np.pi/180.)*earth_rad  # the radius of the great circle by latitude
	gc      = np.pi*2.*h                                      # the circumference of the greate circle
	grid_dx = gc/720.0        								  # horizontal grid cell distance 
	# calculate vertical grid cell distance; this is the same for all lats 
	grid_dy = earth_rad*2.0*np.pi/720.0
	# ----- grid size calcs -----
   

 	# --- subtract 85th percentile from IVT ----# 

 	# create grid; 85th perc. or 100 kg/m/s, whichever is gt. 	# From Guan et. al 2015 
 	p85 = calc_85thp(hr_time_str)
  	p85 = np.where(p85 > 100., p85, 100.)


	# Convert to binary. 1 is above min value
	threshold_array = np.where(ivt > p85, 1, 0)

	# Labeled components. given values of 1:n
	# Each 'blob' identified is assigned a label 
	label_array, num_labels = ndimage.measurements.label(threshold_array)

	# List of the label values 
	labels = np.array(range(num_labels + 1))

	# List size, in pixels, of components 
	sizes = ndimage.sum(label_array, label_array, range(num_labels + 1))
	   
	# Labels to keep (corresponding w/ objects gt. than size threshold)
	keep_label = labels[sizes > size_mask]

	# Label for AR. Used to name/log for later 	
	OBJECT_ID = 1.0

	# loop through labels 
	for label in keep_label:	
		#--------------------------------------------------#
		#  Test Flags; set to true if passing 
		#--------------------------------------------------#
		tc_ivt          = False      # Mean IVT 
		tc_lgh          = False      # Length
		tc_wdt          = False      # Width ratio
		tc_lwr          = False      # Length/width ratio
		tc_merivt       = False      # meridonal ivt
		tc_ecc          = False      # eccentricity > min 
		AR_BASE_FLAG    = False
		# -------------------------------------------------- #


		# Initialize AR object 
		AR_blob = AR_Object()

		# Indicies of blob of interest; same data, for label_indices_alt each index pair is a (1,2) array
		label_indices     = np.where(label_array == label)
		label_indices_alt = np.argwhere(label_array == label)        # indices in [(a,b), ....] fmt


		#----- test if blob is w/in the Pacific -------# 
		pacific        = map(tuple, list(pacific_region()))	   
		blob_reg       = map(tuple, list(label_indices_alt))
		intersect      = list(set(pacific).intersection(blob_reg))
	
		if len(intersect) > 40:  # At least 40 px need to be in Pacific to be considered 
			pass          # if it is in N pacific, proceed throug rest of loop
		else:
			continue      # otherwise, go to the next blob 


		# -------- Centerline Calculations -------- #		
		# Find the center and related properties of the object 
		center            = Centerline(label_indices_alt, label_indices, ivt, lats_mesh)
		length            = center.path_length
		blob_earth_area   = np.sum(grid_dx[label_indices]*grid_dy)
		width             = blob_earth_area/length           # possible divide by zero
		# ----------------------------------------- #		

		
		# -------- Start and End points -----------# 
		start_lat,start_lon   = indices_to_lat_lon(tuple(center.start))
		end_lat,end_lon       = indices_to_lat_lon(tuple(center.end))
		# decide which is which 
		_lons = np.linspace(-180.,179.5,720.)
		lons = np.where (_lons < 0.0, _lons + 180. , _lons - 180.)

		sl = np.where(lons == float(start_lon)) # get index of point along number line 
		el = np.where(lons == float(end_lon))   # get index of point along number line 
		
		if el < sl:  # end lon is further west than start lon; reverse condition
			start_lat,start_lon   = indices_to_lat_lon(tuple(center.end))
			end_lat,end_lon       = indices_to_lat_lon(tuple(center.start))
		# else:
			# do nothing
		# ----------------------------------------#

		

		# --- Great Cirlce Length beteween start and end---------------#
		gc_distance = haversine(start_lat,start_lon,end_lat,end_lon)



		#-------- WIND DIR AND OBJECT RELATIONSHIP HOOKS GO HERE --------#
		v_wgt_mn, u_wgt_mn = ScaleBYq(vgrd[:, label_indices[0], label_indices[1]], ugrd[:, label_indices[0],label_indices[1]], spfh[:, label_indices[0],label_indices[1]])
		# create a grid to store these values on; to be saved as output
		v_wgt_mn_grd = np.zeros_like(ivt)	
		v_wgt_mn_grd[label_indices] = v_wgt_mn
		u_wgt_mn_grd = np.zeros_like(ivt)             # We create these to save as output 
		u_wgt_mn_grd[label_indices] = u_wgt_mn        # We create these to save as output 
		wind_dir_mean      = uv2deg(np.mean(v_wgt_mn),np.mean(u_wgt_mn))   #coverts u and v wind to a direction in degrees
		wind_dir_var       = CircVar(v_wgt_mn, u_wgt_mn)
		wind_speed         = np.hypot(np.mean(v_wgt_mn),np.mean(u_wgt_mn))

		#### -- calculate meridional component of IVT -----### 
		#   meridional
		# 	   |   /  
 		#      |  /IVT
		#      | / 
		#  phi |/    
		#       --------- zonal

		mean_IVT = np.mean(ivt[label_indices])
		meridional_IVT = mean_IVT*np.cos(wind_dir_mean*np.pi/180.)
		#---------------------------------------------------------------#



		# ----- Test if AR makes landfall ----------- # 
		Lloc = center.landfall_location
		if len(Lloc) > 0:
			Lloc_flag  = 'True'
#			Lloc_pt    = str(map(indices_to_lat_lon, Lloc))
			Lloc_lat,Lloc_lon   = indices_to_lat_lon(Lloc[0])    # maybe there are 2 lflocs... lets just pick one.
			
			#print str(map(indices_to_lat_lon, Lloc))
		else:
			Lloc_flag  = 'False'
			#Lloc_pt    = 'None'
			Lloc_lat = 'None'	
			Lloc_lon = 'None'

		# ----- Remove all elements except blob of interest ------- # 
		# ----- Region Props Algorithm ------------ # 
		sub_array          = np.where(label_array == label, 1, 0)
		blob               = regionprops(sub_array, ivt)[0]  # set to 0; only 1 region        
		#blob_dir_corrected = np.abs(np.abs(blob.orientation/np.pi * 180.)- 90.)
		phi = np.abs(blob.orientation/np.pi * 180.)
                
                if (phi >= 0) and (phi <= 90.):
                    blob_phi = 90. - phi
                
                if (phi <= 0) and (phi >= -90.):
                    blob_phi = 90 + abs(phi)
                
                
		#--------------------------------------------------#
		#  Test Flags; set to true if passing 
		#--------------------------------------------------#
		if mean_IVT > ivt_min:
			tc_ivt     = True      #mean IVT of object

		if length > 2000.:
			tc_lgh     = True      # lenght

		if width < 1000.:         
			tc_wdt     = True      # width

		if length/width > 2.:
			tc_lwr     = True      # lw ratio

		if meridional_IVT > 50.:
			tc_merivt  = True      # lw ratio

		if tc_ivt + tc_lgh + tc_wdt + tc_lwr + tc_merivt == 5:
			AR_BASE_FLAG = True
		# -------------------------------------------------- #



		# ------------------------------------------------------ #
		# ------ Put things in dictionary to write out --------- #
		# ------------------- Metadata ------------------------- # 	

		AR_blob.hr_time_str  			 = ivt_timeslice['hr_time_str']
		AR_blob.OBJECT_ID       		 = OBJECT_ID
		AR_blob.filename        		 = ivt_timeslice['filename']
		AR_blob.landfalling   			 = Lloc_flag
		AR_blob.landfall_lat     		 = Lloc_lat
		AR_blob.landfall_lon     		 = Lloc_lon
		AR_blob.object_length  			 = length
		AR_blob.object_width    	         = width
		AR_blob.length_to_width 	         = length/width    # possible divide by zero
		AR_blob.eccentricity    		 = blob.eccentricity 
		AR_blob.mean_IVT 			 = mean_IVT		#blob.mean_intensity... should be the same but not sure.
		AR_blob.meridional_IVT 		         = meridional_IVT		#blob.mean_intensity... should be the same but not sure.
                AR_blob.object_orientation_direction     = str(blob_phi)
		AR_blob.wind_dir_mean                = str(wind_dir_mean)
		AR_blob.wind_dir_var                 = str(wind_dir_var)
		AR_blob.wind_speed                   = str(wind_speed)
		AR_blob.end_lat                      = str(end_lat)
		AR_blob.end_lon                      = str(end_lon)
		AR_blob.start_lat                    = str(start_lat)
		AR_blob.start_lon                    = str(start_lon)
		AR_blob.gc_distance     	     = str(gc_distance)
		AR_blob.AR_BASE_FLAG                 = str(AR_BASE_FLAG)
   		AR_blob.Make_Db()

		# -------  Create Output Files for Saving ------------ # 
		AR_blob.path = center.path

		if AR_BASE_FLAG == True:
			AR_blob.Save_File(label_indices, ivt, center.path, v_wgt_mn_grd, u_wgt_mn_grd)

		OBJECT_ID += 1.0

if __name__=='__main__':
	filename = '/Users/will/Desktop/AR-Detection/data/pgbhnl.gdas.19960516-19960520.nc'
	print filename
	ivt_timeslice = calc_ivt(filename, 0)
	p85 = blob_tester(ivt_timeslice)


