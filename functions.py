from netCDF4 import Dataset
import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from skimage.graph import route_through_array
from scipy.ndimage import filters, morphology, measurements


def setEndCol0(array):
	nrow, ncol = array.shape
	array[:,0] = 0
	array[:,ncol-1] = 0
	return array

def least_cost(array):
	#start, end are tuples of indices (row,col)
	nrow,ncol = array.shape
	start = (nrow/2,0)
	end = (nrow/2, ncol-1)
	indices, cost = route_through_array(array, start, end)
	indices = np.array(indices).T
	path = np.zeros_like(array)
	path[indices[0], indices[1]] = 1
	return path, cost

def costFX(array, percentile_array):
	return np.where(array > percentile_array, 0, 10)


##thsi does not work correctly, at the moment...
def geoToCoord(lat,lon):
	ncrow = -np.asarray(range(-90,91))
	nccol = np.asarray(range(-180,180))
	pixlat = np.asarray(range(0,181))
	pixlon = np.asarray(range(0,360))
	row = pixlat[np.where(ncrow == lat)]
	col = pixlon[np.where(nccol == lon)]
	return col[0], row[0]

def ReadNc(ncfilename, varname):
	#
	ds = Dataset(ncfilename) # dataset of 1/0 years worth of global cfsr precipitable water, 1 hr every 5 days
	var = ds.variables[varname][:]
	try:
		lon = ds.variables['lon_0'][:] 
		lat = ds.variables['lat_0'][:]
	except KeyError:
		try:
			lon = ds.variables['Lon_0'][:] 
			lat = ds.variables['Lat_0'][:]
		except KeyError:
			print 'check your lat/lon var names'
	ds.close()
	return var, lat, lon


def basePlot(ncfilename, varname):
	#string, string
	# try:
	ncfile = Dataset(ncfilename)
	lons = ncfile.variables['lon_0'][:]
	lats = ncfile.variables['lat_0'][:]
	# except
	# 	print 'lons and or lats not found within ncfile'
	# try:
	lon_0 = -140
	lat_0 = 20
  	#lat_0 = lat_0, lon_0=lon_0
	var = ncfile.variables[varname][:]
	m = Basemap(projection='cyl', resolution='c')
	# ,urcrnrlat=60, urcrnrlon=-80, llcrnrlat=-20, llcrnrlon=-280
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = m(lon, lat)
	m.drawparallels(np.arange(-80., 81., 20.), labels=[1,0,0,0], fontsize=10)
	m.drawmeridians(np.arange(-180., 181., 20.), labels=[0,0,0,1], fontsize=10)
	# # Add Coastlines, States, and Country Boundaries
	m.drawcoastlines()
	m.drawstates()
	m.drawcountries()
	cs = m.pcolor(xi,yi,var,latlon=True)
	cbar = m.colorbar(cs, location='bottom', pad="10%")
	# cbar.set_label(tmax_units)
	plt.title('')
	plt.show()

def countPix(array, percentile_array):
	pix = np.greater(array, percentile_array)
	print len(var[pix])

def subset(array, height, width):
	# height/width are tuples of indices in numpy format
	# change later to be geographic, or have it call a converter function
	return array[height[0]:height[1], width[0]:width[1]]


def unsub(array, orig_array, height, width):
	out_arr = np.zeros_like(orig_array)
	out_arr[height[0]:height[1], width[0]:width[1]] = array
	return out_arr


def subhelper(thing, var, h, w):
	dic[var]=subset(thing[var], h, w)
	

def ishow(array):
	plt.imshow(array)
	plt.colorbar()
	plt.show()

def mknames(x):
	y = len(str(x))
	z = str(x)	
	if y == 1:
		return 'percentile00'+z
	if y == 2:
	  	return 'percentile0'+z
	if y == 3:
		return 'percentile'+z


#this is my attempt at solving the common problem of having a 
#series, and wanting to do something for each n and n+1
#takes image, assigns percentile class for each grid cell
#i.e 5th<X<10th <- value 


def iter(prob_dic, array, h, w):
	klist = prob_dic.keys()
	####I/O
	sub = subset(array, h, w)
	res = np.zeros_like(sub)
	####

	####internal functions 
	def foxx(K):
		#sub and fname come from above
		value = 1 
		out_arr = np.zeros_like(sub)
		A = prob_dic[K]
		out_arr[np.where(sub<A)] = value
		return out_arr

	res = sum(list(map(foxx, klist)))
	return res


def spatial_filter(skeleton,image_array,loops,iterations):
	# image_array = np.where(image_array > threshold, 1, 0)
	im = np.zeros_like(image_array)
	s = morphology.generate_binary_structure(2,2)
	for i in range(loops):
		grow = morphology.binary_dilation(skeleton, iterations = iterations, structure = s)
		skeleton = grow
		im = np.where((grow !=0)&(image_array<7), image_array , -1)
		
	return im


def calcIVT(fname):
	u = 'UGRD_P0_L100_GLL0'
	v = 'VGRD_P0_L100_GLL0'
	iso = 'lv_ISBL0'
	spfh = 'SPFH_P0_L100_GLL0'
	ds = Dataset(fname)
	out = np.zeros_like(ds[u][0])
	pressure = ds.variables[iso][:]

	for j in range(10,37):
		dp = pressure[j] - pressure[j-1]
		tvw = (ds[u][j]**2 + ds[v][j]**2)**.5
		out = out + (tvw * ds[spfh][j] * dp)/9.81

	return out 


def write_result(string):	
	f = open('results.csv','a')
	f.write(string) #Give your csv text here.
	## Python will convert \n to os.linesep
	f.close()







