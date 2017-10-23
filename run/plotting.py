import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import rcParams
import glob
import pandas as pd
import sqlite3


################ READ DATABASE ################
db = 'Atmospheric_River.db'
conn = sqlite3.connect(db)
cur  = conn.cursor()

# 'Rivers' is the name of the table 
# Command fetches the column names of the sqlite db 
cur.execute("PRAGMA table_info('Rivers')")
clist = cur.fetchall()


# function to pull the 
def PullCol(col, d, **kwargs):
	# change this so that we can stick in more SQL conditions 
	kwargs.get('condition', None)
	cur.execute("SELECT "+col+" FROM "+d)# WHERE landfalling='True'")  
	field = cur.fetchall()
	field = map( lambda X: reduce(list, X), field)
	return field

# Create a dictionary object to pass things into 
dic = {}
for i in clist:
	dic[str(i[1])] = PullCol(str(i[1]), 'Rivers')


# Create a pandas dataframe from the dictionary we have created
# I know already that my date column is called 'hr_time_str';
df = pd.DataFrame(data=dic)
df.set_index(pd.to_datetime(df['hr_time_str'], format='%Y-%m-%d_%H'), inplace=True)
del df['hr_time_str']

# Count up things 
df['Year'] = df.index.year
df['Month'] = df.index.month
df['Count'] = 1


#image = np.load('../data/outfiles/1996-05-16_18_OBJECT_1.0.npy')
#ilist = glob.glob('../data/outfiles/*npy')



def LocateRecord(date, ID):
	# return a row of pd dataframe by date and object ID
	# date 'yyyy-mm-dd hh:mm:ss'
	# id   '1.0' , etc.
	record = df[(df.index == date) & (df.OBJECT_ID == ID)]
	return record 


def PlotAR(record):
	# INPUT: 1 row of pandas df. returned by record locator 

	############################# 	############################# 
	# parse record to get date-time, other things 
	# Create filename string
	time_str = str(record.index.values[0])[0:10]
	hr       = str(record.index.values[0])[11:13]
	ID       = str(record.OBJECT_ID.values[0])
	direc    = '../data/outfiles/'
	fname    = time_str + '_' + hr + '_OBJECT_' + ID + '.npy'

	############################# 	############################# 

	# construct text box for plotting 
	mean_IVT       = str(np.round(float(record.mean_IVT.values[0]), 2))
	object_length  = str(np.round(float(record.object_length.values[0]), 2))
	wind_dir_mean  = str(np.round(float(record.wind_dir_mean.values[0]), 2))
	textstr = 'mean_IVT: ' + mean_IVT + '\n' + 'object_length: ' + object_length + '\n' + 'wind_dir_mean: ' + wind_dir_mean 



	############################# 	############################# 
	# Create Plot 

	image = np.load(direc+fname)
	land  = np.load('../data/land.npy') # Check me 
	v    = ma.masked_equal(image[1], 0)
	u    = ma.masked_equal(image[2], 0)
	ivt  = ma.masked_equal(image[0], 0)
	path = ma.masked_equal(image[3], 0)*5.0
	land = ma.masked_equal(land, 0)

	fig, ax = plt.subplots(1,1)
	# ax.imshow(land, alpha = .5)
	# ax.imshow(ivt, alpha=.5)
	# ax.imshow(path, alpha=.5)
	#-------------------------------------------------------------------------#
	# 1.c Plot Wind Barbs
	#-------------------------------------------------------------------------#
	yy      = np.arange(0, v.shape[0], 3)
	xx      = np.arange(0, v.shape[1], 3)
	points  = np.meshgrid(yy, xx)                                 
	xi      = np.arange(0, v.shape[0])
	yi      = np.arange(0, v.shape[1])
	xii, yii = np.meshgrid(yi,xi)
	# place a text box in upper left in axes coords

	levels = np.linspace(100., 2000., 20)#39
	# fig.suptitle()
	ax.imshow(land, alpha = 1.0, cmap='Greys',zorder=1)	

	# cbar = plt.colorbar(ivtclr)
	# cbar.ax.set_ylabel('IVT kgm-1s-1')
	ivtclr = ax.contourf(ivt, cmap='magma_r',alpha = .8, levels=levels, zorder=2)
	ax.quiver(xii[points], yii[points], u[points], v[points], scale = 800, pivot='mid', width =0.001, color='black', alpha=.5, zorder=3) 
	ax.imshow(path, alpha = 1.0, interpolation='None', cmap='gray', zorder=4)



	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5, zorder=5)
	ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
	verticalalignment='top', bbox=props)
	return fig, ax 

# 	plt.show()
# for i in ilist:
# 	fig,ax = PlotAR(i)
# 	plt.show()


record = LocateRecord('1996-05-16 12:00:00', '1.0')
fig,ax = PlotAR(record)

