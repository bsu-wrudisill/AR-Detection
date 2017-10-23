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
cur = conn.cursor()


def PullCol(col, db):
	cur.execute("SELECT "+col+" FROM "+db)# WHERE landfalling='True'")  
	field = cur.fetchall()
	field = map( lambda X: reduce(list, X), field)
	return field






hr_time_str   = PullCol('hr_time_str', 'Rivers')
mean_IVT      = PullCol('mean_IVT', 'Rivers')
mean_wind_dir = PullCol('wind_dir_mean', 'Rivers')


dic = {'Date':hr_time_str, 'IVT':mean_IVT, 'mean_wind_dir':mean_wind_dir}
df = pd.DataFrame(data=dic)
df.set_index(pd.to_datetime(df['Date'], format='%Y-%m-%d_%H'), inplace=True)
del df['Date']
df['Year'] = df.index.year
df['Month'] = df.index.month
df['Count'] = 1




#image = np.load('../data/outfiles/1996-05-16_18_OBJECT_1.0.npy')

ilist = glob.glob('../data/outfiles/*npy')

def PlotAR(fname):
	image = np.load(fname)
	land  = np.load('../data/land.npy')
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
	yy      = np.arange(0, v.shape[0], 5)
	xx      = np.arange(0, v.shape[1], 5)
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

	textstr = str(i)
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5, zorder=5)
	ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
	verticalalignment='top', bbox=props)
	return fig, ax 

# 	plt.show()

for i in ilist:
	fig,ax = PlotAR(i)
	plt.show()

