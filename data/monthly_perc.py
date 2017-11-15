import numpy as np
from netCDF4 import Dataset
import glob
import sys
from multiprocessing import Pool
import gc


gc.enable()


def percentile(n):
    '''
    Helper function to calc a given percentile 
    '''
    def percentile_(x):
        return np.percentile(x, n, axis=0)
    percentile_.__name__ = 'percentile_%s' % n
    return percentile_


def calc_ivt(filename, time):
# calculate IVT and other things...
    ds                      = Dataset(filename, format='NETCDF4_CLASSIC')        
#lons                    = ds.variables['lon_0'][:]

#	lats                    = ds.variables['lat_0'][:]
#	_lons_mesh, lats_mesh   = np.meshgrid(lons, lats)
#	lons_mesh               = np.where (_lons_mesh > 180.0, _lons_mesh - 360, _lons_mesh) 
#	hr_time                 = format_date(ds.variables['initial_time0_hours'][time])
#	hr_time_str             = hr_time.strftime("%Y-%m-%d_%H")
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
    ds.close()
    gc.collect()

#U and V winds respectivl
#	vgrd =  ds.variables['VGRD_P0_L100_GLL0'][time,:,:,:]
#	ugrd =  ds.variables['UGRD_P0_L100_GLL0'][time,:,:,:]
#	spfh =  ds.variables['SPFH_P0_L100_GLL0'][time,:,:,:]			 
    return ivt_integral


month = sys.argv[1]
#month ='01'

files =glob.glob('/home/wrudisill/scratch/AR-Detection/data/ncfiles/pgbhnl.gdas.????'+month+'??-*.nc')

def calc_ivt_wrapper(filename):
    print 'working on .....' + filename
    ds    = Dataset(filename, format='NETCDF4_CLASSIC')        
    tlen  = len(ds.variables['initial_time0_hours'])
    store_1file  = np.zeros((tlen-1, 361, 720))
    ds.close()
    for i in range(tlen-1):
        store_1file[i,:,:] = calc_ivt(filename, i)

    return store_1file



prcntl85 = percentile(85)

flen = len(files)
chunk = map(int, np.linspace(0,flen,500)

p = Pool(12)
arlist = map(calc_ivt_wrapper, files)

# np array to contain evert timestep for every file in the month;;; might be very big...
master = np.concatenate(arlist, axis=0)
np.save('./monthly_stats/master_'+month+'.npy', master)

gc.collect()


#out = np.zeros((361,720))

#out = prcntl85(master)

#for i in range(master.shape[0]):
#    for j in range(master.shape[1]):
#        out[i,j] = prcntl85(master[i,j])


# calculate percentile 
# out = prcntl85(master)
#np.save('./monthly_stats/IVT_'+month+'_85thp.npy', out)








