import numpy as np
from netCDF4 import Dataset
import glob
import sys
from multiprocessing import Pool
import gc
import time 

gc.enable()

#this is super ghetto...
flag = int(time.time())


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
    del spfh
    del vgrd
    del phi
    del tV
    del ivt
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
    k = 10 
    #store_1file  = np.zeros((tlen-1, 361, 720))
    store_1file  = np.zeros((k, 361, 720))

    ds.close()
    #for i in range(tlen-1):
    #    store_1file[i,:,:] = calc_ivt(filename, i)
    
    rlist = list(np.random.randint(0,tlen-1, k))
    for i in range(k-1):
        store_1file[i,:,:] = calc_ivt(filename, rlist[i])
    

    return store_1file



prcntl85 = percentile(85)

flen = len(files)


# create a random list of files to compute 
iterations     = 50
number_of_files = 100
#index        = list(np.random.randint(0,flen-1, number_of_files))
#random_files = [files[x] for x in index]
#p = Pool(4)
#arlist = p.map(calc_ivt_wrapper,random_files)
# np array to contain evert timestep for every file in the month;;; might be very big...
#master = np.concatenate(arlist, axis=0)
#out = np.percentile(master, 85, axis=0)
#np.save('./monthly_stats/IVT_'+month+'_'+str(flag)+'_85thp.npy', out)


for k in range(iterations):
   
    index        = list(np.random.randint(0,flen-1, number_of_files))
    random_files = [files[x] for x in index]
   
    p = Pool(4)

    arlist = p.map(calc_ivt_wrapper,random_files)
    p.close()
    p.join
    # np array to contain evert timestep for every file in the month;;; might be very big...
    master = np.concatenate(arlist, axis=0)
    
    #np.save('./monthly_stats/master_'+month+'_01.npy', master)
    out = np.percentile(master, 85, axis=0)
    
    # save percentile
    np.save('./monthly_stats/IVT_'+month+'_'+str(k)+'_85thp.npy', out)

    del master
    del out
    del arlist 
    gc.collect()    






