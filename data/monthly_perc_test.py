import numpy as np 
from multiprocessing import Pool
import sys,os
import glob
import itertools
from netCDF4 import Dataset




#month = sys.argv[1]
month ='01'

files =glob.glob('/home/wrudisill/scratch/AR-Detection/data/ncfiles/pgbhnl.gdas.????'+month+'??-*.nc')
N = len(files)



def calc_ivt_wrapper(ix, jx):

# ---------- Description ------------------
# Retuns the calc IVT function for a grid point, mapping thru the 
# time dimension
# -----------------------------------------

    def _calc_ivt(filename, time, ix, jx):
# calculate IVT and other things...
        ds                      = Dataset(filename, format='NETCDF4_CLASSIC')        
        ivt =  np.zeros(20)
        g = 9.81 
        dp = 25.00 #pa

# loop through pressure levels, calculate IVT 
        for i in range(20):
            spfh = ds.variables['SPFH_P0_L100_GLL0'][time,i,ix,jx]
            vgrd = ds.variables['VGRD_P0_L100_GLL0'][time,i,ix,jx]
            ugrd = ds.variables['UGRD_P0_L100_GLL0'][time,i,ix,jx]
            phi =np.arctan2(vgrd,ugrd)
            tV = (ugrd**2)*(vgrd**2)**.5
            ivt[i] = spfh*tV/g*dp

# sum through pressure levels 
        ivt_integral = np.sum(ivt)  # returns 1 pt.
        ds.close()
        return ivt_integral	


    def calc_ivt(filename, ix, jx):
# finds time dimesnsion and loops thru time

        print 'working on .....' + filename
        ds    = Dataset(filename, format='NETCDF4_CLASSIC')        
        tlen  = len(ds.variables['initial_time0_hours'])
        store_1px  = np.zeros(tlen-1)
        ds.close()

        for i in range(tlen-1):
            value = _calc_ivt(filename, i, ix, jx)
            print value.shape
            store_1px[i] = value
# print _calc_ivt(filename, i, ix, jx)
        return store_1px

    def mapper_calc_ivt(filename):
        return calc_ivt(filename, ix, jx)

# returns function mapper_calc_ivt, where ix and jx are set
# the time dimesion is figured out internally
# returns a length t array 
    return mapper_calc_ivt


grid = np.zeros((361,720))

# ugh... 
for ix in range(361):
    for jx in range(720):
        calc_ivt = calc_ivt_wrapper(ix,jx)  #define calc_ivt function
#p   = Pool(8)                 
        out = map(calc_ivt, files[0:2])       #loop thru files; this makes an n_files*timesteps array
        zdim        = np.concatenate(out, axis=0)
        grid[ix,jx] = np.percentile(zdim, 85)


np.save('test85.npy', grid)






# out_arr[ind[j]:ind[j+1]] = np.percentile(chunk, 85, axis=0)
# del chunk
