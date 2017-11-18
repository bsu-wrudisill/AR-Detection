import numpy as np
import dask.array as da
import dask 

foo = np.load('master_01_01.npy')
y = da.from_array(foo, chunks=100)
std = y.std(axis=0).compute()





print std

