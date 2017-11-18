#/usr/bin/python
import numpy as np
import glob

# Script loops thru 85th p files and calculates the mean of a stack 


def zero_pad(mo):
    if mo < 10:
        mostr = '0'+str(int(mo))
    else:
        mostr = str(int(mo))
    return mostr


for j in range(1,13):
    month = zero_pad(j)
    files = glob.glob('./data/monthly_stats/IVT_'+month+'*.npy')
    flen = len(files)
    master = np.ones((flen, 361,720))
    k = 0 
    for i in files:
        master[k, :, :] = np.load(i)
        k = k+1
        mn = np.mean(master, axis=0)
        np.save('./data/Month_'+month+'.npy', mn)

