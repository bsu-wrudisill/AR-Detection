import atmDetect
from multiprocessing import Pool
from glob import glob
import numpy as np

#enable automatic garbage collection
np.seterr(all='print')

#Multiprocessing Pool
p = Pool(28)

# File List
file_list = glob('/home/wrudisill/scratch/Find_ARs/data/new_ivt_files/IVT*')


# Map Function 
p.map(atmDetect.FindAR_Wrapper, file_list)
p.close()
p.join()






