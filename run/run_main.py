import atmDetect
from multiprocessing import Pool
from glob import glob
import gc


#enable automatic garbage collection
gc.enable()

#Multiprocessing Pool
p = Pool(14)

# File List
file_list = glob('/home/wrudisill/scratch/Find_ARs/data/new_ivt_files/IVT_2010*')

# Map Function 
p.map(atmDetect.FindAR_Wrapper, file_list)

p.close()
p.join()


print 'Done'






