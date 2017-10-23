from reorg import *
from multiprocessing import Pool
from glob import glob
import numpy as np

#enable automatic garbage collection
np.seterr(all='print')

#Multiprocessing Pool
# p = Pool(28)

# File List
# file_list = glob('/home/wrudisill/scratch/Find_ARs/data/new_ivt_files/IVT*')


# # Map Function 
# p.map(atmDetect.FindAR_Wrapper, file_list)
# p.close()
# p.join()

fh = logging.FileHandler("testing.log")
logger = logging.getLogger(__name__)
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
logging.basicConfig(format=FORMAT)
logger.setLevel(level=logging.DEBUG)
logger.addHandler(fh)



foo = []


logger.info('Start the log')

def wrapper(time):
	try:
		ivt_timeslice = calc_ivt('../data/pgbhnl.gdas.19960516-19960520.nc', time)
		bar = blob_tester(ivt_timeslice)
		logger.debug('completed time %s', time)
		return bar

	except Exception as e:
#		logger.findCaller()
		logger.error(e, exc_info=True)
		# logger.debug('something didnt work \n %s', e)


p = Pool(1)
a = range(4)
out = p.map(wrapper, a)
p.close()
p.join


logger.info('Done')


# bar = np.zeros((361,720))
# for i in out:
# 	bar+=i[0].path


# land       = np.load('../data/land.npy')

