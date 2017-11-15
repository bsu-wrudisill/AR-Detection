from reorg import *
from multiprocessing import Pool
from glob import glob
import numpy as np
import sys
#enable automatic garbage collection
np.seterr(all='print')



yr = sys.argv[1]


#yr = '1979'
string = '/home/wrudisill/scratch/AR-Detection/data/ncfiles/pgbhnl.gdas.YEAR*.nc'

flist = glob(string.replace('YEAR',yr))


# Setup logger 
fh_string = 'AR_detect_YEAR.log'
fh = logging.FileHandler(fh_string.replace('YEAR',yr))
logger = logging.getLogger('log')
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
logging.basicConfig(format=FORMAT)
logger.setLevel(level=logging.DEBUG)
logger.addHandler(fh)


logger.info('Start the log')

for f in flist:
	logger.info('will process '+f)

# Wrapper function for processing 
def wrapper(filename):
        # Try, except block
	try:
                time_length = GrepTimeDim(filename)
                print time_length
                for time in range(time_length):
                        ivt_timeslice = calc_ivt(filename, time)
                        blob_tester(ivt_timeslice)
                        logger.debug('completed' + filename + '__' + str(time))

        # Princt traceback information to logger if there is an error
	except Exception as e:
#		logger.findCaller()
		logger.error(e, exc_info=True)
		# logger.debug('something didnt work \n %s', e)

print flist 
p = Pool(8)
p.map(wrapper, flist)
p.close()
p.join

# Finish logger
logger.info('Done')


# bar = np.zeros((361,720))
# for i in out:
# 	bar+=i[0].path
# land       = np.load('../data/land.npy')

