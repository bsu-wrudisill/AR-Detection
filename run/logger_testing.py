import logging


#log_name = datetime.now().strftime('%Y-%m-%d_%H:%M:%S') + '_AR_Detect.log'

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger('foo')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
logger.setLevel(logging.INFO)
fh = logging.FileHandler("testing.log")
fh.setFormatter(formatter)
logger.addHandler(fh)
logger.info('Start the log')



for i in range(10):
    try:
        i/0.0
        logger.info('completed_%s', i)
    except Exception as E:
        logger.info(E)rm

logger.info('Done.')
