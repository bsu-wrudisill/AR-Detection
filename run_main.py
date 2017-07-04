import atmDetect

#if __name__ == '__main__':
#    main()

path = '/home/wrudisill/scratch/Find_ARs/ivt_files/pgbhnl.gdas.20000201-20000205.nc'


for i in range(0,20):
        atmDetect.FindAR(path, i)
        print 'done with %s' %i


