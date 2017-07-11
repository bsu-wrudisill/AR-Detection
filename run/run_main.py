import atmDetect

#if __name__ == '__main__':
#    main()

#path = '/home/wrudisill/scratch/Find_ARs/ivt_files/pgbhnl.gdas.20000201-20000205.nc'

path = 'foo.nc'

for i in range(0,1):
        atmDetect.FindAR(path, i)
        print 'done with %s' %i


