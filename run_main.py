import atmDetect


#if __name__ == '__main__':
#    main()

path = '/home/wrudisill/scratch/Find_ARs/ivt_files/pgbhnl.gdas.20000201-20000205.nc'

#pgbhnl.gdas.20100601-20100605.nc'

ivt_min = 250                     # Minimum IVT value in kg/ms to be retained
size_mask  = 1000                  # Min Grid cell size of object
cell_to_km = 50                   # km
Results_dictionary = {}           # output dictionary


for i in range(0,20):
        atmDetect.FindAR(path, i, ivt_min = 250, size_mask = 1000, cell_to_km = 50)

#        print 'done with %s' %i


