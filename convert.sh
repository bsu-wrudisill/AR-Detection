#/bin/sh

for i in $(cat difflist);
do
    ncl_convert2nc /home/wrudisill/scratch/Find_ARs/CFS_Reanalysis_Partial_RDA/pgbhnl.gdas.$i.grb2 -o input_files/.

done
