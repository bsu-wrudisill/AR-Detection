#!/bin/bash

# ---------------------- Description ----------------------------------- #

# 1) Runs AR Detection code. Supply 'submit_python.sh' with year argumnet.
# 2) Moves Files to /mnt/ while the jobs are running 

# ---------------------- Description ----------------------------------- #



#sbatch submit_python.sh 1979 | cut -d ' ' -f 4 >> .catch.jobs
#sbatch submit_python.sh 1980 | cut -d ' ' -f 4 >> .catch.jobs
#sbatch submit_python.sh 1981 | cut -d ' ' -f 4 >> .catch.jobs
#sbatch submit_python.sh 1983 | cut -d ' ' -f 4 >> .catch.jobs
#sbatch submit_python.sh 1984 | cut -d ' ' -f 4 >> .catch.jobs
#sbatch submit_python.sh 1985 | cut -d ' ' -f 4 >> .catch.jobs


sbatch submit_python.sh 1986 | cut -d ' ' -f 4 >> .catch.jobs
sbatch submit_python.sh 1987 | cut -d ' ' -f 4 >> .catch.jobs
sbatch submit_python.sh 1988 | cut -d ' ' -f 4 >> .catch.jobs
sbatch submit_python.sh 1989 | cut -d ' ' -f 4 >> .catch.jobs
sbatch submit_python.sh 1990 | cut -d ' ' -f 4 >> .catch.jobs
sbatch submit_python.sh 1991 | cut -d ' ' -f 4 >> .catch.jobs
sbatch submit_python.sh 1992 | cut -d ' ' -f 4 >> .catch.jobs
sbatch submit_python.sh 1993 | cut -d ' ' -f 4 >> .catch.jobs
sbatch submit_python.sh 1994 | cut -d ' ' -f 4 >> .catch.jobs
sbatch submit_python.sh 1995 | cut -d ' ' -f 4 >> .catch.jobs
sbatch submit_python.sh 1996 | cut -d ' ' -f 4 >> .catch.jobs


# --- ----- We Need to Move files that get created to /mnt/... but copute nodes to not have /mnt access -------##
# Description:
# 1) Counts the numer of jobs submitted and writes a file of the job ID
# 2) Queries slurm for the list of submitted and counts the numer of matches with the jobs submitted
# 3) Moves files to /mnt/
# 4) When the match count == 0 (none of the submitted jobs are active), the while loop ends, and we stop trying to move files 


# numer of jobs submitted:
COUNTER=$( cat .catch.jobs | wc -l )

# grep jobs from QUEUE
while [ $COUNTER -gt 0 ]; do 
    TMP_COUNTER=0
    # RECALL: > writes over what is previously in the file
    squeue -u wrudisill > .catch.info.tmp                    # get job id's from queue
    job_id=$(grep "$job_id" .catch.info.tmp |  sed "s/^ *//" | cut -d ' ' -f1)  #jobs in queue
    for line in $job_id; do 
	match=$( cat .catch.jobs | grep $line )
	if ! [ -z $match ]; then
	    TMP_COUNTER=$(( $TMP_COUNTER + 1 ))
	fi
    done  #end for
    COUNTER=$TMP_COUNTER    

    #--------- MOVE FILES --------------#
    destination='/mnt/selway/data/data_03/will_r/AR_Detection_Outfiles_0/.'
    opath='/home/wrudisill/scratch/AR-Detection/data/outfiles'
    files=$(ls -U $opath/*.npy 2> /dev/null | wc -l)
    if [ "$files" != "0" ]
    then
	mv -v $opath/*.npy $destination >> transfer_log
    else
	echo 'no files' >> transfer_log
    fi
    #-----------------------------------#

    sleep 10  #wait 10 seconds 
    done


#clean up
rm .catch.jobs
rm .catch.info.tmp
















