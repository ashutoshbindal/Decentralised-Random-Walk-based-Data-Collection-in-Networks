#!/bin/sh
#PBS -N simulations
### Set the project name, your department dc by default
#PBS -P cse
### Request email when job begins and ends
#PBS -m bea
### Specify email address to use for notification.
#PBS -M mcs162663@iitd.ac.in
####
#PBS -l select=1:ncpus=7
### Specify "wallclock time" required for this job, hhh:mm:ss
#PBS -l walltime=150:00:00

#PBS -l software=
# After job starts, must goto working directory. 
# $PBS_O_WORKDIR is the directory from where the job is fired. 
echo "==============================="
echo $PBS_JOBID
cat $PBS_NODEFILE
echo "==============================="
cd $PBS_O_WORKDIR
#job 
time ./run.sh
#NOTE
# The job line is an example : users need to change it to suit their applications
# The PBS select statement picks n nodes each having m free processors
# OpenMPI needs more options such as $PBS_NODEFILE
