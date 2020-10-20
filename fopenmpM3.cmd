#!/bin/bash
#PBS -o logfile.log
#PBS -e errorfile.err
#PBS -l walltime=16:00:00
#PBS -l select=1:ncpus=20
tpdir=`echo $PBS_JOBID | cut -f 1 -d .`
tempdir=$HOME/scratch/job$tpdir
mkdir -p $tempdir
cd $tempdir
cp -R $PBS_O_WORKDIR/* .
export OMP_NUM_THREADS=20
make
make clean
./run > output.txt
rm run
zip -rm result.zip plot* residual.dat restart.dat
mv ../job$tpdir/result.zip  $PBS_O_WORKDIR/.
mv ../job$tpdir $PBS_O_WORKDIR/.
