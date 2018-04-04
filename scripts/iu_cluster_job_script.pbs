#! /bin/bash
#PBS -l	nodes=4:ppn=16,walltime=48:00:00
#PBS -q	gpu
#PBS -m	ae
#PBS -o	out.log
#PBS -e	err.log

# below are the modules you will need to compile the code on bigred2 (see README)
# uncomment the following 3 lines to load the modules at the time of job launch
module swap PrgEnv-cray PrgEnv-gnu
module load boost/1.65.0
module load gsl

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=16
# -d refers to number of cores. this should match ppn in Line 2.
time aprun -n 4 -d 16 ./np_electrostatics_lab -a 2.6775 -b 14.28 -e 2 -E 78.5 -V -60 -v 1 -g 1082 -m 6 -t 0.001 -s 10000 -p 100 -f 10 -M 6 -T 0.001 -k 0.0025 -q 0.001 -L 5 -l 5 -S 10000000 -P 100000 -F 100 -X 10000 -U 1000 -Y 500000 -W 1000000 -B 0.025