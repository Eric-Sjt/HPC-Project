#!/bin/bash
#BSUB -J Explicit
#BSUB -q debug
#BSUB -m 'r13n45'
#BSUB -e %j.err
#BSUB -o %j.out
#BSUB -n 1


module purge
module load intel/2018.4
module load mpi/intel/2018.4

/usr/bin/time -f "run time is: %e" mpirun -np 1 ./explicit.out > $LSB_JOBID.log 2>&1
