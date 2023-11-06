#!/bin/bash
#SBATCH -N 8
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -J ukmerc-celegans-8nodes
#SBATCH --mail-user=yl3722@cornell.edu
#SBATCH --mail-type=ALL
#SBATCH -t 00:05:00

# GASNET_BACKTRACE=1 upcxx-run -shared-heap 16G -n 2 -N 1  ./ukmerc /pscratch/sd/y/yfli03/ELBA_dataset/acinetobacter_baumannii/reads.fa
# upcxx-run -shared-heap 16G -n 16 -N 2  ./ukmerc /pscratch/sd/y/yfli03/ELBA_dataset/acinetobacter_baumannii/reads.fa

#OMP settings
export OMP_NUM_THREADS=16
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

module load upcxx contrib

upcxx-srun -shared-heap 32G -n 64 -N 8 -c 32 --cpu_bind=cores -- ./ukmerc /pscratch/sd/y/yfli03/ELBA_dataset/celeg40x/reads.fa | tee "s4r1b100000t.log"

# upcxx-run -shared-heap 16G -n 8 -N 1 ./ukmerc /pscratch/sd/y/yfli03/ELBA_dataset/acinetobacter_baumannii/reads.fa | tee acine.log
