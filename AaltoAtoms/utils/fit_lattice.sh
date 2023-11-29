#!/bin/bash
#SBATCH --time=04:00:00      # 4 hours
#SBATCH --mem=1000M   # 1G of memory
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH -o model_LDOS.out

export OMP_PROC_BIND=true
echo 'Running on :'$HOSTNAME
srun python3 fit_lattice.py
