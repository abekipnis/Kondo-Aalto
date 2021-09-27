#!/bin/bash
#SBATCH --time=04:00:00      # 4 hours
#SBATCH --mem=1000M   # 1G of memory

srun python3 fit_lattice.py
