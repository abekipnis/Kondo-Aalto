#!/bin/bash
#SBATCH --time=04:00:00      # 4 hours
#SBATCH --mem=1000M   # 1G of memory
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -o model_LDOS.out

export OMP_PROC_BIND=true
echo 'Running on :'$HOSTNAME
#srun python3 fit_lattice.py --emin=-0.02 --emax=0.02 --n_es=5 --ngridpoints=100
srun python3 test_scattering_model.py --emin=-0.067 --emax=0.3 --n_es=3 --ngridpoints=30 --path="test/Createc2_210812.170231.dat"
