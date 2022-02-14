#!/bin/bash
#SBATCH --time=08:00:00      # 4 hours
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH -o model_LDOS.out

export OMP_PROC_BIND=true
echo 'Running on :'$HOSTNAME
#srun python3 fit_lattice.py --emin=-0.02 --emax=0.02 --n_es=5 --ngridpoints=100


# run this with sbatch test_scattering_model.sh

srun python3 test_scattering_model.py
# srun python3 test_scattering_model.py --emin=-0.066 --emax=0.3 --n_es=10 --ngridpoints=100 --path="test/Createc2_210812.170231.dat"
