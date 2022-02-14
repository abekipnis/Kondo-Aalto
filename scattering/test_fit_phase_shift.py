import pdb
import sys

sys.path.append('/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/Kondo data analysis')
sys.path.append("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis")
sys.path.append("/home/kipnisa1/python/Small-Kondo-Corrals/Kondo data analysis")
sys.path.append("/home/kipnisa1/python/Small-Kondo-Corrals/")
import scipy
from numpy import array
import pdb
import numpy as np
import scattering_simulation as ss
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import argparse
import warnings
from time import time
import socket

""" for example, run from the command line like:
python3 test_scattering_model.py --emin=-0.02 --emax=0.02 --n_es=5 --ngridpoints=100
"""

# TODO: save the files without space
# TODO: scatter atom positions over LDOS
# TODO: calculate spectrum along with LDOS maps
# TODO: add in Kondo phase shift for occupied corrals
# TODO: optimize speed
# TODO: longer runs
# TODO: figure out infinity / black dot scattering, etc. derive math better
# TODO: compare point spectra with unoccupied corrals
# TODO: compare CH maps with unoccupied corrals
# TODO: compare with line spectrum from empty corrals

def get_args(parser):
	parser.add_argument("--emin", type=float, default=-0.05)
	parser.add_argument("--emax", type=float, default=0.1)
	parser.add_argument("--n_es", type=int, default=20)
	parser.add_argument("--ngridpoints", type=int, default=20)

	# "/m/phys/project/asp/labdata/Createc_new/STMDATA/Ag/Small Kondo corrals/"
	d_dat = "../test/Createc2_210816.170832.dat"
	d_linespec_dir = "Ag 2021-08-16 2p5 nm radius empty/3p8 nm pm100mV line/"

	parser.add_argument("--path", type=str, default=d_dat)
	parser.add_argument("--linespec_dir", type=str, default=d_linespec_dir)

	args = parser.parse_args()
	return args

if __name__=="__main__":
	host = socket.gethostname()
	print("Running on host: ", host)

	warnings.filterwarnings("ignore")
	args = get_args(argparse.ArgumentParser())
	path = args.path								# the path to the .dat file
	linespec_dir = args.linespec_dir				# path to folder with spectra
	emin = args.emin								# min of energy range for model
	emax = args.emax								# max of energy range for model
	n_es = args.n_es								# number of points in energy range for model
	ngridpoints = args.ngridpoints					# number of grid points (in 1 dimension)

	localdir = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/"
	tritondir = "/scratch/work/kipnisa1/Small Kondo corrals/"

	on_triton = "triton" in host
	dir = tritondir if on_triton else localdir
	linespec_dir = dir + linespec_dir

	if on_triton:
		pdb.set_trace = lambda: 1

	ssim = ss.ScatteringSimulator(linespec_dir, path)
	ssim.fit_scattering_phase_shift(n_bias=20, n_spectra=20)
	exit(0)
