from test_scattering_model import replicate_spectra, get_args
import socket, argparse
import scattering_model as sm
import warnings
import pdb
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
	if emin < sm.E_0.magnitude/sm.electron_charge.magnitude:
		print("minimum energy below surface state onset!")
		exit(0)
	localdir = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/"
	tritondir = "/scratch/work/kipnisa1/Small Kondo corrals/"

	on_triton = "triton" in host
	dir = tritondir if on_triton else localdir
	linespec_dir = dir + linespec_dir

	# deactivate pdb.set_trace function debugging calls if on triton
	if on_triton:
		pdb.set_trace = lambda: 1

	replicate_spectra(linespec_dir, path)
