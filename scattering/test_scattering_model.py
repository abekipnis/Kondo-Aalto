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
	ssim.fit_scattering_phase_shift()
	ssim.simulate_and_save_line_spectrum()
	exit(0)

	t = time()

	erange = np.arange(emin, emax, (emax-emin)/n_es)
	print("Getting spectra for erange: ", erange)
	nmxyrange = c.pix_to_nm(np.arange(0,c.xPix, c.xPix/ngridpoints))
	l = sm.spectrum_along_line(c.pix_to_nm(atoms_g), erange)
	np.save("%s_line_spectrum.npy" %(fname_head), l)
	plt.imshow(np.flipud(np.array(l).T), extent=(0, 10, min(erange), max(erange)), aspect="auto")
	plt.show()
	pdb.set_trace()

	args = [c.pix_to_nm(atoms_g), nmxyrange, erange]
	spectrum = np.array(sm.get_spectra(*args))

	for i, e in enumerate(erange):
		plt.close();
		spectrum_data = spectrum[i,:,:]
		np.save("%s_spectrum_%1.2lf.npy" %(fname_head,e), spectrum_data)
		extent = [0,c.pix_to_nm(c.xPix),0,c.pix_to_nm(c.xPix)]
		plt.imshow(np.rot90(spectrum_data.T), extent=extent);
		x,y = np.array(list(map(c.pix_to_nm,c.gauss_fit_locs)))
		plt.scatter(x,y)
		plt.colorbar()
		plt.savefig("%s_spectrum_%1.2lf.png" %(fname_head,e))
	plt.close();

	fig, (ax1) = plt.subplots(1, 1)
	img = ax1.imshow(np.flipud(np.rot90(spectrum[0,:,:])),
					animated=True,
					extent=[0,c.pix_to_nm(c.xPix),0,c.pix_to_nm(c.xPix)])
	ax1.set_xlabel("nm")
	ax1.set_ylabel("nm")

	ax1.set_title("LDOS(V,r)")

	plt.suptitle("Test scattering model", y=0.95)
	def updatefig(i):
		d = spectrum[i,:,:]

		img.set_array(np.flipud(np.rot90(d)))

		return img,
	anim = animation.FuncAnimation(fig, updatefig, frames=len(erange), interval=50, blit=True) #interval in ms
	try:
		try:
			anim.save('%s_cube_movie.mp4' %(fname_head), writer="ffmpeg", fps=28)
		except:
			anim.save('%s_cube_movie.mp4' %(fname_head), fps=28)
	except:
		print("could not save animation :(")
	plt.show()

	plt.plot(erange, spectrum[:,5,5]);
	plt.savefig("spectrum_test.png")


	# pdb.set_trace()
	# s = sm.get_spectrum_at_middle(c.pix_to_nm(atoms_g), erange)
	# print("it took %1.2lf seconds to get point spectrum" %(time()-t))
	# plt.plot(erange, s);
	# plt.savefig("point_spectrum_test.png")
	# plt.show()

	# l_data = np.rot90(np.array(l))
	# np.save("line_spectrum_test_%s.npy" %(args.path.split("/")[-1]), l_data)
	# plt.imshow(l_data);
	# plt.savefig("line_spectrum_test_%s.png" %(args.path.split("/")[-1]))


	# set this up so it runs command line in chunks ??
	# to avoid the problem where it goes over the 15 minute time limit
	# or just run this using sbatch instead of sh
	# shcommand = """
	# #!/bin/bash
	# #SBATCH --time=04:00:00      # 4 hours
	# #SBATCH --mem=1000M   # 1G of memory
	# #SBATCH --cpus-per-task=4
	# #SBATCH --mem=8G
	# #SBATCH -o model_LDOS.out
	#
	# export OMP_PROC_BIND=true
	# echo 'Running on :'$HOSTNAME
	# srun python3 test_scattering_model.py --emin=-0.066 --emax=0.3 --n_es=10 --ngridpoints=100 --path="test/Createc2_210812.170231.dat"
	# """

	# naive fit from maximum points
	# c.r_n, c.c_n = c.nsphere_fit(atoms_n)

	# better fit from gaussian fits to atoms
	# c.r_g, c.c_g = c.nsphere_fit(atoms_g)

	# c.compare_fits()
	# atompoints, angle, offseta, offsetb, latt = c.fit_lattice(niter=5)
