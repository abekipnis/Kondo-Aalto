from find_atom_positions import CircCorralData
from numpy import array
import pdb
import numpy as np
import scattering_model
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import argparse
import multiprocessing
import warnings
from time import time
import datetime

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

if __name__=="__main__":
	warnings.filterwarnings("ignore")
	parser = argparse.ArgumentParser()
	parser.add_argument("--emin", type=float, default=-0.05)
	parser.add_argument("--emax", type=float, default=0.1)
	parser.add_argument("--n_es", type=int, default=20)
	parser.add_argument("--ngridpoints", type=int, default=20)
	parser.add_argument("--path", type=str, default="test/Createc2_210813.102220.dat")
	args = parser.parse_args()

	if args.emin < scattering_model.E_0.magnitude/scattering_model.electron_charge.magnitude:
		print("minimum energy below surface state onset! ")
		exit(0)

	# "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-12 4p5 nm radius/Createc2_210812.165018.dat"

	now_string = str(datetime.datetime.today())
	datfile = os.path.split(args.path)[-1].strip('.dat')

	fname_head = "%s_%s" %(datfile, now_string)

	c = CircCorralData(args.path, args.path.split("/")[-1])
	c.subtract_plane()
	c.get_region_centroids(diamond_size=5, sigmaclip=2)

	# the box size to fit atom positions
	box_size_nm = 1
	box_size_pix = int(c.nm_to_pix(box_size_nm))
	c.fit_atom_pos_gauss(box_size=box_size_pix)
	c.corral = True
	c.occupied = True

	atoms_n, center_atom_loc = c.remove_central_atom(array(c.centroids))
	atoms_g, center_atom_loc = c.remove_central_atom(c.gauss_fit_locs.T)

	# naive fit from maximum points
	c.r_n, c.c_n = c.nsphere_fit(atoms_n)

	# better fit from gaussian fits to atoms
	c.r_g, c.c_g = c.nsphere_fit(atoms_g)

	# c.compare_fits()
	# atompoints, angle, offseta, offsetb, latt = c.fit_lattice(niter=5)

	erange = np.arange(args.emin, args.emax, (args.emax-args.emin)/args.n_es)
	print("Getting spectra for erange: ", erange)
	nmxyrange = c.pix_to_nm(np.arange(0,c.xPix, c.xPix/args.ngridpoints))
	# this takes too long if using the generated lattice from the fit
	# better to use lattice generated from numpy mesh
	#spectra = scattering_model.gs(atompoints, latt, erange, c.c_g)
	#plt.plot(erange, spectra); plt.imshow()
	#get_spectra(atom_locs (in nm), n_sites (i.e. box size in pixels), r (radius in nm), erange)

	t = time()

	args = [c.pix_to_nm(atoms_g), nmxyrange, erange]

	spectrum = np.array(scattering_model.get_spectra(*args))

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
	# s = scattering_model.get_spectrum_at_middle(c.pix_to_nm(atoms_g), erange)
	# print("it took %1.2lf seconds to get point spectrum" %(time()-t))
	# plt.plot(erange, s);
	# plt.savefig("point_spectrum_test.png")
	# plt.show()

	# l = scattering_model.spectrum_along_line(c.pix_to_nm(atoms_g), erange)
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
