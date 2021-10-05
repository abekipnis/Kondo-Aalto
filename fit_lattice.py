from find_atom_positions import CircCorralData
from numpy import array
import pdb
import numpy as np
import scattering_model
import matplotlib.pyplot as plt
import argparse
import multiprocessing
from time import time
if __name__=="__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--emin", type=float, default=-0.1)
	parser.add_argument("--emax", type=float, default=0.1)
	parser.add_argument("--n_es", type=int, default=2)
	parser.add_argument("--ngridpoints", type=int, default=10)
	parser.add_argument("--path", type=str, default="test/Createc2_210813.102220.dat")
	args = parser.parse_args()

	# "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-12 4p5 nm radius/Createc2_210812.165018.dat"

	""" for example, run from the command line like:
	python3 fit_lattice.py --emin=-0.02 --emax=0.02 --n_es=5 --ngridpoints=100
	"""

	#c = CircCorralData("test/Createc2_210811.092547.dat","test/Createc2_210811.092547.dat")
	print(args.path)
	c = CircCorralData(args.path, "Createc2_210812.165018.dat")

	c.subtract_plane()
	c.get_region_centroids(diamond_size=5, sigmaclip=2)

	# the box size to fit atom positions
	box_size_nm = 1
	box_size_pix = int(c.nm_to_pix(box_size_nm))
	c.fit_atom_pos_gauss(box_size=box_size_pix)
	c.corral = True
	c.occupied = True

	atoms_n = c.remove_central_atom(array(c.centroids))
	atoms_g = c.remove_central_atom(c.gauss_fit_locs.T)

	# naive fit from maximum points
	c.r_n, c.c_n = c.nsphere_fit(atoms_n)

	# better fit from gaussian fits to atoms
	c.r_g, c.c_g = c.nsphere_fit(atoms_g)

	c.compare_fits()
	atompoints, angle, offseta, offsetb, latt = c.fit_lattice(niter=5)
	atompoints = atoms_g
	erange = np.arange(args.emin, args.emax, (args.emax-args.emin)/args.n_es)

	nmxyrange = c.pix_to_nm(np.arange(0,c.xPix, c.xPix/args.ngridpoints))
	# this takes too long if using the generated lattice from the fit
	# better to use lattice generated from numpy mesh
	#spectra = scattering_model.gs(atompoints, latt, erange, c.c_g)
	#plt.plot(erange, spectra); plt.imshow()
	#get_spectra(atom_locs (in nm), n_sites (i.e. box size in pixels), r (radius in nm), erange)
	plt.close()

	t = time()
	s = scattering_model.get_spectrum_at_middle(c.pix_to_nm(atompoints), erange)
	print("it took %1.2lf seconds to get point spectrum" %(time()-t))
	plt.plot(erange, s);
	plt.show()

	pdb.set_trace()

	l = scattering_model.spectrum_along_line(c.pix_to_nm(atompoints), erange)

	spectrum = scattering_model.get_spectra(c.pix_to_nm(atompoints), nmxyrange, erange)
	spectrum = np.array(spectrum)
	plt.plot(erange, s); plt.show()
	for i, e in enumerate(erange):
		plt.close();
		plt.imshow(spectrum[i,:,:], extent=[0,c.pix_to_nm(c.xPix),0,c.pix_to_nm(c.xPix)]);
		plt.savefig("spectrum_testi_%1.2lf.png" %(e))
	plt.close();

	plt.plot(erange, spectrum[:,5,5]);
	plt.savefig("spectrum_test.png")
