from find_atom_positions import CircCorralData
from numpy import array
import pdb
import numpy as np
import scattering_model
import matplotlib.pyplot as plt
#c = CircCorralData("test/Createc2_210811.092547.dat","test/Createc2_210811.092547.dat")
c = CircCorralData("test/Createc2_210813.102220.dat","Createc2_210813.102220")

print (c.pix_to_nm(c.xPix))

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
atompoints, angle, offseta, offsetb, latt = c.fit_lattice(niter=20)
erange = np.arange(-0.020, 0.020, 0.005)

# this takes way too long if using the generated lattice from the fit
# better to use lattice generated from numpy mesh
#spectra = scattering_model.gs(atompoints, latt, erange, c.c_g)
#plt.plot(erange, spectra); plt.imshow()
pdb.set_trace()
#get_spectra(atom_locs (in nm), n_sites (i.e. box size in pixels), r (radius in nm), erange)
spectrum = scattering_model.get_spectra(atompoints/100, 10, 4, erange) #radius doesn't matter much here, just has to be larger than the radius of the actual atoms (can hard code this, given we have c.r)
spectrum = np.array(spectrum)
for i, e in enumerate(erange):
	plt.close(); plt.imshow(spectrum[i,:,:]); plt.savefig("spectrum_testi_%1.2lf.png" %(e))
plt.close(); plt.plot(erange, spectrum[:,5,5]); plt.savefig("spectrum_test.png")

