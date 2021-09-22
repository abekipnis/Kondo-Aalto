from find_atom_positions import CircCorralData
from numpy import array

c = CircCorralData("test/Createc2_210811.092547.dat","test/Createc2_210811.092547.dat")
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

c.fit_lattice()
