from find_atom_positions import CircCorralData
c = CircCorralData("Createc2_210811.092547.dat","Createc2_210811.092547.dat")
c.subtract_plane()
c.get_region_centroids(diamond_size=5, sigmaclip=2)
# the box size to fit atom positions
box_size_nm = 1.5
box_size_pix = int(c.nm_to_pix(box_size_nm))
c.fit_atom_pos_gauss(box_size=box_size_pix)
c.corral = True
c.occupied = True
c.fit_lattice()
