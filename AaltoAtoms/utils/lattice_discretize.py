import numpy as np
import pdb
from sklearn.neighbors import NearestNeighbors

def lattice_discretize(points: list,
                       npix: int,
                       len_nm: float,
                       add_hcp_sites: bool=False,
                       rotation:float=0) -> list:
    """
    snap points to adsorption sites on Ag(111) surface
    """
    trigrid = np.array([[],[]])

    a = 2.887e-10; # m (value for Ag)
    h_dist = 1e9*a*npix/len_nm; # Horizontal distance, pixel units
    v_dist = np.sqrt(h_dist**2-(h_dist/2)**2); # Vertical distance

    x_lim = 2*npix;
    y_lim = 2*npix;
    y_curr = 0
    xx = 0

    #displacement
    disp = 0
    while y_curr < y_lim:
        if disp == 0:
            xx = np.arange(0, x_lim, h_dist)
            yy = np.ones(len(xx))*y_curr
            disp = 1
        else:
            xx = [h_dist/2, x_lim, h_dist]
            yy = np.ones(len(xx))*y_curr
            disp = 0

        trigrid = np.concatenate([trigrid, [xx, yy]], axis=-1)
        y_curr = y_curr + v_dist

    if add_hcp_sites:
        pdb.set_trace()
        hcpgrid = (np.array(trigrid).T + np.array([0.5*h_dist, v_dist/3]).T).T
        trigrid = np.concatenate([trigrid, hcpgrid], axis=-1)

    trigrid = np.array(trigrid).T

    # find the nearest point in the lattice
    neigh = NearestNeighbors(n_neighbors=6)
    neigh.fit(trigrid)

    # get lattice points closest to atoms
    distances, idxs = neigh.kneighbors(points, 1)

    perfect_circle = trigrid[idxs]
    return perfect_circle.reshape(len(perfect_circle),2)
