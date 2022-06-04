import numpy as np
from os import path
from .find_atom_positions import CircCorral
try:
    from AMRL import Createc_Controller
except Exception as e:
    print(e)
    print("Createc_Controller is None")
    print(path.basename(__file__) + " will not work")
    Createc_Controller = None
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
import random

def get_latman_dest(circle, C, safety_distance_nm):
    safe_dist_px = C.nm_to_pix(safe_dist_nm)
    safe_radius = circle[0] - safe_dist_px

    # randomly sample in a circle of this radius and convert to cartesian coords
    r = safe_radius * np.sqrt(random.random())
    theta = random.random()*2*np.pi
    x = circle[1][0] + r * np.cos(theta)
    y = circle[1][0] + r * np.sin(theta)

    # destination for the atom
    dest = [x, y]
    return dest

def get_atom_that_moved_the_most(old_centroids, CC):
    import matplotlib.pyplot as plt

    Cim = np.array(CC.stm.scandata(1,4))
    zconst = float(CC.stm.getparam('ZPiezoConst'))
    nmx = nmy = CC.get_len_nm()
    C = CircCorral(Cim, zconst, nmx, nmy)
    C.subtract_plane()
    C.get_region_centroids(percentile=98, edge_cutoff=0.1, show=False)

    sorted_new = np.array(sorted(C.centroids, key=lambda x:x[0]**2+x[1]**2))
    sorted_old = np.array(sorted(old_centroids, key=lambda x:x[0]**2+x[1]**2))

    cost_matrix = cdist(sorted_old, sorted_new)
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    sorted_new = sorted_new[col_ind,:]
    sorted_old = sorted_old[row_ind,:]

    norm = np.linalg.norm(sorted_old-sorted_new, axis=-1)
    largest_move_index = np.argmax(norm)

    if C.pix_to_nm(np.max(norm)) < 0.2:
        print("atom did not move")
        return False, False

    central_atom_loc = sorted_new[largest_move_index]

    original_circle = C.nsphere_fit(np.delete(sorted_new,largest_move_index, axis=0))

    plt.imshow(C.im)
    plt.scatter(*central_atom_loc)
    plt.show()
    return central_atom_loc, original_circle
