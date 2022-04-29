import numpy as np
from .find_atom_positions import CircCorral, CircCorralData
def get_perfect_circle_positions(radius_nm: float,
                                 C: CircCorral,
                                 showfig: bool=True):
    import matplotlib.pyplot as plt

    """

    """

    n_wall_atoms = len(C.centroids) - 1
    center = C.get_central_atom(C.centroids)

    radius = 7

    theta_vals = np.arange(0, 2*np.pi, 2*np.pi/(n_wall_atoms))
    perfect_circle = np.array([[radius*np.cos(theta), radius*np.sin(theta)]
                    for theta in theta_vals])
    perfect_circle = np.concatenate((perfect_circle, [[0,0]]))
    perfect_circle = list(map(C.nm_to_pix, perfect_circle))
    perfect_circle += np.array(center)

    assert(len(perfect_circle)==len(C.centroids))

    if showfig:
        plt.imshow(C.im);
        plt.scatter(*center)
        plt.scatter(*(np.array(perfect_circle).T))
        plt.scatter(*np.array(C.centroids).T)
        plt.show()

    return perfect_circle, center
