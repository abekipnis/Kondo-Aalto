# Kabsch algorithm for the best rotation between two sets to minimize distance
from scipy.spatial.transform import Rotation
import numpy as np
from os import path
from .find_atom_positions import CircCorral, CircCorralData, clockwiseangle
try:
    from AMRL import Createc_Controller
    from AMRL.Environment.get_atom_coordinate import pixel_to_nm

except Exception as e:
    print(e)
    print("Createc_Controller is None")
    print(path.basename(__file__) + " will not work")

    Createc_Controller = None
from scipy.spatial.transform import Rotation

def minimize_manipulation_distance(CC: Createc_Controller,
                                   C: CircCorral,
                                   perfect_circle: list,
                                   center: list,
                                   dispmat: list,
                                   showfig: bool=True
                                   ) -> list:
    """

    """

    offset_nm = CC.get_xy_nm()
    len_nm = CC.get_len_nm()
    # create 3D arrays to align vectors using algorithm
    # subtract the center
    centroids_3D = np.array(C.centroids) - np.array(center)

    # need to order centroids in same order as perfect circle to make Kabsch algo work
    refvec = [-1, 0]
    clockwise_angle = lambda point: clockwiseangle(point, [0,0], refvec)
    centroids_3D = np.array(sorted(centroids_3D, key=clockwise_angle))

    perfect_circle_3D = np.array(perfect_circle) - np.array(center)

    centroids_3D = np.append(centroids_3D.T, [np.zeros(len(C.centroids)).T], axis=0).T
    perfect_circle_3D = np.append(perfect_circle_3D.T, [np.zeros(len(C.centroids)).T], axis=0).T

    R, rmsd_ = Rotation.align_vectors(centroids_3D, perfect_circle_3D)

    applied_rotation = R.apply(perfect_circle_3D)[:,0:2].T
    applied_rotation = [x + center for x in applied_rotation.T]

    centroids_3D = np.array([c + center for c in centroids_3D[:,0:2]])

    if showfig:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8,8))
        plt.imshow(dispmat);

        plt.scatter(*np.array(centroids_3D).T, label="atom positions", s=4)
        plt.scatter(*np.array(applied_rotation).T, label="perfect circle", s=50)
        plt.legend()

    for n,c in enumerate(centroids_3D):
        #print(*centroids_3D[:][n], applied_rotation[n])
        initial = pixel_to_nm(centroids_3D[:,][n], dispmat, offset_nm, [len_nm, len_nm])

        final = pixel_to_nm(np.array(applied_rotation[n]).T, dispmat[:,:,0], offset_nm, [len_nm, len_nm])

        dx, dy = (np.array(applied_rotation[n]) - np.array(centroids_3D[n]))

        if showfig:
            import matplotlib.pyplot as plt
            plt.text(*c,'%d' %(n), color='blue')
            plt.text(*(applied_rotation[n]-[4,0]),'%d' %(n), color='orange')

            plt.arrow(*centroids_3D[n], dx, dy, color="white")
    if showfig:
        plt.show()
    return applied_rotation, centroids_3D
