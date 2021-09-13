import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
def a(r, k, d0, a0):
    return np.sqrt(2/(np.pi*k*r))*np.exp(np.pi/4j)*((a0*np.exp(2j*d0)-1)/2j)*np.exp((k*r*1j))

def at(r, k):
    return np.sqrt(2/(np.pi*k*r))*np.exp(k*r*1j-np.pi/4)


n_atoms = 8
radius = 2.5
atom_locs = np.array([[radius * np.cos(n*2*np.pi/n_atoms),
                    radius * np.sin(n*2*np.pi/n_atoms)]
                    for n in range(1,n_atoms+1)])
A = np.zeros((n_atoms,n_atoms))
atom_locs[0]
for n in range(n_atoms):
    for m in range(n_atoms):
        A[n][m] = a(np.linalg.norm(atom_locs[n]-atom_locs[m]),1,1,1)
A

a_vec = []
at_vec = [at() for ]

splt.scatter(*np.array(atom_locs).T)
a(1,1,-np.inf,0)
plt.plot(np.cos(2*np.pi*np.arange(0,1,0.01)))


np.exp(-2*1j*np.inf)
