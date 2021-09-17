import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
import pint
import cmath
import numpy as np

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

hbar = 1.0545718e-34 * ureg.joule * ureg.second
m_electron = 9.109e-31 * ureg.kg
electron_charge = 1.6e-19 * ureg.coulomb

def a(r, k, d0, a0):
    return (2/(np.pi*k*r))**0.5*np.exp(np.pi/4j)*((a0*np.exp(2j*d0)-1)/2j)*np.exp((k*r*1j))

def at(r, k):
    return (2/(np.pi*k*r))**0.5*np.exp(k*r*1j-np.pi/4)

def E(k,m_e, E0):
    E = E0 + hbar**2*k**2/(2*m_e)
    return E

def k(E, m_e, E0):
    K = (2*m_e*(E-E0)/(hbar**2))**.5#a
    return K

def LDOS_at_point(x,y):
    aT = [at(Q_(np.linalg.norm([x,y]-atomloc.magnitude),"nm"),k_tip.to("1/nm")) for atomloc in atom_locs]
    a0 = [a(Q_(np.linalg.norm([x,y]-atomloc.magnitude),"nm"),k_tip.to("1/nm"),1,1) for atomloc in atom_locs]
    ld = np.dot(np.dot(np.array(aT),np.linalg.inv(np.ones(n_atoms) -A)),a0).real
    return ld

# if __name__=="__main__":
n_atoms = 8
radius = 2.5
atom_locs = Q_(np.asarray([[radius * np.cos(n*2*np.pi/n_atoms),
                    radius * np.sin(n*2*np.pi/n_atoms)]
                    for n in range(1,n_atoms+1)]), 'nm')

A = np.zeros((n_atoms,n_atoms))
E = Q_(0.005,"volt")*electron_charge
m_e = 0.4*m_electron
E_0 = Q_(-0.067, "volt")*electron_charge
k_tip = k(E, m_e, E_0)
for n in range(n_atoms):
    for m in range(n_atoms):
        A[n][m] = a(Q_(np.linalg.norm(
                        atom_locs[n].magnitude-atom_locs[m].magnitude)
                        ,"nm"),
                        k_tip.to("1/nm"),
                        1, # d0
                        1) # a0
        if n==m:
            A[n][m] = 1



a0 = np.zeros(n_atoms)
aT = np.zeros(n_atoms)
n_sites = 20
m = Q_(np.asarray(np.linspace(-radius, radius, n_sites)),"nm")
X, Y = np.meshgrid(m, m)
LDOS = np.zeros((n_sites,n_sites))

for n in range(n_sites):
    for m in range(n_sites):
        LDOS[n][m] = LDOS_at_point(X[n][m].magnitude,Y[n][m].magnitude)
plt.imshow(LDOS)
(delayed(LDOS_at_point)(X[n][m].magnitude, Y[n][m].magnitude) )


# create a grid over the space
plt.scatter(LDOS)

for a in n_atoms():

a_vec = []
at_vec = [at() for ]

a(1,1,-np.inf,0)
plt.plot(np.cos(2*np.pi*np.arange(0,1,0.01)))


np.exp(-2*1j*np.inf)
