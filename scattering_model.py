import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
import pint
import cmath
import numpy as np
import multiprocessing
from multiprocessing import Pool
import pdb
from time import time
from multiprocessing import Pool, freeze_support
import matplotlib.pyplot as plt
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

hbar = 1.0545718e-34 * ureg.joule * ureg.second
m_electron = 9.109e-31 * ureg.kg
electron_charge = 1.6e-19 * ureg.coulomb
m_e = 0.4*m_electron
E_0 = Q_(-0.067, "volt")*electron_charge

def a(r, k, d0, a0):
    return (2/(np.pi*k*r))**0.5*np.exp(np.pi/4.*1j)*((a0*np.exp(2j*d0)-1)/2j)*np.exp((k*r*1j))

def at(r, k):
    return (2/(np.pi*k*r))**0.5*np.exp(k*r*1j-1j*np.pi/4.)

def E(k,m_e, E0):
    E = E0 + hbar**2*k**2/(2*m_e)
    return E

def k(E, m_e, E0):
    K = (2*m_e*(E-E0)/(hbar**2))**.5#a
    return K

def LDOS_at_point(x, y, A, k_tip, atom_locs):
    n_atoms = len(atom_locs)
    aT = [at(Q_(np.linalg.norm([x,y]-atomloc.magnitude),"nm"),k_tip.to("1/nm")) for atomloc in atom_locs]
    a0 = [a(Q_(np.linalg.norm([x,y]-atomloc.magnitude),"nm"),k_tip.to("1/nm"),1,1) for atomloc in atom_locs]
    ld = 2*np.dot(np.dot(np.array(aT),np.linalg.inv(np.ones(n_atoms) -A)),a0).real
    return ld

def get_atom_locs(n_atoms, radius): #radius in nm
    atom_locs = Q_(np.asarray([[radius * np.cos(n*2*np.pi/n_atoms),
                        radius * np.sin(n*2*np.pi/n_atoms)]
                        for n in range(1,n_atoms+1)]), 'nm')
    return atom_locs

def create_A_matrix(n_atoms, atom_locs, k_tip):
    A = np.zeros((n_atoms,n_atoms))
    delta0 = np.pi/4
    alpha0 = 0
    for n in range(n_atoms):
        for m in range(n_atoms):
            if n==m:
                A[n][m] = 1
            else:
                A[n][m] = a(Q_(
                    np.linalg.norm(atom_locs[n].magnitude-atom_locs[m].magnitude)
                                ,"nm"),
                    k_tip.to("1/nm"),
                    delta0, # d0
                    alpha0) # a0
    return A

def calc_LDOS(atom_locs, nmxyrange, k_tip):
    n_atoms = len(atom_locs)
    atom_locs -= np.mean(atom_locs, axis=0)
    a0 = np.zeros(n_atoms)
    aT = np.zeros(n_atoms)
    m = Q_(np.asarray(nmxyrange),"nm")
    n_sites = len(nmxyrange)
    X, Y = np.meshgrid(m, m)
    LDOS = np.zeros((n_sites,n_sites))
    A = create_A_matrix(n_atoms, atom_locs, k_tip)
    for n in range(n_sites):
        for m in range(n_sites):
            LDOS[n][m] = LDOS_at_point(X[n][m].magnitude, Y[n][m].magnitude,A, k_tip, atom_locs)
    plt.imshow(LDOS)
    plt.colorbar()
    plt.show()
    return LDOS


def c_LDOS(atom_locs, latt_sites, k_tip):
    n_atoms = len(atom_locs)
    a0 = np.zeros(n_atoms)
    aT = np.zeros(n_atoms)
    m = Q_(latt_sites,"nm")
    LDOS = np.zeros(len(latt_sites))
    A = create_A_matrix(n_atoms, atom_locs, k_tip)

    for n0, n in enumerate(latt_sites):
        LDOS[n0] = LDOS_at_point(n[0], n[1], A, k_tip, atom_locs)

    pdb.set_trace()
    plt.scatter(*np.array(latt_sites).T, c=LDOS)
    plt.colorbar()
    plt.show()
    return LDOS

def gs(atom_locs, latt_sites, erange, spectrumpt):
    s = []
    # where in the lattice array do we need to index into to get the
    # spectrum at the point we want?
    speclatidx = np.argmin([np.linalg.norm(spectrumpt-l) for l in latt_sites])
    atom_locs = Q_(atom_locs, "nm")
    for e in erange:
        E = Q_(e,"volt")*electron_charge
        k_tip = k(E, m_e, E_0)
        LDOS = c_LDOS(atom_locs, latt_sites, k_tip)
        s.append(LDOS[speclatidx])
    return s

def get_LDOS(e, atom_locs, nmxyrange):
    ts = time()
    print(multiprocessing.current_process())
    E = Q_(e,"volt")*electron_charge
    k_tip = k(E, m_e, E_0)
    LDOS = calc_LDOS(atom_locs, nmxyrange, k_tip)
    print(time()-ts)
    return LDOS

def get_spectra(atom_locs, nmxyrange, erange):
    atom_locs = Q_(atom_locs, "nm")
    s = []
    print(multiprocessing.cpu_count())

    p = [(e, atom_locs, nmxyrange) for e in erange]
    with Pool(5) as pool:
        s = pool.starmap(get_LDOS, p)
    return s
