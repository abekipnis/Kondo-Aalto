import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
import pint
from pint import set_application_registry
import cmath
import numpy as np
import multiprocessing
import pdb
from time import time
from multiprocessing import Pool, freeze_support
from multiprocessing.pool import ThreadPool
import multiprocessing.pool
from itertools import repeat
# from numba import njit
import warnings
warnings.filterwarnings("ignore")

# TODO: create main function
# TODO: create class and methods, refactor
# TODO: read from line spectra, location of spectra & .dat file etc.
# TODO: save line spectra data
# TODO: adjust line spectra colorbar limits
# TODO: parallelize line spectra calculation

ureg = pint.UnitRegistry()
set_application_registry(ureg)
Q_ = ureg.Quantity

hbar = 1.0545718e-34 * ureg.joule * ureg.second
m_electron = 9.109e-31 * ureg.kg
electron_charge = 1.6e-19 * ureg.coulomb
m_e = 0.4*m_electron
E_0 = Q_(-0.067, "volt")*electron_charge

def a(r, k, d0, a0):
    """
    Example:
        a(Q_(0.1,"meter"), k(Q_(-0.065, "volt")*electron_charge, m_e, Q_(-0.067, "volt")*electron_charge), 1, 1).to_reduced_units()

    Parameters:
    ___________
    r: ureg.Quantity in units of length
    k: ureg.Quantity in units of 1/length
    d0: float, unitless
    a0: float, unitless

    Returns:
    ________
    a: dimensionless imaginary number representing the amplitude at a distance
            r from a scattering center with phase shift d0
    """
    return (2/(np.pi*k*r))**0.5*np.exp(np.pi/4.*1j)*((a0*np.exp(2j*d0)-1)/2j)*np.exp((k*r*1j))

def at(r, k):
    """
    Parameters:
    ___________

    Returns:
    ________
    """
    return (2/(np.pi*k*r))**0.5*np.exp(k*r*1j-1j*np.pi/4.)

def E(k,m_e, E0):
    """
    Parameters:
    ___________

    Returns:
    ________
    """
    E = E0 + hbar**2*k**2/(2*m_e)
    return E

def k(E, m_e, E0):
    """
    Example, Fermi wavelength can be calculated from:
    # Fermi wavelength of Ag(111) ~ 7.5 nm
    # 2*np.pi/(k(Q_(0, "volt")*electron_charge, m_e, Q_(-0.067, "volt")*electron_charge).to("1/nm"))

    Parameters:
    ___________
    E:
    m_e:
    E_0:

    Returns:
    ________
    k:
    """
    K = (2*m_e*(E-E0)/(hbar**2))**.5#a
    return K

#
# # checking that things work
# np.arange(0, 2*np.pi, 0.3)
#
# for d in np.arange(0,2*np.pi, 0.3):
#
# plt.plot([a(Q_(x,"meter"), k(Q_(-0.065, "volt")*electron_charge, m_e, Q_(-0.067, "volt")*electron_charge), 10*np.pi, 10).magnitude.imag for x in np.arange(1e-9, 100e-9, 1e-10)])
# plt.plot([a(Q_(x,"meter"), k(Q_(-0.065, "volt")*electron_charge, m_e, Q_(-0.067, "volt")*electron_charge), 10*np.pi, 10).magnitude.imag for x in np.arange(1e-9, 100e-9, 1e-10)])


def LDOS_at_point(x, y, A, k_tip, atom_locs, n_atoms):
    """


    Parameters:
    ___________
    x: array
    y: array
    A:
    k_tip:
    atom_locs:
    n_atoms: int


    Returns:
    ________
    """
    delta0 = np.pi/4. # why use this?
    alpha0 = 0
    aT = [at(Q_(np.linalg.norm([x,y]-atomloc.magnitude),"nm"), k_tip.to("1/nm")) for atomloc in atom_locs]
    a0 = [a(Q_(np.linalg.norm([x,y]-atomloc.magnitude),"nm"), k_tip.to("1/nm"),delta0, alpha0) for atomloc in atom_locs]

    # aT = [at(np.linalg.norm([x,y]-atomloc),k_tip) for atomloc in atom_locs]
    # a0 = [a(np.linalg.norm([x,y]-atomloc),k_tip,delta0,alpha0) for atomloc in atom_locs]
    ones = np.ones(n_atoms)
    ld = 2*np.dot(np.dot(np.array(aT),np.linalg.inv(ones-A)),a0).real
    return ld

def get_atom_locs(n_atoms, radius): #radius in nm
    """


    Parameters:
    ___________
    n_atoms:
    radius:


    Returns:
    ________
    """
    atom_locs = np.asarray([[radius * np.cos(n*2*np.pi/n_atoms),
                        radius * np.sin(n*2*np.pi/n_atoms)]
                        for n in range(1,n_atoms+1)])
    return atom_locs

def create_A_matrix(n_atoms, atom_locs, k_tip, delta_0=np.pi/4., alpha0=0):
    """
    Create the atom-atom scattering matrix which stays constant for each energy

    Parameters:
    ___________
    n_atoms: int
    atom_locs:
    k_tip:
    delta_0: float
    alpha0: float

    Returns:
    ________
    A: n_atoms*n_atoms array
    """
    A = np.zeros((n_atoms,n_atoms))
    for n in range(n_atoms):
        for m in range(n_atoms):
            if n==m:
                # this number should be infinity (i.e. black dot scatterer)
                # but np.inf will return not interesting data?
                # so we have to use the biggest number possible (?)
                A[n][m] = 1000000
            else:
                A[n][m] = a(Q_(
                    np.linalg.norm(atom_locs[n].magnitude-atom_locs[m].magnitude),"nm"),
                    k_tip.to("1/nm"),
                    delta0, # d0
                    alpha0) # a0
    return A

# @njit(parallel=True)
def calc_LDOS(atom_locs, nmxyrange, k_tip, n_atoms):
    """


    Parameters:
    ___________



    Returns:
    ________


    """
    a0 = np.zeros(n_atoms)
    aT = np.zeros(n_atoms)
    m = Q_(np.asarray(nmxyrange),"nm")
    n_sites = len(nmxyrange)
    print("number of grid points in x", n_sites)
    X, Y = np.meshgrid(m, m)
    LDOS = np.zeros((n_sites,n_sites))
    A = create_A_matrix(n_atoms, atom_locs, k_tip)

    # input parameters for LDOS_at_point()

    p = zip(X.flatten().magnitude,
            Y.flatten().magnitude,
            repeat(A),
            repeat(k_tip),
            repeat(atom_locs),
            repeat(n_atoms))

    # p = np.array([[(X[n][m].magnitude, Y[n][m].magnitude, A, k_tip , atom_locs, n_atoms) for n in range(n_sites)] for m in range(n_sites)]).reshape(n_sites*n_sites,6)
    with Pool(5) as pool:
        LDOS = pool.starmap(LDOS_at_point,p)

    return LDOS

def c_LDOS(atom_locs, latt_sites, k_tip):
    """


    Parameters:
    ___________



    Returns:
    ________
    """
    n_atoms = len(atom_locs)
    a0 = np.zeros(n_atoms)
    aT = np.zeros(n_atoms)
    m = latt_sites #nm
    LDOS = np.zeros(len(latt_sites))
    A = create_A_matrix(n_atoms, atom_locs, k_tip)

    for n0, n in enumerate(latt_sites):
        LDOS[n0] = LDOS_at_point(n[0], n[1], A, k_tip, atom_locs)

    plt.scatter(*np.array(latt_sites).T, c=LDOS)
    plt.colorbar()
    plt.show()
    return LDOS

def gs(atom_locs, latt_sites, erange, spectrumpt):
    """


    Parameters:
    ___________



    Returns:
    ________
    """
    s = []
    # where in the lattice array do we need to index into to get the
    # spectrum at the point we want?
    speclatidx = np.argmin([np.linalg.norm(spectrumpt-l) for l in latt_sites])
    atom_locs = atom_locs
    for e in erange:
        E = Q_(e,"volt")*electron_charge
        k_tip = k(E, m_e, E_0)
        LDOS = c_LDOS(atom_locs, latt_sites, k_tip)
        s.append(LDOS[speclatidx])
    return s

def get_LDOS(e, atom_locs, nmxyrange, n_atoms):
    """


    Parameters:
    ___________



    Returns:
    ________
    """
    ts = time()
    print(multiprocessing.current_process())
    E = Q_(e,"volt")*electron_charge
    k_tip = k(E, m_e, E_0)
    LDOS = calc_LDOS(atom_locs, nmxyrange, k_tip, n_atoms)
    print("time to get LDOS map for E %1.3lf:"  %(E.magnitude), time()-ts)
    return np.array(LDOS).reshape(len(nmxyrange),len(nmxyrange))

def get_spectrum_at_middle(atom_locs, erange):
    """


    Parameters:
    ___________



    Returns:
    ________
    """
    n_atoms = len(atom_locs)

    atom_locs -= np.mean(atom_locs, axis=0)
    atom_locs = Q_(atom_locs, "nm")

    LDOS_spectrum = []
    for e in erange:
        E = Q_(e,"volt")*electron_charge
        k_tip = k(E, m_e, E_0)
        A = create_A_matrix(n_atoms, atom_locs, k_tip)
        LDOS_spectrum.append(LDOS_at_point(0, 0, A, k_tip, atom_locs, n_atoms))
    return LDOS_spectrum

def spectrum_along_line(atom_locs, erange):
    """


    Parameters:
    ___________



    Returns:
    ________
    """
    n_atoms = len(atom_locs)

    atom_locs -= np.mean(atom_locs, axis=0)
    atom_locs = Q_(atom_locs, "nm")

    minx = min([a[0] for a in atom_locs])
    maxx = max([a[0] for a in atom_locs])
    line_spectrum = []
    for l in np.arange(minx.magnitude,maxx.magnitude,(maxx.magnitude-minx.magnitude)/100):
        spectrum = []
        for e in erange:
            E = Q_(e,"volt")*electron_charge
            k_tip = k(E, m_e, E_0)
            A = create_A_matrix(n_atoms, atom_locs, k_tip)
            spectrum.append(LDOS_at_point(l, 0, A, k_tip, atom_locs, n_atoms))
        line_spectrum.append(spectrum)
    return line_spectrum

def get_spectra(atom_locs, nmxyrange, erange):
    """
    Parameters:
    ___________



    Returns:
    ________
    """
    s = []
    atom_locs = Q_(atom_locs, "nm")
    n_atoms = len(atom_locs)
    p = [(e, atom_locs, nmxyrange,n_atoms) for e in erange]
    with ThreadPool(5) as pool:
        s = pool.starmap(get_LDOS, p)
    return s

# TODO: put this all into a class (scattering model)
# TODO: try to simulate real line spectra data (i.e. get real spectra positions in x,y)
# TODO: enable ability to adjust colorscale in simulated line spectra data
# TODO: create line spectrum class object
