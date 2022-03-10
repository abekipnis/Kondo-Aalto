import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
from scipy.special import jv, yn
import pint
from pint import set_application_registry
import cmath
import numpy as np
import multiprocessing
import pdb
import math
from time import time
from multiprocessing import Pool, freeze_support
from multiprocessing.pool import ThreadPool
import multiprocessing.pool
from itertools import repeat
# from numba import njit
import warnings
warnings.filterwarnings("ignore")
class bcolors:
    OK = '\033[92m' #GREEN
    WARNING = '\033[93m' #YELLOW
    FAIL = '\033[91m' #RED
    RESET = '\033[0m' #RESET COLOR
    BLUE = '\u001b[34m'
    PURPLE = "\u001B[35m"
    CYAN = "\u001B[36m"
# TODO: create main function
# TODO: create class and methods, refactor
# TODO: read from line spectra, location of spectra & .dat file etc.
# TODO: save line spectra data
# TODO: adjust line spectra colorbar limits
# TODO: parallelize line spectra calculation
# TODO: try to simulate real line spectra data (i.e. get real spectra positions in x,y)
# TODO: enable ability to adjust colorscale in simulated line spectra data
# TODO: create line spectrum class object

ureg = pint.UnitRegistry()
set_application_registry(ureg)
Q_ = ureg.Quantity

hbar = 1.0545718e-34 * ureg.joule * ureg.second
m_electron = 9.109e-31 * ureg.kg
electron_charge = 1.6e-19 * ureg.coulomb


class ScatteringModel():
    def __init__(self, atom_locs, ):
        self.m_e = 0.42*m_electron
        self.E_0 = Q_(-0.067, "volt")*electron_charge
        self.atom_locs = atom_locs
        return

    def a(self):
        return

    def at(self):
        return

    def A(self):
        return

    def g0ret(self):
        return

    def si(self):
        return

    def line_spectrum_at_points(self):
        return

    def get_LDOS(self):
        return

def a(r, k, d0, a0):
    """

    a0 and d0 make up η = a0 + id0, the complex number representing the
    abdosrptive channel. (s-wave shift for 2D electron scattering off atom on surface)

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
    r:
    k:

    Returns:
    ________
    """
    return (2/(np.pi*k*r))**0.5*np.exp(k*r*1j-1j*np.pi/4.)

def si(e, a, ib):
    return 1j*4*hbar**2/m_e*(np.exp(2j*(a+1j*ib))-1)

def E(k,m_e, E0):
    """
    Parameters:
    ___________
    k: ureg.Quantity in units of 1/length
    m_e: ureg.Quantity in units of mass
    E0: ureg.Quantity in units of ?

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
    K = (2*m_e*(E-E0)/(hbar**2))**.5
    return K

# plt.plot(np.arange(1e-9, 100e-9, 1e-10),[a(Q_(x,"meter"), k(Q_(-0.065, "volt")*electron_charge, m_e, Q_(-0.067, "volt")*electron_charge), np.pi/2, 0).magnitude.real for x in np.arange(1e-9, 100e-9, 1e-10)])
# plt.plot([a(Q_(x,"meter"), k(Q_(-0.065, "volt")*electron_charge, m_e, Q_(-0.067, "volt")*electron_charge), 10*np.pi, 10).magnitude.imag for x in np.arange(1e-9, 100e-9, 1e-10)])


def LDOS_at_point(x, y, A, kt, atom_locs, n_atoms, d0, alpha0):
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
    ds = [Q_(np.linalg.norm([x,y]-al.magnitude),"nm") for al in atom_locs]
    aT, a0 = np.array([[at(d, kt), a(d, kt, d0, alpha0)] for d in ds]).T
    ones = np.ones(n_atoms)

    # In Fiete2003 its reported in equation 9 as Im(G(r,r,ε))/π
    ld = 2*np.dot(np.dot(np.array(aT),np.linalg.inv(ones-A)),a0).real
    return ld

def get_atom_locs(n_atoms, radius):
    """


    Parameters:
    ___________
    n_atoms: int
    radius: float (in nm)


    Returns:
    ________
    """
    atom_locs = np.asarray([[radius * np.cos(n*2*np.pi/n_atoms),
                        radius * np.sin(n*2*np.pi/n_atoms)]
                        for n in range(1,n_atoms+1)])
    return atom_locs

def create_A_matrix(n_atoms, atom_locs, k_tip, delta0=1.36, alpha0=0):
    """
    Create the atom-atom scattering matrix which stays constant for each energy

    Parameters:
    ___________
    n_atoms: int
    atom_locs:
    k_tip:
    delta0: float
    alpha0: float

    Returns:
    ________
    A: n_atoms*n_atoms array
    """
    A = np.zeros((n_atoms,n_atoms))
    for n in range(n_atoms):
        for m in range(n_atoms):
            # if n==m:
            #     # this number should be infinity (i.e. black dot scatterer)
            #     # but np.inf will return not interesting data?
            #     # so we have to use the biggest number possible (?)
            #     # Fiete2003 equation 18


            # or Crommie1995, equation 5
                # A[n][m] =  1 #- a(Q_(
                    # np.linalg.norm(atom_locs[n].magnitude-atom_locs[m].magnitude),"nm"),
                    # k_tip.to("1/nm"),
                    # delta0, # d0
                    # alpha0) # a0
            # else:
            A[n][m] = at(Q_(
                np.linalg.norm(atom_locs[n].magnitude-atom_locs[m].magnitude),"nm"),
                k_tip.to("1/nm"),
                delta0, # d0
                    alpha0) # a0



            # A[n][m] = si(e, a, ib)*
    return A

def g0ret(rp, r, k):
    dr = np.linalg.norm(rp-r)
    return -1j*(m_e/(2*hbar**2))*(jv(0,k*dr )+1j*yn(0, k*dr))

# @njit(parallel=True)
def calc_LDOS(atom_locs, nmxyrange, k_tip, n_atoms):
    """


    Parameters:
    ___________
    atom_locs:
    nmxyrange:
    k_tip:
    n_atoms: int


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
    # use of zip() & repeat() saves memory by not creating big array every time
    kt = k_tip.to("1/nm")

    p = zip(X.flatten().magnitude,
            Y.flatten().magnitude,
            repeat(A),
            repeat(kt),
            repeat(atom_locs),
            repeat(n_atoms))

    with Pool(5) as pool:
        LDOS = pool.starmap(LDOS_at_point,p)

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
        LDOS_spectrum.append(LDOS_at_point(0, 0, A, k_tip.to("1/nm"), atom_locs, n_atoms))
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

def line_spectrum_at_points(points, atom_locs, erange, delta0=1.36, alpha0=1e16):
    """
    Used for fitting the delta0 scattering phase shift and attenuation alpha0
    for Ag quantum corrals on Ag(111) in order to subtract background spectrum
    for more accurate determination of the Kondo temperature in the quantum confined
    system of the corral.

    Parameters:
    ___________
    points: list (units: nm)
    atom_locs: list (units: nm)
    erange: list (units: V)
    delta: float (units: radians)
    alpha0: float

    Returns:
    ________
    line_spectrum:
    """
    n_atoms = len(atom_locs)
    n_bias = len(erange)
    n_pts = len(points)

    print(bcolors.OK + "Calculating spectra at %d locs, %d energies" %(n_pts, n_bias))
    print(bcolors.CYAN + "delta0: %1.2lf, alpha0: %1.2lf" %(delta0, alpha0) + bcolors.RESET)

    atom_locs = Q_(atom_locs, "nm")

    line_spectrum = []

    for e in erange:
        E = Q_(e,"volt")*electron_charge
        k_tip = k(E, m_e, E_0)
        A = create_A_matrix(n_atoms, atom_locs, k_tip)

        print(bcolors.BLUE + "Calculating LDOS for %d pts, at %1.3lf mV" %(n_pts, e) + bcolors.RESET)

        args = zip(points.T[0],
                    points.T[1],
                    repeat(A),
                    repeat(k_tip),
                    repeat(atom_locs),
                    repeat(n_atoms),
                    repeat(delta0),
                    repeat(alpha0))
        with Pool() as p:
            spec = p.starmap(LDOS_at_point, args)

        line_spectrum.append(spec)
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


# def c_LDOS(atom_locs, latt_sites, k_tip, delta0=np.pi/4, alpha0=0):
#     """
#
#
#     Parameters:
#     ___________
#     atom_locs:
#     latt_sites:
#     k_tip:
#
#
#     Returns:
#     ________
#     """
#     n_atoms = len(atom_locs)
#     a0 = np.zeros(n_atoms)
#     aT = np.zeros(n_atoms)
#     m = latt_sites #nm
#     LDOS = np.zeros(len(latt_sites))
#     A = create_A_matrix(n_atoms, atom_locs, k_tip)
#
#     for n0, n in enumerate(latt_sites):
#         LDOS[n0] = LDOS_at_point(n[0], n[1], A, k_tip, atom_locs, delta0=np.pi/4, alpha0=0)
#
#     plt.scatter(*np.array(latt_sites).T, c=LDOS)
#     plt.colorbar()
#     plt.show()
#     return LDOS
