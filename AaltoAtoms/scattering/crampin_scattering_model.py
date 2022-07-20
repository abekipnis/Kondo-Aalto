from scipy.special import hankel1
from scipy.constants import hbar, electron_mass, elementary_charge
import numpy as np
delta = 1j*np.inf
mstar = 0.42*electron_mass # *me


# look also at Braun 2002, Kliewer 2001
def t(delta):
    # Crampin2004 equation 2.3.
    # relationship between t and phase shift delta
    return -2*hbar**2/mstar*np.exp(1j*delta)*np.sin(delta)

def rho_s(r, E):
    # equation 2.4
    return -1/np.pi*G(r,r,E).imag

def k(E):
    r = np.sqrt(2*mstar*(E+0.067)*elementary_charge/(hbar**2))
    return r

def G0(r, rprime, E):
    # equation 2.1
    # free space propagation amplitude of elecron moving in 2D with energy E
    return -1j*mstar/(2*hbar**2)*hankel1(0, k(E)*np.linalg.norm(r-rprime))

def G(r, rprime, E):
    # equation 2.9
    return G0(r, rprime, E) + \
        np.nansum([
            np.nansum([
            G0(r, Ri, E)*tau_matrix(rs,E)*G0(Rj, rprime, E)
                for i, Ri in enumerate(rs)])
            for j,Rj in enumerate(rs)])

def tau_matrix(rs,E):
    G0prime_matrix = np.zeros([len(rs),len(rs)])
    for i, ri in rs:
        for j, rj in rs:
            G0prime_matrix[i,j]=(1-d(i,j))*G0(rs[i],rs[j],E)
    t_matrix = np.ones([len(rs),len(rs)])*t(delta)
    np.nansum([t_matrix, G0prime_matrix], axis=0)
    tau_matrix = np.linalg.inv(np.nansum([t_matrix,-G0prime_matrix], axis=0))
    return tau_matrix


def d(i,j):
    return 1 if i==j else 0

def tbar(i, j):
    return t(delta)*d(i,j)

def G0primebar(i,j,E):
    ri = rs[i]
    rj = rs[j]
    return (1-d(i,j))*G0(ri, rj, E)

def taubar(G0primeb):
    # equation 2.9
    return np.linalg.inv(np.eye(len(rs))*t(delta)-G0primeb)

rs = np.array([[0,0],[1,0]])

np.matmul(tau_matrix(rs,0), G0(rs[0], [0.1,0.1], 0))
[G0([0.1,0.1], rs[0], 0)*tau_matrix(rs,0)*G0(rs[0], [0.1,0.1], 0).T]
np.nanprod([G0([0.1,0.1], rs[0], 0), tau_matrix(rs,0).all(), G0(rs[0], [0.1,0.1], 0).T], axis=0)

r = np.array([1,1])
rprime = np.array([0,1])
G(r, rprime, 0)

G0(1,0,0)

rho_s(np.array([0.5,1]), -1)
