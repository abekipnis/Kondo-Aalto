from scipy.special import hankel1
from scipy.constants import hbar, electron_mass, elementary_charge
import numpy as np
import matplotlib.pyplot as plt
import os
delta = 1j*np.inf
mstar = 0.42*electron_mass # *me

alpha = 0.43 # Braun2002
delta = 0.24*np.pi # Braun2002

# look also at Braun 2002, Kliewer 2001
def t(alpha, delta):
    # Crampin2004 equation 2.3.
    # relationship between t and phase shift delta
    return (alpha*np.exp(2*1j*delta)-1)/2

def rho_s(r, E):
    # equation 2.4
    return -1/np.pi*G(r,r,E).imag

def k(E):
    r = np.sqrt(2*mstar*(E+0.067 + 0.000001j)*elementary_charge/(hbar**2))
    return r

def ai(r, rprime, E):
    # equation 2.1
    # free space propagation amplitude of elecron moving in 2D with energy E
    norm = np.linalg.norm(r-rprime)
    #print(r,rprime, norm, k(E)*norm)
    return hankel1(0, k(E)*norm)

# %% Define the atom locations
atom_locs = np.array([[1e-9,0], [2e-9,1e-9]])

# Define the initial scattering from the tip to these atom locations
aT = np.array([ai(atom_loc, tip_pos, E) for atom_loc in atom_locs ])

E = -0.02


# %% Define the A matrix for a given energy
A = np.zeros((len(atom_locs), len(atom_locs)))
for i,a in enumerate(atom_locs):
    for j, b in enumerate(atom_locs):
        if i!=j:
            A[i][j] = ai(a,b,E)*t(alpha,delta)

# one liner to get A matrix for given energy
np.array([[ai(a,b,E)*t(alpha,delta)  if i !=j else 0 for (j, b) in enumerate(atom_locs)  ] for  (i,a) in enumerate(atom_locs) ])

# %% Recreate a spectrum at the tip location
tip_pos = np.array([0,0])
erange =  np.linspace(-0.067,100, 200)
plt.plot(erange, [np.dot(np.dot(aT,np.linalg.inv(1-np.array([[ai(a,b,E)*t(alpha,delta)  if i !=j else 0 for (j, b) in enumerate(atom_locs)] for  (i,a) in enumerate(atom_locs) ]))),np.array([ai(atom_loc, tip_pos, E) for atom_loc in atom_locs ])).real for E in erange])

# %% Create a dIdV map for the given atom positions and energy
A = np.zeros((len(atom_locs), len(atom_locs)))
for i,a in enumerate(atom_locs):
    for j, b in enumerate(atom_locs):
        if i!=j:
            A[i][j] = ai(a,b,E)*t(alpha,delta)
E = 0
nx, ny = 200, 200
minx, maxx, miny, maxy = min(xvals)-1e-9, max(xvals) + 1e-9,  min(yvals)-1e-9, max(yvals) + 1e-9
map = np.zeros((nx,ny))
X, Y = np.meshgrid(np.linspace(minx, maxx, nx),np.linspace(miny, maxy, ny))
for n, nx in enumerate(X):
    for m, ny in enumerate(nx):
        x, y = X[n,m], Y[n,m]
        aT = np.array([ai(atom_loc, [x,y], E) for atom_loc in atom_locs])
        map[m, n ] = np.dot(np.dot(aT,np.linalg.inv(1-A)),np.array([ai(atom_loc, [x,y], E) for atom_loc in atom_locs])).real
plt.imshow(np.flipud(np.rot90(map)), extent=[minx, maxx,maxy, miny])
plt.scatter(*atom_locs.T)

# %% simulate line spectrum
#tip_positions = np.array([[x,4e-9] for x in np.linspace(0,9e-9, 10)])
lspec_dist = np.linalg.norm(tip_positions[0] - tip_positions[1])
spectra = []
erange = np.linspace(-0.067,0.08, 100)
for tp in tip_positions:
    spectra.append([np.dot(np.dot(np.array([ai(atom_loc, tp, E) for atom_loc in atom_locs]), np.linalg.inv(1-np.array([[ai(a,b,E)*t(alpha,delta)  if i!=j else 0 for (j, b) in enumerate(atom_locs)] for  (i,a) in enumerate(atom_locs) ]))),np.array([ai(atom_loc, tp, E) for atom_loc in atom_locs ])).real for E in erange])
spectra = np.array(spectra)
plt.imshow(np.flipud(spectra.T), aspect='auto', extent=[0, lspec_dist,erange[0],  erange[-1]])

# %% Test loading data for a single corral and line spectrum across whole thing
from AaltoAtoms import CircCorralData, Spec
C = CircCorralData(r"Y:\labdata\Createc\STMDATA\Ag(111)\2022-03 Co Kondo corrals\04-11 Ag Co\A220411.160322.dat",'test')
C.get_region_centroids()
atom_locs = C.pix_to_nm(np.array(C.centroids))*1e-9
xvals = [x[0] for x in atom_locs]
yvals = [x[1] for x in atom_locs]
tip_positions = np.array([[x,np.mean(yvals)] for x in np.linspace(min(xvals)-1e-9,max(xvals)+1e-9, 10)])
plt.scatter(*atom_locs.T)
plt.scatter(*tip_positions.T)
np.linalg.norm(atom_locs[0] -atom_locs[1])

# %% Test loading data for a dat file and line spectrum points
dir = r"Y:\labdata\Createc\STMDATA\Ag(111)\2022-03 Co Kondo corrals\04-11 Ag Co"
dat_vert = ["A220411.141241.dat", "141923.L"]
dat = dat_vert[0]
vert = dat_vert[1]

C = CircCorralData(os.path.join(dir,dat),'test')
C.get_region_centroids()
atom_locs = C.pix_to_nm(np.array(C.centroids))*1e-9
xvals = [x[0] for x in atom_locs]
yvals = [x[1] for x in atom_locs]

files = os.listdir(dir)
files = sorted([f for f in files if f[-4:] =="VERT"])
files = [f for f in files if vert in f]
specs = [Spec(os.path.join(dir,f)) for f in files]
xlocs = [s.XPos_nm for s in specs]
ylocs = [s.YPos_nm for s in specs]
x_nm = np.round(C.image_file.size[0]/10.)

xlocs = np.array(xlocs)-C.image_file.offset[0]/10.+x_nm/2.
ylocs = np.array(ylocs)-C.image_file.offset[1]/10.
tip_positions = np.array([xlocs, ylocs]).T*1e-9
plt.scatter(*atom_locs.T)
plt.scatter(*tip_positions.T)

# %% Compare experimental spectra with simulation
plt.plot(specs[len(specs)//2].bias_mv*0.001, specs[-1].dIdV/np.mean(specs[0].dIdV))
plt.plot(erange, np.flipud(spectra)[-1]/np.nanmean(np.flipud(spectra)[-1]))

 #%% Compare the spectra between theory and experiment at one point
 def compare_spectra(spectra_idx: int, atom_locs: list, alpha: float, delta: float):
    erange = specs[spectra_idx].bias_mv[::4]/1000
    x, y = xlocs[spectra_idx]*1e-9, ylocs[spectra_idx]*1e-9
    dIdV =  [np.dot(np.dot(np.array([ai(atom_loc, [x,y], e) for atom_loc in atom_locs]), np.linalg.inv(1-np.array([[ai(a,b,e)*t(alpha,delta)  if i!=j else 0 for (j, b) in enumerate(atom_locs)] for  (i,a) in enumerate(atom_locs) ]))),np.array([ai(atom_loc, [x,y], e) for atom_loc in atom_locs ])).real for e in erange]
#    plt.plot(erange,dIdV/np.nanmean(dIdV))
#    plt.plot(erange, specs[spectra_idx].dIdV[::4]/np.mean(specs[spectra_idx].dIdV))
    norm = np.linalg.norm(np.subtract(dIdV/np.nanmean(dIdV), specs[spectra_idx].dIdV[::4]/np.mean(specs[0].dIdV), out=np.zeros(len(dIdV)), where=~np.isnan(dIdV)))
    return norm
import scipy

compare_spectra(0, atom_locs, -1.4, 0.57)
plt.ylim(-2,2)

scipy.optimize.minimize(lambda x: compare_spectra(0, atom_locs, *x), [1,1])
compare_spectra(0, atom_locs, alpha, delta)
