import os
dir = r"Z:\Documents\AaltoAtoms\data\04-11 Ag Co"
dir = r"Y:\labdata\Createc\STMDATA\Ag(111)\2022-03 Co Kondo corrals\04-11 Ag Co"
linespectra = [
    ["A220411.141241.dat", "141923.L"],
    ["A220411.145437.dat", "150029.L"],
    ["A220411.153007.dat", "154106.L"],
    ["A220411.160322.dat", "161806.L"],
    ["A220411.165133.dat", "165929.L"],
    ["A220411.173719.dat", "174224.L"],
    ["A220411.183528.dat", "184121.L"],
    ["A220411.193017.dat", "193446.L"],
    ["A220411.200858.dat", "201552.L"],
    ["A220411.211746.dat", "205321.L"],
    ["A220411.214417.dat", "215940.L"],
    ["A220411.224105.dat", "222845.L"],
    ["A220411.233446.dat", "233844.L"],
    ["A220412.010237.dat", "010708.L"]
]

os.path.join(dir, linespectra[0][0])
from AaltoAtoms import CircCorralData, Spec
import numpy as np
import matplotlib.pyplot as plt
pj = lambda f: os.path.join(dir, f)
data = []

for l in linespectra:
    C = CircCorralData(os.path.join(dir, l[0]), l[0])
    C.subtract_plane()
    C.get_region_centroids()
    C.occupied = True
    C.get_corral_radius(1)

    from AaltoAtoms.utils.visualizations import show_line_spectrum
    files = os.listdir(dir)
    files = sorted([f for f in files if f[-4:] =="VERT"])
    spec_timestamp = l[1]
    files = [f for f in files if spec_timestamp in f]
    specs = [Spec(pj(f)) for f in files]
    normed = [(s.dIdV/np.mean(s.dIdV))[np.argmin(np.abs(6.7-s.bias_mv))] for s in specs]

    rad = C.pix_to_nm(C.r)
    data.append([rad, normed, specs])
    #plt.plot()
data[1][1]
data[1][1] + data[1][0]
[d[2][0].biasVoltage for d in data]
[d[2][0].current[0] for d in data]
[d[2][0].FBLogiset for d in data]
[d[2][0].LockinAmpl for d in data]
data[0][2][len(data[0][2])//2].fit_fano(-5,10, savefig=False)

# looking at equations in Limot2005: surface-state localization at adatoms
def na(e, gamma=10):
    # delta_b and delta_s were fit parameters in Limot2005
    delta_b = 700 # mV
    delta_s = 400 # mV

    # delta_e is the imaginary part of the self energy
    # and it consists of contributions from the coupling of the adsorbate
    # level to the bulk and surface states
    delta_e = delta_b + delta_s
    ea = 500 # renormalized eigenenergy of the adsorbate level
    e0 = -67 # surface state onset (mV) extracted from bare spectrum and equation
    def lambda_e(e, gamma, e0, delta_s):
        return delta_s*np.log((e-e0)**2 + (gamma/2)**2)/(2*np.pi)
    return 1/np.pi*delta_e/((e-ea-lambda_e(e, gamma, e0, delta_s))**2+delta_e**2)

[plt.plot(bias_mv,na(bias_mv, g)) for g in np.arange(7,12)]

bias_mv = data[0][2][len(data[0][2])//2].bias_mv
plt.plot(bias_mv,0.5+ np.arctan((bias_mv+67)/(7/2))/np.pi)
plt.plot(data[0][2][len(data[0][2])//2].bias_mv,data[0][2][len(data[0][2])//2].dIdV)


data[0][2][0].__dict__.keys()
for d in data:
    normed = [(s.dIdV/(s.current[0]))[np.argmin(np.abs(0-s.bias_mv))]*100 for s in d[2]]
    plt.plot( normed + d[0])

[plt.plot(np.array(d[1]) + d[0]) for d in data]
