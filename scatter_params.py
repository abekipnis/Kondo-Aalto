import os, scipy
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from collections import namedtuple
from AaltoAtoms import CircCorralData, Spec, analyze_data, plot_Ag_Co_corrals, fit_and_plot_functional_curve
import AaltoAtoms
from AaltoAtoms.Kondo_data_analysis.analyze_data import basepath
from multiprocessing import Pool

fields =  ['datfile', 'height_percentile','vertfile','marker1','marker2', "dataclipmin", "dataclipmax", "fit_order", "edge_cutoff", "chan"]
corralspectrum = namedtuple('corralspectrum', fields, defaults=(None,)*len(fields))
# .dat file and corresponding .VERT file for central Co atom fitting
Co_Ag_corrals = [
    corralspectrum("04-11 Ag Co\A220412.010237.dat", 97, "04-11 Ag Co\A220412.010418.VERT", -5, 40, -10, 60, None ),
    corralspectrum("04-11 Ag Co\A220411.141643.dat", 99, "04-11 Ag Co\A220411.141923.L0017.VERT", -8, 17, -30, 64, 3),
    corralspectrum("04-11 Ag Co\A220411.145437.dat", 99, "04-11 Ag Co\A220411.145852.VERT", -2, 15),
    corralspectrum("04-11 Ag Co\A220411.153007.dat", 99, "04-11 Ag Co\A220411.153513.VERT", -9, 15, -10, 20, 3),
    corralspectrum("04-11 Ag Co\A220411.153007.dat", 99, "04-11 Ag Co\A220411.154106.L0017.VERT", -8, 12, -10, 21, 3),
    corralspectrum("04-11 Ag Co\A220411.161126.dat", 99, "04-11 Ag Co\A220411.161806.L0017.VERT", 0, 11, -5, 20, 3),
    corralspectrum("04-11 Ag Co\A220411.165133.dat", 99, "04-11 Ag Co\A220411.165336.VERT", -13, 50, -25, 60, 3),
    corralspectrum("04-11 Ag Co\A220411.183528.dat", 99, "04-11 Ag Co\A220411.183838.VERT", -13, 30),
    corralspectrum("04-11 Ag Co\A220411.173719.dat", 99, "04-11 Ag Co\A220411.174011.VERT", -12, 15, -25, 50, 3),
    corralspectrum("04-11 Ag Co\A220411.193017.dat", 99, "04-11 Ag Co\A220411.193232.VERT", -7, 25),
    corralspectrum("04-11 Ag Co\A220411.200858.dat", 99, "04-11 Ag Co\A220411.201104.VERT", -7, 25),
    corralspectrum("04-11 Ag Co\A220411.204552.dat", 99, "04-11 Ag Co\A220411.204741.VERT", -8, 17),
    corralspectrum("04-11 Ag Co\A220411.214417.dat", 99, "04-11 Ag Co\A220411.214626.VERT", -4, 15, -10, 22, 3),
    corralspectrum("04-11 Ag Co\A220411.222442.dat", 99, "04-11 Ag Co\A220411.222845.L0016.VERT", -20, 15, -20, 22, 3),
    corralspectrum("04-11 Ag Co\A220411.233446.dat", 99, "04-11 Ag Co\A220411.233625.VERT", -15, 19, -20, 22, 3),
    corralspectrum("04-11 Ag Co\A220411.233446.dat", 99, "04-11 Ag Co\A220411.233844.L0015.VERT", -15, 14, -20, 22, 3),
    corralspectrum("04-06 6nm Ag walls\A220407.155505.dat", 99, "04-06 6nm Ag walls\A220407.160008.L0071.VERT", -5, 20)
]

Co_Co_corrals = [
    corralspectrum("03-28/A220328.182736.dat", 96, "03-28/A220328.183230.L0015.VERT", -30, 30, None, None, None, 0.3), # 2.9 nm
    corralspectrum("03-24/1/A220323.115455.dat", 97, "03-24/1/A220323.115740.VERT", -30, 19),# 3.34 nm
    corralspectrum("4-12 Co-Co/A220412.231331.dat", 99, "4-12 Co-Co/A220412.231512.VERT", 1, 10, -10,20,3,None), #4.42 nm
    corralspectrum("03-26/A220325.102511.dat", 99, "03-26/A220325.102649.VERT", -5,20),# 4.5 nm
    corralspectrum("04-04/A220404.160833.dat", 99, "04-04/A220404.160935.VERT", -10, 15, -20, 45, 3),
    corralspectrum("04-05 6nm Co/A220405.181958.dat", 99, "04-05 6nm Co/A220405.182713.VERT", -10, 20), # 6 nm
    corralspectrum("04-05 6nm Co/A220406.092143.dat", 99, "04-05 6nm Co/A220406.092807.L0010.VERT", -13, 20),
    corralspectrum("4-12 Co-Co/A220412.115834.dat", 99, "4-12 Co-Co/A220412.120450.VERT", -5, 14),
    corralspectrum("4-12 Co-Co/A220412.132510.dat", 99, "4-12 Co-Co/A220412.132822.L0017.VERT", -10, 31),
    corralspectrum("4-12 Co-Co/A220412.154921.dat", 99, "4-12 Co-Co/A220412.160413.VERT", 0, 50, 0, 50, 3),
    corralspectrum("4-12 Co-Co/A220412.163342.dat", 97, "4-12 Co-Co/A220412.163523.VERT", -25, 31),
    corralspectrum("4-12 Co-Co/A220412.183001.dat", 99, "4-12 Co-Co/A220412.184138.L0025.VERT", -10, 20, None, None, None, 2.5),
    corralspectrum("4-12 Co-Co/A220412.223356.dat", 99, "4-12 Co-Co/A220412.223556.VERT", -18, 15, -70, 60, 3, None, None),
    corralspectrum("4-12 Co-Co/A220412.224607.dat", 99, "4-12 Co-Co/A220412.224743.VERT", -12, 35, -20, 45, 3, None),
    corralspectrum("4-12 Co-Co/A220412.225804.dat", 99, "4-12 Co-Co/A220412.230002.VERT", -10, 20),
    corralspectrum("4-12 Co-Co/A220412.233728.dat", 99, "4-12 Co-Co/A220412.233926.VERT", -20, 40, -20, 60, 3, None, None),
    corralspectrum("4-12 Co-Co/A220412.235220.dat", 98, "4-12 Co-Co/A220412.235427.VERT", -10, 12, -10, 30, 3, None, None),
    corralspectrum("4-12 Co-Co/A220413.103724.dat", 95, "4-12 Co-Co/A220413.103854.VERT", -20, 40, -30, 60, 3, None, None),
    corralspectrum("4-12 Co-Co/A220413.105651.dat", 96, "4-12 Co-Co/A220413.105932.VERT", -13, 31),
    corralspectrum("4-12 Co-Co/A220413.113527.dat", 99, "4-12 Co-Co/A220413.113654.VERT", -13, 20),
    corralspectrum("03-29/A220329.113009.dat", 99, "03-29/A220329.113106.VERT", -15, 15, -20, 20, 3, None, None), # 7.42 nm
    corralspectrum("04-01/A220401.010342.dat", 99, "04-01/A220401.011026.L0125.VERT",0, 50, 0, 50, 3), # 8 nm
    corralspectrum("04-15 Co-Co\A220415.180247.dat", 99, "04-15 Co-Co\A220415.180404.VERT", 0, 12, None, None, None, None, 0), # 8.3 nm
    corralspectrum("04-14 Co Co\A220414.194555.dat", 99, "04-14 Co Co\A220414.195507.VERT", -2, 25, -5, 30, 3, None, 0), #
    corralspectrum("04-14 Co Co\A220414.201501.dat", 98, "04-14 Co Co\A220414.201655.VERT", -2, 11, -5, 20, 3, None, 0), #
    corralspectrum("04-14 Co Co\A220414.202552.dat", 97, "04-14 Co Co\A220414.202911.VERT", -20, 60,  -45, 60, 3, None, 0), #
    corralspectrum("04-14 Co Co\A220414.204203.dat", 99, "04-14 Co Co\A220414.204346.VERT", -10, 30, -40, 70, None, None, 0), #
    corralspectrum("04-14 Co Co\A220414.205921.dat", 99, "04-14 Co Co\A220414.210247.VERT", 0, 13, -10, 60, 1, None, 0), #
    corralspectrum("04-14 Co Co\A220414.212955.dat", 99, "04-14 Co Co\A220414.213310.VERT", -8, 25, None, None, None, None, 0) , #
    corralspectrum("04-17 Co Co\position dependence experiment\A220417.213810.dat", 99, "04-17 Co Co\position dependence experiment\A220417.214221.VERT", -18, 40, -25, 60, 3, None, 0)  #
    ]

c = Co_Co_corrals[-5]

S = Spec(os.path.join(basepath, c.vertfile))
S.clip_data(-45,60)
S.remove_background(3)
r = S.fit_fano(marker1=-20, marker2=60, type_fit="default")

C = CircCorralData(os.path.join(basepath, c.datfile), c.datfile)
C.occupied = True
C.corral = True
C.subtract_plane()
C.get_region_centroids(percentile=98, edge_cutoff=0.01)
radius = C.get_corral_radius(1.5, savefig=False)


plt.figure(figsize=(12,6))
for c in Co_Co_data:
    plt.plot(c[3], c[2]/c[2][0]+c[0])
plt.show()

# with Pool() as P:
#     ret = P.starmap(analyze_data, Co_Ag_corrals)
Co_Co_data = np.array(analyze_data(Co_Co_corrals))
Co_Ag_data = np.array(analyze_data(Co_Ag_corrals))


plt.figure(figsize=(9,6))
d = plot_Ag_Co_corrals(dist_cutoff_nm=0.1)
all_Co_Ag_data = np.concatenate([Co_Ag_data[:,0:2].T, np.array([d.radius, d.w])], axis=-1)
plt.scatter(*Co_Co_data[:,0:2].T, c='orange')
plt.scatter(*all_Co_Ag_data, c='blue')

bounds = {
    'Js': (0,1),
    'Jd': (0,1),
    'd1': (-np.pi, np.pi),
    'd2': (-np.pi, np.pi),
    'alpha': (0, 2),
    'A': (0, 20),
    'k': (0,2)
}

p0 = {
    'Js': 0.53,
    'Jd': 0.21,
    'd1': -0.27,
    'd2': -0.24,
    'alpha': 0.88,
    'A': 3.2,
    'k': 0.83
}

p0 = [p0[l] for l in list(p0.keys())]
bounds = np.array([bounds[b] for b in list(bounds.keys())]).T

fit_and_plot_functional_curve(*all_Co_Ag_data, bounds=bounds, p0=p0)
fit_and_plot_functional_curve(*Co_Co_data[:,0:2].T)


c = Co_Ag_corrals[6]

S = Spec(os.path.join(basepath, c.vertfile))
S.clip_data(-25,50)
S.remove_background(3)
r = S.fit_fano(marker1=-12, marker2=35, type_fit="default")

C = CircCorralData(os.path.join(basepath, c.datfile), c.datfile, chan=0)
C.occupied = True
C.corral = True
C.subtract_plane()
C.get_region_centroids(percentile=98, edge_cutoff=0.01)
radius = C.get_corral_radius(1.5, savefig=False)

plt.figure(figsize=(12,6))
for c in Co_Ag_data:
    plt.plot(c[3], c[2]/c[2][-1]+c[0])

#plt.xlim(-20,20)



# plt.ylim(0,30)
# d.to_csv("/Users/akipnis/Desktop/Kondo_width_scatter.txt")

# plt.savefig("/Users/akipnis/Desktop/Li_fit.png")
# plt.axvline(2.88)h
# for k in np.arange(0.6,2,0.2):
#     plt.plot(rng,10.*np.array([w(x=x,k=k) for x in rng]),label=k)
# plt.plot(rng,np.array([w(x=x,k=k) for x in rng]),label=k)
# plt.legend()
# plt.scatter(d["radius"],d["a"])
# plt.scatter(d["radius"], d["q"])

#further work is needed. fitted value varies slightly from position to position
#ideally want 1 point for each corral
