import os, scipy
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from collections import namedtuple
from AaltoAtoms import CircCorralData, Spec

basepath = "Y:/labdata/Createc/STMDATA/Ag(111)/2022-03 Co Kondo corrals"
fields =  ['datfile', 'height_percentile','vertfile','marker1','marker2', "dataclipmin", "dataclipmax", "fit_order"]
corralspectrum = namedtuple('corralspectrum', fields, defaults=(None,)*len(fields))
# .dat file and corresponding .VERT file for central Co atom fitting
Co_Co_corrals = [
                corralspectrum("03-28/A220328.182736.dat", 99, "03-28/A220328.183230.L0015.VERT", -30, 30), # 2.9 nm
                corralspectrum("03-24/1/A220323.115455.dat", 98, "03-24/1/A220323.115740.VERT", -40,20),# 3.34 nm
                corralspectrum("03-26/A220325.102511.dat", 99, "03-26/A220325.102649.VERT", -5,20),# 4.5 nm
                corralspectrum("03-29/A220329.113009.dat", 99, "03-29/A220329.113106.VERT", -15, 15), # 7.42 nm
                corralspectrum("04-01/A220401.010342.dat", 99, "04-01/A220401.011026.L0125.VERT",0, 50, 0, 50, 3), # 8 nm
                corralspectrum("04-04/A220404.160833.dat", 99, "04-04/A220404.160935.VERT", -10, 20),
                corralspectrum("04-05 6nm Co/A220405.181958.dat", 99, "04-05 6nm Co/A220405.182713.VERT", -10, 20), # 6 nm
                corralspectrum("04-05 6nm Co/A220406.092143.dat", 99, "04-05 6nm Co/A220406.092807.L0010.VERT", -13, 20),
                corralspectrum("4-12 Co-Co/A220412.115834.dat", 99, "4-12 Co-Co/A220412.120450.VERT", -5, 14),
                corralspectrum("4-12 Co-Co/A220412.132510.dat", 99, "4-12 Co-Co/A220412.132822.L0017.VERT", -10, 31),
                corralspectrum("4-12 Co-Co/A220412.154921.dat", 99, "4-12 Co-Co/A220412.160413.VERT", -13, 20),
                corralspectrum("4-12 Co-Co/A220412.163342.dat", 99, "4-12 Co-Co/A220412.163523.VERT", -25, 31),
                corralspectrum("4-12 Co-Co/A220412.183001.dat", 99, "4-12 Co-Co/A220412.184138.L0025.VERT", -10, 20),
                corralspectrum("4-12 Co-Co/A220412.223356.dat", 99, "4-12 Co-Co/A220412.223556.VERT", -13, 20),
                corralspectrum("4-12 Co-Co/A220412.224607.dat", 99, "4-12 Co-Co/A220412.224743.VERT", -13, 31),
                corralspectrum("4-12 Co-Co/A220412.224607.dat", 99, "4-12 Co-Co/A220412.224743.VERT", -20, 20),
                corralspectrum("4-12 Co-Co/A220412.225804.dat", 99, "4-12 Co-Co/A220412.230002.VERT", -10, 20),
                corralspectrum("4-12 Co-Co/A220412.231331.dat", 99, "4-12 Co-Co/A220412.231512.VERT", 0, 21),
                corralspectrum("4-12 Co-Co/A220412.233728.dat", 99, "4-12 Co-Co/A220412.233926.VERT", -13, 31),
                corralspectrum("4-12 Co-Co/A220412.235220.dat", 99, "4-12 Co-Co/A220412.235427.VERT", -13, 31),
                corralspectrum("4-12 Co-Co/A220413.103724.dat", 99, "4-12 Co-Co/A220413.103854.VERT", -13, 31),
                corralspectrum("4-12 Co-Co/A220413.105651.dat", 99, "4-12 Co-Co/A220413.105932.VERT", -13, 31),
                corralspectrum("4-12 Co-Co/A220413.113527.dat", 99, "4-12 Co-Co/A220413.113654.VERT", -13, 20),
                ]
Co_Co_corrals[0]

S = Spec(os.path.join(basepath, c.vertfile))
S.clip_data(-0,50)
S.remove_background(3)
r = S.fit_fano(marker1=-0, marker2=50, type_fit="default")

Co_Co_data = []
for n, c in enumerate([Co_Co_corrals[4]]):
    # c = list(Co_Co_corrals.keys())[0]

    # get the width from the spectrum
    S = Spec(os.path.join(basepath, c.vertfile))
    r = S.fit_fano(marker1=c.marker1, marker2=c.marker2)
    width = r[0][1]

    # get the radius from the topography
    C = CircCorralData(os.path.join(basepath, c.datfile), c.datfile)
    C.occupied = True
    C.corral = True
    C.subtract_plane()
    C.get_region_centroids(percentile=98, edge_cutoff=0.01)
    radius = C.get_corral_radius(1.5, savefig=False)

    Co_Co_data.append([radius, width, S.dIdV, S.bias_mv, S, C])
sslen(Co_Co_data)
S.bias_offset

plt.figure(figsize=(12,6))
for c in Co_Co_data:
    plt.plot(c[3], c[2]/c[2][0]+c[0])
    #plt.plot(c[3], c[2]/c[2][np.argmin(abs(c[3]-7.5))]+c[0])
#plt.xlim(-30,30)


def plot_Ag_Co_corrals():
    dir = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/spatial extent Kondo plots/width comparison"
    Ag_Co_files = os.listdir(dir)
    Ag_Co_files = [f for f in Ag_Co_files if f[-9:] =="width.txt" or f[-11:]=="default.txt"]
    d = pd.DataFrame()

    for f in Ag_Co_files:
        data = pd.read_csv(os.path.join(dir,f), delimiter="\t", index_col=False)
        # data = data[data["radius"]<2.5]
        qcutoff = 1.5
        dcutoff = 0.5
        cut_data = data[np.abs(data["q"])<qcutoff]
        cut_data = cut_data[cut_data["dist"]<dcutoff]
        plt.scatter(cut_data["radius"],cut_data["w"],
                    s=50/cut_data["dist"], alpha=0.2)

                    #s=100*(cut_data["a"]/np.max(cut_data["a"]
        # plt.scatter(data["radius"],data["w"],s=100*(data["a"]/np.max(data["a"])))

        d = pd.concat([d,cut_data])
    return d

# plt.scatter(cut_data["dist"],cut_data["w"])
# plt.xlim(2,3)
#
# plt.scatter(np.sqrt((data["x(nm)"]-data["x(nm)"][0])**2+(data["y(nm)"]-data["y(nm)"][0])**2),data["a"])
# plt.scatter(data["y(nm)"],data["w"])
# d.reset_index()

def fit_and_plot_functional_curve(d):
    radii = set(d["radius"])

    def w(x=1, jb=0.53, js=0.21, d1=-0.27, d2=-0.24, alpha=0.88, A=3.2, k=0.83, D=4480):
        rhob=0.27
        rhos0=0.125
        return D*np.exp(-1./(jb*rhob+js*rhos0*(1+A*np.cos(2*k*x+d1+d2)/(k*x)**alpha)))
    rng = np.arange(2.5,4.8,0.01)

    # plt.plot(np.cos(np.arange(0,2*np.pi,0.01)))

    params = scipy.optimize.curve_fit(w, d.dropna()["radius"], d.dropna()["w"],
                                      p0=(0.530, 0.210, -0.27,-0.24,0.88, 3.2),
                                      bounds=np.array([
                                              # (4000,6000), #D
                                              (0,1), #Jb
                                              (0,1), #Js
                                              # (0,1), #k
                                              (-3.14,3.14), #d1
                                              (-3.14,3.14), #d2
                                              (0.87,0.89), #alpha
                                              (0,10) #A
                                              ]).T,
                                      maxfev=6000 )
    rng = np.arange(min(list(radii)),10,0.01)
    plt.plot(rng, np.array([w(x,*params[0]) for x in rng]), label="Co/Ag corrals (our data)")
    # plt.plot(rng, np.array([w(x=x, d1=1.5, D=4000) for x in rng]), label="Fit changed" )

    rng = np.arange(2.88,10,0.01)

    plt.plot(rng, 2*np.array([w(x=x,D=4480) for x in rng]), label="Co/Co corrals (Li et. al.)" )
    plt.fill_between(rng,2*np.array([w(x=x,D=4480+620) for x in rng]),2*np.array([w(x=x,D=4480-620) for x in rng]),facecolor=['orange'], alpha=0.5)
    plt.hlines(13.414,2,10,linestyle="dashed", label="Isolated Co")

    print("jb: %lf meV" %params[0][0])
    print("js %lf" %params[0][1])
    print("d1 %lf" %params[0][2])
    print("d2 %lf" %params[0][3])
    print("alpha %lf" %params[0][4])
    print("A %lf" %params[0][5])
    # print("A %lf"% params[0][6])
    # print("A %lf" %params[0][7])

    plt.xlabel("corral radius (nm)")
    plt.ylabel("FWHM (mV)")
    plt.legend(fontsize="small")
    plt.xlim(2.0, 10)


plt.figure(figsize=(9,6))
#d = plot_Ag_Co_corrals()
plt.scatter(*np.array(Co_Co_data)[:,0:2].T, c='orange')
fit_and_plot_functional_curve(d)


Co_Ag_corrals = {
    corralspectrum()

}
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
