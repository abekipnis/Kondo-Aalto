# %% codecell
import os, scipy
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
import matplotlib
from AaltoAtoms import CircCorralData, Spec, analyze_data, get_old_Ag_Co_corrals, fit_and_plot_functional_curve
from AaltoAtoms.Kondo_data_analysis.analyze_data import show_current_param_fit_result, plot_radial_width_dependence
from AaltoAtoms import show_waterfall
from AaltoAtoms.utils import labellines
import AaltoAtoms
from AaltoAtoms.Kondo_data_analysis.analyze_data import basepath
from multiprocessing import Pool
import pickle
from importlib import reload
# .dat file and corresponding .VERT file for central Co atom fitting
import data_array
reload(data_array)
from data_array import Co_Co_corrals, Co_Ag_corrals
#
import matplotlib

c = Co_Co_corrals[-10]


def show_current_param_fit_result(c):
    matplotlib.rcParams.update({'font.size': 10})
    S = Spec(os.path.join(basepath, c.vertfile))
    S.clip_data(c.dataclipmin, c.dataclipmax)

    C = CircCorralData(os.path.join(basepath, c.datfile), c.datfile, c.chan)
    C.occupied = True
    C.corral = True
    C.subtract_plane()
    C.get_region_centroids(percentile=c.height_percentile, edge_cutoff=c.edge_cutoff)
    radius = C.get_corral_radius(1.5, savefig=False, showfig=False)

    S.remove_background(c.fit_order)
    type_fit = c.type_fit if c.type_fit is not None else "default"
    r = S.fit_fano(marker1=c.marker1, marker2=c.marker2,
                   type_fit=type_fit,
                   showfig=True,
                   q_fixed_val=np.nan,
                   actual_radius=radius)
    return r

res = show_current_param_fit_result(c)

C = CircCorralData(os.path.join(basepath, c.datfile), c.datfile)
C.occupied = True
C.corral = True
C.subtract_plane()
C.get_region_centroids(percentile=97, edge_cutoff=0.01)
radius = C.get_corral_radius(1.5, savefig=False)
2*np.pi*radius/len(C.centroids)
plt.imshow(C.im)

#
# c = Co_Ag_corrals[6]
#
# S = Spec(os.path.join(basepath, c.vertfile))
# S.clip_data(-25,50)
# S.remove_background(3)
# r = S.fit_fano(marker1=-12, marker2=35, type_fit="default")
#
# C = CircCorralData(os.path.join(basepath, c.datfile), c.datfile, chan=0)
# C.occupied = True
# C.corral = True
# C.subtract_plane()
# C.get_region_centroids(percentile=98, edge_cutoff=0.01)
# radius = C.get_corral_radius(1.5, savefig=False)

def create_waterfall():
    dir = r"Y:\labdata\Createc\STMDATA\Ag(111)\2022-03 Co Kondo corrals\04-11 Ag Co"
    wfall = {"A220411.133438.dat": "A220411.134351.L0013.VERT",
    "A220411.141241.dat": "A220411.141923.L0017.VERT",
    "A220411.145437.dat": "A220411.145852.VERT",
    "A220411.153007.dat": "A220411.153513.VERT",
    "A220411.161126.dat": "A220411.161806.L0017.VERT",
    "A220411.165133.dat": "A220411.165336.VERT",
    "A220411.173719.dat": "A220411.174011.VERT",
    "A220411.183528.dat": "A220411.183838.VERT",
    "A220411.193017.dat": "A220411.193232.VERT",
    "A220411.200858.dat": "A220411.201104.VERT",
    "A220411.204552.dat": "A220411.204741.VERT",
    "A220411.215004.dat": "A220411.215940.L0016.VERT",
    #"A220411.222442.dat": "A220411.222845.L0017.VERT",
    #"A220411.233446.dat": "A220411.233625.VERT",
    "A220412.010237.dat": "A220412.010418.VERT"}
    cache = []
    colors = plt.cm.copper(np.linspace(0, 1, len(wfall)))

    for n, dat in enumerate(list(wfall.keys())):
        C = CircCorralData(os.path.join(dir, dat),"")
        C.occupied = True
        C.corral = True
        C.subtract_plane()
        C.get_region_centroids(percentile=99)
        radius = C.get_corral_radius(1.5, savefig=False, showfig=False)

        S = Spec(os.path.join(dir,wfall[dat]))
        norm_val = S.dIdV[0]
        norm_val = S.dIdV[np.argmin(np.abs(7.5-S.bias_mv))]

        #print(wfall[dat])
        #plt.plot(S.dIdV); plt.show()

        cache.append([S.bias_mv, S.dIdV/norm_val + len(wfall.keys())-n*1.05, colors[n], radius])
    return cache

def show_waterfall(cache: list, bias_idx: int=3, dIdV_idx: int =2) -> None:
    plt.figure(figsize=(8,8))

    colors = plt.cm.copper(np.linspace(0, 1, len(cache)))

    for n, c in enumerate(cache):
        r = c[3] #c[5].pix_to_nm(c[5].r)
        plt.plot(c[bias_idx], c[dIdV_idx], color=colors[n], linewidth=4.5, label="%1.1lf nm" %r)

    # plt.text(47, cache[0][1][0] + 0.2,"%1.1lf nm" %(np.round(cache[0][3],1))) # 11 nm
    # plt.text(47, cache[len(cache)//2][1][0] - 0.6,"%1.1lf nm" %(np.round(cache[len(cache)//2][3],1)))
    # plt.text(47, cache[-1][1][0] - 0.9,"%1.1lf nm" %(np.round(cache[-1][3],1))) # 3.6 nm
    plt.yticks([])
    plt.xlabel("Bias (mV)")
    plt.ylabel(r"$dI/dV$ (a.u.)")
    plt.gcf().axes[0].tick_params(direction="in")

    xvals = [50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50 ]
    labellines.labelLines(plt.gca().get_lines(), align=False,fontsize=12, xvals=xvals)

    plt.xlim(-80,80)
    plt.xticks([-80, -60, -40, -20, 0, 20, 40, 60, 80])

    # plt.legend()
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\Co-Ag-spectrum-waterfall.pdf")
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\Co-Ag-spectrum-waterfall.png")
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\Co-Ag-spectrum-waterfall.svg")

    plt.show()

cache = create_waterfall()
matplotlib.rcParams.update({'font.size': 22})
show_waterfall(cache, 0, 1)

if __name__=="__main__":
    matplotlib.rcParams.update({'font.size': 12})

    Co_Co_data = np.array(analyze_data(Co_Co_corrals, showfig=True))
    Co_Ag_data = np.array(analyze_data(Co_Ag_corrals, showfig=True))
    Co_Ag_data = list(sorted(Co_Ag_data, key=lambda x: -x[0]))
    Co_Co_data = list(sorted(Co_Co_data, key=lambda x: -x[0]))

    with open(r'\\home.org.aalto.fi\kipnisa1\data\Documents\Kondo corrals fits\Co_Co_data.pickle', "wb") as handle:
        pickle.dump(Co_Co_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(r'\\home.org.aalto.fi\kipnisa1\data\Documents\Kondo corrals fits\fits.pickle', "rb") as handle:
        Co_Co_data = pickle.load(handle)


    plt.scatter([c[0] for c in Co_Co_data],  [c[1] for c in Co_Co_data])
    plt.ylim(0, 20)

    atol = 0.2
    radius = 4.5
    dataset = Co_Co_data

    def imshow_dIdV_vs_r():
        f = plt.figure(figsize= (12,12))
        #[c for c in dataset if len(c3)==502]
        # center things with bias offset index
        plt.imshow([list(reversed(c[2]/c[2][np.argmin(abs(c[3]-4))])) for c in dataset if len(c[3])==502 and c[3][0]==80 and c[3][-1]==-80], aspect=20,interpolation=None)
        plt.show()

    imshow_dIdV_vs_r()
    #investigate_radius_range(5.5, 0.2, Co_Ag_data)
    def investigate_radius_range(radius: float, tol: float, dataset: list):
        # get all the corral data with radius within tol of radius
        d = [c for c in dataset if np.isclose([c[0]], [radius], atol=atol)]

        # show the spectra for these corrals
        [plt.plot(c[3], c[2], label="%d atoms" %len(c[5].centroids)) for c in d]
        plt.legend()

        lockinampl = [c[4].LockinAmpl for c in d]
        biasvoltage = [c[4].biasVoltage for c in d]
        setpoint_current = [c[4].FBLogiset for c in d]

    plt.xlim(-20, 20)

    plt.xticks([2,3,4,5,6,7,8])


    len(cache)


    AaltoAtoms.utils.visualizations.create_waterfall()
    plt.plot(Co_Co_data[0][3],Co_Co_data[0][2])
    show_waterfall(Co_Co_data, )

    matplotlib.rcParams.update({'font.size': 12})

    plt.scatter([c[0] for c in Co_Ag_data],  [c[1] for c in Co_Ag_data])

    bounds = {
        'Js': (0,2),
        'Jd': (0,2),
        'd1': (-np.pi, np.pi),
        'd2': (-np.pi, np.pi),
        'alpha': (0, 2),
        'A': (0, 20),
        'k': (0,2)
    }

    p0 = {
        'Js': 0.5,
        'Jd': 0.1,
        'd1': -0.27,
        'd2': -0.24,
        'alpha': 0.88,
        'A': 3.2,
        'k': 0.83
    }


    p0 = [p0[l] for l in list(p0.keys())]
    bounds = np.array([bounds[b] for b in list(bounds.keys())]).T

    d = get_old_Ag_Co_corrals(dist_cutoff_nm=0.1)
    all_Co_Ag_data = np.concatenate([np.array(Co_Ag_data)[:,0:2].T, np.array([d.radius, d.w])], axis=-1)

    fig = plt.figure(figsize=(6,6))

    plt.scatter(*all_Co_Ag_data, label="Ag walls")
    plt.scatter(*np.array([c[0:2] for c in Co_Co_data]).T, label="Co walls")
#    all_Co_Ag_data = plot_radial_width_dependence(Co_Ag_data)
    plt.legend()
    kwargs = {"bounds": bounds, "p0": p0, "show_Li_fit": False, "show_isolated_Co": False}


    f = np.concatenate([np.array([c[0:2] for c in Co_Co_data]).T, all_Co_Ag_data], axis=1)
    fit_and_plot_functional_curve(*f,**kwargs)

    plt.gcf().axes[0].tick_params(direction="in")
    plt.xlabel("Corral radius (nm)")
    plt.ylabel("Central Co atom Kondo resonance width (mV)")
    plt.legend(fontsize="small")


    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\Co-Ag-w-r-fit.pdf")
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\Co-Ag-w-r-fit.png")
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\Co-Ag-w-r-fit.svg")
    fit_and_plot_functional_curve(*np.array(Co_Co_data)[:,0:2].T)

#plt.xlim(-20,20)

#[[c[0], c[-1].file] for c in Co_Ag_data]
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
