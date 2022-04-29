# %% codecell
import os, scipy
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
import matplotlib
from AaltoAtoms import CircCorralData, Spec, analyze_data, get_old_Ag_Co_corrals, fit_and_plot_functional_curve
from AaltoAtoms.Kondo_data_analysis.analyze_data import show_current_param_fit_result, plot_radial_width_dependence, show_waterfall
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
c = Co_Co_corrals[1]


show_current_param_fit_result(c)

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

if __name__=="__main__":
    matplotlib.rcParams.update({'font.size': 12})

    Co_Co_data = np.array(analyze_data(Co_Co_corrals, showfig=True))
    with open(r'\\home.org.aalto.fi\kipnisa1\data\Documents\Kondo corrals fits\Co_Co_data.pickle', "wb") as handle:
        pickle.dump(Co_Co_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(r'\\home.org.aalto.fi\kipnisa1\data\Documents\Kondo corrals fits\fits.pickle', "rb") as handle:
        Co_Co_data = pickle.load(handle)

    Co_Co_data = list(sorted(Co_Co_data, key=lambda x: -x[0]))

    plt.scatter([c[0] for c in Co_Co_data],  [c[1] for c in Co_Co_data])
    plt.ylim(0, 20)

    np.isclose([Co_Co_data[0][0]],[8], atol=0.4)
    help(np.isclose)

    atol = 0.2
    radius = 4.5
    dataset = Co_Co_data

    investigate_radius_range(4, 0.2, Co_Ag_data)
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

    show_waterfall(Co_Co_data)

    matplotlib.rcParams.update({'font.size': 12})

    Co_Ag_data = np.array(analyze_data(Co_Ag_corrals, showfig=True))
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
    all_Co_Ag_data = plot_radial_width_dependence()
    fit_and_plot_functional_curve(*all_Co_Ag_data, bounds=bounds, p0=p0)
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
