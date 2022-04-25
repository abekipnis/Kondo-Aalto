from .read_vertfile import Spec
from ..find_atom_positions import CircCorral, CircCorralData
import pandas as pd
import numpy as np
import os
import scipy
import matplotlib.pyplot as plt

basepath = "Y:/labdata/Createc/STMDATA/Ag(111)/2022-03 Co Kondo corrals"

def analyze_data(corrals: list, showfig: bool=False) -> list:
    data = []
    for n, c in enumerate(corrals):
        # get the radius from the topography
        C = CircCorralData(os.path.join(basepath, c.datfile), c.datfile, c.chan)
        C.occupied = True
        C.corral = True
        C.subtract_plane()
        C.get_region_centroids(percentile=c.height_percentile, edge_cutoff=c.edge_cutoff)
        radius = C.get_corral_radius(1.5, savefig=False, showfig=False)
        C.calculate_wall_density()

        # get the width from the spectrum
        S = Spec(os.path.join(basepath, c.vertfile))
        dIdV = S.dIdV
        bias_mv = S.bias_mv + S.bias_offset
        S.clip_data(c.dataclipmin, c.dataclipmax)
        S.remove_background(c.fit_order)

        r = S.fit_fano(marker1=c.marker1, marker2=c.marker2,
                       showfig=showfig, actual_radius=radius,type_fit="wtimes")
        width = r[0][1] #e0, w, q, a, b, c

        r_message =  "radius: %1.1lf nm, " %(radius)
        w_message = "width %1.1lf mV " %(width)
        print(c.datfile, c.vertfile, r_message, w_message)

        data.append([radius, width, dIdV, bias_mv, S, C])
    return data

def get_old_Ag_Co_corrals(dist_cutoff_nm: float):
    dir = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/spatial extent Kondo plots/width comparison"
    dir = r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\Kondo corrals fit data"
    Ag_Co_files = os.listdir(dir)
    Ag_Co_files = [f for f in Ag_Co_files if f[-9:] =="width.txt" or f[-11:]=="default.txt"]
    d = pd.DataFrame()

    for f in Ag_Co_files:
        data = pd.read_csv(os.path.join(dir,f), delimiter="\t", index_col=False)
        # data = data[data["radius"]<2.5]
        qcutoff = 1.5
        cut_data = data[np.abs(data["q"])<qcutoff]
        cut_data = cut_data[cut_data["dist"]<dist_cutoff_nm]
#        plt.scatter(cut_data["radius"],cut_data["w"],
#                    s=50/cut_data["dist"], alpha=0.2)

                    #s=100*(cut_data["a"]/np.max(cut_data["a"]
        # plt.scatter(data["radius"],data["w"],s=100*(data["a"]/np.max(data["a"])))

        d = pd.concat([d,cut_data])
    return d

def fit_and_plot_functional_curve(radius_array, width_array, bounds=None, p0=None):
    def w(x=1, jb=0.53, js=0.21, d1=-0.27, d2=-0.24, alpha=0.88, A=3.2, k=0.83, D=4480):
        rhob=0.27
        rhos0=0.125
        return D*np.exp(-1./(jb*rhob+js*rhos0*(1+A*np.cos(2*k*x+d1+d2)/(k*x)**alpha)))
    rng = np.arange(2.5,4.8,0.01)

    if p0 is None:
        p0 = (0.530, 0.210, -0.27,-0.24,0.88, 3.2, 0.83)

    if bounds is None:
        bounds = np.array([(0,1), #Jb
                (0,1), #Js
                (-3.14,3.14), #d1
                (-3.14,3.14), #d2
                (0,2), #alpha
                (0,10), #A
                (0,2) #k
                ]).T

    params = scipy.optimize.curve_fit(w, radius_array, width_array,
                                      p0,
                                      bounds=bounds,
                                      maxfev=6000 )
    rng = np.arange(min(list(radius_array)),10,0.01)
    plt.plot(rng, np.array([w(x,*params[0]) for x in rng]), label="Co/Ag corrals (our data)")
    # plt.plot(rng, np.array([w(x=x, d1=1.5, D=4000) for x in rng]), label="Fit changed" )

    rng = np.arange(2.88,10,0.01)

    plt.plot(rng, 2*np.array([w(x=x,D=4480) for x in rng]), label="Co/Co corrals (Li et. al.)" )
    plt.fill_between(rng,2*np.array([w(x=x,D=4480+620) for x in rng]),2*np.array([w(x=x,D=4480-620) for x in rng]),facecolor=['orange'], alpha=0.5)
    plt.hlines(13.414,2,10,linestyle="dashed", label="Isolated Co")

    print("jb: %lf meV" %params[0][0])
    print("js %lf mev" %params[0][1])
    print("d1 %lf" %params[0][2])
    print("d2 %lf" %params[0][3])
    print("alpha %lf" %params[0][4])
    print("A %lf mV" %params[0][5])
    print("k %lf nm^-1"% params[0][6])

    plt.xlabel("corral radius (nm)")
    plt.ylabel("FWHM (mV)")
    plt.legend(fontsize="small")
    plt.xlim(2.0, 10)
