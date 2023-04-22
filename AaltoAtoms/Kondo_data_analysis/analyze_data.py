from .read_vertfile import Spec
from ..utils.find_atom_positions import CircCorral, CircCorralData
import pandas as pd
import numpy as np
import os
import scipy
import matplotlib

import matplotlib.pyplot as plt
basepath = r"/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/data"


def show_current_param_fit_result(c):
    matplotlib.rcParams.update({'font.size': 10})
    S = Spec(os.path.join(basepath, c.vertfile))
    S.clip_data(c.dataclipmin, c.dataclipmax)

    C = CircCorralData(os.path.join(basepath, c.datfile), c.datfile, c.chan)
    C.occupied = True
    C.corral = True
    C.subtract_plane()
    kwargs = {'percentile':c.height_percentile, 'edge_cutoff':c.edge_cutoff}
    C.get_region_centroids(**kwargs)
    radius = C.get_corral_radius(1.5, savefig=False, showfig=False)

    S.remove_background(c.fit_order)
    type_fit = c.type_fit if c.type_fit is not None else "default"

    r = S.fit_fano(marker1=c.marker1, marker2=c.marker2,
                   type_fit=type_fit,
                   showfig=True,
                   q_fixed_val=0.3,
                   actual_radius=radius)


def show_waterfall(Co_Co_data):

    def is_standard(c):
        t1 = int(c[4].biasVoltage) == 80
        t2 = int(c[4].FBLogiset) == 1000
        #t3 = int(np.round(c[4].bias_mv[0])) == 80
        return (t1 and t2)

    plt.figure(figsize=(8,8))
    matplotlib.rcParams.update({'font.size': 22})

    # counter for labeling spectra
    counter = 0

    # look at possible color schemes by running help(plt.cm) in interactive python
    colors = plt.cm.copper(np.linspace(0, 1, len(Co_Co_data)))

    # get the value of dIdV at the norm_mv bias
    norm_mv = 7.5
    ns = sum([is_standard(c) for c in Co_Co_data])

    for n,c in enumerate(Co_Co_data):
        if is_standard(c):
            if n !=23 and n!=19 and n!=15 and n!=16 and n!=12:

                # get the value of dIdV at the norm_mv bias
                norm = c[2][np.argmin(np.abs(norm_mv-c[4].bias_mv))]
                plt.plot(c[3], c[2]/norm+ns-counter, color=colors[n], linewidth=5)
                #plt.text(c[3][0], c[2][0]/norm+c[0], n)

                if counter==0:
                    plt.text(50,  c[2][0]/norm +ns-counter+ 0.2, "%1.1lf nm" %(c[0]))
                if counter==2:
                    plt.text(50,  c[2][0]/norm +ns-counter- 0.6, "%1.1lf nm" %(c[0]))
                if counter==3:
                    plt.text(50,  c[2][0]/norm+ns-counter + 0.4, "%1.1lf nm" %(c[0]))
                if counter==6:
                    plt.text(50,  c[2][0]/norm +ns-counter+ 0.2, "%1.1lf nm" %(c[0]))
                if counter==8:
                    plt.text(50,  c[2][0]/norm+ns-counter - 0.25, "%1.1lf nm" %(c[0]))

                counter += 1
    plt.yticks([])
    plt.xlim(-80, 80)
    plt.xlabel("Bias (mV)")
    plt.ylabel(r"$dI/dV$ (a.u.)")
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\Co-Co-spectrum-waterfall.pdf")
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\Co-Co-spectrum-waterfall.png")
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\Co-Co-spectrum-waterfall.svg")


def plot_radial_width_dependence(Co_Ag_data):
    plt.figure(figsize=(9,6))
    d = get_old_Ag_Co_corrals(dist_cutoff_nm=0.1)
    all_Co_Ag_data = np.concatenate([Co_Ag_data[:,0:2].T, np.array([d.radius, d.w])], axis=-1)
    plt.scatter(*np.array(Co_Co_data)[:,0:2].T, c='orange', label="Co walls")
    plt.scatter(*all_Co_Ag_data, c='blue', label="Ag walls")
    plt.xlabel("Corral radius (nm)")
    plt.ylabel("Kondo resonance width (w=FWHM) (mV)")
    plt.legend()
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\w_radius_dependence.png")
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\w_radius_dependence.pdf")
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\w_radius_dependence.svg")
    return all_Co_Ag_data

import pdb
def analyze_data(corrals: list, # list of corralspectrum objects in data_array.py
                 showfig: bool=False, # show the fit figures as they are created
                 savefig: bool=False,
                 fit_type="default",
                 e0_fixed_val=np.nan) -> list:
    """
        Process arrays of corralspectrum objects (defined in data_array.py file)
        - for each element in corralspectrum array:
            - Save files to data directory for storing on github if not there already
            - Get corral radius and calculate wall density
            - Perform background subtractions, etc, get Kondo width from spectrum
            - Append data to array to be further processed
    """
    data = []
    import shutil
    import platform

    basepath = alldatapath
    # enumerate over items in the list
    for n, c in enumerate(corrals):
        # for migrating dat files to data\ folder
        source_dat = os.path.join(basepath, str(c.datfile))
        datfiledir = os.path.dirname(c.datfile)

        dir = os.path.join(alldatapath, datfiledir)
        dir_exists = os.path.exists(dir)

        if not dir_exists:
            os.makedirs(dir)
        dest_path = os.path.join(alldatapath, datfiledir)
        dest_dat = os.path.join(dest_path, os.path.basename(c.datfile))

        if not os.path.exists(dest_dat):
            print("saving to " + dest_dat)
            shutil.copyfile(source_dat, dest_dat)

        # for migrating vert files to data\ folder
        source_vert = os.path.join(basepath, str(c.vertfile))
        vertfiledir = os.path.dirname(c.vertfile)
        dir = os.path.join(alldatapath, vertfiledir)
        dir_exists = os.path.exists(dir)

        if not dir_exists:
            os.makedirs(dir) #need to use makedirs instead of mkdir

        dest_path = os.path.join(alldatapath, vertfiledir)
        dest_vert = os.path.join(dest_path, os.path.basename(c.vertfile))
        if not os.path.exists(dest_vert):
            print("saving to " + dest_vert)
            shutil.copyfile(source_vert, dest_vert)

        print("ANALYZING (zero-indexed) ELEMENT #%d of %d IN ARRAY:" %(n, len(corrals)))

        # get radius from the topography
        C = CircCorralData(source_dat, c.datfile, c.chan)
        C.occupied = True
        C.corral = True
        C.subtract_plane()
        kwargs = {'percentile':c.height_percentile, 'edge_cutoff':c.edge_cutoff}
        C.get_region_centroids(**kwargs)
        radius = C.get_corral_radius(1.5, savefig=False, showfig=False)
        C.calculate_wall_density()

        # get Kondo width from spectrum
        S = Spec(source_vert)
        dIdV = S.dIdV
        bias_mv = S.bias_mv + S.bias_offset
        S.clip_data(c.dataclipmin, c.dataclipmax)
        S.remove_background(c.fit_order)
        type_fit = "default" if c.type_fit is None else c.type_fit

        # if the e0 fixed val is less than marker1
        if e0_fixed_val < c.marker1: #np.abs(c.marker1-e0_fixed_val)<0.5 and"
            # need c._replace() since corralspectrum is an immutable namedtuple
            c = c._replace(marker1=c.marker1 - np.abs(c.marker1-e0_fixed_val) - 2)
        r = S.fit_fano(marker1=c.marker1, marker2=c.marker2,
                       showfig=showfig, savefig=savefig,
                       actual_radius=radius, type_fit=type_fit,
                       e0_fixed_val=e0_fixed_val)
        width = r[0][1] #e0, w, q, a, b, c

        r_message =  "radius: %1.1lf nm, " %(radius)
        w_message = "width %1.1lf mV " %(width)
        print(c.datfile, c.vertfile, r_message, w_message)

        data.append([radius, width, dIdV, bias_mv, S, C, c, r])
    return data


def get_old_Ag_Co_corrals(dist_cutoff_nm: float):
    dir = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/spatial extent Kondo plots/width comparison"
    dir = r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\Kondo corrals fit data"
    Ag_Co_files = os.listdir(dir)
    #Ag_Co_files = [f for f in Ag_Co_files if f[-9:] =="width.txt" or f[-11:]=="default.txt"]
    fname_is_relevant = lambda f: "kondo_fit_param_width" in f and "covs" not in f and ".pdf" not in f
    Ag_Co_files = [f for f in Ag_Co_files if fname_is_relevant(f)]
    d = pd.DataFrame()

    for f in Ag_Co_files:
        print("reading from file %s" %(f))
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

def w(r, jb=0.53, js=0.21, d1=-0.27, A=3.2, k=0.83):
    rhob = 0.27 # 1/eV
    rhos0 = 0.125 # 1/eV
    D = 4480
    return D*np.exp(-1./(jb*rhob+js*rhos0*(1+A*np.cos(2*k*r+d1))))

def fit_and_plot_functional_curve(radius_array: list,
                                  width_array: list,
                                  bounds: list=None,
                                  p0: list=None,
                                  sigma: list=None,
                                  show_Li_fit=True,
                                  show_isolated_Co=True):

    param_dict = {'Jb':'meV',
              'Js':'meV',
              'd1':'',
              'k':'nm^-1',
              }
    if p0 is None:
        p0 = (0.530, 0.210, -0.27, -0.24, 3.2, 0.83)

    if bounds is None:
        bounds = np.array([(0,1), #Jb
                (0,1), #Js
                (-3.14,3.14), #d1
                (0,10), #A
                (0,2), #k
                ]).T

    params, pcov = scipy.optimize.curve_fit(w,
                                      radius_array,
                                      width_array,
                                      p0,
                                      bounds=bounds,
                                      maxfev=10000,
                                      sigma=sigma,
                                      absolute_sigma=False,
                                      loss='cauchy' )
    rng = np.arange(min(list(radius_array)),max(radius_array),0.01)
    line = plt.plot(rng, np.array([w(x,*params) for x in rng]))#, label="Co/Ag corrals (our data)")
    # plt.plot(rng, np.array([w(x=x, d1=1.5, D=4000) for x in rng]), label="Fit changed" )

    if show_Li_fit:
        plt.plot(rng, 2*np.array([w(x=x,D=4480) for x in rng]), label="Co/Co corrals (Li et. al.)" )
        plt.fill_between(rng,2*np.array([w(x=x,D=4480+620) for x in rng]),2*np.array([w(x=x,D=4480-620) for x in rng]),facecolor=['orange'], alpha=0.5)

    if show_isolated_Co:
        plt.hlines(13.414,2,10,linestyle="dashed", label="Isolated Co")

    #https://stackoverflow.com/questions/14581358/getting-standard-errors-on-fitted-parameters-using-the-optimize-leastsq-method-i
    error = []
    print(params)
    for i in range(len(params)):
        try:
            error.append(np.sqrt(pcov[i][i]))
        except:
            error.append(0.00)
    print(error)
    for n, p in enumerate(list(param_dict.keys())):
        print("%s: %lf pm %lf %s" %(p, params[n], error[n], param_dict[p]))
    return params, pcov, line

    #m = max(radius_array)
    #plt.xlim(2.5, )
