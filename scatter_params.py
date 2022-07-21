# %% codecell
# option-shift-enter to Hydrogen - run cell
import os, scipy, pickle
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize
import matplotlib
from AaltoAtoms import CircCorralData, Spec, analyze_data, get_old_Ag_Co_corrals, fit_and_plot_functional_curve
help(fit_and_plot_functional_curve)
from AaltoAtoms.Kondo_data_analysis.analyze_data import show_current_param_fit_result, plot_radial_width_dependence
from AaltoAtoms import show_waterfall, imshow_dIdV_vs_r
from AaltoAtoms.utils import labellines
from AaltoAtoms.Kondo_data_analysis.analyze_data import basepath
from AaltoAtoms.Kondo_data_analysis.analyze_data import w
from importlib import reload
from scipy.signal import decimate
from matplotlib.colors import LogNorm
from scipy.interpolate import interp2d
from itertools import combinations, repeat
import data_array
reload(data_array)

from AaltoAtoms.Kondo_data_analysis import read_vertfile
#read_vertfile = reload(read_vertfile)
#Spec = read_vertfile.Spec

from data_array import Co_Co_corrals, Co_Ag_corrals, corralspectrum

# Create figure 2 (i.e. waterfalls of Co corrals and Ag corrals)
def get_corralspec(datfile: str, l: list) -> list:
    """
    Given datfile string & list of corralspectrum objs (i.e. imported from data_array.py)
    Return the corralspectrum object in the list which has that vertfile as a field
    """
    def is_in(c: corralspectrum, file: str):
        boolean = file in c.datfile
        return boolean
    return [c for c in l if is_in(c, datfile)][0]


def get_total_fit_dIdV(c: corralspectrum) -> list:
    """
    Used in create_waterfall scripts to overlay 3rd order fits w Fano fit
    """
    source_vert = os.path.join(basepath, str(c.vertfile))

    # get Kondo width from spectrum
    S = Spec(source_vert)
    dIdV = S.dIdV
    bias_mv = S.bias_mv + S.bias_offset
    S.clip_data(c.dataclipmin, c.dataclipmax)
    background = S.remove_background(c.fit_order)
    type_fit = "default" if c.type_fit is None else c.type_fit

    r = S.fit_fano(marker1=c.marker1, marker2=c.marker2,
                   showfig=False, actual_radius=None, type_fit=type_fit)

    def in_range(p, m1, m2):
        return p>c.marker1 and p<c.marker2

    enum = enumerate(S.bias_mv)

    # fitting background different than fitting Fano resonance in some cases
    m1, m2 = c.marker1, c.marker2
    n, p = np.array([(n, p) for n, p in enum if in_range(p, m1, m2)]).T

    try: # sometimes there is no (3rd order) background subtraction
        q = np.array(background)[list(map(int, n))]
    except:
        # in these cases we add back nothing
        q = np.zeros(len(n))

    plt.plot(bias_mv, dIdV)
    lin = r[0][4]*np.array([b[1] for b in r[-2]])

    fit_dIdV = r[-1] + q
    fit_bias = [b[1] for b in r[-2]]
    plt.plot(fit_bias, fit_dIdV)
    width = r[0][1] #e0, w, q, a, b, c

    return fit_bias, fit_dIdV


def create_Co_waterfall() -> list:
    """
    Create waterfall of Co Co corrals

    Returns:
    cache: list w elements [bias, didv, color, radius, fit_bias, fit_didv]
    """
    cache = []
    colors = plt.cm.copper(np.linspace(0, 1, len(Co_Co_data)))

    # omit these indices from Co_Co_corrals
    omit_n = [2, 5, 8, 11,13, 17, 22, 23, 27, 28,20,31]

    for n, c in enumerate(Co_Co_data):
        # only show those spectra from -80 to 80 mV
        if np.isclose(80, c[3][0], 0.1):
            #only show those spectra we've decided not to omit
            if n not in omit_n:
                bias = c[3]
                norm_val = c[2][np.argmin([d+7.5 for d in c[2]])]
                dIdV = c[2]/norm_val + c[0]
                radius = c[0]
                d = get_corralspec(c[5].label, Co_Co_corrals)
                fit_bias, fit_dIdV = get_total_fit_dIdV(d)
                fit_dIdV = fit_dIdV/norm_val + c[0]
                cache.append([bias, dIdV, colors[n], radius, fit_bias, fit_dIdV])
    return cache


def show_Co_waterfall(cache: list) -> None:
    """
    Show waterfall of Co Co corrals

    Parameters:
        cache: list w/elements [bias, didv, color, radius, fit_bias, fit_didv]
    """
    plt.figure(figsize=(2, 3))
    matplotlib.rcParams.update({'font.size': 7})
    for n,c in enumerate(cache):
        bias = c[0]
        dIdV = c[1]
        color = c[2]
        r = c[3]
        fit_bias = c[4]
        fit_dIdV = c[5]
        args = bias, dIdV
        kwargs = {"color":color,
                  'label':"%1.1lf nm" %(r),
                  'linewidth':4.5,
                  'alpha':0.8}
        plt.plot(*args, **kwargs)

        plt.plot(fit_bias, fit_dIdV, color="red", alpha=1, linewidth=4)

    # instead of using a bulky legend we can label the lines individually
    # from AaltoAtoms.utils.labellines
    lines = plt.gca().get_lines()
    xvals = map(int,np.linspace(-50,50, len(lines)//2))# [50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50 ]

    labellines.labelLines(lines, align=False, fontsize=7, xvals=xvals)
    plt.yticks([])
    plt.xlabel("Bias (mV)")
    plt.ylabel(r"$dI/dV$ (a.u.)")
    plt.xlim(-80,80)
    plt.axvline(7.3)
    plt.xticks([-80, -60, -40, -20, 0, 20, 40, 60, 80])
    plt.gcf().axes[0].tick_params(direction="in")

    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\images\Co-Co-spectrum-waterfall.pdf")
    #return lines


def create_Ag_cache() -> list:
    """
        Create array with information to plot dIdv for Ag corrals.
    """
    dir = r"Y:\labdata\Createc\STMDATA\Ag(111)\2022-03 Co Kondo corrals\04-11 Ag Co"
    wfall = {"A220411.133438.dat": "A220411.134351.L0013.VERT",
            "A220411.141643.dat": "A220411.141923.L0017.VERT",
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
        dIdV =  S.dIdV/norm_val + len(wfall.keys())-n*1.05
        bias = S.bias_mv
        c = get_corralspec(dat, Co_Ag_corrals)
        fit_bias, fit_dIdV = get_total_fit_dIdV(c)
        fit_dIdV = fit_dIdV/norm_val + len(wfall.keys())-n*1.05
        cache.append([bias, dIdV, colors[n], radius, fit_bias, fit_dIdV])
    return cache

def show_Ag_waterfall(cache: list, bias_idx: int=3, dIdV_idx: int =2) -> None:
    """
    Plot the data created by create_waterfall

    Parameters:
        cache: list of [bias, dIdV + offset, color, radius, fit_bias, fit_dIdV]
        bias_idx:
        dIdV_ idx:

    Returns:
        None
    """
    plt.figure(figsize=(2, 3))
    matplotlib.rcParams.update({'font.size': 7})
    colors = plt.cm.copper(np.linspace(0, 1, len(cache)))

    for n, c in enumerate(cache):
        r = c[3] #c[5].pix_to_nm(c[5].r)
        args = [c[bias_idx], c[dIdV_idx]]
        kwargs = {'color':colors[n],
                  'linewidth':4.5,
                  'label':"%1.1lf nm" %r,
                  'alpha': 0.8}
        plt.plot(*args, **kwargs)

        # also plot the fit on top
        plt.plot(c[-2],c[-1], color="red", alpha=1, linewidth=4 )

    # plt.text(47, cache[0][1][0] + 0.2,"%1.1lf nm" %(np.round(cache[0][3],1))) # 11 nm
    # plt.text(47, cache[len(cache)//2][1][0] - 0.6,"%1.1lf nm" %(np.round(cache[len(cache)//2][3],1)))
    # plt.text(47, cache[-1][1][0] - 0.9,"%1.1lf nm" %(np.round(cache[-1][3],1))) # 3.6 nm
    plt.yticks([])
    plt.xlabel("Bias (mV)")
    plt.ylabel(r"$dI/dV$ (a.u.)")
    plt.gcf().axes[0].tick_params(direction="in")
    plt.axvline(7.3)

    xvals = [50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50 ]
    labellines.labelLines(plt.gca().get_lines(), align=False,fontsize=12, xvals=xvals)

    plt.xlim(-80,80)
    plt.xticks([-80, -60, -40, -20, 0, 20, 40, 60, 80])

    # # plt.legend()
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\images\Co-Ag-spectrum-waterfall.pdf")
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\images\Co-Ag-spectrum-waterfall.png")
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\images\Co-Ag-spectrum-waterfall.svg")

    plt.show()


def show_current_param_fit_result(c: corralspectrum) -> None:
    """
    See how current setup for single corralspectrum works for fit.

    Parameters:
        c: corralspectrum
    Returns:
        None
    """
    basepath = "Z:/Documents/AaltoAtoms/data"
    matplotlib.rcParams.update({'font.size': 10})
    S = Spec(os.path.join(basepath, c.vertfile))
    S.clip_data(c.dataclipmin, c.dataclipmax)
    args = [os.path.join(basepath, str(c.datfile)), c.datfile, c.chan]
    C = CircCorralData(*args)
    C.occupied = True
    C.corral = True
    C.subtract_plane()
    kwargs = {'percentile':c.height_percentile, 'edge_cutoff':c.edge_cutoff}
    C.get_region_centroids(**kwargs)
    radius = C.get_corral_radius(1.5, savefig=False, showfig=True)

    S.remove_background(c.fit_order)
    type_fit = c.type_fit if c.type_fit is not None else "default"
    r = S.fit_fano(marker1=c.marker1, marker2=c.marker2,
                   type_fit=type_fit,
                   showfig=True,
                   e0_fixed_val=np.nan,
                   actual_radius=radius)
    return r

Co_Co_data_loc = os.path.join('data','Co_Co_data.pickle')
Co_Ag_data_loc = os.path.join('data','Co_Ag_data.pickle')
def save_data(e0_fixed_val: float = np.nan) -> None:
    """
    Analyze data defined by corraspectrum arrays in data_array
    Save to location defined in Co_Co_data_loc and Co_Ag_data_loc
    """
    e0_fixed_str = str(e0_fixed_val) if not np.isnan(e0_fixed_val) else ""

    if not np.isnan(e0_fixed_val):
        CC_data_loc = Co_Co_data_loc.strip(".pickle") + "_e0_" + e0_fixed_str+ '.pickle'
        CA_data_loc = Co_Ag_data_loc.strip(".pickle") + "_e0_" + e0_fixed_str+ '.pickle'
    else:
        CC_data_loc = Co_Co_data_loc
        CA_data_loc = Co_Ag_data_loc
    global Co_Co_data
    Co_Co_data = np.array(analyze_data(Co_Co_corrals,
                                        showfig=False,
                                        savefig=False,
                                        e0_fixed_val=e0_fixed_val,
                                        fit_type="default"))
    global Co_Ag_data
    Co_Ag_data = np.array(analyze_data(Co_Ag_corrals,
                                        showfig=False,
                                        savefig=False,
                                        e0_fixed_val=e0_fixed_val,
                                        fit_type="default"))

    # sort by the radius from smallest to largest
    Co_Ag_data = list(sorted(Co_Ag_data, key=lambda x: -x[0]))
    Co_Co_data = list(sorted(Co_Co_data, key=lambda x: -x[0]))

    # save the Cobalt-Cobalt corrals in different pickle files
    with open(CC_data_loc, "wb") as handle:
        pickle.dump(Co_Co_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(CA_data_loc, "wb") as handle:
        pickle.dump(Co_Ag_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

def load_data() -> None:
    """
    Load data saved in Co_Co_data_loc
    """

    global Co_Co_data
    with open(Co_Co_data_loc, "rb") as handle:
        Co_Co_data = pickle.load(handle)

    global Co_Ag_data
    with open(Co_Ag_data_loc, "rb") as handle:
        Co_Ag_data = pickle.load(handle)
#%% load the data from previously saved pickle to save time
load_data()

# %% Look at the mean wall density (atoms/ nm) of corrals
#np.mean(np.append([c[-3].density_per_nm for c in Co_Ag_data],[c[-3].density_per_nm for c in Co_Co_data]))

# %% rerun all data analysis (i.e. Fano fits from data_array.py)
save_data()

# %%
e0range = np.arange(0, 10, 0.25)
for e0 in e0range:
    save_data(e0_fixed_val=7)

# %% Test how params in data_array work to fit Fano to specific spectrum
c = Co_Co_corrals[-1] # select an object from data_array to test
res = show_current_param_fit_result(c) # see how fit works in this case

# %% Create figure 3: i.e. radius-dependence of width
Co_waterfall_cache = create_Co_waterfall()
show_Co_waterfall(Co_waterfall_cache)
Ag_waterfall_cache = create_Ag_cache()
show_Ag_waterfall(Ag_waterfall_cache, 0, 1)

# %% Show the interpolated spectra (i.e. SI figs 3, 4 )
def show_interpolated_spectra(data: str, interpolate:bool=False):
    imshow_dIdV_vs_r(data, downsample=False,
                    interpolate=interpolate, nm_step=0.1,
                    mV_step=0.1,  norm_mV=-67, norm_to_one=True)


fortum for e & insurance
register to telia.


show_interpolated_spectra(Co_Ag_data, interpolate=True)
plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\images\interpolated_Co-Ag_corral_spectra.png", dpi=600)
plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\images\interpolated_Co-Ag_corral_spectra.svg")
plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\images\interpolated_Co-Ag_corral_spectra.pdf")

show_interpolated_spectra(Co_Co_data, interpolate=True)
plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\images\interpolated_Co-Co_corral_spectra.png", dpi=600)
plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\images\interpolated_Co-Co_corral_spectra.svg")
plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\images\interpolated_Co-Co_corral_spectra.pdf")
# %%

# investigate_radius_range(5.5, 0.2, Co_Ag_data)
def investigate_radius_range(radius: float,
                             atol: float,
                             dataset: list) -> list:
    """
    radius: nm
    atol: nm (tolerance around which will get corrals w this radius)
    dataset: list of corralspectra namedtuple objs defined in data_array.py
    """
    # get all the corral data with radius within tol of radius
    d = [c for c in dataset if np.isclose([c[0]], [radius], atol=atol)]

    # show the spectra for these corrals
    s = lambda c: (len(c[5].centroids), c[1], c[0])
    label = lambda c: {'label':"%d atoms, w=%1.1lf mV, r=%1.2lf" %s(c)}
    [plt.plot(c[3], c[2], **label(c)) for c in d]
    #plt.legend()

    lockinampl = [c[4].LockinAmpl for c in d]
    biasvoltage = [c[4].biasVoltage for c in d]
    setpoint_current = [c[4].FBLogiset for c in d]

    return d

retr = investigate_radius_range(4,0.5, Co_Ag_data)

# %% Define functions ot make manuscript Figure 3, make manuscript Figure 3

def w(r: list, jb: float, js: float, d1: float, k: float) -> list:
    """
        Phenomenological model for Kondo resonance width as a function of
        corral radius
    """
    rhob = 0.27 # 1/eV
    rhos0 = 0.125 # 1/eV
    D = 4480 # meV (?)
    A = 1
    return D*np.exp(-1.0/(jb*rhob+js*rhos0*(1+A*np.cos(2*k*r+d1))))


def fit_and_plot_functional_curve(radius_array: list,
                                  width_array: list,
                                  bounds: list=None,
                                  p0: list=None,
                                  sigma: list=None,
                                  show_Li_fit=True,
                                  show_isolated_Co=True) -> tuple:
    param_dict = {'Jb':'meV',
              'Js':'meV',
              'd1':'',
              'k':'nm^-1',
              }
    params, pcov = scipy.optimize.curve_fit(w,
                                      radius_array,
                                      width_array,
                                      p0,
                                      bounds=bounds,
                                      maxfev=100000,
                                      absolute_sigma=True,
                        #              loss="linear",
                        #              f_scale=5,
                        #              method='trf'
                        )
    rng = np.arange(min(list(radius_array)),max(radius_array),0.05)
    plt.plot(rng, np.array([w(x,*params) for x in rng]))#, label="Co/Ag corrals (our data)")
    # plt.plot(rng, np.array([w(x=x, d1=1.5, D=4000) for x in rng]), label="Fit changed" )

    #https://stackoverflow.com/questions/14581358/getting-standard-errors-on-fitted-parameters-using-the-optimize-leastsq-method-i
    error = []
    for i in range(len(params)):
        try:
            error.append(np.sqrt(pcov[i][i]))
        except:
            error.append(0.00)
    for n, p in enumerate(list(param_dict.keys())):
        print("%s: %lf pm %lf %s" %(p, params[n], error[n], param_dict[p]))
    return params, pcov

bounds = {
    'Jb': (0.,0.8),
    'Js': (0,0.7),
    'd1': (0, 2*np.pi),
    'k': (0,2),
}

p0 = {
    'Jb': 0.5,
    'Js': 0.1,
    'd1': 0.27,
    'k': 0.83,
}

p0 = [p0[l] for l in list(p0.keys())]
bounds = np.array([bounds[b] for b in list(bounds.keys())]).T
cm = 1/2.54  # centimeters in inches
fig = plt.figure(figsize=(9*cm,6.6*cm))
ax = fig.add_subplot(111)
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
# subplots((nrows, ncols, index))
ax1 = fig.add_subplot(122)
ax2 = fig.add_subplot(121)#, sharey=ax1)
#fig, (ax1,ax2) = plt.subplots(1,2,, sharey=True) #  fig size for publ. (80 mm x 80 mm)
ax.set_xlabel("r (nm)", size=7)
all_data = np.concatenate([Co_Ag_data, Co_Co_data], axis=0)

# use scatter or errobar ? errorbars are very small for fits
# plt.scatter(*np.array([c[0:2] for c in Co_Ag_data]).T, label="Ag walls", marker="D")
# plt.scatter(*np.array([c[0:2] for c in Co_Co_data]).T, label="Co walls", marker="s")
msize = 2
capsize = 4
elinewidth = 1
kwargs = {'fmt':"o",
          'markersize':msize,
          'solid_capstyle':'butt',
          'capsize':capsize,
          'elinewidth':elinewidth,
          'capthick':elinewidth}
r = np.array([c[0:2] for c in Co_Co_data]).T
ax1.errorbar(*r, yerr=[c[7][1][1][1] for c in Co_Co_data], c="orange", label="Co walls", **kwargs)
r = np.array([c[0:2] for c in Co_Ag_data]).T
ax2.errorbar(*r, yerr=[c[7][1][1][1] for c in Co_Ag_data], c='b',label="Ag walls", **kwargs)
font_kwargs = {'family':'arial', 'style':'normal', 'size':7, 'weight':'normal', 'stretch':'normal'}
_font = matplotlib.font_manager.FontProperties(**font_kwargs)
# ax1.legend(fontsize=7, )
# ax2.legend(fontsize=7, )
kwargs = {"bounds": bounds, "p0": p0, "show_Li_fit": False, "show_isolated_Co": False}
Co = np.array([c[0:2] for c in Co_Co_data]).T
Ag = np.array([c[0:2] for c in Co_Ag_data]).T
f = np.concatenate([Co, Ag], axis=1)
plt.gcf().axes[0].tick_params(direction="in")
plt.gcf().axes[1].tick_params(direction="in")

#https://stackoverflow.com/questions/14581358/getting-standard-errors-on-fitted-parameters-using-the-optimize-leastsq-method-i
w_sigma = [c[7][1][1][1] for c in all_data]
kwargs["sigma"] = w_sigma
nsamples = 200
alpha = 0.005
r_step = 0.05
print("Ag results")
plt.sca(ax2) # set current axis
params, pcov, = fit_and_plot_functional_curve(*Ag,**kwargs)
kwargs = {'mean': params, 'cov': pcov, 'allow_singular': False}
multivar = scipy.stats.multivariate_normal(**kwargs)
rmin = min(Ag[0])
rmax = max(Ag[0])
for i in range(nsamples):
    gauss_sampled_vals = multivar.rvs()
    radius_range = np.arange(rmin, rmax, r_step)
    args = radius_range, w(radius_range, *gauss_sampled_vals)
    ax2.plot(*args, alpha=alpha, color='blue', zorder=0)
plt.gcf().axes[0].tick_params(direction="in")
plt.gcf().axes[1].tick_params(direction="in")
print("Co results")
plt.sca(ax1) # set current axis so fit_and_plot_functional_curve plots on it
kwargs = {"bounds": bounds, "p0": p0, "show_Li_fit": False, "show_isolated_Co": False}
kwargs["sigma"] = w_sigma
params, pcov, = fit_and_plot_functional_curve(*Co,**kwargs)
fit_and_plot_functional_curve(*Co,**kwargs)
kwargs = {'mean':params, 'cov':pcov, 'allow_singular':False}
multivar = scipy.stats.multivariate_normal(**kwargs)
rmin = min(Co[0])
rmax = max(Co[0])
for i in range(nsamples):
    gauss_sampled_vals = multivar.rvs()
    radius_range = np.arange(rmin, rmax, r_step)
    args = radius_range, w(radius_range, *gauss_sampled_vals)
    ax1.plot(*args, alpha=alpha, color='orange', zorder=0)
ax1.set_ylim(2,25)
ax2.set_ylim(2,25)
xticks = ax2.set_xticks(np.arange(2,12,1), fontsize=7,fontproperties=_font)
xticks = ax1.set_xticks(np.arange(2,9,1), fontsize=7,fontproperties=_font)
yticks = ax1.set_yticks(np.arange(5,25,5), fontsize=7, fontproperties=_font)
ylabel = ax2.set_ylabel(r"$\Gamma_0$ (mV)")
ylabel.set_font_properties(_font)
ax2.tick_params(axis='y', which='both', bottom=False, top=False,
                                        labelbottom=False, right=False,
                                        reset=True, labelsize=7)
ax1.tick_params(axis='y', which='both', bottom=False, left=False,
                                        labelleft=False, right=False,
                                        labelright=False)
ax1.tick_params(axis='x',  labelsize=7, direction='in')
ax2.tick_params(axis='x',  labelsize=7, direction='in')
plt.tight_layout()
plt.subplots_adjust(left=None,
                    bottom=None,
                    right=None,
                    top=None,
                    wspace=0.05,
                    hspace=None)
plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\images\Fig3-python-output.png", dpi=600, bbox_inches="tight")
plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\images\Fig3-python-output.pdf", bbox_inches="tight")
#[c[4].fname for c in Co_Co_data if c[1] >17]
# %%

fit_and_plot_functional_curve(*Co,**kwargs)

plt.xlabel("Corral radius (nm)")
plt.ylabel("Central Co atom Kondo resonance width (mV)")
plt.legend(fontsize="small")

# plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\Co-Ag-w-r-fit.svg")
# plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\Co-Ag-w-r-fit.pdf")
# plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\Co-Ag-w-r-fit.png")

fit_and_plot_functional_curve(*np.array(Co_Ag_data)[:,0:2].T)

np.mean([c[-1][0][0] for c in Co_Co_data])

plt.scatter([c[0] for c in Co_Ag_data], [c[-1][0][0] for c in Co_Ag_data])
plt.scatter([c[0] for c in Co_Co_data], [c[-1][0][0] for c in Co_Co_data])
plt.legend(["Ag", "Co"])

plt.scatter([c[0] for c in Co_Ag_data], [c[-1][0][2] for c in Co_Ag_data])
plt.scatter([c[0] for c in Co_Co_data], [c[-1][0][2] for c in Co_Co_data])



# %%
command = '''start matlab -nosplash -nodesktop -r "cd('Z:\Documents\AaltoAtoms\AaltoAtoms\MATLAB\eigenmode_solvers'); plot_eigenspectra(2,10,1)" -logfile log.txt'''
os.system(command)

# %%
