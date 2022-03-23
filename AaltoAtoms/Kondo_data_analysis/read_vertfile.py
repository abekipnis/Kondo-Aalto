import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib.pyplot import figure

from dataclasses import dataclass
from datetime import datetime
from mpl_toolkits.axes_grid1 import make_axes_locatable
from operator import sub
from scipy import optimize
import scipy.signal
from scipy.interpolate import interp1d

import os
import pdb
import traceback
import pathlib
import inspect
import sys
import time
from itertools import product
import re
import createc
from cmath import sqrt
# import matplotlib
# font = {'family' : 'normal',
#         'weight' : 'normal',
#         'size'   : 14}
# import matplotlib
# matplotlib.rc('font', **font)

kb = 8.617333262145e-5 #eV/K

### TO DO:
# - add locations in nm to the axes of Fano fit line spectrum
# - plot the Fano fit parameters as function of position
# - compare spatial dependence of Fano fit parameters to theory
# - calculate broadening of higher energy states w/ Lorentzian fits
# - fit the surface state onset to equation for surf state electron self-energy (lifetime)
# - save the line spectrum file
# - add the standard deviations into the line spectrum plot
# - put fit values in better location in plot - DONE
# - plot the rest of the fit outside the selection range

class LoadTimeMeta(type):
    base_time = time.perf_counter()
    #time.perf_counter_ns() does not work in some python versions

    def __new__(mcs, name, bases, namespace):
        # print(mcs, name, bases, namespace)
        namespace["__class_load_time__"] = time.perf_counter() - LoadTimeMeta.base_time
        return super().__new__(mcs, name, bases, namespace)

class Spec(metaclass=LoadTimeMeta):
    def __init__(self, fname):
        self.fname = fname
        f = open(fname,"rb")
        d = f.readlines()
        self.a = [str(b).split('\\t') for b in d]
        for n, d in enumerate(self.a):
            if "\\r\\n'" in d:
                d.remove("\\r\\n'")
            if "b'DATA\\r\\n'" in d:
                data_idx = n

        #getting the "mystery line" with some data on it
        # only compatible with STM Version 4 or greater (?)
        n = [n for n,b in enumerate(self.a) if "b'DATA\\r\\n'" in b][0] + 1
        l = self.a[n][0].split(" ")
        l.remove("b'")
        l = [m for m in l if m != '']
        l[-1] = l[-1].strip("\\r\\n\'")
        l = [float(m) for m in l]

        self.NPoints = l[0]
        self.VertPosX_DAC = l[1]
        self.VertPosY_DAC = l[2]
        self.ChannelList = l[3]
        self.Channelcount = l[4]
        self.OutchannelList = l[5]
        self.OutChannelCount = l[6]
        self.XPos_nm = l[7]
        self.YPos_nm = l[8]
        self.data = self.a[data_idx+2:]
        self.T_ADC2 = float([b[0].split("=")[1] for b in self.a[0:n-2] if b[0].split("=")[0]=="b'T_ADC2[K]"][0].strip("\\r\\n'"))
        self.T_ADC3 = float([b[0].split("=")[1] for b in self.a[0:n-2] if b[0].split("=")[0]=="b'T_ADC3[K]"][0].strip("\\r\\n'"))

        # auxadc6 is STM head temperature, auxadc7 is cryostat LHe temp
        self.T_AUXADC6 = float([b[0].split("=")[1] for b in self.a[0:n-2] if b[0].split("=")[0]=="b'T_AUXADC6[K]"][0].strip("\\r\\n'"))
        self.T_AUXADC7 = float([b[0].split("=")[1] for b in self.a[0:n-2] if b[0].split("=")[0]=="b'T_AUXADC7[K]"][0].strip("\\r\\n'"))

        self.LockinAmpl = float([b[0].split("=")[1] for b in self.a[0:n-2] if b[0].split("=")[0]=="b'LockinAmpl"][0].strip("\\r\\n'"))
        self.LockinFreq = float([b[0].split("=")[1] for b in self.a[0:n-2] if b[0].split("=")[0]=="b'LockinFreq"][0].strip("\\r\\n'"))
        self.FBLogiset = float([b[0].split("=")[1] for b in self.a[0:n-2] if b[0].split("=")[0]=="b'FBLogIset"][0].strip("\\r\\n'")) # units: pA
        # print([b[0].split("=")[1] for b in self.a[0:n-2] if b[0].split("=")[0]=="b'BiasVoltage / BiasVolt.[mV]"])
        try:
            self.biasVoltage = float([b[0].split("=")[1] for b in self.a[0:n-2] if b[0].split("=")[0]=="b'BiasVoltage / BiasVolt.[mV]"][0].strip("\\r\\n'").replace(",",'.')) # mV
        except:
            pdb.set_trace()
        try:
            self.bias_mv = np.array([float(d[1].replace(",",".")) for d in self.data][19:])
        except:
            pdb.set_trace()
        # try:
        #     self.current = np.array([float(d[4].replace(",",".")) for d in self.data][19:])
        # except:
        self.current = np.array([float(d[4].replace(",",".")) for d in self.data][19:])
        self.dIdV = np.array([float(d[5].replace(",",".")) for d in self.data][19:])

        try:
            self.bias_offset = interp1d(self.current, self.bias_mv)(0)
            self.bias_mv -= self.bias_offset
        except ValueError as ve:
            print("could not calculate bias offset")
        #correcting the bias offset


    def get_corral_radius(self):
        rdf = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/circle fit plots/circle_fits.txt"
        radii_data = np.loadtxt(rdf, delimiter=",", dtype=object, converters={1: float})

        # pdb.set_trace()

        # specs = [Spec(os.path.join(dpath,f)) for f in files]

        #cff = correctly formatted files
        cff = [str(pathlib.PureWindowsPath(r[0])) for r in radii_data]
        have_radius = self.fname in cff
        if have_radius:
            idx = cff.index(self.fname)
            datfile = radii_data[idx][0]
            # datimg = createc.DAT_IMG(os.path.join(dpath,datfile))
            radius = radii_data[idx][1]
            center_x = float(radii_data[idx][2])
            center_y = float(radii_data[idx][3])
            # plot_line(datimg, specs, [center_x,center_y])
        else:
            pdb.set_trace()
            radius = np.nan

    def fit_fano(self, marker1: float = 0, marker2: float = 0,
                 savefig: bool = True, showfig: bool = True, e0_fixed_val = np.nan, w = np.nan, type_fit: str = "default", init_vals=None, actual_radius=None, dist_to_Co=None) -> list:
        # app.update()

        # TODO: implement a 'quit' function i.e. if we want to stop in the middle
        if marker1==0 and marker2==0: # then create / get new markers
            fig , ax = plt.subplots()
            line = plt.plot(self.bias_mv, self.dIdV)
            lines = [line]
            markers = []
            for n, l in enumerate(lines):
                c = "black" #l[0].get_color()
                # initialize markers at the following two locations
                for b in [-20, 20]:
                    mini = np.argmin(np.abs(np.array(self.bias_mv)-b))
                    ix = l[0].get_xdata()[mini]
                    iy = l[0].get_ydata()[mini]
                    d = DraggableMarker(ax=ax, lines=l,
                                        initx=ix, inity=iy, color=c,
                                        marker=">", dir=dir)
                    markers.append(d)
            plt.title(os.path.split(self.fname)[-1])
            if showfig:
                plt.show()

            marker_vals = sorted([m.marker[0].get_xydata()[0] for m in markers],
                                 key=lambda x: x[0])
            marker1 = marker_vals[0][0]
            marker2 = marker_vals[1][0]

        # e0, w, q, a, b, c
        fixed_vals = [e0_fixed_val, np.nan, np.nan, np.nan, np.nan, np.nan]
        if type_fit == 'default':
            popt, pcov, sb, fit_dIdV, p0 = fit_data_fixed_vals(self.bias_mv, self.dIdV, marker1, marker2, fixed_vals)
        elif type_fit == 'wtimes': #using the width*residual as thing to minimize
            popt, pcov, sb, fit_dIdV, p0 = fit_data_w_times_residual(self.bias_mv, self.dIdV, marker1,marker2, fixed_vals, init_vals=init_vals)
        try:
            fig = figure(figsize=(8.5,6.6)) #width, height (inches)
            a1 = plt.subplot(2,1,1)
            trans = plt.gcf().transFigure
            # pdb.set_trace()
            # if len(pcov.shape)<=1 and np.any(np.isnan(pcov)):
            #     pcov = np.eye(pcov.shape[0])*np.zeros(pcov.shape[0])
            # c = tuple(np.array(list(zip(popt, pcov.diagonal()))).flatten())
            # c = tuple(np.array(list(zip(popt, np.zeros(len(popt))))).flatten())
            c = popt
            if dist_to_Co == None:
                dist_to_Co=0
            if actual_radius==None:
                actual_radius=0
            c = list(tuple([marker1]+[marker2]+
                        [self.biasVoltage]+
                        [self.FBLogiset/1000.0]+
                        [self.bias_offset] + list(c)+p0+[actual_radius]+[dist_to_Co]))
            c.insert(7, c[6]*1e-3/kb)
            c = tuple(c)
            plt.text(0.05, 0.1, "Fit Range:\n     min: %1.2lf mV\n"
                                "     max: %1.2lf mV\n\n"
                                "Spectrum acq. params:\n"
                                "     Bias: %1.0lf mV\n"
                                "     Current setpoint: %1.0lf nA\n"
                                "     Bias offset: %1.3lf mV\n\n"
                                "Fit equation:\n"
                                r"    $\frac{dI}{dV}\sim a(\frac{(q+\epsilon)^2}{1+\epsilon^2}) + bV + c$" "\n"
                                r"     $\epsilon = \frac{V-\epsilon_0}{w}$" "\n"
                                r"    $w=2\gamma$, $\gamma=$HWHM" "\n\n"
                                "Fit parameters:\n"
                                r"    $\epsilon_0$" ': %1.2lf ' r"$mV$" '\n'#$\pm$%1.2lf mV\n'
                                '    w: %1.2lf $mV$, '#\pm$%1.2lf mV, '
                                r"$T_K=$" '%1.2lf $K$\n'
                                '    q: %1.2lf\n'#$\pm$%1.2lf\n'
                                '    a: %1.2lf\n'#$\pm$%1.2lf\n'
                                '    b: %1.2lf\n'#$\pm$%1.2lf\n'
                                '    c: %1.2lf\n\n'#$\pm$%1.2lf\n\n'
                                r"Fit initial guess ($\epsilon_0$,w,q,a,b,c)" "\n"
                                "     %1.2lf, %1.2lf, %1.2lf\n     %1.2lf, %1.2lf, %1.2lf\n\n"
                                "Corral radius: %1.3lf nm\n"
                                "Distance to Co: %1.3lf nm\n"
                                %(c),
                                transform=trans)
            a1.plot(self.bias_mv, self.dIdV)
            a1.plot(sb, fit_dIdV, "b-")
            a1.set_ylim(min(self.dIdV), max(self.dIdV))
            f = fano(sb, *popt)
            # f = fix_T(T)(sb, *popt)
            # plt.plot(bias, fix_T(T)(bias, *popt), "black")

            a1.plot(self.bias_mv, fano(self.bias_mv, *popt), 'r--')
            a1.plot(self.bias_mv, fano(self.bias_mv, *popt)-popt[4]*self.bias_mv,"m--")#-0.2*(min(self.dIdV)-popt[4]),"m--")
            # plt.plot(sb, fano(sb, *popt1),'go', markersize=2)

            residY, residtot = residual(fit_dIdV, f) #NORMALIZE THIS BY # OF FIT POINTS

            plt.subplots_adjust(left=0.4)
            a1.set_xlabel(r"Bias ($mV$)")
            a1.set_ylabel(r"$dI/dV$ (a.u.)")

            t1 = os.path.split(self.fname.split(dpath)[-1])
            t = t1[-1]
            a1.set_title("\n".join(list(t1)))
            a1.legend(["data",r'fit data',"model","model - linear background"])
            a1.axvline(marker1, color="b")
            a1.axvline(marker2, color="b")
            # if savefig:
            #     plt.savefig(file.split(".VERT")[0]+"%s_fano_fit.png" %(t))
            # if showfig:
            #     plt.show()
            a2 = plt.subplot(2, 2, 3)
            a2.plot(sb, residY)
            a2.set_xlabel(r"Bias ($mV$)")
            a2.legend(["Fit residuals"])
            a2.set_ylabel(r"$dI/dV$")

            a3 = plt.subplot(2,2,4)
            a3.hist(residY)
            a3.set_xlabel("Residual histogram")
            a3.yaxis.tick_right()
            a3.set_ylabel("Counts")
            a3.set_xlabel(r"$dI/dV$")
            a3.legend(["Residual histogram"])
            a3.set_xlim(min(residY),max(residY))

            a3.yaxis.set_label_position("right")
            # plt.tight_layout()
            if savefig:
                # pdb.set_trace()
                #fp = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/spatial extent Kondo plots/width comparison"
                fp = os.path.dirname(self.fname)
                fig_path = os.path.join(fp,"%s_fit_residual.pdf" %(t.strip(".VERT")))

                # fig_path = os.path.join(os.path.split(self.fname)[0],"%s_fit_residual.pdf" %(t.strip(".VERT")))
                plt.savefig(fig_path)
            if showfig:
                plt.show()

            # plt.close()
        except Exception:
            plt.close()
            # app.update()
            print(traceback.format_exc())
            # messagebox.showerror("showerror","could not plot %s" %(self.fname))
            return [popt, pcov, marker1, marker2, self.XPos_nm, self.YPos_nm, 0]
        return [popt, pcov, marker1, marker2, self.XPos_nm, self.YPos_nm, residtot]

def plot_lrs():
    fig, (ax1,ax2) = plt.subplots(figsize=(8,8), nrows=2)
    for n,f in enumerate(lrs):
        s = Spec(f)
        # s.get_corral_radius()
        # print(s.__dict__["fname"])
        p = np.array(s.dIdV)/s.dIdV[np.argmin(np.abs(s.bias_mv-8))]
        # print(p)
        ax2.plot(s.__dict__["bias_mv"], p)#, label="%1.2lf nm radius" %(radii[n]))
    ax2.set_xlim([-100,100])
    for n, f in enumerate(srs):
        s = Spec(f)
        p = np.array(s.dIdV)/s.dIdV[np.argmin(np.abs(s.bias_mv-8))]
        # print(s.__dict__["fname"], np.max(s.__dict__["bias_mv"]))

        ax1.plot(s.__dict__["bias_mv"],p)
    # type(s.__dict__["bias_mv"])
    ax1.legend(["2.5nm","3.8nm","4.5nm"], loc="upper right")
    ax2.legend(["2.5nm","3.8nm","4.5nm"])
    plt.xlabel("Bias (mV)")
    ax1.set_ylabel("d$\it{I}$/d$\it{V}$ (a.u.)")
    plt.ylabel("d$\it{I}$/d$\it{V}$ (a.u.)")
    # plt.title("dI/dV spectrum of Kondo resonance on Co at corral center")
    plt.savefig("/Users/akipnis/Desktop/spectra.pdf")
    plt.show()

class DraggableMarker():
    """
    Marker class to create boundaries to stitch PSD plot
    Adapted from here: https://stackoverflow.com/questions/43982250/draggable-markers-in-matplotlib
    to use multiple markers at the same time.
    """

    def __init__(self, ax=None, lines=None, initx=0, inity=0, color='red',marker="o", dir=""):
        if ax == None:
            self.ax = plt.gca()
        else:
            self.ax=ax
        if lines==None:
            self.lines=self.ax.lines
        else:
            self.lines=lines
        self.lines = self.lines[:]
        self.tx =  [self.ax.text(0,0,"") for l in self.lines]
        self.marker = [self.ax.plot([initx],[inity], marker=marker, color=color)[0]  for l in self.lines]

        self.draggable=False
        self.dir = dir
        self.c1 = self.ax.figure.canvas.mpl_connect("button_press_event", self.click)
        self.c2 = self.ax.figure.canvas.mpl_connect("button_release_event", self.release)
        self.c3 = self.ax.figure.canvas.mpl_connect("motion_notify_event", self.drag)

    def click(self,event):
        #left click
        if event.button==1:
            #only want to click and drag the marker that we are hovering around / over
            #calculate distance in log space
            mouse_dist_to_marker = np.sqrt((event.xdata-self.marker[0].get_xdata()[0])**2+(event.ydata-self.marker[0].get_ydata()[0])**2)
            if mouse_dist_to_marker<10:
                self.draggable=True
        self.ax.figure.canvas.draw_idle()

    def drag(self, event):
        #move the marker to the nearest point
        if self.draggable:
            self.update(event)
            self.ax.figure.canvas.draw_idle()

    def release(self,event):
        #stop dragging once we release the moust
        self.draggable = False
        # plt.savefig(self.dir+"/PSD_w_markers.pdf")

    def update(self, event):
        for i, line in enumerate(self.lines):
            x,y = self.get_closest(line, event.xdata)
            self.tx[i].set_position((x,y))
            self.tx[i].set_text("x:%1.2lf\ny:%.1E" %(x,y))
            self.marker[i].set_data([x],[y])

    def get_closest(self,line, mx):
        x,y = line.get_data()
        mini = np.argmin(np.abs(x-mx))
        return x[mini], y[mini]

def plot_fano_fit_line(f):
    d1 = pd.read_csv(f, delimiter="\t", index_col=False)
    try:
        covs = pd.read_csv(f.split(".txt")[0] + "_covs.txt",
                                        delimiter="\t", index_col=False)
    except Exception as e:
        print(e)
    #total length of line in nm:
    len_nm = np.linalg.norm(np.array(d1.iloc[0][["x(nm)","y(nm)"]])-np.array(d1.iloc[-1][["x(nm)","y(nm)"]]))
    dists = [np.linalg.norm(np.array(d1.iloc[0][["x(nm)","y(nm)"]])-np.array(d1.iloc[i][["x(nm)","y(nm)"]])) for i in range(len(d1))]

    fs = d1["file"]
    try:
        specs = [Spec(f) for f in fs]
    except:
        specs = [Spec(os.path.join(dpath,f)) for f in fs]
    specdata = np.flipud(np.rot90([s.dIdV for s in specs]))
    biasdata = list(reversed(specs[0].bias_mv))

    # what are the indices of the line?
    ns = fs.str.extract(r"L([0-9]*).VERT").astype(float)
    d = np.array(d1[["e0","w","q","a","b","c"]])
    im = np.array([fano(biasdata,*n) for n in d]).T
    fig, ((ax1, ax2, ax3, ax8), (ax5, ax6, ax7, ax9 )) = plt.subplots(2,4,figsize=(16,6), sharex="col")
    cax = ax1.matshow(np.flipud(im))
    im = ax1.get_images()
    extent =  im[0].get_extent()
    aspect = 1
    ax1.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

    nbias = np.array(biasdata).shape[0]
    no_labels = 5 # how many labels to see on axis 5
    step_bias = int(nbias / (no_labels - 1)) # step between consecutive labels
    bias_positions = np.arange(0,nbias,step_bias)
    bias_labels = list(map(round,list(reversed(biasdata[::step_bias]))))

    ax1.set_yticks(bias_positions)
    ax1.set_yticklabels(["%d" %(int(b)) for b in bias_labels])

    # ax1.set_xticks(np.arange(0, len(ns), 3))
    # pdb.set_trace()
    # nlabels = 6
    # step_nm = len(dists)/(nlabels-1)
    # step_positions = np.arange(0, len(dists), step_nm)
    # pos_labels = list(reversed(np.array(dists)[::step_nm]))
    # pos_labels = ["%1.1lf" %(p) for p in pos_labels]
    # ax1.set_xticks(step_positions)
    # ax1.set_xticklabels(pos_labels)


    # TODO: Fix labeling of indexing (check 2.5nm pm 20mV line spectra in presentation)
    # TODO: account for if there are some data where fit doesn't work
    # TODO: make one of the x axis labels in nanometers

    # ax1.set_xticklabels(["%d" %(int(d[0])) for d in np.array(ns)[0::3]])

    ax1.tick_params(axis="x", bottom=False, top=True, labelbottom=False, labeltop=True)

    #ax2.errorbar(dists, d1["e0"], yerr=covs["e0"], fmt='o')
    ax2.scatter(dists, d1["e0"])

    mine0 = min(d1["e0"])
    maxe0 = max(d1["e0"])
    #pdb.set_trace()
    e0_10pct =  0.1*(maxe0-mine0)
    if not np.isnan(e0_10pct):
        ax2.set_ylim([mine0-e0_10pct, maxe0+e0_10pct])
    # assert not np.isnan(e0_10pct)
    # assert not np.isnan(mine0)
    # assert not np.isnan(maxe0)
    ax2.set_xlim([0,len_nm])
    ax2.set_ylabel("mV")
    ax6.set_xlabel("Distance (nm)")

    # ax3.errorbar(dists, d1["w"], yerr=covs["w"], fmt='o')
    ax3.scatter(dists, d1["w"])

    maxw = max(d1["w"])
    minw = min(d1["w"])
    w_10pct =  0.1*(maxw-minw)
    # assert not np.isnan(maxw)
    # assert not np.isnan(minw)
    ax3.set_ylim(minw-w_10pct, maxw+w_10pct)
    ax7.set_xlabel("Distance (nm)")
    ax3.set_ylabel("mV")


    ax4 = ax3.twinx()
    # ax4.errorbar(dists, d1["w"]/kb*1e-3, yerr=covs["w"]/kb*1e-3, fmt='o')
    ax4.scatter(dists, d1["w"]/kb*1e-3)
    ax4.set_ylim((min(d1["w"])-w_10pct)/kb*1e-3, (max(d1["w"]+w_10pct))/kb*1e-3)

    # expected Kondo temperature for Co on Ag111 is 92K
    # width = 2k_BT_K

    ax5.matshow(specdata)
    nbias = len(biasdata)
    no_labels = 5 # how many labels to see on axis 5
    step_bias = int(nbias / (no_labels - 1)) # step between consecutive labels
    bias_positions = np.arange(0,nbias,step_bias)
    bias_labels = list(map(round,list(reversed(biasdata[::step_bias]))))
    ax5.set_yticks(bias_positions)
    ax5.set_yticklabels(["%d" %(int(b)) for b in bias_labels])

    # ax5.set_xticklabels(np.array(['']+dists).astype(s tr))
    # if all markers are the same
    ax5.axhline(y=nbias-np.argmin(np.abs(biasdata-d1["marker1"].iloc[0])), color="r")
    ax5.axhline(y=nbias-np.argmin(np.abs(biasdata-d1["marker2"].iloc[0])), color="r")

    # if markers are different
    for n, d in enumerate(d1["marker1"].values):
        ax5.scatter(ns[0][n]-min(ns[0]), nbias-np.argmin(np.abs(biasdata-d1["marker2"].iloc[n])),color="r", s=10)
        ax5.scatter(ns[0][n]-min(ns[0]), nbias-np.argmin(np.abs(biasdata-d)),color="r", s=10)

    im = ax5.get_images()
    extent =  im[0].get_extent()
    ax5.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
    ax5.tick_params(axis="x", bottom=True, top=False, labelbottom=True, labeltop=False)

    # ax6.errorbar(dists, d1["q"], yerr=covs["q"], fmt='o')
    ax6.scatter(dists, d1["q"])

    q10pct = 0.1*(max(d1["q"])-min(d1["q"]))
    ax6.set_ylim([min(d1["q"]) - q10pct, max(d1["q"]) + q10pct])
    ax7.scatter(dists, d1["resid"])

    a1, a2 = fs[0].split(os.path.sep)[-3:-1]
    fig.suptitle("Fano fits to Kondo resonance on corral central atom: %s\%s" %(a1, a2))
    # add red lines across where the fit happens

    ax8.scatter(dists,d1["a"])
    ax9.scatter(dists,d1["b"])
    plt.subplots_adjust(wspace=0.345)
    # have to adjust axes limits since sigma blows up
    # ax2.set_ylim([0, 15])
    # ax3.set_ylim([0, 15])
    # ax4.set_ylim([0, 15/kb*1e-3])
    # ax6.set_ylim([-0.5, 2])

    # set titles for the subplots
    ax1.set_title("Fano fit")
    ax2.set_title(r"$\epsilon_0$ (eV)")
    ax3.set_title(r"Width, $T_K$")
    ax5.set_title("Raw data")
    ax6.set_title("q")
    ax7.set_title("Residuals (au)")
    ax9.set_title("b")
    ax8.set_title("a")

    ax5.set_xlabel("Index")
    ax5.set_ylabel("bias (mV)")
    ax4.set_ylabel("Kelvin")
    ax1.set_ylabel("bias (mV)")

    plt.savefig(f.split(".txt")[0]+".pdf")
    # plt.show()

def fermi_dirac(e, mu, T):
    dist = 1/(np.exp((np.array(e)-mu)/(kb*T*1e3)) + 1)
    return dist

def fano(V, e0, w, q, a, b, c):
    """
    Parameters
    __________

    Returns
    _______
    """
    def eps(V, e0, w):
        return (np.array(V)-e0)/(w/2)
    fit_func = a*((q + eps(V, e0, w))**2/(1+eps(V, e0, w)**2))+ b*np.array(V) + c
    return fit_func

def frota(V, e0, gamma_f, a,b,c):
    return a*(sqrt(1j*gammaf/(V-e0+1j*gammaf))) + b*V + c

def fano_t_broad(T, V, e0, w, q, a, b, c):
    """
    Parameters
    __________

    Returns
    _______
    """
    # padding so there are no edge effects from convolution included in fit
    dV = np.abs(V[0] - V[1])
    padmVs = np.abs(V[0]-V[-1])/2 # padding (in mV) on either side of fit range
    pad_idcs = int(padmVs/dV) #padding in array indices based on dV
    V_pad_minus = [min(V)-dV*(i+1) for i in list(range(pad_idcs))]
    V_pad_plus = [max(V)+dV*(i+1) for i in list(reversed(list(range(pad_idcs))))]
    V_padded = np.concatenate([V_pad_plus, V, V_pad_minus])
    assert([V_padded[n]<V_padded[n+1] for n, v in enumerate(V_padded[0:-2])])
    mu = 0 # Fermi energy
    conv = np.diff(fermi_dirac(V_padded, mu, T))
    fit = fano(V_padded, e0, w, q, a, b, c)
    # this is a list with edge effects - we only want the middle of this list
    # function over the space, including thermal broadening, w/no edge effects
    return np.convolve(fit, conv, mode="same")[pad_idcs:-pad_idcs]

def residual(data, fit):
    ld = len(data)
    r = [(data[i]-fit[i]) for i in range(ld)]
    # pdb.set_trace()
    return r, np.sqrt(sum([a*a for a in r]))/ld #normalized by # of data points

def fit_data(bias: np.array, dIdV: np.array, marker1: float, marker2: float):
    """
    Fits the data in dIdV(bias) to a fano lineshape using linear regression
    with reasonable range for the fit bounds.

    Parameters
    __________
    bias: array (in mV)
    dIdV: array
    marker1: float (in mV)
    marker2: float (in mV)

    Returns
    _______

    popt:
    pcov:
    sb, fit_dIdV
    """
    # data for which we are fitting the Fano function
    smallbias = [(n,b) for (n,b) in enumerate(bias) if b>=marker1 and b<=marker2]
    nsb, sb = np.array(smallbias).T
    fit_dIdV = np.array(dIdV)[[int(n) for n in nsb]]

    # initial guess for e0, w, q, a, b, c,
    b0 = (fit_dIdV[-1]-fit_dIdV[0])/(sb[-1]-sb[0])
    e00 = sb[np.argmin(scipy.signal.detrend(fit_dIdV))]
    p0 = [e00, 4, 1, 1, b0, np.mean(fit_dIdV)]

    # bounds for e0, w, q, a, b, c
    bounds = np.array([[min(sb),max(sb)],                   # e0
                        [0,50],                             # w
                        [-np.inf,np.inf],                   # q
                        [0, max(fit_dIdV)],                 # a
                        [-np.inf,np.inf],                   # b
                        [-np.inf,np.inf]]).T                # c
    # fix temperature while fitting other Fano parameters
    def fix_T(T):
        return lambda V, e0, w, q, a, b, c: fano_t_broad(T, V, e0, w, q, a, b, c)

    try:
        # popt, pcov = optimize.curve_fit(fix_T(T), sb, fit_dIdV, p0=p0, bounds=bounds)
        popt, pcov = optimize.curve_fit(fano, sb, fit_dIdV, p0=p0, bounds=bounds)
        return popt, pcov, sb, fit_dIdV

    except RuntimeError as e:
        print(e)
        n = np.ones((len(p0)))*np.nan
        return n, n, sb, fit_dIdV

def lorentz(e, e0, gamma):
    return 1/np.pi*(gamma/2)/((x-x0)**2+(gamma/2)**2)

def fit_data_fixed_vals(bias, dIdV, marker1, marker2, fixed_vals):
    """
    Parameters
    __________

    Returns
    _______
    """
    # data for which we are fitting the Fano function
    smallbias = [(n,b) for (n,b) in enumerate(bias) if b>=marker1 and b<=marker2]
    nsb, sb = np.array(smallbias).T
    fit_dIdV = np.array(dIdV)[[int(n) for n in nsb]]

    # initial guess for e0, w, q, a, b, c,
    b0 = (fit_dIdV[-1]-fit_dIdV[0])/(sb[-1]-sb[0])
    e00 = sb[np.argmin(scipy.signal.detrend(fit_dIdV))]

    p0 = [e00, 4, 0, 1, b0, np.mean(fit_dIdV)]
    hold = [0 if np.isnan(fixed_vals[n]) else 1 for n in range(len(p0)) ]

    p1 = [p if np.isnan(fixed_vals[n]) else fixed_vals[n] for n,p in enumerate(p0) ]
    p0 = [p for n,p in enumerate(p0) if np.isnan(fixed_vals[n]) ]

    # bounds for e0, w, q, a, b, c
    bounds = np.array([ [min(sb),max(sb)],     #none                  # e0
                        [0,50],                             # w
                        [-np.inf,np.inf],                   # q
                        [0, max(fit_dIdV)],                 # a
                        [-np.inf,np.inf],                   # b
                        [-np.inf,np.inf]]).T                # c
    bounds = np.array([b for n,b in enumerate(bounds.T) if np.isnan(fixed_vals[n])]).T

    # https://stackoverflow.com/questions/31705327/scipy-optimize-curve-fit-setting-a-fixed-parameter
    # fix some while fitting other Fano parameters
    def wrapper(V, *args):
        wrapperName = 'fano(V,'
        for i in range(0,len(hold)):
            if hold[i]:
                wrapperName +=str(p1[i])
            else:
                if i%2==0:
                    wrapperName += 'args['+str(i-sum(hold))+']'
                else:
                    wrapperName+='args['+str(i-sum(hold))+']'
            if i<len(hold):
                wrapperName+=','
        wrapperName+=')'
        return eval(wrapperName)
    try:
        # popt, pcov = optimize.curve_fit(fix_T(T), sb, fit_dIdV, p0=p0, bounds=bounds)
        popt, pcov = optimize.curve_fit(wrapper, sb, fit_dIdV, p0=p0, bounds=bounds)
        for i in range(0,len(hold)):
            if hold[i]:
                popt = np.insert(popt, i, p1[i])
                pcov = np.insert(np.insert(pcov, i, 0, axis=1), i, 0, axis=0)
        return popt, pcov, sb, fit_dIdV, p0
    except (RuntimeError, IndexError) as e:
        print(e)
        n = np.ones((len(p0)))*np.nan
        return n, n, sb, fit_dIdV, p0

def fit_data_w_times_residual(bias, dIdV, marker1, marker2, fixed_vals, init_vals=None, scale=True):
    """
    Parameters
    __________

    Returns
    _______
    """
    assert(marker1<marker2)
    # data for which we are fitting the Fano function
    # print(locals().keys())
    # print(locals.get(locals().keys()[0]))
    # print(inspect.getargvalues(inspect.currentframe()))

    smallbias = [(n,b) for (n,b) in enumerate(bias) if b>=marker1 and b<=marker2]
    nsb, sb = np.array(smallbias).T
    fit_dIdV = np.array(dIdV)[[int(n) for n in nsb]]

    # initial guess for e0, w, q, a, b, c,
    b0 = (fit_dIdV[-1]-fit_dIdV[0])/(sb[-1]-sb[0])
    e00 = sb[np.argmin(scipy.signal.detrend(fit_dIdV))]

    if init_vals is None:
        p0 = [e00, 4, 0.0001, 0.5, b0, np.mean(fit_dIdV)]
    else:
        p0 = init_vals
    hold = [0 if np.isnan(fixed_vals[n]) else 1 for n in range(len(p0)) ]

    p1 = [p if np.isnan(fixed_vals[n]) else fixed_vals[n] for n,p in enumerate(p0) ]
    p0 = [p for n,p in enumerate(p0) if np.isnan(fixed_vals[n]) ]

    # bounds for e0, w, q, a, b, c
    bounds = np.array([ [min(sb),max(sb)],     #none                  # e0
                        [2, 30],                             # w
                        [-np.inf,np.inf],                   # q
                        # [-0.2,0.2],
                        [0, max(fit_dIdV)],                 # a
                        [-np.inf,np.inf],                   # b
                        [-np.inf,np.inf]]).T                # c
    bounds = np.array([b for n,b in enumerate(bounds.T) if np.isnan(fixed_vals[n])]).T
    # bounds = [(b[0], b[1]) for b in bounds.T]
    # https://stackoverflow.com/questions/31705327/scipy-optimize-curve-fit-setting-a-fixed-parameter
    def wrapper(V, *args):
        wrapperName = 'fano(V,'
        for i in range(0,len(hold)):
            if hold[i]:
                wrapperName +=str(p1[i])
            else:
                if i%2==0:
                    wrapperName += 'args['+str(i-sum(hold))+']'
                else:
                    wrapperName+='args['+str(i-sum(hold))+']'
            if i<len(hold):
                wrapperName+=','
        wrapperName+=')'
        return eval(wrapperName)

    def objective_function(p, bias, dIdV):
        #parameters are [e0, w, q, a, b, c] so p[1] is w
        of = residual(dIdV, wrapper(bias, *p))[1]*abs(p[2])
        return of #residual(data, fit)

    try:
        # x_scale = [8,2,1,1,1,200]
        lsq_args = [objective_function]
        args = ([b[1] for b in smallbias], fit_dIdV,)
        lsq_kwargs = {"x0": p0, "args": args, "bounds": bounds, "max_nfev":4000}
        if scale:
            x_scale = np.abs(p0)
            # x_cale will complain `x_scale must be 'jac' or array_like with positive numbers`
            x_scale = [b for n,b in enumerate(x_scale) if np.isnan(fixed_vals[n])]
            # print(x_scale, fixed_vals,p0)
            lsq_kwargs["x_scale"] = x_scale

        res = optimize.least_squares(*lsq_args, **lsq_kwargs)
        # res = optimize.least_squares(objective_function, x0=p0,
        #                             args=,
        #                             bounds=bounds, max_nfev=4000,#,
        #                            x_scale=x_scale)#, ftol=3e-16, xtol=3e-16, gtol=3e-16)
        popt = res.x
        pcov = res.jac

        for i in range(0,len(hold)):
            if hold[i]:
                popt = np.insert(popt, i, p1[i])
                pcov = np.insert(np.insert(pcov, i, 0, axis=1), i, 0, axis=0)
        return popt, pcov, sb, fit_dIdV, p0

    except (RuntimeError, IndexError, ValueError) as e:
        print(e)
        n = np.ones((len(p0)))*np.nan
        return n, n, sb, fit_dIdV, p0

def save_fano_fits(files: list, opts: list, covs: list, m1: list, m2: list, path: str, xs: list, ys: list, resid: list, radii: list, dists: list):
    """
    Parameters
    __________

    Returns
    _______
    """
    with open(path, "w") as f:
        f.write("file\tradius\tdist\te0\tw\tq\ta\tb\tc\tmarker1\tmarker2\tx(nm)\ty(nm)\tresid\n")
        for nf, file in enumerate(files):
            # f.write("%s\t" %(os.path.split(file)[-1]))
            f.write("%s\t" %(file))
            f.write("%lf\t" %(radii[nf]))
            f.write("%lf\t" %(dists[nf]))
            try:
                for e in opts[nf]:
                    f.write(str("%1.2lf" %(e))+"\t")
                f.write(str("%1.2lf\t" %(m1[nf])))
                f.write(str("%1.2lf\t" %(m2[nf])))
                f.write(str("%1.2lf\t" %(xs[nf])))
                f.write(str("%1.2lf\t" %(ys[nf])))
                f.write(str("%1.2lf\t" %(resid[nf])))
            except:
                print("could not fit")
            # for c in covs[nf]:
            #     f.write(str(c)+"\t")
            f.write("\n")
    f.close()
    # need to save the 6x6x(# of lines) 3D array of covariances
    # idx = pd.MultiIndex.from_product([range(s) for s in np.array(covs).shape],
    #                                   names=["line_idx","param1","param2"])
    # df = pd.DataFrame({'covs': np.array(covs).flatten()}, index=idx)['covs']
    # df = df.unstack(level='line_idx').swaplevel().sort_index()

    # alternatively could only save 6x(# of lines) 2D array of std deviations
    with open(path.split(".txt")[0] + "_covs.txt","w") as f:
        f.write("file\te0\tw\tq\ta\tb\tc\n")
        for nf, file in enumerate(files):
            f.write("%s\t" %(file))
            try:
                for e in covs[nf].diagonal():
                    f.write("%1.2lf\t" %(e))
            except:
                for n in range(6):
                    f.write("0.0\t")
            f.write("\n")
        f.close()
    # np.savetxt(path.split(".txt")[0] + "_covs.txt", pd.DataFrame(covs[0]).values, fmt="%1.2lf")


def plot_line(image, specs, center):
    """

    Parameters
    __________


    Returns
    _______
    """
    # some_file.py
    # import sys
    # from matplotlib.widgets import Slider, Button
    #
    # # insert at 1, 0 is the script path (or '' in REPL)
    # sys.path.insert(1, '/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/Kondo data analysis')
    # from read_vertfile import Spec
    # import numpy as np
    # import createc

    # from scipy.interpolate import interp1d

    # str(pathlib.Path().resolve())
    # image = createc.DAT_IMG("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-06/Createc2_210806.140146.dat")

    dIdVs = [s.dIdV for s in specs]
    xlocs = [s.XPos_nm for s in specs]
    ylocs = [s.YPos_nm for s in specs]

    dist = np.sqrt(np.sum([(xlocs[0]-xlocs[-1])**2, (ylocs[0]-ylocs[-1])**2]))*10.

    def subtract_plane(im, xPix, yPix):
        X1, X2 = np.mgrid[:xPix, :yPix]
        nxny = xPix*yPix
        X = np.hstack((np.reshape(X1, (nxny, 1)), np.reshape(X2, (nxny, 1))))
        X = np.hstack((np.ones((nxny, 1)), X))
        YY = np.reshape(im, (nxny,1))
        theta = np.dot(np.dot(np.linalg.pinv(np.dot(X.transpose(), X)), X. transpose()), YY)
        plane = np.reshape(np.dot(X, theta), (xPix, yPix))
        im -= plane

    x_nm = np.round(image.size[0]/10.)
    y_nm = np.round(image.size[1]/10.)

    im = image._crop_img(image.img_array_list[2][:][:])
    zconst = image.zPiezoConst # angstroms/V

    DAC_V_conv = interp1d([-524288, 524287],[-10,10])

    im = DAC_V_conv(im)*zconst
    image.meta['fblogiset']
    subtract_plane(im, im.shape[0], im.shape[1])
    im = im-np.min(im)
    ###
    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 20}

    matplotlib.rc('font', **font)
    fig, (ax1,ax2) = plt.subplots(ncols=2, figsize=(12,6))

    fwd = np.array(dIdVs)[:,0:int(len(dIdVs[0])/2)].T
    bwd = reversed(np.array(dIdVs)[:,int(len(dIdVs[0])/2+1):].T)

    mean = np.mean(np.array(list(zip(fwd,bwd))),axis=1)
    mean = [m/m[-1] for m in mean]

    # making aspect larger makes image skinnier
    im1 = ax1.matshow(mean,
                      extent=[0,dist/10.,-100,100],
                      aspect="auto",
                      cmap=plt.get_cmap("plasma"))

    # np.array(dIdVs)[:,0:int(len(dIdVs[0])/2)].shape
    ax1.set_xlabel("nm")
    ax1.xaxis.set_ticks_position("bottom")
    ax1.set_ylabel(r"Bias ($mV$)")

    # ax1.set_title("$dI/dV$ spectra")
    # ax1.hlines(-67,0,25,"r")
    # divider = make_axes_locatable(ax1)
    # divider.__dict__.keys()
    # divider.__dict__['_axes'].__dict__

    # cax = divider.append_axes('right', size='5%', pad=0.05)
    width = 0.02
    height = 0.38
    vertical_position = 0.32
    horizontal_position = 0.1
    # axColor = plt.axes([horizontal_position, vertical_position, width,height])
    # cbar = fig.colorbar(im1, orientation='vertical',shrink=0.25)#,fraction=0.046, pad=0.04,shrink=0.25)

    # cbar.__dict__['ax']
    # cbar.__dict__.keys()
    # cax.xaxis.set_ticks_position("default")

    # cbar.set_ticks([])
    ax1.set_box_aspect(1)
    # get_aspect(ax1)

    im2 = ax2.imshow(im, extent=[0,x_nm,y_nm,0],aspect="equal")
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(im2, cax=cax, orientation='vertical')#,fraction=0.046, pad=0.04)
    cbar.set_label("nm")
    # cbar.set_ticks([])
    ax2.plot(np.array(xlocs)-image.offset[0]/10+x_nm/2,np.array(ylocs)-image.offset[1]/10,"r")
    # pdb.set_trace()
    ax2.scatter(center[0], center[1])
    # ax2.set_title("Topography")
    ax2.set_xlabel("nm")

    # fig.suptitle("Line spectrum from Ag(111) step edge")
    plt.tight_layout()

    look_at_scale = False
    # if changing color scale necessary:
    if look_at_scale:
        c_max = "%1.2lf"
        nmnd = np.nanmin(dIdVs)
        nmxd = np.nanmax(dIdVs)
        ax_cmax  = plt.axes([0.25, 0.15, 0.65, 0.03])
        ax_cmin  = plt.axes([0.25, 0.05, 0.65, 0.03])
        # plt.savefig("/Users/akipnis/Desktop/step_edge_line_spectrum.pdf")
        s_cmax = Slider(ax_cmax, 'max', nmnd, nmxd, valfmt=c_max)
        s_cmin = Slider(ax_cmin, 'min', nmnd, nmxd, valfmt=c_max)

        def update(val, s=None):
            _cmax = s_cmax.val
            _cmin = s_cmin.val
            im1.set_clim(_cmin,_cmax)
            plt.draw()
        s_cmax.on_changed(update)
        s_cmin.on_changed(update)

    # plt.plot(mean[:,int(len(mean[0])/2)])
    # plt.show()
    plt.text(  # position text relative to Figure
        0.025, 0.95, '(a)',
        ha='left', va='top',
        transform=fig.transFigure
    )
    plt.text(  # position text relative to Figure
        0.485, 0.95, '(b)',
        ha='left', va='top',
        transform=fig.transFigure
    )

    plt.show()

srs = []
dpath = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/"
# srs.append(dpath + "Ag 2021-08-13 3p8 nm radius/3p8nm pm20mV line/Createc2_210813.105240.L0029.VERT")
srs.append(dpath + "Ag 2021-08-13 2p5 nm radius/2p5nm radius pm 20mV line spectrum/Createc2_210813.165557.L0029.VERT")

srs.append(dpath + "Ag 2021-08-13 3p8 nm radius/3p8nm pm20mV line/other spectrum/Createc2_210813.103843.VERT")

# srs.append(dpath + "Ag 2021-08-13 3p8 nm radius/pm 20mV spectrum on Co/Createc2_210813.104416.VERT")
# srs.append(dpath + "Ag 2021-07-29 corral built/Line spectrum across corral/Createc2_210729.180059.L0026.VERT")
srs.append(dpath + "Ag 2021-08-12 4p5 nm radius/pm 20mV line/Createc2_210812.171526.L0036.VERT")
# srs.append(dpath + "Ag 2021-08-13 2p5 nm radius/pm20 mV on Co/Createc2_210813.163836.VERT")

# large range spectra
lrs = []
lrs.append(dpath + "Ag 2021-08-13 2p5 nm radius/300mV to -200mV line/Createc2_210813.231403.L0031.VERT")
lrs.append(dpath + "Ag 2021-08-11/3p8 nm radius line spectra pm100mV/Createc2_210811.113827.L0029.VERT")

lrs.append(dpath + "Ag 2021-08-12 4p5 nm radius/Createc2_210812.163415.VERT")
# plot_lrs()
if __name__=="__main__":
    from application import Application
    import tkinter as tk
    # radii data file

    # small range spectra

    # lrs.append(dpath + "Ag 2021-08-09 2p5 nm radius/100mV spectrum on Co/Createc2_210809.162718.VERT")
    # lrs.append(dpath + "Ag 2021-08-12 4p5 nm radius/4p5 nm line spectrum pm100mV/Createc2_210812.154131.L0036.VERT")
    # lrs.append(dpath+ "Ag 2021-08-09 2p5 nm radius/pm 100mV line spectrum across corral/Createc2_210809.154356.L0025.VERT")
    # lrs.append(dpath + "Ag 2021-08-13 2p5 nm radius/pm 100 mV 2p5 nm radius line spectrum/Createc2_210813.173235.L0030.VERT")

    root = tk.Tk()
    app = Application(master=root)
    app.mainloop()
