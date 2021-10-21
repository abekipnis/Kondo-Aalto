import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib.pyplot import figure

from dataclasses import dataclass
from datetime import datetime

from scipy import optimize
import scipy.signal
from scipy.interpolate import interp1d

import os
import pdb
import traceback
import sys
import time
from itertools import product
import regex as re

import tkinter as tk
from tkinter import messagebox
from tkinter.filedialog import askopenfilenames, asksaveasfile, StringVar, OptionMenu

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
    base_time = time.perf_counter_ns()

    def __new__(mcs, name, bases, namespace):
        print(mcs, name, bases, namespace)
        namespace["__class_load_time__"] = time.perf_counter_ns() - LoadTimeMeta.base_time
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

        self.FBLogiset = float([b[0].split("=")[1] for b in self.a[0:n-2] if b[0].split("=")[0]=="b'FBLogIset"][0].strip("\\r\\n'")) # units: pA
        self.biasVoltage = float([b[0].split("=")[1] for b in self.a[0:n-2] if b[0].split("=")[0]=="b'BiasVoltage / BiasVolt.[mV]"][0].strip("\\r\\n'")) # mV
        self.bias_mv = np.array([float(d[1]) for d in self.data][19:])
        self.current = np.array([float(d[4]) for d in self.data][19:])
        self.dIdV = np.array([float(d[5]) for d in self.data][19:])

        self.bias_offset = interp1d(self.current, self.bias_mv)(0)

        #correcting the bias offset
        self.bias_mv -= self.bias_offset

    def fit_fano(self, marker1: float = 0, marker2: float = 0,
                 savefig: bool = True, showfig: bool = True, e0 = np.nan, w = np.nan) -> list:
        app.update()

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
            plt.show()

            marker_vals = sorted([m.marker[0].get_xydata()[0] for m in markers],
                                 key=lambda x: x[0])
            marker1 = marker_vals[0][0]
            marker2 = marker_vals[1][0]

        # e0, w, q, a, b, c
        fixed_vals = [e0, w, np.nan, np.nan, np.nan, np.nan]
        popt, pcov, sb, fit_dIdV = fit_data_fixed_vals(self.bias_mv, self.dIdV, marker1, marker2, fixed_vals)
        try:
            fig = figure(figsize=(8.5,6.6)) #width, height (inches)
            a1 = plt.subplot(2,1,1)
            trans = plt.gcf().transFigure
            c = tuple(np.array(list(zip(popt,pcov.diagonal()))).flatten())
            plt.text(0.05,0.2,
                r"$\frac{dI}{dV}\sim a\frac{(q+\epsilon)^2}{1+\epsilon^2} + bV + c$"+"\n"
                r"$\epsilon = \frac{V-\epsilon_0}{w}$"+"\n"
                r"$w=\gamma/2$, $\gamma=$linewidth"+ "\n\n"
                'e0: %1.2lf$\pm$%1.2lf mV\n\n'
                'w: %1.2lf$\pm$%1.2lf mV\n\n'
                'q: %1.2lf$\pm$%1.2lf\n\n'
                'a: %1.2lf$\pm$%1.2lf\n\n'
                'b: %1.2lf$\pm$%1.2lf\n\n'
                'c: %1.2lf$\pm$%1.2lf\n\n'
                 %(c),
                transform=trans)

            plt.text(0.05, 0.8, "fit range\nmin: %1.2lf mV\n"
                                "max: %1.2lf mV\n\n"
                                "Bias setpoint: %1.2lf mV\n"
                                "Current setpoint: %1.2lf pA\n"
                                "Bias offset: %1.3lf mV"
                                %(marker1,
                                marker2,
                                self.biasVoltage,
                                self.FBLogiset,
                                self.bias_offset),
                                transform=trans)
            a1.plot(self.bias_mv, self.dIdV)
            a1.plot(sb, fit_dIdV,"b-")
            a1.set_ylim(min(self.dIdV), max(self.dIdV))
            f = fano(sb, *popt)

            # f = fix_T(T)(sb, *popt)
            # plt.plot(bias, fix_T(T)(bias, *popt), "black")

            a1.plot(self.bias_mv, fano(self.bias_mv, *popt), 'r--')
            # plt.plot(sb, fano(sb, *popt1),'go', markersize=2)

            residY, residtot = residual(fit_dIdV, f) #NORMALIZE THIS BY # OF FIT POINTS

            plt.subplots_adjust(left=0.4)
            a1.set_xlabel("Bias (mV)")
            t1 = os.path.split(self.fname.split(dpath)[-1])
            t = t1[-1]
            a1.set_title("\n".join(list(t1)))
            a1.legend(["data",r'fit data',"model"])
            # if savefig:
            #     plt.savefig(file.split(".VERT")[0]+"%s_fano_fit.png" %(t))
            # if showfig:
            #     plt.show()
            a2 = plt.subplot(2, 2, 3)
            a2.plot(sb, residY)
            a2.set_xlabel("fit residuals")

            a3 = plt.subplot(2,2,4)
            a3.hist(residY)
            a3.set_xlabel("residual histogram")
            # plt.tight_layout()
            if savefig:
                plt.savefig(file.split(".VERT")[0]+"%s_fit_residual.png" %(t))
            if showfig:
                plt.show()

            plt.close()
        except Exception:
            plt.close()
            app.update()
            print(traceback.format_exc())
            messagebox.showerror("showerror","could not plot %s" %(self.fname))
            return [popt, pcov, marker1, marker2, self.XPos_nm, self.YPos_nm, 0]
        return [popt, pcov, marker1, marker2, self.XPos_nm, self.YPos_nm, residtot]


def plot_lrs():
    for n,f in enumerate(files):
        bias, dIdV, current, a, XPos_nm, YPos_nm, T_AUXADC6 = get_spec_data(f)
        p = np.array(dIdV)/current[0]
        plt.plot(bias, p, label="%1.2lf nm radius" %(radii[n]))
    plt.legend()
    plt.xlabel("Bias (mV)")
    plt.ylabel("d$\it{I}$/d$\it{V}$ (a.u.)")
    plt.title("dI/dV spectrum of Kondo resonance on Co at corral center")
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
    covs = pd.read_csv(f.split(".txt")[0] + "_covs.txt",
                                        delimiter="\t", index_col=False)
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
    # pdb.set_trace()
    ax1.set_yticks(bias_positions)
    ax1.set_yticklabels(["%d" %(int(b)) for b in bias_labels])

    # ax1.set_xticks(np.arange(0, len(ns), 3))
    pdb.set_trace()
    nlabels = 6
    step_nm = int(len(dists)/(nlabels-1))
    step_positions = np.arange(0, len(dists), step_nm)
    pos_labels = list(reversed(np.array(dists)[::step_nm]))
    pos_labels = ["%1.1lf" %(p) for p in pos_labels]
    ax1.set_xticks(step_positions)
    ax1.set_xticklabels(pos_labels)


    # TODO: Fix labeling of indexing (check 2.5nm pm 20mV line spectra in presentation)
    # TODO: account for if there are some data where fit doesn't work
    # TODO: make one of the x axis labels in nanometers

    # ax1.set_xticklabels(["%d" %(int(d[0])) for d in np.array(ns)[0::3]])

    ax1.tick_params(axis="x", bottom=False, top=True, labelbottom=False, labeltop=True)

    ax2.errorbar(dists, d1["e0"], yerr=covs["e0"], fmt='o')
    e0_10pct =  0.1*(max(d1["e0"])-min(d1["e0"]))
    ax2.set_xlim([0,len_nm])
    ax2.set_ylim([min(d1["e0"])-e0_10pct, max(d1["e0"])+e0_10pct])
    ax2.set_ylabel("mV")
    ax6.set_xlabel("Distance (nm)")

    ax3.errorbar(dists, d1["w"], yerr=covs["w"], fmt='o')
    w_10pct =  0.1*(max(d1["w"])-min(d1["w"]))
    ax3.set_ylim(min(d1["w"])-w_10pct, max(d1["w"]+w_10pct))
    ax7.set_xlabel("Distance (nm)")
    ax3.set_ylabel("mV")


    ax4 = ax3.twinx()
    ax4.errorbar(dists, d1["w"]/kb*1e-3, yerr=covs["w"]/kb*1e-3, fmt='o')
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

    ax6.errorbar(dists, d1["q"], yerr=covs["q"], fmt='o')
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
    plt.show()

def fermi_dirac(e, mu, T):
    dist = 1/(np.exp((np.array(e)-mu)/(kb*T*1e3)) + 1)
    return dist

def fano(V, e0, w, q, a, b, c):
    def eps(V, e0, w):
        return (np.array(V)-e0)/w
    fit_func = a*((q + eps(V, e0, w))**2/(1+eps(V, e0, w)**2))+ b*np.array(V) + c
    return fit_func

def fano_t_broad(T, V, e0, w, q, a, b, c):
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


def fit_data_fixed_vals(bias, dIdV, marker1, marker2, fixed_vals):
    # data for which we are fitting the Fano function
    smallbias = [(n,b) for (n,b) in enumerate(bias) if b>=marker1 and b<=marker2]
    nsb, sb = np.array(smallbias).T
    fit_dIdV = np.array(dIdV)[[int(n) for n in nsb]]

    # initial guess for e0, w, q, a, b, c,
    b0 = (fit_dIdV[-1]-fit_dIdV[0])/(sb[-1]-sb[0])
    e00 = sb[np.argmin(scipy.signal.detrend(fit_dIdV))]

    p0 = [e00, 4, 1, 1, b0, np.mean(fit_dIdV)]
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
        # print(wrapperName)
        return eval(wrapperName)
    # fix some while fitting other Fano parameters
    # def fix_e(e):
    #     return lambda V, w, q, a, b, c: fano(V, e, w, q, a, b, c)
        # return lambda **args: fano(args[0], [f if not np.isnan(f) else args[1:][n-sum(np.isnan(fixed_vals[:n]))] for n, f in enumerate(fixed_vals)])
    try:
        # popt, pcov = optimize.curve_fit(fix_T(T), sb, fit_dIdV, p0=p0, bounds=bounds)
        popt, pcov = optimize.curve_fit(wrapper, sb, fit_dIdV, p0=p0, bounds=bounds)
        for i in range(0,len(hold)):
            if hold[i]:
                popt = np.insert(popt, i, p1[i])
                pcov = np.insert(np.insert(pcov, i, 0, axis=1), i, 0, axis=0)
        return popt, pcov, sb, fit_dIdV

    except RuntimeError as e:
        print(e)
        n = np.ones((len(p0)))*np.nan
        return n, n, sb, fit_dIdV

def fit_data(bias, dIdV, marker1, marker2):
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

def residual(data, fit):
    ld = len(data)
    r = [(data[i]-fit[i]) for i in range(ld)]
    return r, np.sqrt(sum([a*a for a in r]))/ld #normalized by # of data points

def save_fano_fits(files: list, opts: list, covs: list, m1: list, m2: list, path: str, xs: list, ys: list, resid: list):
    with open(path, "w") as f:
        f.write("file\te0\tw\tq\ta\tb\tc\tmarker1\tmarker2\tx(nm)\ty(nm)\tresid\n")
        for nf, file in enumerate(files):
            # f.write("%s\t" %(os.path.split(file)[-1]))
            f.write("%s\t" %(file))
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

class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.master.title('Fano fitter')
        self.pack()
        self.d = pd.read_excel("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Small Kondo corral inventory - OneDrive copy.xlsx", header=2)

        self.create_widgets()

    def create_widgets(self):
        self.small = tk.Button(self)
        self.small["text"] = "Analyze small range spectra"
        self.small["command"] = self.analyze_small_range

        self.small.pack(side="top")

        self.large = tk.Button(self)
        self.large["text"] = "Analyze large range spectra"
        self.large["command"] = self.analyze_large_range
        # self.large.grid(row = 1, column = 2, pady = 2)

        self.large.pack(side="top")

        self.choose = tk.Button(self)
        self.choose["text"] = "Choose file(s) to analyze"
        self.choose["command"] = self.open_files
        self.choose.pack(side="top")

        self.keep_fit_bounds = tk.IntVar()
        self.c1 = tk.Checkbutton(self.master, text='Keep same fit bounds',
                            variable=self.keep_fit_bounds, onvalue=1, offvalue=0)
        # self.c1.select()
        self.c1.pack()

        self.save_figures = tk.IntVar()
        self.c2 = tk.Checkbutton(self.master, text='Save figures',
                            variable=self.save_figures, onvalue=1, offvalue=0)
        self.c2.pack()

        self.var3 = tk.IntVar()
        self.c3 = tk.Checkbutton(self.master, text='Create line spectra',
                            variable=self.var3, onvalue=1, offvalue=0)
        # self.c3.select()
        self.c3.pack()

        self.quit = tk.Button(self, text="QUIT", fg="red",
                              command=self.master.destroy)

        self.variable = StringVar(self.master)
        self.opts = self.d["First spectrum filename & path"].dropna().apply(lambda f: f.split("\\")[-1])
        self.opts = self.opts[self.d["Emax"].dropna().index]
        self.variable.set(self.opts.iloc[0]) #default value
        self.w = OptionMenu(self.master, self.variable, *self.opts)
        self.w.pack()

        self.analyze_line = tk.Button(self)
        self.analyze_line["text"] = "Analyze selected line spectrum"
        self.analyze_line["command"] = self.analyze_selected_line
        self.analyze_line.pack(side='top')

    def analyze_selected_line(self):
        v = self.opts==self.variable.get()
        i = self.opts[v].index
        id = self.d.iloc[i]
        Lmin = int(id["Lmin"].values[0])
        Lmax = int(id["Lmax"].values[0])
        Emax = id["Emax"].values[0]
        Emin = id["Emin"].values[0]

        E0f = id["E0_fixed"].values[0]

        Lrange = list(range(Lmin,Lmax+1)) #inclusive

        p = id["First spectrum filename & path"]

        # re.sub(match, sub, string)
        files = [re.sub("0001","%04d" %(L), p.values[0]) for L in Lrange]
        files = [f.replace("\\","*").replace("*",os.path.sep) for f in files]

        specs = [Spec(os.path.join(dpath,f)) for f in files]
        d = []
        for s in specs:
            d.append(s.fit_fano(savefig=self.save_figures.get(),
                                        marker1=Emin, marker2=Emax, e0=E0f, showfig=False))
            self.update()

        opts, covs, m1, m2, xs, ys, resids = np.array(d).T
        path = asksaveasfile(parent=root,
                             initialfile=files[0].split("/")[-1][0:-9]+".txt").name
        save_fano_fits(files, opts, covs, m1, m2, path, xs, ys, resids)
        self.update()

        plot_fano_fit_line(path)
        self.update()


    def analyze_files(self, files):
        opts = []; covs = []; m1s = []; m2s = []; xs = []; ys = []; resids = []
        for n,f in enumerate(files):
            # TODO: decide fit range based on "central" spectrum
            # TODO: plot associated .dat file and spectrum locations
            if n==0 or self.keep_fit_bounds.get()==0:
                s = Spec(f)
                print("Spectrum load time (ms): ", s.__class_load_time__/1e6)

                popt, pcov, marker1, marker2, x, y, r = s.fit_fano(savefig=self.save_figures.get())
            elif self.keep_fit_bounds.get()==0:
                s = Spec(f)
                print("Spectrum load time (ms): ", s.__class_load_time__/1e6)

                popt, pcov, marker1, marker2, x, y, r = s.fit_fano(savefig=self.save_figures.get(),
                                                                marker1=marker1, marker2=marker2)
            else: #use the same markers as previously
                try:
                    s = Spec(f)
                    popt, pcov, marker1, marker2, x, y, r = s.fit_fano(
                                    marker1=marker1, marker2=marker2,
                                    savefig=self.save_figures.get(),
                                    showfig=False)
                except:
                    print(traceback.format_exc())

                    print("could not fit to %s" %(os.path.split(f)[-1]))
            opts.append(popt)
            covs.append(pcov)
            m1s.append(marker1)
            m2s.append(marker2)
            xs.append(x)
            ys.append(y)
            resids.append(r)
        return [opts, covs, m1s, m2s, xs, ys, resids]

    def analyze_large_range(self):
        opts, covs, m1, m2, xs, ys, resids = self.analyze_files(lrs)
        if self.var3.get():
            path = asksaveasfile(parent=root,
                                 # defaultextension=["txt", "*.txt"],
                                 initialfile="test.txt").name
            self.update()
            save_fano_fits(lrs, opts, covs,  m1, m2, "", xs, ys, resids)

    def analyze_small_range(self):
        opts, covs, m1, m2, xs, ys,resids = self.analyze_files(srs)
        if self.var3.get():
            path = asksaveasfile(parent=root,
                                 # defaultextension=["txt", "*.txt"],
                                 initialfile="test.txt").name
            self.update()
            save_fano_fits(srs, opts, covs,  m1, m2, "", xs, ys, resids)

    def open_files(self):
        filenames = askopenfilenames(parent=self.master)
        self.update()
        opts, covs, m1, m2, xs, ys, resids = self.analyze_files(filenames)
        try:
            if self.var3.get():
                path = asksaveasfile(parent=root,
                                     # defaultextension=["txt", "*.txt"],
                                     initialfile=filenames[0].split("/")[-1][0:-9]+".txt").name
                self.update()
                save_fano_fits(filenames, opts, covs, m1, m2, path, xs, ys, resids)
                plot_fano_fit_line(path)
        except (AttributeError, UnboundLocalError) as e:
            print("did not get path")
            print(e)

if __name__=="__main__":
    # small range spectra
    srs = []
    dpath = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/"
    srs.append(dpath + "Ag 2021-08-13 3p8 nm radius/3p8nm pm20mV line/Createc2_210813.105240.L0029.VERT")
    srs.append(dpath + "Ag 2021-08-13 3p8 nm radius/pm 20mV spectrum on Co/Createc2_210813.104416.VERT")
    srs.append(dpath + "Ag 2021-07-29 corral built/Line spectrum across corral/Createc2_210729.180059.L0026.VERT")
    srs.append(dpath + "Ag 2021-08-13 2p5 nm radius/2p5nm radius pm 20mV line spectrum/Createc2_210813.165557.L0029.VERT")
    srs.append(dpath + "Ag 2021-08-12 4p5 nm radius/pm 20mV line/Createc2_210812.171526.L0036.VERT")
    srs.append(dpath + "Ag 2021-08-13 3p8 nm radius/3p8nm pm20mV line/Createc2_210813.103843.VERT")

    srs.append(dpath + "Ag 2021-08-13 2p5 nm radius/pm20 mV on Co/Createc2_210813.163836.VERT")

    # large range spectra
    lrs = []
    lrs.append(dpath + "Ag 2021-08-12 4p5 nm radius/4p5 nm line spectrum pm100mV/Createc2_210812.154131.L0036.VERT")
    lrs.append(dpath+ "Ag 2021-08-09 2p5 nm radius/pm 100mV line spectrum across corral/Createc2_210809.154356.L0025.VERT")
    lrs.append(dpath + "Ag 2021-08-12 4p5 nm radius/Createc2_210812.163415.VERT")
    lrs.append(dpath + "Ag 2021-08-11/3p8 nm radius line spectra pm100mV/Createc2_210811.113827.L0029.VERT")
    lrs.append(dpath + "Ag 2021-08-09 2p5 nm radius/100mV spectrum on Co/Createc2_210809.162718.VERT")
    lrs.append(dpath + "Ag 2021-08-13 2p5 nm radius/pm 100 mV 2p5 nm radius line spectrum/Createc2_210813.173235.L0030.VERT")
    lrs.append(dpath + "Ag 2021-08-13 2p5 nm radius/300mV to -200mV line/Createc2_210813.231403.L0031.VERT")

    root = tk.Tk()
    app = Application(master=root)
    app.mainloop()
