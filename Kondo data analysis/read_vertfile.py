import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
import pdb
from scipy import optimize
import os
import matplotlib.transforms as transforms
from itertools import product
import tkinter as tk
from tkinter.filedialog import askopenfilenames, asksaveasfile
import pandas as pd
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

@dataclass #adds __init__() function, etc.
class Spectrum:
    file: str
    NPoints: int
    XPos_nm: float
    YPos_nm: float
    bias_mv: list
    current: list
    dIdV: list

def get_spec_data(fname: str):
    f = open(fname,"rb")
    d = f.readlines()
    a = [str(b).split('\\t') for b in d]
    for n, d in enumerate(a):
        if "\\r\\n'" in d:
            d.remove("\\r\\n'")
        if "b'DATA\\r\\n'" in d:
            data_idx = n

    #getting the "mystery line" with some data on it
    # only compatible with STM Version 4 or greater (?)
    n = [n for n,b in enumerate(a) if "b'DATA\\r\\n'" in b][0] + 1
    l = a[n][0].split(" ")
    l.remove("b'")
    l = [m for m in l if m != '']
    l[-1] = l[-1].strip("\\r\\n\'")
    l = [float(m) for m in l]

    NPoints = l[0]
    VertPosX_DAC = l[1]
    VertPosY_DAC = l[2]
    ChannelList = l[3]
    Channelcount = l[4]
    OutchannelList = l[5]
    OutChannelCount = l[6]
    XPos_nm = l[7]
    YPos_nm = l[8]
    data = a[data_idx+2:]
    bias_mv = [float(d[1]) for d in data][19:]
    current = [float(d[4]) for d in data][19:]
    dIdV = [float(d[5]) for d in data][19:]
    return bias_mv, dIdV, current, a, XPos_nm, YPos_nm

def plot_large_range_spectra():
    for n,f in enumerate(files):
        bias, dIdV, current, a = get_spec_data(f)
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
    covs = pd.read_csv(f.split(".txt")[0] + "_covs.txt", delimiter="\t", index_col=False)
    #total length of line in nm:
    len_nm = np.linalg.norm(np.array(d1.iloc[0][["x(nm)","y(nm)"]])-np.array(d1.iloc[-1][["x(nm)","y(nm)"]]))
    dists = [np.linalg.norm(np.array(d1.iloc[0][["x(nm)","y(nm)"]])-np.array(d1.iloc[i][["x(nm)","y(nm)"]])) for i in range(len(d1))]

    fs = d1["file"]
    specdata = np.flipud(np.rot90([get_spec_data(f)[1] for f in fs]))
    biasdata = list(reversed(get_spec_data(fs[0])[0]))
    ns = fs.str.extract(r"L([0-9]*).VERT").astype(float)
    d = np.array(d1[["e0","w","q","a","b","c"]])
    bias = np.arange(min(biasdata),max(biasdata),0.2)
    im = np.array([fano(bias,*n) for n in d]).T
    # plt.figure()
    fig, ((ax1, ax2, ax3), (ax5, ax6, ax7)) = plt.subplots(2,3,figsize=(16,6), sharex="col")
    cax = ax1.matshow(np.flipud(im)) #aspect=1/len(d),]extent=[0, len_nm, -100,100]
    # fig.colorbar(cax)
    im = ax1.get_images()
    extent =  im[0].get_extent()
    aspect = 1
    ax1.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

    nbias = bias.shape[0]
    no_labels = 7 # how many labels to see on axis 5
    step_bias = int(nbias / (no_labels - 1)) # step between consecutive labels
    bias_positions = np.arange(0,nbias,step_bias)
    bias_labels = list(reversed(bias[::step_bias]))
    ax1.set_yticks(bias_positions)
    ax1.set_yticklabels(["%1.2lf" %(b) for b in bias_labels])

    ax1.set_xticks(np.arange(0, len(ns), 3))
    ax1.set_xticklabels(["%d" %(int(d[0])) for d in np.array(ns)[0::3]])

    ax1.tick_params(axis="x", bottom=False, top=True, labelbottom=False, labeltop=True)
    # [tick.set_label1("a") for tick in ax1.xaxis.get_major_ticks()]
    # [tick.set_label1("a") for tick in ax1.xaxis.get_major_ticks()]

    ax2.errorbar(dists, d1["e0"], yerr=covs["e0"], fmt='o')
    ax2.set_xlim([0,len_nm])
    ax2.set_ylabel("mV")
    ax6.set_xlabel("Distance (nm)")

    ax3.errorbar(dists, d1["w"], yerr=covs["w"], fmt='o')
    ax7.set_xlabel("Distance (nm)")
    ax3.set_ylabel("mV")

    ax4 = ax3.twinx()
    ax4.errorbar(dists, d1["w"]/kb*1e-3, yerr=covs["w"]/kb*1e-3, fmt='o')
    # expected Kondo temperature for Co on Ag111 is 92K
    # width = 2k_BT_K

    ax5.matshow(specdata)
    nbias = len(biasdata)
    no_labels = 7 # how many labels to see on axis 5
    step_bias = int(nbias / (no_labels - 1)) # step between consecutive labels
    bias_positions = np.arange(0,nbias,step_bias)
    bias_labels = list(reversed(biasdata[::step_bias]))
    ax5.set_yticks(bias_positions)
    ax5.set_yticklabels(["%1.2lf" %(b) for b in bias_labels])

    # if all the markers are the same
    ax5.axhline(y=nbias-np.argmin(np.abs(biasdata-d1["marker1"].iloc[0])), color="r")
    ax5.axhline(y=nbias-np.argmin(np.abs(biasdata-d1["marker2"].iloc[0])), color="r")

    # if the markers are different
    for n, d in enumerate(d1["marker1"].values):
        ax5.scatter(ns[0][n]-min(ns[0]), nbias-np.argmin(np.abs(biasdata-d1["marker2"].iloc[n])),color="r")
        ax5.scatter(ns[0][n]-min(ns[0]), nbias-np.argmin(np.abs(biasdata-d)),color="r")

    im = ax5.get_images()
    extent =  im[0].get_extent()
    ax5.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
    ax5.tick_params(axis="x", bottom=True, top=False, labelbottom=True, labeltop=False)

    ax6.errorbar(dists, d1["q"], yerr=covs["q"], fmt='o')
    ax7.scatter(dists, d1["resid"])

    a1, a2 = fs[0].split("/")[5:][1:3]
    fig.suptitle("Fano fits to Kondo resonance on corral central atom: %s.%s" %(a1, a2))
    # add red lines across where the fit happens

    # ax2.set_ylim([0,20 ])
    # ax3.set_ylim([0,20])
    # ax4.set_ylim([0,20/kb*1e-3])
    # ax6.set_ylim([-0.5, 2])

    # set titles for the subplots
    ax1.set_title("Fano fit")
    ax2.set_title(r"$E_0$ (eV)")
    ax3.set_title(r"Width, $T_K$")
    ax5.set_title("Raw data")
    ax6.set_title("Q")
    ax7.set_title("Residuals (au)")

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

def fano_w_therm_broadening(V, e0, w, q, a, b, c):
    # padding so there are no edge effects from convolution included in fit
    dV = np.abs(V[0] - V[1])
    padmVs = 20 # amount of padding (in mV) on either side of the fit range
    pad_idcs = int(padmVs/dV) #padding in array indices based on dV
    V_pad_minus = [min(V)-dV*i for i in list(range(pad_idcs))]
    V_pad_plus = [max(V)+dV*i for i in list(reversed(list(range(pad_idcs))))]
    V_padded = np.concatenate([V_pad_plus, V, V_pad_minus])

    T = 5 # Kelvin
    mu = 0 # Fermi energy
    conv = np.diff(fermi_dirac(V_padded, mu, T))
    fit = fano(V_padded, e0, w, q, a, b, c)

    # this is a list with edge effects - we only want the middle of this list
    thermally_broadened = np.convolve(fit, conv, mode="same")

    # function over the space, including thermal broadening, w/no edge effects
    return thermally_broadened[pad_idcs:-pad_idcs]

def fit_fano(file: str, marker1: float = 0, marker2: float = 0,
             savefig: bool = True, showfig: bool = True) -> list:
    bias, dIdV, current, a, xpos, ypos = get_spec_data(file)

    # implement thermal broadening integral here?
    if marker1==0 and marker2==0: # then create / get new markers
        fig , ax = plt.subplots()
        line = plt.plot(bias, dIdV)
        lines = [line]
        markers = []
        for n, l in enumerate(lines):
            c = "black" #l[0].get_color()
            # initialize the markers at the following two locations
            for b in [-20, 20]:
                mini = np.argmin(np.abs(np.array(bias)-b))
                ix = l[0].get_xdata()[mini]
                iy = l[0].get_ydata()[mini]
                d = DraggableMarker(ax=ax, lines=l,
                                    initx=ix, inity=iy, color=c,
                                    marker=">", dir=dir)
                markers.append(d)
        plt.title(os.path.split(file)[-1])
        plt.show()

        marker_vals = sorted([m.marker[0].get_xydata()[0] for m in markers], key=lambda x: x[0])
        marker1 = marker_vals[0][0]
        marker2 = marker_vals[1][0]

    # data for which we are fitting the Fano function
    smallbias = [(n,b) for (n,b) in enumerate(bias) if b>=marker1 and b<=marker2]
    nsb, sb = np.array(smallbias).T
    fit_dIdV = np.array(dIdV)[[int(n) for n in nsb]]

    # initial guess for e0, w, q, a, b, c,
    p0 = [8, 6, 1, 1, 0, np.mean(fit_dIdV)]

    # bounds for e0, w, q, a, b, c
    bounds = np.array([[-20,16],
                        [0,20],
                        [-1,3],
                        [-np.inf,np.inf],
                        [-np.inf,np.inf],
                        [-np.inf,np.inf]]).T

    popt, pcov = optimize.curve_fit(fano_w_therm_broadening, sb, fit_dIdV, p0=p0, bounds=bounds)
    fig, ax = plt.subplots()

    # trying to plot the 'boundaries' of the fit using the Hessian diagonal
    # since there are 6 parameters w/+-sd, we have 2**6 = 64 parameter sets
    # t = np.prod(list(product(pcov.diagonal(), [1, -1])),axis=1).reshape(pcov.shape[0],2)
    # a = np.sum([list(product([popt[n]], t[n])) for n in range(len(popt))], axis=2)
    # y = list(product(*a))
    # sort parameter sets by residual to data, pick two with highest residual
    # l = list(sorted(y, key=lambda x: residual(fit_dIdV, fano(sb, *x))))[-2:]

    # plt.fill_between(sb, fano(sb, *l[0]), fano(sb, *l[1]))


    # trans = transforms.blended_transform_factory(
    #     ax.transAxes, ax.transAxes) #ax.transData
    trans = plt.gcf().transFigure

    plt.text(0.05,0.2,
        r"$\frac{dI}{dV}\sim a\frac{(q+\epsilon)^2}{1+\epsilon^2} + bV + c$"+"\n"
        r"$\epsilon = \frac{V-\epsilon_0}{w}$"+"\n"
        r"$w=\gamma/2$, $\gamma=$linewidth"+ "\n\n"
        'e0: %1.2lf$\pm$%1.2lf mV\n\n'
        'w: %1.2lf$\pm$%1.2lf mV\n\n'
        'q: %1.2lf$\pm$%1.2lf\n\n'
        'a: %1.2lf$\pm$%1.2lf\n\n'
        'b: %1.2lf$\pm$%1.2lf\n\n'
        'c: %1.2lf$\pm$%1.2lf'
         %(tuple(np.array(list(zip(popt,pcov.diagonal()))).flatten())),
        transform=trans)

    plt.plot(bias, dIdV)
    plt.plot(sb, fit_dIdV,"b-")
    f = fano(sb, *popt)
    plt.plot(sb, f,'r--')

    plt.subplots_adjust(left=0.4)
    plt.xlabel("Bias (mV)")
    t = os.path.split(file)[-1]
    plt.title(t)
    plt.legend(["best fit",r'fit data',"data"])

    if savefig:
        plt.savefig(file.split(".VERT")[0]+"%s_fano_fit.png" %(t))
    if showfig:
        plt.show()

    plt.close()
    resid = residual(fit_dIdV, f)
    return [popt, pcov, marker1, marker2, xpos, ypos, resid]

def residual(data, fit):
    return np.sqrt(sum([(data[i]-fit[i])**2 for i in range(len(data))]))

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
            for e in covs[nf].diagonal():
                f.write("%1.2lf\t" %(e))
            f.write("\n")
        f.close()
    # np.savetxt(path.split(".txt")[0] + "_covs.txt", pd.DataFrame(covs[0]).values, fmt="%1.2lf")

class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.master.title('Fano fitter')
        self.pack()
        self.create_widgets()

    def create_widgets(self):
        self.small = tk.Button(self)
        self.small["text"] = "Analyze small range spectra"
        self.small["command"] = self.analyze_small_range
        # self.small.grid(row = 0, column = 0, pady = 2, sticky="W")

        self.small.pack(side="top")

        self.large = tk.Button(self)
        self.large["text"] = "Analyze large range spectra"
        self.large["command"] = self.analyze_small_range
        # self.large.grid(row = 1, column = 2, pady = 2)

        self.large.pack(side="top")

        self.choose = tk.Button(self)
        self.choose["text"] = "Choose file(s) to analyze"
        self.choose["command"] = self.open_files
        self.choose.pack(side="top")

        self.var1 = tk.IntVar()
        # self.choose.grid(row = 0, column = 0, pady = 2)
        self.c1 = tk.Checkbutton(self.master, text='Keep same fit bounds',variable=self.var1, onvalue=1, offvalue=0)
        self.c1.select()
        self.c1.pack()


        self.var2 = tk.IntVar()
        # self.choose.grid(row = 0, column = 0, pady = 2)
        self.c2 = tk.Checkbutton(self.master, text='Save figures',variable=self.var2, onvalue=1, offvalue=0)
        self.c2.pack()

        self.var3 = tk.IntVar()
        # self.choose.grid(row = 0, column = 0, pady = 2)
        self.c3 = tk.Checkbutton(self.master, text='Create line spectra',variable=self.var3, onvalue=1, offvalue=0)
        self.c3.select()
        self.c3.pack()

        self.quit = tk.Button(self, text="QUIT", fg="red",
                              command=self.master.destroy)
        # self.quit.pack(side="bottom")

    def analyze_files(self, files):
        opts = []; covs = []; m1s = []; m2s = []; xs = []; ys = []; resids = []
        for n,f in enumerate(files):
            if n==0 or self.var1.get()==0:
                popt, pcov, marker1, marker2, x, y, r = fit_fano(f, savefig=self.var2.get())
            elif self.var1.get()==0:
                popt, pcov, marker1, marker2, x, y,r = fit_fano(f, savefig=self.var2.get(),
                                                                marker1=marker1, marker2=marker2)
            else: #use the same markers as previously
                try:
                    popt, pcov, marker1, marker2, x, y, r = fit_fano(f,
                                    marker1=marker1, marker2=marker2,
                                    savefig=self.var2.get(),
                                    showfig=False)
                except:
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
        opts, covs, m1, m2, xs, ys, resids = self.analyze_files(large_range_spectra)
        save_fano_fits(large_range_spectra, opts, covs,  m1, m2, xs, ys, resids)

    def analyze_small_range(self):
        opts, covs, m1, m2, xs, ys,resids = self.analyze_files(small_range_spectra)
        save_fano_fits(small_range_spectra, opts, covs,  m1, m2, xs, ys, resids)

    def open_files(self):
        filenames = askopenfilenames(parent=self.master)
        # self.master.destroy()
        self.update()
        opts, covs, m1, m2, xs, ys,resids = self.analyze_files(filenames)
        try:
            path = asksaveasfile(parent=root,
                                 # defaultextension=["txt", "*.txt"],
                                 initialfile="test.txt").name
            self.update()
            save_fano_fits(filenames, opts, covs, m1, m2, path, xs, ys, resids)
            if self.var3.get():
                plot_fano_fit_line(path)
        except (AttributeError, UnboundLocalError) as e:
            print("did not get path")
            print(e)

if __name__=="__main__":
    small_range_spectra = []
    small_range_spectra.append("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-12 4p5 nm radius/pm 20mV line/Createc2_210812.171526.L0036.VERT")
    small_range_spectra.append("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-13 3p8 nm radius/3p8nm pm20mV line/Createc2_210813.103843.VERT")
    small_range_spectra.append("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-13 2p5 nm radius/pm20 mV on Co/Createc2_210813.163836.VERT")

    ## these are the locations of the +-100mV spectra on the corrals
    large_range_spectra = []
    large_range_spectra.append("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-12 4p5 nm radius/Createc2_210812.163415.VERT")
    large_range_spectra.append("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-11/3p8 nm radius line spectra pm100mV/Createc2_210811.113827.L0029.VERT")
    large_range_spectra.append("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-09 2p5 nm radius/100mV spectrum on Co/Createc2_210809.162718.VERT")

    root = tk.Tk()
    app = Application(master=root)
    app.mainloop()
