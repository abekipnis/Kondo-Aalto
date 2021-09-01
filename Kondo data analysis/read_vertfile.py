import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy import optimize

files = []
files.append("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-12 4p5 nm radius/Createc2_210812.163415.VERT")
files.append("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-11/3p8 nm radius line spectra pm100mV/Createc2_210811.113827.L0029.VERT")
files.append("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-09 2p5 nm radius/100mV spectrum on Co/Createc2_210809.162718.VERT")

radii = [4.5, 3.8, 2.5]

def get_spec_data(fname):
    f = open(fname,"rb")
    d = f.readlines()
    a = [str(b).split('\\t') for b in d]
    for n, d in enumerate(a):
        if "\\r\\n'" in d:
            d.remove("\\r\\n'")
        if "b'DATA\\r\\n'" in d:
            data_idx = n

    data = a[data_idx+2:-1]
    bias_mv = [float(d[1]) for d in data][20:-1]
    current = [float(d[4]) for d in data][20:-1]
    dIdV = [float(d[5]) for d in data][20:-1]
    return bias_mv, dIdV, current, a

def plot_large_range_spectra():
    for n,f in enumerate(files):
        bias, dIdV, current, a = get_spec_data(f)
        # pdb.set_trace()
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

def fano_function(V, e0, w, q, a, b, c):
    def eps(V, e0, w):
        return (np.array(V)-e0)/w
    return a*((q + eps(V, e0, w))**2/(1+eps(V, e0, w)**2))+ b* V + c

def fit_fano(file: str, marker1: float = 0, marker2: float = 0) -> list:
    bias, dIdV, current, a = get_spec_data(file)

    if marker1==0 and marker2==0: # then create / get new markers
        fig , ax = plt.subplots()
        line = plt.plot(bias, dIdV)
        lines = [line]
        markers = []
        for n, l in enumerate(lines):
            c = l[0].get_color()
            for b in [1, -1]:
                ix = l[0].get_xdata()[b]
                iy = l[0].get_ydata()[b]
                d = DraggableMarker(ax=ax, lines=l,
                                    initx=ix, inity=iy, color=c, marker=">", dir=dir)
                markers.append(d)
        plt.show()

        marker_vals = sorted([m.marker[0].get_xydata()[0] for m in markers], key=lambda x: x[0])
        marker1 = marker_vals[0][0]
        marker2 = marker_vals[1][0]

    smallbias = [(n,b) for (n,b) in enumerate(bias) if b>=marker1 and b<=marker2]
    nsb, sb = np.array(smallbias).T
    fit_dIdV = np.array(dIdV)[[int(n) for n in nsb]]

    # initial guess for e0, w, q
    p0 = [8, 6, 1, 1, 0, np.mean(fit_dIdV)]
    bounds = [(0,16),(0,10),(-1,1), (), (), ()]
    popt, pcov = optimize.curve_fit(fano_function, sb, fit_dIdV, p0=p0)

    plt.plot(sb, fano_function(sb, *popt))
    plt.plot(sb, fit_dIdV)
    plt.show()
    return [popt, pcov, marker1, marker2]

def fit_fano_functions_to_Co_atoms_at_corral_centers(files: list) -> list:
    opts = []; covs = []
    popt, pcov, marker1, marker2 = fit_fano(files[0])
    opts.append(popt)
    covs.append(pcov)
    popt, pcov, *_ = fit_fano(files[1], marker1, marker2)
    opts.append(popt)
    covs.append(pcov)
    popt, pcov, *_ = fit_fano(files[2], marker1, marker2)
    opts.append(popt)
    covs.append(pcov)
    return [opts, covs]

opts, covs = fit_fano_functions_to_Co_atoms_at_corral_centers(files)
pdb.set_trace()

plot_large_range_spectra()
#divide by current at start

files = []
files.append("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-12 4p5 nm radius/pm 20mV line/Createc2_210812.171526.L0036.VERT")
files.append("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-13 3p8 nm radius/3p8nm pm20mV line/Createc2_210813.103843.VERT")
files.append("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-13 2p5 nm radius/pm20 mV on Co/Createc2_210813.163836.VERT")
plot_large_range_spectra()
