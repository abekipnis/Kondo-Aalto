import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy import optimize
import os
import matplotlib.transforms as transforms
from itertools import product
import tkinter as tk
from tkinter.filedialog import askopenfilenames
from fano_fit_line_spectrum import plot_fano_fit_line

radii = [4.5, 3.8, 2.5]

def get_spec_data(fname: str):
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

def plot_fano_fit_line(f):
    d1 = pd.read_csv(f, delimiter="\t", index_col=False)
    fs = d1["file"]
    ns = fs.str.extract(r"L([0-9]*).VERT").astype(float)
    d = np.array(d1[d1.columns[1:]])
    im = [fano(np.arange(-100,100,0.2),*n) for n in d]
    im = np.array(im).T
    plt.imshow(im, aspect=1/len(d), extent=[min(ns), max(ns), -100,100])
    plt.xlabel("line index")
    plt.ylabel("bias (mV)")
    plt.colorbar()
    plt.title("Fano fits to Kondo resonance on corral central atom")
    plt.show()

def fano(V, e0, w, q, a, b, c):
    def eps(V, e0, w):
        return (np.array(V)-e0)/w
    return a*((q + eps(V, e0, w))**2/(1+eps(V, e0, w)**2))+ b* V + c

def fit_fano(file: str, marker1: float = 0, marker2: float = 0, savefig: bool = True) -> list:
    bias, dIdV, current, a = get_spec_data(file)

    # implement thermal broadening integral here?
    if marker1==0 and marker2==0: # then create / get new markers
        fig , ax = plt.subplots()
        line = plt.plot(bias, dIdV)
        lines = [line]
        markers = []
        for n, l in enumerate(lines):
            c = l[0].get_color()
            for b in [-20, 20]:
                mini = np.argmin(np.abs(np.array(bias)-b))
                ix = l[0].get_xdata()[mini]
                iy = l[0].get_ydata()[mini]
                d = DraggableMarker(ax=ax, lines=l,
                                    initx=ix, inity=iy, color=c, marker=">", dir=dir)
                markers.append(d)
        plt.title(os.path.split(file)[-1])
        plt.show()

        marker_vals = sorted([m.marker[0].get_xydata()[0] for m in markers], key=lambda x: x[0])
        marker1 = marker_vals[0][0]
        marker2 = marker_vals[1][0]

    smallbias = [(n,b) for (n,b) in enumerate(bias) if b>=marker1 and b<=marker2]
    nsb, sb = np.array(smallbias).T
    fit_dIdV = np.array(dIdV)[[int(n) for n in nsb]]

    # initial guess for e0, w, q, a, b, c,
    p0 = [8, 6, 1, 1, 0, np.mean(fit_dIdV)]
    bounds = np.array([[0,16],[0,10],[-1,1], [-np.inf,np.inf], [-np.inf,np.inf], [-np.inf,np.inf]]).T

    popt, pcov = optimize.curve_fit(fano, sb, fit_dIdV, p0=p0, bounds=bounds)
    fig, ax = plt.subplots()

    # since there are 6 parameters w/+-sd, we have 2**6 = 64 parameter sets
    t = np.prod(list(product(pcov.diagonal(), [1, -1])),axis=1).reshape(pcov.shape[0],2)
    a = np.sum([list(product([popt[n]], t[n])) for n in range(len(popt))], axis=2)
    y = list(product(*a))
    # sort parameter sets by residual to data, pick two with highest residual
    l = list(sorted(y, key=lambda x: residual(fit_dIdV, fano(sb, *x))))[-2:]
    plt.fill_between(sb, fano(sb, *l[0]), fano(sb, *l[1]))

    plt.plot(sb, fano(sb, *popt),'r--')

    trans = transforms.blended_transform_factory(
        ax.transAxes, ax.transAxes) #ax.transData
    plt.text(0.2,0.75,
        'e0: %1.2lf$\pm$%1.2lf\tw: %1.2lf$\pm$%1.2lf\tq: %1.2lf$\pm$%1.2lf'
        '\n a: %1.2lf$\pm$%1.2lf\tb: %1.2lf$\pm$%1.2lf\tc: %1.2lf$\pm$%1.2lf'
         %(tuple(np.array(list(zip(popt, pcov.diagonal()))).flatten())),
        transform=trans)

    # except Exception as e:
        # print(e)
    plt.plot(bias, dIdV)
    plt.plot(sb, fit_dIdV,"b-")
    t = os.path.split(file)[-1]
    plt.title(t)
    plt.legend(["best fit",r'$\sigma$',"data"])
    if savefig:
        plt.savefig("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/Kondo data analysis/%s_fano_fit.png" %(t))

    plt.show()
    return [popt, pcov, marker1, marker2]

def residual(data, fit):
    return np.sqrt(sum([(data[i]-fit[i])**2 for i in range(len(data))]))

def save_fano_fits(files: list, opts: list, covs: list, markers: list):
    log = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/Kondo data analysis/fano_fit_data.txt"
    with open(log, "w") as f:
        f.write("file\te0\tw\tq\ta\tb\tc\tmarker1\tmarker2\n")
        for nf, file in enumerate(files):
            f.write("%s\t" %(os.path.split(file)[-1]))
            try:
                for e in opts[nf]:
                    f.write(str("%1.2lf" %(e))+"\t")
            except:
                print("could not fit")
            # for c in covs[nf]:
            #     f.write(str(c)+"\t")
            f.write("\n")
    f.close()

class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
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
        self.c1.pack()


        self.var2 = tk.IntVar()
        # self.choose.grid(row = 0, column = 0, pady = 2)
        self.c2 = tk.Checkbutton(self.master, text='Save figures',variable=self.var2, onvalue=1, offvalue=0)
        self.c2.pack()

        self.quit = tk.Button(self, text="QUIT", fg="red",
                              command=self.master.destroy)
        # self.quit.pack(side="bottom")

    def analyze_files(self, files):
        opts = []; covs = []; m1s = []; m2s = []
        for n,f in enumerate(files):
            # define the markers based on the first (-10, 24)
            if n==0 or self.var1.get()==0:
                popt, pcov, marker1, marker2 = fit_fano(f, savefig=self.var2.get())
            else: #use the same markers as previously
                try:
                    popt, pcov, *_ = fit_fano(f,
                                    marker1=marker1, marker2=marker2, savefig=self.var2.get())
                except:
                    print("could not fit to %s" %(os.path.split(f)[-1]))
            opts.append(popt)
            covs.append(pcov)
            m1s.append(marker1)
            m2s.append(marker2)
        return [opts, covs, m1s, m2s]

    def analyze_large_range(self):
        opts, covs = self.analyze_files(large_range_spectra)
        save_fano_fits(large_range_spectra, opts, covs)

    def analyze_small_range(self):
        opts, covs = self.analyze_files(small_range_spectra)
        save_fano_fits(small_range_spectra, opts, covs)

    def open_files(self):
        filenames = askopenfilenames(parent=self.master)
        opts, covs = self.analyze_files(filenames)
        save_fano_fits(filenames, opts, covs)

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
