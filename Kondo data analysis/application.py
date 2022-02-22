import tkinter as tk
from tkinter import messagebox
from tkinter.filedialog import askopenfilenames, asksaveasfile, StringVar, OptionMenu
import pandas as pd
import numpy as np
import pathlib
import re, os
from read_vertfile import Spec, save_fano_fits
import createc
rdf = """/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/circle fit plots/circle_fits.txt"""
dpath = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/"

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


        radii_data = np.loadtxt(rdf, delimiter=",",
                                dtype=object,
                                converters={1: float})
        cff = [str(pathlib.PureWindowsPath(r[0])) for r in radii_data]

        self.variable = StringVar(self.master)
        self.opts = self.d["First spectrum filename & path"].dropna().apply(lambda f: f.split("\\")[-1])
        self.opts = self.opts[self.d["Emax"].dropna().index]
        self.variable.set(self.opts.iloc[0]) #default value
        idcs = [cff.index(v) for v in self.d.iloc[self.opts.index]["Scan file name & path"].values]

        # removes duplicates while preserving order
        vals = list(dict.fromkeys(self.opts.values))

        # reading radii from radii_data array from file
        radii = radii_data[idcs,1]

        #removes duplicates while preserving order
        radii = list(dict.fromkeys(radii))

        print(radii)

        self.w = OptionMenu(self.master, self.variable, *vals)
        self.w.pack()

        self.analyze_line = tk.Button(self)
        self.analyze_line["text"] = "Analyze selected line spectrum"

        # need lambda here otherwise function will run at app start
        self.analyze_line["command"] = lambda: self.analyze_selected_line(self.variable.get())
        self.analyze_line.pack(side='top')

        self.analyze_all = tk.Button(self)
        self.analyze_all["text"] = "Analyze all line spectra"
        self.analyze_all["command"] = self.analyze_all_spectra
        self.analyze_all.pack(side='top')

        self.reload_excel = tk.Button(self)
        self.reload_excel["text"] = "Reload Excel"
        self.reload_excel["command"] = self.reload_excel
        self.reload_excel.pack(side='top')

    def reload_excel(self):
        self.d = pd.read_excel("/Users/akipnis/Desktop/"
                                "Aalto Atomic Scale Physics/"
                                "Small Kondo corral inventory - OneDrive copy.xlsx",
                                header=2)

    def analyze_all_spectra(self):
        for o in list(set(self.opts.values)):
            print("ANALYZING %s" %(o))
            self.analyze_selected_line(o)
        return

    def analyze_selected_line(self,u):
        """
        Parameters
        __________

        Returns
        _______
        """
        v = self.opts==u
        i = self.opts[v].index
        id = self.d.iloc[i]


        # sometimes there are two lines in the excel sheet
        # for example, if the fit bounds need to change between spectra
        # in that case, save the files with an addendum saying which
        # n out of the N from this line spectra it is
        for n in range(len(id["First spectrum filename & path"])):
            Lmin = int(id["Lmin"].values[n])
            Lmax = int(id["Lmax"].values[n])
            Emax = id["Emax"].values[n]
            Emin = id["Emin"].values[n]

            E0f = id["E0_fixed"].values[n]

            Lrange = list(range(Lmin,Lmax+1)) #inclusive

            p = id["First spectrum filename & path"]

            # re.sub(match, sub, string)
            files = [re.sub("0001","%04d" %(L), p.values[0]) for L in Lrange]
            files = [f.replace("\\","*").replace("*",os.path.sep) for f in files]

            radii_data = np.loadtxt(rdf, delimiter=",",
                                    dtype=object,
                                    converters={1: float})

            specs = [Spec(os.path.join(dpath,f)) for f in files]

            #cff = correctly formatted files
            cff = [str(pathlib.PureWindowsPath(r[0])) for r in radii_data]
            have_radius = id["Scan file name & path"].values[n] in cff
            if have_radius:
                idx = cff.index(id["Scan file name & path"].values[n])
                datfile = radii_data[idx][0]
                datimg = createc.DAT_IMG(os.path.join(dpath,datfile))
                radius = radii_data[idx][1]
                center_x = float(radii_data[idx][2])
                center_y = float(radii_data[idx][3])
                # plot_line(datimg, specs, [center_x,center_y])
            else:
                print('could not find corral radius for this spectrum')
                pdb.set_trace()
                radius = np.nan
            d_w = [] #with width*residual
            d_d = [] #default fitting
            radii = np.repeat(radius, len(specs))

            xlocs = [s.XPos_nm for s in specs]
            ylocs = [s.YPos_nm for s in specs]


            xds = np.array(xlocs)-datimg.offset[0]/10.+np.round(datimg.size[0]/10.)/2 - center_x
            yds = np.array(ylocs) - datimg.offset[1]/10 - center_y
            dists = list(map(np.sqrt, (xds**2+ yds**2)))
            # pdb.set_trace()
            for ns, s in enumerate(specs):
                old = s.fit_fano(savefig=self.save_figures.get(),
                                            marker1=Emin, marker2=Emax,
                                            e0_fixed_val=np.nan, showfig=False,
                                            actual_radius=radius,
                                            dist_to_Co=dists[ns])
                d_d.append(old)

                # use the best guess on the old fit for the new fit
                # rescaling the values
                new = s.fit_fano(savefig=self.save_figures.get(),
                                            marker1=Emin, marker2=Emax,
                                            e0_fixed_val=np.nan,
                                            showfig=False, type_fit="wtimes",
                                            init_vals=np.array(old[0]).T,
                                            actual_radius=radius,
                                            dist_to_Co=dists[ns])
                d_w.append(new)
                self.update()

            for d in [[d_w,"width"], [d_d,"default"]]:
                opts, covs, m1, m2, xs, ys, resids = np.array(d[0],dtype=object).T
                # pdb.set_trace()
                # path = asksaveasfile(parent=root,
                #                      initialfile=files[0].split("/")[-1][0:-9]+"_kondo_fit_param_%s.txt" %(d[1])).name
                bpath = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/spatial extent Kondo plots/width comparison/"
                p = files[0].split("/")[-1][0:-9]+"_%d_kondo_fit_param_%s.txt" %(n,d[1])
                path = bpath+p
                # pdb.set_trace()

                save_fano_fits(files, opts, covs, m1, m2, path, xs, ys, resids, radii, dists)
                self.update()

                # plot_fano_fit_line(path)
                self.update()

    def analyze_files(self, files):
        """
        Parameters
        __________

        Returns
        _______
        """
        opts = []; covs = []; m1s = []; m2s = []; xs = []; ys = []; resids = []
        for n,f in enumerate(files):
            # TODO: decide fit range based on "central" spectrum
            # TODO: plot associated .dat file and spectrum locations
            if n==0 or self.keep_fit_bounds.get()==0:
                s = Spec(f)
                # pdb.set_trace()
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
                # path = asksaveasfile(parent=root,
                #                      # defaultextension=["txt", "*.txt"],
                #                      initialfile=filenames[0].split("/")[-1][0:-9]+".txt").name
                self.update()
                save_fano_fits(filenames, opts, covs, m1, m2, path, xs, ys, resids)
                plot_fano_fit_line(path)
        except (AttributeError, UnboundLocalError) as e:
            print("did not get path")
            print(e)
