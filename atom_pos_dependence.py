import os, scipy, datetime
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import createc
from scipy import optimize
from collections import namedtuple
from AaltoAtoms import CircCorralData, Spec

basepath = "Y:/labdata/Createc/STMDATA/Ag(111)/2022-03 Co Kondo corrals"
fields =  ['datfile', 'height_percentile','vertfile','marker1','marker2', "dataclipmin", "dataclipmax", "fit_order", "edge_cutoff", "chan"]
corralspectrum = namedtuple('corralspectrum', fields, defaults=(None,)*len(fields))
# .dat file and corresponding .VERT file for central Co atom fitting

dir = r"Y:\labdata\Createc\STMDATA\Ag(111)\2022-03 Co Kondo corrals\04-17 Co Co\position dependence experiment"
#dir = r"Y:\labdata\Createc\STMDATA\Ag(111)\2022-03 Co Kondo corrals\04-15 Co-Co"
def show_atom_position_dependence(dir):
    files = os.listdir(dir)

    # clean the files of .jpegs, etc.
    files = [f for f in files if (f[-4:]==".dat" or f[-5:]==".VERT")]
    file_times = [os.path.getmtime(os.path.join(dir,f)) for f in files]
    file_times = [datetime.datetime.fromtimestamp(t) for t in file_times]

    def get_corresponding_vert(file, files, dir):
        time = os.path.getmtime(os.path.join(dir,file))
        time = datetime.datetime.fromtimestamp(time)

        timediffs = [abs(time-t)+datetime.timedelta(10*(t<time)) for t in file_times]
        min_idx = np.argsort(timediffs)[1]
        return files[min_idx]

    dat_vert_dict = {}
    for dat in files:
        if ".dat" in dat:
            vert = get_corresponding_vert(dat, files, dir)
            dat_vert_dict[dat] = vert

    dats = list(dat_vert_dict.keys())
    S = Spec(os.path.join(dir,dat_vert_dict[dats[-1]]))
    datpath = os.path.join(dir,dats[-1])
    C = CircCorralData(datpath, "", chan=0)

    C.subtract_plane()
    C.get_region_centroids(percentile=99)
    C.occupied = True
    C.get_corral_radius(1, savefig=False, showfig=False)

    image = createc.DAT_IMG(datpath)
    x_nm = np.round(image.size[0]/10.)
    y_nm = np.round(image.size[1]/10.)
    plt.figure(1)
    plt.imshow(C.im, extent=[0,x_nm,y_nm,0],aspect="equal")

    specs = []
    for d in dats:
        S = Spec(os.path.join(dir,dat_vert_dict[d]))
        r = S.fit_fano(marker1=-5, marker2=10, showfig=False)
        width = r[0][1]

        datpath = os.path.join(dir,d)
        C = CircCorralData(datpath, "", chan=0)
        C.subtract_plane()
        C.get_region_centroids(percentile=99, show=False)
        C.occupied = True
        C.get_corral_radius(1, savefig=False, showfig=False)

        xlocs = [S.XPos_nm]
        ylocs = [S.YPos_nm]

        x_pix = C.nm_to_pix(np.array(xlocs)-image.offset[0]/10.+x_nm/2)
        y_pix = C.nm_to_pix(np.array(ylocs)-image.offset[1]/10.)

        dist_to_center_pix = C.c - [x_pix[0], y_pix[0]]
        dist_to_center_nm = C.pix_to_nm(np.linalg.norm(dist_to_center_pix))

        specs.append([S, dist_to_center_nm, width])
        plt.figure(1)
        plt.scatter(np.array(xlocs)-image.offset[0]/10.+x_nm/2,np.array(ylocs)-image.offset[1]/10.)
    plt.show()
    biasmin = min(specs[0][0].bias_mv)
    biasmax = max(specs[0][0].bias_mv)
    specs = sorted(specs, key=lambda x: x[1])
    plt.plot(specs[0][0].dIdV)
    plt.imshow([s[0].dIdV/s[0].dIdV[np.argmin(abs(0-s[0].bias_mv))] for s in specs],
                extent=[biasmin, biasmax,0,max(np.array(specs)[:,1])],
                aspect="auto",)
    plt.show()
    return specs

specs = show_atom_position_dependence(dir)
biasmin = min(specs[0][0].bias_mv)
biasmax = max(specs[0][0].bias_mv)
plt.imshow([s[0].dIdV/s[0].dIdV[np.argmin(abs(0-s[0].bias_mv))] for s in specs],
            extent=[biasmin, biasmax,0,max(np.array(specs)[:,1])],
            aspect="auto",)
plt.scatter(np.array(specs)[:,1], np.array(specs)[:,2])
