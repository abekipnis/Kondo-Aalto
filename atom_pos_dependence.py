import os, scipy, datetime
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import createc
from scipy import optimize
from collections import namedtuple
from AaltoAtoms import CircCorralData, Spec, clockwiseangle
import pickle
from sklearn import cluster

basepath = "Y:/labdata/Createc/STMDATA/Ag(111)/2022-03 Co Kondo corrals"
fields =  ['datfile', 'height_percentile','vertfile','marker1','marker2', "dataclipmin", "dataclipmax", "fit_order", "edge_cutoff", "chan"]
corralspectrum = namedtuple('corralspectrum', fields, defaults=(None,)*len(fields))
# .dat file and corresponding .VERT file for central Co atom fitting

#lots of data points, medium-ish corral, can't see eigenmodes
dir = r"Y:\labdata\Createc\STMDATA\Ag(111)\2022-03 Co Kondo corrals\04-17 Co Co\position dependence experiment"

# small experiment with only 22 data points
dir = r"Y:\labdata\Createc\STMDATA\Ag(111)\2022-03 Co Kondo corrals\04-08\7nm loose corral - center atom position"

# 233 data points ! good experiment, big corral
dir = r"Y:\labdata\Createc\STMDATA\Ag(111)\2022-03 Co Kondo corrals\04-15 Co-Co"

import pdb
def show_atom_position_dependence(dir, fit_fano=False):
    pickle_file = os.path.join(dir, "spec_dat_dict.pickle")
    file_exists = False
    if os.path.exists(pickle_file):
        file_exists = True
        print("loading pickled data")
        with open(pickle_file, "rb") as handle:
            C, specs = pickle.load(handle)
        print("loaded pickled data")
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

    datpath = os.path.join(dir,dats[0])

    C = CircCorralData(datpath, "", chan=0)
    C.subtract_plane()
    C.get_region_centroids(percentile=99)
    C.occupied = True
    C.get_corral_radius(1, savefig=False, showfig=False)

    image = createc.DAT_IMG(datpath)
    x_nm = np.round(image.size[0]/10.)
    y_nm = np.round(image.size[1]/10.)

    plt.figure(1)
    fig, (ax1, ax2) = plt.subplots(1,2)
    ax2.imshow(C.im, extent=[0,x_nm,y_nm,0],aspect="equal")

    if not os.path.exists(pickle_file):

        specs = []
        for d in dats:
            S = Spec(os.path.join(dir, dat_vert_dict[d]))
            #S.clip_data(-25,50)
            #S.remove_background(3)
            width = False
            if fit_fano:
                r = S.fit_fano(marker1=-5, marker2=10, showfig=False, savefig=False, type_fit="wtimes")
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

            angle_from_center = clockwiseangle([x_pix[0], y_pix[0]], C.c, [1,0])

            if dist_to_center_nm < C.pix_to_nm(C.r):
                specs.append([S, dist_to_center_nm, width, angle_from_center, C])
        #        plt.figure(1)
                xpos = np.array(xlocs)-image.offset[0]/10.+x_nm/2
                ypos = np.array(ylocs)-image.offset[1]/10.
                ax2.scatter(xpos, ypos, color='red')
        #plt.show()
    biasmin = min(specs[0][0].bias_mv)
    biasmax = max(specs[0][0].bias_mv)
    specs = sorted(specs, key=lambda x: x[1])
    #plt.plot(specs[0][0].dIdV)
    norm = [s[0].dIdV[np.argmin(abs(7.5-s[0].bias_mv))] for s in specs]
    d = [list(reversed(s[0].dIdV/norm[n])) for n, s in enumerate(specs)]

    # ax1.matshow(d,
    #             extent=[biasmin, biasmax,0,max(np.array(specs)[:,1])],
    #             aspect="auto",)

    X = specs[0][0].bias_mv
    Y = [s[1] for s in specs]
    Z = [s[0].dIdV for s in specs]

    ax1.pcolormesh(X ,Y, Z)

    ax1.xaxis.set_ticks_position("bottom")
    ax2.yaxis.set_ticks_position("right")
    ax2.set_ylabel('nm')
    ax1.set_xlabel("mV")
    ax1.set_ylabel("distance from center (nm)")
    plt.tight_layout()
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\7nm_position dependence.pdf")
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\7nm_position dependence.png")
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\7nm_position dependence.svg")

    plt.show()
    if not file_exists:
        with open(pickle_file, "wb") as handle:
            pickle.dump([C,specs], handle, protocol=pickle.HIGHEST_PROTOCOL)

    return specs

specs = show_atom_position_dependence(dir)
len(specs)
from sklearn.cluster import MeanShift
allspecs = np.array([s[0].dIdV for s in specs])


kmeans = MeanShift(n_clusters=2).fit(allspecs)
kmeans.labels_
cat = np.argwhere(kmeans.labels_==1).flatten()
cat
plt.pcolormesh(X,np.array(Y)[cat], allspecs[cat].reshape(len(cat),-1))


sfor c in specs:
    X, Y = np.array(c[-1].centroids).T
    m = np.argmin(np.linalg.norm(c[-1].centroids - c[-1].c, axis=-1))

    X = np.delete(X, m)
    Y = np.delete(Y, m)

    X -= X[-1]
    Y -= Y[-1]
    plt.scatter(X, Y)
#    plt.scatter(*c[-1].c)

plt.hlines(np.mean([s[0].current[0] for s in specs]), 0, 300)
plt.scatter([s[0].current[100] for s in specs], [s[0].current[-1] for s in specs])


def narrow_search(specs):
    mean_current0 = np.mean([s[0].current[0] for s in specs])
#if s[0].current[0] > mean_current0 and
    specs_n = [s for s in specs if np.isclose(1,s[3],0.1)]
    X = specs_n[0][0].bias_mv
    Y =  [s[1] for s in specs_n]
    Z = [s[0].dIdV/s[0].dIdV[np.argmin(abs(-30-s[0].bias_mv))] for s in specs_n]
    Z = (np.array(Z).T-np.min(Z, axis=-1))/(np.max(Z, axis=-1)- np.min(Z, axis=-1))


    Zx = np.array([np.polyval(np.polyfit(X,np.array(z),2), X) for z in np.array(Z).T])
    Z -=Zx.T
    plt.figure(figsize=(9,9))
    plt.pcolormesh(X, Y, Z.T)

narrow_search(specs)

normed_s = [(s[0].dIdV-min(s[0].dIdV))/(max(s[0].dIdV)-min(s[0].dIdV)) for s in specs]
plt.hist([s[0].dIdV[np.argmin(abs(5-s[0].bias_mv))] for s in specs], bins=50)
subset=  [s for s in specs  if np.isclose(2,s[-1],0.5)]

plt.plot(subset[7][0].bias_mv, subset[7][0].dIdV)
plt.plot(subset[8][0].dIdV)


plt.scatter([s[1] for s in specs],[s[-1] for s in specs])
plt.plot(np.array(specs)[:,1],np.array(specs)[:,2])


biasmin = min(specs[0][0].bias_mv)
biasmax = max(specs[0][0].bias_mv)

plt.imshow([list(reversed(s[0].dIdV/s[0].dIdV[0])) for s in specs],
            extent=[biasmin, biasmax,0,max(np.array(specs)[:,1])],
            aspect="auto",)
# plt.scatter(np.array(specs)[:,1], np.array(specs)[:,2])
