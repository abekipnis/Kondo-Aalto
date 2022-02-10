"""
This code is for plotting the spatial decay of the Fano resonance fit parameters
as a function of distance from the central Co atom in Co/Ag111 quantum corrals.
"""
import matplotlib.pyplot as plt
import numpy as np
import pdb
import os
import re
import createc
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Rectangle
import matplotlib
import pandas as pd

def subtract_plane(im):
    xPix, yPix = im.shape
    X1, X2 = np.mgrid[:xPix, :yPix]
    nxny = xPix*yPix
    X = np.hstack((np.reshape(X1, (nxny, 1)), np.reshape(X2, (nxny, 1))))
    X = np.hstack((np.ones((nxny, 1)), X))
    YY = np.reshape(im, (nxny,1))
    theta = np.dot(np.dot(np.linalg.pinv(np.dot(X.transpose(), X)), X. transpose()), YY)
    plane = np.reshape(np.dot(X, theta), (xPix, yPix))
    im -= plane
    return im

path = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/grid analysis"

## TODO:
"""
-- make the figure ready for publication by:
---- adding either scale bars or axes
---- making sure that the locations are correctly indexed
---- making a colorbar for the height
----
"""
from skimage.draw import line
from itertools import repeat
import scipy.optimize

def Lorentzian(x, amp1, cen1, wid1):
    cen1 = 0
    return (amp1*wid1**2/((x-cen1)**2+wid1**2))

def show_radial_decay(im, Xs, Ys, title,disc):
    """
    Takes array of parameters from Fano fit,
    x and y values in nm of locations of pixels
    and plots the spatial decay from the central pixel.

    Try to recreate figure 3 from Knorr2002, spatial decay of q)

    Parameters
    __________
    im: Numpy array
    Xs:
    Ys: y values in nm of locations of pixels in im
    title: string, which parameter from the Fano fit this is
    disc: string, differentiating between fits

    Returns
    _______
    None
    """
    # disc = discriminator
    # to differentiate between files where the fit parameters have changed

    d_r = []
    ctr = np.array(list(map(int,np.array(im.shape)/2)))
    ctr = [Xs[ctr[0]][ctr[1]], Ys[ctr[0]][ctr[1]]]
    for nx, i in enumerate(im):
        for ny, j in enumerate(i):

            d = np.sqrt((Xs[ny][nx]-ctr[0])**2 + (Ys[ny][nx]-ctr[1])**2)
            d_r.append([d, j])
    # plt.close()
    x, y = np.array(d_r).T
    plt.figure(title)
    plt.scatter(x,y);
    if title=="a" or title=="w":
        x = x[~np.isnan(y)]
        y = y[~np.isnan(y)]
        if title=="a":
            p0 = [400, 0, 1]
        else:
            p0 = [12.5, 0, 1]
        popt, pcov = scipy.optimize.curve_fit(Lorentzian, x, y, p0=p0,maxfev=4000)
        plt.plot(x, Lorentzian(x, *popt), color="black")
    if title=="q":
        dir = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/spatial extent Kondo plots/width comparison"
        files = os.listdir(dir)

        files = [f for f in files if f[-9:] =="width.txt" or f[-11:]=="default.txt"]
        str(os.path.join(dir,files[0]))

        d = pd.DataFrame()

        # plt.figure(figsize=(9,6))
        for f in files:
            data = pd.read_csv(os.path.join(dir,f), delimiter="\t", index_col=False)
            # data = data[data["radius"]<2.5]
            qcutoff = 1.5
            dcutoff = 0.5
            # cut_data = data[np.abs(data["q"])<qcutoff]
            # cut_data = cut_data[cut_data["dist"]<dcutoff]
            plt.scatter(data["dist"], data["w"], c=data["radius"],cmap="summer")
            plt.clim(2.4,4.8)
        plt.ylim(-0.2, 1.2)
    plt.title(title)

    plt.savefig(os.path.join(path, "%s_%s_decay.pdf" %(disc, title)))

def plot_grid_fit_params(files, xmin, xmax, ymin, ymax):
    topo_file = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-13 2p5 nm radius/grid/Createc2_210814.214635.specgrid.dat"
    img = createc.DAT_IMG(topo_file)

    limit_files = [f for f in files if "limits" in f]
    dat_files = [f for f in files if "limits" not in f]
    pdb.set_trace()
    assert(len(limit_files)>0 and len(dat_files)>0 and len(dat_files)==len(limit_files))

    # img.meta['specgriddx']
    # img.meta['specgriddy']
    # img.meta['specgridnx']
    # img.meta['specgridny']
    # img.meta['specgrid_centerx']
    # img.meta['specgrid_centery']

    X,Y = np.meshgrid(np.arange(0,int(img.meta['specgriddx'])*int(img.meta['specgridnx']),int(img.meta['specgriddx'])),
                      np.arange(0,int(img.meta['specgriddy'])*int(img.meta['specgridny']),int(img.meta['specgriddy'])))
    X=X+2.5
    Y=Y+2.5

    xnm = X[xmin:xmax,ymin:ymax]/img.img_pixels.x*img.size.x/10
    ynm = Y[xmin:xmax,ymin:ymax]/img.img_pixels.y*img.size.y/10

    topo = img.img_array_list[0]
    topo = subtract_plane(topo.copy())

    fig = plt.figure("all_params",figsize=(13,6))
    axes = []

    patt = re.compile(r'(?:specgrid(\S*)[0-9]{4})')

    pdb.set_trace()

    #sort the files so we always see the parameter spreads in the same order
    # TODO: make it so it goes e, w, q, a, b, c  from top to bottom L --> R
    dat_files = sorted(dat_files, key=lambda x: x.replace("$","z")[-5])
    for n,f in enumerate(dat_files):
        this_limit_file = [d for d in limit_files if f[0:-4] in d][0]
        limits = np.loadtxt(os.path.join(path,this_limit_file))
        this_dat_file = np.loadtxt(os.path.join(path,f))
        s = f.split("_")
        title = [s[-2].replace("$\\epsilon",r"$\epsilon$") if "$\\epsilon" in s[-2] else s[-1].strip(".txt")][0]


        diff = re.search(patt, f).group(1)

        show_radial_decay(this_dat_file, xnm, ynm, title, diff)

        plt.figure("all_params")
        n=n+4
        ax = plt.subplot2grid((2,5), (n%2,int(n/2)))
        im = ax.imshow(this_dat_file,
                vmin=limits[0], vmax=limits[1], aspect="equal")

        axes.append(ax)
        divider =  make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')
        ax.set_title(title)
        # cax = fig.add_axes([0.27, 0.8, 0.5, 0.05])
        # fig.colorbar(ax=ax)

    ax = plt.subplot2grid((2,5),(0,0), colspan=2, rowspan=2)
    axes[0].sharex(axes[1])
    axes[2].sharex(axes[3])
    axes[4].sharex(axes[5])

    axes[1].sharey(axes[3])
    axes[3].sharey(axes[5])
    axes[0].sharey(axes[2])
    axes[2].sharey(axes[4])

    for a in axes[0:5]:
        a.set_yticks([])
        a.set_xticks([])

    # axes[1].set_ylabel([1])
    # axes[2].set_ylabel([2])
    # axes[3].set_ylabel([3])

    ax.imshow(topo,extent=[0,img.size.x/10.,0,img.size.y/10])
    # ax.scatter(xnm,ynm,s=0.1)
    ax.add_patch(Rectangle((xnm[0][0], ynm[0][0]), xnm[0][-1]-xnm[0][0], ynm[-1][-1]-ynm[0][0],
                 edgecolor = 'red',
                 facecolor = 'blue',
                 linestyle='dashed',
                 fill=False,
                 lw=5))
    ax.set_title("Topography")
    ax.set_xlabel("nm")
    ax.set_ylabel("nm")
    plt.tight_layout(pad=0, w_pad=0, h_pad=-1)
    # ax.get_xaxis_transform().inverted().transform([30,70])

    # plt.suptitle("Fitting Kondo resonance on central Co atom to Fano function")
    # plt.savefig("/Users/akipnis/Desktop/grid_kondo_fit.pdf")
    plt.show()

if __name__=="__main__":
    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 16}

    matplotlib.rc('font', **font)
    path = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/grid analysis"
    d = os.listdir(path)
    files = sorted([f for f in d if ".txt" in f and "trial" in f])
    plot_grid_fit_params(files, 32,63,33,62)


# dir = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/spatial extent Kondo plots/width comparison"
# files = os.listdir(dir)
#
# files = [f for f in files if f[-9:] =="width.txt" or f[-11:]=="default.txt"]
# str(os.path.join(dir,files[0]))
#
# d = pd.DataFrame()
#
# plt.figure(figsize=(9,6))
# for f in files:
#     data = pd.read_csv(os.path.join(dir,f), delimiter="\t", index_col=False)
#     # data = data[data["radius"]<2.5]
#     qcutoff = 1.5
#     dcutoff = 0.5
#     # cut_data = data[np.abs(data["q"])<qcutoff]
#     # cut_data = cut_data[cut_data["dist"]<dcutoff]
#     plt.scatter(data["dist"], data["w"], c=data["radius"],cmap="summer")
#     plt.clim(2.4,4.8)
    # plt.scatter(cut_data["radius"],cut_data["w"],
    #             s=50/cut_data["dist"], alpha=0.2)
