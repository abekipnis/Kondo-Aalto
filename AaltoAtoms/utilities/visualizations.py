import numpy as np
import matplotlib.pyplot as plt
import os, socket, sys, createc, pdb, matplotlib
from operator import sub
from .find_atom_positions import CircCorralData
from .Kondo_data_analysis.read_vertfile import Spec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import interp1d

from AaltoAtoms import CircCorralData
from AaltoAtoms import Spec

def create_waterfall():
    dir = r"Y:\labdata\Createc\STMDATA\Ag(111)\2022-03 Co Kondo corrals\04-11 Ag Co"
    wfall = {"A220411.133438.dat": "A220411.134351.L0013.VERT",
    "A220411.141241.dat": "A220411.141923.L0017.VERT",
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

        cache.append([S.bias_mv, S.dIdV/norm_val + len(wfall.keys())-n, colors[n], radius])
    return cache


#help(plt.cm)
def show_waterfall(cache):
    plt.figure(figsize=(8,8))
    colors = plt.cm.copper(np.linspace(0, 1, len(cache)))
    for n, c in enumerate(cache):
        plt.plot(c[0], c[1], color=colors[n], linewidth=5)

    plt.text(47, cache[0][1][0] + 0.2,"%1.1lf nm" %(np.round(cache[0][3],1))) # 11 nm
    plt.text(47, cache[len(cache)//2][1][0] - 0.6,"%1.1lf nm" %(np.round(cache[len(cache)//2][3],1)))
    plt.text(47, cache[-1][1][0] - 0.9,"%1.1lf nm" %(np.round(cache[-1][3],1))) # 3.6 nm
    plt.xlim(-80,80)
    plt.yticks([])
    plt.xlabel("Bias (mV)")
    plt.ylabel(r"$dI/dV$ (a.u.)")
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\Co-Ag-spectrum-waterfall.pdf")
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\Co-Ag-spectrum-waterfall.png")
    plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\Co-Ag-spectrum-waterfall.svg")

    plt.show()

def show_line_spectrum(datfile, line_directory, spectrum_timestamp, biasmin, biasmax):
    p = line_directory
    imf = datfile
    image = createc.DAT_IMG(imf)

    c = CircCorralData(imf, imf.split("/")[-1])
    files = os.listdir(p)
    files = sorted([f for f in files if f[-4:] =="VERT"])

    files = [f for f in files if spectrum_timestamp in f]

    # create lambda function to quickly join file paths
    pj = lambda f: os.path.join(p, f)
    specs = [Spec(pj(f)) for f in files]


    dIdVs = [s.dIdV/np.mean(s.dIdV) for s in specs]
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

    im = image._crop_img(image.img_array_list[0][:][:])
    zconst = image.zPiezoConst # angstroms/V

    # DAC conversion factor from Createc documentation from DAC to volts
    DAC_V_conv = interp1d([-524288, 524287],[-10,10])

    # then convert from volts to angstroms
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
    im1 = ax1.matshow(np.array(dIdVs).T,
                      extent=[0,dist/10., biasmin, biasmax],
                      aspect="auto",
                      cmap=plt.get_cmap("plasma")) #making aspect larger makes image skinnier
    ax1.set_xlabel("nm")
    ax1.xaxis.set_ticks_position("bottom")
    ax1.set_ylabel("mV")
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

    im2 = ax2.imshow(im, extent=[0,x_nm,y_nm,0],aspect="equal")
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(im2, cax=cax, orientation='vertical')#,fraction=0.046, pad=0.04)
    cbar.set_label("nm")
    ax2.plot(np.array(xlocs)-image.offset[0]/10+x_nm/2,np.array(ylocs)-image.offset[1]/10,"r")
    c.nm_to_pix(np.array(xlocs)-image.offset[0]/10+x_nm/2)

    ax2.set_xlabel("nm")
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

    # save the line spectrum plot in the same directory with the .dat file
    # save as pdf
    file_name = os.path.basename(imf).strip('.dat') + '_line_spectrum_%s.pdf' %(spectrum_timestamp)
    plt.savefig(os.path.join(os.path.dirname(imf), file_name))

    # save as png
    file_name = os.path.basename(imf).strip('.dat') + '_line_spectrum_%s.png' %(spectrum_timestamp)
    plt.savefig(os.path.join(os.path.dirname(imf), file_name))
    plt.show()

if __name__=="__main__":

    cache = create_waterfall()
    matplotlib.rcParams.update({'font.size': 22})
    show_waterfall(cache)
