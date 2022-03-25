import numpy as np
import matplotlib.pyplot as plt
import os, socket, sys, createc, pdb, matplotlib
from operator import sub
from .find_atom_positions import CircCorralData
from .Kondo_data_analysis.read_vertfile import Spec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import interp1d

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
    file_name = os.path.basename(imf).strip('.dat') + '_line_spectrum_%s.pdf' %(spectrum_timestamp)
    plt.savefig(os.path.join(os.path.dirname(imf), file_name))
    plt.show()
