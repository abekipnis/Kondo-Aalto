
import numpy as np

import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
from pdb import set_trace
import scipy.signal
import plotly.offline
# %pylab inline
import plotly.graph_objects as go
import pandas as pd
import gif
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import createc
import pdb

def animate_cube(cube_array, filename, cut=True, mn=0, sd=0, interval=75, cmap='hot', ptx=0, pty=0, ptx2=0, pty2=0):
    '''
    animates a python cube for quick visualisation.

    INPUT:
        cube_array  : name of 3D numpy array that needs to be animated.
        cut         : trims pixels off of the images edge to remove edge detector effects.
                      Default = True as 0 returns empty array.
        mn          : mean of the cube | Used for contrast
        sd          : std of the cube  | Used for contrast
        interval    : #of ms between each frame.
        cmap        : colormap. Default='hot'

    OUTPUT:
        animated window going through the cube.

    '''
    # cube_array = scipy.signal.detrend(cube_array,axis=2)

    # fig = plt.figure(figsize=(12,8))

    #for the 2.5nm radius empty corral these are at the middle / node
    # ptx, pty = 40, 40 # plot the LDOS from this pixel on the right
    # ptx2, pty2 = 30, 30
    nmx, nmy = get_im_size(filename)
    _, xpix, ypix = cube_array.shape

    def pix_to_nm(pix):
        assert xpix==ypix and nmx==nmy
        return pix * nmx/xpix

    def nm_to_pix(nm):
        assert xpix==ypix and nmx==nmy
        return nm*xpix/nmx*10

    cut = 3 # cut these points out from teh beginning and end of the pixel array

    fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 1]})
    # ax1 = plt.subplot2grid((1,3), (0, 1), colspan = 1)
    # ax2 = plt.subplot2grid((1,3), (0, 2), colspan = 2)
    mn = np.mean(cube_array[-1][cut:-cut,cut:-cut])
    sd = np.std(cube_array[-1][cut:-cut,cut:-cut])
    img = ax1.imshow(cube_array[-1][cut:-cut,cut:-cut],
                    animated=True,
                    cmap=cmap,
                    vmax=mn+3*sd,
                    vmin=mn-3*sd,
                    extent=[0, nmx/10.,0, nmy/10.])
    ax1.set_xlabel("nm")
    ax1.set_ylabel("nm")
    # s1 and s2 are PathCollection objects
    # can use s1.get_sizes() to look at the default size
    s1 = ax1.scatter([ptx],[pty], s=20)
    s2 = ax1.scatter([ptx2],[pty2], s=20) #default size is 36
    cbar = fig.colorbar(img, ax=ax1, shrink=0.6)

    title = ax1.text(0.5,0.85, 'V= %1.2lf mV' %(specvz3[0,0]),
                bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                transform=ax1.transAxes, ha="center")

    # cut off some points from end of spectra because dIdV saturates at the end
    # probably due to not setting the bias to the starting point of the spectra
    # or some other effect
    # makes it hard to see the details in the spectrum
    # we have to either manually change the axes limits or just cut off part of the data in the visualisationquit()

    offset = 10

    px, py, px2, py2 = map(int,map(nm_to_pix,[ptx, pty, ptx2, pty2]))
    ax2.plot(specvz3[:,-1][offset:-1],cube_array[:,px, py][offset:-1])
    ax2.plot(specvz3[:,-1][offset:-1],cube_array[:,px2, py2][offset:-1])

    ax2.yaxis.tick_right()
    ax2.set_xlabel("Bias (mV)")
    ax2.set_ylabel("dI/dV (a.u)")
    ax1.set_title("LDOS(V,r)")
    cbar.set_ticks([])
    ax2.set_adjustable('box')
    ax2.set_aspect(0.1)
    def updatefig(i):
        i = 520-i #to do it in reverse
        d = cube_array[i][cut:-cut,cut:-cut]

        title.set_text('V= %1.2lf mV' %(specvz3[i,0]))
        mn = np.mean(d)
        sd = np.std(d)
        ax1.images[0].colorbar.remove()
        img.set_array(d)#, )
        img.set_clim(vmax=mn+3*sd, vmin=mn-3*sd)
        #
        cbar = fig.colorbar(img, ax=ax1, shrink=0.6)
        cbar.set_ticks([])
        cbar.update_normal(img)
        # cbar.clim(vmin=d.min(),vmax=d.max())
        # cbar.draw_all()
        # ax2.remove()
        ax2.clear()
        ax2.plot(specvz3[:,0][offset:-1],cube_array[:,px, py][offset:-1])
        ax2.plot(specvz3[:,0][offset:-1],cube_array[:,px2, py2][offset:-1])

        ax2.axvline(specvz3[i,0], c='r')
        title2 = ax2.text(0.5,1.1, "dI/dV point spectra from grid", #bbox={'facecolor':'w', 'alpha':1, 'pad':5}
                     transform=ax2.transAxes, ha="center")
        ax2.set_xlabel("Bias (mV)")
        ax2.set_ylabel("dI/dV (a.u)")
        ax2.set_yticks([])
        return img,
    plt.suptitle("2.5 nm radius occupied corral", y=0.95)

    # fig.tight_layout()
    ani = animation.FuncAnimation(fig, updatefig, frames=cube_array.shape[0], interval=interval, blit=True)
    # plt.show()
    ani.save(filename+'_cube_movie.mp4', writer="ffmpeg", fps=28)


def get_im_size(gridfile):
    f = open(gridfile+".dat","rb")
    d = f.readlines()
    a = [str(b).split('\\t') for b in d]
    ls = [float(b.split("=")[-1].rstrip(" \\r\\n'")) for b in np.array(a[0:600]).T[0] if "Length" in b]
    return ls #in angstroms


if __name__ == "__main__":
    # read size of image from .specgrid.dat file
    dir = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/"
    filename= dir +r"Ag 2021-08-13 2p5 nm radius/grid/Createc2_210814.214635.specgrid"
    # datfile = filename
    # filename = dir +r'Ag 2021-08-16 2p5 nm radius empty/Createc2_210816.223358.specgrid'
    # filename = dir+r'/Ag 2021-08-12 4p5 nm radius/grid/Createc2_210813.001749.specgrid'

    # In this section the specgrid parameters are loaded and printed below
    f = open(filename, "rb")
    a = np.fromfile(f, dtype=np.uint32,count=256)
    f.close

    f = open(filename, "rb")
    b = np.fromfile(f, dtype=np.float32,count=256)
    f.close

    version = a[0]
    nx = a[1]
    ny = a[2]
    dx=a[3]
    dy=a[4]
    specxgrid=a[5]
    specygrid=a[6]
    vertpoints=a[7]
    vertmandelay=a[8]
    vertmangain=a[9]
    biasvoltage=b[10]
    tunnelcurrent=b[11]
    imagedatasize=a[12]
    specgriddatasize=a[13]
    specgridchan=a[14]
    specgridchannelselectval=a[15]
    specgriddatasize64=np.int64(a[17])
    specgriddatasize64=(specgriddatasize64 << 32) + a[16]
    xstart=a[18]
    xend=a[19]
    ystart=a[20]
    yend=a[21]
    specgridchannelselectval2=a[22]
    specgridnx=a[23]
    specgridny=a[24]
    specgriddx=a[25]
    specgriddy=a[26]
    specgridcenterx=a[27]
    specgridcentery = a[28]

    count3=vertpoints*3

    specvz = np.fromfile(f, dtype=np.float32,count=count3)
    specvz3 = specvz.reshape(vertpoints,3)
    data = np.fromfile(f, dtype=np.float32)
    f.close

    a, b = int(nx/specgriddx), int(ny/specgriddy)
    try:
        specdata = data.reshape(a,b,len(specvz3),int(len(data)/a/b/len(specvz3)))
    except:
        a, b = xend, yend
        specdata = data.reshape(a,b,len(specvz3),int(len(data)/a/b/len(specvz3)))
    animate_cube(specdata[:,:,:,1].T, filename, cut=False, ptx=3.3, pty=3.6, ptx2=3.3, pty2=4)
