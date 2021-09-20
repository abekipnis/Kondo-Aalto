
import numpy as np

import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import pdb
import scipy.signal
import plotly.offline
# %pylab inline
import plotly.graph_objects as go
import pandas as pd
import gif
import matplotlib.animation as animation

filename=r"S:\Createc_new\STMDATA\Ag\Small Kondo corrals\Ag 2021-08-13 2p5 nm radius\grid\Createc2_210814.214635.specgrid"
filename=r'/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-16 2p5 nm radius empty/Createc2_210816.223358.specgrid'
# filename = r'/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-12 4p5 nm radius/grid/Createc2_210813.001749.specgrid'

# read size of image from /Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-16 2p5 nm radius empty/Createc2_210816.223358.specgrid.dat

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

# pdb.set_trace()
specdata = data.reshape(85,85,len(specvz3),int(len(data)/85/85/len(specvz3)))

from matplotlib.animation import FuncAnimation

def animate_cube(cube_array, cut=True, mn=0, sd=0, interval=75, cmap='hot'):
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

    ptx, pty = 40, 40 # plot the LDOS from this pixel on the right
    ptx2, pty2 = 30, 30

    cut = 3 # cut these points out from teh beginning and end of the pixel array

    fig, (ax1, ax2) = plt.subplots(1, 2)# gridspec_kw={'width_ratios': [2, 1]})
    # ax1 = plt.subplot2grid((1,3), (0, 1), colspan = 1)
    # ax2 = plt.subplot2grid((1,3), (0, 2), colspan = 2)
    mn = np.mean(cube_array[-1][cut:-cut,cut:-cut])
    sd = np.std(cube_array[-1][cut:-cut,cut:-cut])
    img = ax1.imshow(cube_array[-1][cut:-cut,cut:-cut], animated=True,  cmap=cmap, vmax=mn+3*sd, vmin=mn-3*sd)

    ax1.scatter([ptx-cut],[pty-cut],)
    ax1.scatter([ptx2-cut],[pty2-cut],)

    cbar = fig.colorbar(img, ax=ax1, shrink=0.6)

    title = ax1.text(0.5,0.85, 'V= %1.2lf mV' %(specvz3[0,0]), bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                 transform=ax1.transAxes, ha="center")

    ax2.plot(specvz3[:,-1],cube_array[:,ptx, pty])
    ax2.plot(specvz3[:,-1],cube_array[:,ptx2, pty2])


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
        ax2.plot(specvz3[:,0],cube_array[:,ptx, pty])
        ax2.plot(specvz3[:,0],cube_array[:,ptx2, pty2])

        ax2.axvline(specvz3[i,0], c='r')
        title2 = ax2.text(0.5,1.1, "dI/dV point spectrum at corral center", #bbox={'facecolor':'w', 'alpha':1, 'pad':5}
                     transform=ax2.transAxes, ha="center")
        ax2.set_xlabel("Bias (mV)")
        ax2.set_ylabel("dI/dV (a.u)")
        ax2.set_yticks([])
        return img,
    plt.suptitle("3.8 nm radius occupied corral", y=0.8)

    # fig.tight_layout()
    ani = animation.FuncAnimation(fig, updatefig, frames=cube_array.shape[0], interval=interval, blit=True)
    # plt.show()
    ani.save(filename+'_cube_movie.mp4', writer="ffmpeg", fps=28)

animate_cube(specdata[:,:,:,1].T, False)
