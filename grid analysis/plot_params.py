import matplotlib.pyplot as plt
import numpy as np
import os
import createc
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Rectangle

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


import matplotlib
if __name__=="__main__":
    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 16}

    matplotlib.rc('font', **font)

    path = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/grid analysis"
    d = os.listdir(path)
    files = sorted([f for f in d if ".txt" in f])

    topo_file = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-13 2p5 nm radius/grid/Createc2_210814.214635.specgrid.dat"
    img = createc.DAT_IMG(topo_file)
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

    topo = img.img_array_list[0]
    topo = subtract_plane(topo.copy())

    fig = plt.figure(figsize=(13,6))
    axes = []
    for n,f in enumerate(files[::2]):
        limits = np.loadtxt(os.path.join(path,files[n*2+1]))
        n=n+4
        ax = plt.subplot2grid((2,5), (n%2,int(n/2)))
        im = ax.imshow(np.loadtxt(os.path.join(path,f)),
                vmin=limits[0], vmax=limits[1], aspect="equal")
                #[n%2][int(n/2)]
        axes.append(ax)
        divider =  make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')
        s = f.split("_")
        title = [s[-2].replace("$\\epsilon",r"$\epsilon$") if "$\\epsilon" in s[-2] else s[-1].strip(".txt")][0]
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
    # ax._axes.__dict__
    # ax.__dict__
    ax.imshow(topo,extent=[0,img.size.x/10.,0,img.size.y/10])
    # ax.scatter(X[30:70,30:70],Y[30:70,30:70],s=0.1)
    # img.img_pixels.x
    xnm = X[32:63,33:62]/img.img_pixels.x*img.size.x/10
    ynm = Y[32:63,33:62]/img.img_pixels.y*img.size.y/10
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

    plt.suptitle("Fitting Kondo resonance on central Co atom to Fano function")
    # plt.savefig("/Users/akipnis/Desktop/grid_kondo_fit.pdf")
    plt.show()
