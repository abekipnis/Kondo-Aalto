import numpy as np
import matplotlib.pyplot as plt
from pdb import set_trace
import scipy.signal
import pandas as pd
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import importlib
read_vertfile = importlib.import_module("Kondo data analysis.read_vertfile")
import createc
import pdb
import numpy.ma as ma
from matplotlib.widgets import Slider, Button

from multiprocessing import Pool, freeze_support
from find_atom_positions import CircCorralData
class Grid:
    def __init__(self, file):
        self.file = file
        # Specgrid parameters are loaded and printed below
        f = open(file, "rb")
        a = np.fromfile(f, dtype=np.uint32,count=256)
        f.close
        f = open(file, "rb")
        b = np.fromfile(f, dtype=np.float32,count=256)
        f.close

        self.version = a[0]
        self.nx = a[1]
        self.ny = a[2]
        self.dx=a[3]
        self.dy=a[4]
        self.specxgrid=a[5]
        self.specygrid=a[6]
        self.vertpoints=a[7]
        self.vertmandelay=a[8]
        self.vertmangain=a[9]
        self.biasvoltage=b[10]
        self.tunnelcurrent=b[11]
        self.imagedatasize=a[12]
        self.specgriddatasize=a[13]
        self.specgridchan=a[14]
        self.specgridchannelselectval=a[15]
        self.specgriddatasize64=np.int64(a[17])
        self.specgriddatasize64=(self.specgriddatasize64 << 32) + a[16]
        self.xstart=a[18]
        self.xend=a[19]
        self.ystart=a[20]
        self.yend=a[21]
        self.specgridchannelselectval2=a[22]
        self.specgridnx=a[23]
        self.specgridny=a[24]
        self.specgriddx=a[25]
        self.specgriddy=a[26]
        self.specgridcenterx=a[27]
        self.specgridcentery = a[28]

        self.count3=self.vertpoints*3

        self.specvz = np.fromfile(f, dtype=np.float32,count=self.count3)
        self.specvz3 = self.specvz.reshape(self.vertpoints,3)
        self.data = np.fromfile(f, dtype=np.float32)
        f.close

        self.a, self.b = int(self.nx/self.specgriddx), int(self.ny/self.specgriddy)
        try:
            self.specdata = data.reshape(self.a,self.b,len(self.specvz3),int(len(self.data)/self.a/self.b/len(self.specvz3)))
        except:
            self.a, self.b = self.xend, self.yend
            self.specdata = self.data.reshape(self.a,self.b,len(self.specvz3),int(len(self.data)/self.a/self.b/len(self.specvz3)))
        self.cube_array = self.specdata[:,:,:,1].T

        _, self.xpix, self.ypix = self.cube_array.shape
        self.nmx, self.nmy = self.get_im_size()


    def get_im_size(self):
        f = open(self.file+".dat","rb")
        d = f.readlines()
        a = [str(b).split('\\t') for b in d]
        ls = [float(b.split("=")[-1].rstrip(" \\r\\n'")) for b in np.array(a[0:600]).T[0] if "Length" in b]
        return ls #in angstroms

    def fit_Fano_to_grid_data(self, xpixmin, xpixmax, ypixmin,ypixmax, offset):

        xpmx_nm = self.pix_to_nm(xpixmax)
        xpmn_nm = self.pix_to_nm(xpixmin)
        p = np.concatenate([    [   (self.specvz3[:,0][offset:-1],
                                    self.cube_array[:,i,j][offset:-1],
                                    -15,20)
                                    for i in range(xpixmin,xpixmax)]
                            for j in range(ypixmin,ypixmax)])

        with Pool() as pool:
            L = pool.starmap(read_vertfile.fit_data, p)

        nx = xpixmax - xpixmin
        ny = ypixmax - ypixmin
        Lm = np.array(L).reshape(nx,ny,4) # popt, pcov, sb, fit_dIdV
        m = ma.masked_array(Lm[:,:,0], mask=np.any(pd.isnull(Lm),axis=-1))
        c_max = "%1.2lf"
        labels = [r"$\epsilon_0$","w","q","a","b","c"]
        for n, l in enumerate(labels):
            ax = plt.subplot(111)
            plt.subplots_adjust(left=0.25, bottom=0.25)
            dat = np.concatenate(m.data.flatten()).reshape(nx,ny,6)[:,:,n]
            pdb.set_trace()
            img = ax.imshow(dat,
                       extent=[xpmn_nm/10., xpmx_nm/10.,xpmn_nm/10., xpmx_nm/10.,]);
            plt.title(l)
            plt.xlabel("nm")
            plt.ylabel("nm")
            cb = plt.colorbar(img);
            axcolor = 'lightgoldenrodyellow'
            ax_cmax  = plt.axes([0.25, 0.15, 0.65, 0.03])
            ax_cmin  = plt.axes([0.25, 0.05, 0.65, 0.03])

            nmnd = np.nanmin(dat)
            nmxd = np.nanmax(dat)
            if l=="q":
                nmnd=-5
                nmxd=5
            s_cmax = Slider(ax_cmax, 'max', nmnd, nmxd, valfmt=c_max)
            s_cmin = Slider(ax_cmin, 'min', nmnd, nmxd, valfmt=c_max)
            def update(val, s=None):
                _cmax = s_cmax.val
                _cmin = s_cmin.val
                img.set_clim(_cmin,_cmax)
                plt.draw()
            s_cmax.on_changed(update)
            s_cmin.on_changed(update)

            plt.show()

            plt.savefig(l+".png")

    def pix_to_nm(self, pix):
        assert self.xpix==self.ypix and self.nmx==self.nmy
        return pix * self.nmx/self.xpix

    def nm_to_pix(self, nm):
        assert self.xpix==self.ypix and self.nmx==self.nmy
        return nm*self.xpix/self.nmx*10

    def animate_cube(self, cut=True, mn=0, sd=0, interval=75, cmap='hot', ptx=0, pty=0, ptx2=0, pty2=0):
        '''
        animates a cube for visualisation.

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
        offset = 20 # because first N (20) points of spectrum are at 20mV


        # defining range of pixels over which to plot Fano fit in grid
        # xpixmin = ypixmin = 30
        # xpixmax = ypixmax = 70
        # self.fit_Fano_to_grid_data(xpixmin, xpixmax, ypixmin, ypixmax, offset)

        #for the 2.5nm radius empty corral these are at the middle / node
        # ptx, pty = 40, 40 # plot the LDOS from this pixel on the right
        # ptx2, pty2 = 30, 30

        cut = 3 # cut these points out from the beginning and end of the pixel array

        fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 1]})
        mn = np.mean(self.cube_array[-1][cut:-cut,cut:-cut])
        sd = np.std(self.cube_array[-1][cut:-cut,cut:-cut])
        img = ax1.imshow(self.cube_array[-1][cut:-cut,cut:-cut],
                        animated=True,
                        cmap=cmap,
                        vmax=mn+3*sd,
                        vmin=mn-3*sd,
                        extent=[0, self.nmx/10.,0, self.nmy/10.])
        ax1.set_xlabel("nm")
        ax1.set_ylabel("nm")

        # s1 and s2 are PathCollection objects
        # s1.get_sizes() looks at default s (default size is 36)
        s1 = ax1.scatter([ptx],[pty], s=20)
        s2 = ax1.scatter([ptx2],[pty2], s=20)
        cbar = fig.colorbar(img, ax=ax1, shrink=0.6)

        title = ax1.text(0.5,0.85, 'V= %1.2lf mV' %(self.specvz3[0,0]),
                    bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                    transform=ax1.transAxes, ha="center")

        # cut off some points from end of spectra because dIdV saturates at the end
        # probably due to not setting the bias to the starting point of the spectra
        # or some other effect
        # makes it hard to see the details in the spectrum
        # we have to either manually change the axes limits or just cut off part of the data in the visualisation

        offset = 10

        px, py, px2, py2 = map(int,map(self.nm_to_pix,[ptx, pty, ptx2, pty2]))
        ax2.plot(self.specvz3[:,-1][offset:-1], self.cube_array[:,px, py][offset:-1])
        ax2.plot(self.specvz3[:,-1][offset:-1], self.cube_array[:,px2, py2][offset:-1])

        ax2.yaxis.tick_right()
        ax2.set_xlabel("Bias (mV)")
        ax2.set_ylabel("dI/dV (a.u)")
        ax1.set_title("LDOS(V,r)")
        cbar.set_ticks([])
        ax2.set_adjustable('box')
        ax2.set_aspect(0.1)
        def updatefig(i):
            i = 520-i #to do it in reverse
            d = self.cube_array[i][cut:-cut,cut:-cut]

            title.set_text('V= %1.2lf mV' %(self.specvz3[i,0]))
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
            ax2.plot(self.specvz3[:,0][offset:-1], self.cube_array[:,px, py][offset:-1])
            ax2.plot(self.specvz3[:,0][offset:-1], self.cube_array[:,px2, py2][offset:-1])

            ax2.axvline(self.specvz3[i,0], c='r')
            title2 = ax2.text(0.5,1.1, "dI/dV point spectra from grid", #bbox={'facecolor':'w', 'alpha':1, 'pad':5}
                         transform=ax2.transAxes, ha="center")
            ax2.set_xlabel("Bias (mV)")
            ax2.set_ylabel("dI/dV (a.u)")
            ax2.set_yticks([])
            return img,
        plt.suptitle("2.5 nm radius occupied corral", y=0.95)

        ani = animation.FuncAnimation(fig, updatefig, frames=self.cube_array.shape[0], interval=interval, blit=True)
        plt.show()
        ani.save(self.file+'_cube_movie.mp4', writer="ffmpeg", fps=28)



if __name__ == "__main__":
    # read size of image from .specgrid.dat file
    dir = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/"
    filename= dir +r"Ag 2021-08-13 2p5 nm radius/grid/Createc2_210814.214635.specgrid"

    g = Grid(filename)
    # filename = dir +r'Ag 2021-08-16 2p5 nm radius empty/Createc2_210816.223358.specgrid'
    # filename = dir+r'/Ag 2021-08-12 4p5 nm radius/grid/Createc2_210813.001749.specgrid'


    freeze_support()
    g.animate_cube(cut=False, ptx=3.3, pty=3.6, ptx2=3.3, pty2=4)
