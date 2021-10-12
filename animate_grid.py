import numpy as np
import matplotlib.pyplot as plt
from pdb import set_trace
import scipy.signal
import pandas as pd
import matplotlib.animation as animation
# from matplotlib.animation import FuncAnimation
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
            self.specdata = self.data.reshape(self.a,self.b,len(self.specvz3),int(len(self.data)/self.a/self.b/len(self.specvz3)))
        except:
            self.a, self.b = self.xend, self.yend
            self.specdata = self.data.reshape(self.a,self.b,len(self.specvz3),int(len(self.data)/self.a/self.b/len(self.specvz3)))
        self.cube_array = self.specdata[:,:,:,1].T

        _, self.xpix, self.ypix = self.cube_array.shape
        self.nmx, self.nmy = self.get_im_size()
        self.offset = 20 # because first N (20) points of spectrum are at 20mV

    def get_im_size(self):
        f = open(self.file+".dat","rb")
        d = f.readlines()
        a = [str(b).split('\\t') for b in d]
        ls = [float(b.split("=")[-1].rstrip(" \\r\\n'")) for b in np.array(a[0:600]).T[0] if "Length" in b]
        return ls #in angstroms

    def fit_Fano_to_grid_data(self, xpixmin, xpixmax, ypixmin, ypixmax):

        xpmx_nm = self.pix_to_nm(xpixmax)
        xpmn_nm = self.pix_to_nm(xpixmin)

        # fitting the "normal" way, with no parameters fixed (only bounded)
        # p = np.concatenate([    [   (self.specvz3[:,0][offset:-1],
        #                             self.cube_array[:,i,j][offset:-1],
        #                             -15,20)
        #                             for i in range(xpixmin,xpixmax)]
        #                     for j in range(ypixmin,ypixmax)])
        # with Pool() as pool:
        #     L = pool.starmap(read_vertfile.fit_data, p)

        # fitting with params fixed. np.nan means 'unfixed', otherwise fixed at given value
        # [e0, w, q, a, b, c]

        p_fixed = [8.3, np.nan, np.nan, np.nan, -50, np.nan] #fix e0 at 8.3 mV
        pe0 = np.concatenate([    [   (self.specvz3[:,0][self.offset:-1],
                                    self.cube_array[:,i,j][self.offset:-1],
                                    -15, 20, p_fixed)
                                    for i in range(xpixmin,xpixmax)]
                            for j in range(ypixmin,ypixmax)])

        print(read_vertfile.fit_data_fixed_vals(*pe0[int((xpixmax-xpixmin)/2)+(xpixmax-xpixmin)*int((ypixmax-ypixmin)/2)]))

        with Pool() as pool:
            L = pool.starmap(read_vertfile.fit_data_fixed_vals, pe0)

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
            # pdb.set_trace()
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

            # save the data
            f = datetime.today().strftime('%Y-%m-%d') + "_" + l + ".dat"
            d = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/grid analysis"
            f = os.path.join(d,f)
            np.savetxt(data, f)

            plt.show()

            # plt.savefig(l+".png")

    def pix_to_nm(self, pix):
        assert self.xpix==self.ypix and self.nmx==self.nmy
        return pix * self.nmx/self.xpix

    def nm_to_pix(self, nm):
        assert self.xpix==self.ypix and self.nmx==self.nmy
        return nm*self.xpix/self.nmx*10

    def animate_cube(self, mn=0, sd=0, interval=75, cmap='hot', plotpoints=[[0,0],[0,0]]):
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
            plotpoints  : [[x1,y1], [x2,y2]]

        OUTPUT:
            animated window going through the cube.

        '''
        class PauseAnimation:
            def __init__(self, g, plotpoints):
                self.g = g
                gsk = {'width_ratios': [2, 1]}
                self.pts = plotpoints
                self.fig, (self.ax1, self.ax2) = plt.subplots(1, 2, gridspec_kw=gsk)



                mn = np.mean(self.g.cube_array[-1])
                sd = np.std(self.g.cube_array[-1])
                self.img = self.ax1.imshow(self.g.cube_array[-1],
                                animated=True,
                                cmap=cmap,
                                vmax=mn+3*sd,
                                vmin=mn-3*sd,
                                extent=[0, g.nmx/10.,0, g.nmy/10.])
                self.ax1.set_xlabel("nm")
                self.ax1.set_ylabel("nm")

                # s is a PathCollection object
                # s.get_sizes() looks at default size (36)
                for p in plotpoints:
                    s = self.ax1.scatter([p[0]],[p[1]], s=15)
                self.cbar = self.fig.colorbar(self.img, ax=self.ax1, shrink=0.6)

                self.title = self.ax1.text(0.5,0.85, 'V= %1.2lf mV' %(self.g.specvz3[0,0]),
                            bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                            transform=self.ax1.transAxes, ha="center")

                # cut off points from end of spectra because dIdV saturates at end
                # probably due to not setting bias to the starting point of spectra
                # or some other effect
                # makes it hard to see the details in the spectrum
                # we have to either manually change the axes limits or just cut off part of the data in the visualisation

                self.pts = np.array(list(map(int,map(self.g.nm_to_pix, np.array(plotpoints).flatten())))).reshape(len(plotpoints),2)
                # pdb.set_trace()
                self.pts[:,0] = self.g.xpix-self.pts[:,0]
                self.pts[:,1] = self.g.ypix-self.pts[:,1]
                for p in self.pts:
                    self.ax2.plot(self.g.specvz3[:,-1][self.g.offset:-1], self.g.cube_array[:, p[1], p[0]][self.g.offset:-1])

                self.ax2.yaxis.tick_right()
                self.ax2.set_xlabel("Bias (mV)")
                self.ax2.set_ylabel("dI/dV (a.u)")
                self.ax1.set_title("LDOS(V,r)")
                self.cbar.set_ticks([])
                self.ax2.set_adjustable('box')
                self.ax2.set_aspect(0.1)
                self.paused = False


                plt.suptitle("2.5 nm radius occupied corral", y=0.95)

                self.animation = animation.FuncAnimation(self.fig, self.updatefig, frames=g.cube_array.shape[0], interval=interval, blit=True)
                self.fig.canvas.mpl_connect('button_press_event', self.toggle_pausefig)

                # self.animation.save(g.file+'_cube_movie.mp4', writer="ffmpeg", fps=28)

            def toggle_pausefig(self, *args, **kwargs):
                if self.paused:
                    self.animation.resume()
                else:
                    self.animation.pause()
                self.paused = not self.paused

            def updatefig(self, i):
                i = len(self.g.specvz3)-1-i #to do it in reverse
                d = self.g.cube_array[i]

                self.title.set_text('V= %1.2lf mV' %(self.g.specvz3[i,0]))
                mn = np.mean(d)
                sd = np.std(d)
                self.ax1.images[0].colorbar.remove()
                self.img.set_array(d)
                self.img.set_clim(vmax=mn+6*sd, vmin=mn-3*sd)


                self.cbar = self.fig.colorbar(self.img, ax=self.ax1, shrink=0.6)
                self.cbar.set_ticks([])
                self.cbar.update_normal(self.img)

                self.ax2.clear()
                for p in self.pts:
                    self.ax2.plot(self.g.specvz3[:,0][self.g.offset:-1], self.g.cube_array[:,p[1], p[0]][self.g.offset:-1])

                self.ax2.axvline(self.g.specvz3[i,0], c='r')
                self.title2 = self.ax2.text(
                0.5,1.1, "dI/dV point spectra from grid",
                             transform=self.ax2.transAxes, ha="center")
                             #bbox={'facecolor':'w', 'alpha':1, 'pad':5}
                self.ax2.set_xlabel("Bias (mV)")
                self.ax2.set_ylabel("dI/dV (a.u)")
                self.ax2.set_yticks([])
                return self.img,

        P = PauseAnimation(self, plotpoints)


if __name__ == "__main__":
    # read size of image from .specgrid.dat file
    dir = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/"

    # there are three successful grids from the first data set

    #this is the grid with a Cobalt in the center
    filename= dir +r"Ag 2021-08-13 2p5 nm radius/grid/Createc2_210814.214635.specgrid"
    g = Grid(filename)
    # range of pixels over which to plot Fano fit in grid
    xpixmin = ypixmin = 30
    xpixmax = ypixmax = 70

    # TODO: turn this into nm to more easily change

    g.fit_Fano_to_grid_data(xpixmin, xpixmax, ypixmin, ypixmax)

    g.animate_cube(plotpoints=[[3.3, 3.6],[3.3, 3.8], [3.3, 4.0], [3.3, 4.2], [3.3, 4.4], [3.3, 4.6], [3.3, 4.8]])
    plt.show()


    # filename = dir +r'Ag 2021-08-16 2p5 nm radius empty/Createc2_210816.223358.specgrid'
    # g = Grid(filename)
    # g.animate_cube(plotpoints=[[4.58, 4.36],[4.58, 4.86], [4.58, 5.36],[4.58,5.86],[4.58, 6.36],[4.58, 6.86]])
    # plt.show()

    #not sure whats happening with this grid - can't see the data when this code runs
    filename = dir+r'/Ag 2021-08-12 4p5 nm radius/grid/Createc2_210813.001749.specgrid'
    g = Grid(filename)

    g.animate_cube(plotpoints=[[1.58, 1.36]])

    plt.show()

    freeze_support()
