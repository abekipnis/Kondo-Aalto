# some_file.py
import sys
from matplotlib.widgets import Slider, Button

# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/Kondo data analysis')
from read_vertfile import Spec
import numpy as np
import createc
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from operator import sub
from scipy.interpolate import interp1d
# def get_aspect(ax):
#     # Total figure size
#     figW, figH = ax.get_figure().get_size_inches()
#     # Axis size on figure
#     _, _, w, h = ax.get_position().bounds
#     # Ratio of display units
#     disp_ratio = (figH * h) / (figW * w)
#     # Ratio of data units
#     # Negative over negative because of the order of subtraction
#     data_ratio = sub(*ax.get_ylim()) / sub(*ax.get_xlim())
#
#     return disp_ratio / data_ratio

import os
p = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-13 3p8 nm radius/3p8nm pm20mV line/"
image = createc.DAT_IMG("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-13 3p8 nm radius/Createc2_210813.102220.dat")

files = os.listdir(p)
files = sorted([f for f in files if f[-4:] =="VERT"])

specs = [Spec(p+f) for f in files][2:-5]
specs[0].__dict__.keys()
specs[0].__dict__['fname']

image.__dict__['bias']
image.__dict__['current']


specs[5].__dict__['LockinAmpl']
specs[6].__dict__['FBLogiset']
specs[5].__dict__['biasVoltage']
specs[5].__dict__['bias_offset']
specs[5].__dict__['LockinFreq']

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
zconst

DAC_V_conv = interp1d([-524288, 524287],[-10,10])

im = DAC_V_conv(im)*zconst
image.meta['fblogiset']
subtract_plane(im, im.shape[0], im.shape[1])
im = im-np.min(im)
import matplotlib

###
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 20}

matplotlib.rc('font', **font)
fig, (ax1,ax2) = plt.subplots(ncols=2, figsize=(12,6))
im1 = ax1.matshow(np.array(dIdVs).T,extent=[0,dist/10.,-20,20], aspect="auto",cmap=plt.get_cmap("plasma")) #making aspect larger makes image skinnier
ax1.set_xlabel("nm")
ax1.xaxis.set_ticks_position("bottom")
ax1.set_ylabel(r"Bias ($mV$)")
# ax1.set_title("$dI/dV$ spectra")
# ax1.hlines(-67,0,25,"r")
# divider = make_axes_locatable(ax1)
# divider.__dict__.keys()
# divider.__dict__['_axes'].__dict__

# cax = divider.append_axes('right', size='5%', pad=0.05)
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
# get_aspect(ax1)

im2 = ax2.imshow(im, extent=[0,x_nm,y_nm,0],aspect="equal")
divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad=0.05)
cbar = fig.colorbar(im2, cax=cax, orientation='vertical')#,fraction=0.046, pad=0.04)
cbar.set_label("nm")
# cbar.set_ticks([])
ax2.plot(np.array(xlocs)-image.offset[0]/10+x_nm/2,np.array(ylocs)-image.offset[1]/10,"r")
# ax2.set_title("Topography")
ax2.set_xlabel("nm")

# fig.suptitle("Line spectrum from Ag(111) step edge")
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
plt.text(  # position text relative to Figure
    0.025, 0.95, '(a)',
    ha='left', va='top',
    transform=fig.transFigure
)
plt.text(  # position text relative to Figure
    0.475, 0.95, '(b)',
    ha='left', va='top',
    transform=fig.transFigure
)
# plt.show()
plt.savefig("/Users/akipnis/Desktop/Createc2_210813.102220_3p8_pm20mV.pdf")
# image.offset[0]/10
