from scipy import *
import numpy as np
import itertools
import scipy.constants, scipy.special
import matplotlib.pyplot as plt

import matplotlib
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)

radius1 = 2*1e-9 #meters
radius2 = 5*1e-9 #meters
mstar = 0.4*scipy.constants.electron_mass
radius_range = np.arange(radius1,radius2+0.01*radius1,0.01*radius1)

theta = np.arange(0,2*np.pi,0.01*2*np.pi)

def get_azimuthal_function(l):
    ret = np.exp(1j*l*theta)

    if l<0:
        ret*=np.exp(1j*np.pi/2/np.abs(l))
    return ret

def get_radial_function(n,l,radius):
    k_nl = scipy.special.jn_zeros(n,1)/radius
    # plt.plot(np.arange(0,radius1,0.01*radius1),scipy.special.jv(1,k_nl*np.arange(0,radius1,0.01*radius1)))
    return scipy.special.jv(n,k_nl*np.arange(0,radius,0.01*radius))
# plt.subplot(projection="polar")
# plt.pcolormesh(theta,np.arange(0,radius,0.01*radius), (mesh**2).real)

def plot(n,l,radius):
    X, Y = np.meshgrid(get_azimuthal_function(l),get_radial_function(n,l,radius))

    from mpl_toolkits.mplot3d import Axes3D
    # fig = plt.figure()
    # ax = Axes3D(fig)
    # plt.subplot(projection="polar")
    ret = X*Y
    # plt.pcolormesh(theta,np.arange(0,radius,0.01*radius), (ret**2).real)
    return ret


def get_modes(mstar: float, E0_V: float, radius_range_nm: list, n_modes: int):
    h_ = scipy.constants.hbar
    e_ = scipy.constants.elementary_charge
    E = lambda r: h_**2*(scipy.special.jn_zeros(0,n_modes)/r)**2/(2*mstar)/e_-E0_V
    return np.array([E(r) for r in radius_range_nm])


def show_particle_in_box_eigenmodes_vs_corral_radius():
    fig, (ax1,ax3) = plt.subplots(figsize=(11,8), nrows=2, sharex=True)
    ax1.plot(radius_range*1e9,get_modes(mstar, -67, radius_range,3))
    ylim = ax1.get_ylim()
    # ax1.vlines(2.5,-10,10,"r")
    # ax1.vlines(3.8,-10,10,"r")
    # ax1.vlines(4.5,-10,10,"r")
    ax1.set_ylim(ylim)
    # plt.axvline(2.5, c="b")

    for n in range(0,3):
        print(n)
        print(list(range(-n,n+1)))
        for l in range(-n,n+1):
            print("\t",l)
            # ax1.add_hline(2.5)
            mesh = plot(n,l,radius1)

            # These are in unitless percentages of the figure size. (0,0 is bottom left)
            left, bottom, width, height = [0.65+.095*l, 0.6+0.12*n, 0.1, 0.1]
            ax2 = fig.add_axes([left, bottom, width, height],projection="polar")
            ax2.pcolormesh(theta,np.arange(0,radius1,0.01*radius1), (mesh**2).real)
            for r_label in ax2.get_xticklabels():
                r_label.set_text("")
            ax2.set_xticklabels([])
            ax2.set_yticklabels([])
            ax2.__dict__
            ax2.bbox
            # for ax2 in axes:
            ax2.tick_params(color='green', labelcolor='green')
            for spine in ax2.spines.values():
                spine.set_edgecolor('green')

    # left, bottom, width, height = [0.4, 0.3, 0.1, 0.1]
    # ax2 = fig.add_axes([left, bottom, width, height],projection="polar")
    # ax2.pcolormesh(theta,np.arange(0,radius1,0.01*radius1), (plot(2,1,radius1)**2).real)
    # for r_label in ax2.get_xticklabels():
    #     r_label.set_text("")
    # ax2.set_xticklabels([])
    # ax2.set_yticklabels([])

    # ax1.set_xlabel(r"Circular corral radius ($nm$)")
    ax1.set_ylabel(r"Eigenenergy ($eV$)")
    # ax1.text(2.2,  1.6,"n=0")
    # ax1.text(2.2,  0.7,"n=1")
    # ax1.text(2.2,  0.2,"n=2")
    ax1.legend(["n=0", "n=1","n=2"], loc="lower right")
    ax1.set_xlim([2,5])

    # fig.suptitle("Circular corral eigenmodes from particle in a box model\nUsing $m^*=$%1.2lf, $E_0$=-67 mV" %(mstar/scipy.constants.electron_mass))
    # plt.savefig("/Users/akipnis/Desktop/circular_corral_eigenmodes.pdf")


    ###################################
    # occupied modes
    occ_data = np.loadtxt("/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/Corral eigenmode solvers (MATLAB)/Occupied_eigvals.txt",delimiter=",")
    r_occ_data = occ_data[0]
    eval_occ_data = occ_data[1]

    #####
    # fig, ax2 = plt.subplots(figsize=(9,6))
    ax3.plot(radius_range*1e9,[scipy.constants.hbar**2*(scipy.special.jn_zeros(0,1)/r)**2/(2*mstar)/scipy.constants.elementary_charge-0.067 for r in radius_range])
    ax3.plot(r_occ_data,eval_occ_data,'b-')

    ylim = ax3.get_ylim()
    xlim = ax3.get_xlim()

    ax3.vlines(2.5,-10,eval_occ_data[np.argmin(abs(r_occ_data-2.5))],"r")
    ax3.hlines(eval_occ_data[np.argmin(abs(r_occ_data-2.5))],2,2.5,"r")
    ax3.text(2.5-0.15,eval_occ_data[np.argmin(abs(r_occ_data-2.5))]+0.03,r"$E_{0,r=2.5nm}=%1.2lfeV$" %(eval_occ_data[np.argmin(abs(r_occ_data-2.5))]) )

    ax3.vlines(3.8,-10,eval_occ_data[np.argmin(abs(r_occ_data-3.8))],"r")
    ax3.hlines(eval_occ_data[np.argmin(abs(r_occ_data-3.8))],2,3.8,"r")
    ax3.text(3.8-0.2,eval_occ_data[np.argmin(abs(r_occ_data-3.8))]+0.02,r"$E_{0,r=3.8nm}=%1.2lfeV$" %(eval_occ_data[np.argmin(abs(r_occ_data-3.8))]) )

    ax3.vlines(4.5,-10,eval_occ_data[np.argmin(abs(r_occ_data-4.5))],"r")
    ax3.hlines(eval_occ_data[np.argmin(abs(r_occ_data-4.5))],2,4.5,"r")
    ax3.text(4.5-0.2,eval_occ_data[np.argmin(abs(r_occ_data-4.5))]+0.02,r"$E_{0,r=4.5nm}=%1.2lfeV$" %(eval_occ_data[np.argmin(abs(r_occ_data-4.5))]) )

    ax3.set_ylim(ylim)
    ax3.set_xlim([2,5])

    ax3.set_xlabel(r"Circular corral radius ($nm$)")
    ax3.set_ylabel(r"First eigenmode energy ($eV$)")

    ax3.text(0.025, 0.95,"(a)", transform=fig.transFigure)
    ax3.text(0.025, 0.525,"(b)", transform=fig.transFigure)

    ax3.legend(["Unoccupied corrals", "Occupied corrals","Radii of interest"])

    # fig.suptitle("Circular corral ground state from particle in a box model\nUsing $m^*=$%1.2lf$m_e$, $E_0$=-67 mV" %(mstar/scipy.constants.electron_mass))
    plt.tight_layout()
    fig.subplots_adjust(hspace=0)
    # plt.savefig("/Users   /akipnis/Desktop/n-0_energies.pdf")
