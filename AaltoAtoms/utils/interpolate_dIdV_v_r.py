from AaltoAtoms.utils.particle_in_a_box import get_modes, mstar

def imshow_dIdV_vs_r(dataset: list,
                    downsample: bool=False,
                    interpolate:bool=False,
                    norm_mV: float = -75,
                    mV_step:float = 1,
                    nm_step: float = 0.5,
                    enforce_conformity: bool = True,
                    norm_to_one: bool = False,
                    cmap_str="plasma") -> None:
    """

    """

    from scipy.interpolate import griddata
    import itertools
    from itertools import repeat
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    from matplotlib.ticker import MaxNLocator
    import matplotlib.font_manager
    cm = 1/2.54  # centimeters in inches
    fig = plt.figure(figsize=(8*cm,8*cm))
    #plt.imshow([list(reversed(c[2]/c[2][np.argmin(abs(c[3]-7))])) for c in dataset if len(c[3])==502 and c[3][0]==80 and c[3][-1]==-80], aspect=20, interpolation=None, extent=[-80,80, 2, 10])
    if enforce_conformity:
        conforms = lambda x: len(x[3])==502 and x[3][0]==80 and x[3][-1]==-80
        ds = [c for c in dataset if conforms(c)]
    else:
        ds = dataset
    # some data sets have 1002 points instead of 502 like the majority
    # so we have to downsample them to match the length
    # in order to plot them with the other ones in this 3D plot
    if downsample:
        d2 = [d for d in dataset if len(d[3])==1002 and d[3][0]==80]
        dIdV_downsampled = [np.append(d[2][::2],d[2][-1]) for d in d2]
        bias_downsampled = [np.append(d[3][::2],d[3][-1]) for d in d2]
        norm_vals_downsampled = [c[np.argmin(abs(bias_downsampled[n]-norm_mV))] for n,c in enumerate(dIdV_downsampled)]
        Z_downsampled = [d/norm_vals_downsampled[n] for n, d in enumerate(dIdV_downsampled)]

    norm_vals = [c[2][np.argmin(abs(c[3]-norm_mV))] for c in ds]
    X = ds[0][3] # bias (has to be same for all spectra!)
    Z = [c[2]/norm_vals[n] for n,c in enumerate(ds)] # dI/dV signal
    Y = [c[0] for c in ds] # radii

    if norm_to_one:
        Z = [(z - min(z))/(max(z) - min(z)) for z in Z]


    if downsample and len(Z_downsampled) != 0:
        Y = np.append(Y, [c[0] for c in d2])
        Z = np.append(Z, Z_downsampled, axis=0)

    if interpolate:
        # plot the spectra on this new grid
        xstep = mV_step # "pixel size" in in mV
        ystep = nm_step # "pixel size" in nm
        xnew = np.arange(min(X), max(X)+xstep, xstep)
        ynew = np.arange(min(Y), max(Y)+ystep, ystep)
        Xn, Yn = np.meshgrid(xnew, ynew)

        # center things with bias index
        create_pt = lambda x: list(zip(x[3]-x[4].bias_offset, repeat(x[0])))
        pts = np.array([np.array(create_pt(d)) for d in ds])
        pts = pts.reshape(-1,2)

        Z = np.array(Z)
        levels = MaxNLocator(nbins=15).tick_values(Z.min(), Z.max())
        args = np.array(pts), Z.reshape(-1, 2).flatten(), (Xn, Yn)
        gridZ = griddata(*args, method="linear")

        plt.pcolormesh(Xn, Yn, gridZ, cmap=plt.get_cmap(cmap_str))
    else:
        plt.pcolormesh(X, Y, Z, shading="gouraud")#, norm=LogNorm(), vmin=np.array(Z).min(), vmax=np.array(Z).max())
    r_range = np.arange(min(Y), max(Y), 0.1)

    # to get particle in a box eigenmodes analytically
    e0, e1, e2, e3 = get_modes(mstar, 0.067, r_range*1e-9, 4).T

    plt.plot(e0*1e3, r_range, color="red")
    plt.plot(e1*1e3, r_range, color="red")
    plt.plot(e2*1e3, r_range, color="red")
    plt.plot(e3*1e3, r_range, color="red")

    plt.xlim(-80,80)

    from matplotlib.font_manager import FontProperties
    plt.rcParams['font.sans-serif'] = ['Arial']
    # kwargs = {'family':'sans-serif', 'size':7, 'sans-serif': 'Arial'}
    font = {'family' : 'sans-serif','sans-serif': 'Arial',
        'size'   : 7}

    matplotlib.rc('font', **font)
    # _font = FontProperties(**kwargs)
    xlab = plt.xlabel("Bias (mV)")
    ylab = plt.ylabel("Corral radius (nm)")

    # plt.gca().set_xticks([-60, -40, -20, 0,-20,-40,-60], fontname="Arial")
    # plt.gca().set_yticks([4,6,8,10])
    # plt.gca().set_xticklabels([-60, -40, -20, 0,-20,-40,-60], fontproperties=_font)
    # plt.gca().set_yticklabels([4,6,8,10], fontproperties=_font)
    # xticks = plt.xticks(fontproperties=_font)
    # yticks = plt.yticks(fontproperties=_font)
    # xlab.set_font_properties(_font)
    # ylab.set_font_properties(_font)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$dI/dV$ (a.u.)', rotation=270)
    cbar.ax.set_yticks([])
    #plt.tight_layout()
    # # or pull them from MATLAB code
# f = r"\\home.org.aalto.fi\kipnisa1\data\Documents\AaltoAtoms\AaltoAtoms\MATLAB\eigenmode_solvers\eigs-0.6nm-20potential.txt"
# data = []
# with open(f, "r") as handle:
#     lines = handle.readlines()
#     for line in lines:
#         ld = line.split(',')
#         #ld = [fl(l) for l in ld]
#         data.append(ld)

# def convert_float(val):
#     try:
#         return float(val)
#     except ValueError:
#         return np.nan
#
# data = pd.DataFrame(pd.to_numeric(data, errors ='coerce'))
# plt.plot( data[1].apply(lambda x: convert_float(x))*1e3-67, data[0].astype(float))
# plt.plot( data[3].apply(lambda x: convert_float(x))*1e3-67, data[0].astype(float))
# plt.plot( data[4].apply(lambda x: convert_float(x))*1e3-67, data[0].astype(float))
