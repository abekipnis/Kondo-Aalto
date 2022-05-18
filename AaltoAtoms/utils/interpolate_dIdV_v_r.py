import numpy as np
import matplotlib.pyplot as plt
def imshow_dIdV_vs_r(dataset: list,
                    downsample: bool=False,
                    interpolate:bool=False,
                    norm_mV: float = -75,
                    mV_step:float = 1,
                    nm_step: float = 0.5) -> None:
    """

    """

    from scipy.interpolate import griddata
    import itertools
    from itertools import repeat

    fig = plt.figure(figsize=(8,8))
    #plt.imshow([list(reversed(c[2]/c[2][np.argmin(abs(c[3]-7))])) for c in dataset if len(c[3])==502 and c[3][0]==80 and c[3][-1]==-80], aspect=20, interpolation=None, extent=[-80,80, 2, 10])
    conforms = lambda x: len(x[3])==502 and x[3][0]==80 and x[3][-1]==-80
    ds = [c for c in dataset if conforms(c)]

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
    X = ds[0][3] # bias
    Z = [c[2]/norm_vals[n] for n,c in enumerate(ds)]
    Y = [c[0] for c in ds] # radii

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
        #pts = list(itertools.product(X,Y))
        gridZ = griddata(np.array(pts), np.array(Z).reshape(-1, 2).flatten(), (Xn, Yn), method="linear")
        plt.pcolormesh(Xn, Yn, gridZ)
    else:
        plt.pcolormesh(X, Y, Z, shading="gouraud")#, norm=LogNorm(), vmin=np.array(Z).min(), vmax=np.array(Z).max())
    plt.xlabel("Bias (mV)")
    plt.ylabel("Corral radius (nm)")
    plt.show()
