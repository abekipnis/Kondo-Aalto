# %% Imports

from AaltoAtoms import CircCorralData, Spec
import os
from importlib import reload
import data_array
reload(data_array)
from data_array import H2Pc_corrals_edge, H2Pc_corrals_center
from AaltoAtoms.utils.interpolate_dIdV_v_r import imshow_dIdV_vs_r
import numpy as np
import matplotlib.pyplot as plt

# %% Load data
basepath = r"Y:\labdata\Createc\STMDATA\Ag(111)\2022-05-16 Pc depositions"
H2Pc_data = []
#for p in H2Pc_corrals_center:
for p in H2Pc_corrals_edge:

    S = Spec(os.path.join(basepath, p.vertfile))
    args = [os.path.join(basepath, p.datfile), p.datfile, None]
    C = CircCorralData(*args)
    kwargs = {'percentile':p.height_percentile, 'edge_cutoff':p.edge_cutoff}
    C.get_region_centroids(**kwargs)

    C.occupied = True
    C.corral = True
    radius = C.get_corral_radius(1.5, savefig=False, showfig=True)

    width = None

    dIdV, bias_mv = S.dIdV, S.bias_mv

    H2Pc_data.append([radius, width, dIdV, bias_mv, S, C, None, None])

H2Pc_data = sorted(H2Pc_data, key=lambda x: x[0])
# %%

for p in H2Pc_data:
    plt.plot(p[3], p[2]/p[2][np.argmin(np.abs(p[3]+79))] + p[0], label="%1.1lf" %(p[0]))
plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\images\H2Pc_corrals_edge-line.png")

plt.legend()
imshow_dIdV_vs_r(H2Pc_data, enforce_conformity=False, interpolate=True, mV_step=0.1, nm_step=0.1, norm_to_one=True, norm_mV=-79, cmap_str="plasma")
# plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\images\H2Pc_corrals_center.png")
# plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\images\H2Pc_corrals_center.pdf")
# plt.savefig(r"C:\Users\kipnisa1\Dropbox\papers-in-progress\Small Kondo corrals\images\H2Pc_corrals_center.svg")


#plt.matshow([c[4].dIdV for c in H2Pc_data], aspect=11)

#[c[3][0]==80 for c in H2Pc_data]
