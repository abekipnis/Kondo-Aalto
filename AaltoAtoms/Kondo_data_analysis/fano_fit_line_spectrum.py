import numpy as np
import matplotlib.pyplot as plt
from read_vertfile import fano, plot_fano_fit_line

if __name__=='__main__':
    f = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/Kondo data analysis/fano_fit_data_line_spectrum_2p5nm_radius.txt"
    plot_fano_fit_line(f)

# np.arange(1,2,.5)
