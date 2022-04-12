import numpy as np
import matplotlib.pyplot as plt
from pdb import set_trace
import scipy.signal
import pandas as pd
import matplotlib.animation as animation
# from matplotlib.animation import FuncAnimation
import importlib

import createc
import pdb
import os
import numpy.ma as ma
from matplotlib.widgets import Slider, Button
from datetime import datetime
from multiprocessing import Pool, freeze_support
from AaltoAtoms.find_atom_positions import CircCorralData
from AaltoAtoms.grid_analysis.animate_grid import Grid
from AaltoAtoms.Kondo_data_analysis import read_vertfile

if __name__ == "__main__":
    # read size of image from .specgrid.dat file
    filename = r"Y:\labdata\Createc\STMDATA\Ag(111)\2022-03 Co Kondo corrals\04-06 6nm Ag walls\temp.specgrid"
    g = Grid(filename)
    t = g.get_topo()
    plt.imshow(t)
    plt.show()
    g.correlate_topo_to_ZBC()

    g.animate_cube(plotpoints=[[12.6,0.05],[13.28,4.69],[5,14.6], [8.93,2.38], [11.1, 12.57]], title="Ag island on Nb110")


    # correlate 0 bias conductance to height
    plt.show()

    pdb.set_trace()
