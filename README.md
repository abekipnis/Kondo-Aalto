# AaltoAtoms Python package
<img src="https://github.com/abekipnis/Atoms/blob/master/etc/logo.png" alt="drawing" style="width:100px;" align="right"/>
For semi-automated atom manipulation experiments and analysis of Createc STM-AFM dI/dV spectrum and topography files.


## Example outputs using functions / notebooks available here:
  - Line spectra over Co adatoms in quantum corrals on Ag(111)
<img src="https://github.com/abekipnis/Atoms/blob/master/etc/Co_Kondo.gif" width="250">

  - Line spectra over H2Pc molecules in quantum corrals on Ag(111)
<img src="https://github.com/abekipnis/Atoms/blob/master/etc/H2Pc_GIF.gif" width="250">

Example notebooks for data analysis and for automated experimentation:


## Data analysis:
    - line_spectrum_plot
        Display spectroscopic data at various spatial positions
        Along with the related topography data

    - measure_corral_radius_from_scan
        Load a topography scan from the instrument software
        Analyze the topography to extract the corral radius in nm

    - fit_fano_resonance_to_spectrum
        Load a .VERT (spectroscopy) file
        Fit a Fano resonance using the given initial conditions and bounds
        Extract width w, center of lineshape E0, asymmetry parameter q, etc.

## Automated experiments:
    - freeform_manipulation_GUI
        Detect atom positions from a topography scan
        Click to program manipulation paths in the order labeled in the image
        Shift-enter to run procedure and move atoms in order of programming
        Stop if there is a suspected tip change to prevent destroying structure

    - expand_or_contract_corral
        Detect atom positions from topography scan
        Fit atom positions to circle
        Input the size of new corral
        Move wall atoms to new locations for larger/smaller size corral


## Miscellaneous scripts:

    animate_grid:
    - for reading grid spectra and saving movie file
    - also takes fit_fano functions from read_vertfile to fit Fano resonance to grid spectra over occupied corrals

    find_atom_positions
    - contains CircCorralData class, main function that analyzes all .dat files

    fit_lattice
    - script for testing CircCorralData functions (fitting Gaussians to atom positions, fitting circle to corral walls, fitting lattice to atom positions)

    play_latfile
    - for converting .LAT files into .wav's for auditory discrimination of successful/unsuccessful lateral manipulations.

    scattering_model
    - for simulating LDOS rho(r, E) given a set of atom positions and a lattice

## Installation
Use `pip install .` or `pip3 install .` when inside the main directory.

## Reproducibility and data availability
To reproduce manuscript figures from https://arxiv.org/abs/2312.03364, check Jupyter notebook `Create Figures.ipynb` and `Phenomenological model fit.ipynb` as well as MATLAB codes under AaltoAtoms/MATLAB/Main figure scripts.

Data for the manuscript can be found at https://zenodo.org/records/10214707. Some manual editing of file paths may be necessary. For analyses done in Python, tabulated data is found in data_array.py which is used for subsequent analysis routines.
