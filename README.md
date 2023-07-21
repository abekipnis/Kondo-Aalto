# AaltoAtoms Python analysis package
<img src="https://github.com/abekipnis/Atoms/blob/master/logo.png" alt="drawing" style="width:100px;" align="right"/>
For doing data analysis of Createc STM-AFM spectrum and topography files for measurement of the Kondo resonance on magnetic Co atoms inside circular quantum corrals on Ag(111).

Example notebooks for data analysis and for automated experimentation:

    Data analysis: 
    - line_spectrum_plot
        Display spectroscopic data at various spatial positions

    - measure_corral_radius_from_scan
        Load a topography scan from the instrument software and analyze the topography to extract the corral radius in nm

    - fit_fano_resonance_to_spectrum
        Load a .VERT (spectroscopy) file and fit a Fano resonance using the given initial conditions

    Automated experiments:
    - freeform_manipulation_GUI
        Detect atom positions from a topography scan, click to program manipulation paths, shift-enter to run procedure

    - expand_or_contract_corral
        Detect atom positions from topography scan, fit atom positions to circle, determine size of new

Example usage found in, i.e. measure_corral_radius_from_scan.ipynb

Several scripts within:

animate_grid
- for reading grid spectra and saving movie file
- also takes fit_fano functions from read_vertfile to fit Fano resonance to grid spectra over occupied corrals

find_atom_positions
- contains CircCorralData class, main function that analyzes all .dat files in data inventory Excel sheet

fit_lattice
- script for testing CircCorralData functions (fitting Gaussians to atom positions, fitting circle to corral walls, fitting lattice to atom positions)

play_latfile
- for converting .LAT files into .wav's for auditory discrimination of successful/unsuccessful lateral manipulations.

scattering_model
- for simulating LDOS rho(r, E) given a set of atom positions and a lattice

Kondo data analysis:
read_vertfile
- TKinter tool to select files
- Read a spectrum (or a series of spectra) from a .VERT file output by the Createc software
- Perform a fit of a Fano function over the Kondo resonance around a user-chosen range of values
- Get extracted parameters for the Fano resonance
    - linewidth w
    - the center of the Fano lineshape E0
    - the Fano asymmetry parameter q, which can be compared with other values from literature.

## Installation
Use `pip install .` when inside the main directory to install the packages listed in requirements.txt.
