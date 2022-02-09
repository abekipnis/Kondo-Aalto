clear all
close all
clc

% Script to compare against figure 12 (dIdV maps) in Fiete

%% Initializations

% Natural constants, surface state dispersion parameters
InitializeGlobals('Cu')

% Geometry:
% Major axis/radius of the ellipse (atom center points)
a = 71.3e-10;
% a = 6.38641e-9;
% Findings: an ångström makes a difference in the patterns

% Resolution: res-by-res [pixels]
res = 256;

% Number of (corral) atoms
Natoms = 36;

% Eccentricity of the corral (0 results in a circle)
ecc = 0.5;

% Possible scaling of the corral with respect to ls
ellipseScale = 0.8;

% Side length of the scan square [m]
ls = 2*a/ellipseScale;

% Scatterer types: 1 for Kondo
scaType = 1;
focType = 1;

% Energy (10 meV as in the paper)
Energies = 10e-3;

% Plots and saving: folder name
foldern = [datestr(now, 'yyyymmdd'), 'Fietest'];

% File name base;
filen1 = ['withFocus -E0 a', num2str(a, '%10.5e')];
filen0 = ['noFocus -E0 a', num2str(a, '%10.5e')];

%% Geometries

cps1 = MakeGeometry(res, Natoms, ecc, ellipseScale, scaType, 1, focType);
cps0 = MakeGeometry(res, Natoms, ecc, ellipseScale, scaType, 0, focType);
 
% % sanity check plots if needed
% figure;
% imshow(Cps2Grid(cps1, res, ls))
% pause

% figure;
% imshow(Cps2Grid(cps0, res, ls))
% pause

%% Computations
% See ComputeLDOS for details

% [~, LDOS_norm1] = ComputeLDOS_Cos(cps1, Energies, ls, res);
% [~, LDOS_norm0] = ComputeLDOS_Cos(cps0, Energies, ls, res);

[~, LDOS_norm1, invError1] = ComputeLDOS(cps1, Energies, ls, res);

[~, LDOS_norm0, invError0] = ComputeLDOS(cps0, Energies, ls, res);


%% Plots

PlotAndSaveLDOS(LDOS_norm0, Energies, foldern, filen0, 0, ls);
PlotAndSaveLDOS(LDOS_norm1, Energies, foldern, filen1, 0, ls);

% Mirage: normalisation as in the paper
LDOS_mirage = abs(LDOS_norm1 - LDOS_norm0);
PlotAndSaveLDOS(LDOS_mirage, Energies, foldern, 'mirage', 0, ls);

% % Errors
% PlotAndSaveLDOS(invError1, Energies, foldern, [filen1, 'error'], 0, ls);
% PlotAndSaveLDOS(invError0, Energies, foldern, [filen0, 'error'], 0, ls);
