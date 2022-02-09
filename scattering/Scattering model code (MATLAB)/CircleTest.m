% clear all
close all
clc

% To test the consistency of the improved code with fig.4 from Fiete
% 27.04.2020

%% Initializations

% Natural constants and surface parameters
InitializeGlobals('Cu')
% InitializeGlobals('Ag')
global E0
% Geometry:
% Major axis, [m]
major = 88.7e-10;
% major = 4.5e-9;
% Possible scaling of the corral with respect to ls
ellipseScale = 0.75;

% Side length of the scan square [m]
ls = 2*major/ellipseScale;

% Resolution: res-by-res [pixels]
res = 256;

% Number of (corral) atoms
Natoms = 60;

% Adatom type: 1 for Kondo, 0 for normal
type = 1;

% Eccentricity of the corral (0 results in a circle)
ecc = 0;

% Add focus or not
addFocus = 0;

% Energy range
E_min = 420e-3 + E0;
E_max = 480e-3 + E0;
N_E = 40;

Energies = linspace(E_min, E_max, N_E);

N_E = length(Energies);
% Plots and saving: folder name
% foldern = '200428_CircleTest_Nanjiang_fig1';
foldern = [datestr(now, 'yyyymmdd'), '_CircleTest_NE', num2str(N_E)];

% File name base;
filen = 'Circle';

% midlines or not
midlines = 0;

%% Geometry
cps = MakeGeometry(res, Natoms, ecc, ellipseScale, type, addFocus, type);

%% Computations

% LDOS at focus point
%     focus = floor([0;ellipseScale*res/2*ecc] + res/2);

[LDOS, LDOS_norm] = ComputeLDOS(cps, Energies, ls, res);

%     LDOS_focus(ii,jj) = LDOS_norm(focus(1), focus(2));
Topo = LDOS2Topo(LDOS,Energies);

%% Plots

PlotAndSaveLDOS(LDOS_norm, Energies, foldern, [filen '_LDOS'], midlines, ls);
PlotAndSaveLDOS(Topo, Energies, foldern, [filen '_topo'], midlines, ls);

% Pick only the values in Fiete
picks = ceil(linspace(1,N_E,4));
LDOS_picks = LDOS(:,:,picks);
Topo_picks = Topo(:,:,picks);
E_picks = Energies(picks);

PlotAndSaveLDOS(Topo_picks, E_picks, foldern, [filen '_topo_picks'], 1, ls);

