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
major = 71.4e-10 + 3e-10;

% Possible scaling of the corral with respect to ls
% (adjust to a smaller value if the code crashes due to Cps2Grid)
ellipseScale = 0.75;

% Side length of the scan square [m]
ls = 2*major/ellipseScale;

% Resolution: res-by-res [pixels]
res = 256;

% Number of (corral) atoms
Natoms = 48;

% Adatom type: 1 for Kondo, 0 for normal
type = 0;

% Eccentricity of the corral (0 results in a circle)
ecc = 0;

% Add focus or not
addFocus = 0;

% Energy range
E_min = E0;
E_max = 300e-3;
N_E = 1500;

Energies = linspace(E_min, E_max, N_E);
N_E = length(Energies);

% Plots and saving: folder name
foldern = [datestr(now, 'yyyymmdd'), ...
    'CrampinCircleTest_nonRoundedCps_NE', num2str(N_E), ...
    'res', num2str(res)...
    '_LD_',...];%,...
    ...'phase_bd'];
    'aPlus3e-10'];

% File name base;
filen = 'Circle';

% midlines or not
midlines = 0;

%% Geometry
cps = MakeGeometry(res, Natoms, ecc, ellipseScale, type, addFocus, type);

% Use lattice discretization (or not)
cps = LatticeDiscretize(cps, res, ls);

% % Plots for sanity checking
% figure;imshow(Cps2Grid(cps, res, ls))
% pause;
% figure;imshow(Cps2Grid(cps, res, ls))
% pause;

%% Computations

% Center point coordinate
center = floor([res/2;res/2]);

% One LDOS map at Fermi energy
[LDOS, LDOS_norm] = ComputeLDOS(cps, 0, ls, res);

% Point LDOS at center point
[~, LDOS_point] = ComputePointLDOS(cps, Energies, ls, res, center);
    

%% Plots

% Fermi energy map
PlotAndSaveLDOS(LDOS_norm, 0, foldern, [filen '_LDOS'], midlines, ls);
% PlotAndSaveLDOS(-LDOS_norm, 0, foldern, [filen '_LDOS_negative'], midlines, ls);

%% Plotting the spectrum
figure;
plot(Energies, LDOS_point, 'k')
xlim([-0.5, 0.3])
xlabel('E - E_{F} (eV)')
ylabel('LDOS, a.u.')
title(['N_E = ', num2str(N_E)]);
set(gca, 'FontSize', 14);

% Filename and saving: 
filenFull = strcat(filen, ...
    'centerSpectrum',...
    'N_E', num2str(N_E), ...
    'res', num2str(res),...
    '.png');
saveas(gcf, [pwd '\' foldern '\' filenFull])
