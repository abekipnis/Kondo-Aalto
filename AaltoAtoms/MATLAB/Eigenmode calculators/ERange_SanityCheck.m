clear all
close all
clc

%%% Sanity check the effect of ERange in ComputeEigenmodes

% Initialize parameters
InitializeGlobals('Ag')

r = 10;

global E0
global ms

% First with default energy range
[res, model] = ComputeEigenmodes(r, r, "plotAll",false, "E0_local", E0, "ms_local", ms);

% Then with larger energy range
[resE0, modelE0] = ComputeEigenmodes(r, r,"energyRange", 200e-3,...
    "plotAll",false, "E0_local", E0, "ms_local", ms);

%% The eigenenergies should not be different

diff = res.Eigenvalues - resE0.Eigenvalues(1:length(res.Eigenvalues));
figure;plot(diff, 'o')

% Differences close to eps, maximum error around 2e-15. So all is fine in this part

%% Next step: line spectra from the different results

NP = 100;
NAtoms = 30;
WD = WeissDistance(r, NAtoms, 0);

Line1 = ComputeLineSpectra(r,r, NP, res, 'PlotLine',true, ...
    "E0_local", E0, "ms_local", ms,...
    "EBroad",5e-3,...
    "WeissDistance", WD, ...
    'DecayFactor',0.5);
title('EMax = 115meV')

Line2 = ComputeLineSpectra(r,r, NP, resE0, 'PlotLine',true, ...
    "E0_local", E0, "ms_local", ms,...
    "EBroad", 5e-3,...
    "WeissDistance", WD, ...
    'DecayFactor',0.5);
title('EMax = 200meV')

figure;
histogram((Line2-Line1)/median(median(Line2)))
set(gca, 'YScale', 'log')