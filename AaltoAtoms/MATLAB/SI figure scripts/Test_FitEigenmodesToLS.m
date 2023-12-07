clear all
close all
clc

%%%% FIT EIGENMODES TO LINE SPECTRA DATA

%% Imports and initialisations

% Paths to import functions and eigenmode calculators
addpath('C:\Users\aaprom1\Downloads\Temporary Matlab\Matlab\Import functions\Plot VERT with nm')
addpath('C:\Users\aaprom1\Downloads\Temporary Matlab\Matlab\Koralleja ja heijastuksia\Eigenmode approach')

% Surface state dispersion parameters
InitializeGlobals('Ag')

% Load data files
load("CorralData_Abe_SLO.mat", "H2Pc_corralsEdge_datfiles", ...
    "H2Pc_corrals_radii",...
    "H2Pc_corralsEdge_vertfiles",...
    "WallAtoms",...
    "dataPath")

% Sort the data with corral radius
[H2Pc_radii, H2Pc_ind] = sort(H2Pc_corrals_radii);
H2Pc_corralsEdge_vertfiles = H2Pc_corralsEdge_vertfiles(H2Pc_ind);


%% Data imports
% FIRST RUN OF A FILE: SET THE FOLLOWING 
% Based on H2Pc_*_radii
corralIndex = 1;

folderStr = num2str(H2Pc_radii(corralIndex), 3);

% Index of edge spectrum
edgeIndex = str2double(H2Pc_corralsEdge_vertfiles{corralIndex}(end-7:end-5));

% Get and set the number of wall atoms
NWallAtoms = WallAtoms(corralIndex);

% % Get the line spec files with savepath
% [files, path] = SavePath(['SLO line spectra fits/' folderStr ' nm line']);
% load(['SLO line spectra fits/' folderStr ' nm line.mat'])
 
% % Do the import once and save the data
% [Data, header] = CustomVERTImport('files', files, 'path', path,...
%     'normalisation', 'raw');
% save(['SLO line spectra fits/' folderStr, ' nm line.mat'], "Data", "header", "NWallAtoms", '-append')

% % CONSECUTIVE RUNS: uncomment here to load existing variables
load(['SLO line spectra fits/' folderStr ' nm line.mat'])
% load(['SLO line spectra fits/' folderStr '_nm_line_fit.mat'], 'corralIndex', 'edgeIndex', 'NWallAtoms')

% Initialise parameters
Bias = Data.Bias;
dIdV = Data.dIdV;
LineLength = Data.LineLength;
NP = Data.NP;
NE = Data.NE;

% Line spectrum points
Distance = linspace(-LineLength/2, LineLength/2, NP);

% Bias offset correction
Bias = Bias - BiasOffset('file', files{corralIndex}, 'path', path);

% Make a mesh for the fitting
[X, E] = meshgrid(Distance, Bias);

% Normalise the dIdV data
dIdV = (1./dIdV(NE, :)).*dIdV; % Normalise to the minimum energy
dIdV = dIdV./max(max(dIdV)); % normalise to maximum value

% Sanity check plot
DataNorm = 1;

hEigs(1) = figure;
s = surf(X, E, dIdV.*DataNorm); 
view(2)
s.EdgeColor = 'none';
title('Normalised data')
axis tight

%% Use CV tools to extract corral dimensions (or just set the radius)

a = H2Pc_radii(corralIndex);
b = a;

%% Fit initialization

% (global variables are fine for a single fit)
global ms
global m_e
global E0

msLocal = ms;
E0Local = E0;

% Constant offset for data
offset = dIdV(end,1);

% Average distance between corral atoms, with the center atom accounted for
r0 = 0.8;
WD = WeissDistance(a, NWallAtoms, r0);

% Make a handle to CompLineFit_Gen
LineFitHandle = @(params, Line) CompLineFit_Gen('a',params(1)*a,'b',params(1)*b,...
    'EBroad',params(2),...
    'msOrig',msLocal,...
    'E0Orig',E0Local,...
    'E0Coeff', params(6),...
    'CenterPc',true,...
    'atomPotential',params(3),...
    'atomRadius',r0,...
    'MALine',Line,...
    'ERange',Bias.*1e-3,...
    'DecayFactor',params(4),...
    'WeissDistance',WD,...
    'StepHeight',params(5),...
    'MALineOffset', params(7),...
    'EExtension', 50e-3)...
    + offset;

% Parameters: radius correction, EBroad (mV), atomPotential, 
% decay factor, step height, E0 coefficient, spatial offset
% Initial fit parameters, lower and upper bounds
par0 = [0.96, 0.3, 0.5, 0.3, 0.15, 0.97, 0.1];
lb = [0.9, 0, 0.2, 0, 0, 0.5, -0.5];
ub = [1.2, 5, 0.6, 1.5, 1, 1, 0.5];

% Sanity check the initial parameters with a line spectra output
LineFit0 = LineFitHandle(par0, Distance);

hEigs(2) = figure;
s = surf(X,E,LineFit0');
view(2)
axis tight
s.EdgeColor = 'none';
title('Initial guess')

%% Fit the line spectra

% Limit points out from the center
centerPoints = floor(NP*r0/LineLength);

fitDistance = [Distance(1:floor(NP/2 - centerPoints)), ...
    Distance(ceil(NP/2 + centerPoints):end)];
fitdIdV = [dIdV(:,1:floor(NP/2 - centerPoints)), dIdV(:,ceil(NP/2 + centerPoints):end)];

% Run the fit
parOpt = lsqcurvefit(LineFitHandle, par0, fitDistance, fitdIdV'.*DataNorm, lb, ub);

% Line spectra with the fitted parameters
LineOpt = LineFitHandle(parOpt, Distance)';

%% Plot results

% Data heat map
hEigs(3) = figure;
imagesc([Distance(1) Distance(end)], [Bias(1), Bias(end)], ...
    dIdV)
set(gca, 'YDir', 'normal')
xtickformat('%.2f')
% xticklabels({'Focus', 'Center', 'Focus'})
title('Data')
xticks([Distance(1), 0, Distance(end)])
xlabel('Major axis, nm')
ylabel('E, meV');

c = colorbar;
set(c, 'YTickLabel', [])
ylabel(c, {'dI/dV, a.u.'})
set(gca, 'FontSize', 14)


% Fit heat map
hEigs(4) = figure;
imagesc([Distance(1) Distance(end)], [Bias(1), Bias(end)], ...
    LineOpt)
set(gca, 'YDir', 'normal')
xtickformat('%.2f')
% xticklabels({'Focus', 'Center', 'Focus'})
% title(['r = ', num2str(parOpt(1)*a, 3), 'nm, '...
%     'V_A = ', num2str(parOpt(3), 3), ' eV, ', ...
%     '\Gamma = ', num2str(parOpt(2), 3), ' meV'])
title(['Fit (r = ', folderStr, ' nm)'])
xticks([Distance(1), 0, Distance(end)])
xlabel('Major axis, nm')
ylabel('E, meV');

c = colorbar;
set(c, 'YTickLabel', [])
ylabel(c, {'LDOS, a.u.'})
set(gca, 'FontSize', 14)



%% Point spectra: fits and data

figure;

for jj = 1:NP
    pointIndex = jj;
    PointOpt = LineOpt(:,pointIndex);

    plot(Bias, dIdV(:,pointIndex), 'o')
    hold on
    plot(Bias, PointOpt, '--')
    plot([0,0], ylim, 'b--')
    hold off
    legend('Data', 'Fit', 'Location','southeast')
    xlabel('Bias, mV')
    ylabel('dI/dV, a.u.')
    set(gca, 'FontSize', 14)
    title(['Point ', num2str(pointIndex)])
    pause
end

%% Difference map

% Data minus Fit heat map
hEigs(5) = figure;
imagesc([Distance(1) Distance(end)], [Bias(1), Bias(end)], ...
    dIdV - LineOpt)
set(gca, 'YDir', 'normal')
xtickformat('%.2f')
% xticklabels({'Focus', 'Center', 'Focus'})
title('Data - fit')
xticks([-a, 0, a])
xlabel('Major axis, nm')
ylabel('E, meV');

c = colorbar;
% set(c, 'YTickLabel', [])
ylabel(c, {'LDOS, a.u.'})
set(gca, 'FontSize', 14)

% Data minus fit, edge point
hEigs(6) = figure;
plot(Bias, dIdV(:,edgeIndex) - LineOpt(:,edgeIndex))
hold on;
plot([0,0], ylim, 'b--')
hold off
xlabel('Bias, mV')
ylabel('dI/dV, a.u.')
set(gca, 'FontSize', 14)
title(['Data - fit, point ', num2str(edgeIndex)])

%% Extract and save the eigenenergies

% Use ComputeEigenmodes
[res, ~] = ComputeEigenmodes(a*parOpt(1), a*parOpt(1),...
    'atomPotential', parOpt(3),...
    'atomRadius',r0,...
    'centerAtom',true, ...
    'centerPc',true,...
    'plotAll',false,...
    'E0_local',parOpt(6)*E0Local);
EigEns = res.Eigenvalues;

% Append the eigenenergies to saved data
save(['SLO line spectra fits/' folderStr ' nm line.mat'], 'EigEns', '-append');


%% Save the fit variables (and the script as well!)

% Save the variables
save(['SLO line spectra fits/' folderStr, '_nm_line_fit.mat']);

% Prepare a folder for the results and plots
mkdir(['Plots/Eigenmode fits/', folderStr, 'nm'])

% Save the figures
savefig(hEigs, ['Plots/Eigenmode fits/', folderStr, 'nm/EigenmodeFit_', folderStr, 'nm.fig'])
clear('hEigs', 's')

% % Make a directory for the fit scripts
% mkdir('SLO line spectra fits/Line spectra fit scripts')

% Write the script (adjust the file name as needed for different versions)
% Does not name the script properly unless the it is run in its entirety (F5)
FileNameAndLocation=[mfilename('fullpath')];
currentFile = [FileNameAndLocation '.m'];
backup = sprintf(['%s_', folderStr, 'nm.m'], FileNameAndLocation);
copyfile(currentFile, backup)
movefile(backup, './SLO line spectra fits/Line spectra fit scripts')