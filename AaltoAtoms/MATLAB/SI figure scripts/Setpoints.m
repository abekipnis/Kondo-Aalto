clear all
close all
clc

%%%% Setpoints/initial currents as a function of corral radius

% Set paths to data import and eigenmode calculation functions
addpath('C:\Users\aaprom1\Downloads\Temporary Matlab\Matlab\Import functions\Plot VERT with nm')
addpath('C:\Users\aaprom1\Downloads\Temporary Matlab\Matlab\Koralleja ja heijastuksia\Eigenmode approach')

% Load data files
load("CorralData_Abe_SLO.mat", "H2Pc_corralsEdge_datfiles", ...
    "H2Pc_corrals_radii",...
    "H2Pc_corralsEdge_vertfiles",...
    "WallAtoms",...
    "dataPath")

% Sort the data with corral radius
[H2Pc_radii, H2Pc_ind] = sort(H2Pc_corrals_radii);
H2Pc_corralsEdge_vertfiles = H2Pc_corralsEdge_vertfiles(H2Pc_ind);

% Filter out the 3nm corral
H2Pc_radii(1) = [];
H2Pc_corralsEdge_vertfiles(1) = [];
WallAtoms(1) = [];

%% Get all the data

% Nr of radii
NR = length(H2Pc_radii);
% (run from the correct folder or fix later)
lineNames = ls("SLO line spectra fits\*line.mat");

StartEdgeCurrents = zeros(1,NR);
StartCenterCurrents = zeros(1,NR);

% Let's use the bias 78 mV for consistency!
StartBias = zeros(1,NR);

% Loop over radii
for ii = 1:NR
    % Load spectrum fits
    load(['SLO line spectra fits\', lineNames(ii,:)],...
        'files', 'path')
    
    % Import data and correct bias offset
    Data = CustomVERTImport("files",files, "path",path,...
        'correctBiasOffset', true);
    
    % Damn. Indices might be off
    edgeIndex = str2double(H2Pc_corralsEdge_vertfiles{ii}(end-7:end-5));
    centerIndex = ceil(size(Data.I,2)/2); 
    
    % Find the index with bias closest to 78 mV
    [~, BiasInd] = min(abs(Data.Bias - 78));
    StartBias(ii) = Data.Bias(BiasInd);

    StartEdgeCurrents(ii) = Data.I(1, edgeIndex);
    StartCenterCurrents(ii) = Data.I(1, centerIndex);

    % One sanity check plot
    if ii == 5
        figure;plot(Data.Bias, Data.I, 'o')
    end
    
end

%% Plot currents at 78mV

% Currents at 78mV
hPo(1) = figure;
scatter(H2Pc_radii, StartEdgeCurrents.*1e9, 'MarkerEdgeColor','none', 'MarkerFaceColor','flat')
hold on
scatter(H2Pc_radii, StartCenterCurrents.*1e9, 'MarkerEdgeColor','none', 'MarkerFaceColor','flat')
xlabel('Corral radius (nm)')
ylabel('Current at 78 mV (nA)')
% title('Molecule edge spectra')
legend('Edge', 'Center')
axis tight
set(gca, 'FontSize', 12)

% Currents at 80 mV 
% figure;
% scatter(H2Pc_radii, StartCenterCurrents.*1e9, 'MarkerEdgeColor','none', 'MarkerFaceColor','flat')
% xlabel('Corral radius (nm)')
% ylabel('Current at 80 mV (nA)')
% title('Molecule center spectra')
% axis tight
% set(gca, 'FontSize', 12)

% Conductance plots at 78 mV
hPo(2) = figure;
scatter(H2Pc_radii, StartEdgeCurrents.*1e9./(StartBias.*1e-3), 'MarkerEdgeColor','none', 'MarkerFaceColor','flat')
hold on
scatter(H2Pc_radii, StartCenterCurrents.*1e9./(StartBias.*1e-3), 'MarkerEdgeColor','none', 'MarkerFaceColor','flat')
xlabel('Corral radius (nm)')
ylabel('Conductance at 78 mV (nS)')
% title('Molecule edge spectra')
legend('Edge', 'Center')
axis tight
set(gca, 'FontSize', 12)

% Save figures
savefig(hPo, 'Plots/Setpoints/Setpoints_r.fig')
