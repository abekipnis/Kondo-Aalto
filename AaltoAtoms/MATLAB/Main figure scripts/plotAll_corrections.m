clear all
close all
clc

%% Import all the spectra

% Load data
load("CorralData_Abe_SLO.mat"); 

% Add import functions to path
addpath('C:\Users\aaprom1\Downloads\Temporary Matlab\Matlab\Import functions\Plot VERT with nm')

% Define a path to data folders
dataPath = 'P:/asp/labdata/Createc/STMDATA/Ag(111)/2022-05-16 Pc depositions/';

% Sort the H2Pc data based on corral radii
[H2Pc_radii, H2Pc_ind] = sort(H2Pc_corrals_radii);

% % Adjust the edge and center indices (sorted by corral radius)
% newEdgeIndices = [0,11,23,47,26,16,23,14,17,15,19,18,17,19,47];
% newCenterIndices = [0,17,30,40,30,20,19,18,20,19,25,24,24,25,52];

% Load the corrected edge ('split-lobe') and center spectra
load("H2Pc_edge_center_newIndices.mat", "H2Pc_corralsCenter_vertfiles_newInds",...
    "H2Pc_corralsEdge_vertfiles_newInds")

EdgeVerts = H2Pc_corralsEdge_vertfiles_newInds;
CenterVerts = H2Pc_corralsCenter_vertfiles_newInds;

% Filter out the 3.6nm , 6.1nm and 7.4nm corrals
exInds = [1, 7, 10];
H2Pc_radii(exInds) = [];
EdgeVerts(exInds) = [];
CenterVerts(exInds) = [];
newEdgeIndices(exInds) = [];
newCenterIndices(exInds) = [];

% Number of corrals remaining
NP = size(H2Pc_radii,1);

% The save flag and plot folder initialisation
savePlots = false;
% mkdir('Plots')

% Point-wise biases to get everything consistent (besides lock-in and FD broadening)
% Edge and center data have the same biases
Biases = cell(1,NP);
EdgedIdV = cell(1,NP);
CenterdIdV = cell(1,NP);

for kk = 1:NP
    EdgeData = CustomVERTImport('files',strtrim(EdgeVerts{kk}), ...
        'path',dataPath, 'normalisation', 'avg',...
        'correctBiasOffset',true);
    Biases{kk} = EdgeData.Bias;
    EdgedIdV{kk} = EdgeData.dIdV;

    CenterData = CustomVERTImport('files',strtrim(CenterVerts{kk}), ...
        'path',dataPath, 'normalisation', 'avg',...
        'correctBiasOffset',true);
    CenterdIdV{kk} =  CenterData.dIdV;

end

% % Save data for easier access later
% save('H2Pc_edge_center_newIndices.mat', 'EdgeData', 'CenterData', ...
%    'EdgedIdV', 'CenterdIdV',...
%     'Biases', 'dataPath', 'newCenterIndices', 'newEdgeIndices', 'H2Pc_corralsEdge_vertfiles_newInds',...
%     'H2Pc_corralsCenter_vertfiles_newInds');

%% Edge data plot

% Opens a figure for some reason
fullColors = colormap('cool');
cLabel = 'Corral radius, nm';

% Displacement between spectra at varying radii
radDisp = 0.8;

% Normalisation for generating ticks to plot
normRadii = (H2Pc_radii - min(H2Pc_radii))./(max(H2Pc_radii)-min(H2Pc_radii));

h(1) = figure;
colormap cool
hold on
for ii  = 1:NP
    plot(Biases{ii}, ...
        EdgedIdV{ii}+...+ ii*radDisp,...%
        H2Pc_radii(ii)*radDisp, ...
        '-', ...
        "Color", fullColors(round(ii*256/NP), :),...
        'LineWidth',1.5)
end
% hold on
% a = ylim; 
axis tight
hold on
plot([0,0], ylim, 'k--')
title('Edge spectra')
xlabel('Bias, mV')
ylabel('dI/dV, a.u.')

c = colorbar('Ticks', normRadii, ...
    'TickLabels', arrayfun(@(x) sprintf('%.2f',x),H2Pc_radii,'un',0),...
    'FontSize',12);
c.Label.String = cLabel;
c.Label.VerticalAlignment = 'top';
c.Label.FontSize = 14;
xlim([-80,80])
set(gcf, 'Position', [680   240   674   758])
set(gca, 'FontSize', 14)

% Save if flag is set
if savePlots
    exportgraphics(gca, 'Plots/EdgeSpectra_corrected.pdf', ...
        'ContentType', 'vector')
end

%% Center data plot
h(2) = figure;
colormap cool
hold on
for ii  = 1:NP
    plot(Biases{ii}, ...
        CenterdIdV{ii} +...+ ii*radDisp,...%
        H2Pc_radii(ii)*radDisp, ...
        '-', ...
        "Color", fullColors(round(ii*256/NP), :),...
        'LineWidth',1.5)
end
% hold on
% a = ylim; 
axis tight
hold on
plot([0,0], ylim, 'k--')
title('Center spectra')
% h.CurrentAxes.LineWidth = 1.5;
% lines = findobj(gcf,'Type','Line');
% for i = 1:numel(lines)
%   lines(i).LineWidth = 1.5;
% end
xlabel('Bias, mV')
ylabel('dI/dV, a.u.')
% pause
% end
%     title(['r: ', num2str(H2Pc_radii(ii))])
%     pause;
c = colorbar('Ticks', normRadii, ...
    'TickLabels', arrayfun(@(x) sprintf('%.2f',x),H2Pc_radii,'un',0),...
    'FontSize',12);
%             c = colorbar('Ticks', linspace(0,1, NP), 'TickLabels', linspace(1,NP, NP));
c.Label.String = cLabel;
c.Label.VerticalAlignment = 'top';
c.Label.FontSize = 14;
xlim([-80,80])
set(gcf, 'Position', [680   240   674   758])
set(gca, 'FontSize', 14)

if savePlots
    exportgraphics(gca, 'Plots/CenterSpectra_corrected.pdf', ...
        'ContentType', 'vector')
end


%% Edge spectra minus the center spectra

fullColors = colormap('cool');
% fullColors = colormap('winter');

TL = [min(H2Pc_radii), max(H2Pc_radii)];
cLabel = 'Corral radius, nm';
% radDisp = 0.4;
normRadii = (H2Pc_radii - min(H2Pc_radii))./(max(H2Pc_radii)-min(H2Pc_radii));


h(3) = figure;
colormap cool
hold on
for ii  = 1:NP
    plot(Biases{ii}, ...
        EdgedIdV{ii} - CenterdIdV{ii} +...+ ii*radDisp,...
        H2Pc_radii(ii)*radDisp, ...
        '-', ...
        "Color", fullColors(round(ii*256/NP), :),...
        'LineWidth',1.5)
end
% hold on
% a = ylim; 
axis tight
plot([0,0], ylim, 'k--')
% title(['Edge - center, r = ', num2str(H2Pc_radii(ii)), 'nm'])
% hold off
xlabel('Bias, mV')
ylabel('dI/dV, a.u.')
% pause
% end
title('Edge - Center')
%     title(['r: ', num2str(H2Pc_radii(ii))])
%     pause;
c = colorbar('Ticks', normRadii, ...
    'TickLabels', arrayfun(@(x) sprintf('%.2f',x),H2Pc_radii,'un',0),...
    'FontSize',12);
%             c = colorbar('Ticks', linspace(0,1, NP), 'TickLabels', linspace(1,NP, NP));
c.Label.String = cLabel;
c.Label.VerticalAlignment = 'top';
c.Label.FontSize = 14;
set(gca, 'FontSize', 14)
xlim([-80,80])
set(gcf, 'Position', [680   240   674   758])

if savePlots
    exportgraphics(gca, 'Plots/EdgeMinusCenter_corrected.pdf', 'ContentType','vector')
end


%% Edge spectra divided by center spectra

h(4) = figure;
colormap cool
hold on
for ii  = 1:NP
    plot(Biases{ii}, ...
        EdgedIdV{ii}./CenterdIdV{ii}+...+ ii*radDisp,...
        H2Pc_radii(ii)*radDisp, ...
        '-', ...
        "Color", fullColors(round(ii*256/NP), :),...
        'LineWidth',1.5)
end
axis tight
hold on;plot([0,0], ylim, 'k--')
xlabel('Bias, mV')
ylabel('dI/dV, a.u.')
title('Edge/Center')
%     title(['r: ', num2str(H2Pc_radii(ii))])
%     pause;
% c = colorbar('Ticks', [0 1], 'TickLabels', arrayfun(@(x) sprintf('%.2f',x),TL,'un',0));
c = colorbar('Ticks', normRadii, ...
    'TickLabels', arrayfun(@(x) sprintf('%.2f',x),H2Pc_radii,'un',0),...
    'FontSize',12);
%             c = colorbar('Ticks', linspace(0,1, NP), 'TickLabels', linspace(1,NP, NP));
c.Label.String = cLabel;
c.Label.VerticalAlignment = 'top';
c.Label.FontSize = 14;
set(gca, 'FontSize', 14)
xlim([-80,80])
set(gcf, 'Position', [680   240   674   758])

if savePlots
    exportgraphics(gca, 'Plots/EdgeByCenter_corrected.pdf', ...
        'ContentType', 'vector')
    savefig(h, 'Plots/CenterSubtractions_corrected.fig')
end
