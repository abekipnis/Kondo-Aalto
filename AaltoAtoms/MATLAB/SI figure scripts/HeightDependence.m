clear all
close all
clc

%% Add import tools to path 

addpath('C:\Users\aaprom1\Downloads\Temporary Matlab\Matlab\Import functions\Plot VERT with nm')

%% Get paths to data files

% SavePath('Path_heightDependence_full')
load("Path_heightDependence_full.mat")

% Drop duplicates (in decreasing index order)
files(17) = []; %.163259
files(13) = []; %.162625
files(11) = []; %.162239
files(6) = []; %.162017

% Get the z-offsets from headers
zOffsets = GetHeaderParam('ZOffset_A', 'files', files, 'path', path);

%% Quick plots of the VERT files

% Heat map plot
% PlotVertFunction('files', files, 'path', path, ...
%     'normalisation', 'avg',...
%     'multiPlotType', 'heatMap', ...
%     'xAxisType', 'custom',...
%     'customXLimits', [0, max(zOffsets)], ...
%     'customXLabel', 'z-offset (Å)');

% Stack plot
[h,data,~] = PlotVertFunction('files', files, 'path', path, ...
    'normalisation', 'avg',...
    'multiPlotType', 'stack', ...
    'stackDisp', 0.15,...
    'xAxisType', 'custom',...
    'customXLimits', [0, max(zOffsets)], ...
    'customXLabel', 'z-offset (Å)');

close(h)

%% Further analysis

NP = data.NP;
NE = data.NE;
dIdV = data.dIdV;

% Let's get point-wise biases
Biases = zeros(NE,NP);

for kk = 1:NP
    biasData = CustomVERTImport('files',files{kk}, 'path',path, 'normalisation', 'avg',...
        'correctBiasOffset',true);
    Biases(:,kk) = biasData.Bias;
end

%% Plot the stack properly with bias offsets taken into account
hStack = figure;

% Plot color and stack displacement setup
fullColors = colormap('cool');
stackDisp = 0.15;

hold on
for jj = 1:NP
    plot(Biases(:,jj), dIdV(:,jj) + stackDisp*jj, ...
        "Color", fullColors(round(jj*256/NP), :))
end
axis tight
xlabel('Bias (mV)')
ylabel('d{\itI}/d{\itV} (arb. units)')

% TL = [0, max(zOffsets)];
TL = zOffsets;
cLabel = 'Z-offset (Å)';

% Adjust the colorbar
% c = colorbar('Ticks', [0,1], 'TickLabels', arrayfun(@(x) sprintf('%.1f',x),TL,'un',0));
c = colorbar('Ticks', linspace(0,1,NP), 'TickLabels', arrayfun(@(x) sprintf('%.1f',x),TL,'un',0));
c.Label.String = cLabel;
% c.Label.VerticalAlignment = 'bottom';
c.Label.VerticalAlignment = 'top';
xlim([-78, 78])
set(gca, 'FontSize', 14)

% Save the figure
savefig(hStack, 'Plots/Height dependence/HDSpectra_BiasOffsetFixed.fig')


%% Plot the currents

hCurrents(1) = figure;
xlabel('Bias (mV)')
ylabel('Current (nA)')
yline(0, 'k--')
hold on
for ii = 1:NP
    % 
    % if ii>14
    %     plot(data.Bias-1.7, data.I(:,ii))
    % else
    %     plot(data.Bias, data.I(:,ii))
    % end
    plot(Biases(:,ii), data.I(:,ii).*1e9)
    hold on
    % plot([0,0], ylim, 'k--')
    % hold off
    % xlim([-4,4])
    % ylim([-5e-11, 5e-11])
    % pause
end
axis tight
plot([0,0], ylim, 'k--')
set(gca, 'FontSize', 12)

% Starting currents as a function of z offset
% Get the starting conductances at 78 mV (bias offset)
StartCurrents = data.I(2,:);
StartCurrents(1:11) = data.I(7,1:11);

hCurrents(2) = figure;
scatter(zOffsets, StartCurrents.*1e9, 'MarkerEdgeColor','none', 'MarkerFaceColor', 'flat')
xlabel('Z offsets (Å)')
ylabel('Current at 78 mV (nA)')
set(gca, 'FontSize', 12)

% Save current figures
savefig(hCurrents, 'Plots/Height dependence/HDCurrents.fig')


%% Fit Fano lineshapes

Bias = data.Bias.*1e-3; % Back to eV
Biases = Biases.*1e-3;

% Params: a, q, eps0, w, b, c
Fano = @(params, V) params(1).*((params(2) + (V- params(3))./params(4)).^2)...
    ./(1+ ((V- params(3))./params(4)).^2) + params(5).*V + params(6);

% Initial guesses, lower and upper bounds
inits = [0.015, ...
    7.5, ...
    -0.001, ...
    0.013, ...
    1.5, ...
    0.8];

ub = [0.05,...
    10,...
    3e-3,...
    30e-3,...
    5,...
    1.0];

lb = [1e-4,...
    3.5,...
    -6e-3,...
    10e-3,...
    0.2,...
    0];

% Fit energy range
vmin = -0.03;
% [~, indVMin] = min(abs(Bias - vmin));

vmax = 0.03;
% [~, indVMax] = min(abs(Bias - vmax));

fitParams = zeros(6, NP);

% Loop over spectra
for jj = 1:NP
    % Initialise the bias indices
    [~, indVMin] = min(abs(Biases(:,jj) - vmin));
    [~, indVMax] = min(abs(Biases(:,jj) - vmax));
    
    % Perform the fit
    [fitP, resNormFit] = lsqcurvefit(Fano, inits, ...
        Biases(indVMax:indVMin,jj), dIdV(indVMax:indVMin, jj),...
        lb, ub);

    % Store fit parameters
    fitParams(:,jj) = fitP;

    % Plot the data, initial guess and fit result
    hFanoFitsWithData(jj) = figure;
    plot(Biases(indVMax:indVMin,jj).*1e3, dIdV(indVMax:indVMin, jj), 'bo')
    hold on
    plot(Biases(indVMax:indVMin,jj).*1e3, Fano(fitP, Biases(indVMax:indVMin,jj)), 'r--')
    plot(Biases(indVMax:indVMin,jj).*1e3, Fano(inits, Biases(indVMax:indVMin,jj)), 'm--')
    hold off
    title(['Z-offset ', num2str(zOffsets(jj), 3), ' Å'])
    xlabel('Bias (mV)')
    ylabel('dI/dV (arb. units)')
    legend('Data', 'Fit', 'Starting parameters', 'Location','best')
    set(gca, 'FontSize', 12)

    % Return to the stack plot and overlay the fits on the data
    figure(hStack)
    hold on
    hf = plot(Biases(indVMax:indVMin,jj).*1e3, Fano(fitParams(:,jj), Biases(indVMax:indVMin,jj)) + jj*stackDisp, 'r-');
    legend(hf, 'Fano fit')

end

% Save the figures with Fano fits
savefig(hFanoFitsWithData, 'Plots/Height dependence/FanoFitPlotsWithData.fig')
savefig(hStack, 'Plots/Height dependence/HDSpectra_BiasOffsetFixed.fig')

%% Plot Fano widths
hFano(1) = figure;
plot(zOffsets, fitParams(4,:).*1e3, 'o', 'MarkerEdgeColor','none', 'MarkerFaceColor','b')
% axis tight
xlabel('Z offsets (Å)')
ylabel('Fano width (mV)')
title('w')
set(gca, 'FontSize', 12)

%% Plot the eps0s
hFano(2) = figure;
plot(zOffsets, fitParams(3,:).*1e3, 'o', 'MarkerEdgeColor','none', 'MarkerFaceColor','b')
% axis tight
xlabel('Z offsets (Å)')
ylabel('\epsilon_0 (mV)')
title('\epsilon_0')
set(gca, 'FontSize', 12)

%% Plot the q
hFano(3) = figure;
plot(zOffsets, fitParams(2,:), 'o', 'MarkerEdgeColor','none', 'MarkerFaceColor','b')
% axis tight
xlabel('Z offsets (Å)')
ylabel('q')
title('q')
set(gca, 'FontSize', 12)

% Save the Fano parameter plots
savefig(hFano, 'Plots/Height dependence/FanoFitsToHD_parOpt.fig')
