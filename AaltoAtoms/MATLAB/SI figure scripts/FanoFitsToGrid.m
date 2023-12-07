clear all
close all
clc
format longE

%%% Fit grid spectra with Fano curves and plot statistics
%% Add import tools to path

addpath('C:\LocalUserData\User-data\aaprom1\Matlab\Import functions\Plot VERT with nm')


%% Import data (saved to a .mat file)

% Initialize path
% SavePath('isolatedKondoGridPath')

% Import data
% load('isolatedKondoGridPath.mat')
% 
% [header, V, ~, fullData, dIdV] = ImportSpecGrid('file', files, 'path', path, ...
%     'normalisation', 'median');
% save('IsolatedKondoGrid.mat')

load("IsolatedKondoGrid.mat");


%% Sanity check: plot point spectra
% figure;
% for ii = 1:32
%     for jj = 1:32
%         plot(V, squeeze(dIdV(ii,jj,:)), 'o')
%         xlabel('Bias, V')
%         ylabel('dI/dV, a.u.')
%         title(['x = ' num2str(ii) ', y = ' num2str(jj)])
%         pause
%     end
% end

% Exemplary coordinates
% x12, y 21
% x 14, y 22
% x 23, y 8 
% x 28 y 17
% x 31, 16

%% Fano fitting: initializations

% Scan frame size [nm]
x = [0, 1.9];
y = x;

close all

% Anonymous Fano fit function
% Parameters in order: a, q, eps0, w, b, c
% NOTE: w is HWHM, in paper Gamma0 is FWHM
Fano = @(params, V) params(1).*((params(2) + (V- params(3))./params(4)).^2)...
    ./(1+ ((V- params(3))./params(4)).^2) + params(5).*V + params(6);

% Initial fit parameters and their upper and lower bounds
inits = [0.015, ...
    7.5, ...
    -0.002, ...
    0.018, ...
    1.5, ...
    0.7];

ub = [0.05,...
    10,...
    3e-3,...
    30e-3,...
    5,...
    1.0];

lb = [1e-4,...
    3.5,...
    -5e-3,...
    13e-3,...
    0.2,...
    0];

% Energy window for fitting and corresponding bias vector indices
vmin = -0.03;
[~, indVMin] = min(abs(V - vmin));

vmax = 0.06;
[~, indVMax] = min(abs(V - vmax));

% Sanity check initial parameters with an example point
% x = 12; y = 21;
% 
% figure;
% % plot(V, squeeze(dIdV(12,21,:)), 'o')
% plot(V, squeeze(dIdV(x,y,:)), 'o')
% hold on
% plot(V, Fano(inits, V), 'r-')
% xlabel('Bias, V')
% ylabel('dI/dV, a.u.')

% resNorm = sum((squeeze(dIdV(x,y,indVMax:indVMin)) - Fano(inits, V(indVMax:indVMin))).^2);
% Residual norm threshold for fitting: ~15? 10?
resThresh = 11;

% Fit parameter and residual norm tensors
fitParams = zeros(header.xend, header.yend, length(inits));
resNorms = zeros(header.xend, header.yend);

%% Fano fits: run the grid
% Float nonsense: change the type of V
V = double(V);
figure;

for ii = 1:header.xend
    for jj = 1:header.yend
        % Compute the residual norm for Fano fit with initial guess params
        resNorm = sum((squeeze(dIdV(ii,jj,indVMax:indVMin)) - Fano(inits, V(indVMax:indVMin))).^2);

        % If residual norm is smaller than threshold, proceed to fit
        if resNorm < resThresh
            % Fit the guy with the thing in the range
            [fitPara, resNormFit] = lsqcurvefit(Fano, inits, ...
                V(indVMax:indVMin), squeeze(dIdV(ii,jj,indVMax:indVMin)),...
                lb, ub);

            % Save fit results
            fitParams(ii,jj,:) = fitPara;
            resNorms(ii,jj) = resNormFit;
            
            % % plot if neurotic
            % plot(V, squeeze(dIdV(12,21,:)), 'o')
            % plot(V, squeeze(dIdV(ii,jj,:)), 'o')
            % hold on
            % plot(V, Fano(fitPara, V), 'r-')
            % xlabel('Bias, V')
            % ylabel('dI/dV, a.u.')
            % % pause
            % hold off
        end
    end
end

%% Plot the fit parameters
% Unit adjustments for eps and w (eV to meV)
fitParams(:,:,3) = fitParams(:,:,3).*1e3;
fitParams(:,:,4) = fitParams(:,:,4).*1e3;

% Params: a, q, eps0, w, b, c
parNames = {'a', 'q', 'epsilon_0', 'w', 'b', 'c'};
parUnits = {'a.u.', 'a.u.', 'meV', 'meV', '1/V', 'a.u.'};

% Path to plot folder
plotsave = 'Plots/Isolated grid/';

% Plot the fit results
for kk = 1:length(inits)
    
    % % Subplot approach
    % figure;
    % subplot(1,2,1)
    % imagesc(fitParams(:,:,kk))
    % axis square
    % set(gca, 'YDir', 'normal')
    % subplot(1,2,2)
    % histogram(nonzeros(fitParams(:,:,kk)))
    % xlabel(parUnits{kk})
    % axis square
    % title(parNames{kk})
    currentFitParam = fitParams(:,:,kk);
    
    % Separate frames
    figure;
    % subplot(1,2,1)
    imagesc(x, y, currentFitParam)
    axis square
    set(gca, 'YDir', 'normal')
    xlabel('x, nm')
    ylabel('y, nm')

    c = colorbar;
    c.Label.String = [parNames{kk}, ', ' parUnits{kk}];
    c.Label.VerticalAlignment = "top";
    c.Label.FontSize = 12;
    
    % Find smallest non-zero value for color
    clim_min = min(nonzeros(fitParams(:,:,kk)));
    clim_max = max(nonzeros(fitParams(:,:,kk)));
    clim([clim_min, clim_max]);

    set(gca, 'Fontsize', 12)
    
    % Save the figure
    exportgraphics(gca, [plotsave 'IsolatedGrid_' parNames{kk} '.pdf'], ...
        'ContentType','vector',...
        'BackgroundColor', 'none')
    
    % Plot the fit parameter histogram
    figure;
    histogram(nonzeros(fitParams(:,:,kk)))
    axis square
    xlabel(parUnits{kk})
    ylabel('Counts')
    title(parNames{kk})
    set(gca, 'Fontsize', 12)
    
    % Save the figure
    exportgraphics(gca, [plotsave 'IsolatedGrid_' parNames{kk} '_hist.pdf'], ...
        'ContentType','vector',...
        'BackgroundColor', 'none')
end

%% Plot residual norms

figure;
subplot(1,2,1)
imagesc(resNorms)
set(gca, 'YDir', 'normal')
subplot(1,2,2)
histogram(nonzeros(resNorms))
title('Fit residual norms')

%% Plot the dIdV at zero bias as well

[~, indZB] = min(abs(V));

figure;
imagesc(x,y, dIdV(:,:,indZB))
set(gca, 'YDir', 'normal')
axis square
xlabel('x, nm')
ylabel('y, nm')

c = colorbar;
c.Label.String = 'dI/dV, a.u.';
c.Label.VerticalAlignment = "top";
c.Label.FontSize = 12;

% % Find smallest non-zero value for color
% clim_min = min(nonzeros(fitParams(:,:,kk)));
% clim_max = max(nonzeros(fitParams(:,:,kk)));
% clim([clim_min, clim_max]);

set(gca, 'Fontsize', 12)

% Save the figure
exportgraphics(gca, [plotsave 'IsolatedGrid_zeroBiasMap_crop8.pdf'], ...
    'ContentType','vector',...
    'BackgroundColor', 'none')
