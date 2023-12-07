clear all
close all
clc

%%%% Calculate DOS at Fermi in the corral center as a function of corral 
% radius 
% Eigenmode calculation parameters from fits to Co-Co corral line spectra

%% Initialize paths, dispersion parameters and data

% Initialise the paths to import tools and eigenmode calculation functions
addpath('C:\Users\aaprom1\Downloads\Temporary Matlab\Matlab\Import functions\Plot VERT with nm')
addpath('C:\Users\aaprom1\Downloads\Temporary Matlab\Matlab\Koralleja ja heijastuksia\Eigenmode approach')

% Initialize dispersion parameters (global variables converted to locals later)
InitializeGlobals('Ag')
global E0
global ms

% Load and sort data and fit files
load("H2Pcs\CorralData_Abe_SLO.mat", ...
    'Co_Co_radii',... %vector with the corral radii
    'Co_Co_vertfiles'); % cell array with VERT files

[CoCo_radii, CoCo_ind] = sort(Co_Co_radii);
CoCo_verts = Co_Co_vertfiles(CoCo_ind);

% Get the files that actually have been fitted
fitFiles = ls("Co Co eigenmode fits\*_line_fit.mat");
NFits = size(fitFiles, 1);

%% Initialize parameters

% Number of simulated corral radii
NP = 70;

% Corral radius vector
rs = linspace(2.5, 11, NP);

% Numbers of wall atoms
Ns = ceil(linspace(12, 30, NP));

% Central atom radius (if used)
% r0 = 0.25;

% Number of simulated bias points (even number needed to include E_F)
NE = 150;

% Energy range cap ([meV])
EMax = 100;

% Bias vectors
% Bias100 = linspace(-100, EMax, NE+1);
Bias150 = linspace(-100, EMax +50, NE+1);

% DOS(E_F) vectors
% DosAtFermi100 = zeros(NP,1);
DosAtFermi150 = zeros(NP,NFits);
DosAtFermiOcc = zeros(NP,NFits);
DosAtFermiOccNorm = DosAtFermiOcc;

% Index of Fermi energy in the bias vectors
% [~, indEF100] = min(abs(Bias100));
[~, indEF150] = min(abs(Bias150));

% Number of simulated line spectrum points along the radius
NP_spec = 60;

% Average fit parameters
parOptAvg = zeros(1,7);

% Loop over the fit files to extract average fit parameters
for jj = 1:NFits
    
    % Load the fit parameters
    load(['Co Co eigenmode fits/' strtrim(fitFiles(jj,:))],...
        'parOpt', 'r0');
    fixedParams = parOpt;

    % Central atom potential
    % V0 = 1.0;
    V0 = parOpt(3);

    % Rolling average for optimal parameters
    if jj ==1
        parOptAvg = parOpt;
    else
        parOptAvg = (parOptAvg + parOpt)./2;
    end

    % "Localize" the global parameters
    ms_local = ms;
    E0_local = E0;
    
    % Initialise the line spectrum function
    % Params: r, WD, HMax, CenterAtom
    % FixedParams: r correction, EBroad, (V0), Decay factor, Step height,
    % E0 coefficient, major axis line offset
    % (note the setting for normaliseLine must be false for this to work)
    LineHandle = @(params, Line, Bias) CompLineFit_Gen(...
        'a',fixedParams(1)*params(1),...
        'b',fixedParams(1)*params(1),...
        'EBroad',fixedParams(2),...
        'msOrig',ms_local,...
        'E0Orig',E0_local,...
        'E0Coeff',fixedParams(6),...
        'CenterAtom',params(4),...
        'MALine',Line,...
        'atomPotential',V0,...
        'atomRadius',r0,... % has worked reasonably well with Co
        'ERange',Bias'.*1e-3,...
        'DecayFactor',fixedParams(4),...
        'WeissDistance',params(2),...
        'StepHeight',fixedParams(5),...
        'plotAll', false,...
        'plotLine', false,...
        'HMax', params(3),...
        'MALineOffset', fixedParams(7),...
        'EExtension', 50e-3,...
        'normaliseLine', false);
    
    % Generate the LDOS(E_F, r)-curves
    tic
    parfor ii = 1:length(rs)
        % Initialise points for line spec
        Distance = linspace(-rs(ii), rs(ii), NP_spec +1);
        
        % Calculate the line spectrum
        % Line100 = LineFitHandle([rs(ii), rs(ii), ...
        %     rs(ii)/20], ...
        %     Distance, Bias100);
    
        Line150 = LineHandle([rs(ii), WeissDistance(rs(ii), Ns(ii), r0), ...
            rs(ii)/20, false], ...
            Distance, Bias150);
    
        LineOcc = LineHandle([rs(ii), WeissDistance(rs(ii), Ns(ii), r0), ...
            rs(ii)/20, true], ...
            Distance, Bias150);
        
        % ylim([-100, EMax])
    
        % Extract the DOS at Fermi at corral center
        DosAtFermiOcc(ii,jj) = mean(LineOcc(NP_spec/2+1,indEF150-5:indEF150+5));
        DosAtFermi150(ii,jj) = mean(Line150(NP_spec/2+1,indEF150-5:indEF150+5));
    
    end
    toc
    
    % Plot the DOS as a function of corral radius

    % % Sanity check plot with different energy ranges
    % figure;
    % % plot(rs, DosAtFermi100, '-o')
    % hold on
    % plot(rs, DosAtFermi150, '-o')
    % % title(['Occupied, V_0 = ', num2str(V0), ' eV, r_0 = ', num2str(r0), ' nm'])
    % title('Empty')
    % xlabel('r, nm')
    % ylabel('DOS(E_F), a.u.')
    % legend('Emax 100mV', 'Emax 150mV')
    % set(gca,"FontSize", 12)
    
    % Normalise occupied DOS
    DosAtFermiOccNorm(:,jj) = (DosAtFermiOcc(:,jj)-min(DosAtFermiOcc(:,jj)))./(max(DosAtFermiOcc(:,jj))- min(DosAtFermiOcc(:,jj))).*max(DosAtFermi150(:,jj));
    
    h(jj) = figure;
    plot(rs, DosAtFermi150(:,jj), '-', 'LineWidth', 1.5)
    hold on
    plot(rs, DosAtFermiOccNorm(:,jj), '-', 'LineWidth', 1.5)
    % title(['Occupied, V_0 = ', num2str(V0), ' eV, r_0 = ', num2str(r0), ' nm'])
    title('E_{max} = 150 meV')
    xlabel('r (nm)')
    ylabel('DOS(E_F) (arb. units)')
    legend('Empty', ['V = ', num2str(V0, 2),' eV, r_0 = ', ...
        num2str(r0, 2), ' nm'], 'Location','best')
    set(gca,"FontSize", 12)

    % pause
end

% Save the variables after running the loop for the first time
% save('AllDosesAtFermi_CoCo')

%% Plot all the curves together
% If the prior loop has already been run, import the data
load('AllDosesAtFermi_CoCo.mat')

h(jj+1) = figure;
h_leg(1:NFits) = plot(rs, DosAtFermi150, 'b-', 'LineWidth', 1.5);
hold on
h_leg((NFits + 1):2*NFits) = plot(rs, DosAtFermiOccNorm, 'r-', 'LineWidth', 1.5);
axis tight
title('E_{max} = 150 meV')
xlabel('r (nm)')
ylabel('DOS(E_F) (arb. units)')
legend(h_leg([1, (NFits +1)]), 'Empty', 'With Co', 'Location','best')
set(gca,"FontSize", 12)

%% Generate the LDOS(r)-curve with averaged fit parameters
fixedParams = parOptAvg;

DosAtFermi150Avg = zeros(NP,1);
DosAtFermiOccAvg = zeros(NP,1);

% Central atom potential
% V0 = 1.0;
V0 = parOptAvg(3);
% 
% ms_local = ms;
% E0_local = E0;

% Initialise the line spectra function
% Params: as before
LineHandle = @(params, Line, Bias) CompLineFit_Gen(...
    'a',fixedParams(1)*params(1),...
    'b',fixedParams(1)*params(1),...
    'EBroad',fixedParams(2),...
    'msOrig',ms_local,...
    'E0Orig',E0_local,...
    'E0Coeff',fixedParams(6),...
    'CenterAtom',params(4),...
    'MALine',Line,...
    'atomPotential',V0,...
    'atomRadius',r0,... % has worked reasonably well with Co
    'ERange',Bias'.*1e-3,...
    'DecayFactor',fixedParams(4),...
    'WeissDistance',params(2),...
    'StepHeight',fixedParams(5),...
    'plotAll', false,...
    'plotLine', false,...
    'HMax', params(3),...
    'MALineOffset', fixedParams(7),...
    'EExtension', 50e-3,...
    'normaliseLine', false);

tic
parfor ii = 1:length(rs)
    % Initialise points for line spectra
    Distance = linspace(-rs(ii), rs(ii), NP_spec +1);
    
    % % Calculate the line spectrum
    % Line100 = LineFitHandle([rs(ii), rs(ii), ...
    %     rs(ii)/20], ...
    %     Distance, Bias100);

    Line150 = LineHandle([rs(ii), WeissDistance(rs(ii), Ns(ii), r0), ...
        rs(ii)/20, false], ...
        Distance, Bias150);

    LineOcc = LineHandle([rs(ii), WeissDistance(rs(ii), Ns(ii), r0), ...
        rs(ii)/20, true], ...
        Distance, Bias150);
    
    % Extract the DOS at Fermi at center
    DosAtFermiOccAvg(ii) = mean(LineOcc(NP_spec/2+1,indEF150-5:indEF150+5));
    DosAtFermi150Avg(ii) = mean(Line150(NP_spec/2+1,indEF150-5:indEF150+5));

end
toc

% Normalise occupied DOS
DosAtFermiOccNormAvg = (DosAtFermiOccAvg-min(DosAtFermiOccAvg))./(max(DosAtFermiOccAvg)- min(DosAtFermiOccAvg)).*max(DosAtFermi150Avg);

% Plot the DOS at Fermi
hAvg = figure;
plot(rs, DosAtFermi150Avg, '-', 'LineWidth', 1.5)
hold on
plot(rs, DosAtFermiOccNormAvg, '-', 'LineWidth', 1.5)
axis tight
title('Co-Co: Averaged fit parameters')
xlabel('r (nm)')
ylabel('DOS(E_F) (arb. units)')
legend('Empty', ['V = ', num2str(V0, 2),' eV, r_0 = ', ...
    num2str(r0, 2), ' nm'], 'Location','best')
set(gca,"FontSize", 12)

%% Save parameters and figures

save('AllDosesAtFermi_CoCo.mat', 'DosAtFermi150', 'DosAtFermiOcc', 'DosAtFermiOccNorm', ...
    "DosAtFermi150Avg", "DosAtFermiOccNormAvg", "DosAtFermiOccAvg", '-append')
savefig(h, 'Co Co eigenmode fits/AllDosesAtFermi_CoCo.fig')
exportgraphics(h(jj+1), 'Co Co eigenmode fits/LDOS_center_EF_CoCo.pdf', ...
        'ContentType', 'vector', 'BackgroundColor', 'none')