function LineSpec = ComputeLineSpectra(a,b,NP,results, options)
    arguments
        a double
        b double
        NP double
        results 
        
        options.MALine 
        options.EBroad double = 10e-3;
        options.ERange (1,2) double = [-85e-3, 115e-3];
        options.NE double = 300;
        options.PlotLine logical = false;
        options.XLim (1,2) double
        options.YLim (1,2) double
        options.Cassini logical = false;
        options.WeissDistance double = 3;
        options.DecayFactor double = 0.5;
        options.ms_local
        options.E0_local
    end
    
%COMPUTELINESPECTRA Generates a simulated line spectrum from PDE solutions
%   Mandatory inputs:
%       - a: major axis length, nm
%       - b: minor axis length, nm
%       - NP: number of xy-points in the line spectra
%       - results: an EigenResults struct from solvepdeeig
%   Optional inputs as name-value pairs:
%       - MALine: Major axis line as a (zero-centered) vector
%       - EBroad: Broadening of confined mode peaks, eV [double]
%       - ERange: Energy range for line spectrum, eV [(1,2) double]
%       - NE: Number of energy points in the line spectrum [double]
%       - PlotLine: logical flag for plotting the line spectrum
%       - XLim: input of xlim function for plotting [(1,2) double]
%       - YLim: input of ylim function for plotting [(1,2) double]
%       - Cassini: flag for using the Cassini oval geometry
%       - WeissDistance: the average distance between atoms in a corral,
%       used in the determination of energy broadening
%       - DecayFactor: Determines an energy-dependent component of confined
%       mode energy broadening.
%       - ms_local: local effective mass of surface state electrons
%       - E0_local: local surface state onset energy, eV
%
%   Output:
%       - LineSpec: The line spectrum as a (NE, NP) -matrix
% TO DO:
%   - WeissDistance should be a flag instead of the numerical value
%   (marginal speedup for not counting it many times, but it's inconvenient
%   as is)

% HARK! Knives have been used to force this to run in parallel loops (the
% global parameters were probably a bad idea to begin with)
% Initialise stuff
ms = options.ms_local;
% if isempty(ms)
%     global ms
% end

E0 = options.E0_local;
% if isempty(E0)
%     global E0
% end
% global m_e
m_e = 5.685600000000000e-12;
hbar = 6.58210000000000e-16;

% Eccentricity for ellipses vs. Cassini ovals (and circles)
if options.Cassini || a == b
    ecc = 1;
else
    ecc = sqrt(1-(b/a)^2);
end

Energies = results.Eigenvalues;

% % Sanity check plot
% figure;
% pdemesh(model);
% hold on
% plot(mesh.Nodes(1,LPoints), mesh.Nodes(2,LPoints), 'or', 'MarkerFaceColor','g')

% Define the major axis points
if isfield(options, 'MALine')
    MALine = options.MALine;
else
    MALine = linspace(-a*ecc, a*ecc, NP);
end

% y-positions along which to interpolate the solutions (all 0)
yq = zeros(1,NP);

% Interpolate the eigenmode solutions along the major axis
if length(Energies) == 1
    uLine = interpolateSolution(results, MALine, yq);
else
    uLine = interpolateSolution(results, MALine, yq, 1:length(Energies));
end

% Broaden and square the wave functions
% % Legacy formulation of a Lorentzian peak
% Broadening = @(E, Ec, width) ...
%     (E-options.ERange(1) + options.DecayFactor*options.EBroad).^options.DecayExponent./(pi.*width).*...
%     (width.^2./((E - Ec).^2 + width.^2));
% % Gaussian broadening with energy-dependent decay
% Broadening = @(E, Ec, width) ...
%     (E-min(options.ERange) + options.DecayFactor*options.EBroad).^options.DecayExponent...
%     ./(width.*sqrt(2*pi)).*exp((-1.*(E - Ec).^2)./(2*width.^2));

% Broadening according to Weiss et al.: Gaussian with an energy-dependent
% width component
Const = options.DecayFactor*hbar*sqrt(2/ms)/(options.WeissDistance*1e-9);

% Linear correction to Weiss et al.
GWidth = @(E) options.EBroad + Const.*sqrt(abs(E - E0));

Broadening = @(E, Ec) ...
    1./(GWidth(E).*sqrt(2*pi)).*exp((-1.*(E - Ec).^2)./(2.*GWidth(E).^2));

% Square and broaden the eigenmodes
LDOSLine = @(E, x) uLine(x,:).^2*Broadening(E, Energies + E0);

% Evaluate the line spectrum function in the given energy window
Ens = linspace(options.ERange(1), options.ERange(2), options.NE);
LineSpec = LDOSLine(Ens, 1:NP);

% Plot the line spectrum if needed
if options.PlotLine
    figure;
    imagesc([MALine(1) MALine(end)], [Ens(1), Ens(end)].*1e3, ...
        LineSpec')
    set(gca, 'YDir', 'normal')
    xtickformat('%.2f')
    % xticklabels({'Focus', 'Center', 'Focus'})
    title(['m^* = ', num2str(ms./m_e, 3), 'm_e'])
    
    % xticks based on eccentricity
    if ecc == 0
        xticks([-a, 0, a])
    else
        xticks([-a*ecc, 0, a*ecc])
    end
    
    xlabel('Major axis, nm')
    ylabel('E, meV');

    if isfield(options, 'XLim')
        xlim(options.XLim);
    end

    if isfield(options, 'YLim')
        ylim(options.YLim);
    end

    c = colorbar;
    set(c, 'YTickLabel', [])
    ylabel(c, {'LDOS, a.u.'})
    set(gca, 'FontSize', 14)
end


end

