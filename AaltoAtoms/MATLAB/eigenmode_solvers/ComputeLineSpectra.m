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

    end
    
%COMPUTELINESPECTRA Generates a simulated line spectrum from PDE solutions
%   Mandatory inputs:
%       - a: major axis length, nm
%       - b: minor axis length, nm
%       - results: an EigenResults struct from solvepdeeig
%       - model: the PDEModel struct corresponding to results 
%   Optional inputs as name-value pairs:
%       - EBroad: Broadening of Lorentzian peaks, eV [double]
%       - ERange: Energy range for line spectrum, eV [(1,2) double]
%       - NE: Number of energy points in the line spectrum [double]
%       - XLim: input of xlim function [(1,2) double]
%       - YLim: input of ylim function [(1,2) double]
%   Output:
%       - LineSpec: The line spectrum as a matrix

% Initialise stuff
global ms
global E0
global m_e

ecc = sqrt(1-(b/a)^2);
% mesh = model.Mesh;
% u = results.Eigenvectors;
Energies = results.Eigenvalues;
% % Get the mesh points along the major axis
% LPoints = findNodes(mesh, 'box', [-a, a], ...
%     [-options.LineWidth, options.LineWidth]);
% 
% % Sanity check plot
% figure;
% pdemesh(model);
% hold on
% plot(mesh.Nodes(1,LPoints), mesh.Nodes(2,LPoints), 'or', 'MarkerFaceColor','g')
% 
% % Sort the x-points
% [~, sortInd] = sort(mesh.Nodes(1,LPoints));
% LPoints = LPoints(sortInd);
% 
% % Eigenfunction values on the line
% LineVac = u(LPoints,:);

% Define the major axis points
if isfield(options, 'MALine')
    MALine = options.MALine;
else
    MALine = linspace(-a*ecc, a*ecc, NP);
end

yq = zeros(1,NP);

uLine = interpolateSolution(results, MALine, yq, 1:length(Energies));

% Broaden and square the wave functions
% imaginary width component makes the spectra strange 
% HARDCODED DECAY OF PEAK INTENSITIES WITH ENERGY
Broadening = @(E, Ec, width) (E-options.ERange(1) + 3*options.EBroad).^-1.5./(pi.*width).*...
    (width.^2./((E - Ec).^2 + width.^2));

LDOSLine = @(E, x) uLine(x,:).^2*Broadening(E, Energies + E0, ...
    options.EBroad);


% % Manoharan formulation: matrix dimension issues
% LDOSLine = @(E, x) -1/pi.*imag(...
%     uLine(x,:).^2./(E - Energies + E0 + 1i*options.EBroad));
% Imaginary formulation does not generate decaying peak intensities

% Plot the line spectrum somehow
Ens = linspace(options.ERange(1), options.ERange(2), options.NE);

LineSpec = LDOSLine(Ens, 1:NP);
% figure;plot(Ens.*1e3, LDOSLine(Ens, 35:40))

% Major axis in nm
% MALine = mesh.Nodes(1,LPoints);
if options.PlotLine
    figure;
    imagesc([MALine(1) MALine(end)], [Ens(1), Ens(end)].*1e3, ...
        LineSpec')
    set(gca, 'YDir', 'normal')
    xtickformat('%.2f')
    % xticklabels({'Focus', 'Center', 'Focus'})
    title(['m^* = ', num2str(ms./m_e, 3), 'm_e'])
    xticks([-a*ecc, 0, a*ecc])
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

