function ComputePointSpectrum(a,b, results,model, point, options)
    arguments
        a double
        b double
        results 
        model 
        point (1,2) double
        
%         options.LineWidth double = 0.05;
        options.EBroad double = 10e-3;
        options.ERange (1,2) double = [-85e-3, 115e-3];
        options.NE double = 300;
    end
    
%COMPUTELINESPECTRA Generates a simulated line spectrum from PDE solutions
%   Mandatory inputs:
%       - a: major axis length, nm
%       - b: minor axis length, nm
%       - results: an EigenResults struct from solvepdeeig
%       - model: the PDEModel struct corresponding to results
%       - point: Where to compute the point spectrum
%   Optional inputs as name-value pairs:
%       - EBroad: Broadening of Lorentzian peaks, eV [double]
%       - ERange: Energy range for line spectrum, eV [(1,2) double]
%       - NE: Number of energy points in the line spectrum [double]
%       - XLim: Input for xlim in eV

% Initialise stuff
global ms
global E0
global m_e

% ecc = sqrt(1-(b/a)^2);
mesh = model.Mesh;
u = results.Eigenvectors;
Energies = results.Eigenvalues;

% Get the mesh point closest to asked point
LPoint = findNodes(mesh, 'nearest', point');

% Sanity check plot
figure;
pdemesh(model);
hold on
plot(mesh.Nodes(1,LPoint), mesh.Nodes(2,LPoint), 'or', 'MarkerFaceColor','g')

% % Sort the x-points
% [~, sortInd] = sort(mesh.Nodes(1,LPoint));
% LPoint = LPoint(sortInd);

% Eigenfunction values on the line
EigAtPoint = u(LPoint,:);

% Broaden and square the wave functions

Broadening = @(E, Ec, width) 1./(pi.*width).*...
    (width.^2./((E - Ec).^2 + width.^2));

LDOSPoint = @(E) EigAtPoint.^2*Broadening(E, Energies + E0, options.EBroad);

% Plot the line spectrum somehow
Ens = linspace(options.ERange(1), options.ERange(2), options.NE);
% figure;plot(Ens.*1e3, LDOSLine(Ens, 35:40))

% MALine = mesh.Nodes(1,LPoint);

figure;
% imagesc([min(MALine) max(MALine)], [Ens(1), Ens(end)].*1e3, ...
%     LDOSPoint(Ens)')
plot(Ens.*1e3, LDOSPoint(Ens))
% set(gca, 'YDir', 'normal')
axis tight
xtickformat('%.2f')
% xticklabels({'Focus', 'Center', 'Focus'})
title(['m^* = ', num2str(ms./m_e, 3), 'm_e'])
% xticks([-a*ecc, 0, a*ecc])
ylabel('LDOS, a.u.')
xlabel('E, meV');

if isfield(options, 'XLim')
    xlim(options.XLim);
end

set(gca, 'FontSize', 14)


end

