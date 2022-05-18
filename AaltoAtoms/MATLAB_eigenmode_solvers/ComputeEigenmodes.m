function [results, model] = ComputeEigenmodes(a,b, options)
    arguments
        a double
        b double
        options.focusAtom (1,1) logical = true;
        options.atomPotential double = 0.9;
        options.atomRadius double = 0.6;
        options.energyRange double = 115e-3;
        options.HMax double = 0.3;
        options.plotAll (1,1) logical = true;
        options.plotEnergyTol (1,1) double = 40e-3;
        options.savePlots (1,1) logical = false;
        options.plotPath string = '';
    end

%COMPUTEEIGENMODES Solves the hard-walled ellipse eigenmodes by FEM
%   Necessary inputs:
%       - a: major axis length, nm
%       - b: minor axis length, nm
%   Optional inputs in name-value pairs:
%       - focusAtom: adds an atom into the focus point [Logical]
%       - atomPotential: Potential of the focus atom [double], eV
%       - atomRadius: Radius of the focus atom [double], nm
%       - energyRange: energy window for solutions [double], eV
%       - HMax: FEM mesh refinement parameter [double], nm
%       - plotAll: flag for plotting eigenmodes, eigenenergies, mesh
%       - plotEnergyTol: window around E_F for plotting eigenmodes
%       - savePlots: flag for saving plots
%       - plotPath: path to folder for saving figures [string]
%   Outputs:
%       - results: struct containing eigenmodes and eigenenergies 
%       - model: the FEM model created by the PDE toolbox

global hbar
global ms
global E0

% Let's compute eccentricity as well
ecc = sqrt(1-(b/a)^2);

% Schrödinger constant in natural units
d = 2*ms/hbar^2*(1e-9)^2;

% plot output folder initialisation
if options.savePlots
    foldern = options.plotPath;
    mkdir(foldern);
end

% Initialise PDE model
model = createpde();

% very strange geometry initialization with pdetools
xc = 0;
yc = 0;
angEllipse = 0;

Ellipse = [4 xc yc a b angEllipse]';

% An occupied focus point
focii = [ecc*a, 0; -ecc*a, 0];

% dl = decsg(Ellipse);
if options.focusAtom
    xc = focii(1,1);
    yc = focii(1,2);

    atomCircle = [1 xc yc options.atomRadius 0 0]';
    gd = [Ellipse, atomCircle];
    %... or, alternatively
    dl = decsg(gd);
else
    dl = decsg(Ellipse);
end

% Sanity check geometry plot
if options.plotAll
    figure;
    pdegplot(dl, 'EdgeLabels','on','FaceLabels','on')
end

geometryFromEdges(model, dl);

%% Boundary condition: Dirichlet (h*u = r)
applyBoundaryCondition(model,'dirichlet',...
    'Edge',1:model.Geometry.NumEdges,...
    'u',0);

% Equation coefficients: should adjust c and/or d
% Equation form for eigenvalues: -nabla*(c*nabla u) + au = lambda*d*u
if options.focusAtom
    specifyCoefficients(model, 'm',0,...
                           'd',d,...
                           'c',1,...
                           'a',0,...
                           'f',0,...
                           'Face', 2);
                       
    specifyCoefficients(model, 'm',0,...
                           'd',d,...
                           'c',1,...
                           'a',d*options.atomPotential,...
                           'f',0,...
                           'Face', 1);  
else
    specifyCoefficients(model, 'm',0,...
                           'd',d,...
                           'c',1,...
                           'a',0,...
                           'f',0);
end


% Eigenvalue range: [E0, slightly above E_F] in eV
r = [0, abs(E0)+ options.energyRange];

% Generate meshes
generateMesh(model,'Hmax',options.HMax);

if options.plotAll
    figure;
    pdeplot(model);
end

%% Solve the eigenvalue pde
% The evalc supresses consol output
[~, results] = evalc('solvepdeeig(model, r)');


%% Whole lotta plots
if options.plotAll
    Etol = options.plotEnergyTol;
    u = results.Eigenvectors;
    Energies = results.Eigenvalues;
    
    indF = find(abs(Energies + E0) < Etol);
    
    % Eigenenergy plots
    figure;
    plot((Energies + E0).*1e3, 'o', ...
        'color', 'b',...
        'MarkerFaceColor', 'b');
    yline(0, '--')
    ylim((r+E0).*1e3)
    xlabel('Eigenmode index')
    ylabel('E, meV');
    set(gca, 'FontSize', 14)
    
    
    
    % Eigenmode plots
    for jj = 1:length(indF)
        figure;
        pdeplot(model, 'XYData', u(:,indF(jj)))
        hold on
        plot(focii(:,1), focii(:,2), 'ko', 'MarkerSize', 30)
        axis equal
        colormap jet
        c = colorbar;
        set(c, 'YTickLabel', [])
        t = ylabel(c, {'\Psi'});
        set(t, 'Rotation', 0, 'HorizontalAlignment', 'left')
        title(['E = ', num2str((Energies(indF(jj)) + E0)*1e3), ' meV'])
        xlabel('x, nm')
        ylabel('y, nm')
        set(gca, 'FontSize', 14)
        
        if options.savePlots
            filen = ['Eigenmode_vac_', num2str(indF(jj)), ...
                '_a_', strrep(num2str(a, 3), '.', 'p'), 'nm' ...
                '_b_', strrep(num2str(b, 3), '.', 'p'), ...
                'nm.pdf'];
            exportgraphics(gcf, [pwd '\' foldern '\' filen])
        end
    end

end

