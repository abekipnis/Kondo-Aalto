function [results, model] = ComputeEigenmodes(a,b, options)
    arguments
        a double
        b double
        options.geometryType string = 'ellipse';
        options.focusAtom (1,1) logical = false;
        options.centerAtom (1,1) logical = false;
        options.centerPc (1,1) logical = false;
        options.atomPotential double = 0.9;
        options.atomRadius double = 0.6;
        options.energyRange double = 115e-3;
        options.HMax double = 0.3;
        options.Hedge double = 0.05;
        options.plotAll (1,1) logical = false;
        options.plotEnergyTol (1,1) double = 40e-3;
        options.savePlots (1,1) logical = false;
        options.plotPath string = '';
        options.ms_local double = 2.2742e-12; % Default: Ag(111)!!
        options.E0_local double = -0.0670;    % Default: Ag(111)!!
    end

%COMPUTEEIGENMODES Solves the eigenmodes of a hard-walled geometry by FEM
%   Necessary inputs:
%       - a: major axis length, nm
%       - b: minor axis length, nm (in a circle a = b)
%   Optional inputs in name-value pairs:
%       - geometryType: either 'ellipse' or 'Cassini' for a Cassini oval
%       - focusAtom: adds an atom into an ellipse focus point [logical]
%       - centerAtom: adds an atom into the corral center [logical]
%       - centerPc: adds a 'Pc molecule' into the corral center [logical]
%       - atomPotential: Potential of the added atom/molecule [double], eV
%       - atomRadius: Radius of the added atom/molecule [double], nm
%       - energyRange: upper energy limit (E_F = 0) for computing solutions [double], eV
%       - HMax: FEM mesh refinement parameter [double], nm
%       - Hedge: FEM mesh refinement parameter for edge areas [double], nm?
%       - plotAll: flag for plotting eigenmodes, eigenenergies and mesh
%       - plotEnergyTol: window around E_F for plotting eigenmodes
%       [double], eV
%       - savePlots: flag for saving plots
%       - plotPath: path to folder for saving figures [string]
%       - ms_local: effective mass of surface state electrons
%       (default Ag(111))
%       - E0_local: surface state onset energy, eV (default Ag(111))
%
%   Outputs:
%       - results: struct containing eigenmodes and eigenenergies 
%       - model: the FEM model created by the PDE toolbox
% FIXME:
% - Would make sense to determine HMax relative to other dimensions

%% Initialize the geometry
% HARK! Knives have been used to force this to run in parallel loops (the
% global parameters were probably a bad idea to begin with)

ms = options.ms_local;
% if isempty(ms)
%     global ms
% end

E0 = options.E0_local;
% if isempty(E0)
%     global E0
% end

hbar = 6.582100000000000e-16;

% global hbar
% global ms
% global E0

if isempty(ms) && isempty(E0) && isempty(hbar)
    error('Dispersion parameters not set: use InitializeGlobals() to specify the surface!')
end

% Let's compute eccentricity as well
ecc = sqrt(1-(b/a)^2);

% Schrödinger constant in natural units
d = 2*ms/hbar^2*(1e-9)^2;

% Plot output folder initialisation
if options.savePlots
    foldern = options.plotPath;
    mkdir(foldern);
end

% Initialise PDE model
model = createpde();

% Very strange geometry initialization with pdetools
switch options.geometryType
    case 'ellipse'        
        xc = 0;
        yc = 0;
        angEllipse = 0;

        Ellipse = [4 xc yc a b angEllipse]';

        % An occupied focus point
        focii = [ecc*a, 0; -ecc*a, 0];

        % dl = decsg(Ellipse);
        if options.focusAtom || options.centerAtom
            
            % Only put the atom in the focus point: center is the default
            if options.focusAtom
                xc = focii(1,1);
                yc = focii(1,2);
            end
            
            atomCircle = [1 xc yc options.atomRadius 0 0]';
            gd = [Ellipse, atomCircle];
            %... or, alternatively
            dl = decsg(gd);

        elseif options.centerPc
            % Parametrise a Pc molecule? a cross should suffice :D
            rect_hor = [3 4 ...
                -options.atomRadius options.atomRadius options.atomRadius -options.atomRadius ...
                -options.atomRadius/3 -options.atomRadius/3 options.atomRadius/3 options.atomRadius/3]';
            
            rect_vert = [3 4 ...
                -options.atomRadius/3 options.atomRadius/3 options.atomRadius/3 -options.atomRadius/3 ...
                -options.atomRadius -options.atomRadius options.atomRadius options.atomRadius]';
            
            % Pad the Ellipse geometry vector to match the rectangles
            Ellipse = [Ellipse;0;0;0;0];
            gd = [Ellipse, rect_vert, rect_hor];

            ns = char('ellipse', 'rect1', 'rect2'); ns = ns';
            sf = 'ellipse+rect1+rect2';
            
            [dl, bt] = decsg(gd, sf, ns);

            % Remove the edges inside the PC
            [dl,~] = csgdel(dl,bt, [4 7 10 13]);

        else
            dl = decsg(Ellipse);
        end

        g1 = geometryFromEdges(model, dl);
        
        % Sanity check geometry plot
        if options.plotAll
            figure;
            pdegplot(dl, 'EdgeLabels','on','FaceLabels','on')
        end

    case 'Cassini'
        % A Cassini oval geometry
        N = 300;
        thetas = linspace(0, 2*pi*(1-1/N), N);

        % Parametric expression
        M = @(t) 2*b^2*cos(2.*t) + 2*sqrt(-b^4 + a^4 + b^4*cos(2.*t).^2);

        X = @(t) sqrt(M(t)./2).*cos(t);
        Y = @(t) sqrt(M(t)./2).*sin(t);
        
        % Create a mesh
        pg = polyshape(X(thetas), Y(thetas));
        
        % Defining focii
        focii = [a, 0; -a, 0];

        if options.focusAtom || options.centerAtom
            % Adding an atom into the focus point:
            
%             xc = focii(1,1);
%             yc = focii(1,2);
            if options.focusAtom
                circleX = @(t) options.atomRadius.*cos(t) + a;
            else
                circleX = @(t) options.atomRadius.*cos(t);
            end

            circleY = @(t) options.atomRadius.*sin(t);
            
            pg = addboundary(pg, circleX(thetas), circleY(thetas));
                    
            tr = triangulation(pg);
            g1 = model.geometryFromMesh(tr.Points', tr.ConnectivityList');
            
            addFace(g1, 2);

        else
            tr = triangulation(pg);
            model.geometryFromMesh(tr.Points', tr.ConnectivityList');
        end
%         
%         tr = triangulation(pg);
%         
%         model.geometryFromMesh(tr.Points', tr.ConnectivityList');

    otherwise
        disp('Geometry not supported: please implement or use another one')
        return
end


%% Boundary condition: Dirichlet (h*u = r)
% Hard wall approximation: wave functions go to zero on the corral edge
applyBoundaryCondition(model, 'dirichlet',...
    'Edge', 1:model.Geometry.NumEdges,...
    'u', 0);
% Dirichlet boundary condition on the focus atom??

% Generate mesh
% generateMesh(model,'Hmax',options.HMax,...
%     'Hedge', {[1,2,3,4], 0.05, [5,6,7,8], options.HMax});

generateMesh(model,'Hmax',options.HMax);

% Equation coefficients: should adjust c and/or d
% Equation form for eigenvalues: -nabla*(c*nabla u) + au = lambda*d*u

if options.focusAtom || options.centerAtom || options.centerPc

    % First confirm the face IDs
    fID_el = getFaceID(model, [0.99*a,0]);
    
    if options.focusAtom
        fID_atom = getFaceID(model, [focii(1,1), 0]);
    else
        fID_atom = getFaceID(model, [0, 0]);
    end
    % FIXME: doesn't work for very small corrals where origin and focus are
    % overlapping
    
    specifyCoefficients(model, 'm',0,...
                           'd',d,...
                           'c',1,...
                           'a',0,...
                           'f',0,...
                           'Face', fID_el);
                       
    specifyCoefficients(model, 'm',0,...
                           'd',d,...
                           'c',1,...
                           'a',d*options.atomPotential,...
                           'f',0,...
                           'Face', fID_atom);

    % Get the edge IDs of the "guy" in the center or focus point
    eid_guy = faceEdges(g1, fID_atom);
    % eid_else = faceEdges(g1, fID_el);

    % Refine mesh
    generateMesh(model,'Hmax',options.HMax, ...
        'Hedge', {eid_guy, options.Hedge});

else
    specifyCoefficients(model, 'm',0,...
                           'd',d,...
                           'c',1,...
                           'a',0,...
                           'f',0);
end

% Eigenvalue range: [E0, slightly above E_F] in eV
r = [0, abs(E0) + options.energyRange];

% Geometry plot with different areas
if options.plotAll
    figure;
    pdeplot(model);
end

%% Solve the eigenvalue PDE
% The evalc suppresses consol output
[~, results] = evalc('solvepdeeig(model, r)');
% [~, results] = solvepdeeig(model, r);


%% Plot results
if options.plotAll
    Etol = options.plotEnergyTol;
    u = results.Eigenvectors;
    Energies = results.Eigenvalues;
    
    % Index of Fermi energy
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
        
        % Save the plots if needed
        if options.savePlots
            filen = ['Eigenmode_vac_', num2str(indF(jj)), ...
                '_a_', strrep(num2str(a, 3), '.', 'p'), 'nm' ...
                '_b_', strrep(num2str(b, 3), '.', 'p'), ...
                'nm.pdf'];
            exportgraphics(gcf, [pwd '\' foldern '\' filen])
        end
    end

end

