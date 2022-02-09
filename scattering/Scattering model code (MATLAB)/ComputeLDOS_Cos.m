function [LDOS,LDOS_norm] = ComputeLDOS_Cos(cps, Energies, ls, res)
%COMPUTELDOS Computes the LDOS by approximating Bessel with Cos()^2
%   For reference, see Crommie et al 1993 (Science vol. 262)
%   Inputs:
%       - cps: the pixel coordinates of the scatterers
%       - Energies: energies at which the LDOS are computed (in eV)
%       - ls: side length of the scan area (in m)
%       - res: pixel resolution (square scans)
%   Outputs:
%       - LDOS: the LDOSes calculated at the specified energies, a.u.
%       - LDOS_norm: LDOS normalized with the surface state background

% Initializations
global rho_surf
global hbar
global ms
global E0
global E_F

% Parallelization and globals do not mix
localGlobals = [hbar, ms, E0, E_F];

% Initializing LDOS arrays
LDOS = zeros(res,res, length(Energies));
LDOS_norm = LDOS;

% Initialize Green's function matrices
Natoms = size(cps, 1);

% Generate grid from cps
skipGrid = Cps2Grid(cps, res, ls);

for Ee = 1:length(Energies)
    
    tic;
    % Current energy value
    EC = Energies(Ee);
    disp(strcat('Starting computations at ', num2str(EC), 'eV'))

    % Phase: use as fitting parameter or try different options.
    phase = 85/180*pi;

    % Start parallel loop to compute LDOS at different points
    parfor x = 1:res
        for y = 1:res
            
            % First things first: no calculations too close to scattering points
            if skipGrid(x,y) >= 1 % Index consistency checked
                continue
            end
                
            % LDOS through Crommie's paper: the cos^2 approximation
            k = k_disp(Ee, localGlobals);
            
            % Distances from atoms
            rhoD = pdist([x y; cps(:,1:2)]);
            rho = rhoD(1:Natoms).*ls./res;
            LDOS(x,y, Ee) = sum(1./(k.*rho).*(cos(k.*rho - pi/4 + phase).^2 ...
                - cos(k.*rho - pi/4)));
            % Sum dimensions? A consistent implementation would use only
            % functions or only matrices, so .* should be best.
        end
        disp(strcat('LDOS computed at x=', num2str(x)));
        
    end
    toc;
    
    % LDOS normalisation to surface state DOS
    LDOS_norm(:,:,Ee) = LDOS(:,:,Ee);%./rho_surf;
end

end

% UTILITY FUNCTIONS: G0 and k_disp cannot contain global variables

% k-value from surface state dispersion
function [k] = k_disp(E, localGlobals)
    hbar = localGlobals(1);
    ms = localGlobals(2);
    E0 = localGlobals(3);
    E_F = localGlobals(4);
    
    % Eq 1: E0 should be negative
    k = abs(sqrt((E - E_F - E0).*2*ms))/hbar;
end


% s and phase_Kondo are outside parfor loop, so global variables are fine

% Phase shift from Kondo impurity: Eq 27
function ph_K = phase_Kondo(E)
    global delta_bg
    global delta_primes
    global Gamma
    global eps0
    
    ph_K = delta_bg + 1i*delta_primes + atan((E - eps0)./(Gamma/2));
end


