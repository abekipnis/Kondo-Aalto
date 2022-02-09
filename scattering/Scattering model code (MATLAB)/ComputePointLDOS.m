function [LDOS,LDOS_norm] = ComputePointLDOS(cps, Energies, ls, res, point)
%COMPUTELDOS Computes the LDOS by the Green's function method
%   For reference, see Fiete et al 2003 (Rev. Mod. Phys 75)
%   Inputs:
%       - cps: the pixel coordinates of the scatterers
%       - Energies: energies at which the LDOS are computed (in eV)
%       - ls: side length of the scan area (in m)
%       - res: pixel resolution (square scans)
%       - point: pixel point in which to calculate the LDOS
%   Outputs:
%       - LDOS: the LDOSes calculated at the specified energies, a.u.
%       - LDOS_norm: LDOS normalized with the surface state background


global rho_surf
global phase_n
global phase_bd
global hbar
global ms
global E0
global E_F

localGlobals = [hbar, ms, E0, E_F];

% Initializing LDOS vectors
LDOS = zeros(length(Energies), 1);
LDOS_norm = LDOS;

% Initialize Greens function matrices
Natoms = size(cps, 1);
G0Mat = zeros(Natoms,Natoms);

% Matrix A
I = eye(Natoms);

for Ee = 1:length(Energies)
    
    tic;
    % Current energy
    EC = Energies(Ee);
    disp(strcat('Starting computations at ', num2str(EC), 'eV'))

    % Initialize G0Mat and A before the scan point loop
    for col = 1:Natoms 
        for row = (col+1):Natoms
            G0Mat(row, col) = G0(cps(row,1:2), cps(col,1:2), EC, ls, res, localGlobals);
            G0Mat(col, row) = G0Mat(row,col);
        end
    end
    
    % (Initialize s-vector if several types of scatterers have been added)
    sVec = cps(:,3)'.*s(phase_Kondo(EC)) + ~cps(:,3)'.*s(phase_bd);

    % A simple matrix operation should be fast enough for a small matrix
    A = I - sVec*G0Mat;
    % Is the Kronecker delta really well described by an identity matrix?

    % Then, the fun begins. Initialize the G0 vector
    G0Vec = zeros(Natoms,1);

    % TO DO: implement without for loops
    for jj = 1:Natoms
        G0Vec(jj) = G0([point(1) point(2)], cps(jj,1:2), EC, ls, res, localGlobals);
    end

    GVec = A\G0Vec;     % Numbers changed

    LDOS(Ee) = -1/pi*imag(localGlobals(2)/(2*localGlobals(1)^2)*(...
        G0([point(1) point(2)], [point(1) point(2)], EC, ls, res, localGlobals) + ...
        sum(sVec'.*G0Vec.*GVec)));
    % Sum dimensions? A consistent implementation would use only
    % functions or only matrices, so .* should be best.

    disp(strcat('LDOS computed at x=', num2str(point(1)), 'y=', num2str(point(2))));    
    toc;
    
    % LDOS normalisation
    % Simple way: relative to rho_surf
    LDOS(Ee) = LDOS(Ee) + rho_surf;
    LDOS_norm(Ee) = LDOS(Ee)./rho_surf;
end




end


% UTILITY FUNCTIONS:

% Outgoing free Green's function
function [G0] = G0(rp, r, E, ls, res, localGlobals)
    % Outgoing free Green's function: page 938 in Fiete
    G0 = -1i*(...ms/(2*hbar^2)*(...
        besselj(0, k_disp(E, localGlobals).*norm(rp - r).*ls/res) + ...
        1i*bessely(0, k_disp(E, localGlobals).*norm(rp - r).*ls/res));
end

% k-value from surface state dispersion
function [k] = k_disp(E, localGlobals)
    hbar = localGlobals(1);
    ms = localGlobals(2);
    E0 = localGlobals(3);
    E_F = localGlobals(4);
    
    % Eq 1: E0 should be negative
    k = sqrt((E - E0).*2*ms)/hbar;
end

% s-wave scattering: Eq 17
function [s] = s(phase)
%     global hbar    
%     global ms
    
%     s = 4*1i*hbar^2/ms.*(exp(2*1i*phase) - 1);
    s = 4*1i.*(exp(2*1i*phase) - 1);

end

% Phase shift from Kondo impurity: Eq 27
function ph_K = phase_Kondo(E)
    global delta_bg
    global delta_primes
    global Gamma
    global eps0
    
    ph_K = delta_bg + 1i*delta_primes + atan((E - eps0)./(Gamma/2));
end
