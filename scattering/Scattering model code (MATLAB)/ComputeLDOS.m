function [LDOS,LDOS_norm, invError] = ComputeLDOS(cps, Energies, ls, res)
%COMPUTELDOS Computes the LDOS by the Green's function method
%   For reference, see Fiete et al 2003 (Rev. Mod. Phys 75)
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
global phase_bd
global hbar
global ms
global E0
global E_F

% Parallelization and globals do not mix
localGlobals = [hbar, ms, E0, E_F];

% Initializing LDOS arrays
LDOS = zeros(res,res, length(Energies));
LDOS_norm = LDOS;
invError = LDOS;

% Initialize Green's function matrices
Natoms = size(cps, 1);

G0Mat = zeros(Natoms,Natoms);

I = eye(Natoms);

% Generate grid from cps
skipGrid = Cps2Grid(cps, res, ls);

for Ee = 1:length(Energies)
    
    tic;
    
    % Current energy value
    EC = Energies(Ee);
    disp(strcat('Starting computations at ', num2str(EC), 'eV'))

    % Initialize G0Mat, SVec and A before the scan point loop (see utility
    % functions below)
    for col = 1:Natoms 
        % bessely diverges for r=0: no computations along the diagonal
        for row = (col+1):Natoms
            G0Mat(row, col) = G0(cps(row,1:2), cps(col,1:2), ...
                EC, ls, res, localGlobals);
            
            % Symmetry of Bessel functions
            G0Mat(col, row) = G0Mat(row,col);
        end
    end
    
    % Initialize s-vector if several types of scatterers have been added
    sVec = cps(:,3)'.*s(phase_Kondo(EC)) + ~cps(:,3)'.*s(phase_bd);

    % A simple matrix operation should be fast enough for a small matrix
    A = I - sVec*G0Mat;
    % s*G0: prefactors cancel to unity
    
    % Is the Kronecker delta properly described by an identity matrix?
    
    % Start parallel loop to compute LDOS at different points
    parfor x = 1:res
        for y = 1:res
            
            % First things first: no calculations too close to scattering points
            if skipGrid(x,y) >= 1 % Index consistency checked
                continue
            end
                
            % Then, the fun begins. Initialize the G0 vector
            G0Vec = zeros(Natoms,1);

            % TO DO: implement without for loops
            for jj = 1:Natoms
                G0Vec(jj) = G0([x y], cps(jj,1:2), EC, ls, res, localGlobals);
            end
            
            % Eq 18
            GVec = A\G0Vec;
            % The relative error: abs(A*GVec - G0Vec)./G0Vec
            invError(x,y, Ee) = mean(abs((A*GVec - G0Vec)./G0Vec));
%             disp(['Mean relative error: ', ...
%                 num2str(mean(abs((A*GVec - G0Vec)./G0Vec)))])
            % Relative error in matrix inversion ~3e-15: good enough?
            % Add some error and see what happens? Jack squat, but still
            
            % LDOS through Eq 15
            LDOS(x,y, Ee) = -1/pi*imag(localGlobals(2)/(2*localGlobals(1)^2)*...
                (G0([x y], [x y], EC, ls, res, localGlobals) + ...
                sum(sVec'.*G0Vec.*GVec)));
            
%           Looking into the origins of negative LDOS
%             disp(['Signum sum for G0Vec*GVec: ', ...
%                 num2str(sum(sign(sum(real(G0Vec.*GVec)))))]);
% %             plot(x*y, sum(sign(real(G0Vec.*GVec))), 'x k');
% 
%             disp(['Signum sum for real part of G0Vec: ', ...
%                 num2str(sum(sign(real(G0Vec))))]);
%             
%             if LDOS(x,y,Ee) <0
%                 if imag(G0([x y], [x y], EC, ls, res, localGlobals))>0
%                     disp('Negative LDOS from G0!')
%                 end
%             
%                 if imag(sum(sVec'.*G0Vec.*GVec)) > 0
%                     disp('Negative LDOS from sVec.*G0Vec.*GVec!')
%                     
%                     disp(['Signum sum: ', ...
%                         num2str(sum(sign(imag(sVec'.*G0Vec.*GVec))))]);
%                     
% %                     disp(['Signum sum for sVec: ', ...
% %                         num2str(sum(sign(imag(sVec))))]);
% %                     disp(['Real components of sVec: ', num2str(sum(abs(real(sVec))))])
%                     disp(['Signum sum for real part of G0Vec: ', ...
%                         num2str(sum(sign(real(G0Vec))))]);
%                     
%                     disp(['Signum sum for GVec: ', ...
%                         num2str(sum(sign(imag(GVec))))]);
%                     
%                     disp(['Signum sum for G0Vec*GVec: ', ...
%                         num2str(sum(sign(real(G0Vec.*GVec))))]);                    
%                     
%                 end
%             end
            % remaining prefactor from G0: ms/hbar^2
            % Sum dimensions? A consistent implementation would use only
            % functions or only matrices, so .* should be best.
        end
        disp(strcat('LDOS computed at x=', num2str(x)));
        
    end
    toc;
    
    % LDOS normalisation to surface state DOS
    % Let's interpret the Green's function as a perturbation to the surface
    % state DOS background: should solve negative LDOS in most cases
    LDOS(:,:,Ee) = LDOS(:,:,Ee) + rho_surf;
    LDOS_norm(:,:,Ee) = LDOS(:,:,Ee)./rho_surf;
end

end


% UTILITY FUNCTIONS: G0 and k_disp cannot contain global variables

% Outgoing free Green's function
function [G0] = G0(rp, r, E, ls, res, localGlobals)
    % Stupid re-initialization to make parfor work
    hbar = localGlobals(1);
    ms = localGlobals(2);

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


% s and phase_Kondo are outside parfor loop, so global variables are fine

% s-wave scattering: Eq 17
function [s] = s(phase)
    global hbar    
    global ms
    
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


