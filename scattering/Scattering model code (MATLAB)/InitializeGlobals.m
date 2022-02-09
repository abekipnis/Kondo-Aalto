function InitializeGlobals(surface)
%INITIALIZEGLOBALS Initializes the necessary constants as global variables
%   The following functions require these initializations:
%       - MakeGeometry.m
%       - ComputeLDOS.m
%       - PlotAndSaveLDOS.m
%   Input:
%       - surface: string with either 'Ag' or 'Cu'

% Natural constants
global E0
global E_F
global hbar
global ms
global rho_surf
global delta_bg
global delta_primes
global Gamma
global eps0
global phase_n
global phase_bd
global a

% possible additions:
% global scale % as in 2a/ls
% global atomDiam

% Fermi energy: let's see if this matters at all
E_F = 0;

hbar = 6.5821e-16;  % eVs

% Electron mass
m_e = 5.6856e-12;   % eV/((m/s)^2)

% Surface specifics
switch surface
    case 'Cu'
        E0 = -450e-3;       % eV
        ms = 0.38*m_e;      % eV/((m/s)^2)
        a = 2.546e-10;      % m
    case 'Ag'
        E0 = -67e-3;        % eV
        ms = 0.4*m_e;       % eV/((m/s)^2)
        a = 2.887e-10;      % m
    otherwise
        disp([surface ' surface state parameters unknown: using Cu(111) parameters instead'])
        E0 = -450e-3;        % eV
        ms = 0.38*m_e;      % eV/((m/s)^2)
end

% Density of surface states above E0
rho_surf = ms/(pi*hbar^2);

% Kondo phase shift for Cu(111) - who knows how these work for Ag(111)
delta_bg = pi/4;        % rad
delta_primes = 1.5;     % no dimension
Gamma = 9e-3;           % eV
eps0 = E_F - 1e-3;      % eV

% "Normal" phase shifts
phase_n = pi/2;

% For Fe on Cu(111), for instance
phase_bd = 1i*inf;

end

