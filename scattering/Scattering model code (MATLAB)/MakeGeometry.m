function [cps] = MakeGeometry(res, Natoms, ecc, scale, scattererType, addFocus, focusType)
%GEOMETRY Creates a geometry in which the computations will be performed
%   Inputs:
%       - res: scan resolution, total number of pixels defined as (res*res)
%       - Natoms: how many atoms one wants the structure to consist of
%       - ecc: Eccentricity of the ellipse. 0 results in circle, max 1.
%       - scale: a scaling factor to specify the major axis wrt scan window
%       - scattererType: 0 for black disc, 1 for Kondo
%       - addFocus: if true, adds an atom to the focus point of the ellipse
%       - focusType: 0 for black disc, 1 for Kondo
%   Outputs:
%       - cps: scatterer point coordinates in pixel units, and scatterer
%           types
%   TO DO:
%       - Variable number of arguments
%       - Lattice discretization: shift adatoms to allowed sites
%       - Hardcode scale: useless parameter, replace with major axis length

cps = ones(Natoms, 3);

% How to make the interatom distance uniform? Some elliptical integral
cPoints = linspace(0,2*pi*(1-1/Natoms),Natoms);

% Ellipse parametrised through eccentricity
ellipse = [sqrt(1-ecc^2)*cos(cPoints); sin(cPoints)];

% Scale and move properly
corral = scale*res/2*ellipse + res/2;
% TO DO: check for bad indices

% Removing in-built rounding:
% cps(:,1:2) = floor(corral');
cps(:,1:2) = corral';

% Scatterer type allocation
cps(:,3) = scattererType.*ones(Natoms,1);

% Add an atom into a focus point (or middle if ecc == 0)
if addFocus
    focus = [res/2;scale*res/2*ecc + res/2;focusType];
    cps(length(cps) + 1,:) = focus;
end

end

