function [Topo] = LDOS2Topo(LDOS, Energies)
%LDOS2Topo Integrate LDOS over energy to get a constant current map 
%   Inputs:
%       - LDOS: the LDOS stack to be integrated (xres,yres,Energies)
%       - Energies: the energies corresponding to the LDOS maps.
%   Outputs:
%       - Topo: A (resx, resy, Energies-1) topography map stack
%       - Topo_norm: Normalized topographies
%   To do:
%       - Negative LDOSes will be a problem

global E_F

% First determine if Energies spans over EF: no sense integrating if not
if sign(min(Energies-E_F)) == sign(max(Energies-E_F))
    error('The energy range must include E_F: go back and compute!')
end

% Find the energy closest to Fermi
indF = find(Energies - E_F == min(abs(Energies - E_F)),1);

% Compute the integrals
Topo = zeros(size(LDOS));

% Negative energies present?
if indF > 1
    % Direction: maybe add a minus sign afterwards
    range = linspace(indF,1,indF);
    Topo(:,:, range) = cumtrapz(Energies(range), LDOS(:,:,range), 3);
    Topo(:,:,indF:end) = cumtrapz(Energies(indF:end), LDOS(:,:,indF:end), 3);
else
    Topo = cumtrapz(Energies, LDOS, 3);
end


end

