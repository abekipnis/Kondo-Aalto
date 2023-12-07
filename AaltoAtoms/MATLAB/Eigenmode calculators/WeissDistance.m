function [W] = WeissDistance(r, N, r0)
%WEISSDISTANCE Computes the "Weiss distance", or average nearest neighbor
% distance for circular corral wall atoms/molecules
%   For details see https://arxiv.org/pdf/2304.06571.pdf
% Inputs:
%   - r: radius of the corral
%   - N: number of atoms in the corral
%   - r0: radius of a central atom (set to 0 if considering an empty corral)
% 
% Output:
%   - W: average distance between corral atoms

% Initialize a circle with N points
thetas = linspace(0, 2*pi*(1-1/N), N);
circle = [r*cos(thetas); r*sin(thetas)]';

% get the distances with pdist
dists = pdist(circle);
dists = dists(1:N-1);

% Elaborate handling of direct scattering in case of central
% atoms/molecules
if r0 > 0
    % Number of atoms "covered" by a central atom of radius r0
    N_out = ceil((2*pi - 4*acos(r0/r))/(2*pi)*N);

    [~,indRem] = maxk(dists, N_out);
    
    dists(indRem) = r/2;
end

% Average
W = mean(dists);

end