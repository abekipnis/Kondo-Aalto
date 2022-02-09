function [cps_LD, trigrid] = LatticeDiscretize(cps, res, ls)
%LATTICEDISCRETIZE Snaps the cps points to triangular lattice points
%   Inputs:
%       - cps: (x[pix],y[pix],type)-matrix containing impurity sites and types
%       - res: pixel resolution of simulated scan
%       - ls: square scan window side length, [m]
%   Outputs:
%       - cps_LD: Lattice-discretized version of the cps input
%       - trigrid: the triangular lattice used to generate cps_LD

%% Initializations
global a
cps_LD = cps;

%% Make a triangular lattice
% This part mostly copy-pasted from https://se.mathworks.com/matlabcentral/answers/474193-how-to-generate-points-in-triangular-lattice-pattern
h_dist = a*res/ls; % Horizontal distance, pixel units
v_dist = sqrt(h_dist^2-(h_dist/2)^2); % Vertical distance
x_lim = res;
y_lim = res;

% Generate grid
trigrid = [];
y_current = 0;
xx = 0;
displacement = 0;
while y_current < y_lim
    if displacement == 0
        xx = [0:h_dist:x_lim]';
        yy = ones(length(xx), 1)*y_current;
        displacement = 1;
    else
        xx = [h_dist/2:h_dist:x_lim]';
        yy = ones(length(xx), 1)*y_current;
        displacement = 0;
    end
    trigrid = [trigrid; [xx,yy]];
    y_current = y_current + v_dist;
end

%% CPS nearest neighbors

% Find the nearest point in the triangular lattice
NN_IDs = knnsearch(trigrid, cps(:,1:2));

cps_LD(:,1:2) = trigrid(NN_IDs,:);
   

end

