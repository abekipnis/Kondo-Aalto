function [skipGrid] = Cps2Grid(cps, res, ls)
%CPS2GRID Converts a list of impurity points to a res*res pixel grid
%   Inputs:
%       - cps: (x[pix],y[pix],type)-matrix containing impurity sites and types
%       - res: pixel resolution of simulated scan
%       - ls: square scan window side length, [m]
%   Output:
%       - skipGrid: a (res*res) binary map, with circular areas around the 
%               impurity points set to 1
%
% Note that the impurity type is not used in this function.

% Removing in-built rounding from cps and adding it here:
cps(:,1:2) = round(cps(:,1:2));

Natoms = size(cps,1);
skipGrid = zeros(res, res);

% Atom diameter in meters and pixels
atomDiam_m = 14e-10;
atomDiam = res/ls*atomDiam_m;

% In pixels the diameter must be odd
if mod(floor(atomDiam),2) == 0
    atomDiam_pix = ceil(atomDiam);
else
    atomDiam_pix = floor(atomDiam);
end

% Produce the circle mask to be added on impurity points
[columnsInImage, rowsInImage] = meshgrid(1:atomDiam_pix, 1:atomDiam_pix);
circleMask = (rowsInImage - ceil(atomDiam_pix/2)).^2 ...
    + (columnsInImage - ceil(atomDiam_pix/2)).^2 <= (atomDiam_pix/2).^2;

% Mask every atom
for jj = 1:Natoms
    side = floor(atomDiam_pix/2);
    skipGrid((cps(jj,1)-side):(cps(jj,1)+side),...
        (cps(jj,2)-side):(cps(jj,2)+side))...
        = skipGrid((cps(jj,1)-side):(cps(jj,1)+side),...
        (cps(jj,2)-side):(cps(jj,2)+side))+ circleMask;
end

end

