function xStart_m = getMyCustomStartPoints(x0_v, alphaMax, nPoints)
% Usage: xStart_m = getMyCustomStartPoints(x0_v, alphaMax, nPoints)
%
% getMyCustomStartPoints creates a random cloud of start points, uniformly
% distributed in the cone of (half-) angle <alphaMax>, centered on the
% directions specified by <x0_v> on each sphere implied by <x0_v>. This is
% useful for multistart searches over fixel orientations on nFixels unit 
% spheres, for example.
% 
% Input arguments:
%
%   x0_v: (2*nFixels) x 1 array of estimated fixel orientations. The nth 
%           fixel has estimated orientation given by theta = x0_v(nFixels+n)
%           and phi = x0_v(n), in radians.
%   alphaMax: (half-) cone angle in each unit sphere, in radians. All 
%           random start points are within alphaMax of x0_v in each sphere.
%   nPoints: number of start points returned. 
%
% Output variables:
%
%   xStart_m: nPoints x (2*nFixels) array of start points for multistart
%           search. The nth startpoint is on the nth row.
%
% History:
%   20240722: created (AWA).
%
% getMyCustomStartPoints calls the functions:
%   randDir_m = getRandomSpherePoints(theta0, phi0, alphaMax, nPoints)
%   vuPointsOnSphere3(theta_v, phi_v, color_s, title_s)

% Preliminaries:
nFixels = length(x0_v)/2;

% Generate start points for each fixel:
xStart_m = zeros(nPoints, 2*nFixels);
% figure
for fIndex = 1:nFixels
    % Add random angle (in a cone) to get start point for search:
    theta0 = x0_v(nFixels + fIndex);
    phi0 = x0_v(fIndex);
    randDir_m = getRandomSpherePoints(theta0, phi0, alphaMax, nPoints);
    thetaStart_v = acos(randDir_m(:, 3));
    phiStart_v = atan2(randDir_m(:, 2), randDir_m(:, 1));
    % Add randomly rotated axis:
    xStart_m(:, nFixels + fIndex) = thetaStart_v;   % theta values.
    xStart_m(:, fIndex) = phiStart_v;               % phi values.
    % Show start points for each subFixel:
    title_s = ['Fixel ', num2str(fIndex), ' start points'];
    % subplot(1, nFixels, fIndex)
    % vuPointsOnSphere3(thetaStart_v, phiStart_v, 'r', title_s);
end
