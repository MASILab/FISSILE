function randDir_m = getRandomSpherePoints(theta0, phi0, alphaMax, nPoints)
% Usage: randDir_m = getRandomSpherePoints(theta0, phi0, alphaMax, nPoints)
%
% getRandomSpherePoints finds the location of nPoints randomly selected 
% on the unit sphere within the (half) cone angle alphaMax from the center
% point at (theta0, phi0). The random points are uniformly sampled on the
% spherical surface.
%
% Input parameters:
%
%   theta0: polar angle of the cone axis on the unit sphere (in radians).
%   phi0: azimuthal angle of the cone axis on the unit sphere (in radians).
%   alphaMax: half cone angle, i.e., points at the intersection of the cone
%       and the unit sphere are separated from the cone axis by angle
%       alphaMax at the origin. For now, it's assumed that alphaMax <= pi/2.
%       Units are radians.
%   nPoints: number of randomly (and uniformly) sampled points on the
%       spherical surface patch inside the cone.
%
% Output variable:
%
%   randDir_m: nPoints x 3 array of coordinates of the random points. Each
%       row gives the [x, y, z] coordinates of one randomly sampled point.
%
% See Response Functions II, section H for details.
% 
% History:
%   20240611: created (AWA).
%
% getRandomSpherePoints calls the functions:
%   rot_m = yRotMatrix(angle)
%   rot_m = zRotMatrix(angle)

% Trap for bad input values:
if (alphaMax > pi/2 || alphaMax < 0)
    error(['alphaMax = ', num2str(alphaMax), ' is not a legal value.'])
end

% Get nPoints random values of sin(alpha), cos(alpha), and beta (for polar
% angle alpha and azimuthal angle beta on the unit sphere, prior to rotation):
%rng(1);     % Reset random number generator.
sinAlpha_v = sin(alphaMax) * rand(1, nPoints);  % Uniform on [0, sin(alphaMax)].
cosAlpha_v = sqrt(1 - sinAlpha_v.^2);
beta_v = 2*pi * rand(1, nPoints);   % Uniform on [0, 2*pi].

% Make initial coordinate array--one point in each column:
coord_m = [sinAlpha_v.*cos(beta_v); sinAlpha_v.*sin(beta_v); cosAlpha_v];

% Rotate by theta0 around y axis, then by phi0 around z axis:
coordPrime_m = zRotMatrix(phi0) * yRotMatrix(theta0) * coord_m;
randDir_m = coordPrime_m.';
