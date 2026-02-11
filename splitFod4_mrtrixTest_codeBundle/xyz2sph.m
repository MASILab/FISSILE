function [theta, phi] = xyz2sph(x, y, z)
% Transform Cartesian to spherical coordinates.
% theta: angle from +z axis
% phi: counterclockwise angle in xy-plane from +x axis
% angles are in radians
% YS

[phi, el, ~] = cart2sph(x, y, z);
theta = pi/2 - el;
end
