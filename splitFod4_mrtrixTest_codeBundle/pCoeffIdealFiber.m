function pCoeff_v = pCoeffIdealFiber(theta, phi, maxOrder)
% Usage: pCoeff_v = pCoeffIdealFiber(theta, phi, maxOrder)
%
% pCoeffIdealFiber calculates the spherical harmonic expansion
% coefficients for an ideal fiber (two delta functions, 
% symmetric around the origin) at orientation (theta, phi). 
% The coefficients are in standard order (in a row vector).
%
% Input arguments:
%
%   theta: polar angle of fixel axis (in radians).
%   phi: azimuthal angle of fixel axis (in radians).
%   maxOrder: maximum order of spherical harmonics representing the FOD
%       (assumed to be an even integer). 
%
% Output variables:
%
% pCoeff_v: 1 by (maxOrder+1)^2 row vector of FOD coefficients, in
%   'standard' format--the order (l) increases from left to right and for 
%   each order, the degree m increases from -l to l. Note that odd-order 
%   coefficients are included (for compatibility with other routines), but 
%   are identically zero. Hence the order of the coefficients in the vector
%   is
%   (l,m) = (0,0), (1,-1), (1,0), (1,1), (2,-2), (2,-1), ..., (maxOrder, maxOrder)
%   Note also that the coefficients are generally complex. 
% 
% pCoeffIdealFiber calls the functions:
%   ynm_m = sphericalHarm(theta_v, phi_v, maxOrder)

% Find indices of even orders:
evenIndex_v = [];
for order = 0:2:maxOrder
    evenIndex_v = [evenIndex_v, (order^2 + 1):(order+1)^2];
end

% Find complex conjugate of spherical harmonics at (theta, phi):
coeff_v = conj(sphericalHarm(theta, phi, maxOrder));

% Get even coefficients for the plm:
pCoeff_v = zeros(1, (maxOrder+1)^2);
pCoeff_v(evenIndex_v) = coeff_v(evenIndex_v);
