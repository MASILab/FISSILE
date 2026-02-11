function P_v = evalP2(theta_v, phi_v, pCoeff_v)
% Usage: P_v = evalP2(theta_v, phi_v, pCoeff_v)
%
% evalPX evaluates P(theta_v, phi_v) reconstructed from the 
% SHD with coefficients pCoeff_v (coefficients in standard
% order, in a row vector). 
%
% Input variables are
%
%   theta_v:  nTheta x 1 column vector of theta values (in radians).
%   phi_v:    nTheta x 1 column vector of phi values (in radians).
%   pCoeff_v: 1 x nTerms row vector of FOD coefficients, in standard order.
%
% Output variables are:
%
%   P_v:    nTheta x 1 vector of FOD values. The FOD value for 
%                   (theta_v(n), phi_v(n))
%           is in P_v(n).
%
% evalPX calls the functions:
%   ynm_m = sphericalHarm(theta_v, phi_v, maxOrder)

% Preliminaries:
maxOrder = sqrt(length(pCoeff_v)) - 1;

% Calculate spherical harmonics at the desired orientation. ynm_m has one
% row for each (theta, phi) pair and nTerms columns (standard order):
ynm_m = sphericalHarm(theta_v, phi_v, maxOrder); 
P_v = real(ynm_m * pCoeff_v.');
