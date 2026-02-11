function P = evalP2(ori_v, pCoeff_v)
% Usage: P = evalP2(ori_v, pCoeff_v)
%
% evalP evaluates P(theta, phi) reconstructed from the 
% SHD with coefficients pCoeff_v (coefficients in standard
% order, in a row vector). evalP is written to be called by an optimization
% function (e.g., fminsearch) to find the maximum of P(theta, phi), 
% for example. The orientation vector, ori_v, has elements 
% [theta, phi].
%
% evalP calls the functions:
%   ynm_m = sphericalHarm(theta_v, phi_v, maxOrder)

% Preliminaries:
maxOrder = sqrt(length(pCoeff_v)) - 1;
theta = ori_v(1);
phi = ori_v(2);

% Calculate spherical harmonics at the desired orientation:
ynm_v = sphericalHarm(theta, phi, maxOrder); % Row vector.
P = real(sum(ynm_v .* pCoeff_v));
