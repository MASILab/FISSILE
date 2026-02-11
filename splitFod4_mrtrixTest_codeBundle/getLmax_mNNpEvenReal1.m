function L = getLmax_mNNpEvenReal1(nTerms)
% Usage: L = getLmax_mNNpEvenReal1(nTerms)
%
% getLmax_mNNpEvenReal1 takes the number of coefficients, nTerms, of the SH
% representation of an FOD or DW signal, in mNNpEvenReal1 format, and
% returns the maximum order, L. See Response Functions II (H.9), for
% example, for details. Use the relation
%
%   nTerms = L^2/2 + 3*L/2 + 1
%
% or 
% 
%   L^2 + 3*L + 2*(1-nTerms) = 0.

% Solve the quadratic:
L = (-3 + sqrt(9 + 4*2*(nTerms-1)))/2;
