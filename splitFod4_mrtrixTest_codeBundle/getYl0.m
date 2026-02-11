function Yl0_m = getYl0(lMax, theta_v)
% Usage: Yl0_m = getYl0(lMax, theta_v)
%
% getYl0 takes a (column) vector of theta (polar angle) values and returns
% the corresponding values of Yl0 for l = 0:2:lMax. The product of Yl0_m
% with (the column vector) sl0Prime_v gives S'(theta), the signal at
% polar angles theta_v in the FF.
%
% Input arguments;
%   lMax: maximum order (assumed to be even).
%   theta_v: nTheta x 1 column vector of polar angles (in radians).
%
% Output arguments:
%   Yl0_m: nTheta x (lMax/2 + 1) array of Yl0 values. The nth row
%       corresponds to the polar angle theta_v(n) and the pth column
%       corresponds to l = 2*(p-1).
%
% getYl0 is useful for calculating the response function S'(theta) in the
% FF.
%
% History:
%   2024/05/08: Created (AWA).
%
% getYl0 calls the functions:

% Setup variables:
l_v = 0:2:lMax;
nTheta = length(theta_v);

% Loop over even values of l:
Yl0_m = zeros(nTheta, length(l_v));
for lIndex = 1:length(l_v)
    % Get current l value:
    l = l_v(lIndex);
    % Get the values of Pl0, the Legendre polynomial for order l, degree m:
    Plm_m = legendre(l, cos(theta_v));   % Each row corresponds to one value of m.
    Yl0_v = sqrt((2*l+1)/(4*pi)) * Plm_m(1,:);     % Save m=0 values.
    Yl0_m(:, lIndex) = Yl0_v.';
end
