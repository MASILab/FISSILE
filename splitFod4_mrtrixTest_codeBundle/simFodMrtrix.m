function plmLabMrtrix_v = simFodMrtrix(theta_v, phi_v, f_v, vuFlag)
% Usage: plmLabMrtrix_v = simFodMrtrix(theta_v, phi_v, f_v, vuFlag)
%
% simFodMrtrix constructs the summation of <nFixels> ideal FODs (an "ideal" 
% FOD is a delta function, truncated to the maximum order, lMax).
%
% Input:
%
%   theta_v: nFixels x 1 column vector of polar angles (in radians). The
%           polar angle of the nth fixel is theta_v(n).
%   phi_v:  nFixels x 1 column vector of azimuthal angles (in radians). The
%           azimuthal angle of the nth fixel is phi_v(n).
%   f_v:    nFixels x 1 column vector of (relaxation-weighted) volume
%           fractions. The volume fraction of the nth fixel is f_v(n).
%   vuFlag: flag for displaying FOD. If vuFlag == 1, the net FOD is 
%           displayed.
%
% Output:
%
%   plmLabMrtrix_v: nTermsMRtrix x 1 column vector of the net FOD, 
%           written in the MRtrix3 format, where
%                   nTermsMRtrix = (lMax + 1)*(lMax + 2)/2.
%
% simFodMrtrix calls the functions:
%   pCoeff_v = pCoeffIdealFiber(theta, phi, maxOrder)
%   plmMRtrix3_v = standard2mrtrix(plmAll_v)

% Define parameters:
lMax = 6;   % Maximum spherical harmonic order. 

% Derived parameters:
nFixels = length(theta_v);

% Loop over fixels:
plmAll_v = zeros((lMax+1)^2, 1);    % Initialize net FOD (standard format).
for fIndex = 1:nFixels
    % Get coefficients of current fixel (in standard format):
    pCoeff_v = f_v(fIndex) * pCoeffIdealFiber(theta_v(fIndex), ...
        phi_v(fIndex), lMax);
    plmAll_v = plmAll_v + pCoeff_v.'; 
end

% Display net FOD, if desired:
if (vuFlag == 1)
    alpha = 1;  % Opacity.
    el = 30;    % Viewing angle elevation (degrees).
    az = 45;    % Viewing azimuthal angle (degrees).
    showFad3(plmAll_v, alpha, az, el)
    title('Net FOD')
end

% Convert to MRtrix3 format:
plmLabMrtrix_v = standard2mrtrix(plmAll_v);
