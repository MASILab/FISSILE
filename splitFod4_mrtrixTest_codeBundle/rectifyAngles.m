function xOut_v = rectifyAngles(xIn_v)
% Usage: xOut_v = rectifyAngles(xIn_v)
%
% rectifyAngles flips (and unwinds, if necessary) phi in order to guarantee 
% that 0 <= theta and -2*pi <= phi <= 2*pi. This is useful when a nonlinear
% optimizer increments theta to become negative, but SH require theta >= 0.
%
% rectifyAngles is written for fixel orientation searches, as in
% simFixelDecomp5.m. It is assumed that theta <= pi.
%
% Input:
%
%   xIn_v: (2*nFixels x 1) array of FF polar and azimuthal angles. The
%       orientaion of the nth fixel axis is given by 
%                   theta = xIn_v(nFixels + n)
%       and 
%                   phi = xIn_v(n).
%
% Output:
%
%   xOut_v: (2*nFixels x 1) array of equivalent angles, where theta >= 0.
%
% rectifyAngles calls the function:
%   [theta, phi] = rectifyTheta(theta, phi)

% Get number of fixels:
nFixels = length(xIn_v)/2;

% Return immediately if theta >= 0 for all n:
if all(xIn_v((nFixels+1):(2*nFixels)) >= 0)
    % Nothing to fix, so return:
    xOut_v = xIn_v;
    return
end

% Otherwise, loop through fixels, fixing the theta values, if necessary:
xOut_v = zeros(size(xIn_v));
for fIndex = 1:nFixels
    % Get angles for current fixel:
    theta = xIn_v(nFixels + fIndex);
    phi = xIn_v(fIndex);
    % Trap for theta < 0:
    if (theta < 0)
        [theta, phi] = rectifyTheta(theta, phi);
    end
    % Unwind phi, if necessary:
    if (abs(phi) > 2*pi)
        phi = phi - sign(phi)*floor(abs(phi)/(2*pi)) * 2*pi;
    end
    % Set ouput values:
    xOut_v(fIndex) = phi;
    xOut_v(nFixels + fIndex) = theta;
end % End of fIndex FOR.
