function gridP_m = shCoeffToGridFun(coeff_v, realFlag, nNegFlag, minFlag, sh_3d, nTheta, nPhi)
% usage: gridP_m = shCoeffToGridFun(coeff_v, realFlag, nNegFlag, minFlag, sh_3d, nTheta, nPhi)
%        (nNegFlag...optional)
% reconstruct the spherical function on a grid of points based on its SH
% coefficents.
% nNegFlag: =1 to set all negative values of gridP_m to zeros.
% minusminFlag: =1 to minus min(gridP), use for plot QBALL ODF
% 
% call:
%   gridThetaPhi
%   sphericalHarm_xh4
%
%
% by XH
% last update 08/06/08

nTerms = length(coeff_v);
order = sqrt(nTerms)-1;

if ~exist('realFlag', 'var') || isempty(realFlag)
    realFlag = 0;
end
if ~exist('nNegFlag', 'var') || isempty(nNegFlag)
    nNegFlag = 1;
end
if ~exist('minFlag', 'var') || isempty(minFlag)
    minFlag = 0;
end

if ~exist('sh_3d', 'var') || (isempty(sh_3d))
    if ~exist('nTheta', 'var') || (isempty(nTheta)) || ~exist('nPhi', 'var') || (isempty(nPhi))
        nTheta = 80;
        nPhi = 160;
    end   
    [gridTheta_v, gridPhi_v] = gridThetaPhi(nTheta, nPhi);
    sh_3d = sphericalHarm_xh4(gridTheta_v, gridPhi_v, order, 1, 0, realFlag);
else
    nTheta = size(sh_3d,1);
    nPhi = size(sh_3d,2);
end

gridP_m = zeros(nTheta, nPhi);
for term = 1:nTerms
    gridP_m = gridP_m + real(coeff_v(term) * squeeze(sh_3d(:,:,term)));
end

if nNegFlag == 1
    gridP_m(gridP_m<0) = 0;
end

if minFlag == 1
    gridP_m = gridP_m - min(gridP_m(:));
end

