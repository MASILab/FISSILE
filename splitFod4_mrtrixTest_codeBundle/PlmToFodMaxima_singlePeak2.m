function fodMax_v = PlmToFodMaxima_singlePeak2(ProbCoeff_v, realFlag, nSamples, nTheta, nPhi, flipFlag)
% usage: fodMax_v = PlmToFodMaxima_singlePeak2(ProbCoeff_v, realFlag, nSamples, nTheta, nPhi, flipFlag)
% find the magnitude and orientation of the single-peak FOD function.
% 
% input: ProbCoeff_v -- SH coefficients of the FOD function.
%        realFlag    -- 1 if ProbCoeff_v is real, 0 otherwise
%        nSamples -- # of evenly distributed sampling points
%        nTheta -- # of theta values in range of [0, pi]
%        nPhi -- # of phi values in range of [0, 2*pi]
%        flipFlag -- if 1, the FOD peak will be flipped to the upper half
%                    spherical (fodMax_v(2)<=pi/2)
%
% output: fodMax_v -- 1 by 3 vector: [magnitue, theta, phi]
%
% by XH
%
% History:
% last modified 12/09/08 (XH)
% 12/29/2020: replaced readHardDirs_xh2 with getDirectionsWrapper(n) (AWA).

order = sqrt(length(ProbCoeff_v)) - 1;

if ~exist('realFlag', 'var') || isempty(realFlag)
    realFlag = 0;
end

if ~exist('nSamples', 'var') || (isempty(nSamples))
    if ~exist('nTheta', 'var') || (isempty(nTheta)) || ~exist('nPhi', 'var') || (isempty(nPhi))
        nSamples = 1002;
    else      
        [gridTheta_v, gridPhi_v, gridPhi_m, gridTheta_m] = gridThetaPhi(nTheta, nPhi);
        sh_3d = sphericalHarm_xh4(gridTheta_v, gridPhi_v, order, 1, 0, realFlag);
        nNegFlag = 1;
        minFlag = 0;
        gridP_m = shCoeffToGridFun(ProbCoeff_v, realFlag, nNegFlag, minFlag, sh_3d);
        [maxMag, maxId] = max(gridP_m(:));
        fodMax_v = [maxMag, gridTheta_m(maxId), gridPhi_m(maxId)];
        return
    end
end

if ~exist('flipFlag', 'var') || isempty(flipFlag)
    flipFlag = 0;
end

%uniDirs_m = readHardDirs_xh2(nSamples);
uniDirs_m = getDirectionsWrapper(nSamples);
[uniTheta_v, uniPhi_v] = xyz2sph(uniDirs_m(:, 1), uniDirs_m(:, 2), uniDirs_m(:, 3));
% sample the spherical function
gridFlag = 0;
evenOrderFitFlag = 1;
% evenOrderFitFlag = 0;
uniDesign_m = sphericalHarm_xh5(uniTheta_v, uniPhi_v, order, gridFlag, evenOrderFitFlag, realFlag);
fod_v = real(uniDesign_m * ProbCoeff_v);
[maxMag, maxId] = max(fod_v);
fodMax_v = [maxMag, uniTheta_v(maxId), uniPhi_v(maxId)];

% flip to the upper hemisphere if necessary
if flipFlag == 1
   if fodMax_v(2) > pi/2
       fodMax_v(2) = pi-fodMax_v(2);
       if fodMax_v(3) < pi
           fodMax_v(3) = fodMax_v(3) + pi;
       else
           fodMax_v(3) = fodMax_v(3) - pi;
       end
   end
end