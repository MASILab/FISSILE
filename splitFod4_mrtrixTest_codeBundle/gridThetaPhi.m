function [gridTheta_v, gridPhi_v, gridPhi_m, gridTheta_m] = gridThetaPhi(nTheta, nPhi)
% by XH 041408

gridTheta_v = pi * (0:(nTheta-1)).' / (nTheta-1);
gridPhi_v = 2*pi * (0:(nPhi-1)).' / (nPhi-1);

if nargout > 2
    [gridPhi_m, gridTheta_m] = meshgrid(gridPhi_v.', gridTheta_v);
end