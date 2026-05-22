function rotated_fod = rotate_nnMpEvenRealRotation(euler_m, fod)
%rotate_nnMpEvenRealRotation Rotate an FOD in nnMpEvenRotation format
%
% The input parameters are:
%
%   euler_m: 1 by 3 matrix of Euler angles ([alpha, beta, gamma]) 
%       giving the orientation of the lab frame,
%       relative to the fiber frame, or more generally from the initial
%       frame to the frame described by the Euler angles.
%   fod: the Nterms ((l//2+1)*(l+1)) x N matrix of FODs to rotate
%
% Output variables are:
%
%   rotated_fod : Nterms x N rotated FOD
%
% nnMpEvenRotationMatrix calls the functions:
%   nnMpEvenRotationMatrix
%   index2lm_mNNpEvenReal

% Find lMax
Nterms = size(fod, 1);
lMaxTrial = 0;
while (getNumCoeffs(lMaxTrial, 1, 'mNNpEvenReal') < Nterms)
    lMaxTrial = lMaxTrial + 2;  % Assume lMax is even.
end
lMax = lMaxTrial;

[l_values, m_values, ~] = index2lm_mNNpEvenReal(1:Nterms);
R = zeros(Nterms, Nterms);
for l = 0:2:lMax
    R_l = nnMpEvenRotationMatrix_v2(euler_m, l);
    R(l_values == l, l_values == l) = R_l;
end
rotated_fod = fod' * R;
rotated_fod = rotated_fod';

end