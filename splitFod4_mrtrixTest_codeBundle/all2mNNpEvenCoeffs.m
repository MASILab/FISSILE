function mNNpEvenCoeff_v = all2mNNpEvenCoeffs(allMCoeff_v)
% Usage: mNNpEvenCoeff_v = all2mNNpEvenCoeffs(allMCoeff_v)
%
% all2mNNpEvenCoeffs transforms a complex vector of SH coefficients [in 
% standard order: (p=0,m=0), (1,-1), (1,0), (1,1), etc.] describing
% a real function on a sphere to a vector of real SH coefficients for
% p=0:2:pMax, m>=0 (non-negative m) in the "non-negative m" order: (p=0,m=0), 
% (2,0), real part of (2,1), imaginary part of (2,1), real part of (2,2), 
% imaginary part of (2,2), etc. Note this gives 2p+1 (real-valued) 
% coefficients for each order, p (derived from the 2p+1
% complex coefficients in the input array). Note the total number of 
% coefficients is
%       nCoeff = Sum{2p+1, p=0:2:pMax} = (pMax/2 + 1)^2
% The simplification to non-negative m is due to the symmetry in m for 
% coefficients of a real-valued function:
%   c(p,-m) = (-1)^m * conj(c(p,m))
% which follows from the same symmetry for the spherical harmonics.
% The output is a column vector, regardless of whether the input is.

% Find pMax:
pMax = sqrt(length(allMCoeff_v)) - 1;
% Trap for bad input vectors:
if (pMax ~= floor(pMax))
    error('Length of input vector must be a perfect square...')
end

% Build output vector of coefficients:
mNNpEvenCoeff_v = [];
for p=0:2:pMax
    for m=0:p
        % Get complex coefficient for current (p,m):
        cIndex = p^2 + p + m + 1;   % Index in complex vector.
        if (m == 0)
            mNNpEvenCoeff_v = [mNNpEvenCoeff_v; real(allMCoeff_v(cIndex))];
        else
            mNNpEvenCoeff_v = [mNNpEvenCoeff_v; real(allMCoeff_v(cIndex)); ...
                 imag(allMCoeff_v(cIndex))];
        end
    end
end
