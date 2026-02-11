function plmMRtrix3_v = standard2mrtrix(plmAll_v)
% Usage: plmMRtrix3_v = standard2mrtrix(plmAll_v)
%
% standard2mrtrix converts from coefficients of complex spherical
% harmonics (in standard format) to real coefficients of real spherical
% harmonics, in mrtrix3 format. See "Spherical Harmonics notes," section D
% for details.
%
% Input variables:
%
%   plmAll_v: 1 by (lMax+1)^2 row vector of (generally) complex
%               coefficients of a real-valued function with antipodal
%               symmetry. Coefficients are stored in 'standard' order:
%               (l,m) = (0,0), (1,-1), (1,0), (1,1), (2,-1), ...
%
% Output variables:
%
%   plmMRtrix3_v: 1 by (lMax+1)*(lMax+2)/2 row vector of real coefficients
%               of real spherical harmonic functions, as defined for
%               MRtrix3 (see the MRtrix3 docs).
%

% Get parameters:
nTermsStandard = length(plmAll_v);

% Derived parameters:
lMax = sqrt(nTermsStandard) - 1;
nTermsMRtrix3 = (lMax + 1)*(lMax + 2)/2;

% Loop over (l,m) combinations:
plmMRtrix3_v = zeros(1, nTermsMRtrix3);
for l = 0:2:lMax
    for m = -l:l
        % Get indices for current (l,m):
        indexStandard = (l+1)^2 - l + m;
        indexMRtrix3 = l*(l+1)/2 + m + 1;
        % Assign coefficient (see notes [D.14]):
        if (m <  0)
            plmMRtrix3_v(indexMRtrix3) = (-1)^m * sqrt(2) * ...
                imag(plmAll_v(indexStandard));
        elseif (m == 0)
            plmMRtrix3_v(indexMRtrix3) = plmAll_v(indexStandard);
        else
            plmMRtrix3_v(indexMRtrix3) = sqrt(2) * ...
                real(plmAll_v(indexStandard));
        end
    end
end

