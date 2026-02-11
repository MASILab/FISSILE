function plmAll_v = mrtrix2standard(plmMRtrix3_v)
% Usage: plmAll_v = mrtrix2standard(plmMRtrix3_v)
%
% mrtrix2standard converts from coefficients of the real spherical
% harmonics used in MRtrix3 (real coefficients of real spherical harmonics, 
% in MRtrix3 format) to standard complex format. See "Spherical Harmonics 
% notes," section D for details.
%
% Input variables:
%
%   plmMRtrix3_v: 1 by (lMax+1)*(lMax+2)/2 row vector of real coefficients
%               of real spherical harmonic functions, as defined for
%               MRtrix3 (see the MRtrix3 docs).
%
% Output variables:
%
%   plmAll_v: 1 by (lMax+1)^2 row vector of (generally) complex
%               coefficients of a real-valued function with antipodal
%               symmetry. Coefficients are stored in 'standard' order:
%               (l,m) = (0,0), (1,-1), (1,0), (1,1), (2,-1), ...
%
% The inverse function is standard2mrtrix.m.

% Get lMax. Note that nTermsMRtrix3 = (lMax + 1)*(lMax + 2)/2:
nTermsMRtrix3 = length(plmMRtrix3_v);
% Solve quadratic eqn: 
%       lMax^2 + 3*lMax - 2*(nTermsMRtrix3-1) = 0
% or 
%       lMax = (-3 + sqrt(9 + 4*2*(nTermsMRtrix3-1)))/2
%       lMax = (sqrt(1 + 8*nTermsMRtrix3) - 3)/2
lMax = (sqrt(1 + 8*nTermsMRtrix3) - 3)/2;

% Derived parameters:
nTermsStandard = (lMax+1)^2;

% Loop over (l,m) combinations:
plmAll_v = zeros(1, nTermsStandard);
for l = 0:2:lMax
    for m = -l:l
        % Get indices for current (l,m):
        indexStandard = (l+1)^2 - l + m;
        indexStandardOpp = (l+1)^2 - l - m; % Index for m with opposite sign.
        % Note that (l+1)*(l+2)/2 - l + m = (l^2 + l)/2 + 1 + m:
        indexMRtrix3 = l*(l+1)/2 + m + 1; 
        % Assign coefficient (use notes Eq. [D.14] and symmetry [D.11]):
        if (m <  0)
            % MRtrix coefficients for m<0 give imaginary parts of clm.
            % Add imaginary part of clm to plmAll_v:
            plmAll_v(indexStandard) = plmAll_v(indexStandard) + ...
                1i * (-1)^m * plmMRtrix3_v(indexMRtrix3)/sqrt(2);
            plmAll_v(indexStandardOpp) = plmAll_v(indexStandardOpp) + ...
                -1i * plmMRtrix3_v(indexMRtrix3)/sqrt(2);            
        elseif (m == 0)
            % MRtrix coefficient for m=0 is cl0:
             plmAll_v(indexStandard) = plmMRtrix3_v(indexMRtrix3);
        else
            % MRtrix coefficients for m>0 give real parts of clm. Add the
            % real part of clm to plmAll_v:
            plmAll_v(indexStandard) = plmAll_v(indexStandard) + ...
                plmMRtrix3_v(indexMRtrix3)/sqrt(2);
            plmAll_v(indexStandardOpp) = plmAll_v(indexStandardOpp) + ...
                (-1)^m * plmMRtrix3_v(indexMRtrix3)/sqrt(2);            
        end
    end
end

