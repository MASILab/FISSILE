function Q = axialAsymIndex(plm_v, lMax)
% Usage: Q = axialAsymIndex(plm_v, lMax)
%
% axialAsymIndex calculates the axial asymmetry index, Q, for the FOD
% described by the SH coefficients plm_v. See Response Functions II (H.57)
% for details.
%
% Input variables:
%
%   plm_v: nTerms x 1 array of SH coefficients describing the FOD, in any
%       format recognized by uniTrans2standard.m.
%   lMax: scalar maximum order in SH representation.
%
% Output variables:
%
%   Q: scalar axial anisotropy index (Q>=0).
%
% History:
% 2024/08/28: renamed from axialAnisoIndex (AWA).
%
% axialAsymIndex calls the functions:
%   slmAll_m = uniTrans2standard(slm_m, lMax)
%   [l_v, m_v] = index2orderAndDegree(lMax)

% Convert to standard format:
plmAll_v = uniTrans2standard(plm_v, lMax);

% Get standard order indices:
[l_v, m_v] = index2orderAndDegree(lMax);

% Calculate Q:
P2Int = sum(plmAll_v.*conj(plmAll_v));
Pmean2Int = sum(plmAll_v(m_v == 0).^2);
Q = (P2Int - Pmean2Int)/Pmean2Int;
