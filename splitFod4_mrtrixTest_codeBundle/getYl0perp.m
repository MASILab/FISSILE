function yl0_v = getYl0perp(lMax)
% Usage: yl0_v = getYl0perp(lMax)
%
% getYl0perp calculates Yl0 values in the plane perpendicular to the FF
% axis (theta = pi/2). Values are calculated for even values of l and
% returned in mNNpEvenReal1 format.
%
% Input variable:
%
%   lMax: maximum order of Yl0.
%
% Output variable:
%
%   yl0_v: nTerms x 1 array of Yl0(pi/2,0) values, in mNNpEvenReal1 format.
%
% This is useful for constraining the FOD to have zero amplitude
% perpendicular to the FF axis. The value of the FOD is
%
%        P(pi/2, phi) = Sum[pl0 * Yl0(pi/2, 0), {l, 0, lMax}]
%
% assuming axial symmetry of the FOD. Note that this is independent of phi
% because only m=0 terms are included in the sum.
%
% Use the expression in Arfken and Weber (4th ed.), p. 741 for even (l + m).
%
% getYl0perp calls the functions:
%   nTerms = getNumCoeffs(lMax, l0Flag, format_s)
%   [l_v, m_v, isReal_v] = index2lm_mNNpEvenReal(index_v)
%   function val_v = facFac2(n_v)

% Preliminaries:
l0Flag = 1;     % Includes l=0.
nTerms = getNumCoeffs(lMax, l0Flag, 'mNNpEvenReal');
[l_v, m_v, isReal_v] = index2lm_mNNpEvenReal(1:nTerms);

% Get Yl0(pi/2, 0) values. See Response Functions II (E.2b):
yl0_v = zeros(nTerms, 1);   % Initialize vector.
l0_v = l_v(m_v == 0);       % l values for (l,m=0) coefficients.
yl0_v(m_v == 0) = (-1).^(l0_v/2) .* sqrt((2*l0_v + 1)/(4*pi)) .* ...
    factorial(l0_v) ./ facFac2(l0_v).^2;
