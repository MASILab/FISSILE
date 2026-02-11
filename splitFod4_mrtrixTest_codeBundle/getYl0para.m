function yl0_v = getYl0para(lMax)
% Usage: yl0_v = getYl0para(lMax)
%
% getYl0para calculates Yl0 values along the FF axis (theta = 0). 
% Values are calculated for even values of l and returned in 
% mNNpEvenReal1 format.
%
% Input variable:
%
%   lMax: maximum order of Yl0.
%
% Output variable:
%
%   yl0_v: nTerms x 1 array of Yl0(0,0) values, in mNNpEvenReal1 format.
%
% This is useful for constraining the FOD to have maximum amplitude
% along the FF axis. The value of the FOD is
%
%        P(0, 0) = Sum[pl0 * Yl0(0, 0), {l, 0, lMax}]
%
% assuming axial symmetry of the FOD.
%
% Use the expression in Arfken and Weber (4th ed.), p. 740 for m=0.
%
% getYl0perp calls the functions:
%   nTerms = getNumCoeffs(lMax, l0Flag, format_s)
%   [l_v, m_v, isReal_v] = index2lm_mNNpEvenReal(index_v)

% Preliminaries:
l0Flag = 1;     % Includes l=0.
nTerms = getNumCoeffs(lMax, l0Flag, 'mNNpEvenReal');
[l_v, m_v, isReal_v] = index2lm_mNNpEvenReal(1:nTerms);

% Get Yl0(0, 0) values. See Response Functions II (E.2b):
yl0_v = zeros(nTerms, 1);   % Initialize vector.
l0_v = l_v(m_v == 0);       % l values for (l,m=0) coefficients.
yl0_v(m_v == 0) = sqrt((2*l0_v + 1)/(4*pi));
