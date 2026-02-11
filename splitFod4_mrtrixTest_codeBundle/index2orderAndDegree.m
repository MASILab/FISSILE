function [l_v, m_v] = index2orderAndDegree(lMax)
% Usage: [l_v, m_v] = index2orderAndDegree(lMax)
%
% index2lm returns the order and degree of each element of a coefficient
% vector in standard format through the order lMax.
%
% Input argument:
%
%   lMax: maximum order in the standard format vector. Assumed to be even.
%
% Output variables:
%
%   l_v: 1 by (lMax+1)^2 row vector of l values for each element of a 
%       standard format coefficient vector through order lMax.
%   m_v: 1 by (lMax+1)^2 row vector of m values for each element of a 
%       standard format coefficient vector through order lMax.
%
% The inverse function is orderAndDegree2index.m.

% Loop through elements in standard order:
l_v = [];
m_v = [];
for l = 0:lMax
    for m = -l:l
        l_v = [l_v, l];
        m_v = [m_v, m];
    end
end
