function logicValue = iseven(ipValue)
% Usage: logicValue = iseven(ipValue)
%
% iseven tests whether the input value is even. Non-integers
% and odd integers return 0.

% Test whether ipValue is a multiple of 2:
logicValue = (2*floor(ipValue/2) == ipValue);