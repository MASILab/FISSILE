function logicValue = isodd(ipValue)
% Usage: logicValue = isodd(ipValue)
%
% isodd tests whether the input value is odd. Note that
% non-integers and even integers return 0.

% Test whether ipValue is 1 greater than a multiple of 2:
logicValue = (2*floor((ipValue-1)/2) + 1 == ipValue);