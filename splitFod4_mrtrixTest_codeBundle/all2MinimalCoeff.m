function pCoeffMinimal_v = all2MinimalCoeff(pCoeff_v)
% Usage: pCoeffMinimal_v = all2MinimalCoeff(pCoeff_v)
%
% all2MinimalCoeff converts a full (row) vector of spherical harmonic coefficients,
% pCoeff_v, in standard order, to its minimal equivalent representation. The odd 
% orders are dropped as are the negative m coefficients for even orders. 
% The (complex) coefficients of the minimal representation are returned in 
% the shorter vector pCoeffMinimal_v.

% Find maxOrder:
maxOrder = sqrt(length(pCoeff_v)) - 1;

% Loop through orders, dropping the odd order coefficients:
pCoeffMinimal_v = [];
for order = 0:2:maxOrder
    % Keep coefficients for non-negative m only. 'order^2' coefficients
    % for lower orders, 'order' coefficients for negative m:
    startIndex = order^2 + order + 1;   
    stopIndex = (order + 1)^2;
    pCoeffMinimal_v = [pCoeffMinimal_v, ...
        pCoeff_v(startIndex:stopIndex)];
end
    