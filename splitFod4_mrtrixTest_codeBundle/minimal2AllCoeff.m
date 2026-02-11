function pCoeff_v = minimal2AllCoeff(pCoeffMinimal_v)
% Usage: pCoeff_v = minimal2AllCoeff(pCoeffMinimal_v)
%
% minimal2AllCoeff takes a (row) vector of the minimal equivalent 
% representation of an FOD (even order, non-negative m spherical harmonic 
% coefficients) and expands it to include all orders, assuming zero 
% odd-order contribution and all m, assuming the FOD is real, i.e., 
%       p(l,m) = (-1)^m * conj(p(l,-m))
% The output (row) vector, pCoeff_v, includes all orders through maxOrder, 
% which is assumed to be even. 

% Find maxOrder:
maxOrder = 2*(sqrt(length(pCoeffMinimal_v)) - 1);

% Loop through orders, dropping the odd order coefficients:
pCoeff_v = zeros(1, (maxOrder+1)^2);
remaining_v = pCoeffMinimal_v;
for order = 0:2:maxOrder
    % Get coefficients for this order, non-negative m:
    nonNegmCoeff_v = remaining_v(1:(order+1));
    if (order>0)
        % Flip positive m to negative m coefficients:
        negm_v = (-order:-1);
        negmCoeff_v = (-1).^negm_v .* conj(fliplr(nonNegmCoeff_v(2:end)));
        orderCoeff_v = [negmCoeff_v, nonNegmCoeff_v];
    else
        orderCoeff_v = nonNegmCoeff_v;
    end
    % Drop these from remaining coefficient vector:
    remaining_v = remaining_v((order+2):end);
    % Copy coefficients to pCoeff_v:
    startIndex = order^2 + 1;
    stopIndex = (order + 1)^2;
    pCoeff_v(startIndex:stopIndex) = orderCoeff_v;
end
    