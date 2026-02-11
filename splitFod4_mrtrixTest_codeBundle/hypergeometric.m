function F_v = hypergeometric(a, b, c, z_v)
% Usage: F_v = hypergeometric(a, b, c, z_v)
%
%   hypergeometric calculates the hypergeometric function, 2F1(a, b; c; z)
% for
%   1). a or b equal to an integer <= 0 and
%   2). c not zero or a negative integer.
% where z_v is a column vector. The output, F_v has the same size as z_v.
% [See Rose, Elementary theory of angular momentum, p. 53].

% Make sure z_v is a column vector:
if (size(z_v, 2) > 1)
    z_v = z_v.';
end
nz = length(z_v);

% Check conditions:
cLogic = (c ~= 0) & ~((floor(c) == c) & (c < 0));
aLogic = (floor(a) == a) & (a <= 0);
bLogic = (floor(b) == b) & (b <= 0);
% Test c:
if (cLogic == 0)
    error(['c = ', num2str(c), ' is zero or a negative integer'])
end
% Test z and a or b:
if ((aLogic == 0) && (bLogic == 0))
    error(['Neither a = ', num2str(a), ' nor b = ', ...
        num2str(b), ' is an integer <= 0'])
end

% Calculate terms in power series expansion of F:
minInt = min(aLogic * a, bLogic * b);   % Find lowest negative integer.
maxOrder = -minInt;
F_m = zeros(nz, maxOrder);
for order = 1:maxOrder
    F_m(:, order) = ( prod(a:(a+order-1)) * prod(b:(b+order-1)) / prod(c:(c+order-1)) ) * ...
        (z_v.^order) / factorial(order);
end
% F_m   % Show terms, for debugging.
F_v = ones(nz, 1) + sum(F_m, 2);   % Zeroth order plus higher orders.
