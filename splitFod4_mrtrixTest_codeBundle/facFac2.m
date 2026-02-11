function val_v = facFac2(n_v)
% Usage: val_v = facFac2(n_v)
%
%   facFac2 returns the value n!! for each integer, n, in n_v. This is 
% defined as 1 for n = 0 or 1, and otherwise is
%       n!! = n(n-2)(n-4)...2   for even n
%       n!! = n(n-2)(n-4)...1   for odd n
% 

% Trap for any non-integer values:
if any(round(n_v) ~= n_v)
    % n is non-integer:
    error(['Not all input values, ', num2str(n_v), ', are integers'])
end

% Loop through input values:
val_v = ones(size(n_v));    % Initialize output.
for j = 1:length(n_v)
    n = n_v(j);
    % Calculate only for integer n>1:
    if (n > 1)
        % Set up product:
        p_v = n:-2:1;
        val_v(j) = prod(p_v);
    end
end
