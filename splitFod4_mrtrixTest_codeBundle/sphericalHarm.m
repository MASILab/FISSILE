function ynm_m = sphericalHarm(theta_v, phi_v, maxOrder)
% Usage: ynm_m = sphericalHarm(theta_v, phi_v, maxOrder)
%
% sphericalHarm calculates the values of the spherical
% harmonics Ynm(theta, phi) at the points given by the
% elements of theta_v, phi_v. Row p of ynm_m gives the
% values of Ynm at theta_v(p), phi_v(p). Each column of
% ylm_m gives the values for one (n,m) combination, 
% where n and m increase from left to right (the 
% 'standard' order: m increases from -n to +n, then
% n increments).

% Make sure angle vectors are column vectors. Assume they
% have the same size:
[nRows, nCols] = size(theta_v);
if (nCols > 1)
	% Transpose from row to column vectors:
	theta_v = theta_v.';
	phi_v = phi_v.';
	% Make sure input arrays are vectors:
	if (nRows > 1)
		error('sphericalHarm accepts vector arguments only')
	end	
end

% Create ynm_m matrix. Note that matlab's 'legendre' function 
% returns the Associated Legendre polynomials, which are related 
% to the spherical harmonics of the same (n,m) by
% 		Ynm(theta, phi) = k * Pnm(cos(theta)) * exp(i*m*phi)
% where
%	 k = sqrt( (2*n+1) * prod(1:(n-m)) / (4*pi * prod(1:(n+m))))
% The output of legendre(n, x_v) is a matrix where rows correspond
% to m = 0, 1, 2, ..., n and columns correspond to the elements of
% x_v:
ynm_m = [];
for n = 0:maxOrder
    % Trap for illegal value of x=cos(theta):
    if any(imag(cos(theta_v)) ~= 0) || (any(abs(cos(theta_v))) > 1)
        dummy = 1;
    end
    % Proceed with legal values:
	legendre_m = legendre(n, cos(theta_v.'));
	% Calculate -ve m polynomials from +ve index values:
	for m = -n:n
		if (m < 0)
	        negFactor = (-1)^m * prod(1:(n+m)) / prod(1:(n-m));
		    legendre_v = negFactor * legendre_m(-m+1, :);
		else
			legendre_v = legendre_m(m+1, :);
		end
		k = sqrt((2*n+1) * prod(1:(n-m)) / (4*pi * prod(1:(n+m))));
		exp_v = exp(1i*m*phi_v.');
	    sh_v = k * legendre_v .* exp_v;	% Ynm in a row vector.
		% Add this term to ynm matrix:
		ynm_m = [ynm_m, sh_v.'];
	end
end
