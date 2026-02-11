function ynm_m = sphericalHarm_xh4(theta_v, phi_v, maxOrder, gridFlag, evenFlag, realFlag)
% Usage: ynm_m = sphericalHarm_xh4(theta_v, phi_v, maxOrder, gridFlag, evenFlag, realFlag)
%
% sphericalHarm_xh1 calculates the values of the spherical
% harmonics Ynm(theta, phi) at the points given by the
% elements of theta_v, phi_v.
%
% based on sphericalHarm
% gridFlag == 1:
% the output is a 3d matrix. element(p,q) of each page of
% ynm_m gives the values of Ynm at theta_v(p), phi_v(q). each page of ynm_m
% gives the values  for one (n,m) combination, where n and m increase in
% the 'standard' order: m increases from -n to +n, then n increments)
%
% gridFlag == 0 (default):
% sphericalHarm calculates the values of the spherical
% harmonics Ynm(theta, phi) at the points given by the
% elements of theta_v, phi_v. Row p of ynm_m gives the
% values of Ynm at theta_v(p), phi_v(p). Each column of
% ylm_m gives the values for one (n,m) combination,
% where n and m increase from left to right (the
% 'standard' order: m increases from -n to +n, then
% n increments).
%
% ver3: based on ver1, complex fitting
%       add evenFlag, when evenFlag==1, gives only even order items
%
% ver4: based on ver3 and HKJ's sphericalHarmonicRealY.m
%       add input realFlag
%
% by XH
% last update 06/19/08



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

% by default, gridFlag = 0
if ~exist('evenFlag', 'var') || isempty(evenFlag)
    evenFlag = 0;
end

if ~exist('gridFlag', 'var') || isempty(gridFlag)
    gridFlag = 0;
end

if ~exist('realFlag', 'var')
    realFlag = 0;
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
    if (evenFlag == 1) && (mod(n,2) == 1)
        continue
    else
        if realFlag == 0    % complex SH
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
                exp_v = exp(i*m*phi_v.');
                if gridFlag == 0
                    sh_v = k * legendre_v .* exp_v;	% Ynm in a row vector.
                    % Add this term to ynm matrix:
                    ynm_m = [ynm_m, sh_v.'];
                else
                    sh_m = k * legendre_v.' * exp_v;	% Ynm in a matrix.
                    % Add this term to ynm matrix:
                    ynm_m = cat(3, ynm_m, sh_m);
                end
            end

        else    % real SH
            if n == 0   % Set to row vector.
                Plm_m = legendre(n,cos(theta_v))';
            else
                Plm_m = legendre(n,cos(theta_v));
            end

            for m = -n:n
                The_v = sqrt( (2*n+1)/2*prod(1:(n-abs(m)))/prod(1:(n+abs(m))) )*Plm_m(abs(m)+1,:);

                if gridFlag == 0
                    if m > 0
                        ynm_v = The_v'.*cos(m*phi_v)/sqrt(pi);
                    elseif m == 0
                        ynm_v = The_v'/sqrt(2*pi);
                    else
                        ynm_v = The_v'.*sin(-m*phi_v)/sqrt(pi);
                    end
                    ynm_m = [ynm_m ynm_v];

                else
                    if m > 0
                        ynm_v = The_v'*cos(m*phi_v)'/sqrt(pi);
                    elseif m == 0
                        ynm_v = The_v'*ones(1,length(phi_v))/sqrt(2*pi);
                    else
                        ynm_v = The_v'*sin(-m*phi_v)'/sqrt(pi);
                    end
                    ynm_m = cat(3, ynm_m, ynm_v);
                end

            end % m

        end % end of if realFlag

    end % end of (evenFlag == 1) & (mod(n,2) == 1)

end % end of n
