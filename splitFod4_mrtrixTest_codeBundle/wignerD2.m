function D_v = wignerD2(l, mPrime, m, alpha_v, beta_v, gamma_v)
% Usage: D_v = wignerD2(l, mPrime, m, alpha_v, beta_v, gamma_v)
%
%   wignerDX returns the rotation operator for spherical harmonics:
%       Rotated[Ylm] = Sum{D(l, mPrime, m, alpha, beta, gamma) YlmPrime}
% where the sum is over the index mPrime. [See Rose, Elementary theory 
% of angular momentum, pp. 52-3, for details]. The angles should 
% be passed in column vectors, where the nth rotation is specified
% by the Euler angles alpha_v(n), beta_v(n), and gamma_v(n). The output
% vector, D_v, has the same size as alpha_v, etc.
%   Both algorithms for calculating d(beta) from Rose are given. Note 
% that eqn (4.13) avoids the divergence in (4.14) when beta = pi.
%
%   wignerDX calls the functions:
%       F_v = hypergeometric(a, b, c, z_v)

% Preliminaries:
alg_s = 'Rose_4.13';    % 'Rose_4.13' or 'Rose_4.14'.

% Make sure all angle vectors are column vectors:
if (size(alpha_v, 1) < size(alpha_v, 2))
    alpha_v = alpha_v.';
end
if (size(beta_v, 1) < size(beta_v, 2))
    beta_v = beta_v.';
end
if (size(gamma_v, 1) < size(gamma_v, 2))
    gamma_v = gamma_v.';
end

% Find d:
switch alg_s
    case 'Rose_4.14'
        swapFlag = 0;
        if (m > mPrime)
            % Swap m and mPrime, and change sign of beta:
            temp = m;
            m = mPrime;
            mPrime = temp;
            beta_v = -beta_v;
            swapFlag = 1;
        end
        % Find d assuming mPrime >= m:
        F_v = hypergeometric(mPrime-l, -m-l, mPrime-m+1, -tan(beta_v/2).^2);
        d_v = sqrt(factorial(l-m) * factorial(l+mPrime) / (factorial(l+m) * factorial(l-mPrime))) * ...
            ((cos(beta_v/2).^(2*l+m-mPrime) .* (-sin(beta_v/2)).^(mPrime-m)) / factorial(mPrime-m)) ...
            .* F_v;
        if (swapFlag == 1)
            % Swap m and mPrime back to their original values:
            mPrime = m;
            m = temp;
        end
    case 'Rose_4.13'
        % The summation index kappa must satisfy:
        %   kappa >= 0
        %   kappa >= m - mPrime
        %   kappa <= l + m
        %   kappa <= l - mPrime
        % so 
        %   max(0, m - mPrime) <= kappa <= min(l+m, l-mPrime)
        kappa_v = max(0, m-mPrime):min(l+m, l-mPrime);  % Must have at least one element (0).
        dSum_v = zeros(size(beta_v));
        for kappaIndex = 1:length(kappa_v)
            kappa = kappa_v(kappaIndex);
            term_v = (-1)^kappa * cos(beta_v/2).^(2*l+m-mPrime-2*kappa) .* ...
                (-sin(beta_v/2)).^(mPrime-m+2*kappa) / ...
                (factorial(l-mPrime-kappa)*factorial(l+m-kappa)*...
                factorial(kappa+mPrime-m)*factorial(kappa));
            dSum_v = dSum_v + term_v;
        end
        d_v = sqrt(factorial(l+m)*factorial(l-m)*factorial(l+mPrime) * ...
            factorial(l-mPrime)) * dSum_v;
end

% Now calculate D [Rose, eqn. 4.12]:
D_v = exp(-1i*mPrime*alpha_v) .* d_v .* exp(-1i*m*gamma_v);


    
    
    