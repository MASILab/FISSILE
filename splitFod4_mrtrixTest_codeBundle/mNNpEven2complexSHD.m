function plm_v = mNNpEven2complexSHD(slm_v)
% Usage: plm_v = mNNpEven2complexSHD(slm_v)
%
% mNNpEven2complexSHD takes a vector of SH coefficients, slm_v, in a 
% compact representation (describing real-valued functions that are
% symmetric around the origin) and returns the corresponding full SHD
% representation, plm_v, using complex coefficients (all orders 
% and m values). 
%   
% Note that this function does not involve the Descoteaux convention 
% for real SHDs (see all2complexSHD.m for that). 
%
% The input argument is:
%   slm_v: compact representation of a real, symmetric function in terms 
%       of real-valued coefficients for non-negative m (column vector). 
%       If the complex coefficients 
%           p(l,m) = u(l,m)+i*v(l,m), [u and v real]
%       for all (l,m), then for each order, the function is defined by 
%       the the 2l+1 coefficients 
%           u(l,0),u(l,1),v(l,1),u(l,2),v(l,2),...,u(l,l),v(l,l) 
%       This defines the "non-negative m" order. Since the function is
%       symmetric around the origin, only orders l=0,2,4,...,lMax
%       are included. 
%
% The output argument is:
%   plm_v: column vector of complex coefficients in "standard order":
%           (l=0,m=0),(l=1,m=-1),(l=1,m=0),(l=1,m=1), ..., (l=lMax,m=lMax)
%       for all orders l=0,1,...,lMax.

% Get pMax, the solution to p^2 + 3p - 2(Ne-1) = 0 [see (B.41), amended to
% include p=0 term]:
Ne = length(slm_v);
pMax = (-3 + sqrt(9 + 4*2*(Ne-1)))/2;   % pMax assumed to be even.

% Trap for input vector of wrong length:
if (pMax ~= round(pMax))
    error(['Input vector has the wrong length (pMax = ', num2str(pMax), ...
        ')'])
end

% Loop over orders:
for p=0:pMax
    if (p==0)
        plm_v = slm_v(1);
    elseif isodd(p)
        plm_v = [plm_v; zeros(2*p+1, 1)];
    else
        % p is even and >0. Get segment of slm_v for current p:
        nePrev = (p-2)*(p+1)/2 + 1; % Number of elements for lower p.
        index_v = (nePrev+1):(nePrev+2*p+1);
        cReal_v = slm_v(index_v);   % slm_v elements for current p.
        % Loop over m values. For m>0,
        %   p(l,-m) = (-1)^m * [u(l,m)-i*v(l,m)]:
        % while
        %   p(l,0) = u(l,0)
        % and
        %   p(l,m) = u(l,m)+i*v(l,m)
        % Calculate these in order:
        for m=-p:p 
            % Get u and v for current (l,|m|):
            if (m==0)
                u = cReal_v(1);
                v = 0;
            else
                u = cReal_v(2*abs(m));
                v = cReal_v(2*abs(m)+1);
            end
            % Now get complex coefficient:
            if (m<0)
                plm = (-1)^m * (u-1i*v);
            elseif (m==0)
                plm = u;
            else
                % m>0:
                plm = u + 1i*v;
            end
            % Append element for current m:
            plm_v = [plm_v; plm];
        end % End of m FOR.
    end     % End of p IF.
end         % End of p FOR.
        