function slmAll_m = uniTrans2standard(slm_m, lMax)
% Usage: slmAll_m = uniTrans2standard(slm_m, lMax)
%
% uniTrans2standard is a general purpose utility for converting signal 
% coefficients to the standard representation. 
%
% The signal coefficients of the nth fiber are held in the nth column of 
% slm_m. The input coefficients are represented in one of the following 
% formats:
%
%   m0pEven: m=0 (real-valued) coefficients for 
%               l=0:2:lMax (for l0Flag = 1) or 
%               l=2:2:lMax (for l0Flag = 0). 
%           sl0 coefficients are listed in order of increasing l. For 
%           example, 
%               [s(0,0), s(2,0), s(4,0), s(6,0), ...]
%   mNNpEvenReal: m>=0 (real-valued) coefficients for 
%               l=0:2:lMax (for l0Flag = 1) or 
%               l=2:2:lMax (for l0Flag = 0). 
%           slm coefficients are listed in order of increasing l and for
%           each l, increasing m. For example, 
%               [s(0,0), s(2,0), real(s(2,1)), imag(s(2,1)),  
%                real(s(2,2)), image(s(2,2)), ...]
%   mNNpEvenComp: m>=0 (complex-valued) coefficients for 
%               l=0:2:lMax (for l0Flag = 1) or 
%               l=2:2:lMax (for l0Flag = 0). 
%           slm coefficients are listed in order of increasing l and for
%           each l, increasing m. For example, 
%               [s(0,0), s(2,0), s(2,1), s(2,2), s(4,0), ...]
%           Also known as the 'minimal' format (the minimum number of
%           coefficients).
%   standard: complex-valued coefficients for l=0:lMax, m=-l:l, (l=0 is 
%           always included). slm coefficients are listed in order of 
%           increasing l and for each l, increasing m. For example,
%               [s(0,0), s(1,-1), s(1,0), s(1,1), s(2,-2), s(2,-1), s(2,0),
%                s(2,1), s(2,2), s(3,-3), ...]
%
% Note that if l0Flag = 0, the s00 coefficient is set to 0 in the output.
%
% The input parameters are 
%
%   slm_m: nTerms by nFibers matrix of signal coefficients. Coefficients 
%       for the nth fiber are held in the nth column.
%   lMax: maximum order represented in slm_m. Assumed to be even.
%
% The output argument is
%
% slmAll_m: (lMax+1)^2 by nFibers matrix of signal coefficients (in the 
%       same frame as slm_m). Coefficients for the nth fiber are in 
%       column n (in standard format). 
%
% History:
%   2021/06/09: created from showFiberSignals3 (AWA). 
%
% uniTrans2standard calls the functions:
%   [rep_s, l0Flag] = findSlmRep(nTerms, lMax)
%   [aOpt, AisoOpt] = sl02a_v3(ratio_s, s20, s40, s60)
%   slmAll_v  = mNNpEven2complexSHD(slm_v)
%   pCoeff_v = minimal2AllCoeff(pCoeffMinimal_v)
%   w_m = rotFodCoeffMatrix(maxOrder, alpha, beta, gamma)

% Derived parameters:
[nTerms, nFibers] = size(slm_m);

% Convert to standard format, if necessary. First find input representation:
[rep_s, l0Flag] = findSlmRep(nTerms, lMax);
repFlag_s = [rep_s, num2str(l0Flag)];   % Combined rep_s and l0Flag.
% Convert:
slmAll_m = zeros((lMax+1)^2, nFibers);   % Initialize standard format array.
switch repFlag_s
    case {'m0pEven0', 'm0pEven1'}
        l_v = 0:2:lMax;
        m0index_v = l_v.*(l_v-1)/2 + 1; % Indices of m=0 terms in mNNpEvenReal format.
        % Loop through fibers:
        for fIndex = 1:nFibers
            slm_v = slm_m(:, fIndex);
            slmAll_v = zeros(lMax*(lMax+3)/2 + 1, 1);
            % Set l=0 signal, if necessary:
            if (l0Flag == 0)
                % Set missing l=0 coefficient to 0:
                s00 = 0;
                slmAll_v(m0index_v) = [s00; slm_v];
            else
                % l=0 coefficient is included already:
                slmAll_v(m0index_v) = slm_v;
            end           
            % Change to standard format:
            slmAll_m(:, fIndex) = mNNpEven2complexSHD(slmAll_v);
        end
    case {'mNNpEvenReal0', 'mNNpEvenReal1'}
        % Loop through fibers:
        for fIndex = 1:nFibers
            slm_v = slm_m(:, fIndex);
            % Set l=0 signal, if necessary:
            if (l0Flag == 0)
                % Set missing l=0 coefficient to 0:
                s00 = 0;
                slmAll_v = [s00; slm_v];
            else
                % l=0 coefficient is included already:
                slmAll_v = slm_v;
            end           
            % Change to standard format:
            slmAll_m(:, fIndex) = mNNpEven2complexSHD(slmAll_v);
        end
    case {'mNNpEvenComp0', 'mNNpEvenComp1'}
        % Loop through fibers:
        for fIndex = 1:nFibers
            slm_v = slm_m(:, fIndex);
            % Set l=0 signal, if necessary:
            if (l0Flag == 0)
                % Set missing l=0 coefficient to 0:
                s00 = 0;
                slm_v = [s00; slm_v];
            end           
            % Change to standard format. Note that minimal2AllCoeff takes
            % and returns row vectors:
            slmAll_m(:, fIndex) = minimal2AllCoeff(slm_v.').';
        end
    case {'standard0', 'standard1'}
        % No need to transform frames:
        if (l0Flag == 0)
            % Loop through fibers:
            for fIndex = 1:nFibers
                slm_v = slm_m(:, fIndex);
                % Set missing l=0 coefficient to 0:
                s00 = 0;
                slm_v = [s00; slm_v];   % Pre-pend s00.
                % Save it:
                slmAll_m(:, fIndex) = slm_v;
            end
        else
            slmAll_m = slm_m;
        end
end
