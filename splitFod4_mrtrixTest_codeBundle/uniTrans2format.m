function slm_m = uniTrans2format(slmAll_m, opFormat_s)
% Usage: slm_m = uniTrans2format(slmAll_m, opFormat_s)
%
% uniTrans2format is a general purpose utility for converting signal 
% coefficients from standard order to the representation specified by 
% opFormat_s. 
%
% The input parameters are 
%
%   slmAll_m: (lMax+1)^2 by nFibers matrix of signal coefficients.  
%       Coefficients for the nth fiber are held in the nth column.
%   opFormat_s: [pmFormat_s, l0Flag_s], where pmFormat_s specifies the 
%       order and assumptions placed on the slm coefficients (see below) 
%       and l0Flag_s = '1' if the l=p=0 term is included and is '0' 
%       otherwise.
%
% The output argument is
%
% slm_m: nTerms by nFibers matrix of signal coefficients (in the 
%       same frame as slmAll_m). Coefficients for the nth fiber are in 
%       column n in the format specified by opFormat_s.   
%                
% Options for pmFormat_s are the following:
%
%   'm0pEven': m=0 (real-valued) coefficients for 
%               l=0:2:lMax (for l0Flag = 1) or 
%               l=2:2:lMax (for l0Flag = 0). 
%           sl0 coefficients are listed in order of increasing l. For 
%           example, 
%               [s(0,0), s(2,0), s(4,0), s(6,0), ...]
%   'mNNpEvenReal': m>=0 (real-valued) coefficients for 
%               l=0:2:lMax (for l0Flag = 1) or 
%               l=2:2:lMax (for l0Flag = 0). 
%           slm coefficients are listed in order of increasing l and for
%           each l, increasing m. For example, 
%               [s(0,0), s(2,0), real(s(2,1)), imag(s(2,1)),  
%                real(s(2,2)), image(s(2,2)), ...]
%   'mNNpEvenComp': m>=0 (complex-valued) coefficients for 
%               l=0:2:lMax (for l0Flag = 1) or 
%               l=2:2:lMax (for l0Flag = 0). 
%           slm coefficients are listed in order of increasing l and for
%           each l, increasing m. For example, 
%               [s(0,0), s(2,0), s(2,1), s(2,2), s(4,0), ...]
%           Also known informally as the 'minimal' format (the minimum
%            number of coefficients).
%   'standard': complex-valued coefficients for l=0:lMax, m=-l:l, (l=0 is 
%           always included). slm coefficients are listed in order of 
%           increasing l and for each l, increasing m. For example,
%               [s(0,0), s(1,-1), s(1,0), s(1,1), s(2,-2), s(2,-1), s(2,0),
%                s(2,1), s(2,2), s(3,-3), ...]
%           This is the input format and is provided just for completeness.
%
% History:
%   2021/06/09: created from showFiberSignals3 (AWA).
%   2022/06/22: created from uniTrans2format (AWA). 
%
% uniTrans2format calls the functions:
% remove   [rep_s, l0Flag] = findSlmRep(nTerms, lMax)
%   nTerms = getNumCoeffs(lMax, l0Flag, pmFormat_s)
%   mNNpEvenCoeff_v = all2mNNpEvenCoeffs(allMCoeff_v)
% remove  [aOpt, AisoOpt] = sl02a_v3(ratio_s, s20, s40, s60)
%   slmAll_v  = mNNpEven2complexSHD(slm_v)
%   pCoeff_v = minimal2AllCoeff(pCoeffMinimal_v)
%   pCoeffMinimal_v = all2MinimalCoeff(pCoeff_v)
% remove  w_m = rotFodCoeffMatrix(maxOrder, alpha, beta, gamma)

% Derived parameters:
[nStandardTerms, nFibers] = size(slmAll_m);
lMax = sqrt(nStandardTerms) - 1;

% Split information in opFormat_s:
pmFormat_s = opFormat_s(1:(end-1));
l0Flag = str2double(opFormat_s(end));

% Convert to output format, if necessary. First initialize output array:
nTerms = getNumCoeffs(lMax, l0Flag, pmFormat_s);
slm_m = zeros(nTerms, nFibers);   % Initialize standard format array.
switch opFormat_s
    case {'m0pEven0', 'm0pEven1'}
        l_v = (2*(1-l0Flag)):2:lMax;
        m0index_v = l_v.*(l_v+1) + 1; % Indices of m=0 terms in standard format.
        % Loop through fibers:
        for fIndex = 1:nFibers
            slm_m(:, fIndex) = slmAll_m(m0index_v, fIndex);
        end
    case {'mNNpEvenReal0', 'mNNpEvenReal1'}
        % Loop through fibers:
        for fIndex = 1:nFibers
            slmAll_v = slmAll_m(:, fIndex);
            mNNpEvenCoeff_v = all2mNNpEvenCoeffs(slmAll_v);
            % Account for l=0 coefficient:
            if (l0Flag == 0)
                % Drop l=0:
                slm_v = mNNpEvenCoeff_v(2:end);
            else
                % Inlcude l=0 coefficient:
                slm_v = mNNpEvenCoeff_v;
            end           
            % Add to op array:
            slm_m(:, fIndex) = slm_v;
        end
    case {'mNNpEvenComp0', 'mNNpEvenComp1'}
        % Loop through fibers:
        for fIndex = 1:nFibers
            slmAll_v = slmAll_m(:, fIndex);
            pCoeffMinimal_v = all2MinimalCoeff(slmAll_v.').';   % Col. vector.
            % Account for l=0 coefficient:
            if (l0Flag == 0)
                % Drop l=0 coefficient:
                slm_v = pCoeffMinimal_v(2:end);
            else
                slm_v = pCoeffMinimal_v;
            end           
            % Add to op array:
            slm_m(:, fIndex) = slm_v;
        end
    case {'standard0', 'standard1'}
        % No need to transform frames, just account for l=0 coefficient:
        if (l0Flag == 0)
            slm_m = slmAll_m(2:end, :);
        else
            slm_m = slmAll_m;
        end
end
