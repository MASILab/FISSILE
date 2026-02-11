function nTerms = getNumCoeffs(lMax, l0Flag, rep_s)
% Usage: nTerms = getNumCoeffs(lMax, l0Flag, rep_s)
%
% getNumCoeffs returns the number of coefficients required to specify
% the spherical harmonic coefficients of a function of max order lMax
% in the rep_s representation, either including the l=0 term or not.
% lMax is assumed to be even. 
%
% Note that for any lMax>0, nTerms uniquely determines rep_s and l0Flag.
%
% Input arguments are
%
%   lMax: maximum order of slm coefficients (assumed to be even). lMax
%       can be a scalar or array.
%   l0Flag: scalar = 1 if l=0 coefficient is included, 0 otherwise. For 
%       the standard representation, l0Flag=1 always.
%   rep_s: string indicating the representation. Options are
%       'm0pEven', 'mNNpEvenReal', 'mNNpEvenComp', 'standard'.
%
% The output variable is
%
%   nTerms: number of terms in the spherical harmonic representation.
%       nTerms has the same dimensions as lMax (scalar or array). 
%
% The representations are defined using coefficients slm as follows:
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

% Find the number of terms (in each representation) for given lMax, l0Flag,
% and representation:
switch rep_s
    case 'm0pEven'
        nTerms = lMax/2 + l0Flag;  
    case 'mNNpEvenReal'
        nTerms = lMax.*(lMax+3)/2 + l0Flag;
    case 'mNNpEvenComp'
        nTerms = (lMax/2 + 1).^2 + (l0Flag-1); 
    case 'standard'
        nTerms = (lMax+1).^2 + (l0Flag-1);
end
