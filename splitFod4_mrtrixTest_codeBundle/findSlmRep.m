function [rep_s, l0Flag] = findSlmRep(nTerms, lMax)
% Usage: [rep_s, l0Flag] = findSlmRep(nTerms, lMax)
%
% findSlmRep returns the format (representation) of the signal, that is
% which coefficients are explicitly listed and in what order they appear.
%
% Input arguments are:
%   
%   nTerms: number of coefficients.
%   lMax: maximum order of slm coefficients (assumed to be even).
%
% Output variables are:
%   rep_s: string indicating the signal representation. Options are
%       'm0pEven', 'mNNpEvenReal', 'mNNpEvenComp', 'standard'.
%   l0Flag: 1 if l=0 coefficient is included, 0 otherwise. For the
%       standard representation, l0Flag=1 always.
%
% The coefficients are represented in one of the following formats:
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
% findSlmRep calls the functions:
%   nTerms = getNumCoeffs(lMax, l0Flag, rep_s)

% Define parameters:

% Find the number of terms (in each representation) for given lMax and
% l0Flag (added to the variable name as a suffix):
nTerms_m0pEven0 = getNumCoeffs(lMax, 0, 'm0pEven');   
nTerms_m0pEven1 = getNumCoeffs(lMax, 1, 'm0pEven');
nTerms_mNNpEvenReal0 = getNumCoeffs(lMax, 0, 'mNNpEvenReal');
nTerms_mNNpEvenReal1 = getNumCoeffs(lMax, 1, 'mNNpEvenReal');
nTerms_mNNpEvenComp0 = getNumCoeffs(lMax, 0, 'mNNpEvenComp');      
nTerms_mNNpEvenComp1 = getNumCoeffs(lMax, 1, 'mNNpEvenComp');      
nTerms_standard = getNumCoeffs(lMax, 1, 'standard');

% Find representation with matching nTerms.  Note that for any lMax>0, 
% nTerms uniquely determines rep_s and l0Flag:
switch nTerms
    case nTerms_m0pEven0
        rep_s = 'm0pEven';
        l0Flag = 0;
    case nTerms_m0pEven1
        rep_s = 'm0pEven';
        l0Flag = 1;
    case nTerms_mNNpEvenReal0
        rep_s = 'mNNpEvenReal';
        l0Flag = 0;
    case nTerms_mNNpEvenReal1
        rep_s = 'mNNpEvenReal';
        l0Flag = 1;  
    case nTerms_mNNpEvenComp0
        rep_s = 'mNNpEvenComp';
        l0Flag = 0;
    case nTerms_mNNpEvenComp1
        rep_s = 'mNNpEvenComp';
        l0Flag = 1; 
    case nTerms_standard
        rep_s = 'standard';
        l0Flag = 1;
    otherwise
        error('Unknown representation...')
end
