function [l_v, m_v, isReal_v] = index2lm_mNNpEvenReal(index_v)
% Usage: [l_v, m_v, isReal_v] = index2lm_mNNpEvenReal(index_v)
%
% index2lm_mNNpEvenReal takes a row vector of indices and returns the
% corresponding order (l), degree (m), and isReal flag (1 for the
% real part of a coefficient, 0 for the imaginary part) in mNNpEvenReal1 
% format.
%
% The input arguments are
%
% index_v: row vector of indices into a coefficient vector (in mNNpEvenReal
%       format).
%
% Output arguments are
%
% l_v: row vector of l values for each element in index_v (mNNpEvenReal format).
% m_v: row vector of m values for each element in index_v (mNNpEvenReal format). 
% isReal_v: row vector of flags for the real part of a coefficient (1 for
%   real, 0 for imaginary).
%
% index2lm_mNNpEvenReal calls the functions
%   nTerms = getNumCoeffs(lMax, l0Flag, rep_s)

% Find lMax. Increment lMaxTrial until the number of terms (for that lMax 
% in mNNpEvenReal format) exceeds the maximum index in index_v:
lMaxTrial = 0;
while (getNumCoeffs(lMaxTrial, 1, 'mNNpEvenReal') < max(index_v))
    lMaxTrial = lMaxTrial + 2;  % Assume lMax is even.
end
lMax = lMaxTrial;

% Create arrays:
lAll_v = [];
mAll_v = [];
isRealAll_v = [];
% Loop over orders:
for l = 0:2:lMax
    for m = 0:l
        lAll_v = [lAll_v, l];
        mAll_v = [mAll_v, m];
        isRealAll_v = [isRealAll_v, 1];
        % Add imaginary part:
        if (m > 0)
            lAll_v = [lAll_v, l];
            mAll_v = [mAll_v, m];
            isRealAll_v = [isRealAll_v, 0];
        end
    end
end

% Now select the elements of each array using the indices in index_v:
l_v = lAll_v(index_v);
m_v = mAll_v(index_v);
isReal_v = isRealAll_v(index_v);
