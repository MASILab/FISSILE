function acc = angularCorrCoeff(p_v, q_v)
% Usage: acc = angularCorrCoeff(p_v, q_v)
%
% angularCorrCoeff calculates the correlation coefficient
% of two functions, P(theta, phi) and Q(theta, phi), ovar
% all orientations, based on their spherical harmonic 
% decompositions. The coefficients of the spherical harmonic 
% expansions of P and Q are p_v and q_v, respectively 
% (coefficients are given in 'standard order', i.e., 
% p_v(1) is the (0,0) coefficient, p_v(2) is the (1,-1)
% coefficient, and so on). 
%
% ver _xh: 030507 by xh
%          make sure p_v and q_v have same size


% Preliminaries:

[dim1, dim2] = size(p_v);
if dim1 ~= size(q_v,1)
    q_v = q_v.';
end

% Make sure vectors have the same length:
np = length(p_v);
nq = length(q_v);
nCoeffs = min([np, nq]);
p_v = p_v(1:nCoeffs);
q_v = q_v(1:nCoeffs);

% Subtract the mean value from the expansions of P and Q:
zmp_v = p_v(2:nCoeffs);
zmq_v = q_v(2:nCoeffs);

% Calculate norm of each vector:
pNorm = sqrt(sum(zmp_v .* conj(zmp_v)));
qNorm = sqrt(sum(zmq_v .* conj(zmq_v)));

% Calculate normalized vectors and correlation coefficient:
if (pNorm > 0 & qNorm > 0)
    np_v = zmp_v / pNorm;
    nq_v = zmq_v / qNorm;

    % Calculate correlation coefficient:
    acc = real(sum(np_v .* conj(nq_v)));
else
    acc = 0;
end