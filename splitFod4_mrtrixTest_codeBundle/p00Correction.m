function [plmPrimeOut_m, p00Correction_v, p00Correction2_v] = p00Correction(plmPrime_m)
% Usage: [plmPrimeOut_m, p00Correction_v, p00Correction2_v] = p00Correction(plmPrime_m)
%
% p00Correction evaluates the integral of each FOD on the unit sphere, 
% from theta = pi/2 - beta to pi/2 + beta (the 'belt' around the middle
% of the unit sphere). Beta is the elevation angle of the first exremum
% from the midpoint, theta = pi/2 (beta = 0). See Response Functions II, 
% (H.25) for details.
%
% Input variables:
%
%   plmPrime_m: nTerms x nFibers array of FOD coefficients, in mNNpEvenReal1
%               format.
%
% Output variables are 
%
%   plmPrimeOut_m: nTerms x nFibers with corrected p00n values (mean value
%               on the belt around the equator nearly equal to zero).
%   p00Correction_v: value to add to the p00 coefficient to make the belt
%               integral equal zero.
%   p00Correction2_v: correction values (p00Correction_v) after subtracting
%               their mean (to guarantee they sum to zero).
%
% History:
%   2024/06/30: created from evalBeltIntegral3.m (AWA). 
%
% p00Correction calls the functions:
%   L = getLmax_mNNpEvenReal1(nTerms)
%   [l_v, m_v, isReal_v] = index2lm_mNNpEvenReal(index_v)
%   pl0All_v = uniTrans2standard(pl0Prime_v, lMax)
%   P = evalP2(ori_v, pCoeff_v)

% Preliminaries:
[nTerms, nFibers] = size(plmPrime_m);        % In mNNpEvenReal1 format.
lMax = getLmax_mNNpEvenReal1(nTerms);
[l_v, m_v, isReal_v] = index2lm_mNNpEvenReal(1:nTerms);
pl0Prime_m = plmPrime_m(m_v == 0, :); % Only m=0 values contribute to the integral.

% Derived parameters:
pl0All_m = uniTrans2standard(pl0Prime_m, lMax); 

% Find beta (elevation angle of first extremum from theta = pi/2) for each
% fiber:
nTheta = 181;
index90 = 91;
theta_v = (0:(nTheta-1))*pi/180;    % In radians.

% Evaluate each FOD at each theta:
P_m = zeros(nTheta, nFibers);
intVal_v = zeros(1, nFibers);
meanVal_v = zeros(1, nFibers);
p00Correction_v = zeros(1, nFibers);
% Loop over fibers:
for fIndex = 1:nFibers
    P_v = zeros(nTheta, 1);
    % Loop over angles:
    for aIndex = 1:length(theta_v)
        theta = theta_v(aIndex);  % In radians.
        % Use only m=0 terms and phi=0:
        P_v(aIndex) = evalP([theta, 0], pl0All_m(:, fIndex).'); % Row vector input.
    end
    P_m(:, fIndex) = P_v;
    % Find the sign of the 2nd derivative at theta = pi/2.
    % Get central difference approximation to derivatives. Note that 1st and
    % last points are invalid. Get 1st derivative:
    %DP_v = (circshift(P_v, -1) - circshift(P_v, 1))/2;  % Divide by dTheta.
    % Now get the 2nd derivative:
    D2P_v = circshift(P_v, -1) - 2*P_v + circshift(P_v, 1); % / dTheta^2.
    signD2At90 = sign(D2P_v(index90));

    % Find the next two sign changes in D2P_v, moving in the -theta direction:
    index = index90-5;
    while (sign(D2P_v(index))*sign(D2P_v(index-1)) == 1)
        index = index - 1;
    end
    indexL1 = index-1;      % First inflection point.
    index = indexL1-5;
    while (sign(D2P_v(index))*sign(D2P_v(index-1)) == 1)
        index = index - 1;
        % For debugging:
        if (index == 1)
            disp(['Error in p00Correction.m: plmPrime_m = ', ...
                num2str(plmPrime_m)])
        end
    end
    indexL2 = index - 1;    % Second inflection point.

    % Search for the extremum between theta_v(indexL1) and theta_v(indexL2).
    % Flip function if necessary (i.e., if extremum is a maximum, flip it over
    % to find the local minimum). Use only m=0 coefficients, so beta is
    % independent of phi (so set phi=0):
    [thetaOptL, ampL0] = fminbnd(@(x) -signD2At90 * ...
        evalP([x, 0], pl0All_m(:, fIndex).'), theta_v(indexL2), theta_v(indexL1));
    % Get elevation angle:
    beta = pi/2 - thetaOptL;
    disp(['For fiber ', num2str(fIndex), ', Estimated beta = ', ...
        num2str(beta*180/pi), ' degrees'])
    sinBeta_v(fIndex) = sin(beta);
    %ampL = -signD2At90*ampL0;

    % Integrate from theta = pi/2 - beta to pi/2 + beta.
    % Loop through orders:
    lSum = zeros(size(beta));
    for l = 0:2:lMax
        % Get pl0:
        pl0Prime = pl0Prime_m(l/2 + 1, fIndex);
        % Get Legendre function values:
        if l == 0
            Pminus_m = zeros(1, length(beta));    % Zero for order < 0.
        else
            Pminus_m = legendre(l-1, sinBeta_v(fIndex));
        end
        Pplus_m = legendre(l+1, sinBeta_v(fIndex));
        lSum = lSum + (pl0Prime/sqrt(2*l+1))*(Pplus_m(1,:) - Pminus_m(1,:));
    end

    % Get integral values and mean values:
    intVal_v(fIndex) = sqrt(4*pi) * lSum;     % I(beta) of (H.25).
    area = 4*pi * sinBeta_v(fIndex);
    meanVal_v(fIndex) = intVal_v(fIndex) / area;    % Weighted mean.
    % Add this value to p00 to force meanVal to be zero in belt region:
    p00Correction_v(fIndex) = -sqrt(4*pi) * meanVal_v(fIndex); % = -lSum/sinBeta.

    % Perform integral numerically, as a check.
    % Integrate the FOD from theta = pi/2 - beta to pi/2 + beta. Note that the
    % integrand must accept (and return) a row vector. Use (H.18). Note that
    % m ~=0 terms integrate to zero, so just include m=0 terms in integral:
    numFlag = 0;    % 1 to show numerical integral, as a check.
    if (numFlag == 1)
        Ibeta_m0 = 4*pi*integral(@(x) sin(x).*evalP2(x.', zeros(size(x.')), ...
            pl0All_m(:, fIndex).').', thetaOptL, pi/2);
        disp(['   Numerical value of I(beta) = ', num2str(Ibeta_m0)])
        disp(['   Mean value of FOD on unit sphere from thetaOptL to PI/2 = ', ...
            num2str(Ibeta_m0/(4*pi*sinBeta_v(fIndex)))])
        disp(['   Analytic valaue of FOD mean = ', num2str(meanVal_v(fIndex))])
    end

end    % End of fIndex FOR.

% Make sure corrections sum to zero:
p00Correction2_v = p00Correction_v - mean(p00Correction_v);

% Update the plm matrix:
plmPrimeOut_m = plmPrime_m;
plmPrimeOut_m(1, :) = plmPrime_m(1, :) + p00Correction2_v;
