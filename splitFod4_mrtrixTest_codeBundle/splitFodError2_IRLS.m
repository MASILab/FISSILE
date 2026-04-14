function [fracFitError, fitCost, penaltyCost, plmLabModel_v, plmPrime_m] = ...
    splitFodError2_IRLS(ffAxes_v, plmLab_v, testRotFlag, figName_s, W)
% Usage: [fracFitError, fitCost, penaltyCost, plmLabModel_v, plmPrime_m] = ...
%    splitFodError2(ffAxes_v, plmLab_v, testRotFlag, figName_s)
%
% splitFodError finds the error implied by the Fiber Frame (FF) orientations
% specified by ffAxes_v, where asymmetric (m ~= 0) terms in the estimated 
% FF FODs are penalized. The hypothesis is that minimizing the asymmetric
% terms will give the true FF axis orientations. 
% 
% Note that rotation by the Euler angles takes the Z axis of the LF to 
% the FF axis, so eulerAlpha is phi and eulerBeta is theta for the FF axis 
% in the lab frame. Using Euler angles, rather than (theta, phi) avoids
% complications accommodating theta<0 values. 
%
% For details, see "Response Functions II", section H in Notability.
%
% Input arguments are:
%
%   ffAxes_v: (2*Nf) x 1 array of FF axis orientations, defined by the
%       (eulerAlpha, eulerBeta) values that transform the lab frame (LF) 
%       to the fiber frame (FF). The eulerAlpha value for the nth FF axis 
%       ffAxes_v(n) and the eulerBeta value is given by ffAxes_v(Nf + n). 
%   plmLab_v: nTerms x 1 column vector of the net LF FOD coefficients in 
%       mNNpEvenReal1 format. 
%   testRotFlag: 1 to show fiber FODs in FF and rotated to LF.
%
% Output variables are:
%
%   sse: sum of squared errors, i.e., the norm of the error vector (= the
%       difference between plmLab_v and the sum of FF plmPrime coefficients
%       transformed to the lab frame) and the penalties for m~=0 terms.        
%
% History:
%   2024/06/10: created from splitSignalError4.m (AWA).
%
% splitFodErrorX calls the functions:
%   [l_v, m_v, isReal_v] = index2lm_mNNpEvenReal(index_v)
%   Xl_m = nnMpEvenRotationMatrix(euler_m, l)
%   slmAll_m = uniTrans2standard(slm_m, lMax)
%   showFad3(sl0All_v, alpha, az, el)
%   yl0_v = getYl0para(lMax)
%   yl0_v = getYl0perp(lMax)
%   Yl0_m = getYl0(lMax, theta_v)
%   n_m = getDirectionsWrapper(n)

%% Setup parameters.

% Get useful constants:
nTerms = length(plmLab_v);  % mNNpEvenReal1 format. 
Nf = length(ffAxes_v)/2;      % Number of fiber axes.
eta2 = 5;    % Regularization (penalty) term.

% Get derived paramters. First, find lMax using 
%   nTerms = lMax.*(lMax+3)/2 + l0Flag
% or
%   lMax^2 + 3*lMax + 2(l0Flag - nTerms) = 0
% which implies
%   lMax = [-3 + sqrt(9 + 8(nTerms-l0Flag))]/2
% for mNNpEvenReal format:
l0Flag = 1;
lMax = (-3 + sqrt(9 + 8*(nTerms-l0Flag)))/2;

% Set flags:
%testRotFlag = 1;            % 1 to check rotation matrix visually.

% Get 3xNf matrix of [alpha, beta, gamma].' values in col n for the nth FF.
% For each FF, eulerGamma = 0:
euler_m = [ffAxes_v(1:Nf).'; ffAxes_v((Nf+1):(2*Nf)).'; zeros(1, Nf)]; 

% Get index vectors for mNNpEvenReal1 format vectors:
[l_v, m_v, isReal_v] = index2lm_mNNpEvenReal(1:nTerms);

%% Calculate rotation matrix, R_m, from fiber axes.

% Get the transformation matrix, R_m. This is constructed in three steps:
% 1). Xlf_m is the 
%           (2l+1) by (2l+1) 
% matrix for a single FF and l that transforms coefficients (mNNpEvenReal
% format for the current order, l, for both input and output) from the nth  
% FF to the lab frame. Xlf_m is calculated by nnMpEvenRotationMatrix.m using 
% Euler angle definitions of the FF axes.
%
% 2). Xf_m is a nTerms x nTerms matrix that transforms the plmPrime
% coefficients for a given FF to the LF. It is constructed from Xlf_m blocks 
% for different l values, appended down columns.
%
% 3). R_m transforms the plmPrime coefficients for all FFs to the LF (in
% mNNpEvenReal1 format). R_m is therefore 
%               nTerms x (Nf*nTerms)  
% It is constructed from Xf_m blocks for different FFs, appended along rows.
%
% Note that euler_m gives the transformation from the LF to the FF. 
% This is the inverse of the transformation from FF to LF, which we need
% here, so take the opposite Euler angles in the opposite order. 
%
% Loop over FFs and even orders:
R_m = zeros(nTerms, Nf*nTerms);
for fIndex = 1:Nf
    % Get transformation matrix for each FF:
    Xf_m = zeros(nTerms, nTerms); 
    euler_v = euler_m(:, fIndex).'; % Row vector: eulerAlpha, eulerBeta, 0.
    for  l = 0:2:lMax
        % Input to nnMpEvenRotationMatrix is a row vector for one FF. Output
        % is (2l+1)x(2l+1) matrix:
        Xlf_m = nnMpEvenRotationMatrix(-fliplr(euler_v), l);
        % Xf_m is block-diagonal (l terms transform to the same l):
        Xf_m(l_v == l, l_v == l) = Xlf_m;
    end
    R_m(:, (fIndex-1)*nTerms + (1:nTerms)) = Xf_m;
end

% Create constraint matrix for each fiber, one row for each m ~= 0: 
nCon = nTerms - (lMax/2 + 1);   % Number of penalty constraints.
con_m = zeros(nCon, nTerms); 
for l = 2:2:lMax
    for m = 1:l
        % Get row and column index of imaginary component with current
        % (l,m) in mNNpEvenReal1 format--see Response Functions II, (H.9)
        % and (H.10):
        imagRow = (l^2)/2 - l + 2*m;    % m=0 terms are neglected.
        imagCol = (l^2)/2 - l/2 + 2*m + 1;
        % Assign values to penalize real and imag components:
        con_m(imagRow-1, imagCol-1) = m*eta2;   % Real term.
        con_m(imagRow, imagCol) = m*eta2;       % Imag term.
    end
end

% Add constraint matrix (along diagonal) in block-diagonal matrix conAll_m:
conAll_m = zeros(Nf*nCon, Nf*nTerms);
for fIndex = 1:Nf
    conAll_m((fIndex-1)*nCon + (1:nCon), (fIndex-1)*nTerms + (1:nTerms)) = ...
        con_m;
end

% Create net transformation matrix (rotations with constraints):
T_m = [R_m; conAll_m];
% Show rank of T_m, if viewing is switched on:
if (testRotFlag == 1)
    disp(['Rank of T_m = ', num2str(rank(T_m))])
end

%% Solve for FF FOD coefficients.

% Assuming T_m * plmPrime_v = plmLab_v, find plmPrime_v:
plmPrime_v = pinv(T_m) * [plmLab_v; zeros(Nf*nCon, 1)];

% Visualize FF FODs, if desired:
if testRotFlag == 1
    % Show rotation for each FF:
    figure
    set(gcf, 'Name', figName_s)
    alpha = 1;
    az = 0;
    el = 0;
    for fIndex = 1:Nf
        % Get (Nf*nTerms) x 1 vector of plmPrime_v values--only the
        % coefficients for FF=fIndex are non-zero: 
        plmPrimeViz_v = zeros(Nf*nTerms, 1);
        % First get row indices for current fiber in plmPrime_v:
        fIndex_v = (fIndex-1)*nTerms + (1:nTerms);
        plmPrimeViz_v(fIndex_v) = plmPrime_v(fIndex_v); % FOD to visualize.
        % Rotate into LF:
        plmLabTest_v = R_m * plmPrimeViz_v;
        subplot(Nf, 2, (fIndex-1)*2 + 1)
        % Display current fiber FOD in FF:
        plmPrimeAll_v = uniTrans2standard(plmPrime_v(fIndex_v), lMax);
        % Show it:
        showFad3(plmPrimeAll_v, alpha, az, el);
        title(['Fiber ', num2str(fIndex), ' FOD in FF'])
        subplot(Nf, 2, (fIndex-1)*2 + 2)
        % Display FOD for current fiber in LF:
        plmLabTestAll_v = uniTrans2standard(plmLabTest_v, lMax);
        showFad3(plmLabTestAll_v, alpha, az, el);
        title(['Fiber ', num2str(fIndex), ' FOD in LF'])
    end     % End of fIndex FOR.
end         % End of testRotFlag IF.

%% Calculate shape penalty.

% Get cost for each FF:
plmPrime_m = reshape(plmPrime_v, nTerms, Nf); % Return one fiber FOD per column.
yl0_0_v = getYl0para(lMax);
yl0_90_v = getYl0perp(lMax);
% % Get Yl0 at theta = 30 degrees (for l=0:2:lMax):
% Yl0_30_v = getYl0(lMax, pi/6);     
% yl0_30_v = zeros(nTerms, 1);
% yl0_30_v(m_v == 0) = Yl0_30_v;     % Column vector in mNNpEvenReal1 format.
% % Get Yl0 at theta = 45 degrees (for l=0:2:lMax):
% Yl0_45_v = getYl0(lMax, pi/4);     % Yl0 at theta = 45 degrees for l=0:2:lMax.
% yl0_45_v = zeros(nTerms, 1);
% yl0_45_v(m_v == 0) = Yl0_45_v;     % Column vector in mNNpEvenReal1 format.
% % Get Yl0 at theta = 60 degrees (for l=0:2:lMax):
% Yl0_60_v = getYl0(lMax, pi/3);     % Yl0 at theta = 60 degrees for l=0:2:lMax.
% yl0_60_v = zeros(nTerms, 1);
% yl0_60_v(m_v == 0) = Yl0_60_v;     % Column vector in mNNpEvenReal1 format.
Yl0_15_v = getYl0(lMax, 15*pi/180);     % Yl0 at theta = 60 degrees for l=0:2:lMax.
yl0_15_v = zeros(nTerms, 1);
yl0_15_v(m_v == 0) = Yl0_15_v;     % Column vector in mNNpEvenReal1 format.
Yl0_80_v = getYl0(lMax, 80*pi/180);     % Yl0 at theta = 60 degrees for l=0:2:lMax.
yl0_80_v = zeros(nTerms, 1);
yl0_80_v(m_v == 0) = Yl0_80_v;     % Column vector in mNNpEvenReal1 format.
shapeCost_v = zeros(1, Nf);
% Loop over FFs:
for fIndex = 1:Nf
    plmFF_v = plmPrime_m(:, fIndex);
    fod0 = plmFF_v.' * yl0_0_v; % P(0,0).
    %fod30 = plmFF_v.' * yl0_30_v; % P(pi/6,0). % Note: ringing makes this
    %problematic.
    %fod45 = plmFF_v.' * yl0_45_v; % P(pi/4,0).
    % fod60 = plmFF_v.' * yl0_60_v; % P(pi/3,0).
    fod90 = plmFF_v.' * yl0_90_v; % P(pi/2,0).
    fod15 = plmFF_v.' * yl0_15_v;
    fod80 = plmFF_v.' * yl0_80_v;
    % Each difference should be < 0; use a RELU:
    % shapeCost_v(fIndex) = max([fod90 - fod60, fod60 - fod0, 0]); 
    % shapeCost_v(fIndex) = sum(max([fod80 - fod15, fod15 - fod0, fod80 - fod0], 0).^2); 
    shapeCost_v(fIndex) = sum(max([fod90 - fod80, fod80 - fod15, fod15 - fod0, fod90 - fod0, fod80 - fod0], 0).^2); 
end
shapeCost = max(shapeCost_v);    % Take max over fibers.       

%% Calculate error in model fit.

% Calculate error and useful output variables. Start with FOD fit accuracy.
% Calculate fit in angular ('ang') or spherical harmonic ('sh') domain:
fitCost_s = 'ang';  
plmLabModel_v = R_m * plmPrime_v;
switch fitCost_s
    case 'sh'
        fitError_v = plmLabModel_v - plmLab_v;
        fitCost = fitError_v'*fitError_v;
    case 'ang'
        % Get modeled and actual FODs at discrete sample points in space:
        nDirs = 1002;
        n_m = getDirectionsWrapper(nDirs);
        theta_v = acos(n_m(:, 3));
        phi_v = atan2(n_m(:, 2), n_m(:, 1));
        plmLabModelAll_v = uniTrans2standard(plmLabModel_v, lMax);
        fodLabModel_v = evalP2(theta_v, phi_v, plmLabModelAll_v.');
        plmLabAll_v = uniTrans2standard(plmLab_v, lMax);
        fodLab_v = evalP2(theta_v, phi_v, plmLabAll_v.');
        % Compare FODs where they are both above the amplitude threshold.
        % This elminates small lobes due to truncation and noise:
        fodThresh = 0.15 * max(fodLab_v);
        angMask_v = (min(fodLab_v, fodLabModel_v) > fodThresh);
        fitError_v = (fodLabModel_v - fodLab_v).*angMask_v;
        % Get fitCost in anguluar space [see Response Functions II notes, 
        % (H.60)]:
        fitCost = 4*pi * diag(W) * (fitError_v'*fitError_v)/nDirs;
end
% Calculate |m| > 0 penalty:
penalty_v = conAll_m * plmPrime_v;
penaltyCost = penalty_v'*penalty_v;
% Get sign(p00) cost; use RELU:
p00SignCost = max(-plmPrime_m(1,:).*(ones(1, Nf) - sign(plmPrime_m(1,:)))/2);

% Penalty for dominant p20:
pl0Prime_m = plmPrime_m(m_v == 0,:);
mean_pl0 = mean(pl0Prime_m(3:end,:)) + 1e-8;
ratio = pl0Prime_m(2,:) ./ mean_pl0;
penalty_p20 = max([0, ratio - 3])^2;

% Get final SSE value. Penalty for |m| > 0 terms is part of plmPrime_v 
% optimization:
sse = sum(fitCost) + 1 * shapeCost + 2 * penalty_p20 + p00SignCost;
% For debugging:
if (testRotFlag == 1)
    disp(['fitCost = ', num2str(sum(fitCost)/sse), '*SSE, shapeCost = ', ...
        num2str(1*shapeCost/sse), '*SSE, penaltyCost = ', ...
        num2str(2*penalty_p20/sse), '*SSE'])
end

% Return fractional error (includes fit and shape cost):
switch fitCost_s
    case 'sh'
        plmNorm2 = plmLab_v' * plmLab_v;
    case 'ang'
        plmNorm2 = plmLabAll_v' * plmLabAll_v;
end
fracFitError = sqrt(sse/plmNorm2);
