function [xOpt_v, plmPrimeMrtrix_m, solutions_sa, xInit_v] = ...
    splitFod4_mrtrix(plmLabMrtrix_v)
% Usage: [xOpt_v, plmPrimeCor_m, solutions_sa, xInit_v] = ...
%   splitFod4_mrtrix(plmLabMrtrix_v)
%
% splitFod4_mrtrix estimates the fiber frame (FF) axes using Xin Hong's 
% watershed algorithm. These estimates are then used to initialize the 
% multistart search in optFiberFrames.m. If the fracFitError is not below 
% threshold, then the subFixel FOD with the largest axial asymmetry is 
% split iteratively, until the fracFitError goal is met.
%
% Input variables:
%
%   plmLabMrtrix_v: nTermsMrtrix x 1 array of lab frame (LF) FOD coefficients
%               (in MRtrix3 format) for given lMax.
%
% Output variables:
%
%   xOpt_v:     (2*nSubFixels) x 1 array of optimized FF orientations. The
%               phi values for each FF are held in the first nSubFixels 
%               elements; the theta values are in the last nSubFixels elements. 
%               All angles are in radians and the orientation of the nth FF
%               is given by 
%                   (theta, phi) = (xOpt_v(nSubFixels+n), xOpt_v(n)).
%   plmPrimeMrtrix_m: nTermsMrtrix x nFF array of FF FOD coefficients (in
%               MRtrix3 format) for given lMax. The FOD of the nth
%               fiber is in the nth column.
%   solutions_sa:   Structure array of solutions found by multistart or 
%               GlobalSearch optimization. Fields are X (local solution
%               coordinates in the search space), Fval (value of cost
%               function at the solution), ExitFlag (search exit condition)
%               Output (structure holding information about the search
%               process from fmincon), and X0 (initial point that lead to X).
%   xInit_v:    (2*nSubFixels) x 1 array of initial FF orientations. The phi 
%               values for each FF are held in the first nSubFixels elements; 
%               the theta values are in the last nSubFixels elements. All
%               angles are in radians and the orientation of the nth FF is
%               given by 
%                   (theta, phi) = (xInit_v(nSubFixels+n), xInit_v(n)).
%               These are returned to the calling program for diagnostic
%               purposes. 
%
% splitFod4_mrtrix calls the functions:
%   plmAll_v = mrtrix2standard(plmMRtrix3_v)
%   slm_m = uniTrans2format(slmAll_m, opFormat_s)
%   [xOpt_v, plmPrimeCor_m, solutions_sa, xInit_v] = splitFod4(plmLab_v)
%   lMax = getLmax_mNNpEvenReal1(length(plmLab_v))
%   vuPointsOnSphere(thetaOpt_v, phiOpt_v)

% Convert coefficients to mNNpEvenReal1 format:
plmLabAll_v = mrtrix2standard(plmLabMrtrix_v');    % Complex row vector output.
plmLab_v = uniTrans2format(plmLabAll_v.', 'mNNpEvenReal1');

% Get useful parameters:
lMax = sqrt(length(plmLabAll_v)) - 1;
nTermsMrtrix = length(plmLabMrtrix_v);

% Enforce FOD normalization:
plmLab_v(1) = 1/sqrt(4*pi);     % FOD integrates to 1.

% Call splitFod4:
[xOpt_v, plmPrimeCor_m, solutions_sa, xInit_v] = splitFod4(plmLab_v);

%% Evaluate solution.

% Show FF FODs:
testRotFlag1 = 1;   % Show subFixels.
[fracFitError, fitCost, penaltyCost, plmLabModel_v, plmPrime_m] = ...
    splitFodError2(xOpt_v, plmLab_v, testRotFlag1, 'Optimized Solution');
figure
alphaOp = 1;
az = 45;
el = 30;
subplot(1, 2, 1)
showFad3(uniTrans2standard(plmLab_v, lMax), alphaOp, az, el)
title('Input LF FOD')
subplot(1, 2, 2)
plmLabModelAll_v = uniTrans2standard(plmLabModel_v, lMax);
showFad3(plmLabModelAll_v, alphaOp, az, el)
title('Sum over FF FODs')

% Report optimized fixel properties and convert output plm to MRtrix3 format:
nSubFixels = length(xOpt_v)/2;
thetaOpt_v = zeros(nSubFixels, 1);  % xInit_v + optimized dOrient_v.
phiOpt_v = zeros(nSubFixels, 1);  % xInit_v + optimized dOrient_v.
plmPrimeMrtrix_m = zeros(nTermsMrtrix, nSubFixels);
figure
for fIndex = 1:nSubFixels
    thetaEst = xOpt_v(nSubFixels + fIndex);
    phiEst = xOpt_v(fIndex);
    fEst = sqrt(4*pi)*plmPrimeCor_m(1, fIndex); % Use p00 corrected FODs.   
    % Trap for thetaEst < 0:
    if (thetaEst < 0)
        [thetaEst, phiEst] = rectifyTheta(thetaEst, phiEst);
    end
    % Unwind phi, if necessary:
    if (abs(phiEst) > 2*pi)
        phiEst = phiEst - sign(phiEst)*floor(abs(phiEst)/(2*pi)) * 2*pi;
    end
    % Show results:
    disp(['For fixel ', num2str(fIndex), '...'])
    disp(['   Estimated theta = ', num2str(thetaEst*180/pi), ' degrees'])
    disp(['   Estimated phi = ', num2str(phiEst*180/pi), ' degrees'])
    disp(['   Estimated f = ', num2str(fEst)])
    % Save net (rectified) theta, phi:
    thetaOpt_v(fIndex) = thetaEst;   
    phiOpt_v(fIndex) = phiEst; 
    % Show FF FODs:
    subplot(1, nSubFixels, fIndex)
    plmPrimeCorAll_v = uniTrans2standard(plmPrimeCor_m(:, fIndex), lMax);
    showFad3(plmPrimeCorAll_v, alphaOp, az, el);
    title(['SubFixel ', num2str(fIndex)])
    % Convert to MRtrix3 format:
    plmPrimeMrtrix_m(:, fIndex) = standard2mrtrix(plmPrimeCorAll_v)';
end
% Plot the subFixels on the surface of a sphere:
vuPointsOnSphere(thetaOpt_v, phiOpt_v);     % Opens new figure itself.
title('Optimized fixel orientations')


