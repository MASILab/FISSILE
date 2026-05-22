function [xOpt_v, plmPrimeMrtrix_m, solutions_sa, xInit_v] = ...
    splitFod4_mrtrix(plmLabMrtrix_v, fracFitErrorThresh, eta2, max_fibers)
% Usage: [xOpt_v, plmPrimeCor_m, solutions_sa, xInit_v] = ...
%   splitFod4_mrtrix(plmLabMrtrix_v, fracFitErrorThresh, eta2, max_fibers)
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
%   fracFitErrorThresh: 1e-6 for simulation, 5e-2 for actual data.
%   eta2: Axial symmetry constraint. 5 by default.
%   max_fibers: maximum number of fibers to fit. 10 by default.
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

if nargin < 3
    eta2 = 5;
end
if nargin < 4
    max_fibers = 10;
end

% Convert coefficients to mNNpEvenReal1 format:
plmLabAll_v = mrtrix2standard(plmLabMrtrix_v');    % Complex row vector output.
plmLab_v = uniTrans2format(plmLabAll_v.', 'mNNpEvenReal1');

% Get useful parameters:
lMax = sqrt(length(plmLabAll_v)) - 1;
nTermsMrtrix = length(plmLabMrtrix_v);

% Enforce FOD normalization:
plmLab_v = plmLab_v / plmLab_v(1) * 1 / sqrt(4*pi); % FOD integrates to 1.

% Call splitFod4:
[xOpt_v, plmPrimeCor_m, solutions_sa, xInit_v] = splitFod4(plmLab_v, fracFitErrorThresh, eta2, max_fibers);

%% Evaluate solution.

% Show FF FODs:
testRotFlag1 = 0;   % Show subFixels.
[fracFitError, fitCost, penaltyCost, plmLabModel_v, plmPrime_m] = ...
    splitFodError2(xOpt_v, plmLab_v, testRotFlag1, 'Optimized Solution', eta2);
% figure
% alphaOp = 1;
% az = 45;
% el = 30;
% subplot(1, 2, 1)
% showFad3(uniTrans2standard(plmLab_v, lMax), alphaOp, az, el)
% title('Input LF FOD')
% subplot(1, 2, 2)
% plmLabModelAll_v = uniTrans2standard(plmLabModel_v, lMax);
% showFad3(plmLabModelAll_v, alphaOp, az, el)
% title('Sum over FF FODs')

% Report optimized fixel properties and convert output plm to MRtrix3 format:
nSubFixels = length(xOpt_v)/2;
thetaOpt_v = zeros(nSubFixels, 1);  % xInit_v + optimized dOrient_v.
phiOpt_v = zeros(nSubFixels, 1);  % xInit_v + optimized dOrient_v.
plmPrimeMrtrix_m = zeros(nTermsMrtrix, nSubFixels);
plmPrimeMrtrix_lf_m = zeros(nTermsMrtrix, nSubFixels);
% figure
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
    % % Show FF FODs:
    % subplot(1, nSubFixels, fIndex)
    plmPrimeCorAll_v = uniTrans2standard(plmPrimeCor_m(:, fIndex), lMax);
    % showFad3(plmPrimeCorAll_v, alphaOp, az, el);
    % title(['SubFixel ', num2str(fIndex)])
    % Convert to MRtrix3 format:
    plmPrimeMrtrix_m(:, fIndex) = standard2mrtrix(plmPrimeCorAll_v)';
end
% 
% figure
% for fIndex = 1:nSubFixels
%     subplot(1, nSubFixels+1, fIndex)
%     showFad3(plmPrimeStandard_lf_m(:, fIndex), 1, 50, 25);
%     xlim([-0.3, 0.3])
%     ylim([-0.3, 0.3])
%     zlim([-0.3, 0.3])
% end
% subplot(1, nSubFixels+1,  nSubFixels+1)
% showFad3(sum(plmPrimeStandard_lf_m, 2), 1, 50, 25);
% xlim([-0.3, 0.3])
% ylim([-0.3, 0.3])
% zlim([-0.3, 0.3])

% Plot the subFixels on the surface of a sphere:
% vuPointsOnSphere(thetaOpt_v, phiOpt_v);     % Opens new figure itself.
% title('Optimized fixel orientations')


