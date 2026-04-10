function [xOpt_v, plmPrimeCor_m, solutions_sa, xInit_v] = splitFod4(plmLab_v, fracFitErrorThresh)
% Usage: [xOpt_v, plmPrimeCor_m, solutions_sa, xInit_v] = splitFod4(plmLab_v, fracFitErrorThresh)
%
% splitFod4 estimates the FF axes using Xin Hong's watershed algorithm.
% These estimates are used to initialize the multistart search in
% optFiberFrames.m. If the fracFitError is not below threshold, then
% the subFixel FOD with the largest axial asymmetry is split iteratively, 
% until the fracFitError goal is met.
%
% Input variables:
%
%   plmLab_v:   nTerms x 1 array of lab frame (LF) FOD coefficients (in
%               mNNpEvenReal1 format) for given lMax.
%   fracFitErrorThresh: 1e-6 for simulation, 5e-2 for actual data.
%
%
% Output variables:
%
%   xOpt_v:     (2*nSubFixels) x 1 array of optimized FF orientations. The
%               phi values for each FF are held in the first nSubFixels 
%               elements; the theta values are in the last nSubFixels elements. 
%               All angles are in radians and the orientation of the nth FF
%               is given by 
%                   (theta, phi) = (xInit_v(nSubFixels+n), xInit_v(n)).
%   plmPrimeCor_m: nTerms x nFF array of FF FOD coefficients (in
%               mNNpEvenReal1 format) for given lMax. The FOD of the nth
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
% splitFod4 calls
%   L = getLmax_mNNpEvenReal1(nTerms)
%   plmLabAll_v = uniTrans2standard(plmLab_v, lMax)
%   [pCoeff_m, fodMax_m, ACC, oriFodMax_m, nOriFib] = splitFOD3(plmLabAll_v, ...
%       realFlag, volThreshold, angThreshold, plotFlag, flipFlag)
%   [xOpt_v, plmPrimeCor_m, solutions_sa] = optFiberFrames(xInit_v, plmLab_v)
%   Q = axialAsymIndex(plm_v, lMax)

%% Get initial FF orientation estimates.

% Convert plmLab_v to standard format for splitFOD3.m:
lMax = getLmax_mNNpEvenReal1(length(plmLab_v));
plmLabAll_v = uniTrans2standard(plmLab_v, lMax);

% Split FOD peaks:
realFlag = 0;
volThreshold = 0.15;    % Large enough to reject Gibbs ringing.
angThreshold = 30*(pi/180); 
plotFlag = 0;   % 1 to show plots for debugging, 0 to skip.
flipFlag = 1;
% Call Xin Hong's watershed algorithm (note capitalization of 'FOD'):
[pCoeff_m, fodMax_m, ACC, oriFodMax_m, nOriFib] = splitFOD3(plmLabAll_v, ...
    realFlag, volThreshold, angThreshold, plotFlag, flipFlag);
% output:   pCoeff_m -- each col gives the SH coefficients of a component. 
%                      Components are listed in the descending order of 
%                      the magnitude.
%           fodMax_m -- each row is the FOD maximum of a component,
%                       [magnitude, theta, phi]. Rows are listed in the
%                       descending order of the magnitude.
%           ACC -- anguar correlation coefficient btwn the original and
%                  reconstructed FOD
%           oriFodMax_m -- FOD maxima matrix of the original FOD. each row
%                          is a maximum: [magnitude, theta, phi]. rows are
%                          listed in the descending order of the magnitude.
%           nOriFib -- # of peaks of the original FOD

% Get orientation of each fiber frame axis:
thetaInit_v = fodMax_m(:, 2);              % Column vector.
phiInit_v = fodMax_m(:, 3);
xInit_v = [phiInit_v; thetaInit_v];

%% Loop until error in FOD fit falls below threshold.

% Top of loop:
%Qthresh = 1e-6;
f_v = fodMax_m(:, 1);   % Initialize estimated volume fractions.
testRotFlag0 = 0;       % Don't test rotations.
splitFlag = 1;
while (splitFlag == 1)
    % Show initial orientations:
    nSubFixels = length(xInit_v)/2;
    for fIndex = 1:nSubFixels
        % Show initial point for each subFixel:
        disp('Initial point for search:')
        disp(['   For subFixel', num2str(fIndex), ', magnitude = ', ...
            num2str(f_v(fIndex))])
        disp(['   For subFixel', num2str(fIndex), ', thetaInit = ', ...
            num2str(xInit_v(nSubFixels + fIndex)*180/pi), ' degrees'])
        disp(['   For subFixel', num2str(fIndex), ', phiInit = ', ...
            num2str(xInit_v(fIndex)*180/pi), ' degrees'])
    end

    % Use multistart to search for more accurate estimates of FF orientation:
    [xOpt_v, plmPrimeCor_m, solutions_sa] = optFiberFrames(xInit_v, plmLab_v);
    % Get value of cost function:
    [fracFitError, fitCost, penaltyCost, plmLabModel_v, plmPrime_m] = ...
        splitFodError2(xOpt_v, plmLab_v, testRotFlag0, '');
    disp(['For optimized solution, fracFitError = ', num2str(fracFitError)])

    % Check for asymmetric FF FODs:
    Q_v = zeros(1, nSubFixels);
    for fIndex = 1:nSubFixels
        plmPrimeCor_v = plmPrimeCor_m(:, fIndex);
        Q = axialAsymIndex(plmPrimeCor_v, lMax);
        disp(['For subFixel ', num2str(fIndex), ', Q = ', num2str(Q)])
        Q_v(fIndex) = Q;
    end
    Qmax = max(Q_v);

    % Determine whether fit is acceptable. If not, add another subFixel at 
    % the orientation of the fiber with the greatest FOD axialAsymmetry:
    fprintf('fracFitError: %.4f\n', fracFitError);
    fprintf('fracFitErrorThresh: %.4f\n', fracFitErrorThresh);
    if (fracFitError > fracFitErrorThresh)
        % Split FOD with largest axialAsymmetry:
        QmaxIndex = find(Q_v == Qmax, 1, 'first');
        phiQmax = xOpt_v(QmaxIndex);
        thetaQmax = xOpt_v(nSubFixels + QmaxIndex);
        % Add a new FF; append a copy of the current FF axis to xInit_v:
        xInit_v = [xOpt_v(1:nSubFixels); phiQmax; ...
            xOpt_v((nSubFixels+1):(2*nSubFixels)); thetaQmax];
        disp(['  duplicating FF #', num2str(QmaxIndex), ' in xInit_v...'])
        % Divide volume fraction:
        f_v(QmaxIndex) = f_v(QmaxIndex)/2;
        f_v = [f_v; f_v(QmaxIndex)];
        disp(['Repeating optFiberFrames loop with ', ...
            num2str(length(f_v)), ' subFixels...'])
    else
        splitFlag = 0;  % Exit loop.
        disp(['Exiting optFiberFrames loop with ', ...
            num2str(length(f_v)), ' subFixels...'])
    end

end % of splitFlag WHILE.

%% End of main function.
end
