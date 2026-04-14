function [xOpt_v, plmPrimeCor_m, solutions_sa] = optFiberFrames(xInit_v, plmLab_v)
% Usage: [xOpt_v, plmPrimeCor_m, solutions_sa] = optFiberFrames(xInit_v, plmLab_v)
%
% optFiberFrames optimizes the fiber frame (FF) orientations, given the
% voxel-level FOD (in plmLab_v) and realistic constraints on the FF FODs.
% optFiberFrames appears to have the best performance using the
% 'multistart' search algorithm.
%
% Input variables:
%
%   xInit_v:    (2*nSubFixels) x 1 array of initial FF orientations. The phi 
%               values for each FF are held in the first nSubFixels elements; 
%               the theta values are in the last nSubFixels elements. All
%               angles are in radians and the orientation of the nth FF is
%               given by 
%                   (theta, phi) = (xInit_v(nSubFixels+n), xInit_v(n)).
%   plmLab_v:   nTerms x 1 array of lab frame (LF) FOD coefficients (in
%               mNNpEvenReal1 format) for given lMax.
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
%
% optFiberFrames calls the functions:
%   xStart_m = getMyCustomStartPoints(xInit_v, alphaMax, nTrials)
%   xOpt_v = rectifyAngles(xOpt_v)
%   [fracFitError, fitCost, penaltyCost, plmLabModel_v, plmPrime_m] = ...
%       splitFodError2(xOpt_v, plmLab_v, testRotFlag0, ' ')
%   [plmPrimeCor_m, p00Correction_v, p00Correction2_v] = ...
%       p00Correction(plmPrime_m)

% Get nSubFixels:
nSubFixels = length(xInit_v)/2;

% Define constants:
testRotFlag0 = 0;   % For optimization.
testRotFlag1 = 1;   % For visualization.

% Define lower-bound 1D array of angle increments:
delta = 35;     % Half-width of search window for each subfixel [degrees].                                                                                                                                          
xlb_v = zeros(2*nSubFixels, 1);
xlb_v(1:nSubFixels) = -(delta*pi/180) * ones(nSubFixels, 1);   % phi.
xlb_v(nSubFixels + (1:nSubFixels)) = -(delta*pi/180) * ones(nSubFixels, 1);   % theta.
% Define upper-bound 1D array:
xub_v = zeros(2*nSubFixels, 1);
xub_v(1:nSubFixels) = (delta*pi/180) * ones(nSubFixels, 1);
xub_v(nSubFixels + (1:nSubFixels)) = (delta*pi/180) * ones(nSubFixels, 1);
% Define function to optimize:
fitError = @(dOrient_v) splitFodErrorWrapper(xInit_v + dOrient_v, plmLab_v, ...
    testRotFlag0, '');  % For fminsearch (simplex), simulanneal, GlobalSearch.
fitErrorMS = @(x_v) splitFodErrorWrapper(x_v, plmLab_v, ...
    testRotFlag0, '');  % For multistart.
fitErrorSO = @(x_v) splitFodErrorWrapper(x_v(:), plmLab_v, ...
    testRotFlag0, '');  % For multistart.
% Find splitFodError at initial search point:
fracFitErrorInit = splitFodErrorWrapper(xInit_v, plmLab_v, testRotFlag0, '');
disp(['   fractional fit error at initial search point = ', ...
    num2str(fracFitErrorInit)])

% Launch search for the best increment to add to xInit_v (for GlobalSearch)
% or best angles (for multistart):
search_s = 'multistart_irls'; % fminsearch, simulanneal, multisart or GlobalSearch.
tStart = tic;
switch search_s
    case 'fminsearch'       % Simplex search.
        %[xOpt_v, fval] = fminsearch(fitError, zeros(size(xInit_v)));
        options_stc = optimset('PlotFcns', 'optimplotfval', ...
            'MaxFunEvals', 1e4, 'TolX', 1e-6, 'TolFun', 1e-6, ...
            'Display', 'off');
        problem_stc.objective = fitError;
        problem_stc.x0 = zeros(size(xInit_v));
        problem_stc.solver = 'fminsearch';
        problem_stc.options = options_stc;
        [xOpt_v, fval] = fminsearch(problem_stc);
        solutions_sa = [];  % Only returned by multistart and GlobalSearch.
        disp('For simplex search...')
    case 'simulanneal'
        options = optimoptions('simulannealbnd', 'PlotFcns',...
          {@saplotbestf});
        [xOpt_v, fval] = simulannealbnd(fitError, xInit_v, xlb_v, xub_v, ...
            options);
        solutions_sa = [];  % Only returned by multistart and GlobalSearch.        
        disp('For simulated annealing search...')
    case 'multistart'
        % Get number of trial points:
        switch nSubFixels
            case 1
                nTrials = 1;
            case 2
                nTrials = 10;
            otherwise
                nTrials = 60;
        end
        alphaMax =  30 * (pi/180);   % 30 degrees (matches simFixelDecomp8.m).
        xStart_m = getMyCustomStartPoints(xInit_v, alphaMax, nTrials);
        tpoints = CustomStartPointSet(xStart_m);    % MATLAB fn.
        % Redefine lb and ub for multistart:
        lb_v = [-pi*ones(nSubFixels, 1); -(pi/2)*ones(nSubFixels, 1)];
        ub_v = [pi*ones(nSubFixels, 1); (pi/2)*ones(nSubFixels, 1)];
        % Create objects:
        disp(['Running multistart with ', num2str(nTrials), ' trials...'])
        %rng default % For reproducibility
        opts = optimoptions(@fmincon, 'UseParallel', false);
        probObj = createOptimProblem('fmincon', 'x0', zeros(size(xInit_v)), ...
               'objective', fitErrorMS, 'lb', lb_v, 'ub', ub_v, 'options', opts);
        msObj = MultiStart('Display', 'iter');
        [xOpt_v, fval, exitflag, output, solutions_sa] = run(msObj, probObj, tpoints);
        xOpt_v = rectifyAngles(xOpt_v); % Insure that angles are rectified
    case 'multistart_irls'
            % Get number of trial points:
            switch nSubFixels
                case 1
                    nTrials = 1;
                case 2
                    nTrials = 10;
                otherwise
                    nTrials = 60;
            end
            alphaMax =  30 * (pi/180);   % 30 degrees (matches simFixelDecomp8.m).
            xStart_m = getMyCustomStartPoints(xInit_v, alphaMax, nTrials);
            tpoints = CustomStartPointSet(xStart_m);    % MATLAB fn.
            % Redefine lb and ub for multistart:
            lb_v = [-pi*ones(nSubFixels, 1); -(pi/2)*ones(nSubFixels, 1)];
            ub_v = [pi*ones(nSubFixels, 1); (pi/2)*ones(nSubFixels, 1)];
            % Create objects:
            disp(['Running multistart with ', num2str(nTrials), ' trials...'])
            %rng default % For reproducibility
            opts = optimoptions(@fmincon, 'UseParallel', false);
            msObj = MultiStart('Display', 'iter');

            maxIRLSIter = 5;
            tolW = 1e-3;
            W = diag(ones(1002,1));
            xCurr_v = xInit_v;
            
            for IRLSIter = 1:maxIRLSIter
                fitErrorIRLS = @(x_v) splitFodErrorWrapperIRLS(x_v, plmLab_v, ...
                    testRotFlag0, '', W);
                probObj = createOptimProblem('fmincon', 'x0', xCurr_v, ...
                   'objective', fitErrorIRLS, 'lb', lb_v, 'ub', ub_v, 'options', opts);

                [xCurr_v, fval, exitflag, output, solutions_sa] = run(msObj, probObj, tpoints);
                [fracFitError, fitCost, penaltyCost, plmLabModel_v, plmPrime_m] = splitFodError2_IRLS(xCurr_v, plmLab_v, testRotFlag0, '', W);

                % Update W based on fit error
                diagonal_W = ones(size(W, 1)) / max(tolW, fitCost);
                W = diag(diagonal_W);
            end
            xOpt_v = rectifyAngles(xOpt_v); % Insure that angles are rectified.
    case 'surrogateopt'
        % Get number of trial points:
        switch nSubFixels
            case 1
                nTrials = 1;
            case 2
                nTrials = 10;
            otherwise
                nTrials = 60;
        end
    
        alphaMax = 30 * (pi/180);   % 30 degrees (matches simFixelDecomp8.m).
        xStart_m = getMyCustomStartPoints(xInit_v, alphaMax, nTrials);
    
        % Finite bounds required by surrogateopt:
        lb_v = [-pi*ones(nSubFixels, 1); -(pi/2)*ones(nSubFixels, 1)];
        ub_v = [ pi*ones(nSubFixels, 1);  (pi/2)*ones(nSubFixels, 1)];
    
        % Make sure the initial points are inside the bounds.
        xStart_m = min(max(xStart_m, lb_v.'), ub_v.');
    
        % Make sure each row is one point, with the correct number of columns.
        if size(xStart_m, 2) ~= numel(lb_v) && size(xStart_m, 1) == numel(lb_v)
            xStart_m = xStart_m.';
        end        
        xStart_m = min(max(xStart_m, lb_v.'), ub_v.');
        
        opts = optimoptions('surrogateopt', ...
            'Display', 'iter', ...
            'InitialPoints', xStart_m, ...
            'MaxFunctionEvaluations', 5000, ...
            'UseParallel', false);
        
        solutions_sa = [];
        [xOpt_v, fval, exitflag, output] = surrogateopt(fitErrorSO, lb_v, ub_v, opts);
        xOpt_v = xOpt_v(:);
        xOpt_v = rectifyAngles(xOpt_v);
    case 'GlobalSearch'
        constrain_s = 'y';  % Select constrained or unconstrained search.
        if strcmp(constrain_s, 'y')
            % Setup global search:
            nTrialPoints0 = 2000;  % was 25^nSubFixels.
            NumStageOnePoints = max([200, 1000]);  % Default = 200.
            NumTrialPoints = max([1500, nTrialPoints0]);    % Default = 1000.
            MaxWaitCycle = 10;              % Default = 20.
            BasinRadiusFactor = 0.5;        % Default = 0.2.
            DistanceThresholdFactor = 0.75; % Default = 0.75.
            PenaltyThresholdFactor = 0.2;   % Default = 0.2.
            % Define GlobalSearch object:
            gsObj = GlobalSearch('Display', 'iter', ...
                'NumStageOnePoints', NumStageOnePoints, ...
                'NumTrialPoints', NumTrialPoints, ...
                'StartPointsToRun', 'bounds', ...
                'MaxWaitCycle', MaxWaitCycle, ...
                'BasinRadiusFactor', BasinRadiusFactor, ...
                'DistanceThresholdFactor', DistanceThresholdFactor, ...
                'PenaltyThresholdFactor', PenaltyThresholdFactor);
            problemObj = createOptimProblem('fmincon', 'x0', zeros(size(xInit_v)), ...
                'objective', fitError, 'lb', xlb_v, 'ub', xub_v);
            disp(['   running constrained global search with NumStageOnePoints = ', ...
                num2str(NumStageOnePoints), ', NumTrialPoints = ', ...
                num2str(NumTrialPoints)])
        else
            % Setup unconstrained global search:
            gsObj = GlobalSearch('Display', 'iter', 'NumTrialPoints', 1000);
            problemObj = createOptimProblem('fmincon', 'x0', zeros(size(xInit_v)), ...
                'objective', fitError);
            disp('   running unconstrained global search...')
        end
        % Run it:
        [xOpt_v, fval, exitflag, output, solutions_sa] = run(gsObj, problemObj);
        disp('For global search...')
end
tEnd = toc(tStart); % tEnd is the elapsed time in seconds, calculated by toc.
disp(['Elapsed time = ', num2str(tEnd/60), ' minutes'])
fracFitErrorOpt = fval;
disp(['   Fractional error at optimized solution = ', ...
    num2str(fracFitErrorOpt)])

%% Correct p00 coefficients.

% Get solution:
if strcmp(search_s, 'multistart')
    [fracFitError, fitCost, penaltyCost, plmLabModel_v, plmPrime_m] = ...
        splitFodError2(xOpt_v, plmLab_v, testRotFlag0, ' ');
elseif strcmp(search_s, 'surrogateopt')
    [fracFitError, fitCost, penaltyCost, plmLabModel_v, plmPrime_m] = ...
        splitFodError2(xOpt_v, plmLab_v, testRotFlag0, ' ');
else
    [fracFitError, fitCost, penaltyCost, plmLabModel_v, plmPrime_m] = ...
        splitFodError2(xInit_v + xOpt_v, plmLab_v, testRotFlag0, ' ');
end

% Correct p00 terms (this does not affect plmLabModel):
[plmPrimeCor_m, p00Correction_v, p00Correction2_v] = p00Correction(plmPrime_m);

% End of main function.
end

%% Local functions.

function fracFitError = splitFodErrorWrapper(ffAxes_v, plmLab_v, testRotFlag, figName_s)
% Usage: fracFitError = splitFodErrorWrapper(ffAxes_v, plmLab_v, testRotFlag, figName_s)
%
% splitFodErrorWrapper calls the function splitFodError2 and returns only
% the first output variable (for fmincon).
%
% Call splitFodError2:
[fracFitError, fitCost, penaltyCost, plmLabModel_v, plmPrime_m] = ...
    splitFodError2(ffAxes_v, plmLab_v, testRotFlag, figName_s);
end

function fracFitError = splitFodErrorWrapperIRLS(ffAxes_v, plmLab_v, testRotFlag, figName_s, W)
% Usage: fracFitError = splitFodErrorWrapper(ffAxes_v, plmLab_v, testRotFlag, figName_s)
%
% splitFodErrorWrapper calls the function splitFodError2 and returns only
% the first output variable (for fmincon).
%
% Call splitFodError2:
[fracFitError, fitCost, penaltyCost, plmLabModel_v, plmPrime_m] = ...
    splitFodError2_IRLS(ffAxes_v, plmLab_v, testRotFlag, figName_s, W);
end

