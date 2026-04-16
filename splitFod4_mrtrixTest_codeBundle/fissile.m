function fissile(input_file, output_file, fracFitErrorThresh, n_jobs)
%FISSILE Run FISSILE to separate crossing fibers
% FISSILE(input_file, output_file, fracFitErrorThresh, n_jobs) run FISSILE
% on a mat file called input_file containing total_fod, which should
% be n_voxel x n_SH in mrtrix/tournier07 basis. The output file will be a 
% mat file with n_voxel structs containing total_fod, separate_fod, 
% est_theta, est_phi, and est_v. Note that separate_fod is in the frame of 
% the fiber and needs to be rotated to align with voxel space. 
% fracFitErrorThresh determines when to split fibers, 5e-2 for in vivo data
% and 5e-6 for simulations. n_jobs is the number of parfor processes.
%
% Example:
%   matlab -batch "run_fissile('fissile_input.mat', 'fissile_output.mat',
%   5e-2, 10)"

    %% Test sensitivity of FISSILE to volume fraction and angle between fibers
    % addpath(genpath('./splitFod4_mrtrixTest_codeBundle'));
    % addpath(genpath('./uniform_sphere_sampling'));
    % addpath(genpath('./dirichlet-matlab'));


    if ischar(fracFitErrorThresh) || isstring(fracFitErrorThresh)
        fracFitErrorThresh = str2double(fracFitErrorThresh);
    end

    if ischar(n_jobs) || isstring(n_jobs)
        n_jobs = str2num(n_jobs);
    end

    S = load(input_file);
    fprintf('input_file: %s\n', input_file);
    fprintf('output_file: %s\n', output_file);
    fprintf('fracFitErrorThresh: %.4f\n', fracFitErrorThresh);
    fprintf('n_jobs: %d\n', n_jobs);
    
    if ~isfield(S, 'signal')
        error('The input .mat file does not contain a variable named "signal".')
    end
    signal = S.signal;
    n_iter = size(signal, 1);
    data = repmat(struct( ...
        'total_fod', [], ...
        'separate_fod', [], ...
        'rotate_separate_fod', [], ...
        'est_theta', [], ...
        'est_phi', [], ...
        'est_v', []), n_iter, 1);
    
    % Run FISSILE on signal
    total_fod_c = cell(n_iter, 1);
    separate_fod_c = cell(n_iter, 1);
    rotate_separate_fod_c = cell(n_iter, 1);
    est_theta_c = cell(n_iter, 1);
    est_phi_c = cell(n_iter, 1);
    est_v_c = cell(n_iter, 1);

    n_jobs = min(n_jobs, n_iter);
    fprintf('n_voxels: %d\n', n_iter);
    fprintf('Using %d parallel processes\n', n_jobs);
    parpool("Processes", n_jobs);
    parfor i = 1:n_iter
        fprintf('Progress: %3d%% (%d/%d)\n', round(100*(i-1)/size(signal,1)), i, size(signal,1));
    
        % Split into fibers
        [xOpt_v, plmPrimeCor_m, plmPrimeCor_separate_m, solutions_sa, xInit_v] = splitFod4_mrtrix(signal(i,:), fracFitErrorThresh);
    
        total_fod_c{i} = signal(i,:);
        separate_fod_c{i} = plmPrimeCor_m';
        rotate_separate_fod_c{i} = plmPrimeCor_separate_m';
        n_subfixel = size(plmPrimeCor_m, 2);
        est_theta_c{i} = xOpt_v(n_subfixel+1:end)';
        est_phi_c{i} = xOpt_v(1:n_subfixel)';
        est_v_c{i} = sqrt(4*pi)*plmPrimeCor_m(1,:);
    end
    
    data = struct;
    for i = 1:n_iter
        data(i).total_fod = total_fod_c{i};
        data(i).separate_fod = separate_fod_c{i};
        data(i).est_theta = est_theta_c{i};
        data(i).est_phi = est_phi_c{i};
        data(i).est_v = est_v_c{i};
        data(i).rotate_separate_fod = rotate_separate_fod_c{i};
    end

    % Replace "data" with "results" and save
    save(output_file, "data")
end
