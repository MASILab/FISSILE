% splitFod4_mrtrixTest is a script that tests splitFod4_mrtrix using either
%   1. data from an MRtrix3 FOD .mif file, or
%   2. simulated data.
%
% splitFod4_mrtrixTest calls the functions:
%   fod_4d = readMrtrixFod(filePath_s)
%   plmMrtrix3_v = simFodMrtrix(theta_v, phi_v, f_v, vuFodFlag)
%   [xOpt_v, plmPrimeCor_m, solutions_sa, xInit_v] = ...
%       splitFod4_mrtrix(plmMrtrix3_v)

% Set parameters:
source_s = 'sim';   % 'sim' for simulation or 'mif' for MRtrix3 FOD file.
vuFodFlag = 1;      % 1 to view FOD.

%% Get data.

% Select FOD:
switch source_s
    case 'mif'
        % Define MRtrix3 FOD file:
        dataDir_s = '/Users/andera/Data/ChristosNHP/Viking2019_02_02_mrtrix';
        disp(['Reading data from ', dataDir_s, filesep, file_s, '...'])
        file_s = 'seriesAfod.mif';

        % Read FOD data from file:
        fod_4d = readMrtrixFod([dataDir_s, filesep, file_s]);

        % Identify voxel to send to splitFod:
        nFibersSaved = input('Enter # of desired fibers (0 to explore voxels): ');
        switch nFibersSaved
            case 1
                row = 1;
                col = 1;
                slice = 1;
            case 2
                row = 41;
                col = 40;
                slice = 18;
            case 3
                row = 1;
                col = 1;
                slice = 1;
            otherwise
                contFlag = 1;
                figure
                % Set initial slice:
                slice = round(size(fod_4d, 3)/2);
                while (contFlag == 1)
                    % Show p00_m map for current slice:
                    subplot(1, 2, 1)
                    p00_m = squeeze(fod_4d(:, :, slice, 1));
                    imagesc(p00_m)
                    colormap("gray")
                    axis image
                    title(['Slice ', num2str(slice), ' p00 map'])
                    % Prompt for voxel of interest:
                    disp('Click on VOI...')
                    [x, y] = ginput(1);
                    row = round(y);
                    col = round(x);
                    color_3d = cat(3, p00_m, p00_m, p00_m)/max(p00_m(:));
                    color_3d(row, col, :) = [1; 0; 0];
                    image(color_3d)
                    axis image
                    title(['Slice ', num2str(slice), ' p00 map'])
                    % Show FOD:
                    subplot(1, 2, 2)
                    plmMrtrix3_v = squeeze(fod_4d(row, col, slice, :));
                    plmAll_v = mrtrix2standard(plmMrtrix3_v);
                    showFad3(plmAll_v, 1, 45, 30);
                    % Get next voxel location:
                    slice = input('Enter slice number (or 0 to use current voxel): ');
                    % Test for exit:
                    if (slice == 0)
                        contFlag = 0;
                    end
                end % of contFlag WHILE.
        end     % of nFibersSaved SWITCH.
        % If nFibersSaved is not zero, show location and FOD:
        if (nFibersSaved > 0)
            figure
            % Show p00_m map for current slice:
            subplot(1, 2, 1)
            p00_m = squeeze(fod_4d(:, :, slice, 1));
            color_3d = cat(3, p00_m, p00_m, p00_m)/max(p00_m(:));
            color_3d(row, col, :) = [1; 0; 0];
            image(color_3d)
            axis image
            title(['Slice ', num2str(slice), ' p00 map'])
            % Show FOD:
            subplot(1, 2, 2)
            plmMrtrix3_v = squeeze(fod_4d(row, col, slice, :));
            plmAll_v = mrtrix2standard(plmMrtrix3_v);
            showFad3(plmAll_v, 1, 45, 30);
        end
    case 'sim'
        % Define simulated fixel orientations and amplitudes. Due to long
        % run-times, keep number of fixels <= 3:
        disp('Simulating data...')
        nFixelsSim = 3;
        switch nFixelsSim
            case 1
                theta_v = 20*pi/180;   % Polar angle of each fixel.
                phi_v = 45*pi/180;     % Azimuthal angle of each fixel.
            case 2
                theta_v = [0; 20]*pi/180;   % Polar angle of each fixel.
                phi_v = [0; 0]*pi/180;      % Azimuthal angle of each fixel.
            case 3
                theta_v = [0; 20; 60]*pi/180;   % Polar angle of each fixel.
                phi_v = [0; 0; 45]*pi/180;      % Azimuthal angle of each fixel.
        end
        f_v = ones(length(theta_v), 1)/length(theta_v); % Vol. fraction.
        plmMrtrix3_v = simFodMrtrix(theta_v, phi_v, f_v, vuFodFlag);
end

%% Split FOD into fixels.

% Send to splitFod4_mrtrix:
[xOpt_v, plmPrimeCor_m, solutions_sa, xInit_v] = splitFod4_mrtrix(plmMrtrix3_v);
