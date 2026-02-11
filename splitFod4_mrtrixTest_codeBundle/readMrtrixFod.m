function fod_4d = readMrtrixFod(filename_s)
% Usage: fod_4d = readMrtrixFod(filename_s)
%
% readMrtrixFod reads an FOD .mif file into the 4D array fod_4d.
%
% Input:
%   
%   filename_s: 1D character array giving the path and filename of the FOD
%               .mif file calculated by MRtrix3.
%
% Output:
%
%   fod_4d: 4D array with indices [row, col, slice, SH index]. The
%           spherical harmonic indices are for the real spherical
%           harmonics, as defined in the MRtrix documentation.
%
% readMrtrixFod calls the funcitons:
%   image = read_mrtrix (filename)
%   [xzView, yzView, xyView] = vuVolumeOrtho2(image_3d, waitForUserFlag)

% Define flags:
vuVolFlag = 0;  % 1 to view p00 values.
vuFodFlag = 0;  % 1 to view FOD of selected voxels.

% Read data:
image_stc = read_mrtrix(filename_s);
fod_4d = image_stc.data;

% View 3D volume of p00 values:
if (vuVolFlag == 1)
    waitForUserFlag = 1;
    p00_3d = squeeze(fod_4d(:, :, :, 1));
    [xzView, yzView, xyView] = vuVolumeOrtho2(p00_3d, waitForUserFlag);
end

% Convert to standard (complex) format:
if (vuFodFlag == 1)
    contFlag = 1;
    figure
    % Set initial slice:
    slice = round(size(fod_4d, 3)/2);
    while (contFlag == 1)
        % Show p00_m map for current slice:
        subplot(1, 2, 1)
        imagesc(squeeze(fod_4d(:, :, slice, 1)))
        colormap("gray")
        axis image
        title(['Slice ', num2str(slice), ' p00 map'])
        % Prompt for voxel of interest:
        disp('Click on VOI...')
        [x, y] = ginput(1);
        row = round(y);
        col = round(x);
        % Show FOD:
        subplot(1, 2, 2)
        plmMrtrix3_v = squeeze(fod_4d(row, col, slice, :));
        plmAll_v = mrtrix2standard(plmMrtrix3_v);
        showFad3(plmAll_v, 1, 45, 30);
        % Get next voxel location:
        slice = input('Enter slice number (or 0 to exit): ');
        % Test for exit:
        if (slice == 0)
            contFlag = 0;
        end
    end
end