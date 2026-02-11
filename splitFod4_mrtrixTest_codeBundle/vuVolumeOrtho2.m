function [xzView, yzView, xyView] = vuVolumeOrtho2(image_3d, waitForUserFlag)
% Usage: [xzView, yzView, xyView] = vuVolumeOrtho2(image_3d, waitForUserFlag)
%
% vuVolumeOrthoX displays orthogonal views of the 3D image volume given
% in image_3d. Note that y = row, x = col, z = slice. If waitForUserFlag=0,
% then the default orthogonal views are displayed and the program exits.
% 
% The image is scaled so the <normPercentile>-th percentile voxel intensity
% and above is mapped to that single value (to prevent the fat signal, for
% example, from making all the other brain tissue too dark).
% 
% On exit, the program returns the [row, col, slice] coordinates of the
% last voxel selected.

% Preliminaries:
[nRows, nCols, nSlices] = size(image_3d);
normPercentile = 98;    % All voxels higher than this are rendered as 1.

% Define flags:
crossHairFlag = 1;  % 1 to display cross hairs.

% Scale image to normPercentile magnitude:
image_v = image_3d(:);
[sortImage_v, sortIndex_v] = sort(image_v);
nVoxels = nRows*nCols*nSlices;
normVal = sortImage_v(round(nVoxels*normPercentile/100));
lessThanIndex_v = find(image_v <= normVal); % Index in image_v.
greaterThanIndex_v = find(image_v > normVal); % Slower, but straightforward.
normImage_3d = zeros(size(image_3d));
normImage_3d(lessThanIndex_v) = image_v(lessThanIndex_v);
normImage_3d(greaterThanIndex_v) = normVal*ones(size(greaterThanIndex_v));
% Replace image_3d with normImage_3d:
image_3d = normImage_3d;
imageMax = max(image_3d(:));
imageMin = min(image_3d(:));

% Show orthogonal views:
figure
% Make figure name = input argument name (does not work for expressions):
set(gcf, 'Name', inputname(1))  
xyView = round(nSlices/2);    % Index for xy view.
xzView = round(nRows/2);
yzView = round(nCols/2);
contFlag = 1;
position_ca = {};
while (contFlag == 1)
    % Show xy view:
    xyImage_m = squeeze(image_3d(:, :, xyView));
    xyAxes_h = subplot(2, 2, 1);
    % Display image--make scaling consistent across views:
    imagesc(xyImage_m, [imageMin, imageMax])    
    colormap(gray)
    axis image
    title(['XY view = ', num2str(xyView)])
    % Add cross hairs:
    if (crossHairFlag == 1)
        % Each column in cross hair coordinate matrices is a line. Specify
        % horizontal line in column 1:
        xyX_m = [1,     yzView; ...
            nCols, yzView];
        xyY_m = [xzView, 1; ...
            xzView, nRows];
        hold on
        line(xyX_m, xyY_m, 'Color', 'y')
        hold off
    end
    % Show xz view:
    xzImage_m = squeeze(image_3d(xzView, :, :)).';
    xzAxes_h = subplot(2, 2, 3);
    imagesc(xzImage_m, [imageMin, imageMax])
    colormap(gray)
    axis image
    title(['XZ view = ', num2str(xzView)])
    % Add cross hairs:
    if (crossHairFlag == 1)
        % Each column in cross hair coordinate matrices is a line. Specify
        % horizontal line in column 1:
        xzX_m = [1,     yzView; ...
            nCols, yzView];
        xzY_m = [xyView, 1; ...
            xyView, nSlices];
        hold on
        line(xzX_m, xzY_m, 'Color', 'y')
        hold off
    end
    % Show yz view:
    yzImage_m = squeeze(image_3d(:, yzView, :));
    yzAxes_h = subplot(2, 2, 2);
    imagesc(yzImage_m, [imageMin, imageMax])
    colormap(gray)
    axis image
    title(['YZ view = ', num2str(yzView)])
    % Add cross hairs:
    if (crossHairFlag == 1)
        % Each column in cross hair coordinate matrices is a line. Specify 
        % horizontal line in column 1:
        yzX_m = [1,       xyView; ...
            nSlices, xyView];
        yzY_m = [xzView, 1; ...
            xzView, nRows];
        hold on
        line(yzX_m, yzY_m, 'Color', 'y')
        hold off
    end
    % Test for no-user-input:
    if (waitForUserFlag == 0)
        break
    end
    % Prompt for input:
    [x, y] = ginput(1);
    currentAxes_h = gca;
    % Determine the view that was clicked on:
    switch currentAxes_h
        case xyAxes_h
            % User clicked on xyView. Update coordinates:
            xzView = round(y);
            yzView = round(x);
            % Test for exit:
            if (xzView < 1 || xzView > nRows || ...
                    yzView < 1 || yzView > nCols)
                contFlag = 0;
            end
        case xzAxes_h
            % User clicked on xzView. Update coordinates:
            xyView = round(y);
            yzView = round(x);
            % Test for exit:
            if (xyView < 1 || xyView > nSlices || ...
                    yzView < 1 || yzView > nCols)
                contFlag = 0;
            end
        case yzAxes_h
            % User clicked on yzView. Update coordinates:
            xyView = round(x);
            xzView = round(y);
             % Test for exit:
            if (xyView < 1 || xyView > nSlices || ...
                    xzView < 1 || xzView > nRows)
                contFlag = 0;
            end
        otherwise
            % Exit loop:
            contFlag = 0;
    end
    % Display coordinates and value:
    if (contFlag == 1)  % Display only if coordinates are legal.
        position_s = ['Row = ', num2str(xzView), ...
            ', Col = ', num2str(yzView), ', Slice = ', num2str(xyView)];
        value_s = ['  Value = ', num2str(image_3d(xzView, yzView, xyView))];
        subplot(2, 2, 4)
        cla
        text('String', position_s, 'FontSize', 8, 'Units', 'normalized', ...
            'Position', [0, 0.9])
        text('String', value_s, 'FontSize', 8, 'Units', 'normalized', ...
            'Position', [0, 0.8])
        text('String', 'Click outside subplot axes to exit...', ...
            'FontSize', 8, 'Units', 'normalized', ...
            'Position', [0, 0.7])
        axis off
    end
end