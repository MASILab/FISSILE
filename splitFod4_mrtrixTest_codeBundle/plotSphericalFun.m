function [fig, surf_h] = plotSphericalFun(SHCoeff_v, realFlag, zeroFlag, minFlag, title_s, ...
                                fig, subRow, subCol, subId, sh_3d, nTheta, nPhi, ...
                                sizeScale, contourFlag, colorScheme_v, location_v, alphaValue, normalFlag)
% usage: plotSphericalFun(SHCoeff_v, realFlag, zeroFlag, minFlag, title_s, fig, subRow, subCol, subId, sh_3d, ...
%                         nTheta, nPhi, sizeScale, contourFlag, colorScheme_v, location_v, normalFlag)
%        (realFlag... optional)
% plot spherical function surface (eg, FAD, ODF)
% colorScheme_v: [2,-1,3] means xPat_m = yImage_m; yPat_m = -xImage_m;
%                zPat_m = zImage_m
%                [0 0 0] means gray color
% alphaValue: 1-opaque, 0-clear(transparent)

if ~exist('realFlag') | isempty(realFlag)
    realFlag = 0;
end
if ~exist('zeroFlag') || isempty(zeroFlag)
    zeroFlag = 1;
end
if ~exist('minFlag') || isempty(minFlag)
    minFlag = 0;
end
if ~exist('title_s')
    title_s = '';
end
if ~exist('sh_3d')
    sh_3d = [];
end
if ~exist('nTheta') || isempty(sh_3d)
    nTheta = 40;	    % For spherical plots. 80-2.2785 degrees
    nPhi = 40;		    % For spherical plots.
end

if ~exist('sizeScale') || isempty(sizeScale)
    sizeScale = 1;
end

if ~exist('contourFlag') || isempty(contourFlag)
    contourFlag = 0;
end

if ~exist('location_v') || isempty(location_v)
    location_v = zeros(3,1);
end

if ~exist('alphaValue') || isempty(alphaValue)
    alphaValue = 1;
end

if ~exist('normalFlag') || isempty(normalFlag)
    normalFlag = 1;
end


[a,b] = size(SHCoeff_v);
if a==1 || b==1 % coefficent vector
    gridP_m = shCoeffToGridFun(SHCoeff_v, realFlag, zeroFlag, minFlag, sh_3d, nTheta, nPhi);
else    % input is acturally gridP_m
    gridP_m = SHCoeff_v;
    nTheta = a;
    nPhi = b;
    if zeroFlag == 1
        gridP_m(gridP_m<0) = 0;
    end
    if minFlag == 1
        gridP_m = gridP_m - min(gridP_m(:));
    end
end

nTheta = size(gridP_m,1);
nPhi = size(gridP_m,2);
gridTheta_v = pi * (0:(nTheta-1)).' / (nTheta-1);
gridPhi_v = 2*pi * (0:(nPhi-1)).' / (nPhi-1);
[gridPhi_m, gridTheta_m] = meshgrid(gridPhi_v.', gridTheta_v);

if exist('fig', 'var')
    if fig ~= 0
        figure(fig)
    else
        fig = figure;
    end
else
    fig = figure;
end
if exist('subRow') & ~isempty(subRow) & subRow~=0
    subplot(subRow, subCol, subId)
end

maxGridP = max(gridP_m(:));
if (maxGridP > 0)
    if normalFlag == 1
        r_m = sizeScale * gridP_m / maxGridP;
    else
        r_m = sizeScale * gridP_m;
    end
else
    r_m = zeros(size(gridP_m));
end

xImage_m = r_m .* sin(gridTheta_m) .* cos(gridPhi_m);
yImage_m = r_m .* sin(gridTheta_m) .* sin(gridPhi_m);
zImage_m = r_m .* cos(gridTheta_m);

if ~exist('colorScheme_v') | isempty(colorScheme_v)
    xPat_m = xImage_m;
    yPat_m = yImage_m;
    zPat_m = zImage_m;
elseif ~any(colorScheme_v)
    xPat_m = 0.2*ones(size(xImage_m));
    yPat_m = 0.2*ones(size(yImage_m));
    zPat_m = 0.2*ones(size(zImage_m));
else
    % Translate from patient to image coordinates:
    for k = 1:3
        switch k
             case 1
                P_s = 'x';
            case 2
                P_s = 'y';
            case 3
                P_s = 'z';
        end
        if colorScheme_v(k) > 0
            sign_s = '';
        else
            sign_s = '-';
        end
        switch abs(colorScheme_v(k))
            case 1
                I_s = 'x';
            case 2
                I_s = 'y';
            case 3
                I_s = 'z';
        end
        eval([P_s, 'Pat_m = ', sign_s, I_s, 'Image_m;'])        
    end
end   

% Plot it. Set color to radius:
color_3d = cat(3, abs(xPat_m), abs(yPat_m), abs(zPat_m));
color_3d = color_3d / max(color_3d(:)) * 2;

if contourFlag == 1
    surf_h = surfc(location_v(1) + xImage_m, location_v(2) + yImage_m, location_v(3) + zImage_m, color_3d, 'FaceAlpha',alphaValue);
else
    surf_h = surf(location_v(1) + xImage_m, location_v(2) + yImage_m, location_v(3) + zImage_m, color_3d, 'FaceAlpha',alphaValue);
end
if ~isempty(title_s)
    title(title_s)
end
axis equal
% alpha(alphaValue)
