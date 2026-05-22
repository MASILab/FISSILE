function showFad3(plm_v, alpha, az, el)
% Usage: showFad3(plm_v, alpha, az, el)
%
% showFod3 draws the FOD described by the spherical harmonic coefficients 
% in plm_v with opacity alpha, viewed from camera position (az, el).

% Parameters:

% Get dimensions:
nTerms = length(plm_v);

% Derived parameters:
maxOrder = sqrt(nTerms) - 1;

% Find orientation angles:
nTheta = 60;
nPhi = 60;
gridTheta_v = pi * (0:(nTheta-1)).' / (nTheta-1);
gridPhi_v = 2*pi * (0:(nPhi-1)).' / (nPhi-1);
[gridPhi_m, gridTheta_m] = meshgrid(gridPhi_v.', gridTheta_v);

% Find the values of the spherical harmonics at the regularly-spaced
% {gridTheta, gridPhi} points:
sh_3d = [];
for n = 0:maxOrder
    legendre_m = legendre(n, cos(gridTheta_v.')); % .' for n=0 case.
    % Calculate -ve m polynomials from +ve index values:
    for m = -n:n
        if (m < 0)
            negFactor = (-1)^m * prod(1:(n+m)) / prod(1:(n-m));
            legendre_v = negFactor * legendre_m(-m+1, :);
        else
            legendre_v = legendre_m(m+1, :);
        end
        k = sqrt((2*n+1) * prod(1:(n-m)) / (4*pi * prod(1:(n+m))));
        exp_v = exp(1i*m*gridPhi_v.');
        % Put spherical harmonic into matrix form, rows correspond
        % to theta, columns to phi:
        sh_m = k * legendre_v.' * exp_v;	% Y(n,m) in a matrix.
        % Add this term to 3D array of Y(n,m) values. First page is
        % n=0, then (n=1,m=-1), (n=1,m=0), (n=1,m=1), etc:
        sh_3d = cat(3, sh_3d, sh_m);
    end
end

%% Draw FODs.

% Find P(theta, phi):
gridP_m = zeros(nTheta, nPhi);	% Prob. map for 'fiber1'.
for term = 1:nTerms
    gridP_m = gridP_m + ...
        real(plm_v(term) * squeeze(sh_3d(:,:,term)));
end

% Get mask for +ve values:
mask_m = (gridP_m > 0);
gridP_m = mask_m .* gridP_m;

% Normalize to a maximum value of 1:
% maxP = max(gridP_m(:));
% if (maxP > 0)
%     gridP_m = gridP_m / maxP;
% end

% Convert from spherical (r, theta, phi) to cartesian coordinates:
x_m = gridP_m .* sin(gridTheta_m) .* cos(gridPhi_m);
y_m = gridP_m .* sin(gridTheta_m) .* sin(gridPhi_m);
z_m = gridP_m .* cos(gridTheta_m);

% Define color array:
r_m = sqrt(x_m.^2 + y_m.^2 + z_m.^2);
red_m = abs(x_m) ./ r_m;
green_m = abs(y_m) ./ r_m;
blue_m = abs(z_m) ./ r_m;
color_3d = cat(3, red_m, green_m, blue_m);

% Plot it:
surf_h = surf(x_m, y_m, z_m, color_3d, 'EdgeColor', 'none');
set(surf_h, 'FaceLighting', 'phong', 'FaceColor', 'interp', ...
    'AmbientStrength', 0.5, 'FaceAlpha', alpha)
light_h = light('Position', [1, -1, 1], 'Style', 'infinite');
%surf(x_m, y_m, z_m)
%shading flat
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')

% Define view and axes:
%az = 0;
%el = 90;
view(az, el)
