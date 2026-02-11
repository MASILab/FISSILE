function vuPointsOnSphere3(theta_v, phi_v, color_s, title_s)
% Usage: vuPointsOnSphere3(theta_v, phi_v, color_s, title_s)
%
% vuPointsOnSphereX displays points on a spherical surface. The positions of
% the points are given by corresponding elements of theta_v (the polar
% angle in a 1D array) and phi_v (the azimuthal angle in a 1D array). The
% color of the point is given by the single character in color_s (e.g., 'b').
% The plot is given the title in <title_s>.
%
% The sphere is drawn in the current figure (which allows points to be
% added in different colors by multiple calls to vuPointsOnSphere2). 

% Define unit sphere for plotting points on:
nPoints = 30;
[xSphere_m, ySphere_m, zSphere_m] = sphere(nPoints);
sphereColor_m = ones(size(xSphere_m));
faceAlpha = 1.0;  % Transparency (there's a bug in exporting graphics for values <1).
level = 0.4;

% Display sphere:
%figure
surf(xSphere_m, ySphere_m, zSphere_m, sphereColor_m, ...
    'FaceAlpha', faceAlpha)
colormap([level, level, level])
axis image
%axis off
hold on

% Find coordinates on sphere:
xPoints_v = sin(theta_v).*cos(phi_v);
yPoints_v = sin(theta_v).*sin(phi_v);
zPoints_v = cos(theta_v);

% Display rectified v1 points:
liftOff = 1.01; % Bring points slightly above surface.
markerSize = 12;    % In points/3 (i.e., 1/(3*72)").
plot3(liftOff*xPoints_v, liftOff*yPoints_v, ...
    liftOff*zPoints_v, ['.', color_s], 'MarkerSize', markerSize)

% View sphere and points:
view([1,1,1])
xlabel('X')
ylabel('Y')
zlabel('Z')
% Add title:
title(title_s)
hold off

