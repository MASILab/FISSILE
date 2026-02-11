function rot_m = yRotMatrix(angle)
% Usage: rot_m = yRotMatrix(angle)
%
% yRotMatrix calculates the rotation matrix for the given angle (in
% radians):
%
%					cos		0		sin
%			Rot =	0		1		0
%					-sin	0		cos

% Calulate matrix:
rot_m = [cos(angle) 0 sin(angle); 0 1 0; -sin(angle)  0 cos(angle)];