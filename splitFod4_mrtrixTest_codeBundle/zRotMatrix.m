function rot_m = zRotMatrix(angle)
% Usage: rot_m = zRotMatrix(angle)
%
% zRotMatrix calculates the rotation matrix for the given angle (in
% radians):
%
%					cos		-sin	0
%			Rot =	sin		cos		0
%					0		0		1

% Calulate matrix:
rot_m = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];