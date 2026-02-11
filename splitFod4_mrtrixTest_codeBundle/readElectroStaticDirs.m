function n_m = readElectroStaticDirs(n)
% Usage: n_m = readElectroStaticDirs(n)
%
% Read directions from a text file downloaded from the website "Minimal
% Energy Arrangements of Points on a Sphere", by RH Hardin, NJA Sloane, 
% and WD Smith, originally at:
%   http://www2.research.att.com/~njas/electrons/
% but now at
%   http://neilsloane.com/electrons/
% For 3 dimensions, the list goes up to 282 point charges. On the local
% filesystem, the files should be named harDirs_elec_3_<nDirs>.txt. The 
% files list the (x,y,z) coordinates of the point charges, one coordinate
% per line (x1, y1, z1, x2, y2, z2, etc. on separate lines). 

% Set the path:
if ispc
    % Assume the Dell laptop:
    path_s = 'C:\Users\AWA\Documents\matlab\harDti';
else
    % Assume the mbp:
    path_s = ['/Users/andera/Documents/Documents - VUSE73DCFNNF/', ...
            'MATLAB/harDti'];
end
% Get filename:
dirsFile_s = [path_s, filesep, 'harDirs_elec_3_', num2str(n), '.txt'];
dirsFid = fopen(dirsFile_s, 'r', 'native', 'UTF-8');

% Read components:
comps_v = [];
while 1
	line_s = fgetl(dirsFid);
	if ~isstr(line_s)
		break
    end
    if isempty(comps_v)
        comps_v = str2num(line_s(2:end));   % Skip 1st character (UTF-8 code 65279)
    else
        comps_v = [comps_v; str2num(line_s)];
    end
end
fclose(dirsFid);

% Reformat to dirs_m:
n_m = reshape(comps_v, 3, n).';

