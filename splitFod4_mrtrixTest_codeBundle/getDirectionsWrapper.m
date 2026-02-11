function n_m = getDirectionsWrapper(n)
% Usage: n_m = getDirectionsWrapper(n)
%
% getDirectionsWrapper gets a list of directions from any of three sources:
%
%   1. Philips standard sets (low: 6, medium: 15, high: 32, or our '92' set)
%   2. Electrostatic repulsion of n charges in 3D
%   3. Any set written by Mathematica, har4049.nb, for example.
%
% Note that 92 can be taken from any of the 3 sets, but here we use the 
% 92 directions calculated with Mathematica (tesselation of a sphere) and
% included in the Philips sets (added by hand). 
%
% Input argument:
%   n: integer number of diffusion weighting directions.
%   
% Output variable:
%   n_m: (n x 3) matrix of x, y, and z components of unit vectors
%       (approximately) uniformly distributed on a sphere. The components
%       of the nth unit vector are on the nth row of n_m in the order x, y,
%       z.
%
% getDirectionsWrapper calls the functions:
%   n_m = getPhilipsDirs(res_s, overplusFlag)
%   n_m = readElectroStaticDirs(n)
%   n_m = getDirsFromFile(dirsFile_s, tabFlag)

% Find the directions for given n:
switch n
    case {6, 15, 32, 46, 92}
        overplusFlag = 0;
        % Define resolution string:
        switch n
            case 6
                res_s = 'low';
            case 15
                res_s = 'medium';
            case 32
                res_s = 'high';
            case {46, 92}
                res_s = '92';
        end
        n_m = getPhilipsDirs(res_s, overplusFlag);
        % Take the first n directions (for the Philips direction sets, the
        % (n+1)-th row corresponds to the synthetic isotropic weighting. 
        % For 46 directions, take the 1st half of the 92 direction set (the
        % second half are just the opposite of the first half):
        n_m = n_m(1:n, :);  % Take the first n vectors.
    case {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 128, ...
            130, 132, 192, 212, 272, 282}
        n_m = readElectroStaticDirs(n);
    case 23
        % Read harDirs23set2.txt text file:
        tabFlag = 0;
        dirsFile_s = ['/Users/andera/Docs/', ...
            'MATLAB/harDti/harDirs23set2.txt'];
        n_m = getDirsFromFile(dirsFile_s, tabFlag);  
    case {1, 2, 3, 4, 5, 7, 8, 9}
        % Take the first n directions from readElectroStaticDirs(10):
        n10_m = readElectroStaticDirs(10);
        n_m = n10_m(1:n, :);  % Take the first n vectors.
    case {4002}
        % Read text file created by har4049.nb, for example:
        tabFlag = 0;
        dirsFile_s = ['/Users/andera/Docs/', ...
            'Mathematica/harDirs', num2str(n), '.txt'];
        n_m = getDirsFromFile(dirsFile_s, tabFlag);        
    otherwise
        % Read text file created by har4049.nb, for example:
        tabFlag = 1;
        dirsFile_s = ['/Users/andera/Docs/', ...
            'Mathematica/harDirs', num2str(n), '.txt'];
        n_m = getDirsFromFile(dirsFile_s, tabFlag);
end
