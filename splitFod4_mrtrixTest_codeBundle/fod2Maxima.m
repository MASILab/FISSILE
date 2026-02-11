function [fadMax_m, nPoles] = fod2Maxima(gridP_m, nProcFibers, nTheta, nPhi, volFraRatioThreshold)
% usage: [fadMax_m, nPoles] = fod2Maxima(gridP_m, nProcFibers, nTheta, nPhi, volFraRatioThreshold)
%        (nProcFibers... input optional)
% calculate FAD from probability, output the maxima, and the
% number of local maxima.
% note: pancake is presented by 2  maxima with crossing angle of pi/2 
% note: return only the upper lobe of the FAD maxima (theta <= pi/2+)
% except for pancake along z axis where 2 maxima with same phi and theta
% pi/2 away are returned
% if volFraRatioThreshold is not given or empty, return unmodified nPoles;
% otherwise, ignore small lobes (if volFraRatioThreshold == 0, use default
% value (1/5))
%
% call:
%   maxInSpheCoor2
%   ReadMjCmpn
%
% by XH
% last modified 01/01/09


if ~exist('nProcFibers', 'var') || isempty(nProcFibers)
    nProcFibers = 3;
end

if ~exist('nTheta', 'var') || isempty(nTheta)
    nTheta = 80;
    nPhi = 160;
end
gridTheta_v = pi * (0:(nTheta-1)).' / (nTheta-1);
gridPhi_v = 2*pi * (0:(nPhi-1)).' / (nPhi-1);


% get probability maxima
mask_m = (gridP_m > 0);
if mask_m == zeros(size(gridP_m))
    fadMax_m = zeros(nProcFibers,3);
    nPoles = 0;
else
    work_m = gridP_m .* mask_m;
    % find the most probable direction:
    fadMax_m = maxInSpheCoor2(work_m, nProcFibers);
    nPoles = ReadMjCmpn(fadMax_m, 0);
    if nPoles == 0   % if it's a sphere, set all components to zeros
        if nProcFibers > 1
            fadMax_m(2:nProcFibers, :) = zeros(nProcFibers-1, 3);
        end
    else    % convert the indexes of theta, phi for each components
            % into their values
        if mod(nPoles, 1) > 0   % pancake exist
           if fadMax_m(1,2) == -10    
               % pancake along z, pancake gives the largest peak
               % represent the pancake by 2 peaks crossed at pi/2
               fadMax_m(1,2) = 10;
               fadMax_m(3,:) = fadMax_m(2,:);
               fadMax_m(2,:) = fadMax_m(1,:);
               fadMax_m(2,2) = fadMax_m(1,2) + nTheta/2;
           elseif fadMax_m(2,2) == -10    % pancake along z, largest peak in xy plane
               fadMax_m(2,2) = 10;
               fadMax_m(3,:) = fadMax_m(2,:);
               fadMax_m(3,2) = fadMax_m(2,2) + nTheta/2;
           elseif fadMax_m(1,3) == -10
               % pancake in xy plan, pancake gives larges peak
               fadMax_m(1,3) = 10;
               fadMax_m(3,:) = fadMax_m(2,:);
               fadMax_m(2,:) = fadMax_m(1,:);
               fadMax_m(2,3) = fadMax_m(1,3) + nPhi/2;
           elseif fadMax_m(2,3) == -10
               fadMax_m(2,3) = 10;
               fadMax_m(3,:) = fadMax_m(2,:);
               fadMax_m(3,3) = fadMax_m(2,3) + nPhi/2;
           end
        end
        
        fadMax_m(:,2) = gridTheta_v(fadMax_m(:,2));
        fadMax_m(:,3) = gridPhi_v(fadMax_m(:,3));
        
    end
end

% ignore the small lobes, get actual # of lobes
if exist('volFraRatioThreshold', 'var') && ~isempty(volFraRatioThreshold)
    
    if volFraRatioThreshold == 0
        volFraRatioThreshold = 1/5;
    end
    
    if nPoles > 1
        if fadMax_m(2,1)/fadMax_m(1,1) >= volFraRatioThreshold  % fadMax_m is sorted according to the magnitude of each local maxima
            if (size(fadMax_m,1)>2) && (fadMax_m(3,1)/fadMax_m(1,1) < volFraRatioThreshold)
                nPoles = 2;
            end
        else
            nPoles = 1;
        end
    end
end