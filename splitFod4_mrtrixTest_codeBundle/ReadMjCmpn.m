function [nCmpns, cmpn_m] = ReadMjCmpn(cmpn_m, squeezeFlag);

% return the number of actual components, 
% if squeezeFlag = 1 then retain only the non-zero parts

if nargin < 2
    squeezeFlag = 1;
end

% actual number of the major components
nCmpns = length(find(cmpn_m(:,1) ~= 0));
if nCmpns >= 1 
    if cmpn_m(1,2:3) == [-10,-10]
        nCmpns = 0;     % sphere
    elseif any(cmpn_m(:,2) == -10)  % pancake along z
        nCmpns = nCmpns + 0.1;   
    elseif any(cmpn_m(:,3) == -10)  % pancake in xy plane
        nCmpns = nCmpns + 0.2;   
    end
end

if squeezeFlag == 1
    if nCmpns < 1
        cmpn_m = cmpn_m(1,:);
    else
        cmpn_m = cmpn_m(1:nCmpns, :);
    end
end
