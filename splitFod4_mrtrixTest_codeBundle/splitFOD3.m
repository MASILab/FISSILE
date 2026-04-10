function [coeff_m, fodMax_m, ACC, oriFodMax_m, nOriFib] = ...
    splitFOD3(coeff_v, realFlag, volThreshold, angThreshold, plotFlag, flipFlag)
% split FOD into unique fibers
%
% input:    coeff_v -- SH coefficients of the original FOD
%           realFlag -- 1 if coeff_v in real values, 0 for complex values
%           volThreshold -- threshold for the peak magnitude. if the
%                           magnitude ratio of a small peak and the max
%                           peak is less than volThreshold, the small peak 
%                           is ignored/combined into the large peak.
%           angThreshold -- if the crossing angle btwn 2 peaks is smaller
%                           than angThreshold, these 2 peaks are combined.
%           plotFlag -- 1 to plot FODs and gridP maps, for test purpose
%           flipFlag -- if 1, the FOD peak will be flipped to the upper half
%                       spherical (fodMax_v(2)<=pi/2)
%           
% output:   coeff_m -- each col is the SH coefficients of a component. rows
%                      are listed in the descending order of the magnitude.
%           fodMax_m -- each row is the FOD maximum of a component.
%                       [magnitude, theta, phi]. rows are listed in the
%                       descending order of the magnitude.
%           ACC -- anguar correlation coefficient btwn the original and
%                  reconstructed FOD
%           oriFodMax_m -- FOD maxima matrix of the original FOD. each row
%                          is a maximum: [magnitude, theta, phi]. rows are
%                          listed in the descending order of the magnitude.
%           nOriFib -- # of peaks of the original FOD
%
% ver2: add further combining if necessary (previous combining may change
%       the magnitude of the major component). note, may result in smaller
%       fiber numbers than calling PlmToFadMaxima.
%       add fodMaxFlag to choose from 3 ways to find the FOD maximum.
% 
%
% by XH
%
% History:
% last modified 06/30/08 (XH).
% 12/27/2020: changed PlmToFodMaxima_singlePeak to version 2 (AWA).


%% prelim
fodMaxFlag = 2;
switch fodMaxFlag
    case 1
        fodMax_s = ['fodMax_v = PlmToFadMaxima(cf_v, realFlag, 1);'];
    case 2  % fastest, the speed of nDir=1002(default) is about 1/7 of method 1, angular resolution about 7 degree
        nDir = 1002;    % 4049 takes about 4 time of as long as 1002, angular resolution about 3.5 degree
        fodMax_s = ['fodMax_v = PlmToFodMaxima_singlePeak2(cf_v, realFlag, nDir, [], [], flipFlag);'];
    case 3  % about the same speed of method 1. [80, 160] gives angluar resolution of about 3 degree
        fodMax_s = ['fodMax_v = PlmToFodMaxima_singlePeak(cf_v, realFlag, [], 80, 160, flipFlag);'];
end
        

if ~exist('volThreshold', 'var') || isempty(volThreshold)
    volThreshold = 1/5;
end
if ~exist('angThreshold', 'var') || isempty(angThreshold)
    angThreshold = 35*pi/180;
end
if ~exist('plotFlag', 'var') || isempty(plotFlag)
    plotFlag = 0;
end
if ~exist('flipFlag', 'var') || isempty(flipFlag)
    flipFlag = 0;
end

nTheta = 80;
nPhi = 160;
% working on a grid of half theta range, double phi range
nHalfRows = ceil((nTheta)/2);
nHalfCols = ceil((nPhi + 1)/2);
 
if plotFlag
    fodFig = figure;
    mapFig = figure;
    nSubRowFod = 3;
    nSubColFod = 4;
    nSubRowMap = 5;
    nSubColMap = 2;
    normalSizeFlag = 0;

    % plot the original FOD (commented out; AWA)
    plotSphericalFun(coeff_v, realFlag, [],[],[],fodFig, nSubRowFod, ...
         nSubColFod,1,[], [], [], [], [], [], [], [], normalSizeFlag);
    [az,el] = view;
end

%% split the probability map
% calculate the probability, make negative values zero. 
% (AWA) The array gridP_m holds the FOD values for theta on [0,pi] and
% phi on [0, 2*pi]:
nNegFlag = 1;
gridP_m = shCoeffToGridFun(coeff_v, realFlag, nNegFlag, [], [], nTheta, nPhi);

% flip the sign
modGridP_m = -gridP_m;

% make small lobes zero
maxP = min(modGridP_m(:));
fodThres = min(modGridP_m(:)) * volThreshold + max(modGridP_m(:)) * (1 - volThreshold);
%modGridP_m(modGridP_m > fodThres) = 0;
modGridP_m(modGridP_m > maxP * volThreshold) = 0;

% half of theta range, double phi range
work_m = modGridP_m(1:nHalfRows,:);
work_m = [work_m(:, nHalfCols:(nPhi-1)), work_m, work_m(:, 2:(nHalfCols-1))];

% watershed
ws_m = watershed(work_m);

% recover original phi range
modWs_m = ws_m(:, nPhi-nHalfCols+1 : 2*nPhi-nHalfCols);

% combine connected regions on both sides of the line phi=0/2pi
phi0_v = modWs_m(:,1);
zeroIdx_v = find(phi0_v==0);    % theta index where watershed ridge occurs
nZero = length(zeroIdx_v);

if nZero == 0   % no ridge occurs on the phi=0/2pi line, combine 2 sides
    regionIdx1 = modWs_m(1,1);
    regionIdx2 = modWs_m(1, end);
    modWs_m(modWs_m == regionIdx2) = regionIdx1;

elseif nZero < length(phi0_v)   % ridge cross the line, combine segment by segment
    for k = 1:nZero
        zeroIdx = zeroIdx_v(k);
        if zeroIdx > 1
            regionIdx1 = modWs_m(zeroIdx-1, 1);
            regionIdx2 = modWs_m(zeroIdx-1, end);
            modWs_m(modWs_m == regionIdx2) = regionIdx1;
        end
    end
    % combine the last segment
    regionIdx1 = modWs_m(min(zeroIdx+1, nHalfRows),1);
    regionIdx2 = modWs_m(min(zeroIdx+1, nHalfRows), end);
    modWs_m(modWs_m == regionIdx2) = regionIdx1;
end

% plot the maps
if plotFlag
    figure(mapFig);
    for k = 1:8
        subplot(nSubRowMap, nSubColMap, k)
        switch k
            case 1
                map_m = gridP_m;
                title_s = 'original Probability map P';
            case 2
                map_m = modGridP_m;
                title_s = '-(+P), w/o small lobes';
            case 3
                map_m = work_m;
                title_s = 'half theta range, double phi range';
            case 4
                wsP_m = (work_m+1).*(ws_m~=0);
                map_m = wsP_m;
                title_s = 'P w/ watershed';
            case 5
                map_m = label2rgb(ws_m);
                title_s = 'divided regions';
            case 6
                map_m = wsP_m(:, nPhi-nHalfCols+1 : 2*nPhi-nHalfCols);
                title_s = 'P w/ watershed, original phi range';
            case 7
                map_m = label2rgb(ws_m(:, nPhi-nHalfCols+1 : 2*nPhi-nHalfCols));
                title_s = 'divided regions, original phi range';
            case 8
                map_m = label2rgb(modWs_m);
                title_s = 'regions combined by phi';
        end
        if ismember(k, [5,7])
            imshow(map_m);
        else
            imagesc(map_m); axis equal; colorbar
        end
        title(title_s)
        % show half theta range
        if k == 2
            hold on
            line([0 nPhi], [nHalfRows nHalfRows], 'color', 'k')
        end
        % show the original phi range
        if ismember(k, [3,4])
            minRowIdx = -nHalfRows/10;
            maxRowIdx = nHalfRows * 1.1;
            hold on
            line([nPhi-nHalfCols+1 nPhi-nHalfCols+1], [minRowIdx, maxRowIdx], 'color', 'k')
            line([2*nPhi-nHalfCols 2*nPhi-nHalfCols], [minRowIdx, maxRowIdx], 'color', 'k')
        end
    end
end



%% split FOD
label_v = unique(modWs_m(:));
count = length(label_v)-1;   % excluding the ridge value 0
cmbWs_m = modWs_m;

if count == 0 || count == 1   % no need to split
    [oriFodMax_m, nOriFib] = fod2Maxima(gridP_m, 3, nTheta, nPhi, volThreshold);
    coeff_m = coeff_v;
    fodMax_m = oriFodMax_m(1,:);
    ACC = 1;

    % for plotting
    if plotFlag
        nFib = 0;
        sumNewCoef_v = coeff_v;
    end
    
else
    order = sqrt(length(coeff_v))-1;
    evenOrderFitFlag = 1;
    [gridTheta_v, gridPhi_v, gridPhi_m, gridTheta_m] = gridThetaPhi(nTheta, nPhi);
    hlfGridPhi_m = gridPhi_m(1:nHalfRows,:);        % (AWA) theta = 0 through pi/2.
    hlfGridTheta_m = gridTheta_m(1:nHalfRows,:);    % (AWA) theta = 0 through pi/2.
    
    existLabel_v = [];
    coeff_m = [];
    fodMax_m = [];
    nFib = 0;
    for fibIdx = 1:count
        currentLabel = label_v(fibIdx+1);
        mask_m = (modWs_m==currentLabel);
        p_m = gridP_m(1:nHalfRows, :) .* mask_m;

        newCoeff_v = harmonicTransform2_xh2(p_m(:), hlfGridTheta_m(:), hlfGridPhi_m(:), order, realFlag, evenOrderFitFlag);
        if any(newCoeff_v)  % some region may contain no peak
            % find maximum. AWA: eval calculates 
            %       fodMax_v = [magnitude, theta, phi] 
            % of peak:
            cf_v = newCoeff_v;
            eval(fodMax_s)
            
            if nFib == 0
                combineFlag = 0;
            else
                % calculate crossing angle w.r.t. existing fibers
                crsAng_v = angleBetw2vect(fodMax_v(2), fodMax_v(3), fodMax_m(1:nFib,2), fodMax_m(1:nFib,3), 1);
                
                % if the smallest crsAnle is smaller than the threshold, combine
                [minAng, minIdx] = min(crsAng_v);
                if minAng < angThreshold
                    cmbCoeff_v = newCoeff_v + coeff_m(:,minIdx);
                    cf_v = cmbCoeff_v;
                    eval(fodMax_s)

                    coeff_m(:,minIdx) = cmbCoeff_v;
                    fodMax_m(minIdx, :) = fodMax_v;
                    combineFlag = 1;
                    cmbWs_m(cmbWs_m==currentLabel) = existLabel_v(minIdx);
                else
                    combineFlag = 0;
                end

                % if volume fraction too small, combine
                if combineFlag == 0
                    volRatio_v = fodMax_m(:,1)./fodMax_v(1);
                    idx_v = find(volRatio_v > 1);
                    volRatio_v(idx_v) = 1./volRatio_v(idx_v);
                    [minVol, minIdx] = min(volRatio_v);
                    if minVol < volThreshold
                        cmbCoeff_v = newCoeff_v + coeff_m(:,minIdx);
                        cf_v = cmbCoeff_v;
                        eval(fodMax_s)

                        coeff_m(:,minIdx) = cmbCoeff_v;
                        fodMax_m(minIdx, :) = fodMax_v;
                        combineFlag = 1;
                        cmbWs_m(cmbWs_m==currentLabel) = existLabel_v(minIdx);
                    else
                        combineFlag = 0;
                    end
                end
            end

            if combineFlag == 0
                nFib = nFib+1;
                coeff_m(:,nFib) = newCoeff_v;
                fodMax_m(nFib, :) = fodMax_v;
                existLabel_v = [existLabel_v, currentLabel];
            end
        end
    end

    % sort coeff_m, fodMax_m in descending order of magnitude
    [srtVol_v, srtIdx_v] = sort(fodMax_m(:,1), 'descend');
    fodMax_m = fodMax_m(srtIdx_v, :);
    coeff_m = coeff_m(:, srtIdx_v);
    existLabel_v = existLabel_v(srtIdx_v);
    
    % further combining
    if size(fodMax_m, 1) > 1
        combineFlag = 1;
        while combineFlag
            volRatio_v = fodMax_m(end,1)./fodMax_m(1:end-1,1);
            idx_v = find(volRatio_v < volThreshold);
            if ~isempty(idx_v)
                % combine
                currentLabel = existLabel_v(end);
                if length(idx_v) > 1
                    % if there are 1+ major components, combine with the one with smaller angle
                    ang_v = angleBetw2vect(fodMax_m(end,2), fodMax_m(end,3), fodMax_m(idx_v,2), fodMax_m(idx_v,3),1);
                    [minAng, minIdx] = min(ang_v);
                    cmbIdx = idx_v(minIdx);
                else
                    cmbIdx = idx_v;
                end

                cmbCoeff_v = coeff_m(:,end) + coeff_m(:,cmbIdx);
                cf_v = cmbCoeff_v;
                eval(fodMax_s)

                coeff_m(:,cmbIdx) = cmbCoeff_v;
                coeff_m = coeff_m(:, 1:end-1);
                fodMax_m(cmbIdx, :) = fodMax_v;
                fodMax_m = fodMax_m(1:end-1, :);
                nFib = nFib - 1;
                existLabel_v = existLabel_v(1:end-1);
                cmbWs_m(cmbWs_m==currentLabel) = existLabel_v(cmbIdx);
                combineFlag = 1;
            else
                combineFlag = 0;
            end
        end
    end
    
    sumNewCoef_v = sum(coeff_m, 2);

    % calcularte ACC btwn original coeff_v and the reconstructed sum coeff_v
    if nargout > 2 || plotFlag == 1
        ACC = angularCorrCoeff_xh(coeff_v, sumNewCoef_v);
    end

    if nargout > 3 || plotFlag == 1
        [oriFodMax_m, nOriFib] = fod2Maxima(gridP_m, 3, nTheta, nPhi, volThreshold);
    end
end



if plotFlag
    % plot FOD
    % set axes limit for the original FOD
    figure(fodFig);
    subplot(nSubRowFod, nSubColFod, 1);
    if nOriFib >= 2
        volRat = oriFodMax_m(2,1)/oriFodMax_m(1,1);
    else
        volRat = 0;
    end
    title_s = {'original FOD', ['nFib=', num2str(nOriFib)], ['volRat=', num2str(volRat)]};
    title(title_s)
    xlim([maxP -maxP]);ylim([maxP -maxP]);zlim([maxP -maxP])

    % plot each peak
    for k = 1:nFib
        newCoeff_v = coeff_m(:,k);
        plotSphericalFun(newCoeff_v, realFlag, [],[],['fiber',num2str(k)],fodFig,...
            nSubRowFod, nSubColFod, nSubColFod+k, [], [], [], [], [], [], [], [], normalSizeFlag);
        view(az,el);
        xlim([maxP -maxP]);ylim([maxP -maxP]);zlim([maxP -maxP])
    end

    % sum FOD
    plotSphericalFun(sumNewCoef_v, realFlag, [],[],[],fodFig,nSubRowFod, ...
        nSubColFod,2,[], [], [], [], [], [], [], [], normalSizeFlag);
    view(az,el);
    xlim([maxP -maxP]);ylim([maxP -maxP]);zlim([maxP -maxP])
    if nFib > 1
        volRat = srtVol_v(2)/srtVol_v(1);
    else
        volRat = 0;
    end
    title_s = {'reconstructed FOD', ['ACC=',num2str(ACC)], ['volRat=', num2str(volRat)]};
    title(title_s)
    
    % plot combined map
    figure(mapFig)
    subplot(nSubRowMap, nSubColMap, 9)
    imshow(label2rgb(cmbWs_m));
    title('regions combined by crsAng & vol')
    
%     subplot(nSubRowMap, nSubColMap, 10)
%     wsP2_m = (work_m(:, nPhi-nHalfCols+1 : 2*nPhi-nHalfCols)+1) .* (cmbWs_m~=0);
%     imagesc(wsP2_m); axis image; colorbar
%     title('final regions')
    
end

