function maxima_m = maxInSpheCoor2(A_m, numMax)

% based on maxInSpheCoor
% difference: 1. only one output (kick out nPoles)
%             2. when nPoles=0, maxima_m(1,2:3)=[0,0] instead of [1,1]
%                when nPoles=1, maxima_m(1,2:3)~=[0,0]

% maxInSpheCoor find unique maxima in a spherical coordinates defined
% by matrix A_m. numMax is the required number of maxima. nPoles is the
% actual number of unique maxima.
% the first dim of A_m is theta varying [0,pi/2], the 2nd is phi varying [0,2*pi].
% If all the elements of A_m are equal, then nPoles = 0.
% nPoles is the actual number of unique maxima. The output maxima_m is a
% [numMax,3] matrix, each row gives the value and index of the maxima.
% All the rows are sorted by the value in descending order.
% Initial row values are [0 1 1] for those unoccupied rows.
% sphere: first row of maxima_m is set to [radius, -10, -10]. -10 is chosen
% b/c it is not valid index nor in the valid theta, phi range;
% pancake: row of maxima_m is set to [radius,-10,x] where -10 means
% theta could be [0 pi] while phi is same; or [radius,x,-10] where the
% pancake is in x-y plane (phi could be [0 2*pi]
%
% note: sometimes all the maxima appear in the lower half-sphere
%           (theta>pi/2)
%       when theta near pi/2, theta resolution is pi/2 +-
%           thetaGridThreahold; phi resolution is +-1.
%       
%       pancake part not finished yet 
%


% if no 2nd input, look for up to 10 unique maxima.
if nargin < 2
    numMax = 5;
end

% initialize matrix. default max value is zero, default row, col is 1
maxima_m = ones(numMax, 3);
maxima_m(:, 1) = 0;

[nRows, nCols] = size(A_m);
% [nHalfRows, nCols-1] is the size of the actual working matrix
% nHalfRows = ceil(nRows/2) + 1;
thetaGridThreahold = ceil(nRows/40);    % offset from pi/2 for theta
nHalfRows = ceil(nRows/2) + thetaGridThreahold; % sometimes maxima appear lower farther

% get a full-length,double-width matrix:
% the first part ends with col (nCols-1), and the 3rd part begin from col 2,
% b/c the 1st and last col of A_m are the same(when phi equals 0 and 2*pi)
B_m = A_m;
nHalfCols = ceil((nCols + 1)/2);
B_m = [B_m(:, nHalfCols:(nCols-1)), B_m, B_m(:, 2:(nHalfCols-1))];

% get the mask of local maxima
lclMaxMsk_m = imregionalmax(B_m);

% if all the elements of A_m are equal (the mask is a unit matrix):
if lclMaxMsk_m == ones(size(B_m))
    maxima_m(1,:) = [A_m(1,1), -10, -10];

else
    % consider only the [0:pi/2+, 0:2*pi) part (size:[nHalfRows, nCols-1])
    lclMaxMsk_m = lclMaxMsk_m(1:nHalfRows, ...
            (nCols - nHalfCols + 1):(2 * nCols - nHalfCols - 1));
        
    % if maxima occurs at theta==0, then maxium occur for all the phi.
    % keep only the first one.
    if lclMaxMsk_m(1,:) == ones(1,nCols-1)
        lclMaxMsk_m(1,2:end) = zeros;
    end
    
    [nzRow_v, nzCol_v] = find(lclMaxMsk_m > 0);
    nzRow_v = unique(nzRow_v);
    nzCol_v = unique(nzCol_v);
    
    if (length(nzRow_v) == size(lclMaxMsk_m,1)) % one pancake along z
        if (length(nzCol_v) == 1)   % one pancake only
            phiIdx = find(lclMaxMsk_m(1,:), 1); % indix of phi
            maxima_m(1,:) = [A_m(1,1), -10, phiIdx];
            
%         elseif ??? % plus one peak in xy plane
            
        else  % plus multiple peaks
            % ...
        end
        
    elseif length(nzCol_v) == size(lclMaxMsk_m,2) % one pancake in xy plane
        if (length(nzRow_v) == 1) % one pancake only
            thetaIdx = find(lclMaxMsk_m(:, 1), 1); % index of theta
            maxima_m(1,:) = [A_m(nHalfRows,1), thetaIdx, -10];
            
        elseif (lclMaxMsk_m(1,1) == 1)   % plus one peak alone z
            magRod = A_m(1,1);
            thetaIdx = max(nzRow_v(nzRow_v ~= 1));  % index of theta
            magPancake = A_m(nHalfRows,1);
            if magRod >= magPancake
                maxima_m(1:2,:) = [magRod 1 1; magPancake thetaIdx -10];
            else
                maxima_m(1:2,:) = [magPancake thetaIdx -10; magRod 1 1];
            end
            
        else    % plus other peaks
            % ...
        end

    else % no pancake
        % get rid of neighboring same value points
        lclMaxMsk_m = bwmorph(lclMaxMsk_m, 'shrink');

        % get the index of the local max points:
        lclMaxInd_v = find(lclMaxMsk_m == 1);
        % number of the local maxima
        nLclMax = length(lclMaxInd_v);

        % loop till there's no more local maxium, or already got enough maxima,
        % each time find the maxium, eliminates it and all others maxima
        % in the same axis:
        nPoles  = 0;
        while (nLclMax > 0) & (nPoles < numMax)
            % find the largest one among the local maxima as the current point
            lclMax_m = lclMaxMsk_m .* A_m(1:nHalfRows, 1:(nCols-1));
            maxVal = max(lclMax_m(:));
            [row,col] = find(lclMax_m == maxVal);
            row = row(1);
            col = col(1);
            nPoles = nPoles + 1;

            % initialize mask
            crspPnMsk_m = zeros(nHalfRows, nCols-1);
            crspPnMsk_m(row,col) = 1;

            % if there are more than one local maxi, check if some of them in
            % the same axis
            if nLclMax > 1

                % if phi = 0, 2*pi (col = 1,2, nCols-1),
                % pick the 2pi corresponding point
                if (col < 3) | (col == nCols-1)
                    testRow1 = max(1, (row-1));
                    testRow2 = min(nHalfRows, (row+1));
                    crspPnMsk_m(testRow1:testRow2, 1:2) = 1;
                    crspPnMsk_m(testRow1:testRow2, (nCols-1)) = 1;
                end

                % if theta = pi/2 (row = (nHalfRows-thetaGridThreahold-1) : nHalfRows),
                % also pick the pi corresponding point
                if (row >= (nHalfRows-thetaGridThreahold-1))
                    testRow1 = nHalfRows - thetaGridThreahold - 1;
                    testRow2 = nHalfRows;

                    % if phi = pi, pick points with phi = 0 or 2*pi
                    if (col >= floor(nCols/2)) & (col <= ceil(nCols/2+1))
                        crspPnMsk_m(testRow1:testRow2, 1:2) = 1;
                        crspPnMsk_m(testRow1:testRow2, nCols-1) = 1;

                    else    % pick points pi away
                        if col < nHalfCols
                            testCol1 = col + floor(nCols/2-1);
%                             testCol2 = min((nCols-1), (col + ceil(nCols/2)));
                            testCol2 = min((nCols-1), (col + ceil(nCols/2)+1));
                        else
%                             testCol1 = max(1, (col - ceil(nCols/2)));
                            testCol1 = max(1, (col - ceil(nCols/2)-1));
                            testCol2 = col - floor(nCols/2-1);
                        end
                        crspPnMsk_m(testRow1:testRow2, testCol1:testCol2) = 1;
                    end     % end of col
                end     % end of row

            end     % end of nLclMax

            % eliminate current maxi point and all its corresponding points
            crspPnMsk_m = crspPnMsk_m .* lclMaxMsk_m;
            crspPntIdx_v = find(crspPnMsk_m == 1);
            % zero all these points:
            lclMaxMsk_m(crspPntIdx_v) = 0;
            % update the number of the local maxima
            lclMaxInd_v = find(lclMaxMsk_m == 1);
            nLclMax = length(lclMaxInd_v);

            maxima_m(nPoles,:) = [maxVal, row, col];

        end     % end of while nLclMax & nPoles
    end
end     % end of lclMaxMsk_m