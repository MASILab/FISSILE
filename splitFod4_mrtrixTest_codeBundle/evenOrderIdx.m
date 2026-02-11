function idx_v = evenOrderIdx(order, oddFlag)
% function evenOrderIdx gives the indices of even-order items for spherical
% harmonics

if ~exist('oddFlag')
    oddFlag = 0;    % default: even order
end

idx_v = [];
if oddFlag == 1
    kk_v = 1:2:order;
else
    kk_v = 0:2:order;
end

for k = kk_v
    idx_v = [idx_v, (k^2+1):(k+1)^2];
end
