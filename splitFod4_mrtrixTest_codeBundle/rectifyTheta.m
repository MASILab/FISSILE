function [opTheta_v, opPhi_v] = rectifyTheta(ipTheta_v, ipPhi_v)
% Usage: [opTheta_v, opPhi_v] = rectifyTheta(ipTheta_v, ipPhi_v)
%
% rectifyTheta tests the sign of each element of ipTheta_v. For any 
% (theta,phi) pair with theta<0, flip the sign of theta and
% add (or subtract) pi to/from phi. This maintains the location of the
% physical point for the purposes of a search on the unit sphere, while 
% keeping theta>0 (required by sphericalHarm.m).

% If theta<0, flip the sign of theta and subtract pi from phi if phi>=0 or 
% add pi to phi if phi<0 (to keep phi on [-pi,pi]). Start by initializing
% output arrays:
opTheta_v = ipTheta_v;
opPhi_v = ipPhi_v;
% Overwrite values where ipTheta < 0:
negThetaIndex_v = find(ipTheta_v < 0);        
opTheta_v(negThetaIndex_v) = -ipTheta_v(negThetaIndex_v);
opPhi_v(negThetaIndex_v) = ipPhi_v(negThetaIndex_v) - ...
    2*((ipPhi_v(negThetaIndex_v) >= 0) - 0.5) * pi;
end