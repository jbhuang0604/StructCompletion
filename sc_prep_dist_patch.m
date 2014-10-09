function [wDistPatch, wDistImg] = sc_prep_dist_patch(distMap, trgPixPos, iLvl, optS)

% SC_PREP_DIST_PATCH
% 
% Precompute patch weights and summation for patch matching and voting
% 
% Input:
%   - distMap: 
%   - trgPixPos: 
%   - iLvl:
%   - optS:
% Output:
%   - distWPatch:
%   - distWImg:

[imgH, imgW] = size(distMap);

% Compute the patch weights
wDistPatch  = sc_prep_target_patch(distMap, trgPixPos, optS);
wDistPatch = bsxfun(@minus, wDistPatch, wDistPatch(optS.pMidPix,:));
wDistPatch  = optS.wDist(iLvl).^ (- wDistPatch); 

% Sum of the patch weights for all unknown pixels
numUvPix = size(wDistPatch, 2);

wDistImg = zeros(imgH, imgW, 'single');
indMap = reshape(1:imgH*imgW, imgH, imgW);
indPatch  = sc_prep_target_patch(indMap, trgPixPos, optS);

for i = 1: numUvPix
    wDistImg(indPatch(:,i)) = wDistImg(indPatch(:,i)) + wDistPatch(optS.pMidPix,i);
end
 
end