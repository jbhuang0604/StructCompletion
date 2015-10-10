function [wPatchR, wPatchSumImg] = sc_prep_dist_patch(distMap, trgPixPos, optS)

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
wPatchR  = sc_prep_target_patch(distMap, trgPixPos, optS);
wPatchR  = squeeze(wPatchR);

wPatchR  = bsxfun(@minus, wPatchR, wPatchR(optS.pMidPix,:));
wPatchR  = optS.wDist(optS.iLvl).^ (wPatchR);

% Sum of the patch weights for all unknown pixels
numUvPix = size(wPatchR, 2);

wPatchSumImg = zeros(imgH, imgW, 'single');
indMap    = reshape(1:imgH*imgW, imgH, imgW);
indPatch  = sc_prep_target_patch(indMap, trgPixPos, optS);
indPatch  = squeeze(indPatch);
for i = 1: numUvPix
    wPatchSumImg(indPatch(:,i)) = wPatchSumImg(indPatch(:,i)) + wPatchR(:,i);
end
 
end