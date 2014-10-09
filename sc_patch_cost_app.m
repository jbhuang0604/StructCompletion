function [costApp, uvBias] = sc_patch_cost_app(trgPatch, srcPatch, wDistPatch, optS)

% SC_PATCH_COST_APP
%
% Compute the weighted sum of the absolute difference between
% cost between source and target patches

% Input:
%   - trgPatch, srcPatch, wDistPatch, optS
% Output:
%   - costApp
%   - uvBias

% Initialization

numUvValidPix = size(wDistPatch, 2);

% Apply bias correction
if(optS.useBiasCorrection)
    % Mean of source and target patch
    meanTrgPatch = mean(trgPatch, 1);
    meanSrcPatch = mean(srcPatch, 1);
    
    % Compute bias and clamp it to inteval [optS.minBias, optS.maxBias]
    biasPatch = meanTrgPatch - meanSrcPatch;
    biasPatch = sc_clamp(biasPatch, optS.minBias, optS.maxBias);
    
    % Update the UV map for gain and bias
    uvBias = reshape(biasPatch, 3, numUvValidPix);
    
    % Compute patch appearance cost
    srcPatch = bsxfun(@plus, srcPatch, biasPatch);
end
patchDist = trgPatch - srcPatch;

% Weight patch appearance cost by their distances to hole boarder
wDistPatchC = reshape(wDistPatch, optS.pNumPix, 1, numUvValidPix);

% Sum of absolute distance
if(strcmp(optS.costType, 'L1'))
    patchDist = abs(patchDist);
elseif(strcmp(optS.costType, 'L2'))
    patchDist = patchDist.^2;
end

% Apply weights
patchDist = bsxfun(@times, patchDist, wDistPatchC);

patchDist = sum(sum(patchDist, 1),2);
patchDist = reshape(patchDist, 1, numUvValidPix);

% Weight normalization
sumDistWeight = sum(wDistPatch, 1);
costApp = patchDist./sumDistWeight;

end