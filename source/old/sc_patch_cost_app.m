function [costApp, uvBias] = sc_patch_cost_app(trgPatch, srcPatch, optS)

% SC_PATCH_COST_APP
%
% Compute the weighted sum of the absolute difference between
% cost between source and target patches

% Input:
%   - trgPatch:  target patch
%   - srcPatch:  source patch
%   - optS:      parameter
% Output:
%   - costApp:   appearance cost
%   - uvBias     bias for the closet neighbor

% =========================================================================
% Apply bias correction
% =========================================================================

uvBias = [];
if(optS.useBiasCorrection)
    % Mean of source and target patch
    meanTrgPatch = mean(trgPatch, 1);
    meanSrcPatch = mean(srcPatch, 1);
    
    % Compute bias and clamp it to inteval [optS.minBias, optS.maxBias]
    uvBias = meanTrgPatch - meanSrcPatch;
    uvBias = sc_clamp(uvBias, optS.minBias, optS.maxBias);
    
    % Update the UV map for bias
%     uvBias = reshape(biasPatch, 3, numUvValidPix);
    
    % Apply the bias correction to source patch
    srcPatch = bsxfun(@plus, srcPatch, uvBias);
end

% =========================================================================
% Compute weighted distance
% =========================================================================
% Compute distance
patchDist = trgPatch - srcPatch;

% Sum of absolute distance
if(strcmp(optS.costType, 'L1'))
    patchDist = abs(patchDist);
elseif(strcmp(optS.costType, 'L2'))
    patchDist = patchDist.^2;
end

% Apply weights
patchDist = bsxfun(@times, patchDist, optS.wPatch);
costApp  = squeeze(sum(sum(patchDist, 1),2));
end