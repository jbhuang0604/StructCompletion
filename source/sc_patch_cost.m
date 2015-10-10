function [costPatchCand, uvBiasCand] = ...
    sc_patch_cost(trgPatch, srcPatch, modelPlane, uvPlaneIDData, ...
    trgPos, srcTform, srcPosMap, bdPos, optS)

numUvPix      = size(srcPatch, 3);
costPatchCand = zeros(numUvPix, 4, 'single');

% Patch cost - appearance
[costApp, uvBiasCand] = sc_patch_cost_app(trgPatch, srcPatch, optS);

% Patch cost - planar guidance
srcPos = srcTform(:,7:8);
costPlane = sc_patch_cost_plane(modelPlane.mLogLPlaneProb, uvPlaneIDData, trgPos, srcPos);

if(optS.useCoherence) % Patch cost - coherence
    costCoherence = sc_patch_cost_coherence(trgPos, srcPosMap, srcTform, optS);
end

% Patch cost - directional guidance
costDirect = sc_patch_cost_direct(uvPlaneIDData, trgPos, srcPos, modelPlane, optS);

% Patch cost - proximity cost
costProx = sc_patch_cost_prox(srcPos, trgPos, bdPos, optS);

% costReg = double(uvPixUpdateSrc == 1);

% Weighted sum of the appearance cost and guidance cost
costPatchCand(:,1) = costApp;
costPatchCand(:,2) = optS.lambdaPlane * costPlane;

if(optS.useCoherence)
    costPatchCand(:,3) = optS.lambdaCoherence * costCoherence;
end

costPatchCand(:,4) = optS.lambdaDirect * costDirect;
costPatchCand(:,5) = optS.lambdaProx   * costProx;
% costPatchCand(3,:) = optS.lambdaReg    * costReg;

% [To-Do] Adaptive weighting for patch cost
% costPatchCand(2:5,:) = (iLvl/optS.numPyrLvl)*costPatchCand(2:5,:);

end

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

function cost = sc_patch_cost_coherence(trgPos, srcPosMap, srcTform, opt)

% SC_PATCH_COST_COHERENCE: spatial coherence cost
%
% Input
%   - trgPos:    target patch positions [numPix] x [2]
%   - srcPosMap: source patch map [imgH] x [imgW] x [2]
%   - srcPos:    source patch positions [numPix] x [2]
%   - srcTfmG:   source patch geometric transformation [numPix] x [9]
%   - opt:       parameters
% Output:
%   - cost:      spatio-temporal coherence cost

numPix = size(trgPos,1);
[imgH, imgW, ~] = size(srcPosMap);

% initialize coherence cost
cost = zeros(numPix, 1, 'single');

% =========================================================================
% Compute spatial coherence cost
% =========================================================================
for i = 1:4
    % source patch position prediction
    v = opt.propDir(i,:);
    srcPosP  = sc_trans_tform(srcTform, v);
    srcPosP =  srcPosP(:,7:8);
    
    % source patch positions of neighboring target patches
    trgPosN = bsxfun(@plus, trgPos, v);
    trgIndN = uint32(sub2ind([imgH, imgW], trgPosN(:,2), trgPosN(:,1)));
    srcPosN = sc_uvMat_from_uvMap(srcPosMap,  trgIndN);
    
    validSrc = (srcPosN(:,1) ~= 0);
    
    % add cost if the differences are high
    diff           = zeros(numPix, 1, 'single');
    diff(validSrc) = sum(abs(srcPosP(validSrc,:) - srcPosN(validSrc,:)), 2) > 1;
    cost    = cost + diff;
end

% cost = opt.lambdaCoherence*cost;

end

function costPlane = sc_patch_cost_plane(mLogLPlaneProb, uvPlaneIDData, trgPixSub, srcPixSub)

% SC_PATCH_COST_PLANE
%
% Compute planar costs (See Eqn 11 in the paper)

[imgH, imgW, numPlane] = size(mLogLPlaneProb);

srcPixSub     = round(srcPixSub);
uvPlaneIDData = double(uvPlaneIDData);

trgPixIndCur = sub2ind([imgH, imgW, numPlane], trgPixSub(:,2), trgPixSub(:,1), uvPlaneIDData);
srcPixIndCur = sub2ind([imgH, imgW, numPlane], srcPixSub(:,2), srcPixSub(:,1), uvPlaneIDData);

costPlane = mLogLPlaneProb(trgPixIndCur) + mLogLPlaneProb(srcPixIndCur);

end

function costProx = sc_patch_cost_prox(srcPos, trgPos, uvDtBdPixPos, optS)

% SC_PATCH_COST_PROX

% Encourage to sample patches near to the target patch

d = srcPos - trgPos;
d = sqrt(sum(d.^2,2));

d = d/optS.imgSize;
uvDtBdPixPos = uvDtBdPixPos/optS.imgSize;

% costProx = max(0, d - optS.proxThres);
% Shrinkage thresholding
costProx = max(0, d - uvDtBdPixPos - optS.proxThres);
% costProx = max(0, d - optS.proxThres);

end


function costDirect = sc_patch_cost_direct(uvPlaneIDData, trgPos, srcPos, modelPlane, optS)

% SC_PATCH_COST_DIRECT
%
% Compute the directional cost (See Eqn 13 in the paper)
%
% Input
%   -
% Output
%   -

numUvPix = size(trgPos, 1);
costDirect = optS.lambdaDirect*ones(numUvPix, 2, 'single');

for indPlane = 1: modelPlane.numPlane
    % Retrieve the uv pixels that have the current plane label
    uvPlaneIndCur  = uvPlaneIDData == indPlane;
    numPlanePixCur = sum(uvPlaneIndCur);
    
    if(indPlane == modelPlane.numPlane)
        costDirect(uvPlaneIndCur, :) = optS.imgSize*optS.directThres;
    else
        % The rectifying transformation for the plane
        rectMat = modelPlane.rectMat{indPlane};
        h7 = rectMat(3,1);
        h8 = rectMat(3,2);
        
        if(numPlanePixCur~=0)
            trgPosCur = trgPos(uvPlaneIndCur, :) - 1;
            srcPosCur = srcPos(uvPlaneIndCur, :) - 1;
            
            for iTheta = 1:2
                rotMat = modelPlane.rotMat{indPlane, iTheta};
                
                rotRecMat = rotMat;
                rotRecMat(3,1) = h7;    rotRecMat(3,2) = h8;
                rotRecMat = rotRecMat';
                
                trgPosCurRect = cat(2, trgPosCur, ones(numPlanePixCur,1))*rotRecMat;
                trgPosCurRect = bsxfun(@rdivide, trgPosCurRect, trgPosCurRect(:,3));
                
                % Source patch center position in the rectified domain
                srcPosCurRect = cat(2, srcPosCur, ones(numPlanePixCur,1))*rotRecMat;
                srcPosCurRect = bsxfun(@rdivide, srcPosCurRect, srcPosCurRect(:,3));
                
                costDirect(uvPlaneIndCur, iTheta) = abs(srcPosCurRect(:,2) - trgPosCurRect(:,2));
            end
        end
    end
end

costDirect = min(costDirect, [], 2);
costDirect = min(costDirect/optS.imgSize, optS.directThres);

end