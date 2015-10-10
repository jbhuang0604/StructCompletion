function [NNF, nUpdate]= sc_update_NNF(trgPatch, img, NNF, modelPlane, modelReg, optS)

% SC_UPDATE_NNF
%
% Update the nearest neighbor field using the PatchMatch algorithm
%
% Input:
%   - trgPatch
%   - img:
%   - NNF
%   - modelPlane
%   - modelReg
%   - optS
% Output:
%   - NNF
%   - nUpdate

nUpdate = zeros(1,3);

for i = 1:optS.numPassPerIter
    % propagate along four directions
    for iDirect = 1:4
        [NNF, n] = sc_propagate(trgPatch, img, NNF, modelPlane, optS, iDirect);
        nUpdate(1) = nUpdate(1) + n;
    end
    
    % Random sampling
    %     [NNF, n] = sc_regular_search(trgPatch, img, NNF, modelPlane, modelReg, optS);
    %     nUpdate(3) = nUpdate(3) + n;
    
    % Regularity guided sampling
    [NNF, n] = sc_random_search(trgPatch, img, NNF, modelPlane, optS);
    nUpdate(2) = nUpdate(2) + n;
end

end


% PatchMatch - propagate step
function [NNF, nUpdateTotal] = sc_propagate(trgPatch, img, NNF, modelPlane, optS, iDirect)

% SC_PROPAGATE: update the nearest neighbor field using propagation
% Input:
%   - trgPatch:     target patch
%   - img:          image at current level
%   - NNF:          current nearest neighbor field
%   - modelPlane:   planar structure model
%   - optS
%   - indDirection: propagation direction
% Output:
%   - NNF, nUpdateTotal

nUpdateTotal = 0;

% The positions of neighboring pixels
uvPixN = NNF.uvPixN{iDirect};
uvPixActiveInd = uvPixN.validInd;

numUpdatePix = NNF.uvPix.numPix;
while(numUpdatePix ~= 0)
    % Prepare uvPix, uvPixNCur
    uvPix.sub     = NNF.uvPix.sub(uvPixActiveInd, :);  uvPix.ind     = NNF.uvPix.ind(uvPixActiveInd, :);
    uvPixNCur.sub = uvPixN.sub(uvPixActiveInd, :);     uvPixNCur.ind = uvPixN.ind(uvPixActiveInd, :);
    
    uvDtBdPixPosCur = NNF.uvDtBdPixPos(uvPixActiveInd);
    
    trgPatchCur   = trgPatch(:,:, uvPixActiveInd);
    srcPosCur     = NNF.uvTform.data(uvPixActiveInd, 7:8);
    uvCostCur     = NNF.uvCost.data(uvPixActiveInd);
    uvPlaneIDCur  = NNF.uvPlaneID.data(uvPixActiveInd);
    
    % Get srcPosMap for computing coherence cost
    srcPosMapCur  = NNF.uvTform.map(:,:,7:8);
    
    % Active pixel positions
    uvPixActivePos = find(uvPixActiveInd);
    
    % Get candidate uvTform candidates
    uvTformCand = sc_uvMat_from_uvMap(NNF.uvTform.map, uvPixNCur.ind);
    
    % Generate candidate transformation by propagation
    uvTformCand = sc_trans_tform(uvTformCand, optS.propDir(iDirect,:));
    
    % Check if the nearest neighbors are valid source patches
    uvValidSrcInd = sc_check_valid_uv(uvTformCand(:,7:8), NNF.validPix.mask);
    
    % Check if the nearest neighbors are already the same as the existing one
    diff = abs(uvTformCand(:,7:8) - srcPosCur);
    uvValidDistInd = ((diff(:,1) > 1 ) | (diff(:,2) > 1 ));
    
    % Valid pixel indices
    uvValidInd = uvValidSrcInd & uvValidDistInd;
    
    numUvValid = sum(uvValidInd);
    
    if(numUvValid ~= 0)
        uvPixValid.sub = uvPix.sub(uvValidInd,:);
        uvPixValid.ind = uvPix.ind(uvValidInd);
        
        uvDtBdPixPosCur = uvDtBdPixPosCur(uvValidInd);
        
        trgPatchCur    = trgPatchCur(:,:, uvValidInd);
        uvTformCand    = uvTformCand(uvValidInd, :);
        uvCostCur      = uvCostCur(uvValidInd);
        uvPlaneIDCand  = uvPlaneIDCur(uvValidInd);
        
        uvPixUpdatePos = uvPixActivePos(uvValidInd);
        
        % Grab source patches
        srcPatch = sc_prep_source_patch(img, uvTformCand, optS);
        
        % Compute patch matching cost
        [costPatchCandAll, uvBiasCand] = ...
            sc_patch_cost(trgPatchCur, srcPatch, modelPlane, uvPlaneIDCand, ...
            uvPixValid.sub, uvTformCand, srcPosMapCur, uvDtBdPixPosCur, optS);
        costPatchCand = sum(costPatchCandAll, 2);
        
        % Check which one to update
        updateInd = costPatchCand < uvCostCur;
        
        uvPixUpdatePos = uvPixUpdatePos(updateInd);
        numUpdatePix = size(uvPixUpdatePos, 1);
    else
        numUpdatePix = 0;
    end
    
    % Update NNF data
    if(numUpdatePix ~= 0)
        nUpdateTotal = nUpdateTotal + numUpdatePix;
        
        % === Update NNF data ===
        NNF.uvTform.data(uvPixUpdatePos, :) = uvTformCand(updateInd,:);
        NNF.uvCost.data(uvPixUpdatePos)     = costPatchCand(updateInd);
        NNF.uvPlaneID.data(uvPixUpdatePos)  = uvPlaneIDCand(updateInd);
        
        % Apply bias correction
        if(optS.useBiasCorrection)
            NNF.uvBias.data(:,:,uvPixUpdatePos)   = uvBiasCand(:,:,updateInd);
        end
        %         NNF.update.data(uvPixUpdatePos)   = 1; % updateInd;
        
        % === Update NNF map ===
        uvPixValidInd = uvPixValid.ind(updateInd);
        NNF.uvTform.map    = sc_update_uvMap(NNF.uvTform.map, uvTformCand(updateInd,:), uvPixValidInd);
        NNF.uvCost.map     = sc_update_uvMap(NNF.uvCost.map, costPatchCand(updateInd), uvPixValidInd);
        NNF.uvPlaneID.map  = sc_update_uvMap(NNF.uvPlaneID.map, uvPlaneIDCand(updateInd), uvPixValidInd);
        
        %         NNF.update.map  = sc_update_uvMap(NNF.update.map, 1, uvPixValidInd);
        
        % === Update uvPixActiveInd ===
        uvPixNextSub = uvPixValid.sub(updateInd,:);
        uvPixNextSub = bsxfun(@plus, uvPixNextSub, optS.propDir(iDirect,:));
        uvPixNextInd = sub2ind([NNF.imgH, NNF.imgW], uvPixNextSub(:,2), uvPixNextSub(:,1));
        
        updateMap = NNF.uvPix.mask;
        updateMap(uvPixNextInd) = 0;
        uvPixActiveInd = ~updateMap(NNF.uvPix.ind);
        uvPixActiveInd = uvPixActiveInd & uvPixN.validInd;
    end
end

end

% PatchMatch - random search step
function [NNF, nUpdateTotal] = sc_random_search(trgPatch, img, NNF, modelPlane, optS)

% SC_RANDOM_SEARCH: update the nearest neighbor using random sampling

[imgH, imgW, nCh] = size(img);

uvPix    = NNF.uvPix;
numUvPix = size(uvPix.sub, 1);

searchRad = max(imgH, imgW)/2;
nUpdateTotal = 0;

% uvPixActiveInd = true(numUvPix, 1);
iter = 1;
% while(iter <= optS.numRandSample)
while(searchRad > 1)
    iter = iter + 1;
    % Reduce search radius by half
    searchRad = searchRad/2;
    if(searchRad < 1)
        break;
    end
    
    % Get srcPosMap for computing coherence cost
    srcPosMapCur  = NNF.uvTform.map(:,:,7:8);
    
    uvTformCandCur = sc_uvMat_from_uvMap(NNF.uvTform.map, uvPix.ind);
    
    % Draw random samples
    srcPos = uvTformCandCur(:,7:8) + 2*searchRad*(rand(numUvPix, 2) - 0.5);
    
    % Draw plane ID candidate
    uvPlaneIDCand = sc_draw_plane_id(NNF.uvPlaneID.planeProbAcc);
       
    % Estimate the domain transformation
    uvTformCand = sc_src_domain_tform(uvPlaneIDCand, modelPlane, [], srcPos, NNF.uvPix.sub, 1);
    
    % === Reject invalid samples ===
    % Check if the scale of the source patch valid
    uvTformScale    = sc_scale_tform(uvTformCand);
    uvValidScaleInd = (uvTformScale > optS.minScale) & (uvTformScale < optS.maxScale);
    % Check if the souce patch is valid
    uvValidSrcInd   = sc_check_valid_uv(uvTformCand(:,7:8), NNF.validPix.mask);
    
    uvValidInd = uvValidSrcInd & uvValidScaleInd;
    
    % Active pixels
    uvPixActivePos = find(uvValidInd);
    numActPix = size(uvPixActivePos, 1);
    
    if(numActPix~=0)
        % Update
        trgPatchCur      = trgPatch(:,:,uvValidInd);
        uvCostDataCur    = NNF.uvCost.data(uvValidInd);
        uvTformCandCur   = uvTformCand(uvValidInd, :);
        uvPlaneIDCandCur = uvPlaneIDCand(uvValidInd);
        
        uvPixValid.sub   = uvPix.sub(uvValidInd,:);
        uvPixValid.ind   = uvPix.ind(uvValidInd);
        
        uvDtBdPixPosCur = NNF.uvDtBdPixPos(uvValidInd);
        
        % Grab source patches
        srcPatch = sc_prep_source_patch(img, uvTformCandCur, optS);
        
        [costPatchCandAll, uvBiasCand] = ...
            sc_patch_cost(trgPatchCur, srcPatch, modelPlane, uvPlaneIDCandCur, ...
            uvPixValid.sub, uvTformCandCur, srcPosMapCur, uvDtBdPixPosCur, optS);
        costPatchCand = sum(costPatchCandAll, 2);
        % Check which one to update
        updateInd = (costPatchCand < uvCostDataCur);
        nUpdate = sum(updateInd);
        
        if(nUpdate~=0)
            uvPixActivePos = uvPixActivePos(updateInd);
            
            nUpdateTotal = nUpdateTotal + nUpdate;
            
            % === Update NNF data ===
            NNF.uvTform.data(uvPixActivePos, :) = uvTformCandCur(updateInd,:);
            NNF.uvPlaneID.data(uvPixActivePos)  = uvPlaneIDCandCur(updateInd);
            NNF.uvCost.data(uvPixActivePos)     = costPatchCand(updateInd);
            if(optS.useBiasCorrection)
                NNF.uvBias.data(:,:,uvPixActivePos)   = uvBiasCand(:,:,updateInd);
            end
            
            % === Update NNF map ===
            uvPixValidInd     = uvPixValid.ind(updateInd);
            NNF.uvTform.map   = sc_update_uvMap(NNF.uvTform.map, uvTformCandCur(updateInd,:), uvPixValidInd);
            NNF.uvPlaneID.map = sc_update_uvMap(NNF.uvPlaneID.map, uvPlaneIDCandCur(updateInd), uvPixValidInd);
            NNF.uvCost.map    = sc_update_uvMap(NNF.uvCost.map, costPatchCand(updateInd), uvPixValidInd);
        end
    end
end
end