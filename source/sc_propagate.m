function [NNF, nUpdateTotal] = sc_propagate(trgPatch, wDistPatch, img, NNF, modelPlane, optS, indDirection, iLvl)

% SC_PROPAGATE: update the nearest neighbor field using propagation

% Input:
%   - trgPatch, wDistPatch, img, NNF, modelPlane, optS, indDirection
% Output:
%   - NNF, nUpdateTotal

[imgH, imgW, nCh] = size(img);
imgSize = max(imgH, imgW);

nUpdateTotal = 0;

% The positions of neighboring pixels
uvPixN = NNF.uvPixN{indDirection};
uvPixActiveInd = true(1, NNF.uvPix.numUvPix);
uvPixActiveInd = uvPixActiveInd & uvPixN.validInd;

numUpdatePix = NNF.uvPix.numUvPix;
while(numUpdatePix ~= 0)
    % Prepare uvPix, uvPixNCur
    uvPix.sub     = NNF.uvPix.sub(:, uvPixActiveInd); uvPix.ind     = NNF.uvPix.ind(:, uvPixActiveInd);
    uvPixNCur.sub = uvPixN.sub(:, uvPixActiveInd);    uvPixNCur.ind = uvPixN.ind(:, uvPixActiveInd);
    uvDtBdPixPosCur = NNF.uvDtBdPixPos(:, uvPixActiveInd);
    trgPatchCur   = trgPatch(:,:, uvPixActiveInd);
    wDistPatchCur = wDistPatch(:,uvPixActiveInd);
    srcPosCur     = NNF.uvTform.data(7:8, uvPixActiveInd);    
    uvCostCur     = NNF.uvCost.data(:, uvPixActiveInd);
    uvPlaneIDCur  = NNF.uvPlaneID.map(uvPixNCur.ind);
    uvPixUpdateSrc= NNF.uvPixUpdateSrc.data(uvPixActiveInd);
    
    % Active pixel positions
    uvPixActivePos = find(uvPixActiveInd);

    % Get candidate uvTform candidates
    uvTformCand = sc_uvMat_from_uvMap(NNF.uvTform.map, uvPixNCur);
    
    % Generate candidate transformation by propagation
    uvTformCand = sc_trans_tform(uvTformCand, optS.propDir(:,indDirection));
    
    % Check if the nearest neighbors are valid source patches
    uvValidSrcInd = sc_check_valid_uv(uvTformCand(7:8,:), NNF.validPix.mask);
    % Check if the nearest neighbors are already the same as the existing one
    diff = abs(uvTformCand(7:8,:) - srcPosCur);
    uvValidDistInd = ((diff(1,:) > 1 ) | (diff(2,:) > 1 ));
    
    % Valid pixel indices
    uvValidInd = uvValidSrcInd & uvValidDistInd;
    
    numUvValid = sum(uvValidInd);
    
    if(numUvValid ~= 0)
        trgPatchCur    = trgPatchCur(:,:, uvValidInd);
        wDistPatchCur  = wDistPatchCur(:,uvValidInd);
        uvTformCand    = uvTformCand(:, uvValidInd);        
        uvCostCur = uvCostCur(uvValidInd);
        
        uvPixUpdatePos = uvPixActivePos(uvValidInd);
        uvPixValid.sub = uvPix.sub(:,uvValidInd);
        uvPixValid.ind = uvPix.ind(uvValidInd);
        uvPlaneIDCand  = uvPlaneIDCur(uvValidInd);
        
        uvDtBdPixPosCur = uvDtBdPixPosCur(:, uvValidInd);
        
        % Grab source patches
        srcPatch = sc_prep_source_patch(img, uvTformCand, optS);
        
        % Compute patch matching cost
        [costPatchCandAll, uvBiasCand] = ...
            sc_patch_cost(trgPatchCur, srcPatch, wDistPatchCur, modelPlane, uvPlaneIDCand, ...
            uvPixValid.sub, uvTformCand(7:8,:), uvDtBdPixPosCur, zeros(1, numUvValid), imgSize, optS, iLvl);
        costPatchCand = sum(costPatchCandAll, 1);

        % Check which one to update
        updateInd = costPatchCand < uvCostCur;
        
        uvPixUpdatePos = uvPixUpdatePos(updateInd);
        numUpdatePix = size(uvPixUpdatePos, 2);
    else
        numUpdatePix = 0;
    end
    
    % Update NNF data
    if(numUpdatePix ~= 0)
        nUpdateTotal = nUpdateTotal + numUpdatePix;
        
        % === Update NNF data ===
        NNF.uvTform.data(:, uvPixUpdatePos) = uvTformCand(:,updateInd);
        NNF.uvCost.data(uvPixUpdatePos)     = costPatchCand(updateInd);        
        NNF.uvPlaneID.data(uvPixUpdatePos)  = uvPlaneIDCand(updateInd);
        
        % Apply bias correction
        if(optS.useBiasCorrection)
            NNF.uvBias.data(:,uvPixUpdatePos)   = uvBiasCand(:,updateInd);
        end
        NNF.update.data(:,uvPixUpdatePos)   = 1; % updateInd;
        
        % Label as update by propagation
        NNF.uvPixUpdateSrc.data(uvPixUpdatePos) = 3;
        
        % === Update NNF map ===
        NNF.uvTform.map    = sc_update_uvMap(NNF.uvTform.map, uvTformCand(:,updateInd), uvPixValid, updateInd);
        NNF.uvCost.map     = sc_update_uvMap(NNF.uvCost.map, costPatchCand(updateInd), uvPixValid, updateInd);
        NNF.uvPlaneID.map  = sc_update_uvMap(NNF.uvPlaneID.map, uvPlaneIDCand(updateInd), uvPixValid, updateInd);
        if(optS.useBiasCorrection)
            NNF.uvBias.map  = sc_update_uvMap(NNF.uvBias.map, uvBiasCand(:,updateInd), uvPixValid, updateInd);
        end
        NNF.update.map  = sc_update_uvMap(NNF.update.map, 1, uvPixValid, updateInd);
        NNF.uvPixUpdateSrc.map = sc_update_uvMap(NNF.uvPixUpdateSrc.map, 3, uvPixValid, updateInd);
        
        % === Update uvPixActiveInd ===
        uvPixNextSub = uvPixValid.sub(:,updateInd);
        uvPixNextSub(1,:) = uvPixNextSub(1,:) + optS.propDir(1,indDirection);
        uvPixNextSub(2,:) = uvPixNextSub(2,:) + optS.propDir(2,indDirection);
        uvPixNextInd = sub2ind([NNF.imgH, NNF.imgW], uvPixNextSub(2,:), uvPixNextSub(1,:));
        
        updateMap = NNF.uvPix.mask;
        updateMap(uvPixNextInd) = 0;
        uvPixActiveInd = ~updateMap(NNF.uvPix.ind);
        uvPixActiveInd = uvPixActiveInd & uvPixN.validInd;
    end
end

end