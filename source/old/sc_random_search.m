function  [NNF, nUpdateTotal] = sc_random_search(trgPatch, wDistPatch, img, NNF, modelPlane, optS, iLvl)

% SC_RANDOM_SEARCH: update the nearest neighbor using random sampling

[imgH, imgW, nCh] = size(NNF.uvTform.map);

uvPix = NNF.uvPix;
numUvPix = size(uvPix.sub, 2);

searchRad = max(imgH, imgW)/2;
nUpdateTotal = 0;

uvPixActiveInd = true(1, numUvPix);
iter = 1;
while(iter <= optS.numRandSample)
    % while(searchRad > 1)
    iter = iter + 1;
    % Reduce search radius by half
    searchRad = searchRad/2;
    %     if(searchRad < 1)
    %         break;
    %     end
    
    uvTformCandCur = sc_uvMat_from_uvMap(NNF.uvTform.map, uvPix);
    
    % Draw random samples
    srcPos = uvTformCandCur(7:8,:) + 2*searchRad*(rand(2, numUvPix) - 0.5);
    srcPos(1,:) = sc_clamp(srcPos(1,:), optS.pRad+1, imgW - optS.pRad);
    srcPos(2,:) = sc_clamp(srcPos(2,:), optS.pRad+1, imgH - optS.pRad);
    
    % Draw plane ID candidate
    uvPlaneIDCand = sc_draw_plane_id(NNF.uvPlaneID.planeProbAcc);
    
    % Estimate the domain transformation
    uvTformCand = sc_src_domain_tform(uvPlaneIDCand, modelPlane, [], srcPos, NNF.uvPix.sub, 1);
    
    A = isnan(uvTformCand);
    if(sum(A(:))~=0)
        error('uvTformCand contains NaN');
    end
    % === Reject invalid samples ===
    % Check if the scale of the source patch valid
    uvTformScale = sc_scale_tform(uvTformCand);
    uvValidScaleInd = (uvTformScale > optS.minScale) & (uvTformScale < optS.maxScale);
    % Check if the souce patch is valid
    uvValidSrcInd = sc_check_valid_uv(uvTformCand(7:8,:), NNF.validPix.mask);
    % Check if the cost is already low
    uvValidCostInd = NNF.uvCost.data > optS.rsThres;
    
    uvValidInd = uvPixActiveInd & uvValidSrcInd & ...
        uvValidCostInd & uvValidScaleInd;
    
    uvPixActivePos = find(uvValidInd);
    numActPix = size(uvPixActivePos, 2);
    
    if(numActPix~=0)
        % Update
        trgPatchCur      = trgPatch(:,:,uvValidInd);
        wDistPatchCur    = wDistPatch(:, uvValidInd);
        uvCostDataCur    = NNF.uvCost.data(:,uvValidInd);
        uvTformCandCur   = uvTformCand(:, uvValidInd);
        uvPlaneIDCandCur = uvPlaneIDCand(uvValidInd);
        
        uvPixValid.sub = uvPix.sub(:,uvValidInd);
        uvPixValid.ind = uvPix.ind(uvValidInd);
        
        uvDtBdPixPosCur = NNF.uvDtBdPixPos(:, uvValidInd);
        
        % Grab source patches
        srcPatch = sc_prep_source_patch(img, uvTformCandCur, optS);
        
        [costPatchCandAll, uvBiasCand] = ...
            sc_patch_cost(trgPatchCur, srcPatch, wDistPatchCur, modelPlane, uvPlaneIDCandCur, ...
            uvPixValid.sub, uvTformCandCur(7:8,:), uvDtBdPixPosCur, zeros(1, numActPix), imgSize, optS, iLvl);
        costPatchCand = sum(costPatchCandAll, 1);
        % Check which one to update
        updateInd = (costPatchCand < uvCostDataCur);
        nUpdate = sum(updateInd);
        
        if(nUpdate~=0)
            uvPixActivePos = uvPixActivePos(updateInd);
            
            uvPixActiveInd(uvPixActivePos) = 0;
            
            nUpdateTotal = nUpdateTotal + nUpdate;
            
            % === Update NNF data ===
            NNF.uvTform.data(:, uvPixActivePos) = uvTformCandCur(:,updateInd);
            NNF.uvPlaneID.data(uvPixActivePos)  = uvPlaneIDCandCur(updateInd);
            NNF.uvCost.data(uvPixActivePos)     = costPatchCand(updateInd);
            if(optS.useBiasCorrection)
                NNF.uvBias.data(:,uvPixActivePos)   = uvBiasCand(:,updateInd);
            end
            NNF.update.data(uvPixActivePos) = 2;
            
            %
            NNF.uvPixUpdateSrc.data(uvPixActivePos) = 2;
            
            % === Update NNF map ===
            NNF.uvTform.map = sc_update_uvMap(NNF.uvTform.map, uvTformCandCur(:,updateInd), uvPixValid, updateInd);
            NNF.uvPlaneID.map = sc_update_uvMap(NNF.uvPlaneID.map, uvPlaneIDCandCur(updateInd), uvPixValid, updateInd);
            NNF.uvCost.map  = sc_update_uvMap(NNF.uvCost.map, costPatchCand(updateInd), uvPixValid, updateInd);
            if(optS.useBiasCorrection)
                NNF.uvBias.map  = sc_update_uvMap(NNF.uvBias.map, uvBiasCand(:,updateInd), uvPixValid, updateInd);
            end
            NNF.update.map  = sc_update_uvMap(NNF.update.map, 1, uvPixValid, updateInd);
            NNF.uvPixUpdateSrc.map  = sc_update_uvMap(NNF.uvPixUpdateSrc.map, 2, uvPixValid, updateInd);
        end
    end
end


end