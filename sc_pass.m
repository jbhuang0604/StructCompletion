function [img, NNF] = sc_pass(img, holeMask, NNF, wDistPatch, wDistImg, numIterLvl, modelPlaneCur, modelRegCur, imgSize, iLvl, optS, lockImgFlag)

% SC_PASS

% Modified PatchMatch algorithm

for iter = 1 : numIterLvl
    
    % === Compute the patch matching cost at the current level ===
    
    % Prepare target and source patches
    trgPatch = sc_prep_target_patch(img, NNF.uvPix.sub,    optS);
    srcPatch = sc_prep_source_patch(img, NNF.uvTform.data, optS);
    
    % Compute patch matching cost
    [uvCostcur, NNF.uvBias.data] = ...
        sc_patch_cost(trgPatch, srcPatch, wDistPatch, modelPlaneCur, NNF.uvPlaneID.data, ...
        NNF.uvPix.sub, NNF.uvTform.data(7:8,:), NNF.uvDtBdPixPos, NNF.uvPixUpdateSrc.data, imgSize, optS, iLvl);
    NNF.uvCost.data = sum(uvCostcur, 1);
    NNF.uvCost.map = sc_update_uvMap(NNF.uvCost.map, NNF.uvCost.data, NNF.uvPix, true(1, NNF.uvPix.numUvPix));
    
    % Initialize update index map (for early termination)
    NNF.update.data = false(1, NNF.uvPix.numUvPix);
    NNF.update.map  = false(NNF.imgH, NNF.imgW);
    
    % === Update the NNF using the PatchMatch algorithm ===
    [NNF, nUpdate]= sc_update_NNF(trgPatch, wDistPatch, img, NNF, modelPlaneCur, modelRegCur, iLvl, optS);
    avgPatchCost = mean(NNF.uvCost.data, 2);
    
    % === Update the image ===
    if(~lockImgFlag)
        img = sc_voting_update(wDistPatch, img, NNF, NNF.uvPix, holeMask, wDistImg, optS);
    end
    % === Visualizing the progress ===
    visFlag = 1; % Set to 0 for faster computation
    if(visFlag) 
        NNFVis = sc_vis_nnf(NNF);
        figure(1);
        subplot(2,2,1), imshow(img);
        title('Current Image', 'fontsize', 16);
        subplot(2,2,2), imshow(NNFVis.uvTfomMapVis);
        title('Nearest neighbor field', 'fontsize', 16);
        subplot(2,2,3), imshow(NNFVis.uvPlaneIDMapVis);
        title('Plane label (R, G, B, W)', 'fontsize', 16); colormap jet
        subplot(2,2,4), imshow(NNFVis.uvCostMapVis);
        title('Patch matching cost', 'fontsize', 16); colormap jet
    end
    fprintf('    %3d\t%12d\t%12d\t%12d\t%14f\n', iter, nUpdate(1), nUpdate(2), nUpdate(3), avgPatchCost);
end


end