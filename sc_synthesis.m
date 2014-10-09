function [imgPyr, imgPyrNNF] = ...
    sc_synthesis(imgPyr, maskPyr, modelPlane, modelReg, optS)

% SC_SYNTHESIS:
%
% Patch-based synthesis using patchmatch algorithm and structural guidance
%
% Input:
%   - imgPyr:
%   - maskPyr:
%   - modelPlane:
%   - modelReg:
%   - optS:
% Output:
%   - imgPyr
%   - imgPyrNNF
%

imgPyrNNF = cell(optS.numPyrLvl, 1);
NNF = [];
numIterLvl = optS.numIter;
pyrLvl = optS.numPyrLvl: -1 : optS.topLevel;


for iLvl = pyrLvl
    % Initialize level
    holeMask = maskPyr{iLvl};
    [imgH, imgW] = size(holeMask);
    imgSize = max(imgH, imgW);
    
    modelPlaneCur = modelPlane{iLvl};
    modelRegCur   = modelReg{iLvl};
    
    % === Prepare img and NNF for the current level ===
    fprintf('--- Initialize NNF: ');
    tic;
    [img, NNF, wDistPatch, wDistImg] = ...
        sc_init_lvl_nnf(imgPyr{iLvl}, NNF, holeMask, modelPlaneCur, modelRegCur, iLvl, optS);
    tNNF = toc;
    fprintf('done in %6.3f seconds.\n', tNNF);
    
    % Number of iterations at the currect level
    numIterLvl = max(numIterLvl - optS.numIterDec, optS.numIterMin);
    
    fprintf('--- Pass... level: %d, #Iter: %d, #uvPixels: %7d\n', iLvl, numIterLvl, NNF.uvPix.numUvPix);
    fprintf('--- %3s\t%12s\t%12s\t%12s\t%10s\n', 'iter', '#PropUpdate', '#RandUpdate', '#RegUpdate', 'AvgCost');
    
    for iter = 1 : numIterLvl
        
        % === Compute the patch matching cost at the current level ===
        % Prepare target and source patches
        trgPatch = sc_prep_target_patch(img, NNF.uvPix.sub,    optS);
        srcPatch = sc_prep_source_patch(img, NNF.uvTform.data, optS);
        
        % Compute patch matching cost
        [NNF.uvCost.data, NNF.uvBias.data] = ...
            sc_patch_cost(trgPatch, srcPatch, wDistPatch, modelPlaneCur, NNF.uvPlaneID.data, ...
            NNF.uvPix.sub, NNF.uvTform.data(7:8,:), NNF.uvDtBdPixPos, NNF.uvPixUpdateSrc.data, imgSize, optS, iLvl);
        NNF.uvCost.map = sc_update_uvMap(NNF.uvCost.map, NNF.uvCost.data, NNF.uvPix, true(1, NNF.uvPix.numUvPix));
        
        % Initialize update index map (for early termination)
        NNF.update.data = false(1, NNF.uvPix.numUvPix);
        NNF.update.map  = false(NNF.imgH, NNF.imgW);
        
        % === Update the NNF using the PatchMatch algorithm ===
        [NNF, nUpdate]= sc_update_NNF(trgPatch, wDistPatch, img, NNF, modelPlaneCur, modelRegCur, iLvl, optS);
        avgPatchCost = mean(NNF.uvCost.data, 2);
        
        % === Update the image ===
        img = sc_voting_update(wDistPatch, img, NNF, NNF.uvPix, holeMask, wDistImg, optS);
        
        % === Visualizing the progress ===
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
        
        fprintf('    %3d\t%12d\t%12d\t%12d\t%14f\n', iter, nUpdate(1), nUpdate(2), nUpdate(3), avgPatchCost);
    end
    
    imgPyr{iLvl} = img;
    imgPyrNNF{iLvl} = NNF;
end

% Visualization
NNFVis = sc_vis_nnf(NNF);

figure(2);
subplot(1,2,1), imshow(img); title('Completion results before blending');
subplot(1,2,2), imshow(NNFVis.uvTfomMapVis); title('Nearest neighor field');

% To-Do: Poisson blending
% imgFinal = sc_poisson_blend(img, holeMask, distMap, optS);
%
% figure(3),
% subplot(1,3,1); imshow(img);
% subplot(1,3,2); imshow(imgFinal);
% subplot(1,3,3); imshow(imgFinal-img+0.5);

end