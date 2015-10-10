function imgPyr = ...
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

NNF = [];
numIterLvl = optS.numIter;
pyrLvl = optS.numPyrLvl: -1 : optS.topLevel;

% Coarse-to-fine image completion
for iLvl = pyrLvl
    
    % =====================================================================
    % Initialize level
    % =====================================================================
    holeMask = maskPyr{iLvl};
    [imgH, imgW] = size(holeMask);
    
    % Update optS with current level and imgSize
    optS.iLvl = iLvl;
    optS.imgSize = max(imgH, imgW);
    
    % Use the plane model and regularity model at the current level
    modelPlaneCur = modelPlane{iLvl};
    modelRegCur   = modelReg{iLvl};
    
    % Prepare img and NNF for the current level
    fprintf('--- Initialize NNF: ');
    tic;
    [img, NNF] = sc_init_lvl_nnf(imgPyr{iLvl}, NNF, holeMask, ...
        modelPlaneCur, modelRegCur, optS);
    tNNF = toc;
    fprintf('done in %6.3f seconds.\n', tNNF);
    
    % Number of iterations at the currect level
    numIterLvl = max(numIterLvl - optS.numIterDec, optS.numIterMin);
    optS.numIterLvl = numIterLvl;
    
    % =====================================================================
    % Optimization at the current level
    % =====================================================================
    fprintf('--- Pass... level: %d, #Iter: %d, #uvPixels: %7d\n', iLvl, numIterLvl, NNF.uvPix.numPix);
    fprintf('--- %3s\t%12s\t%12s\t%12s\t%10s\n', 'iter', '#PropUpdate', '#RandUpdate', '#RegUpdate', 'AvgCost');
    
    if(optS.iLvl == optS.numPyrLvl)
        lockImgFlag = 1;
        [img, NNF] = sc_pass(img, holeMask, NNF, modelPlaneCur, modelRegCur, optS, lockImgFlag);
        lockImgFlag = 0;
        [img, NNF] = sc_pass(img, holeMask, NNF, modelPlaneCur, modelRegCur, optS, lockImgFlag);
    else
        [img, NNF] = sc_pass(img, holeMask, NNF, modelPlaneCur, modelRegCur, optS, lockImgFlag);
    end
    
    
    % =====================================================================
    % (Optional) visualizing progress
    % =====================================================================
    
    visFlag = 0;
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
    
    % Save the result
    imgPyr{iLvl} = img;
end

% Visualization
% NNFVis = sc_vis_nnf(NNF);

% figure(2);
% subplot(1,2,1), imshow(img); title('Completion results before blending');
% subplot(1,2,2), imshow(NNFVis.uvTfomMapVis); title('Nearest neighor field');

end