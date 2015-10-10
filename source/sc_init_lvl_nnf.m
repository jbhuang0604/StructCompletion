% SC_INIT_LVL_NNF
function [img, NNF] = sc_init_lvl_nnf(img, NNF, holeMask, modelPlaneCur, modelRegCur, optS)
% SC_INIT_LVL_NNF: Initialize the nearest neighbor field for the current level
%
% Input:
%   - img:           Image at the current level
%   - NNF:           Nearest neigbhbor field
%   - holeMask:      Binary hole mask at the current level
%   - modelPlaneCur: planar structure model
%   - modelRegCur:   regularity structure model
%   - optS:          parameter for patch-based synthesis
% Output:
%   - img:           updated image
%   - NNF:           updated nearest neighbor field

if(optS.iLvl == optS.numPyrLvl)
    % Initialize the NNF for the coarest level using random sampling
    NNF = sc_init_nnf(holeMask, modelPlaneCur, modelRegCur, optS);
else
    % Initialize the NNF upsampling of NNF from previous level
    NNF = sc_upsample(holeMask, NNF, modelPlaneCur, modelRegCur, optS);
    
    trgPatch = sc_prep_target_patch(img, NNF.uvPix.sub,    optS);
    srcPatch = sc_prep_source_patch(img, NNF.uvTform.data, optS);
    
    % Compute patch matching cost
    [~, NNF.uvBias.data] = ...
        sc_patch_cost(trgPatch, srcPatch, modelPlaneCur, NNF.uvPlaneID.data, ...
        NNF.uvPix.sub, NNF.uvTform.data, NNF.uvTform.map(:,:,7:8), NNF.uvDtBdPixPos, optS);   
    
    % Synthesize image using the upsampled NNF 
    img = sc_voting(img, NNF, holeMask, optS);
end

end

% SC_INIT_NNF
function NNF = sc_init_nnf(holeMask, modelPlane, modelReg, optS)

% SC_INIT_NNF: Initialize the nearest neighbor field
% Input:
%     - holeMask: the binary mask indicate the hole pixels
%     - modelPlane: model of detected plane structures
%     - modelReg: model of detected reguarity structures
%     - optS: parameters for synthesis process
% Output:
%   Image height and width
%     - NNF.imgH: image height
%     - NNF.imgW: image width
%   Precomputing pixel positions for PatchMatch algorithm
%     - NNF.uvPix: the patch positions containing hole pixels
%     - NNF.uvPixN: the neighboring patch positions of 4 connected pixels
%     - NNF.validPix: the patch positions containing known pixels
%     - NNF.trgPatchInd: the indice to retrieve target patches
%   Initialize NNF components:
%     - NNF.uvTform: the transformation specifying the source patch
%         - NNF.uvTform.data:   numUvPix x 9
%         - NNF.uvTform.map:    H x W x 9
%     - NNF.uvPlaneID: the plane assignment for unknown region
%         - NNF.uvPlaneID.data: numUvPix x 1
%         - NNF.uvPlaneID.map:  H x W
%     - NNF.uvBias: the bias between source and target patches
%         - NNF.uvBias.data:    numUvPix x 3
%         - NNF.uvBias.map:     H x W x 3
%     - NNF.uvGain: the gain between source and target patches
%         - NNF.uvGain.data:    numUvPix x 3
%         - NNF.uvGain.map:     H x W x 3
%     - NNF.uvCost: the matching cost between source and target patches
%         - NNF.uvCost.data:    numUvPix x 3
%         - NNF.uvCost.map:     H x W x 3

% =========================================================================
% Initialize uvPix and validPix
% =========================================================================
% Image size for the nearest neighbor field
[NNF.imgH, NNF.imgW] = size(holeMask);
% Get validPix and uvPix
[NNF.validPix, NNF.uvPix] = sc_init_level(holeMask, optS.pSize, optS.pRad);

% =========================================================================
% Initialize uvPixN
% =========================================================================
% 4-connected pixels for spatial propagation
%     2
%     |
% 1 - x - 3
%     |
%     4
NNF.uvPixN = getUvPixN(NNF, optS);

% =========================================================================
% Initialize indMap
% =========================================================================
NNF.trgPatchInd = getTargPathInd(NNF, optS);

% =========================================================================
% Initialize uvPlaneID
% =========================================================================

% Draw a random sample from U[0,1] for each uv pixel to sample its plane ID
% uvPlaneIDData= sc_draw_plane_id(planeProbAccData);
uvPlaneIDData = modelPlane.numPlane*ones(NNF.uvPix.numPix, 1, 'uint8');

% Save data
NNF.uvPlaneID.data = uvPlaneIDData;
NNF.uvPlaneID.map  = zeros(NNF.imgH, NNF.imgW, 'uint8');
NNF.uvPlaneID.map(NNF.uvPix.ind) = NNF.uvPlaneID.data;
NNF.uvPlaneID.planeProbAcc = sc_prep_plane_prob_acc(modelPlane.planeProb, NNF.uvPix.ind);
NNF.uvPlaneID.mLogLikelihood = -log(NNF.uvPlaneID.planeProbAcc);

% =========================================================================
% Initialize uvTform
% =========================================================================
% Initialize unTform.data with random samples
if(NNF.validPix.numPix)
    % Random sampling
    randInd   = randi(NNF.validPix.numPix, NNF.uvPix.numPix, 1);
    uvRandSub = NNF.validPix.sub(randInd, :);
else
    uvRandSub = (NNF.imgW/2)*ones(NNF.uvPix.numPix, 2);
end
NNF.uvTform.data = sc_src_domain_tform(NNF.uvPlaneID.data, modelPlane, modelReg, uvRandSub, NNF.uvPix.sub, 1);
NNF.uvTform.map = zeros(NNF.imgH, NNF.imgW, 9, 'single');
NNF.uvTform.map = sc_update_uvMap(NNF.uvTform.map, NNF.uvTform.data, NNF.uvPix.ind);

% =========================================================================
% Initialize uvBias
% =========================================================================
NNF.uvBias.data = zeros(1, 3, NNF.uvPix.numPix, 'single');
NNF.uvBias.map  = zeros(NNF.imgH, NNF.imgW, 3, 'single');

% =========================================================================
% Initialize uvCost
% =========================================================================
NNF.uvCost.data = zeros(NNF.uvPix.numPix, 1, 'single');
NNF.uvCost.map  = zeros(NNF.imgH, NNF.imgW, 'single');

% =========================================================================
% Initialize update
% =========================================================================
% NNF.update.data = false(NNF.uvPix.numPix, 1);
% NNF.update.map  = false(NNF.imgH, NNF.imgW);

% =========================================================================
% Initialize distMap, wPatchR, NNF.wPatchR
% =========================================================================

[NNF.distMap, ~] = bwdist(~holeMask, 'euclidean');
[NNF.wPatchR, NNF.wPatchSumImg] = sc_prep_dist_patch(NNF.distMap, NNF.uvPix.sub, optS);

% %% === Initialize uvPixUpdateSrc ===
% NNF.uvPixUpdateSrc.data = zeros(NNF.uvPix.numUvPix);
% NNF.uvPixUpdateSrc.map  = zeros(NNF.imgH, NNF.imgW);

NNF.uvDtBdPixPos = double(NNF.distMap(NNF.uvPix.ind));

end

% SC_UPSAMPLE
function NNF_H = sc_upsample(holeMask, NNF_L, modelPlane, modelReg, optS)
% SC_UPSAMPLE: upsample the nearest neighbor field

% Input:
%     - holeMask: the binary mask indicate the hole pixels
%     - NNF_L: nearest neighbor field of the low resolution image
%     - optS: parameters for synthesis process
% Output:
%   Image height and width
%     - NNF_H.imgH: image height
%     - NNF_H.imgW: image width
%   Precomputing pixel positions for PatchMatch algorithm
%     - NNF_H.uvPix: the patch positions containing hole pixels
%     - NNF_H.uvPixN: the neighboring patch positions of 4 connected pixels
%     - NNF_H.validPix: the patch positions containing known pixels
%     - NNF_H.trgPatchInd: the indice to retrieve target patches
%   Initialize NNF components:
%     - NNF_H.uvTform: the transformation specifying the source patch
%         - NNF.uvTform.data:   9 x numUvPix
%         - NNF.uvTform.map:    H x W x 9
%     - NNF_H.uvPlaneID: the plane assignment for unknown region
%         - NNF.uvPlaneID.data: 1 x numUvPix
%         - NNF.uvPlaneID.map:  H x W
%     - NNF_H.uvBias: the bias between source and target patches
%         - NNF.uvBias.data:    3 x numUvPix
%         - NNF.uvBias.map:     H x W x 3
%     - NNF_H.uvCost: the matching cost between source and target patches
%         - NNF.uvCost.data:    3 x numUvPix
%         - NNF.uvCost.map:     H x W x 3

% =========================================================================
% Initialize uvPix and validPix
% =========================================================================
[NNF_H.imgH, NNF_H.imgW] = size(holeMask);
[NNF_H.validPix, NNF_H.uvPix] = sc_init_level(holeMask, optS.pSize, optS.pRad);

% =========================================================================
% Initialize uvPixN
% =========================================================================
NNF_H.uvPixN = getUvPixN(NNF_H, optS);

% =========================================================================
% Initialize trgPatchInd
% =========================================================================
NNF_H.trgPatchInd = getTargPathInd(NNF_H, optS);

% =========================================================================
% Initialize uvPixL for upsampling
% =========================================================================

imgH_H = NNF_H.imgH;    imgW_H = NNF_H.imgW;
imgH_L = NNF_L.imgH;    imgW_L = NNF_L.imgW;

sX = imgH_L/imgH_H;     sY = imgW_L/imgW_H;
uvPixL.sub = round(NNF_H.uvPix.sub*diag([sX, sY]));
uvPixL.sub(:,1) = sc_clamp(uvPixL.sub(:,1), optS.pRad+1, imgW_L - optS.pRad);
uvPixL.sub(:,2) = sc_clamp(uvPixL.sub(:,2), optS.pRad+1, imgH_L - optS.pRad);
uvPixL.ind = uint32(sub2ind([imgH_L, imgW_L], uvPixL.sub(:,2), uvPixL.sub(:,1)));

% =========================================================================
% Initialize uvPlaneID
% =========================================================================

NNF_H.uvPlaneID.data = sc_uvMat_from_uvMap(NNF_L.uvPlaneID.map, uvPixL.ind);
NNF_H.uvPlaneID.data(NNF_H.uvPlaneID.data==0) = 1;
NNF_H.uvPlaneID.map = zeros(NNF_H.imgH, NNF_H.imgW, 'uint8');
NNF_H.uvPlaneID.map = sc_update_uvMap(NNF_H.uvPlaneID.map, NNF_H.uvPlaneID.data, NNF_H.uvPix.ind);

NNF_H.uvPlaneID.planeProbAcc = sc_prep_plane_prob_acc(modelPlane.planeProb, NNF_H.uvPix.ind);
NNF_H.uvPlaneID.mLogLikelihood = -log(NNF_H.uvPlaneID.planeProbAcc);


% === Initialize uvPixUpdateSrc ===
% NNF_H.uvPixUpdateSrc.map = zeros(NNF_H.imgH, NNF_H.imgW);
% NNF_H.uvPixUpdateSrc.data = sc_uvMat_from_uvMap(NNF_L.uvPixUpdateSrc.map, uvPixL);
% NNF_H.uvPixUpdateSrc.map = sc_update_uvMap(NNF_H.uvPixUpdateSrc.map, NNF_H.uvPixUpdateSrc.data, NNF_H.uvPix, true(1,NNF_H.uvPix.numUvPix));

% =========================================================================
% Update uvTform
% =========================================================================

uvTform_L = sc_uvMat_from_uvMap(NNF_L.uvTform.map, uvPixL.ind);
uvTform_L(:,7:8) = uvTform_L(:,7:8)*diag([1/sX, 1/sY]);

% Refinement
refineVec = NNF_H.uvPix.sub - uvPixL.sub*diag([1/sX, 1/sY]);
uvTform_H = sc_trans_tform(uvTform_L, refineVec);

% Clamp
uvValid_H = sc_check_valid_uv(uvTform_H(:,7:8), NNF_H.validPix.mask);

uvInvalidInd = ~uvValid_H;
nInvalidUv_H = sum(uvInvalidInd);
if(nInvalidUv_H)
    randInd = randi(size(NNF_H.validPix.ind, 1), nInvalidUv_H, 1);
    uvTform_H(uvInvalidInd,7:8) = NNF_H.validPix.sub(randInd, :);
end

% Update uvTform.map
NNF_H.uvTform.data = uvTform_H;
I = reshape(eye(3), 1, 1, 9);
NNF_H.uvTform.map = repmat(I, [imgH_H, imgW_H, 1]);
NNF_H.uvTform.map = sc_update_uvMap(NNF_H.uvTform.map, NNF_H.uvTform.data, NNF_H.uvPix.ind);

% =========================================================================
% Initialize uvBias
% =========================================================================
NNF_H.uvBias.data = sc_uvMat_from_uvMap(NNF_L.uvBias.map, uvPixL.ind);
NNF_H.uvBias.map  = zeros(imgH_H, imgW_H, 3, 'single');
NNF_H.uvBias.map  = sc_update_uvMap(NNF_H.uvBias.map, NNF_H.uvBias.data, NNF_H.uvPix.ind);
NNF_H.uvBias.data = reshape(NNF_H.uvBias.data', [1, 3, NNF_H.uvPix.numPix]);

% =========================================================================
% Initialize uvCost
% =========================================================================
NNF_H.uvCost.map  = zeros(imgH_H, imgW_H, 'single');
NNF_H.uvCost.data = sc_uvMat_from_uvMap(NNF_L.uvCost.map, uvPixL.ind);
NNF_H.uvCost.map = sc_update_uvMap(NNF_H.uvCost.map, NNF_H.uvCost.data, NNF_H.uvPix.ind);


% =========================================================================
% Initialize distMap, wPatchR, wPatchSumImg
% =========================================================================

[NNF_H.distMap, ~] = bwdist(~holeMask, 'euclidean');
[NNF_H.wPatchR, NNF_H.wPatchSumImg] = sc_prep_dist_patch(NNF_H.distMap, NNF_H.uvPix.sub, optS);

NNF_H.uvDtBdPixPos = double(NNF_H.distMap(NNF_H.uvPix.ind));

end

% Utility functions
function planeProbAcc = sc_prep_plane_prob_acc(planeProb, uvPixInd)

numPlane = size(planeProb, 3);
numUvPix = size(uvPixInd, 1);

planeProbAcc = zeros(numUvPix, numPlane + 1, 'single');

% Compute the accumulative probability
for i = 1: numPlane
    planeProbCur = planeProb(:,:,i);
    planeProbAcc(:, i+1) = planeProbCur(uvPixInd);
    % Accumulative probability, starting from 0
    if(i ~= 1)
        planeProbAcc(:, i+1) = planeProbAcc(:, i+1) + planeProbAcc(:, i);
    end
end

end

function [validPix, uvPix] = sc_init_level(mask, psize, prad)

% Functionality:
%   get holePixels, uvPixels, validPixels in the level
%
% holePixels: locations of holes pixels
% uvPixels: locations of patches with holes
% validPixels: locations of patches without holes

% =========================================================================
% Get uvPix: center locations of patches with holes
% =========================================================================
uvMask = imdilate(mask, strel('square', double(psize)));
uvMask([1:prad, end-prad+1:end], :) = 0;
uvMask(:, [1:prad, end-prad+1:end]) = 0;

uvPix = getUvPix(uvMask);

% =========================================================================
% Get validPixels: center locations of patches without holes
% =========================================================================
validMap = ~uvMask;
validMap([1:prad,end-prad+1:end], :) = 0;
validMap(:, [1:prad,end-prad+1:end]) = 0;

validPix = getUvPix(validMap);

end

function uvPix = getUvPix(uvMap)

[rUv, cUv] = find(uvMap);
uvPix.sub    = single(cat(2, cUv, rUv));
uvPix.ind    = uint32(sub2ind(size(uvMap), rUv, cUv));
uvPix.mask   = uvMap;
uvPix.numPix = size(uvPix.ind, 1);

end

function uvPixN = getUvPixN(NNF, optS)

uvPixN = cell(4,1);
for i = 1: 4
    uvPixN{i}.sub = bsxfun(@minus, NNF.uvPix.sub, optS.propDir(i,:));
    uvPixN{i}.ind = uint32(sub2ind([NNF.imgH, NNF.imgW], ...
        uvPixN{i}.sub(:,2), uvPixN{i}.sub(:,1)));
    uvPixN{i}.validInd = NNF.uvPix.mask(uvPixN{i}.ind);
end

end

function trgPatchInd = getTargPathInd(NNF, optS)

pixIndMap = single(reshape(1:NNF.imgH*NNF.imgW, NNF.imgH, NNF.imgW));
indTrgPatch   = sc_prep_target_patch(pixIndMap, NNF.uvPix.sub, optS);
trgPatchInd = cat(2, indTrgPatch, ...
    indTrgPatch+NNF.imgH*NNF.imgW, indTrgPatch + 2*NNF.imgH*NNF.imgW);
trgPatchInd = reshape(trgPatchInd, [optS.pNumPix*3, NNF.uvPix.numPix]);

end