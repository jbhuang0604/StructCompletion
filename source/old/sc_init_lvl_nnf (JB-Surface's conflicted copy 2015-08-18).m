function [img, NNF] = sc_init_lvl_nnf(img, NNF, holeMask, modelPlaneCur, modelRegCur, optS)

% SC_INIT_LVL_NNF

% Initialize the nearest neighbor field for the current level

% Input:
%   - NNF, img, holeMask, modelPlaneCur, modelRegCur, iLvl, optS
% Output:
%   - img
%   - NNF

% Prepare distance weight
% [distMap, idMap] = bwdist(~holeMask, 'euclidean');

if(optS.iLvl == optS.numPyrLvl)
    % Initialize the NNF for the coarest level using random sampling
    NNF = sc_init_nnf(holeMask, modelPlaneCur, modelRegCur, optS);
else
    % Initialize the NNF upsampling of NNF from previous level
    NNF = sc_upsample(holeMask, NNF, modelPlaneCur, modelRegCur, optS);
    [wDistPatch, wDistImg] = sc_prep_dist_patch(distMap, NNF.uvPix.sub, iLvl, optS);
    img = sc_voting(img, NNF, NNF.uvPix, holeMask, wDistPatch, wDistImg, optS);
end

% NNF.uvDtBdPixPos = double(distMap(NNF.uvPix.ind));

end


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

% TO-DO: change to column major

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
% 4-connected pixels
%     2
%     |
% 1 - x - 3
%     |
%     4
NNF.uvPixN = cell(4,1);
for i = 1: 4
    NNF.uvPixN{i}.sub = bsxfun(@minus, NNF.uvPix.sub, optS.propDir(i,:));
    NNF.uvPixN{i}.ind = uint32(sub2ind([NNF.imgH, NNF.imgW], ...
        NNF.uvPixN{i}.sub(:,2), NNF.uvPixN{i}.sub(:,1)));
    NNF.uvPixN{i}.validInd = NNF.uvPix.mask(NNF.uvPixN{i}.ind);
end

% =========================================================================
% Initialize indMap
% =========================================================================
NNF.pixIndMap = single(reshape(1:NNF.imgH*NNF.imgW, NNF.imgH, NNF.imgW));
indTrgPatch   = sc_prep_target_patch(NNF.pixIndMap, NNF.uvPix.sub, optS);
NNF.trgPatchInd = cat(2, indTrgPatch, ...
    indTrgPatch+NNF.imgH*NNF.imgW, indTrgPatch + 2*NNF.imgH*NNF.imgW);

% =========================================================================
% Initialize uvPlaneID
% =========================================================================
% Accumulative plane probability
planeProbAccData = sc_prep_plane_prob_acc(modelPlane.planeProb, NNF.uvPix.ind);

% Draw a random sample from U[0,1] for each uv pixel to sample its plane ID
uvPlaneIDData= sc_draw_plane_id(planeProbAccData);
% uvPlaneIDData = modelPlane.numPlane*ones(NNF.uvPix.numPix, 1, 'uint8');
% uvPlaneIDData = ones(NNF.uvPix.numPix, 1, 'single');

% Save data
NNF.uvPlaneID.data = uvPlaneIDData;
NNF.uvPlaneID.map  = zeros(NNF.imgH, NNF.imgW, 'uint8');
NNF.uvPlaneID.map(NNF.uvPix.ind) = NNF.uvPlaneID.data;
NNF.uvPlaneID.planeProbAcc = planeProbAccData;
NNF.uvPlaneID.mLogLikelihood = -log(planeProbAccData);

% =========================================================================
% Initialize uvTform
% =========================================================================
% Initialize unTform.data with random samples
if(NNF.validPix.numPix)
    % Random sampling
    randInd   = randi(NNF.validPix.numPix, NNF.uvPix.numPix, 1);
    uvRandSub = NNF.validPix.sub(randInd, :);
else
    uvRandSub = (NNF.imgW/2)*ones(2, NNF.uvPix.numUvPix);
end
NNF.uvTform.data = sc_src_domain_tform(NNF.uvPlaneID.data, modelPlane, modelReg, uvRandSub, NNF.uvPix.sub, 1);
NNF.uvTform.map = zeros(NNF.imgH, NNF.imgW, 9);
NNF.uvTform.map = sc_update_uvMap(NNF.uvTform.map, NNF.uvTform.data, NNF.uvPix.ind);

% =========================================================================
% Initialize uvBias
% =========================================================================
NNF.uvBias.data = zeros(NNF.uvPix.numPix, 3, 'single');
NNF.uvBias.map  = zeros(NNF.imgH, NNF.imgW, 3, 'single');

% =========================================================================
% Initialize uvCost
% =========================================================================
NNF.uvCost.data = zeros(NNF.uvPix.numPix, 3, 'single');
NNF.uvCost.map  = zeros(NNF.imgH, NNF.imgW, 3, 'single');

% =========================================================================
% Initialize update
% =========================================================================
NNF.update.data = false(NNF.uvPix.numPix, 1);
NNF.update.map  = false(NNF.imgH, NNF.imgW);

% =========================================================================
% Initialize distMap, wPatchR, NNF.wPatchR
% =========================================================================

[NNF.distMap, ~] = bwdist(~holeMask, 'euclidean');
[NNF.wPatchR, NNF.wPatchSumImg] = sc_prep_dist_patch(NNF.distMap, NNF.uvPix.sub, optS);

% %% === Initialize uvPixUpdateSrc ===
% NNF.uvPixUpdateSrc.data = zeros(1, NNF.uvPix.numUvPix);
% NNF.uvPixUpdateSrc.map  = zeros(NNF.imgH, NNF.imgW);

end

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