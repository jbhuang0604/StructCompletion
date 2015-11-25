% function NNF = sc_init_nnf(holeMask, modelPlane, modelReg, optS)
% 
% % SC_INIT_NNF: Initialize the nearest neighbor field
% % Input:
% %     - holeMask: the binary mask indicate the hole pixels
% %     - modelPlane: model of detected plane structures
% %     - modelReg: model of detected reguarity structures
% %     - optS: parameters for synthesis process
% % Output:
% %   Image height and width
% %     - NNF.imgH: image height
% %     - NNF.imgW: image width
% %   Precomputing pixel positions for PatchMatch algorithm
% %     - NNF.uvPix: the patch positions containing hole pixels
% %     - NNF.uvPixN: the neighboring patch positions of 4 connected pixels
% %     - NNF.validPix: the patch positions containing known pixels
% %     - NNF.trgPatchInd: the indice to retrieve target patches
% %   Initialize NNF components:
% %     - NNF.uvTform: the transformation specifying the source patch
% %         - NNF.uvTform.data:   9 x numUvPix
% %         - NNF.uvTform.map:    H x W x 9
% %     - NNF.uvPlaneID: the plane assignment for unknown region
% %         - NNF.uvPlaneID.data: 1 x numUvPix
% %         - NNF.uvPlaneID.map:  H x W
% %     - NNF.uvBias: the bias between source and target patches
% %         - NNF.uvBias.data:    3 x numUvPix
% %         - NNF.uvBias.map:     H x W x 3
% %     - NNF.uvGain: the gain between source and target patches
% %         - NNF.uvGain.data:    3 x numUvPix
% %         - NNF.uvGain.map:     H x W x 3
% %     - NNF.uvCost: the matching cost between source and target patches
% %         - NNF.uvCost.data:    3 x numUvPix
% %         - NNF.uvCost.map:     H x W x 3
% 
% %% === Initialize uvPix and validPix ===
% [NNF.imgH, NNF.imgW] = size(holeMask);
% [NNF.validPix, NNF.uvPix] = sc_init_level(holeMask, optS.pSize, optS.pRad);
% 
% %% === Initialize uvPixN ===
% NNF.uvPixN = cell(4,1);
% for i = 1: 4
%     NNF.uvPixN{i}.sub(1,:) = NNF.uvPix.sub(1,:) - optS.propDir(1,i);
%     NNF.uvPixN{i}.sub(2,:) = NNF.uvPix.sub(2,:) - optS.propDir(2,i);
%     NNF.uvPixN{i}.ind = sub2ind([NNF.imgH, NNF.imgW], NNF.uvPixN{i}.sub(2,:), NNF.uvPixN{i}.sub(1,:));
%     NNF.uvPixN{i}.validInd = NNF.uvPix.mask(NNF.uvPixN{i}.ind);
% end
% 
% %% === Initialize indMap ===
% NNF.pixIndMap = reshape(1:NNF.imgH*NNF.imgW, NNF.imgH, NNF.imgW);
% indTrgPatch = sc_prep_target_patch(NNF.pixIndMap, NNF.uvPix.sub, optS);
% 
% NNF.trgPatchInd = cat(1, indTrgPatch, ...
%     indTrgPatch+NNF.imgH*NNF.imgW, indTrgPatch + 2*NNF.imgH*NNF.imgW);
% 
% %% === Initialize uvPlaneID ===
% % Accumulative plane probability
% planeProbAccData = sc_prep_plane_prob_acc(modelPlane.planeProb, NNF.uvPix.ind);
% 
% % Draw a random sample from U[0,1] for each uv pixel to sample its plane ID
% % uvPlaneIDData= sc_draw_plane_id(planeProbAccData);
% uvPlaneIDData = modelPlane.numPlane*ones(1, size(NNF.uvPix.ind, 2));
% 
% % Save data
% NNF.uvPlaneID.data = uvPlaneIDData;
% NNF.uvPlaneID.map  = zeros(NNF.imgH, NNF.imgW);
% NNF.uvPlaneID.map(NNF.uvPix.ind) = NNF.uvPlaneID.data;
% NNF.uvPlaneID.planeProbAcc = planeProbAccData;
% NNF.uvPlaneID.mLogLikelihood = -log(planeProbAccData);
% 
% %% === Initialize uvTform ===
% % Initialize unTform.data with random samples
% if(NNF.validPix.numValidPix)
% % Random sampling
% randInd = randi(NNF.validPix.numValidPix, NNF.uvPix.numUvPix, 1);
% uvRandSub = NNF.validPix.sub(:, randInd);
% else
%     uvRandSub = (NNF.imgW/2)*ones(2, NNF.uvPix.numUvPix);
% end
% NNF.uvTform.data = sc_src_domain_tform(NNF.uvPlaneID.data, modelPlane, modelReg, uvRandSub, NNF.uvPix.sub, 1);
% NNF.uvTform.map = zeros(NNF.imgH, NNF.imgW, 9);
% NNF.uvTform.map = sc_update_uvMap(NNF.uvTform.map, NNF.uvTform.data, NNF.uvPix, true(1,NNF.uvPix.numUvPix));
% 
% %% === Initialize uvBias ===
% NNF.uvBias.data = zeros(3, NNF.uvPix.numUvPix);
% NNF.uvBias.map    = im2single(zeros(NNF.imgH, NNF.imgW, 3));
% 
% %% === Initialize uvCost ===
% NNF.uvCost.data = zeros(3, NNF.uvPix.numUvPix);
% NNF.uvCost.map  = im2single(zeros(NNF.imgH, NNF.imgW));
% 
% %% === Initialize update ===
% NNF.update.data = false(1, NNF.uvPix.numUvPix);
% NNF.update.map  = false(NNF.imgH, NNF.imgW);
% 
% %% === Initialize uvPixUpdateSrc ===
% NNF.uvPixUpdateSrc.data = zeros(1, NNF.uvPix.numUvPix);
% NNF.uvPixUpdateSrc.map  = zeros(NNF.imgH, NNF.imgW);
% 
% end