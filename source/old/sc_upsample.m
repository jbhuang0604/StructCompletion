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

%% === Initialize uvPix and validPix ===
[NNF_H.imgH, NNF_H.imgW] = size(holeMask);
[NNF_H.validPix, NNF_H.uvPix] = sc_init_level(holeMask, optS.pSize, optS.pRad);

%% === Initialize uvPixN ===
NNF_H.uvPixN = cell(4,1);
for i = 1: 4
    NNF_H.uvPixN{i}.sub(1,:) = NNF_H.uvPix.sub(1,:) - optS.propDir(1,i);
    NNF_H.uvPixN{i}.sub(2,:) = NNF_H.uvPix.sub(2,:) - optS.propDir(2,i);
    NNF_H.uvPixN{i}.ind = sub2ind([NNF_H.imgH, NNF_H.imgW], NNF_H.uvPixN{i}.sub(2,:), NNF_H.uvPixN{i}.sub(1,:));
    NNF_H.uvPixN{i}.validInd = NNF_H.uvPix.mask(NNF_H.uvPixN{i}.ind);
end

%% === Initialize indMap ===
NNF_H.pixIndMap = reshape(1:NNF_H.imgH*NNF_H.imgW, NNF_H.imgH, NNF_H.imgW);
indTrgPatch = sc_prep_target_patch(NNF_H.pixIndMap, NNF_H.uvPix.sub, optS);

NNF_H.trgPatchInd = cat(1, indTrgPatch, ...
    indTrgPatch+NNF_H.imgH*NNF_H.imgW, indTrgPatch + 2*NNF_H.imgH*NNF_H.imgW);

%% === Initialize uvPixL ===

imgH_H = NNF_H.imgH;    imgW_H = NNF_H.imgW;
imgH_L = NNF_L.imgH;    imgW_L = NNF_L.imgW;

sX = imgH_L/imgH_H;     sY = imgW_L/imgW_H;
uvPixL.sub = round(diag([sX, sY])*NNF_H.uvPix.sub);
uvPixL.sub(1,:) = sc_clamp(uvPixL.sub(1,:), optS.pRad+1, imgW_L - optS.pRad);
uvPixL.sub(2,:) = sc_clamp(uvPixL.sub(2,:), optS.pRad+1, imgH_L - optS.pRad);
uvPixL.ind = sub2ind([imgH_L, imgW_L], uvPixL.sub(2,:), uvPixL.sub(1,:));

%% === Initialize uvPlaneID ===

NNF_H.uvPlaneID.data = sc_uvMat_from_uvMap(NNF_L.uvPlaneID.map, uvPixL); 
NNF_H.uvPlaneID.data(NNF_H.uvPlaneID.data==0) = 1;
NNF_H.uvPlaneID.map = zeros(NNF_H.imgH, NNF_H.imgW);
NNF_H.uvPlaneID.map = sc_update_uvMap(NNF_H.uvPlaneID.map, NNF_H.uvPlaneID.data, NNF_H.uvPix, true(1,NNF_H.uvPix.numUvPix));
NNF_H.uvPlaneID.planeProbAcc = sc_prep_plane_prob_acc(modelPlane.planeProb, NNF_H.uvPix.ind);

%% === Initialize uvPixUpdateSrc ===
NNF_H.uvPixUpdateSrc.map = zeros(NNF_H.imgH, NNF_H.imgW);
NNF_H.uvPixUpdateSrc.data = sc_uvMat_from_uvMap(NNF_L.uvPixUpdateSrc.map, uvPixL);
NNF_H.uvPixUpdateSrc.map = sc_update_uvMap(NNF_H.uvPixUpdateSrc.map, NNF_H.uvPixUpdateSrc.data, NNF_H.uvPix, true(1,NNF_H.uvPix.numUvPix));

%% === Update uvTform ===

uvTform_L = sc_uvMat_from_uvMap(NNF_L.uvTform.map, uvPixL);
uvTform_L(7:8,:) = diag([1/sX, 1/sY])*uvTform_L(7:8,:);

% Refinement
refineVec = NNF_H.uvPix.sub - diag([1/sX, 1/sY])*uvPixL.sub;
uvTform_H = sc_trans_tform(uvTform_L, refineVec);

% Clamp
uvTform_H(7,:) = sc_clamp(uvTform_H(7,:), optS.pRad+1, imgW_H - optS.pRad);
uvTform_H(8,:) = sc_clamp(uvTform_H(8,:), optS.pRad+1, imgH_H - optS.pRad);
uvValid_H = sc_check_valid_uv(uvTform_H(7:8,:), NNF_H.validPix.mask);

uvInvalidInd = ~uvValid_H;
nInvalidUv_H = sum(uvInvalidInd);
if(nInvalidUv_H)
    randInd = randi(size(NNF_H.validPix.ind, 2), nInvalidUv_H, 1);
    uvRand = NNF_H.validPix.sub(:, randInd);
    uvTform_H(7,uvInvalidInd) = uvRand(1,:);
    uvTform_H(8,uvInvalidInd) = uvRand(2,:);
end

% Update uvTform.map
NNF_H.uvTform.data = uvTform_H;
I = reshape(eye(3), 1, 1, 9);
NNF_H.uvTform.map = repmat(I, [imgH_H, imgW_H, 1]);
NNF_H.uvTform.map = sc_update_uvMap(NNF_H.uvTform.map, NNF_H.uvTform.data, NNF_H.uvPix, true(1,NNF_H.uvPix.numUvPix));

%% === Initialize uvBias ===
NNF_H.uvBias.map    = im2single(zeros(imgH_H, imgW_H, 3));
NNF_H.uvBias.data = sc_uvMat_from_uvMap(NNF_L.uvBias.map, uvPixL);
NNF_H.uvBias.map = sc_update_uvMap(NNF_H.uvBias.map, NNF_H.uvBias.data, NNF_H.uvPix, true(1,NNF_H.uvPix.numUvPix));

%% === Initialize uvCost ===
NNF_H.uvCost.map  = im2single(zeros(imgH_H, imgW_H));
NNF_H.uvCost.data = sc_uvMat_from_uvMap(NNF_L.uvCost.map, uvPixL);
NNF_H.uvCost.map = sc_update_uvMap(NNF_H.uvCost.map, NNF_H.uvCost.data, NNF_H.uvPix, true(1,NNF_H.uvPix.numUvPix));

end