function NNFVis= sc_vis_nnf(NNF)

mask = NNF.uvPix.mask;

[imgH, imgW] = size(mask);

bdMask = bwperim(mask);
bdMask = bdMask(:,:,ones(3,1));
mask = mask(:,:,ones(3,1));

%% uvTfomMapVis
NNFVis.uvTfomMapVis = vis_tform_map(NNF.uvTform.map, mask, bdMask);

%% uvCostMapVis

NNFVis.uvCostMapVis = (NNF.uvCost.map - min(NNF.uvCost.map(:)))/max(NNF.uvCost.map(:));

%% uvBiasMapVis

NNFVis.uvBiasMapVis = NNF.uvBias.map + 0.5;
NNFVis.uvBiasMapVis(bdMask) = 1;

%% uvPixUpdateSrc

% NNFVis.uvPixUpdateSrcMap = zeros(imgH, imgW, 3);
% for ch = 1:3
%     NNFVis.uvPixUpdateSrcMap(:,:,ch) = im2double(NNF.uvPixUpdateSrc.map == ch);
% end
% NNFVis.uvGainMapVis = NNF.uvGain.map;
% NNFVis.uvGainMapVis = NNFVis.uvGainMapVis - 0.5;
% NNFVis.uvGainMapVis(bdMask) = 1;

%% uvPlaneIDMapVis
numPlane = size(NNF.uvPlaneID.planeProbAcc, 2) - 1;
NNFVis.uvPlaneIDMapVis = zeros(imgH, imgW, 3);
for i = numPlane:-1:1
    if(i == numPlane)
        fpIDMap = NNF.uvPlaneID.map == i;
        for ch = 1: 3
            NNFVis.uvPlaneIDMapVis(:,:,ch) = NNFVis.uvPlaneIDMapVis(:,:,ch) +  im2double(fpIDMap);
        end
        %         NNFVis.uvPlaneIDMapVis(:,:,1) = fpIDMap;
        %         NNFVis.uvPlaneIDMapVis(:,:,2) = fpIDMap;
        %         NNFVis.uvPlaneIDMapVis(:,:,3) = fpIDMap;
    else
        NNFVis.uvPlaneIDMapVis(:,:,i) = NNFVis.uvPlaneIDMapVis(:,:,i) + ...
            im2double(NNF.uvPlaneID.map == i);
    end
end

% NNFVis.uvPlaneIDMapVis()

NNFVis.uvPlaneIDMapVis = im2double(NNFVis.uvPlaneIDMapVis);
% NNFVis.uvPlaneIDMapVis = NNF.uvPlaneID.map/numPlane;

end

function NNFMapVis = vis_tform_map(NNFMap, mask, bdMask)

[imgH, imgW, ch] = size(mask);

% Initialize NNFMapVis
NNFMapVis = zeros(imgH, imgW, 3, 'single');
[X, Y] = meshgrid(1:imgW, 1:imgH);
X = X/imgW;     Y = Y/imgH;

NNFMapVis(:,:,2) = 0.5;
NNFMapVis(:,:,1) = X;
NNFMapVis(:,:,3) = Y;

% Prepare the visualization of the NNF
NNFMapCurVis = zeros(imgH, imgW, 3, 'single');
NNFMapCurVis(:,:,2) = 0.5;
NNFMapCurVis(:,:,1) = NNFMap(:,:,7)/imgW;
NNFMapCurVis(:,:,3) = NNFMap(:,:,8)/imgH;

NNFMapVis = NNFMapVis.*(1-mask) + NNFMapCurVis.*mask;
NNFMapVis(bdMask) = 1;

end