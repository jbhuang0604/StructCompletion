function img = sc_voting(img, NNF, holeMask, optS)

% SC_VOTING: update the image using NNF
%
% Input:
%   - img:      input image at the current level
%   - NNF       nearest neighbor field
%   - holeMask  hole mask
%   - optS      parameters
% Output:
%   - img

numUvPix = NNF.uvPix.numPix;
[imgH, imgW, nCh] = size(img);

% Prepare source patch
srcPatch = sc_prep_source_patch(img, NNF.uvTform.data, optS);

% Update NNF data
if(optS.useBiasCorrection)
    biasPatch = NNF.uvBias.data;
    srcPatch = bsxfun(@plus, srcPatch, biasPatch);
end

% Target patch index
trgPatchInd = NNF.trgPatchInd;

% Weight by distance to the closest known pixel
wPatchR  = reshape(NNF.wPatchR, optS.pNumPix, 1, numUvPix);
srcPatch = bsxfun(@times, srcPatch, wPatchR);

% Compute weighted average from source patches
srcPatch = reshape(srcPatch, optS.pNumPix*nCh, numUvPix);

% Initialization
imgAcc    = zeros(imgH, imgW, nCh, 'single');
for i = 1 : numUvPix
    imgAcc(trgPatchInd(:,i)) = imgAcc(trgPatchInd(:,i)) + srcPatch(:,i);
end
imgAcc = bsxfun(@rdivide, imgAcc, NNF.wPatchSumImg);

% Merge with known region
% uvMask = NNF.uvPix.mask(:,:,ones(nCh,1));
holeMask = holeMask(:,:,ones(nCh,1));
img(holeMask) = imgAcc(holeMask);
img = sc_clamp(img, 0, 1);

end

