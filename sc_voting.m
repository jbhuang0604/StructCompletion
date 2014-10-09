function img = sc_voting(img, NNF, uvPix, holeMask, wDistPatch, wDistImg, optS)

% SC_VOTING: update the image using NNF
%
% Input: 
%   - wDistPatch
%   - img
%   - NNF
%   - uvPix
%   - holeMask
%   - wDistImg, optS
% Output:
%   - img

[imgH, imgW, nCh] = size(img);

% Prepare source patch
numUvPix = size(uvPix.ind, 2);
uvValid.ind = true(1, numUvPix);    uvValid.pos = 1:numUvPix;
srcPatch = sc_prep_source_patch(img, NNF.uvTform.data, optS);

% Apply bias correction
if(optS.useBiasCorrection)
    biasPatch = reshape(NNF.uvBias.data, 1, nCh, numUvPix);
    srcPatch = bsxfun(@plus, srcPatch, biasPatch);
end
% Weight by distance to the hole border
wDistPatchC = reshape(wDistPatch(optS.pMidPix, :), 1, 1, numUvPix);
srcPatch = bsxfun(@times, srcPatch, wDistPatchC);

% Compute weighted average from source patches
srcPatch = reshape(srcPatch, optS.pNumPix*nCh, numUvPix);
imgAcc = zeros(imgH, imgW, 3, 'single');
for i = 1:numUvPix
    imgAcc(NNF.trgPatchInd(:,i)) = imgAcc(NNF.trgPatchInd(:,i)) + srcPatch(:,i);
end
imgAcc = imgAcc./wDistImg(:,:,ones(1,1,nCh));

% Merge with known regions
holeMask = holeMask(:,:,ones(1,1,nCh));
img(holeMask) = imgAcc(holeMask);

end

