function img = sc_voting_update(img, NNF, holeMask, optS)
% SC_VOTING_UPDATE: Fast voting update, instead of computing the
% weighted average from all source patches, here we only update patches
% that have been updated

[imgH, imgW, nCh] = size(img);

numUvPix = sum(NNF.update.data);
if(numUvPix~=0)
    % Prepare source patch
    uvTform = NNF.uvTform.data(NNF.update.data, :);
    srcPatch = sc_prep_source_patch(img, uvTform, optS);
    
    % Update NNF data
    if(optS.useBiasCorrection)
        biasPatch = NNF.uvBias.data(:, :, NNF.update.data);
    end
    trgPatchInd = NNF.trgPatchInd(:, NNF.update.data);
    
    % Apply bias and gain correction
    if(optS.useBiasCorrection)
        srcPatch = bsxfun(@plus, srcPatch, biasPatch);
    end
    % Weight by distance to the closest known pixel
    wPatchR  = NNF.wPatchR(:, NNF.update.data);
    wPatchR  = reshape(wPatchR, optS.pNumPix, 1, numUvPix);
    srcPatch = bsxfun(@times, srcPatch, wPatchR);
    
    % Compute weighted average from source patches
    srcPatch = reshape(srcPatch, optS.pNumPix*nCh, numUvPix);
    
    % Initialization
    imgAcc    = zeros(imgH, imgW, 3, 'single');
    weightAcc = optS.voteUpdateW*ones(imgH, imgW, 'single');
    
    for i = 1 : numUvPix
        imgAcc(trgPatchInd(:,i)) = imgAcc(trgPatchInd(:,i)) + srcPatch(:,i);
        weightAcc(trgPatchInd(1:optS.pNumPix,i)) = ...
            weightAcc(trgPatchInd(1:optS.pNumPix,i)) + wPatchR(:,i);
    end
    imgAcc = imgAcc + optS.voteUpdateW*img;
    imgAcc = bsxfun(@rdivide, imgAcc, weightAcc);
    
    % Merge with known region
    holeMask = holeMask(:,:,ones(1,1,nCh));
    img(holeMask) = imgAcc(holeMask);
end

end