function img = sc_voting_update(wDistPatch, img, NNF, uvPix, holeMask, wDistImg, optS)

% SC_VOTING_UPDATE
% 
% Fast voting update, instead of computing the weighted average from all
% source patches, here we only update patches that have been updated
%
[imgH, imgW, nCh] = size(img);

numUvPix = sum(NNF.update.data);
if(numUvPix~=0)
    
    % Prepare source patch
    uvTform = NNF.uvTform.data(:, NNF.update.data);
    srcPatch = sc_prep_source_patch(img, uvTform, optS);
    
    % Update NNF data
    if(optS.useBiasCorrection)
        biasPatch = NNF.uvBias.data(:, NNF.update.data);
    end
    wDistPatch = wDistPatch(:, NNF.update.data);
    trgPatchInd = NNF.trgPatchInd(:, NNF.update.data);
    
    % Apply bias and gain correction
    if(optS.useBiasCorrection)
        biasPatch = reshape(biasPatch, 1, nCh, numUvPix);
        srcPatch = bsxfun(@plus, srcPatch, biasPatch);
    end
    % Weight by distance to the closest known pixel
    wDistPatchC = reshape(wDistPatch(optS.pMidPix, :), 1, 1, numUvPix);
    srcPatch = bsxfun(@times, srcPatch, wDistPatchC);
        
    % Compute weighted average from source patches
    srcPatch = reshape(srcPatch, optS.pNumPix*nCh, numUvPix);
    imgAcc = zeros(imgH, imgW, 3, 'single');
    weightAcc = optS.voteUpdateW*ones(imgH, imgW, 'single');
    
    for i = 1:numUvPix
        imgAcc(trgPatchInd(:,i)) = imgAcc(trgPatchInd(:,i)) + srcPatch(:,i);
        weightAcc(trgPatchInd(1:optS.pNumPix,i)) = weightAcc(trgPatchInd(1:optS.pNumPix,i)) + wDistPatch(optS.pMidPix, i);
    end
    imgAcc = imgAcc + optS.voteUpdateW*img;
    imgAcc = imgAcc./weightAcc(:,:,ones(1,1,nCh));
    
    % Merge with known region
    holeMask = holeMask(:,:,ones(1,1,nCh));
    img(holeMask) = imgAcc(holeMask);
    
end

end