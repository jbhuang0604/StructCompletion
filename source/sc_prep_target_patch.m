function trgPatch = sc_prep_target_patch(img, uvPixSub, optS)
% SC_PREP_TARGET_PATCH: Prepare target patches with centers uvPixSub
%
% Input:
%   - img:       input image
%   - uvPixSub:  target patch position
%   - optS:      parameter
% Output:
%   - trgPatch: [pNumPix] x [3] x [numUvPix]

[imgH, imgW, nCh] = size(img);

% Updated target patch extraction
numUvPix = size(uvPixSub, 1);

% Prepare target patch position
uvPixSub    = reshape(uvPixSub', [1, 2, numUvPix]);
trgPatchPos = bsxfun(@plus, optS.refPatchPos(:,1:2), uvPixSub);

% Sample target patch
trgPatchInd = zeros(optS.pNumPix, nCh, numUvPix, 'single');
for i = 1: nCh
    trgPatchInd(:, i, :) = sub2ind([imgH, imgW, nCh], ...
        trgPatchPos(:,2,:), trgPatchPos(:,1,:), ...
        i*ones(optS.pNumPix, 1, numUvPix, 'single'));
end

% Get target patch via indexing
trgPatch = img(trgPatchInd);

end