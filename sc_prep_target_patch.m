function trgPatch = sc_prep_target_patch(img, uvPixSub, optS)

% SC_PREP_TARGET_PATCH
%
% Prepare target patches with centers uvPixSub
%
% Input:
%   - img
%   - uvPixSub 
%   - optS
% Output:
%   - trgPatch

[imgH, imgW, nCh] = size(img);

% Updated target patch extraction
numUvPix = size(uvPixSub, 2);

% Prepare target patch position
refPatchPos = reshape(optS.refPatchPos(1:2,:)', optS.pNumPix, 2, 1);
trgPatchCenter = reshape(uvPixSub, 1, 2, numUvPix);
trgPatchPos = bsxfun(@plus, refPatchPos, trgPatchCenter);

% Sample target patch
img = im2double(img);
trgPatchPos = im2double(trgPatchPos);
trgPatch = mirt2D_mexinterp(img, trgPatchPos(:, 1, :), trgPatchPos(:, 2, :));
if(size(trgPatch, 3) > 1)
    trgPatch = permute(trgPatch, [1, 3, 2]);
end

end