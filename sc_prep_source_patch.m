function srcPatch = sc_prep_source_patch(img, uvTform, optS)

% SC_PREP_SOURCE_PATCH
%
% Prepare source patches according to uvTform
%
% Input:
%   - img
%   - uvTform 
%   - optS
% Output:
%   - srcPatch

numUvPix =  size(uvTform, 2);

% Prepare source patch sampling position
srcPatchPos = zeros(optS.pNumPix, 3, numUvPix);

for i = 1 : optS.pNumPix
    dx = optS.refPatchPos(1,i);
    dy = optS.refPatchPos(2,i);
    
    srcPatchPos(i, 1, :) = uvTform(1,:)*dx + uvTform(4,:)*dy + uvTform(7,:);
    srcPatchPos(i, 2, :) = uvTform(2,:)*dx + uvTform(5,:)*dy + uvTform(8,:);
    srcPatchPos(i, 3, :) = uvTform(3,:)*dx + uvTform(6,:)*dy + uvTform(9,:);
end
srcPatchPos(:, 1:2,:) = bsxfun(@rdivide, srcPatchPos(:, 1:2,:), srcPatchPos(:, 3,:));

% Avoid sample out of boundary positions
srcPatchPos(:,1,:) = sc_clamp(srcPatchPos(:,1,:), 1, size(img,2));
srcPatchPos(:,2,:) = sc_clamp(srcPatchPos(:,2,:), 1, size(img,1));

% Sampling source patch
srcPatch = mirt2D_mexinterp(img, srcPatchPos(:, 1, :), srcPatchPos(:, 2, :));
srcPatch = permute(srcPatch, [1,3,2]);
end