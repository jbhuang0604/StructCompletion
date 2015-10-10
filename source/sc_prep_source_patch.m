function srcPatch = sc_prep_source_patch(img, uvTform, optS)

% SC_PREP_SOURCE_PATCH
%
% Prepare source patches according to uvTform
%
% Input:
%   - img:       input image
%   - uvPixSub:  target patch position
%   - optS:      parameter
% Output:
%   - srcPatch: [pNumPix] x [3] x [numUvPix]

numUvPix =  size(uvTform, 1);

% Prepare source patch sampling position
% srcPatchPos = zeros(optS.pNumPix, 3, numUvPix);

% Get srcPatchPos
c1 = reshape(uvTform(:,1:3)', 1, 3, numUvPix);
c2 = reshape(uvTform(:,4:6)', 1, 3, numUvPix);
c3 = reshape(uvTform(:,7:9)', 1, 3, numUvPix);

% Get the source patch pixel positions
srcPatchPos = bsxfun(@times, optS.refPatchPos(:,1), c1) + ...
              bsxfun(@times, optS.refPatchPos(:,2), c2);
srcPatchPos = bsxfun(@plus, srcPatchPos, c3);

% Convert back to Eucledian coordinate
srcPatchPos = bsxfun(@rdivide, srcPatchPos, srcPatchPos(:, 3, :));

% Grab the color values of source patch using bilinear interpolation
srcPatch = vgg_interp2(img, srcPatchPos(:,1,:), srcPatchPos(:,2,:), 'linear', 0);

% Convert to the format [pNumPix] x [3] x [numUvPix]
srcPatch = permute(srcPatch, [1,3,2]);

% for i = 1 : optS.pNumPix
%     dx = optS.refPatchPos(1,i);
%     dy = optS.refPatchPos(2,i);
%     
%     srcPatchPos(i, 1, :) = uvTform(1,:)*dx + uvTform(4,:)*dy + uvTform(7,:);
%     srcPatchPos(i, 2, :) = uvTform(2,:)*dx + uvTform(5,:)*dy + uvTform(8,:);
%     srcPatchPos(i, 3, :) = uvTform(3,:)*dx + uvTform(6,:)*dy + uvTform(9,:);
% end
% srcPatchPos(:, 1:2,:) = bsxfun(@rdivide, srcPatchPos(:, 1:2,:), srcPatchPos(:, 3,:));

% % Avoid sample out of boundary positions
% % srcPatchPos(:,1,:) = sc_clamp(srcPatchPos(:,1,:), 1, size(img,2));
% % srcPatchPos(:,2,:) = sc_clamp(srcPatchPos(:,2,:), 1, size(img,1));
% 
% % Sampling source patch
% srcPatch = mirt2D_mexinterp(img, srcPatchPos(:, 1, :), srcPatchPos(:, 2, :));


% srcPatch = permute(srcPatch, [1,3,2]);
end