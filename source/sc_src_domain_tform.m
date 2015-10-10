function uvTformData = sc_src_domain_tform(uvPlaneID, modelPlane, modelReg, srcPos, trgPos, sampleRandReg)

% SC_SRC_DOMAIN_TFORM: The function estimates the source patch domain
% transformation based on the source and target position, and the plane parameters
%
% Input:
%     - uvPlaneID:  numUvPix x 1, the plane label of uv pixels
%     - modelPlane: plane model
%     - modelReg:   regularity model
%     - srcPos:     numUvPix x 2, the center position of source patch
%     - trgPos:     numUvPix x 2, the center position of target patch
%     - sampleRandReg:
%         1: random sampling mode
%         0: regular sampling mode
% Output:
%     - uvTformData: numUvPix x 9

numUvPix    = size(srcPos, 1);

uvTformData = zeros(numUvPix, 9, 'single');
I = eye(3);

for indPlane = 1: modelPlane.numPlane
    % The rectifying transformation for the plane
    rectMat = modelPlane.rectMat{indPlane};
    h7 = rectMat(3,1);
    h8 = rectMat(3,2);
    
    % Retrieve the uv pixels that have the current plane label
    uvPlaneIndCur = uvPlaneID == indPlane;
    numPlanePixCur = sum(uvPlaneIndCur);
    
    if(numPlanePixCur)
        % Target patch center position in the rectified domain
        trgPosCur = trgPos(uvPlaneIndCur, :) - 1;
        trgPosCurR = sc_apply_tform_H(trgPosCur, h7, h8);
        
        % Get dRect
        if(sampleRandReg)    % Random sampling
            % Source patch center position in the rectified domain
            srcPosCur = srcPos(uvPlaneIndCur, :) - 1;
            srcPosCurR = sc_apply_tform_H(srcPosCur, h7, h8);
            
            % Displacement vector from target to source position
            dRect = srcPosCurR - trgPosCurR;
        else                 % Regularity guided sampling
            dRect = zeros(numPlanePixCur, 2, 'single');
            
            numDispVecCur = modelReg.numDispVec(indPlane);
            if(numDispVecCur~=0)
                randInd = randi(numDispVecCur, numPlanePixCur, 1);
                dRect = modelReg.dispVec{indPlane}(randInd, :);
            end
        end
        
        % Compute the transformation that maps from target to source
        % (See Eqn 8 in the paper)
        if(size(dRect,2) ~= 0)
            uvTformCur = zeros(numPlanePixCur, 9, 'single');
            uvTformCur(:,[1,4,7]) = bsxfun(@times, dRect(:,1), [h7, h8, 1]);
            uvTformCur(:,[2,5,8]) = bsxfun(@times, dRect(:,2), [h7, h8, 1]);
            
            dTemp = dRect*[h7;h8]; % dTemp = dx*h7 + dy*h8
            uvTformCur(:,[3,6,9]) = bsxfun(@times, dTemp, -[h7, h8, 1]);
            uvTformCur = bsxfun(@plus, uvTformCur, I(:)');
            
            % Apply the offset to cancel out the dependency of the target position
            % (See Eqn 9 in the paper)
            uvTformData(uvPlaneIndCur, :)   = sc_trans_tform(uvTformCur, trgPosCur);
            uvTformData(uvPlaneIndCur, 7:8) = uvTformData(uvPlaneIndCur, 7:8) + 1;
        end
    end
end


end

function y = sc_apply_tform_H(x, h7, h8)
% Apply homography H with third row [h7, h8, 1] to 2D points x

y = x(:,1)*h7 + x(:,2)*h8 + 1;
y = bsxfun(@rdivide, x(:,1:2), y + eps);

end