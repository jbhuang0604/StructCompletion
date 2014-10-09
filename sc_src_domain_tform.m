function uvTformData = sc_src_domain_tform(uvPlaneID, modelPlane, modelReg, srcPos, trgPos, sampleRandReg)

% SC_SRC_DOMAIN_TFORM:
%
% The function estimates the source patch domain transformation based on
% the source and target position, and the plane parameters
%
% Input:
%     - uvPlaneID: 1 x numUvPix, the plane label of uv pixels
%     - modelPlane: plane model
%     - modelReg: regularity model
%     - srcPos: 2 x numUvPix, the center position of source patch
%     - trgPos: 2 x numUvPix, the center position of target patch
%     - sampleRandReg:
%         1: random sampling mode
%         0: regular sampling mode
% Output:
%     - uvTformData: 9 x numUvPix,

%
 
numUvPix = size(srcPos, 2);

uvTformData = zeros(9, numUvPix);
I = eye(3);


for indPlane = 1: modelPlane.numPlane
    % The rectifying transformation for the plane
    rectMat = modelPlane.rectMat{indPlane};
    h7 = rectMat(3,1);
    h8 = rectMat(3,2);
    
    % Retrieve the uv pixels that have the current plane label
    uvPlaneIndCur = uvPlaneID == indPlane;
    numPlanePixCur = sum(uvPlaneIndCur);
    
    % Target patch center position in the rectified domain
    trgPosCur = trgPos(:, uvPlaneIndCur) - 1;
    trgPosCurRect = rectMat*cat(1, trgPosCur, ones(1, numPlanePixCur));
    trgPosCurRect = bsxfun(@rdivide, trgPosCurRect, trgPosCurRect(3,:));
    
    if(sampleRandReg)
        % Random sampling
        
        % Source patch center position in the rectified domain
        srcPosCur = srcPos(:, uvPlaneIndCur) - 1;
        srcPosCurRect = rectMat*cat(1, srcPosCur, ones(1, numPlanePixCur));
        
        validInd = srcPosCurRect(3,:) ~= 0;
        srcPosCurRect(1:2, validInd) = bsxfun(@rdivide, srcPosCurRect(1:2,validInd), srcPosCurRect(3,validInd));
        % Displacement vector from target to source position
        dRect = srcPosCurRect(1:2,:) - trgPosCurRect(1:2,:);        
        dRect(:, ~validInd) = 0;
    else
        % Regularity guided sampling
        dRect = zeros(2, numPlanePixCur);
        
        numDispVecCur = modelReg.numDispVec(indPlane);
        if(numDispVecCur~=0)
            randInd = randi(numDispVecCur, 1, numPlanePixCur);
            dRect = modelReg.dispVec{indPlane}(:, randInd);
        end
    end
    % Compute the transformation that maps from target to source
    % (See Eqn 8 in the paper)
    if(size(dRect,2)~=0)
    uvTformCur = zeros(9, numPlanePixCur, 'single');
    uvTformCur([1,4,7],:) = bsxfun(@times, dRect(1,:), [h7;h8;1]);
    uvTformCur([2,5,8],:) = bsxfun(@times, dRect(2,:), [h7;h8;1]);
    
    dTemp = [h7, h8]*dRect; % dTemp = dx*h7 + dy*h8
    uvTformCur([3,6,9],:) = bsxfun(@times, dTemp, -[h7; h8; 1]);
    uvTformCur = bsxfun(@plus, uvTformCur, I(:));
    
    % Apply the offset to cancel out the dependency of the target position
    % (See Eqn 9 in the paper)
    % HERE there is a NAN occur
    uvTformData(:, uvPlaneIndCur) = sc_trans_tform(uvTformCur, trgPosCur);
    uvTformData(7:8, uvPlaneIndCur) = uvTformData(7:8, uvPlaneIndCur) + 1;
    end
end


end