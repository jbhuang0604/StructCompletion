function modelReg = sc_extract_regularity(img, mask, modelPlane, optA)

% SC_EXTRACT_REGULARITY
%
% Extract regularity model
%
% Output: modelReg
% The regularity model stores the set of dominant displacement vectors in the
% rectified domain for each plane.
%
% modelReg.plane{indPlane}.dispVec
% modelReg.plane{indPlane}.costMapReg


modelReg = [];

%%
img = im2single(img);
imgGray = rgb2gray(img);
[imgH, imgW] = size(imgGray);

%% Matching SIFT features
% === Load image ===
% img(cat(3, mask, mask, mask)) = 1;

[frames, descriptor] = vl_covdet(imgGray, 'EstimateAffineShape', true, 'PeakThreshold', optA.PeakThreshold);

% Select only first 1500 features
if(size(frames, 2) > optA.maxNumLocalFeat)
    fInd = randperm(size(frames, 2));
    fInd = fInd(1:optA.maxNumLocalFeat);
    frames = frames(:, fInd);
    descriptor = descriptor(:, fInd);
end

% === Feature matching ===
kdtree = vl_kdtreebuild(descriptor);
[index, distance] = vl_kdtreequery(kdtree, descriptor, descriptor, 'NumNeighbors', optA.numQueryNN);

% === Save feature matching data ===

featMatchData.index = index;
featMatchData.distance = distance;
featMatchData.K = optA.numQueryNN;
featMatchData.frames = frames;
featMatchData.numFeat = size(frames, 2);

%% Get the displacement vectors in each plane

modelReg.plane = cell(modelPlane.numPlane, 1);

% Load feature matching data
for indPlane = 1 : modelPlane.numPlane
    
    % === compute the rectifying homography ===
    H = eye(3);
    H(3,:) = modelPlane.plane{indPlane}.vLine;
    
    % === get rectified feature positions ===
    framesRect = cat(1, featMatchData.frames(1:2,:), ones(1, featMatchData.numFeat));
    framesRect = H*framesRect;
    framesRect = framesRect./repmat(framesRect(3,:), 3, 1);
    
    % === use 2-NN ===
    nn = featMatchData.index(1:featMatchData.K, :);
    
    targetPosRect = cat(2, framesRect(1:2, nn(1,:)), framesRect(1:2, nn(1,:)), framesRect(1:2, nn(2,:)));
    sourcePosRect = cat(2, framesRect(1:2, nn(2,:)), framesRect(1:2, nn(3,:)), framesRect(1:2, nn(3,:)));
    
    targetPos = cat(2, (featMatchData.frames(1:2, nn(1,:))), ...
        (featMatchData.frames(1:2, nn(1,:))), ...
        (featMatchData.frames(1:2, nn(2,:))));
    sourcePos = cat(2, (featMatchData.frames(1:2, nn(2,:))), ...
        (featMatchData.frames(1:2, nn(3,:))), ...
        (featMatchData.frames(1:2, nn(3,:))));
    targetPos = round(targetPos);
    sourcePos = round(sourcePos);
    
    % === compute feature matching weights ===
    indTargetPos = sub2ind([imgH, imgW], targetPos(2,:), targetPos(1,:));
    indSourcePos = sub2ind([imgH, imgW], sourcePos(2,:), sourcePos(1,:));
    
    planeProb = modelPlane.postProb(:,:,indPlane);
    weightVec = planeProb(indTargetPos).*planeProb(indSourcePos);
    
    % === get the matched features on the plane ===
    validMatchInd = weightVec >= optA.prodProbThres;
    sourcePosRect = sourcePosRect(:, validMatchInd);
    targetPosRect = targetPosRect(:, validMatchInd);
    
    % === compute the displacement vectors in the rectified space
    dispVecRect = sourcePosRect - targetPosRect;
    distDispVecRect = sum(dispVecRect.^2, 1);
    validDispVecInd = distDispVecRect > optA.minDistDispVec;
    dispVecRect = dispVecRect(:, validDispVecInd);
    
    dispVecRect = cat(2, dispVecRect, -dispVecRect);
    dispVecRect = cat(2, dispVecRect, 2*dispVecRect);
    validDispVecInd = (abs(dispVecRect(1,:)) < 1000 & abs(dispVecRect(2,:)) < 1000);
    dispVecRect = dispVecRect(:, validDispVecInd);
    
    [clustCent, point2cluster, clustMembsCell] = MeanShiftCluster(dispVecRect, optA.msBandwidth);
    
    % === Filtering out weak cluster ===
    numClust = length(clustMembsCell);
    
    numMemberInCluster = zeros(numClust, 1);
    for i = 1: numClust
        numMemberInCluster(i) = length(clustMembsCell{i});
    end
    validClusterInd = numMemberInCluster >= optA.msMinNumClusterMember;
    clustCent = clustCent(:, validClusterInd);
    
    % === 
    modelReg.plane{indPlane}.dispVec = clustCent;
    modelReg.plane{indPlane}.numDispVec = size(clustCent, 2);

end