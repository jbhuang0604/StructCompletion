function modelPlane = sc_detect_plane_from_vp(vpData, img, mask, optA)

% SC_DETECT_PLANE_FROM_VP: simple plane detection algorithm
% Input:
%     - vpData: vanishing point data
%     - img: input image
%     - mask: hole mask
% Output:
%     - modelPlane

%%

modelPlane = [];

% === Setting up ===
[imgH, imgW, ch] = size(img);
HfilterX = fspecial('gaussian', [1, optA.filterSize], optA.filterSigma);
HfilterY = HfilterX';
% fspecial('gaussian', optA.filterSize, optA.filterSigma);

img = im2double(img);

% === Supporting lines spatial support estimation ===
shapeInserter = vision.ShapeInserter('Shape', 'Lines','BorderColor', 'White');

for i = 1: vpData.numVP
    % The support lines
    imgLines = zeros(imgH, imgW);
    imgLines = step(shapeInserter, imgLines, int16(round(vpData.vp{i}.lines(:,1:4))));
    % Spatial density estimation via blurring
    imgLinesPosMap = imgLines;
    for k = 1:optA.numFilterIter
        imgLinesPosMap = imfilter(imgLinesPosMap, HfilterX, 'conv', 'replicate');
    end
    for k = 1:optA.numFilterIter
        imgLinesPosMap = imfilter(imgLinesPosMap, HfilterY, 'conv', 'replicate');
    end
    
    % Save results
    modelPlane.vp{i}.imgLines = imgLines;
    modelPlane.vp{i}.imgLinesPosMap = imgLinesPosMap;
end


% === Estimate plane support and plane parameters ===
numPlane = (vpData.numVP)*(vpData.numVP-1)/2;
% Initialize plane data
modelPlane.plane = cell(numPlane, 1);

indPlane = 1;
% A pair of vanishing points forms a plane hypothesis
for i = 1: vpData.numVP - 1
    for j = i+1: vpData.numVP
        % Compute the vanishing line
        modelPlane.plane{indPlane}.vLine = vLineFromTwoVP(vpData.vp{i}.pos, vpData.vp{j}.pos);
        % Element-wise product of two support line density
        modelPlane.plane{indPlane}.imgPlaneProb = modelPlane.vp{i}.imgLinesPosMap.*modelPlane.vp{j}.imgLinesPosMap; % Product of two probability maps
        
        %
        modelPlane.plane{indPlane}.imgPlaneProb(mask) = 0;
        modelPlane.plane{indPlane}.score = sum(modelPlane.plane{indPlane}.imgPlaneProb(:));
        
        modelPlane.plane{indPlane}.imgPlaneProb(mask) = 1e-10;
        modelPlane.plane{indPlane}.sourceVP = [i, j];
        
        indPlane = indPlane + 1;
    end
end


% === Compute rectified rotation parameters ===

for i = 1: numPlane
    for vpInd = 1: 2
        
        linesCurr = vpData.vp{modelPlane.plane{i}.sourceVP(vpInd)}.lines;
        invalidLineInd = linesCurr(:,5) == 0;
        linesCurr = linesCurr(~invalidLineInd,:);
        numLines = size(linesCurr, 1);
        
        vLineCurr = modelPlane.plane{i}.vLine;
        
        % Rectified homography
        H = eye(3);
        H(3,:) = vLineCurr;
        
        linesStart = cat(2, linesCurr(:,1), linesCurr(:,2),ones(numLines, 1))';
        linesEnd = cat(2, linesCurr(:,3), linesCurr(:,4),ones(numLines, 1))';
        
        linesStartRect = H*linesStart;
        linesStartRect = linesStartRect./repmat(linesStartRect(3,:), 3, 1);
        
        linesEndRect = H*linesEnd;
        linesEndRect = linesEndRect./repmat(linesEndRect(3,:), 3, 1);
        
        linesVec = linesStartRect(1:2, :) - linesEndRect(1:2, :);
        linesSign = linesEndRect(2,:) > linesStartRect(2, :);
        linesSign = 2*linesSign - 1;
        linesLength = sqrt(sum(linesVec.^2, 1));
        linesCos = linesSign.*linesVec(1,:)./linesLength; % repmat(linesLength, 2, 1);
        
        theta = acos(linesCos);
        
        % Estimate average theta so that all the supporting lines aligned
        % with the x-axis
        thetaAvg = mean(theta, 2);
        for iter = 1: 5
            thetaDiff = theta - thetaAvg;
            indLargeTheta = thetaDiff > pi/2;
            theta(indLargeTheta) = pi - theta(indLargeTheta);
            
            indSmallTheta = thetaDiff < -pi/2;
            theta(indSmallTheta) = pi + theta(indSmallTheta);
            thetaAvg = mean(theta, 2);
        end
        
        thetaEst = thetaAvg;
        
        modelPlane.plane{i}.rotPar(vpInd) = thetaEst;
    end
end


% === Add a fronto-parallel plane ===

modelPlane.plane{indPlane}.vLine = [0 0 1];
modelPlane.plane{indPlane}.imgPlaneProb = 1e-5*ones(imgH, imgW);
modelPlane.plane{indPlane}.imgPlaneProb(mask) = 0;
modelPlane.plane{indPlane}.score = sum(modelPlane.plane{indPlane}.imgPlaneProb(:));
modelPlane.plane{indPlane}.imgPlaneProb(mask) = 1e-10;
modelPlane.plane{indPlane}.rotPar(1) = 0;
modelPlane.plane{indPlane}.rotPar(2) = 0;

numPlane = numPlane + 1;

modelPlane.numPlane = numPlane;

% === Compute posterior probability ===

planeProb = zeros(imgH, imgW, numPlane);
for i = 1 : numPlane
    planeProb(:,:,i) = modelPlane.plane{i}.imgPlaneProb;
end
planeProbSum = sum(planeProb, 3);
planeProb = (planeProb)./repmat(planeProbSum, [1, 1, numPlane]);

planeProbSum = 1 + numPlane*optA.probConst;
% planeProb = (planeProb + optA.probConst)./repmat(planeProbSum, [1, 1, numPlane]);
planeProb = (planeProb + optA.probConst)/planeProbSum;

modelPlane.postProbHole = planeProb;

% === Propagate posterior probability into the hole region ===

% Get border pixels
borderImg = im2double(mask);
borderImg = imfilter(borderImg, fspecial('gaussian'));
borderImg = (borderImg~=0) & (borderImg~=1);
borderImg = borderImg & ~mask;

holeImg = zeros(imgH, imgW) + inf;
holeImg(borderImg) = 0 ;
[~, neighbors] = vl_imdisttf(single(holeImg)) ;
[v_, u_] = ind2sub([imgH, imgW], neighbors);

[u, v] = meshgrid(1:imgW,1:imgH);

holeProb = zeros(imgH, imgW, numPlane);
for i = 1:numPlane
    planeProbCh = planeProb(:,:,i);
    holeProb(:,:,i) = planeProbCh(neighbors);
end
holeProbSum = sum(holeProb, 3);
holeProb = holeProb./repmat(holeProbSum, [1, 1, numPlane]);

maskPlane = repmat(mask, [1,1,numPlane]);
planeProb = (~maskPlane).*planeProb + maskPlane.*holeProb;

modelPlane.postProb = planeProb;



end

function vLine = vLineFromTwoVP(vp1, vp2)

A = cat(1, vp1, vp2);

[U S V] = svd(A, 0);
vLine = V(:,end);
vLine = vLine/vLine(3); % [h7, h8, 1]

end