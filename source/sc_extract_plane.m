function modelPlane = sc_extract_plane(imgName, img, mask, optA)

% SC_EXTRACT_PLANE:
% 
% Extract plane model from an image
%
% Output: modelPlane
% The model plane has the following data structure
%
% modelPlane
% modelPlane.vpData:                        Vanishing point data
% modelPlane.numPlane                       Number of planes
% modelPlane.plane{indPlane}.vLine          The vanishing line of the plane
% modelPlane.plane{indPlane}.imgPlaneProb;  The planar location density
% modelPlane.plane{indPlane}.sourceVP       Which two VPs form the plane
% modelPlane.plane{indPlane}.rotPar(vpInd)  Rotation parameters for aligning
%                                           two sets of lines with the x-axis                                           x-axis
% modelPlane.plane{indPlane}.postProb       Posteior probability of the plane

modelPlane = [];

%% VP detection
imgPath = 'data';
vpFilePath = 'cache\vpdetection';
vpFileName  = [imgName(1:end-4), '-vanishingpoints.txt'];

if(~exist(fullfile(vpFilePath, 'text', vpFileName), 'file'))
    vpExeFile = 'vpdetection.exe';
    vpDetectCMD = [vpExeFile, ' -indir ', imgPath, ' -infile ', imgName, ' -outdir ', vpFilePath];
    system(vpDetectCMD);
end
% Read vanishing point data
vpData = sc_read_vpdata(fullfile(vpFilePath, 'text', vpFileName));

modelPlane = sc_detect_plane_from_vp(vpData, img, mask, optA);

end

function vpData = sc_read_vpdata(fileName)

% SC_READ_VPDATA: read the data from vanishing point detection algorithm
% Input:
%   - fileName: the txt file containing the pre-computed vanishing point
%   detection code
% Output:
%   - vpData
%   The data structure of vpData
%       - vpData.numVP: number of detected vanishing points
%       - vpData.vp{i}.pos: the vanishing point position in the homogenous coordiante
%       - vpData.vp{i}.score: the score of the vanishing point
%       - vpData.vp{i}.numLines: number of lines supporting the vanishing point
%       - vpData.vp{i}.lines{j}.p1: (x1, y1): starting position
%       - vpData.vp{i}.lines{j}.p2: (x2, y2): ending position
%       - vpData.vp{i}.lines{j}.length: length of the line segment

vpData = [];

% Read data
fid = fopen(fileName);

%% Parse VP positions
temp = fscanf(fid, '%s ', [1 5]);
numVP = 0;
readVPFlag = 1;
VP = [];
while(readVPFlag)
    numVP = numVP + 1;
    vpCurr = fscanf(fid, '%g %g %g %g %g', [5 1]);
    if(~isempty(vpCurr))
        VP(:,numVP) = vpCurr;
    else
        temp = fscanf(fid, '%s ', [1 6]);
        readVPFlag = 0;
    end
end
VP = VP';

vpData.numVP = size(VP, 1);

% Save VP position data
for i = 1: vpData.numVP
    vpData.vp{i}.pos = VP(i, 1:3);
    vpData.vp{i}.score = VP(i, 4);
    vpData.vp{i}.numLines = VP(i, 5);
end

%% Parse each set of line segments for the corresponding VP
for i = 1: vpData.numVP
    numLine = fscanf(fid, '%d ', [1 1]);
    lines = fscanf(fid, '%g %g %g %g %g', [5 numLine]);
    vpData.vp{i}.lines = lines';
end

fclose(fid);
end

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
        %         modelPlane.plane{indPlane}.imgPlaneProb(mask) = 0;
        modelPlane.plane{indPlane}.score = sum(modelPlane.plane{indPlane}.imgPlaneProb(:));
        
        %         modelPlane.plane{indPlane}.imgPlaneProb(mask) = 1e-10;
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
modelPlane.plane{indPlane}.imgPlaneProb = optA.fpPlaneProb*ones(imgH, imgW);
% modelPlane.plane{indPlane}.imgPlaneProb(mask) = 0;
modelPlane.plane{indPlane}.score = sum(modelPlane.plane{indPlane}.imgPlaneProb(:));
% modelPlane.plane{indPlane}.imgPlaneProb(mask) = 0;
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

modelPlane.postProbHole = planeProb;

% === Propagate posterior probability into the hole region ===

    % Get border pixels
%     borderImg = im2double(mask);
%     borderImg = imfilter(borderImg, fspecial('gaussian'));
%     borderImg = (borderImg~=0) & (borderImg~=1);
%     borderImg = borderImg & ~mask;
[distMap, idMap] = bwdist(~mask, 'euclidean');

maskInt = mask;
maskInt(1,:) = 0;   maskInt(end,:) = 0;
maskInt(:,1) = 0;   maskInt(:,end) = 0;

for i = 1:numPlane
    planeProbCh = planeProb(:,:,i);
    planeProb(:,:,i) = planeProbCh(idMap);
    planeProb(:,:,i) = roifill(planeProb(:,:,i), maskInt);
end

planeProbSum = sum(planeProb, 3);
planeProb = (planeProb)./repmat(planeProbSum, [1, 1, numPlane]);

planeProbSum = 1 + numPlane*optA.probConst;
planeProb = (planeProb + optA.probConst)/planeProbSum;

% planeProb(:,:,2) = roifill(planeProb(:,:,2), mask);
% roifill(planeProb(:,:,3), mask);
% roifill(planeProb(:,:,4), mask);

if(0)
    
    % Get border pixels
    borderImg = im2double(mask);
    borderImg = imfilter(borderImg, fspecial('gaussian'));
    borderImg = (borderImg~=0) & (borderImg~=1);
    borderImg = borderImg & ~mask;
    
    holeImg = zeros(imgH, imgW) + inf;
    holeImg(borderImg) = 0;
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
    
end

modelPlane.postProb = planeProb;



end

function vLine = vLineFromTwoVP(vp1, vp2)

A = cat(1, vp1, vp2);

[U S V] = svd(A, 0);
vLine = V(:,end);
vLine = vLine/vLine(3); % [h7, h8, 1]

end