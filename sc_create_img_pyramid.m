function [imgPyr, maskPyr, scaleImgPyr] = sc_create_img_pyramid(img, mask, optS)

% SC_CREAT_IMG_PYRAMID
% 
% Create image pyramid with linear or log scale for coarse to fine image
% completion
% 
% Input: 
%   - img:  Image with hole
%   - mask: Hole mask
%   - optS: options
% Output:  
%   - imgPyr:      Image pyramid
%   - maskPyr:     Mask pyramid
%   - scaleImgPyr: Image dimensions in each level 

% Image size in the high-resolution image
[imgHeight, imgWidth, nCh] = size(img);

% Compute the coarsest image scale
imgSizeMin = min(imgHeight, imgWidth);
coarestScale = optS.coarestImgSize/imgSizeMin;

% Compute the scale in each layer in the image pyramid
if(optS.useLogScale) % use log scale
    scalePyr = 2.^linspace(0, log2(coarestScale), optS.numPyrLvl);
else % use linear scale
    scalePyr = linspace(1, coarestScale, optS.numPyrLvl);
end

% Image size in each layer
imgHPyr = round(imgHeight *scalePyr);
imgWPyr = round(imgWidth  *scalePyr);

% Initialize image pyramid
imgPyr  = cell(optS.numPyrLvl, 1);
maskPyr = cell(optS.numPyrLvl, 1);
scaleImgPyr = cell(optS.numPyrLvl, 1);

% Finest level
imgPyr{1} = img;
maskPyr{1} = logical(mask);
scaleImgPyr{1}.imgScale = 1;
scaleImgPyr{1}.imgSize = [imgHeight, imgWidth];

% Downsampled images
for k = 2: optS.numPyrLvl
    imgPyr{k} = imresize(img, [imgHPyr(k), imgWPyr(k)], optS.resampleKernel);
    maskCur = imresize(mask, [imgHPyr(k), imgWPyr(k)], optS.resampleKernel);
    maskPyr{k} = maskCur ~= 0;
    scaleImgPyr{k}.imgScale = scalePyr(k);
    scaleImgPyr{k}.imgSize  = [imgHPyr(k), imgWPyr(k)];
end

% Get the inital solution
[~, idMap] = bwdist(~maskPyr{optS.numPyrLvl}, 'euclidean');
 
% Intepolate only in the interior to avoid dark values near the image
% borders
maskInt = maskPyr{optS.numPyrLvl};
maskInt(1,:) = 0;   maskInt(end,:) = 0;
maskInt(:,1) = 0;   maskInt(:,end) = 0;
 
for ch = 1: 3
    imgCh = imgPyr{optS.numPyrLvl}(:,:,ch);
    imgCh = imgCh(idMap); 
    imgPyr{optS.numPyrLvl}(:,:,ch) = roifill(imgCh, maskInt);
end

end