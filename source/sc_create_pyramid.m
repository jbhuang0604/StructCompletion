function [imgPyr, maskPyr, scaleImgPyr] = sc_create_pyramid(img, mask, optS)

% SC_CREAT_IMG_PYRAMID
%
% Create image pyramid with linear or log scale for coarse to fine image
% completion
%
% Input:
%   - img:  input image with hole
%   - mask: Hole mask
%   - optS: options
% Output:
%   - imgPyr:      Image pyramid
%   - maskPyr:     Mask pyramid
%   - scaleImgPyr: Image dimensions in each level

% Image size in the high-resolution image
[imgHeight, imgWidth, nCh] = size(img);
img = sc_init_coarsest_level(img, logical(mask));

% =========================================================================
% Create pyramid: scale
% =========================================================================
scaleImgPyr = sc_create_scale_pyramid(imgHeight, imgWidth, optS);

% =========================================================================
% Create pyramid: mask
% =========================================================================
maskPyr = sc_create_image_pyramid(mask, scaleImgPyr, 'mask', optS); 

% =========================================================================
% Create pyramid: image
% =========================================================================
imgPyr = sc_create_image_pyramid(img, scaleImgPyr, 'image', optS); 

% =========================================================================
% Recover image boundary 
% =========================================================================

% nCh = 3;
% for iLvl = 1: optS.numPyrLvl
%     maskCur  = maskPyr{iLvl};
%     bdRegion = maskCur < 0.99 & maskCur > 0.1;
%     bdRegionC = bdRegion(:,:,ones(nCh,1));
%     
%     imgCur   = imgPyr{iLvl};
%     imgCurBd = bsxfun(@rdivide, imgCur, maskCur);
%     imgCur(bdRegionC) = imgCurBd(bdRegionC);
%     imgPyr{iLvl} = imgCur;
% end


% Initialize the coarsest level
imgPyr{optS.numPyrLvl} = sc_init_coarsest_level(imgPyr{optS.numPyrLvl}, maskPyr{optS.numPyrLvl});

% Convert to single type
imgPyr = cellfun(@im2single, imgPyr, 'UniformOutput', false);

end

function img = sc_init_coarsest_level(img, mask)

% Get the inital solution
[~, idMap] = bwdist(~mask, 'euclidean');

% Intepolate only in the interior to avoid dark values near the image borders
maskInt = mask;
maskInt(1,:) = 0;   maskInt(end,:) = 0;
maskInt(:,1) = 0;   maskInt(:,end) = 0;

for ch = 1: 3
    imgCh = img(:,:,ch);
    imgCh = imgCh(idMap);
    img(:,:,ch) = roifill(imgCh, maskInt);
end

end

function imgPyr = sc_create_image_pyramid(img, scaleImgPyr, imageType, optS)

% h = fspecial('gaussian', 5, 1);

% Initialize image pyramid
imgPyr  = cell(optS.numPyrLvl, 1);

% The finest level
imgPyr{1} = img;

%
for iLvl = 2: optS.numPyrLvl
    imgHCurLvl = scaleImgPyr{iLvl}.imgSize(1);
    imgWCurLvl = scaleImgPyr{iLvl}.imgSize(2);
    
    % Previous layer
    imgCur = imgPyr{iLvl - 1};

    % Anti-alising by blurring
%     imgCur   = imfilter(imgCur, h, 'same', 'replicate', 'conv');
      
    % Resampling
    imgPyr{iLvl} = imresize(imgCur,  [imgHCurLvl, imgWCurLvl], optS.resampleKernel);
end


if(strcmp(imageType, 'mask'))
    % Convert resampled masks into logical type
    for iLvl = 1: optS.numPyrLvl
        imgPyr{iLvl} = imgPyr{iLvl} > 0.5;
    end
end

end

function scaleImgPyr = sc_create_scale_pyramid(imgHeight, imgWidth, optS)

% Compute the coarsest image scale
imgSizeMin   = min(imgHeight, imgWidth);
coarestScale = optS.coarestImgSize/imgSizeMin;

% Compute the scale in each layer in the image pyramid
if(optS.useLogScale)      % use log scale
    scalePyr = 2.^linspace(0, log2(coarestScale), optS.numPyrLvl);
else                      % use linear scale
    scalePyr = linspace(1, coarestScale, optS.numPyrLvl);
end

% Image size in each layer
imgHPyr = round(imgHeight *scalePyr);
imgWPyr = round(imgWidth  *scalePyr);

% Initialize scales
scaleImgPyr = cell(optS.numPyrLvl, 1);

% Finest level
scaleImgPyr{1}.imgScale = 1;
scaleImgPyr{1}.imgSize = [imgHeight, imgWidth];

% Downsampled images
for k = 2: optS.numPyrLvl
    scaleImgPyr{k}.imgScale = scalePyr(k);
    scaleImgPyr{k}.imgSize  = [imgHPyr(k), imgWPyr(k)];
end

end