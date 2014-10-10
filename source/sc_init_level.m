function [validPix, uvPix] = sc_init_level(mask, psize, prad)

% Functionality: 
%   get holePixels, uvPixels, validPixels in the level
%
% holePixels: locations of holes pixels
% uvPixels: locations of patches with holes
% validPixels: locations of patches without holes
% mask = (mask);

[imgH, imgW] = size(mask);

% % Get holePixels
% [rH, cH] = find(mask);
% holePix.sub = cat(2, cH, rH)';
% holePix.ind = sub2ind([imgH, imgW], rH, cH)';

% Get  uvPixSub, uvPixInd
uvMaps = imdilate(mask, strel('square', double(psize)));
uvMaps([1:prad, end-prad+1:end], :) = 0;
uvMaps(:, [1:prad, end-prad+1:end]) = 0;
[rUv, cUv] = find(uvMaps);
uvPix.sub = cat(2, cUv, rUv)';
uvPix.ind = sub2ind([imgH, imgW], rUv, cUv)';
uvPix.mask = uvMaps;
uvPix.numUvPix = size(uvPix.ind, 2);

% Get validPixels
validMap = ~uvMaps;
validMap([1:prad,end-prad+1:end], :) = 0;
validMap(:, [1:prad,end-prad+1:end]) = 0;
[rV, cV] = find(validMap);
validPix.sub = cat(2, cV, rV)';
validPix.ind = sub2ind([imgH, imgW], rV, cV)';
validPix.mask = validMap;
validPix.numValidPix = size(validPix.ind, 2);

end