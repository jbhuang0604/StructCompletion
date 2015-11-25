function [validPix, uvPix] = sc_init_level(mask, psize, prad)

% Functionality: 
%   get holePixels, uvPixels, validPixels in the level
%
% holePixels: locations of holes pixels
% uvPixels: locations of patches with holes
% validPixels: locations of patches without holes

% =========================================================================
% Get uvPix: center locations of patches with holes
% =========================================================================
uvMask = imdilate(mask, strel('square', double(psize)));
uvMask([1:prad, end-prad+1:end], :) = 0;
uvMask(:, [1:prad, end-prad+1:end]) = 0;

uvPix = getUvPix(uvMask);

% =========================================================================
% Get validPixels: center locations of patches without holes
% =========================================================================
validMap = ~uvMask;
validMap([1:prad,end-prad+1:end], :) = 0;
validMap(:, [1:prad,end-prad+1:end]) = 0;

validPix = getUvPix(validMap);

end

function uvPix = getUvPix(uvMap)

[rUv, cUv] = find(uvMap);
uvPix.sub    = single(cat(2, cUv, rUv));
uvPix.ind    = uint32(sub2ind(size(uvMap), rUv, cUv));
uvPix.mask   = uvMap;
uvPix.numPix = size(uvPix.ind, 1);

end