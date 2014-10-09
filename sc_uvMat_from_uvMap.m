function [uvMat] = sc_uvMat_from_uvMap(uvMap, uvPix)

% Function: conver uvMap into matrix form

% Extract uv transform from the NNF
numUvPixels = size(uvPix.ind,2);
[imgH, imgW, nCh] = size(uvMap);

uvMat = zeros(nCh, numUvPixels, 'single');

for i = 1: nCh
    uvMat(i,:) = uvMap(uvPix.ind + (i-1)*imgH*imgW);
end

% Indices of valid uv pixels
% uvValid.ind = uvMapMask(uvPix.ind);
% uvValid.pos = find(uvValid.ind);

end