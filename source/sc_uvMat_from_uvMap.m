function uvMat = sc_uvMat_from_uvMap(map, uvPixInd)

% Function: conver uvMap into matrix form

% Extract uv transform from the NNF
[imgH, imgW, nCh] = size(map);

offset = uint32((0:nCh-1)*imgH*imgW);
uvPixInd = bsxfun(@plus, uvPixInd, offset);
uvMat = map(uvPixInd);

% map(uvPixInd) = data;
% 
% 
% 
% uvMat = zeros(nCh, numUvPixels, 'single');
% 
% for i = 1: nCh
%     uvMat(i,:) = map(uvPix.ind + (i-1)*imgH*imgW);
% end

% Indices of valid uv pixels
% uvValid.ind = uvMapMask(uvPix.ind);
% uvValid.pos = find(uvValid.ind);

end