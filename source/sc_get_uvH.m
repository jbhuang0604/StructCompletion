function uvH = sc_get_uvH(uvMap, uvPixels)

% Extract uv transform from the NNF
[imgH, imgW, ch] = size(uvMap);

numUvPixels = size(uvPixels,1);
uvPixelsInd = sub2ind([imgH, imgW], uvPixels(:,2), uvPixels(:,1));

uvH = single(zeros(8, numUvPixels));
for i = 1: ch
    uvMapTemp = uvMap(:,:,i);
    uvH(i,:) = uvMapTemp(uvPixelsInd);
end
uvH = cat(1, uvH, ones(1, numUvPixels));
 
end