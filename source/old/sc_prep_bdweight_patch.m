function distWeightTable = sc_prep_bdweight_patch(distMap, uvPixels, pSize)

[imgH, imgW] = size(distMap);

numUvPixels = size(uvPixels,1);
uvPixelsInd = sub2ind([imgH, imgW], uvPixels(:,2), uvPixels(:,1));




end