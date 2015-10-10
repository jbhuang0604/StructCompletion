function uvValidInd = sc_check_valid_uv(srcPos, validSrcMask)

% Rounding to the closest integer
srcPos   = round(srcPos);
numUvPix = size(srcPos, 1);

% Initialize uvValidInd
uvValidInd = false(numUvPix, 1);

% Check only valid 
validSrcInd = (srcPos(:,1) >= 1) & (srcPos(:,1) <= size(validSrcMask,2)) ...
    & (srcPos(:,2) >= 1) & (srcPos(:,2) <= size(validSrcMask,1));

% uvValidInd
uvInd = sub2ind(size(validSrcMask), srcPos(validSrcInd,2), srcPos(validSrcInd,1));
uvValidInd(validSrcInd) = validSrcMask(uvInd);


end