function costProx = sc_patch_cost_prox(srcPos, trgPos, uvDtBdPixPos, sizeImg, optS)

% SC_PATCH_COST_PROX

% Encourage to sample patches near to the target patch
%

d = srcPos - trgPos;
d = sqrt(sum(d.^2,1));

d = d/sizeImg;
uvDtBdPixPos = uvDtBdPixPos/sizeImg;

% costProx = max(0, d - optS.proxThres);
% Shrinkage thresholding?
costProx = max(0, d - uvDtBdPixPos - optS.proxThres);
% costProx = max(0, d - optS.proxThres);

end