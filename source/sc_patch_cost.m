function [costPatchCand, uvBiasCand] = ...
    sc_patch_cost(trgPatchCur, srcPatch, wDistPatchCur, modelPlane, uvPlaneIDData, ...
    trgPos, srcPos, bdPos, uvPixUpdateSrc, imgSize, optS, iLvl)

numUvPix = size(wDistPatchCur, 2);
costPatchCand = zeros(5, numUvPix);

% Patch cost - appearance
[costApp, uvBiasCand] = sc_patch_cost_app(trgPatchCur, srcPatch, wDistPatchCur, optS);

% Patch cost - planar guidance
costPlane = sc_patch_cost_plane(modelPlane.mLogLPlaneProb, uvPlaneIDData, trgPos, srcPos);

% Patch cost - directional guidance
costDirect = sc_patch_cost_direct(uvPlaneIDData, trgPos, srcPos, modelPlane, imgSize, optS);

% Patch cost - proximity cost
costProx = sc_patch_cost_prox(srcPos, trgPos, bdPos, imgSize, optS);

costReg = double(uvPixUpdateSrc == 1);

% Weighted sum of the appearance cost and guidance cost
costPatchCand(1,:) = costApp;
costPatchCand(2,:) = optS.lambdaPlane  * costPlane;
costPatchCand(3,:) = optS.lambdaReg    * costReg;
costPatchCand(4,:) = optS.lambdaDirect * costDirect;
costPatchCand(5,:) = optS.lambdaProx   * costProx;

% [To-Do] Adaptive weighting for patch cost
% costPatchCand(2:5,:) = (iLvl/optS.numPyrLvl)*costPatchCand(2:5,:);

end
