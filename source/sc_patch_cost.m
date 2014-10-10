function [costPatchCand, uvBiasCand] = ...
    sc_patch_cost(trgPatchCur, srcPatch, wDistPatchCur, modelPlane, uvPlaneIDData, ...
    trgPos, srcPos, bdPos, uvPixUpdateSrc, imgSize, optS, iLvl)

% Patch cost - appearance
[costApp, uvBiasCand] = sc_patch_cost_app(trgPatchCur, srcPatch, wDistPatchCur, optS);

% Patch cost - planar guidance
costPlane = sc_patch_cost_plane(modelPlane.mLogLPlaneProb, uvPlaneIDData, trgPos, srcPos);

% Patch cost - directional guidance
costDirect = sc_patch_cost_direct(uvPlaneIDData, trgPos, srcPos, modelPlane, imgSize, optS);

% Patch cost - proximity cost
costProx = sc_patch_cost_prox(srcPos, trgPos, bdPos, imgSize, optS);

costReg = double(uvPixUpdateSrc == 1);

% Graduately tuned down the structural guidance cost
w = (iLvl)/optS.numPyrLvl;
w = 1;

% Weighted sum of the appearance cost and guidance cost
costPatchCand = costApp + ...
    w*(optS.lambdaPlane  * costPlane + ...
     optS.lambdaReg    * costReg + ...
     optS.lambdaDirect * costDirect + ...
     optS.lambdaProx   * costProx);


end
