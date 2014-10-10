function     [img, NNF, wDistPatch, wDistImg] = sc_init_lvl_nnf(img, NNF, holeMask, modelPlaneCur, modelRegCur, iLvl, optS)

% SC_INIT_LVL_NNF

% Initialize the nearest neighbor field for the current level

% Input: 
%   - NNF, img, holeMask, modelPlaneCur, modelRegCur, iLvl, optS
% Output:
%   - img
%   - NNF 
%   - wDistPatch
%   - wDistImg

% Prepare distance weight
[distMap, idMap] = bwdist(~holeMask, 'euclidean');

if(iLvl == optS.numPyrLvl)
    % Initialize the NNF for the coarest level using random sampling
    NNF = sc_init_nnf(holeMask, modelPlaneCur, modelRegCur, optS);
    [wDistPatch, wDistImg] = sc_prep_dist_patch(distMap, NNF.uvPix.sub, iLvl, optS);
else
    % Initialize the NNF upsampling of NNF from previous level
    NNF = sc_upsample(holeMask, NNF, modelPlaneCur, modelRegCur, optS);
    [wDistPatch, wDistImg] = sc_prep_dist_patch(distMap, NNF.uvPix.sub, iLvl, optS);
    img = sc_voting(img, NNF, NNF.uvPix, holeMask, wDistPatch, wDistImg, optS);
end
 
NNF.uvDtBdPixPos = double(distMap(NNF.uvPix.ind));

end