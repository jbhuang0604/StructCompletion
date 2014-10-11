function [NNF, nUpdate]= sc_update_NNF(trgPatch, wDistPatch, img, NNF, modelPlane, modelReg, iLvl, optS)

% SC_UPDATE_NNF
%
% Update the nearest neighbor field using the PatchMatch algorithm
%
% Input:
%   - trgPatch
%   - wDistPatch
%   - img,
%   - NNF
%   - modelPlane
%   - modelReg
%   - optS
% Output:
%   - NNF
%   - nUpdate

nUpdate = zeros(1,3);

for i = 1:optS.numPassPerIter
    % propagate along four directions
    for iDirect = 1:4
        [NNF, n] = sc_propagate(trgPatch, wDistPatch, img, NNF, modelPlane, optS, iDirect, iLvl);
        nUpdate(1) = nUpdate(1) + n;
    end
    
    if(iLvl > optS.propOnlyLevel)
        % Random sampling
        [NNF, n] = sc_regular_search(trgPatch, wDistPatch, img, NNF, modelPlane, modelReg, optS, iLvl);
        nUpdate(3) = nUpdate(3) + n;
        
        % Regularity guided sampling
        [NNF, n] = sc_random_search(trgPatch, wDistPatch, img, NNF, modelPlane, optS, iLvl);
        nUpdate(2) = nUpdate(2) + n;
    end
end

end