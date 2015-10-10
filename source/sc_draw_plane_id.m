function uvPlaneIDData= sc_draw_plane_id(planeProbAccData)

% SC_DRAW_PLANE_ID
%
% Input: 
%   - planeProbAccData
% Output:
%   - uvPlaneIDData

numUvPix = size(planeProbAccData, 1);
numPlane = size(planeProbAccData, 2) - 1;

randSample = rand(numUvPix, 1);
uvPlaneIDData = zeros(numUvPix, 1, 'uint8');

for indPlane = 1: numPlane
    indSamplePlane = (planeProbAccData(:,indPlane) < randSample ) & ...
        (planeProbAccData(:, indPlane + 1) >= randSample);
    uvPlaneIDData(indSamplePlane) = indPlane;
end

end