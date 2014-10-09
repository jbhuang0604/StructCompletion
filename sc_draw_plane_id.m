function uvPlaneIDData= sc_draw_plane_id(planeProbAccData)

% SC_DRAW_PLANE_ID
%
% Input: 
%   - planeProbAccData
% Output:
%   - uvPlaneIDData

numUvPix = size(planeProbAccData, 2);
numPlane = size(planeProbAccData, 1) - 1;

randSample = rand(1, numUvPix);
uvPlaneIDData = zeros(1, numUvPix, 'uint8');

for indPlane = 1: numPlane
    indSamplePlane = (planeProbAccData(indPlane,:) < randSample ) & ...
        (planeProbAccData(indPlane + 1, :) >= randSample);
    uvPlaneIDData(indSamplePlane) = indPlane;
end

end