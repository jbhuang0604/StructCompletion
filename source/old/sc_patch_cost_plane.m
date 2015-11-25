
function costPlane = sc_patch_cost_plane(mLogLPlaneProb, uvPlaneIDData, trgPixSub, srcPixSub)

% SC_PATCH_COST_PLANE
% 
% Compute planar costs (See Eqn 11 in the paper)

[imgH, imgW, numPlane] = size(mLogLPlaneProb);

srcPixSub     = round(srcPixSub);
uvPlaneIDData = double(uvPlaneIDData);

trgPixIndCur = sub2ind([imgH, imgW, numPlane], trgPixSub(:,2), trgPixSub(:,1), uvPlaneIDData);
srcPixIndCur = sub2ind([imgH, imgW, numPlane], srcPixSub(:,2), srcPixSub(:,1), uvPlaneIDData);

costPlane = mLogLPlaneProb(trgPixIndCur) + mLogLPlaneProb(srcPixIndCur);

end 