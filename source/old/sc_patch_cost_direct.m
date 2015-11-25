function costDirect = sc_patch_cost_direct(uvPlaneIDData, trgPos, srcPos, modelPlane, optS)

% SC_PATCH_COST_DIRECT
%
% Compute the directional cost (See Eqn 13 in the paper)
%
% Input
%   -
% Output
%   -

costDirect = optS.lambdaDirect*ones(2, size(trgPos, 2));

for indPlane = 1: modelPlane.numPlane
    % Retrieve the uv pixels that have the current plane label
    uvPlaneIndCur  = uvPlaneIDData == indPlane;
    numPlanePixCur = sum(uvPlaneIndCur);
    
    if(indPlane == modelPlane.numPlane)
        costDirect(:, uvPlaneIndCur) = optS.sizeImg*optS.directThres;
    else
        % The rectifying transformation for the plane
        rectMat = modelPlane.rectMat{indPlane};
        h7 = rectMat(3,1);
        h8 = rectMat(3,2);
        
        if(numPlanePixCur~=0)
            trgPosCur = trgPos(:, uvPlaneIndCur) - 1;
            srcPosCur = srcPos(:, uvPlaneIndCur) - 1;
            
            for iTheta = 1:2
                rotMat = modelPlane.rotMat{indPlane, iTheta};
                
                rotRecMat = rotMat;
                rotRecMat(3,1) = h7;    rotRecMat(3,2) = h8;
                
                trgPosCurRect = rotRecMat*cat(1, trgPosCur, ones(1, numPlanePixCur));
                trgPosCurRect = bsxfun(@rdivide, trgPosCurRect, trgPosCurRect(3,:));
                
                % Source patch center position in the rectified domain
                srcPosCurRect = rotRecMat*cat(1, srcPosCur, ones(1, numPlanePixCur));
                srcPosCurRect = bsxfun(@rdivide, srcPosCurRect, srcPosCurRect(3,:));
                
                costDirect(iTheta, uvPlaneIndCur) = abs(srcPosCurRect(2,:) - trgPosCurRect(2,:));
            end
        end
    end
end

costDirect = min(costDirect, [], 1);
costDirect = costDirect/sizeImg;
costDirect(costDirect > optS.directThres) = optS.directThres;


end