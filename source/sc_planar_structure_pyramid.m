function [modelPlanePyr, modelRegPyr] = sc_planar_structure_pyramid(scaleImgPyr, modelPlane, modelReg)

%
% SC_PLANAR_STRUCTURE_PYRAMID:
%
% Pre-compute the models of planes and regularity according to the image
% pyramid, rescale parameters according to the image sizes.
%
% Output: 
%   - modelPlanePyr: plane model
%   - modelRegPyr: regularity model

numLevel = length(scaleImgPyr);

modelPlanePyr = cell(numLevel, 1);
modelRegPyr = cell(numLevel, 1);

planePostProb = modelPlane.postProb;

for iLvl = 1: numLevel
    % === Update plane model ===
    scaleImgCur = scaleImgPyr{iLvl}.imgScale;
    sizeImgCur = scaleImgPyr{iLvl}.imgSize;
    
    % Update posterior plane probability
    planePostProbCur = imresize(planePostProb, sizeImgCur, 'bicubic');
    planePostProbCurSum = sum(planePostProbCur, 3);
    planePostProbCur = bsxfun(@rdivide, planePostProbCur, planePostProbCurSum);
    modelPlanePyr{iLvl}.numPlane = modelPlane.numPlane;
    modelPlanePyr{iLvl}.planeProb = planePostProbCur;
    modelPlanePyr{iLvl}.mLogLPlaneProb = -log(planePostProbCur);
    
    % Update rectification matrix and rotation parameters
    for iPlane = 1: modelPlane.numPlane
        H = eye(3);
        vLine = modelPlane.plane{iPlane}.vLine;
        vLine(1:2) = vLine(1:2)/scaleImgCur;
        H(3,:) = vLine;
        modelPlanePyr{iLvl}.rectMat{iPlane} = H;
        modelPlanePyr{iLvl}.rotPar{iPlane} = modelPlane.plane{iPlane}.rotPar;

        % Construct the rotation matrices
        for iTheta = 1:2
            t = modelPlanePyr{iLvl}.rotPar{iPlane}(iTheta);
            Hr = eye(3);
            Hr(1,1) = cos(t);   Hr(1,2) = - sin(t);
            Hr(2,1) = sin(t);   Hr(2,2) =   cos(t);
            modelPlanePyr{iLvl}.rotMat{iPlane, iTheta} = Hr;
        end
    end
     
    % === Update reguarlity model ===
    for iPlane = 1: modelPlane.numPlane
        modelRegPyr{iLvl}.dispVec{iPlane} = scaleImgCur*modelReg.plane{iPlane}.dispVec;
        modelRegPyr{iLvl}.numDispVec(iPlane) = modelReg.plane{iPlane}.numDispVec;
    end
end

end