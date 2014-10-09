function imgCompletion = sc_complete(imgFileName, optA, optS)
% SC_COMPLETE 
%
% Single image completion with user-defined mask
% 
% Input: 
%     - imgFileName: image file name
%     - optA: parameters for analysis
%     - optS: parameters for synthesis
% Output:
%     - imgCompletion: completed result
% 
% Example:
%   imgFileName = '002_input_hole.png';
%   [optA, optS] = sc_init_opt;
%   imgCompletion = sc_complete(imgFileName, optA, optS);
%   
% Disclaimer: 
%   This is a Matlab re-implementation of the paper:
% 
%   Jia-Bin Huang, Sing Bing Kang, Narendra Ahuja, and Johannes Kopf, 
%   Image Completion using Planar Structure Guidance,
%   ACM Transactions on Graphics (Proceedings of SIGGRAPH 2014), 33(4), 2014
%  
%   It is provided for educational/researrch purpose only.
%   Please cite our paper if you found the software useful for your work.
%   
%   
% Contact:   
%   Jia-Bin Huang
%   University of Illinois, Urbana-Champaign 
%   www.jiabinhuang.com  
%   jbhuang0604@gmail.com 
    
%% Temp input
% rand('seed', 1);
imgID = 46; 
imgFileName = [num2str(imgID, '%03d'), '_input_hole.png'];
addpath(genpath('external'));
addpath('external\vlfeat-0.9.19\toolbox\mex\mexw64');

% Option parameters
fprintf('- Initialize parameters \n');
[optA, optS] = sc_init_opt;

%% === Planar structure extraction ===
fprintf('- Extract planar structures \n'); 
tic;
[img, mask, modelPlane, modelReg] = sc_extract_planar_structure(imgFileName, optA);
mask = im2double(mask);
tAnalysis = toc;
fprintf('Done in %6.3f seconds.\n\n', tAnalysis);

%% === Guided image completion ===
% Construct image pyramid for coarse-to-fine image completion
fprintf('- Construct image pyramid: \n'); 
tic; 
% Image
[imgPyr, maskPyr, scaleImgPyr] = sc_create_img_pyramid(img, mask, optS);
% Structure constraints
[modelPlane, modelReg] = sc_planar_structure_pyramid(scaleImgPyr, modelPlane, modelReg);
tImgPyramid = toc;   
fprintf('Done in %6.3f seconds.\n\n', tImgPyramid);
    
% Completion by synthesis
fprintf('- Image completion using planar structure guidance\n'); 
tic; 
[imgPyr, imgPyrNNF] = sc_synthesis(imgPyr, maskPyr, modelPlane, modelReg, optS);
tSynthesis = toc;
fprintf('Synthesis took %6.3f seconds.\n', tSynthesis);
 
% Return the top level 
imgCompletion = imgPyr{optS.topLevel};

% Visualizing completion process 
% To-Do
% sc_vis_synthesis_process(imgCompletePyr, imgPyrNNF, imgPyrCost);

end