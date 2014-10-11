function [optA, optS] = sc_init_opt()

% SC_INIT_OPT
% Initialize parameters for analysis and synthesis

fprintf('- Initialize parameters \n');

optA = [];
optS = [];

%%  Synthesis: guided completion

% === Patch size ===
optS.pSize = 7;                        % Patch size (odd number), use larger patch for more coherent region
optS.pRad  = floor(optS.pSize/2);      % Patch radius
optS.pNumPix = optS.pSize*optS.pSize;  % Number of pixels in a patch
optS.pMidPix = round(optS.pNumPix/2);  % The center of the patch

% === Multi-resolution parameters ===
optS.numPyrLvl = 10;                   % Number of coarse to fine layer
optS.coarestImgSize = 32;              % The size of the smallest image in the pyramid
optS.resampleKernel = 'lanczos3';      % Resampling kernel: 'lanczos3', 'bicubic', 'bilinear'
optS.useLogScale = 1;                  % Use log scales or linear scales for downsampling

% Weighting parameters for patch match, larger weight put more emphasis on
% pixels near to known pixels
optS.wDist = 2.^linspace(1, 1, optS.numPyrLvl);

optS.topLevel = 1;                     % Which level to stop
optS.propOnlyLevel = 0;

% === Parameters for domain transformation ===
% This scale range is used to reject unlikely patch transformation
optS.minScale = 0.25;                  % Mininum patch scale variation
optS.maxScale = 4;                     % Maximum patch scale variation

% === Parameters for photometric compensation ===
optS.minBias = -0.1;                   % Mininum bias compensation
optS.maxBias =  0.1;                   % Maximum bias compensation

% === Number of iterations per level ===
optS.numIter    = 30;                  % The initial iteration number, large hole might require
% more iterations
optS.numIterDec = optS.numIter/optS.numPyrLvl;   % Number of decrements
optS.numIterMin = optS.numIter/optS.numPyrLvl;   % Minimum number of iterations
optS.numPassPerIter = 1;

% === Weighting parameters ===
% To-Do: adaptive parameter selection
optS.lambdaPlane  = 0.05;               % Weighting parameters for planar cost

optS.lambdaDirect = 10;                  % Weighting parameters for directional cost
optS.directThres  = 0.02;

optS.lambdaProx   = 0.2;               % Weighting parameters for proximity cost
optS.proxThres    = 0.0;

optS.lambdaReg    = -0.02;             % Weighting parameters for encouraging
% regularity-guided sampling

% === Method configuration ===
optS.useCoherence  = 1;                % Not implemented yet
optS.useRegGuide   = 1;
optS.usePlaneGuide = 1;
optS.useBiasCorrection = 1;

% === Precomputed patch position in the reference position ===
[X, Y] = meshgrid(-optS.pRad:optS.pRad, -optS.pRad:optS.pRad);
optS.refPatchPos = single(cat(1, X(:)', Y(:)', ones(1, optS.pSize*optS.pSize)));

% === Propagation directions ===
optS.propDir = [1 0; 0 1; -1 0; 0 -1]';

optS.rsThres = 0.05;                   % Random sampling threshold for early termination
optS.voteUpdateW = 20;                 % Smaller number for quick voting update

optS.numRegSample  = 5;                % Number of regularity-guided sampling per iteration
optS.numRandSample = 5;                % Number of coarse-to-fine random sampling per iteration

% [To-Do] robust norm, e.g., huber
optS.costType = 'L1';                  % Patch appearance cost, other option: 'L2'

%% Analysis: extracting planar structure
% === Plane parameters ===
% Detect up to 3 vanishing points, always include dummy VP (froto-parallel plane)
optA.numVP = 4;
% Fixed density for frontal parallel plane [To-Do] should adapt to image size
optA.fpPlaneProb = 1e-5;

% === Regularity parameters ===
% Maximum number of local features
optA.numFeatMatch = 2000;

% Blurring operations for estimating the plane spatial support
optA.filterSize = 100;
optA.filterSigma = 50;
optA.numFilterIter = 10;

% Add the constant to all plane posterior probabilities
optA.probConst = 0.05;

% Parameters for feature extraction and matching
optA.PeakThreshold = 0.04;
optA.maxNumLocalFeat = 1000;
optA.numQueryNN = 3;
optA.maxDispVec = 1000;
optA.minDistDispVec = 100;

% Parameters for the Mean-Shift clustering algorithm
optA.msBandwidth = 20;
optA.msMinNumClusterMember = 5;

% Threshold for determining whether a pair of feature matches is on a plane
optA.prodProbThres = 0.5;
end