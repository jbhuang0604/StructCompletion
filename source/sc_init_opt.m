function [optA, optS] = sc_init_opt()

% SC_INIT_OPT: Initialize parameters for analysis and synthesis
% 

fprintf('- Initialize parameters \n');

optA = [];
optS = [];

%%  Synthesis: guided completion
% =========================================================================
% Patch size
% =========================================================================
optS.pSize = 9;                        % Patch size (odd number), use larger patch for more coherent region
optS.pRad  = floor(optS.pSize/2);      % Patch radius
optS.pNumPix = optS.pSize*optS.pSize;  % Number of pixels in a patch
optS.pMidPix = round(optS.pNumPix/2);  % The center of the patch

% =========================================================================
% Multi-resolution parameters
% =========================================================================
optS.numPyrLvl = 10;                   % Number of coarse to fine layer
optS.coarestImgSize = 32;              % The size of the smallest image in the pyramid
optS.resampleKernel = 'bicubic';       % Resampling kernel: 'lanczos3', 'bicubic', 'bilinear'
optS.useLogScale = 1;                  % Use log scales or linear scales for downsampling

% Weighting parameters for patch match, larger weight put more emphasis on
% pixels near to known pixels
optS.wDist = 2.^linspace(1, 0.0, optS.numPyrLvl);
optS.topLevel = 1;                     % Which level to stop

% Patch weighting kernel
h = fspecial('gaussian', [optS.pSize, optS.pSize], optS.pRad);
h = h(:)/sum(h(:));
optS.wPatch = h;

% =========================================================================
% Parameters for domain transformation and photometric compensation
% =========================================================================
% This scale range is used to reject unlikely patch transformation
optS.minScale = 0.75;                   % Mininum patch scale variation
optS.maxScale = 1.25;                   % Maximum patch scale variation

% Parameters for photometric compensation
optS.minBias = -0.05;                   % Mininum bias compensation
optS.maxBias =  0.05;                   % Maximum bias compensation

% =========================================================================
% Coarse-to-fine iterations
% =========================================================================
% Number of iterations per level
optS.numIter    = 10;                  % The initial iteration number, large hole might require
                                       % more iterations
optS.numIterDec = optS.numIter/optS.numPyrLvl;   % Number of decrements
optS.numIterMin = 3;                   % Minimum number of iterations
optS.numPassPerIter = 1;

% =========================================================================
% Weighting parameters
% =========================================================================
% To-Do: adaptive parameter selection
optS.lambdaCoherence = 1e-2;            % Weighting parameter for coherence

% Planar cost
optS.lambdaPlane  = 5e-2;               % Weighting parameters for planar cost

% Directional cost
optS.lambdaDirect = 5;                  % Weighting parameters for directional cost
optS.directThres  = 0.05;               % Directional cost threshold

optS.lambdaProx   = 0e-2;               % Weighting parameters for proximity cost
optS.proxThres    = 0.25;

optS.lambdaReg    = -0.02;              % Weighting parameters for encouraging
% regularity-guided sampling

% =========================================================================
% Method configuration
% =========================================================================
optS.useCoherence  = 0;                % Not implemented yet
optS.useRegGuide   = 1;
optS.usePlaneGuide = 1;
optS.useBiasCorrection = 1;

% === Precomputed patch position in the reference position ===
[X, Y] = meshgrid(-optS.pRad:optS.pRad, -optS.pRad:optS.pRad);
optS.refPatchPos = single(cat(2, X(:), Y(:), ones(optS.pSize*optS.pSize, 1)));

% === Propagation directions ===
optS.propDir = [1 0; 0 1; -1 0; 0 -1];

% optS.rsThres =  0.05;                  % Random sampling threshold for early termination
% optS.voteUpdateW = 1;                    % Smaller number for quick voting update

optS.numRegSample  = 5;                % Number of regularity-guided sampling per iteration
optS.numRandSample = 5;                % Number of coarse-to-fine random sampling per iteration

% [To-Do] robust norm, e.g., huber
optS.costType = 'L1';                  % Patch appearance cost, other option: 'L2'

%% Analysis: extracting planar structure
% === Plane parameters ===
% Detect up to 3 vanishing points, always include dummy VP (froto-parallel plane)
optA.numVP = 4;

% Fixed density for frontal parallel plane [To-Do] should adapt to image size
optA.fpPlaneProb = 1e-4;

% === Regularity parameters ===
% Maximum number of local features
optA.numFeatMatch = 2000;

% Blurring operations for estimating the plane spatial support
optA.filterSize = 100;
optA.filterSigma = 50;
optA.numFilterIter = 20;

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