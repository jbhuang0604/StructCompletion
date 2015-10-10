function imgBlend = sc_poisson_blend(imgTrg, imgSrc, holeMask)
% SC_POISSON_BLEND: Blend the source and target image using a simple 
% discrete Poisson solver
%
% Input:
%   - imgTrg: original image with hole
%   - imgSrc: synthesized image
%   - holeMask: specify the hole pixels
% Output:
%   - imgBlend: blended image

[imgH, imgW, nCh] = size(imgSrc);

% Initialize the reconstructed image
imgRecon = zeros(imgH, imgW, nCh);

% Independently process each channel
for ch = 1: nCh
    % Source and target image
    imgS = imgSrc(:,:,ch);
    imgT = imgTrg(:,:,ch);
    
    % prepare discrete Poisson equation
    [A, b] = solvePoisson(holeMask, imgS, imgT);

    % solve Poisson equation
    x = A\b;
    imgRecon(:,:,ch) = reshape(x, [imgH, imgW]);
end

% Combined with the known region in the target
holeMaskC = cat(3, holeMask, holeMask, holeMask);
imgBlend = holeMaskC.*imgRecon + ~holeMaskC.*imgTrg;

end

function [A, b] = solvePoisson(holeMask, imgS, imgT)

% Prepare the linear system of equations for Poisson blending

[imgH, imgW] = size(holeMask);
N = imgH*imgW;

% Number of unknown variables
numUnknownPix = sum(holeMask(:));

maxNumInd = 8*numUnknownPix; % Max number of indices for constructing the sparse matrix
maxNumEqn = 4*numUnknownPix; % Max number of equations (4-neighbor)

% 4-neighbors: dx and dy
dx = [1, 0, -1,  0];    dy = [0, 1,  0, -1];

% Initialize (I, J, S), for sparse matrix A where A(I(k), J(k)) = S(k)
I = zeros(maxNumInd, 1);
J = zeros(maxNumInd, 1);
S = zeros(maxNumInd, 1);

% Initialize b
b = zeros(maxNumEqn, 1);

% Precompute unkonwn pixel position
[pi, pj] = find(holeMask == 1);
pind = sub2ind([imgH, imgW], pi, pj);

% Precompute the 4-neighbor of the unkonwn pixel positions
qi = bsxfun(@plus, pi, dy);
qj = bsxfun(@plus, pj, dx);

% Handling cases at image borders
validN = (qi >= 1) & (qi <= imgH) & (qj >= 1) & (qj <= imgW);
qind = zeros(size(validN));
qind(validN) = sub2ind([imgH, imgW], qi(validN), qj(validN));

c = 1; % index counter for matrix A
e = 1; % equation counter

% Set up the matrix A and the vector b
for k = 1: numUnknownPix
    pind_cur = pind(k);
    for n = 1: 4 % 4-neighbor
        if(validN(k, n)) % if the neighbor pixel q lies in the image
            qind_cur = qind(k, n);
            if(holeMask(qind_cur))
                % A(e, pind_cur) = 1;
                I(c) = e;   J(c) = pind_cur;    S(c) = 1;   c = c + 1;
                % A(e, qind_cur) = 1;
                I(c) = e;   J(c) = qind_cur;    S(c) = -1;  c = c + 1;
                
                % gradient constraint
                b(e) = imgS(pind_cur) - imgS(qind_cur);
                e = e + 1;
            else
                % A(e, pind_cur) = 1;
                I(c) = e;   J(c) = pind_cur;    S(c) = 1;   c = c + 1;
                
                % boundary constraint
                b(e) = imgS(pind_cur) - imgS(qind_cur) + imgT(qind_cur);
                e = e + 1;
            end
        end
    end
end

% Clean up unused entries
nEqn = e - 1;
b = b(1:e-1);

nInd = c - 1;
I = I(1:nInd);  J = J(1:nInd);  S = S(1:nInd);

% Construct the sparse matrix A
A = sparse(I, J, S, nEqn, N);

end