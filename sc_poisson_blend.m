function imgBlend = sc_poisson_blend(img, holeMask, distMap, optS)

% SC_POISSON_BLEND
% Blend the source and target image
% Input:
%   - img: completed image
%   - holeMask: hask specify the hole pixels
%   - distMap: distanace transform
%   - optS
% Output:
%   - imgBlend: blended image

[imgH, imgW, nCh] = size(img);
img = sc_clamp(img, 0, 1);

% Compute gradient x and gradient y
hx = [-1, 1];
hy = hx';

imgGx = imfilter(img, hx, 'replicate', 'same');
imgGy = imfilter(img, hy, 'replicate', 'same');

% Specify known region
knowMask = ~holeMask;

% Specify boundary region
bdMask = distMap == 1;

% Weight for source and target region
% weightTrg = 1;

weightSrc = (1e2*knowMask + 1);
weightTrg = holeMask;

% Create the gradient operator
[Gx, Gy] = CreatGradOperator(imgH, imgW);

n = imgH*imgW;
% Id = sparse( 1:n, 1:n, ones(n,1));
Id = sparse( 1:n, 1:n, weightSrc(:));

% Apply the weight
% Gx = weightTrg*Gx(holeMask(:), :);
% Gy = weightTrg*Gy(holeMask(:), :);
% Id = Id(knowMask(:),:);

% Prepare the linear system
A = cat(1, Gx, Gy, Id);
At = A';
AtA = At*A;

imgBlend = zeros(imgH, imgW, nCh);

% Indepedently blend each color channels
for ch = 1: nCh
    % Gradient and color image
    imgGxCh    = imgGx(:,:,ch);
    imgGyCh    = imgGy(:,:,ch);
    imgColorCh = img(:,:,ch);
    
    % Set the gradient at the boundary to be zeros
    imgGxCh(bdMask) = 0;
    imgGyCh(bdMask) = 0;
    
    % Select the source and target data
%     imgGxCh    = imgGxCh(holeMask(:));
%     imgGyCh    = imgGyCh(holeMask(:));
%     imgColorCh = imgColorCh(knowMask(:));
    
    % Prepare the linear system
    b = cat(1, weightTrg(:).*imgGxCh(:), weightTrg(:).*imgGyCh(:), weightSrc(:).*imgColorCh(:));
    
    % Solve the weighted least square
    imgRes = (AtA) \ (At*b);
    imgRes = reshape(imgRes, imgH, imgW);
    imgBlend(:,:,ch) = imgRes;
end

% figure(1),
% subplot(1,3,1); imshow(img);
% subplot(1,3,2); imshow(imgBlend);
% subplot(1,3,3); imshow(imgBlend-img+0.5);

end

function [Gx, Gy] = CreatGradOperator(imgH, imgW)

D = spdiags([-ones(imgH,1), ones(imgH,1)], [0 1], imgH, imgH + 1);
D(:,end) = [];
D(end,end) = 0;

Gy = kron(speye(imgW),D);

D = spdiags([-ones(imgW,1), ones(imgW,1)], [0 1], imgW, imgW + 1);
D(:,end) = [];
D(end,end) = 0;

Gx = kron(D, speye(imgH));


end