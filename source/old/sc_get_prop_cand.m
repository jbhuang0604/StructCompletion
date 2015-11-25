function uvHCurrCand = sc_get_prop_cand(uvHCurr, dx, dy)

% Get propagation neighbor

uvHCurrCand = zeros(size(uvHCurr));

Ht = [1, 0, dx; 0, 1, dy; 0 0 1];
numUvPixels = size(uvHCurr, 2);

% for i = 1: numUvPixels
%     H = reshape(uvHCurr(:,i), 3, 3);
%     Hc = H*Ht;
%     uvHCurrCand(:, i) = Hc(:);
% end
% 
% uvHCurrCand = uvHCurrCand./repmat(uvHCurrCand(9,:), 9, 1);

end