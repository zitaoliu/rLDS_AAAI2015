function [ psdMatrix, eigVecs, eigVals ] = psdProjection( matrix )
%COVCONVERT Summary of this function goes here
%   Detailed explanation goes here

matrix = matrix ./ 2;

matrix = matrix + matrix';

% if sum(sum(isnan(matrix))) > 0 || sum(sum(isinf(matrix))) > 0
%     disp('')
% end

[eigVecs, eigVals] = eig(matrix);

eigValsVec = diag(eigVals);

if sum(eigValsVec < 0) > 0 && abs(min(eigValsVec)) > 1e-15 
    eigVals(eigVals < 0) = eps;
    psdMatrix = eigVecs * eigVals * eigVecs';
else
    psdMatrix = matrix;
end

end

