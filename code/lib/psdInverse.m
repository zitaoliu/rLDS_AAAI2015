function [ invMat ] = psdInverse( eigVecs, eigVals )
%PSDINV Summary of this function goes here
%   Detailed explanation goes here

eigValsVec = diag(eigVals);
idx = eigValsVec > 1e-4;

invMat = eigVecs(:, idx) * diag(1./eigValsVec(idx)) * eigVecs(:, idx)';

end

