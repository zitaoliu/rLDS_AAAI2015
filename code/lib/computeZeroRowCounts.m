function [ count, matrix ] = computeZeroRowCounts( matrix )
%COMPUTEZEROROWCOUNTS Summary of this function goes here
%   Detailed explanation goes here

count = 0;

for rowIdx = 1:size(matrix, 1)
    if sum(abs(matrix(rowIdx, :))) < 1e-6
        count = count + 1;
        matrix(rowIdx, :) = 0;
    end
    
end

end

