function [ minA, fvals ] = istaANorm( A, sdgMaxIter, gamma1, beta, ...
    stepSize, threshold, lambda, invQ, lambda2A )
%ISTANORM Summary of this function goes here
%   Detailed explanation goes here

fold = computeFuncValANorm(A, invQ, lambda, gamma1, beta, lambda2A);

fvals = [fold];

minFuncVal = realmax;
minA = 0;

for i = 1:sdgMaxIter
    
    gradient = invQ*A*gamma1 - invQ*beta + lambda2A * A;
    tmp = A - stepSize .* gradient;
    
    if sum(sum(isnan(tmp))) > 0 || sum(sum(isinf(tmp))) > 0
        break;
    end
    
    A = softThresholderNorm(lambda*stepSize, tmp);
    
    fnew = computeFuncValANorm(A, invQ, lambda, gamma1, beta, lambda2A);
    fvals = [fvals fnew];
    
    if minFuncVal > fnew
        minA = A;
        minFuncVal = fnew;
    end
    
    if abs(fnew-fold)*2/(abs(fold)+abs(fnew)+eps) < threshold
        break;
    end
    
    fold = fnew;
end



end

