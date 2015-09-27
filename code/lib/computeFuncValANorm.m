function [ fval ] = computeFuncValANorm( A, invQ, lambda, gamma1, beta, lambda2A )
%COMPUTEFUNCVALANORM Summary of this function goes here
%   Detailed explanation goes here

fval = 0.5 * trace(gamma1*A'*invQ*A - 2*A'*invQ*beta);

[~, s, ~] = svd(A);

fval = fval + lambda * sum(diag(s));

fval = fval + lambda2A * norm(A, 'fro')^2;

end

