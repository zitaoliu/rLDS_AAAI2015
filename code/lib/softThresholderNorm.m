function [ newX ] = softThresholderNorm( lambda, x )
%SOFTTHRESHOLDERMC Summary of this function goes here
%   Detailed explanation goes here


[u, sigma, v] = svd(x);


tmp = diag(sigma) - lambda;

tmp(tmp < 0) = 0;

newX = u*diag(tmp)*v';



end

