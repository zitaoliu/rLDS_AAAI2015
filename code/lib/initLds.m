function [ init ] = initLds( hiddenStateSize, observSize, initLdsOption, seed )
%INITLDS Summary of this function goes here
%   Detailed explanation goes here

if initLdsOption == 1
   
    rng(seed)
    init.A = rand(hiddenStateSize, hiddenStateSize);
    init.C = rand(observSize, hiddenStateSize);
    init.Q = 0.01*eye(hiddenStateSize);
    init.R = 0.01*eye(observSize);
    
    init.initx = rand(hiddenStateSize, 1);
    init.initV = 0.01*eye(hiddenStateSize);
end


end

