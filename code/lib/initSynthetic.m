function [ config ] = initSynthetic(  )

%% data processing
load('../data/synthetic.mat');
dataMat = data;

config.rawMat = dataMat;
config.T = size(dataMat, 2);
config.numTs = size(dataMat, 1);
config.dataMeanVec = mean(dataMat, 2);
config.dataMat = config.rawMat - repmat(config.dataMeanVec, 1, config.T);

%% initalization config
config.hankelSizeOption = 3;
config.initLdsOption = 1;

%% EM configuration
config.maxIter = 200;

config.trainIdx = 1:200;
config.hiddenStates = [15 20 30];

config.sdgMaxIter = 500;
config.barrier = 1e-4;

% %% Sparse Group Lasso Lds config
% config.lambdaAs = 10.^[3:6];
% config.lambda2As = 10.^[3];
% config.lambda3As = [linspace(6000, 10000, 50) linspace(10^4, 10^6, 50)];
% config.valiRatio = 0.8;
% 
% 
% config.sdgThreshold = 1e-4;

end

