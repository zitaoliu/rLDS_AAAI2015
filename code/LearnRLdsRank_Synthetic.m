clc;
clear;

addAllToPath();

config = initSynthetic();
observSize = size(config.dataMat, 1);

lambda2As = [1000 1000 1000];
lambda3As = [9500 30000 30000];

trainData = config.dataMat(:, config.trainIdx);
hankelSize = chooseHankelSize( config.hankelSizeOption, size(trainData, 2) );

models = cell(length(config.hiddenStates), 1);

for stateIdx = 1:length(config.hiddenStates)
    
    hiddenStateSize = config.hiddenStates(stateIdx);
    init = learn_spectral(hiddenStateSize, hankelSize, trainData);
    
    lambda2A = lambda2As(stateIdx);
    lambda3A = lambda3As(stateIdx);
    
    [models{stateIdx, 1}.A, models{stateIdx, 1}.C, models{stateIdx, 1}.Q, ...
        models{stateIdx, 1}.R, models{stateIdx, 1}.initx, models{stateIdx, 1}.initV, ~] ...
        = learn_kalman_rank(trainData, init.A, init.C, init.Q, init.R, init.initx, ...
        init.initV, lambda3A, lambda2A, config.sdgMaxIter, config.maxIter);
   
end


