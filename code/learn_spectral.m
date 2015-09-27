function [ model ] = learn_spectral( hiddenState, hankelSize, data )
%LEARN_SPECTRAL_L1 Summary of this function goes here
%   Detailed explanation goes here

Y = data;

d = hankelSize;

m   = size(Y,1);
t   = size(Y,2);

% Subtract Mean from Observations
Ymean = mean(Y,2);
Y = Y - Ymean*ones(1,t);

% Build Hankel Matrix
tau = t - d+1;
D = zeros(d*m, tau);

for i = 1:d
    for j = 1:tau
        D((i-1)*m + 1:i*m,j) = Y(:,(j-1) + i);
    end
end

D = D';

% Perform SVD on Hankel Matrix
if size(D,2) < size(D,1)
    [V,S,U] = svd(D,0);
else
    [U,S,V] = svd(D',0);
end

n = hiddenState;

if n <= size(S, 1)
    
    V = V(:,1:n);
    S = S(1:n,1:n);
    U = U(:,1:n);
    
    Xhat = S*V';
    Chat = U(1:m,:);
    
    % Estimate Dynamics Matrix
    S1 = Xhat(:,1:end-1);
    S2 = Xhat(:,2:end);
    
    Ahat = S2*pinv(S1);
    
    % Estimate Evolution and Observation Noise Covariance Matrices
    What = S2-Ahat*S1;
    Qhat = cov(What');
    
    if m <= 100
        Vhat = Y(:,1:size(Xhat, 2)) - Chat*Xhat;
        Rhat = cov(Vhat');
    else
        Rhat = [];
    end
    
    model.initx = Xhat(:,1);
    model.initV = 0.01* eye(n, n);
    model.A = Ahat;
    model.C = Chat;
    model.R = Rhat;
    model.Q = Qhat;
    
else
    model = cell(0);
end

end

