%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEMO - GROUPLASSO
% Plot path for a MKL example
%
% the group Lasso is used after an incomplete cholesky decomposition is
% performed. See e.g.
% Francis R. Bach, Michael I. Jordan. Predictive
% low-rank decomposition for kernel methods. 
% Proceedings of the Twenty-second International Conference on Machine Learning (ICML), 2005.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath utilities
clear all;
seed=1;
randn('state',seed);
rand('state',seed);

m = 4;      % number of groups
n = 1000;   % number of data points
r = 2;      % number of selected groups
noise_std = 2e-1;
kmax = 5;
n0_coeff=1;
alpha_kernel = 2;

SigmaG_sqrt = randn(m,m);
SigmaG = SigmaG_sqrt' * SigmaG_sqrt;

% normalize to unit trace
diagonal = zeros(m,1);
for i=1:m
    trK(i) = SigmaG(i,i) ;
    diagonal(i) = trK(i);
end
SigmaG = diag( 1./diagonal.^.5) * SigmaG * diag( 1./diagonal.^.5);
SigmaG_sqrt =   SigmaG_sqrt * diag( 1./diagonal.^.5);

% generate random functions
% generate W0c -> random coefficients on the eigenfunctions
W0c = zeros((kmax+1)*m,1);
nmax_coeff = kmax+1;
coeff_eig = zeros(nmax_coeff,r);
coeff0_eig = zeros(nmax_coeff,r);
for j=1:r
    if 0
        coeff_eig(:,j) =  (randn(nmax_coeff,1))  .* (1./(1:nmax_coeff)'.^1);
        coeff_cos(1,j) = 0;
        coeff_sin(1,j) = 0;
    else
        coeff_eig(:,j) = (randn(nmax_coeff,1)) .* ( (1:nmax_coeff)' ./ ( n0_coeff + (1:nmax_coeff).^2 )');
    end
    normj = norm( coeff_eig(:,j) );
    coeff0_eig(:,j) = coeff_eig(:,j) / normj;   % normalized to unit norm
    normj = normj / (  ( 1 + 2 * rand) / 3 );
    coeff_eig(:,j) = coeff_eig(:,j) / normj;    % normalized to norm between 1/3 and 1
    W0c((j-1)*(kmax+1)+1:j*(kmax+1)) =  coeff_eig(:,j);
end

[conditions,SXX] = compute_consistency_condition_MKL(kmax,m,SigmaG,alpha_kernel,W0c,1:r);
% consistency condition (see JMLR paper)
conditions

% now generates data
a = 1./4./diag(SigmaG);
b = ones(m,1) * alpha_kernel;
c = ( a.^2 + 2 * a .* b).^.5;

XX = randn(n,m) * SigmaG_sqrt;  % THIS IS THE TRUE DATA
X = zeros(n,m*(kmax+1));        % THIS IS THE DATA PROJECTED ON EIGENFUNCIONS
for j=1:m
    for k=0:kmax
        X(:,k+1+(j-1)*(kmax+1)) = build_eigenfunction(a(j),b(j),k,XX(:,j));
    end
end
Y =  X * W0c + randn(n,1) * noise_std * sqrt(W0c' * SXX * W0c);


% compute least-square estimate and check norm
lambda = 1e-4;
W_LS = (X'*X/n + eye(m*(kmax+1)) * lambda) \ (X'*Y/n);
groups = cell(1,m);
for i=1:m
    groups{i}= (i-1)*(kmax+1)+1:i*(kmax+1);
end
% check that the least square estimate is not too bad...
W2norms(W_LS-W0c,groups)
groups_c = groups;

% NOW BUILD DATA FOR MKL !!!
% PERFORMS INCOMPLETE CHOLESKY TO GO BACK TO INPUT SPACE

lambda = 1e-12;
mmax = 100;
x = [];
groups = cell(1,m);
for j=1:m

    [G,P,mm] = icd_gauss(XX(:,j)',b(j),n*lambda*1e-3,min(n,mmax)); % perform largest ICD to get lambda
    I = P(1:mm);
    [temp,Pi] = sort(P);
    G = G(Pi,1:mm);
    groups{j} = size(x,2)+1:size(x,2)+size(G,2);
    x = [ x , G ];
end



% center
x = x - repmat(mean(x,1),size(x,1),1);
Y = Y - repmat(mean(Y,1),size(Y,1),1);
Sigma = x' * x / n;
Q = x' * Y / n;
W_LS = ( Sigma + 1e-4 * eye(length(Sigma))) \ Q;

alpha = 1e-4;
lambdamax = 0;
lambdamax = max( 2 * W2norms(Q,groups) );
lambdamin = 1e-5 * sum( Q(:) .* W_LS(:) ) / sum( W2norms(W_LS,groups)  );
nlambdas = 200;
lambdas = exp( log(lambdamax) + ( log(lambdamin) - log(lambdamax) ) / nlambdas * ( 0:(nlambdas-1) ));

optparam.tol = 1e-16;
optparam.kmax = 200;
optparam.tol_df = 1e-20;
optparam.display= 0;
W = zeros(size(Sigma,1),1);
svds = [];
for ilambda = 1:length(lambdas)
    lambda = lambdas(ilambda);
    W = newton_grouplasso(W,Sigma,Q,groups,lambda,alpha);
    norms(:,ilambda) = W2norms(W,groups);

end


% normalized to sum of eta = 1
% go back to mus
mus = lambdas ./ sum( norms,1 );
norms = norms ./ repmat(sum(norms,1),m,1);

minlambda = max(find(max(norms)>.999));
[a,b] = max(norms(:,minlambda-1));
norms(:,1:minlambda-1)=0;
norms(b,1:minlambda-1)=1;

% add only part of what's before the beginning of the path
minlambda = max(1,minlambda-3);
maxlambda= max( find( abs(sum( abs(norms(:,2:end)-norms(:,1:end-1)) )) > 1e-4 ));

plot(-log(mus(minlambda:maxlambda)),norms(:,minlambda:maxlambda)','linewidth',2); hold on;
plot(-log(mus(minlambda:maxlambda)),repmat(W2norms(W0c,groups_c)/ sum(W2norms(W0c,groups_c) ) ,1,length(lambdas(minlambda:maxlambda)))',':','linewidth',2); hold off
axis( [ min(-log(mus(minlambda:maxlambda)) ) max(-log(mus(minlambda:maxlambda)) ) 0 max(norms(:))*1.1] );
xlabel('-log(\mu)','fontsize',16);
ylabel('\eta','fontsize',16);
title('group norms of MKL estimate','fontsize',16)

