%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEMO - GROUPLASSO
% Plot path for a group lasso example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
addpath utilities


seed=8;
randn('state',seed);
rand('state',seed);

d = 2;      % dimension of groups
m = 4;      % number of groups
n = 200;   % number of data points
r = 2;      % number of selected groups
noise_std = 2e-1;

% generate random Gaussians with fixed group sizes
dm = d * m;
groups = cell(d,1);
for i=1:m, groups{i} = (i-1)*d+1:i*d; end

SigmaG_sqrt = randn(dm,dm);
SigmaG = SigmaG_sqrt' * SigmaG_sqrt;

% normalize to unit trace
diagonal = zeros(dm,1);
for i=1:m
    trK(i) = trace( SigmaG(groups{i},groups{i}) );
    diagonal(groups{i}) = trK(i);
end
SigmaG = diag( 1./diagonal.^.5) * SigmaG * diag( 1./diagonal.^.5);
SigmaG_sqrt =   SigmaG_sqrt * diag( 1./diagonal.^.5);



% generate sparse weights
W0 = randn(dm,1);
for i=1:r, W0(groups{i}) = W0(groups{i}) / norm( W0(groups{i}) ) * ( 1 + 2 * rand) / 3; end
for i=r+1:m, W0(groups{i}) = 0; end
J = [];
for i=1:r, J = [ J, groups{i} ]; end
Jc = [];
for i=r+1:m, Jc = [ Jc, groups{i} ]; end

X = randn(n,d*m) * SigmaG_sqrt;
Y =  X * W0 + randn(n,1) * noise_std * sqrt(W0' * SigmaG * W0);

% center
X = X - repmat(mean(X,1),size(X,1),1);
Y = Y - repmat(mean(Y,1),size(Y,1),1);

Sigma = X' * X / n;
Q = X' * Y / n;

% compute least-square estimate
W_LS = Sigma \ Q;


% check consistency of estimates
norm(W0-W_LS)/norm(W0)
norm(SigmaG-Sigma,'fro') / norm(SigmaG,'fro')

% check condition for consistency
invSigmaGJJ = inv( SigmaG(J,J) );
W0norm = W0;
for i=1:r,
    W0norm(groups{i}) = W0norm(groups{i}) / norm( W0norm(groups{i}) );
end
for i=r+1:m
    conditions(i) = norm( SigmaG(groups{i},J) * invSigmaGJJ * W0norm(J) );
end
conditions(r+1:m)


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
W = zeros(dm,1);
svds = [];
for ilambda = 1:length(lambdas)
    lambda = lambdas(ilambda);
    W = newton_grouplasso(W,Sigma,Q,groups,lambda,alpha);

    W
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

% add only parto of what's before the beginning of the path (plus joli)
minlambda = max(1,minlambda-4);
maxlambda= max( find( abs(sum( abs(norms(:,2:end)-norms(:,1:end-1)) )) > 1e-4 ));


plot(-log(mus(minlambda:maxlambda)),norms(:,minlambda:maxlambda)','linewidth',2); hold on;
plot(-log(mus(minlambda:maxlambda)),repmat(W2norms(W0,groups)/ sum( W2norms(W0,groups)) ,1,length(lambdas(minlambda:maxlambda)))',':','linewidth',2); hold off
axis( [ min(-log(mus(minlambda:maxlambda)) ) max(-log(mus(minlambda:maxlambda)) ) 0 max(norms(:))*1.1] );
xlabel('-log(\mu)','fontsize',16);
ylabel('\eta','fontsize',16);
title('group norms of group lasso estimate','fontsize',16)

