%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEMO - GROUPLASSO
% Estimate probability of correct pattern selection and RMS
% in two cases (inconsistent and consistent)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
addpath utilities

seed= 8;                    % "well chosen" random seed to get consistency
randn('state',seed);
rand('state',seed);

d = 2;      % dimension of groups
m = 4;      % number of groups
n = 200;    % number of data points
r = 2;      % number of selected groups
noise_std = 2e-1;   % noise standard deviation

nrep = 5;
n = 1000;
nlambdas = 200;


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

% generate sparse weigts
W0 = randn(dm,1);
for i=1:r, W0(groups{i}) = W0(groups{i}) / norm( W0(groups{i}) ) * ( 1 + 2 * rand) / 3; end

for i=r+1:m, W0(groups{i}) = 0; end
J = [];
for i=1:r, J = [ J, groups{i} ]; end
Jc = [];
for i=r+1:m, Jc = [ Jc, groups{i} ]; end

% check consistency condition
invSigmaGJJ = inv( SigmaG(J,J) );
W0norm = W0;
for i=1:r, 
    W0norm(groups{i}) = W0norm(groups{i}) / norm( W0norm(groups{i}) );
end
for i=r+1:m
   conditions(i) = norm( SigmaG(groups{i},J) * invSigmaGJJ * W0norm(J) );
end
conditions = conditions(r+1:m)
% if all less than one, group lasso is consisten
% if not,



% determine start and end of path
Q0 = SigmaG * W0;
lambdamax = 0;
lambdamax = max( 2 * W2norms(Q0,groups) );
lambdamin = 1e-5 * sum( Q0(:) .* W0(:) ) / sum( W2norms(W0,groups)  );
lambdas = exp( log(lambdamax) + ( log(lambdamin) - log(lambdamax) ) / nlambdas * ( 0:(nlambdas-1) ));



fprintf('%d replications \n',nrep);

for irep=1:nrep
    irep
    
    X = randn(n,d*m) * SigmaG_sqrt;
    Y =  X * W0 + randn(n,1) * noise_std * sqrt(W0' * SigmaG * W0);
    % center
    X = X - repmat(mean(X,1),size(X,1),1);
    Y = Y - repmat(mean(Y,1),size(Y,1),1);
    Sigma = X' * X / n;
    Q = X' * Y / n;


    alpha = 1e-4;   % barrier parameter

    optparam.tol = 1e-16;
    optparam.kmax = 100;
    optparam.tol_df = 1e-20;
    optparam.display= 0;
    W = zeros(dm,1);

    for ilambda = 1:length(lambdas)
        lambda = lambdas(ilambda);

        W = newton_grouplasso(W,Sigma,Q,groups,lambda,alpha);
        correct(ilambda,irep) = (...
            sum( abs( ( W2norms(W,groups) >1e-6 ) -   [ones(r,1); zeros(m-r,1)] ), 1 ) == 0 );
        dists(ilambda,irep) = norm(W-W0,'fro')^2;
    end

end


subplot(2,2,1);
plot(-log10(lambdas),squeeze(mean(correct,2))','linewidth',2)

xlabel('-log(\lambda)','fontsize',16);
ylabel('P(correct pattern)','fontsize',16);
title('Consistent','fontsize',16);
axis([ min(-log10(lambdas)), max(-log10(lambdas)), 0 , 1 ]);


subplot(2,2,2);
plot(-log10(lambdas),log10(mean(dists,2))','linewidth',2);
xlabel('-log(\lambda)','fontsize',16);
ylabel('log(RMS)','fontsize',16);
title('Consistent','fontsize',16);





%%%%%%%%%%%%%%%%
% INCONSISTENT %
%%%%%%%%%%%%%%%%


seed= 4;      % "well chosen" random seed to get inconsistency
randn('state',seed);
rand('state',seed);

d = 2;      % dimension of groups
m = 4;      % number of groups
n = 200;    % number of data points
r = 2;      % number of selected groups
noise_std = 2e-1;   % noise standard deviation



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

% generate sparse weigts
W0 = randn(dm,1);
for i=1:r, W0(groups{i}) = W0(groups{i}) / norm( W0(groups{i}) ) * ( 1 + 2 * rand) / 3; end

for i=r+1:m, W0(groups{i}) = 0; end
J = [];
for i=1:r, J = [ J, groups{i} ]; end
Jc = [];
for i=r+1:m, Jc = [ Jc, groups{i} ]; end

% check consistency condition
invSigmaGJJ = inv( SigmaG(J,J) );
W0norm = W0;
for i=1:r, 
    W0norm(groups{i}) = W0norm(groups{i}) / norm( W0norm(groups{i}) );
end
for i=r+1:m
   conditions(i) = norm( SigmaG(groups{i},J) * invSigmaGJJ * W0norm(J) );
end
conditions = conditions(r+1:m)
% if all less than one, group lasso is consisten
% if not,



% determine start and end of path
Q0 = SigmaG * W0;
lambdamax = 0;
lambdamax = max( 2 * W2norms(Q0,groups) );
lambdamin = 1e-5 * sum( Q0(:) .* W0(:) ) / sum( W2norms(W0,groups)  );
lambdas = exp( log(lambdamax) + ( log(lambdamin) - log(lambdamax) ) / nlambdas * ( 0:(nlambdas-1) ));
save lambdas_general_consistent lambdas



fprintf('%d replications \n',nrep);
for irep=1:nrep
    irep
    X = randn(n,d*m) * SigmaG_sqrt;
    Y =  X * W0 + randn(n,1) * noise_std * sqrt(W0' * SigmaG * W0);
    % center
    X = X - repmat(mean(X,1),size(X,1),1);
    Y = Y - repmat(mean(Y,1),size(Y,1),1);
    Sigma = X' * X / n;
    Q = X' * Y / n;


    alpha = 1e-4;

    optparam.tol = 1e-16;
    optparam.kmax = 100;
    optparam.tol_df = 1e-20;
    optparam.display= 0;
    W = zeros(dm,1);

    for ilambda = 1:length(lambdas)
        lambda = lambdas(ilambda);

        barriers = 10.^[0:.5:4];
        for ibar = 1:length(barriers)
            [W,exit_parameters] = minimize_newton(W,@grouplasso,optparam,Sigma,Q,groups,lambda,alpha/barriers(ibar));
        end

        correct(ilambda,irep) = (...
            sum( abs( ( W2norms(W,groups) >1e-6 ) -   [ones(r,1); zeros(m-r,1)] ), 1 ) == 0 );
        dists(ilambda,irep) = norm(W-W0,'fro')^2;

    end

end


subplot(2,2,3);
plot(-log10(lambdas),squeeze(mean(correct,2))','linewidth',2)

xlabel('-log(\lambda)','fontsize',16);
ylabel('P(correct pattern)','fontsize',16);
title('Inconsistent','fontsize',16);
axis([ min(-log10(lambdas)), max(-log10(lambdas)), 0, 1 ]);

subplot(2,2,4);
plot(-log10(lambdas),log10(mean(dists,2))','linewidth',2);
xlabel('-log(\lambda)','fontsize',16);
ylabel('log(RMS)','fontsize',16);
title('Inconsistent','fontsize',16);


