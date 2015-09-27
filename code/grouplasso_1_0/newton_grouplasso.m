function W = newton_grouplasso(W,Sigma,Q,groups,lambda,alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% performs estimation of the grouplasso for parameter lambda
% the loss function which is used is .5 * W' * Sigma * W - Q' * W
%
% INPUT
% Sigma : quadratic form (typically covariance matrix of covariates)
% Q     : vector (linear form, typically covariance between covariates and
%         responses)
% groups: cell array of group indices
% lambda: regularization parameters
% alpha : barrier parameter = 1e-4
%
% see demo scripts to see example of use
%

if isempty(W), W = zeros(length(Q),1); end


% very precise: goal is not speed but precision
% optparam.tol = 1e-16;
optparam.tol = 1e-6;
optparam.kmax = 100;
%optparam.tol_df = 1e-20;
optparam.tol_df = 1e-6;
optparam.display= 0;



barriers = 10.^[0:.5:4];
for ibar = 1:length(barriers)
    [W,exit_parameters] = minimize_newton(W,@grouplasso,optparam,Sigma,Q,groups,lambda,alpha/barriers(ibar));
end