function [bestA, bestC, bestQ, bestR, bestInitx, bestInitV, LL] = ...
    learn_kalman_sparse(data, A, C, Q, R, initx, initV, lambdaA,...
    lambda2A, barrier, max_iter, diagQ, diagR, ARmode, constr_fun, varargin)

% LEARN_KALMAN_SPARSE Find the ML parameters of a stochastic Linear Dynamical System using EM with Group Lasso.
%
% [bestA, bestC, bestQ, bestR, bestInitx, bestInitV, LL] = learn_kalman_sparse(data, A, C, Q, R, initx, initV, lambdaA, lambda2A, barrier) fits
% the parameters which are defined as follows
% 
%   x(t+1) = bestA*x(t) + w(t),  w ~ N(0, bestQ),  x(0) ~ N(bestInitx, bestInitV)
%   y(t)   = bestC*x(t) + v(t),  v ~ N(0, bestR)
% 
% A is the initial value, bestA is the final value, etc. The same for C, Q, R, initx, initV.
%
% lambdaA is the penalty parameter for group lasso
% lambda2A is the penalty parameter for frobenius norm
% barrier is the barrier value for the SOCP
%
% DATA(:,t,l) is the observation vector at time t for sequence l. If the sequences are of
% different lengths, you can pass in a cell array, so DATA{l} is an O*T matrix.
% LL is the "learning curve": a vector of the log lik. values at each iteration.
% LL might go positive, since prob. densities can exceed 1, although this probably
% indicates that something has gone wrong e.g., a variance has collapsed to 0.
%
% There are several optional arguments, that should be passed in the following order.
% learn_kalman_sparse(data, A, C, Q, R, initx, initV, lambdaA, lambda2A, barrier, max_iter, diagQ, diagR, ARmode, constr_fun, varargin)
% MAX_ITER specifies the maximum number of EM iterations (default 10).
% DIAGQ=1 specifies that the Q matrix should be diagonal. (Default 0).
% DIAGR=1 specifies that the R matrix should also be diagonal. (Default 0).
% ARMODE=1 specifies that C=I, R=0. i.e., a Gauss-Markov process. (Default 0).
%
% For details about learning LDS with regularization, see
% - Liu, Z., & Hauskrecht, M..  A Regularized Linear Dynamical System Framework for Multivariate Time Series Analysis. 
% - In Proceedings of the AAAI Conference on Artificial Intelligence. 2015.
% - download link: http://www.zitaoliu.com/download/aaai2015.pdf
%
% The group lasso solver is downloaded from http://www.di.ens.fr/~fbach/grouplasso/
%
% The majority pieces of this code are based the Kalman Filter Matlab
% package(http://www.cs.ubc.ca/~murphyk/Software/Kalman/kalman.html)
%
% For details, see
% - Ghahramani and Hinton, "Parameter Estimation for LDS", U. Toronto tech. report, 1996
% - Digalakis, Rohlicek and Ostendorf, "ML Estimation of a stochastic linear system with the EM
%      algorithm and its application to speech recognition",
%       IEEE Trans. Speech and Audio Proc., 1(4):431--442, 1993.


if nargin < 11, max_iter = 10; end
if nargin < 12, diagQ = 0; end
if nargin < 13, diagR = 0; end
if nargin < 14, ARmode = 0; end
if nargin < 15, constr_fun = []; end

thresh = 1e-4;

if ~iscell(data)
    N = size(data, 3);
    data = num2cell(data, [1 2]); % each elt of the 3rd dim gets its own cell
else
    N = length(data);
end

N = length(data);
ss = size(A, 1);
os = size(C,1);

alpha = zeros(os, os);
Tsum = 0;
for ex = 1:N
    %y = data(:,:,ex);
    y = data{ex};
    %T = length(y);
    T =size(y,2);
    Tsum = Tsum + T;
    alpha_temp = zeros(os, os);
    for t=1:T
        alpha_temp = alpha_temp + y(:,t)*y(:,t)';
    end
    alpha = alpha + alpha_temp;
end

previous_loglik = -inf;
loglik = 0;
converged = 0;
num_iter = 1;
LL = [];

% Convert to inline function as needed.
if ~isempty(constr_fun)
    constr_fun = fcnchk(constr_fun,length(varargin));
end

groupsA = cell(ss,1);
for i=1:ss
    groupsA{i} = (i-1)*ss + 1:i*ss;
end

maxLoglik = -Inf;
bestA = A;
bestC = C;
bestQ = Q;
bestR = R;
bestInitx = initx;
bestInitV = initV;

while ~converged && (num_iter <= max_iter)
    
    %%% E step
    
    delta = zeros(os, ss);
    gamma = zeros(ss, ss);
    gamma1 = zeros(ss, ss);
    gamma2 = zeros(ss, ss);
    beta = zeros(ss, ss);
    P1sum = zeros(ss, ss);
    x1sum = zeros(ss, 1);
    loglik = 0;
    
    for ex = 1:N
        y = data{ex};
        %T = length(y);
        [beta_t, gamma_t, delta_t, gamma1_t, gamma2_t, x1, V1, loglik_t] = ...
            Estep(y, A, C, Q, R, initx, initV, ARmode);
        beta = beta + beta_t;
        gamma = gamma + gamma_t;
        delta = delta + delta_t;
        gamma1 = gamma1 + gamma1_t;
        gamma2 = gamma2 + gamma2_t;
        P1sum = P1sum + V1 + x1*x1';
        x1sum = x1sum + x1;
        loglik = loglik + loglik_t;
    end
    LL = [LL loglik];
    
    num_iter =  num_iter + 1;
    
    if loglik > maxLoglik 
        bestA = A;
        bestC = C;
        bestQ = Q;
        bestR = R;
        bestInitx = initx;
        bestInitV = initV;
        
        maxLoglik = loglik;
    end
    
    
    [gamma1, gamma1EigVecs, gamma1EigVals] = psdProjection(gamma1);
    [~, gammaEigVecs, gammaEigVals] = psdProjection(gamma);
    [gamma2, ~, ~] = psdProjection(gamma2);
    [alpha, ~, ~] = psdProjection(alpha);
    [P1sum, ~, ~] = psdProjection(P1sum);
    
    %%% M step
    
    Tsum1 = Tsum - N;
    A = beta * psdInverse(gamma1EigVecs, gamma1EigVals);

    Q = (gamma2 - A*beta') / Tsum1;
    [Q, QEigVecs, QEigVals] = psdProjection( Q );
    if diagQ
        Q = diag(diag(Q));
    end
    
    
    % group lasso on A here
    if diagQ
        L = sqrt(diag(1./diag(Q)));
        detQ = det(Q);
    else
        
        %% compute the approximated decomposition of inv(Q)
        QEigValsVec = diag(QEigVals);
        qIdx = QEigValsVec > 1e-4;

        detQ = prod(QEigValsVec(qIdx));
        L = QEigVecs(:, qIdx) * sqrt(diag(1./QEigValsVec(qIdx)));

    end
    
    Ldim = size(L, 2);
    
    LI = kron(L', eye(ss));
    H = LI' * kron(eye(Ldim), gamma1) * LI;
    H = H + lambda2A .* eye(size(H));
    b = (L(:)'* kron(eye(Ldim), beta) * LI)';
    
    a = newton_grouplasso(A(:), H, b, groupsA, lambdaA, barrier);
    A = reshape(a, ss, ss)';
    
    [numZeroRows, A] = computeZeroRowCounts(A);
    
    if numZeroRows == size(A, 1)
        break;
    end
    
    if ~ARmode
        C = delta * psdInverse(gammaEigVecs, gammaEigVals);
        R = (alpha - C*delta') / Tsum;
        [R, ~, ~] = psdProjection( R );
        if diagR
            R = diag(diag(R));
        end
        
    end
    
    initx = x1sum / N;
    initV = P1sum/N - initx*initx';
    
    if ~isempty(constr_fun)
        [A,C,Q,R,initx,initV] = feval(constr_fun, A, C, Q, R, initx, initV, varargin{:});
    end
    
    converged = em_converged(loglik, previous_loglik, thresh);
    previous_loglik = loglik;
end



%%%%%%%%%

function [beta, gamma, delta, gamma1, gamma2, x1, V1, loglik] = ...
    Estep(y, A, C, Q, R, initx, initV, ARmode)
%
% Compute the (expected) sufficient statistics for a single Kalman filter sequence.
%

[os, T] = size(y);
ss = length(A);

if ARmode
    xsmooth = y;
    Vsmooth = zeros(ss, ss, T); % no uncertainty about the hidden states
    VVsmooth = zeros(ss, ss, T);
    loglik = 0;
else
    [xsmooth, Vsmooth, VVsmooth, loglik] = kalman_smoother(y, A, C, Q, R, initx, initV);
end

delta = zeros(os, ss);
gamma = zeros(ss, ss);
beta = zeros(ss, ss);
for t=1:T
    delta = delta + y(:,t)*xsmooth(:,t)';
    gamma = gamma + xsmooth(:,t)*xsmooth(:,t)' + Vsmooth(:,:,t);
    if t>1 beta = beta + xsmooth(:,t)*xsmooth(:,t-1)' + VVsmooth(:,:,t); end
end
gamma1 = gamma - xsmooth(:,T)*xsmooth(:,T)' - Vsmooth(:,:,T);
gamma2 = gamma - xsmooth(:,1)*xsmooth(:,1)' - Vsmooth(:,:,1);

x1 = xsmooth(:,1);
V1 = Vsmooth(:,:,1);



