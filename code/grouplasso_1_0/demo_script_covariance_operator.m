%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEMO - covariance operator      
% (1) compute covariance operator on eigenfunctions of the corresponding
% marginal covariance operators for jointly Gaussian variables
% (2) plot eigenfunctions for marginal covariance operator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath utilities
clear all;
seed=1;
randn('state',seed);
rand('state',seed);

m = 3;                  % number of variables
alpha_kernel = 2;       % kernel width
kmax = 30;              % number of eigenfunctions to keep

% generate random covariance matrix with m variables
SigmaG_sqrt = randn(m,m);
SigmaG = SigmaG_sqrt' * SigmaG_sqrt;





% compute all factorials and store them
global factorials;
if isempty(factorials)
    factorials(1)=1;
    factorials(2)=1;
    for k=3:100
        factorials(k) = factorials(k-1) * (k-1);
    end
end


% vectors of parameters
a = 1./4./diag(SigmaG);
b = ones(m,1) * alpha_kernel;
c = ( a.^2 + 2 * a .* b).^.5;
A = a + b + c;
B = b./A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) CORRELATION OPERATOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute correlation between k-th and l-th eigenfunctions for i, j;
C = eye((kmax+1)*m);
for i=1:m
    for j=1:m
        if (i~=j)

            for l=0:kmax
                for k=0:kmax
                    C(k+1 + (i-1)*(kmax+1),l+1 + (j-1)*(kmax+1)) = ( 2^(k+l) * factorials(k+1) * factorials(l+1) * sqrt( c(i) * c(j) / a(i) / a(j) ) )^(-.5) ...
                        / 2 / pi / sqrt( SigmaG(i,i) * SigmaG(j,j) - SigmaG(i,j) * SigmaG(i,j) ) / a(i)^.5 / a(j).^.5 / 2 ...
                        * compute_DklQ(k,l,[.5-.5*a(i)/c(i),0;0,.5-.5*a(j)/c(j)] + 1/4 * ...
                        inv( [ SigmaG(i,i)*c(i), SigmaG(i,j) * c(i)^.5 * c(j)^.5 ; ...
                        SigmaG(i,j) * c(i)^.5 * c(j)^.5 , SigmaG(j,j)*c(j) ] ) );

                end
            end
        end
    end
end


% compute means of each eigenfunction
mu = zeros((kmax+1)*m,1);
for i=1:m
    for l=0:kmax
        if mod(l,2)==0
            k=l/2;
            mu(l+1 + (i-1)*(kmax+1)) =  c(i)^.25 / a(i)^.25 / 2^k / factorials(k+1) * factorials(2*k+1)^.5 *(-1)^k  *...
                ( (a(i)-c(i))/(a(i)+c(i)) )^k * sqrt( 2 * a(i) / ( c(i) + a(i) ) );
        end
    end
end

% compute square root matrices (to compute correlation operator)
sqrtM = eye((kmax+1)*m);
for i=1:m
    mi = mu((i-1)*(kmax+1)+1:i*(kmax+1));
    sqrti = eye(kmax+1)- mi*mi' / (mi'*mi) + mi*mi' * ( 1 - mi'*mi)^-.5 / (mi'*mi);
    sqrtM((i-1)*(kmax+1)+1:i*(kmax+1),(i-1)*(kmax+1)+1:i*(kmax+1))=sqrti;
end

lambdas = zeros((kmax+1)*m,1);  % eigenvalues of the marginal covariance operators
for i=1:m
    lambdas((i-1)*(kmax+1)+1:i*(kmax+1)) = (2*a(i)/A(i))^.5 * B(i).^(0:kmax);
end


CXX = sqrtM * (C - mu*mu')*sqrtM;                           % correlation operator on the eigenfunction basis
SXX = diag(lambdas.^.5) * (C - mu*mu') * diag(lambdas.^.5); % covariance operator on the eigenfunction basis
subplot(1,2,1);

imagesc(CXX)




% compare correlation of operators and original correlation: they are
% equal (conjecture!!)

% normalize to unit trace
diagonal = zeros(m,1);
for i=1:m
    trK(i) = SigmaG(i,i) ;
    diagonal(i) = trK(i);
end
CorrG = diag( 1./diagonal.^.5) * SigmaG * diag( 1./diagonal.^.5);

min(eig(CorrG))
min(eig(CXX))




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) plot eigenfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2);
j = 1;
X = -5:.01:5;       % abcissas for function evaluations
for k=0:5
EIGENFUNCTIONS(k+1,:) = build_eigenfunction(a(j),b(j),k,X);
end
plot(X,EIGENFUNCTIONS')


