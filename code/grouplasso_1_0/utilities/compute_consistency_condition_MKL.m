function [conditions,SXX] = compute_consistency_condition_MKL(kmax,m,SigmaG,alpha_kernel,W0,indJ)

global factorials;
if isempty(factorials)
    factorials(1)=1;
factorials(2)=1;
for k=3:100
    factorials(k) = factorials(k-1) * (k-1);
end
end

    a = 1./4./diag(SigmaG);
    b = ones(m,1) * alpha_kernel;
    c = ( a.^2 + 2 * a .* b).^.5;
    A = a + b + c;
    B = b./A;


    % compute correlation between k-th and l-th eigenfunctions for i=1, j =2;
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

    % compute means
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

    sqrtM = eye((kmax+1)*m);
    for i=1:m
        mi = mu((i-1)*(kmax+1)+1:i*(kmax+1));
        sqrti = eye(kmax+1)- mi*mi' / (mi'*mi) + mi*mi' * ( 1 - mi'*mi)^-.5 / (mi'*mi);
        sqrtM((i-1)*(kmax+1)+1:i*(kmax+1),(i-1)*(kmax+1)+1:i*(kmax+1))=sqrti;
    end

    lambdas = zeros((kmax+1)*m,1);
    for i=1:m
        lambdas((i-1)*(kmax+1)+1:i*(kmax+1)) = (2*a(i)/A(i))^.5 * B(i).^(0:kmax);
    end
    CXX = sqrtM * (C - mu*mu')*sqrtM;
    SXX = diag(lambdas.^.5) * (C - mu*mu') * diag(lambdas.^.5);




    groups = cell(1,m);
    for i=1:m
        groups{i}= (i-1)*(kmax+1)+1:i*(kmax+1);
    end
    J = [];
    for i=indJ, J = [ J, groups{i} ]; end
    Jc = [];
    for i=setdiff(1:m,indJ), Jc = [ Jc, groups{i} ]; end

 


    % check condition
    invSigmaGJJ = inv( SXX(J,J) );
    W0norm = W0;
    for i=indJ,
        W0norm(groups{i}) = W0norm(groups{i}) / norm( W0norm(groups{i}) );
    end
    for i=setdiff(1:m,indJ)
        conditions(i) = norm( SXX(groups{i},J) * invSigmaGJJ * W0norm(J) );
    end
