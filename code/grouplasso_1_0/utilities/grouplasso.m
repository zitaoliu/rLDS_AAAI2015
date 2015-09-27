function [f,grad,hess,lambda2, descent] = grouplasso(W,Sigma,Q,groups,lambda,alpha)
dm = size(Q,1);
SigmaW = Sigma * W;
f = .5 * W' * SigmaW - sum( W .* Q );
s = W2norms(W,groups);

f = f + lambda * sum( s + 2 * alpha * log( 1 + exp( -s / alpha ) ) );

if nargout>1
    % compute gradient
    grad = SigmaW - Q;
    for i=1:length(groups)
        grad(groups{i}) = grad(groups{i}) + lambda * W(groups{i}) ...
             / alpha * phiexp(s(i)/alpha);
             
    end
    
end
if nargout > 2
    % compute hessian
    hess = Sigma;


    for i=1:length(groups)

        hess(groups{i},groups{i}) = hess(groups{i},groups{i}) + lambda * ( ...
        eye(length(groups{i})) / alpha * phiexp(s(i)/alpha) ...
        + W(groups{i})  * W(groups{i})' * ( ...
        1/alpha^3 * phiexp2(s(i)/alpha) ...
    )   );
           
    end
    
    hess = (hess + hess')./2;
    
%     error_flag = 0;
%     eval('[eigVecs, eigVals] = eig(hess);', 'error_flag = 1;');
%     if (error_flag)
%         [eigVecs, eigVals] = eig(hess + 0.001 * eye(size(hess)));
%     else
%        [eigVecs, eigVals] = eig(hess);
%     end
    
    [eigVecs, eigVals] = eig(hess);
    
    hessInv = psdInverse( eigVecs, eigVals );
    
    descent = - hessInv * grad;
    lambda2 = - grad' * descent;
    %if isnan(lambda2), keyboard; end
end

function phi = phiexp(t);

if t>1e-8
    phi = (1-exp(-t))/t/(1+exp(-t));
else
   phi = ( 1+ t/2) / ( 1 + exp(t) ); 
end



function phi = phiexp2(t);

if t>1e-8
    phi = exp(-2*t) / ( 1 + exp(-t) )^2 * (  ( 2 * exp(-t) - (  1 - exp(-2*t)  ) / t ) / t / t );
else
   phi = ( -1/3 - t/3 ) / ( 1 + exp(-t) )^2; 
end