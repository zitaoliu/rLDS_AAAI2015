function D = compute_DklQ(k,l,Q)
global factorials;
k = k+1;
l = l+1;
D = 0;
% check same parity
if mod(k,2)~=mod(l,2), return; end
Qinv = inv(Q);
if mod(k,2)==1
    % both even (k-1 not k!)
    k = (k-1) / 2;
    l = (l-1) / 2;
    temp =  (-1)^(k+l) * pi / sqrt( det(Q) );
    for i=0:min(k,l)
        D =  D + 4^i * ( 1-Qinv(1,1) )^(k-i) * ( 1-Qinv(2,2) )^(l-i) * Qinv(1,2)^(2*i) * ...
            factorials(2*l+1) * factorials(2*k+1) / factorials(2*i+1) / factorials(k-i+1) / factorials(l-i+1);
    end
    D = D * temp;
  
else
    % both odd  (k-1 not k!)
    k = (k-2) / 2;
    l = (l-2) / 2;
    temp =  (-1)^(k+l) * pi / sqrt( det(Q) );
    for i=0:min(k,l)
        D =  D + 2^(2*i+1) * ( 1-Qinv(1,1) )^(k-i) * ( 1-Qinv(2,2) )^(l-i) * Qinv(1,2)^(2*i+1) * ...
            factorials(2*l+1+1) * factorials(2*k+1+1) / factorials(2*i+1+1) / factorials(k-i+1) / factorials(l-i+1);
    end
    D = D * temp;
end
