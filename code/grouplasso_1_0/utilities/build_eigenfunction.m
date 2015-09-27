function f = build_eigenfunction(a,b,k,x)
c = ( a.^2 + 2 * a .* b).^.5;
A = a + b + c;
B = b./A;
lambda = (2*a/A)^.5 * B^k;
f = lambda^.5 * ( c^.5 / a^.5 / 2^k / factorial(k) )^.5 * exp( -(c-a) .* x.*x ) ...
    .* hermite_polynomials(k+1,x * sqrt(2*c) );
