function f = ricatti_diff(Plin,X,Q,nA,nB,kA,kB,muA,muB,kSA,kSB,kI,sIn,dilution_rate,...
    lambda)
n = nA + nB;
growth_vector = [growth_monod(X(n+1),muA,kSA); growth_monod(X(n+2),muB,kSB)];
x_trunc = X(1:n);
f = -Q -2.*Plin.*(growth_vector- dilution_rate)...
    + 1./lambda.*Plin.^2.*growth_vector.^2.*x_trunc.^2;
end