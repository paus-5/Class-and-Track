function f = sF_diff(sF,X,z,P,nA,nB,yA,yB,muA,muB,kSA,kSB,kI,sIn,dilutionRate,...
    lambda)
n = nA + nB;
growth_vector = [growth_monod(X(n+1),muA,kSA); growth_monod(X(n+2),muB,kSB)];
x_trunc = X(1:n);
f = -z -((growth_vector- dilutionRate) - 1./lambda.*(x_trunc.^2).*growth_vector.^2.*P).*sF;
end