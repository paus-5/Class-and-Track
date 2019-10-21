function B = Bx(X,nA,nB,kA,kB,muA,muB,kSA,kSB,kI)
addpath('..//');
n = nA+ nB;
% growthRowVector = [growthMonod(X(n+1),muA,kSA)' growthHaldane(X(n+2),muB,kSB,kI)'];
growth_row_vector = [growthMonod(X(n+1),muA,kSA)' growthMonod(X(n+2),muB,kSB)'];
B = [diag(growth_row_vector.*X(1:n)');
    [-kA; zeros(nB,1)]'.*growth_row_vector.*X(1:n)' ; 
    [kA; -kB]'.*growth_row_vector.*X(1:n)';
    [zeros(nA,1); kB]'.*growth_row_vector.*X(1:n)'];
end