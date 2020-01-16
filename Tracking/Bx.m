function B = Bx(X,nA,nB,kA,kB,muA,muB,kSA,kSB,kI)
addpath('..//');
n = nA+ nB;
% growthRowVector = [growthMonod(X(n+1),muA,kSA)' growthHaldane(X(n+2),muB,kSB,kI)'];
% growth_row_vector = [growth_monod(X(n+1),muA,kSA)' growth_monod(X(n+2),muB,kSB)'];
growth_row_vector = [muA ; muB]';
B = [diag(growth_row_vector.*X(1:n)');
    [-kA; zeros(nB,1)]'.*growth_row_vector.*X(1:n)' ; 
    [kA; -kB]'.*growth_row_vector.*X(1:n)';
    [zeros(nA,1); kB]'.*growth_row_vector.*X(1:n)'];
end