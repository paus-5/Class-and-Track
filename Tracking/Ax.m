function A = Ax(X,nA,nB,kA,kB,muA,muB,kSA,kSB,kI,sIn,dilutionRate)
addpath('..//');
n = nA+ nB;
% growthRowVector = [growthMonod(X(n+1),muA,kSA)' growthHaldane(X(n+2),muB,kSB,kI)'];
growth_row_vector = [growthMonod(X(n+1),muA,kSA)' growthMonod(X(n+2),muB,kSB)'];
A = [diag(growth_row_vector-dilutionRate) zeros(n,3);
    [-kA; zeros(nB,1)]'.*growth_row_vector (sIn/X(n+1) -1)*dilutionRate 0 0; 
    [kA; -kB]'.*growth_row_vector 0 -dilutionRate 0;
    [zeros(nA,1); kB]'.*growth_row_vector zeros(1,2) -dilutionRate] ;
end