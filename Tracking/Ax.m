function A = Ax(X,nA,nB,kA,kB,muA,muB,kSA,kSB,kI,sIn,dilution_rate)
addpath('..//');
n = nA+ nB;
% growthRowVector = [growthMonod(X(n+1),muA,kSA)' growthHaldane(X(n+2),muB,kSB,kI)'];
% growth_row_vector = [growth_monod(X(n+1),muA,kSA)' growth_monod(X(n+2),muB,kSB)'];
growth_row_vector = [muA ; muB]';
A = [diag(growth_row_vector-dilution_rate) zeros(n,3);
    [-kA; zeros(nB,1)]'.*growth_row_vector (sIn/X(n+1) -1)*dilution_rate 0 0; 
    [kA; -kB]'.*growth_row_vector 0 -dilution_rate 0;
    [zeros(nA,1); kB]'.*growth_row_vector zeros(1,2) -dilution_rate] ;
end