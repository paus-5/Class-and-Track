function J = dynamic_jacobian(z,A,nA,nB,kA,kB,muA,muB,kSA,kSB,D)
n = nA+nB;
J = zeros(n+3,n+3);
A1 = A(1:nA,:);
A2 = A(nA+1:n,:);
growthA = growth_monod(z(n+1),muA,kSA);
growthB = growth_monod(z(n+2),muB,kSB);
d_growthA = d_monod(z(n+1),muA,kSA);
d_growthB = d_monod(z(n+2),muB,kSB);
%growth(s,muMax,kS)
%dGrowth(s,muMax,kS)
J(1:n, 1:n) = diag([growthA; growthB].*(ones(n,1) +A*z(1:n))-D*ones(n,1))+diag([growthA; growthB].*z(1:n))*A;
J(1:n,n+1) = z(1:n).*[d_growthA ; zeros(nB,1)].*(ones(n,1) + A*z(1:n));
J(1:n,n+2) = z(1:n).*[zeros(nA,1) ; d_growthB].*(ones(n,1) + A*z(1:n));
J(1:n,n+3) = zeros(n,1);
J(n+1,1:n) = -([ones(nA,1) + A1*z(1:n) ; zeros(nB,1)]'*diag([kA; kB].*[growthA; growthB]))...
    + (diag([kA; kB].*[growthA; growthB])*z(1:n))'*[A1; zeros(nB,n)]  ;
J(n+1,n+1) = -D - [ones(nA,1) + A1*z(1:n) ; zeros(nB,1)]'*([kA; kB].*[d_growthA; zeros(nB,1)].*z(1:n));
J(n+1,n+2:n+3) = zeros(1,2);
J(n+2,1:n) = [(ones(nA,1) + A1*z(1:n)) ;-(ones(nB,1) +A2*z(1:n))]'*diag([kA;kB].*[growthA; growthB])...
    + (diag([kA;kB].*[growthA; growthB])*z(1:n))'*[A1; -A2];
J(n+2,n+1) = [(ones(nA,1) + A1*z(1:n)) ;-(ones(nB,1) +A2*z(1:n))]'*([kA; kB].*[d_growthA; zeros(nB,1)].*z(1:n));
J(n+2,n+2) = -D+ [(ones(nA,1) + A1*z(1:n)) ;-(ones(nB,1) +A2*z(1:n))]'*([kA; kB].*[zeros(nA,1) ; d_growthB].*z(1:n));
J(n+2,n+3) = 0;
J(n+3,1:n) = [zeros(nA,1); ones(nB,1) + A2*z(1:n)]'*diag([kA; kB].*[growthA; growthB])...
    +  (diag([kA;kB].*[growthA; growthB])*z(1:n))'*[zeros(nA,n); A2] ;
J(n+3,n+1) = 0;
J(n+3,n+2) = [zeros(nA,1) ;((ones(nB,1) + A2*z(1:n)))]'*([kA; kB].*[zeros(nA,1) ; d_growthB].*z(1:n));
J(n+3,n+3) = -D;
end