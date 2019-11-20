function y = dynamic(A,nA,nB,kA,kB,muA,muB,kSA,kSB,D,S_in)
n = nA+nB;
dynamicX =  @(t,z) diag(z(1:n))*(diag([growth_monod(z(n+1),muA,kSA); growth_monod(z(n+2),muB,kSB)])*(ones(n,1)+A*z(1:n)) - ones(n,1)*D(t));
dynamicS1 =  @(t,z) (S_in(t) - z(n+1))*D(t) +(diag(-kA)*(ones(nA,1) + A(1:nA,:)*z(1:n)))'*diag(growth_monod(z(n+1),muA,kSA))*z(1:nA);
dynamicS2 =  @(t,z) -z(n+2)*D(t) + (diag([kA; -kB])*(ones(n,1)+A*z(1:n)))'*diag([growth_monod(z(n+1),muA,kSA); growth_monod(z(n+2),muB,kSB)])*z(1:n);
dynamicS3 =  @(t,z) -z(n+3)*D(t) +(diag(kB)*(ones(nB,1) + A((nA+1):n,:)*z(1:n)))'*diag(growth_monod(z(n+2),muB,kSB))*z((nA+1):n);
y  = @(t,z)  [dynamicX(t,z) ; dynamicS1(t,z); dynamicS2(t,z); dynamicS3(t,z)];
end