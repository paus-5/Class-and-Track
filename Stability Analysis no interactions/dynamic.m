function y = dynamic(kA,kB,muA,muB,kSA,kSB,D,S_in)
dynamicX =  @(t,z) diag(z(1:2))*([growth_monod(z(3),muA,kSA); growth_monod(z(4),muB,kSB)]- ones(2,1)*D(t));
dynamicS1 =  @(t,z) (S_in(t) - z(3))*D(t) -kA*growth_monod(z(3),muA,kSA)*z(1);
dynamicS2 =  @(t,z) -z(4)*D(t) + kA*growth_monod(z(3),muA,kSA)*z(1) - kB*growth_monod(z(4),muB,kSB)*z(2);
dynamicS3 =  @(t,z) -z(5)*D(t) + kB*growth_monod(z(4),muB,kSB)*z(2);
y  = @(t,z)  [dynamicX(t,z) ; dynamicS1(t,z); dynamicS2(t,z); dynamicS3(t,z)];
end