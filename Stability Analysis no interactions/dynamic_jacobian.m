function J = dynamic_jacobian(z,kA,kB,muA,muB,kSA,kSB,D)
J = zeros(5,5);
J(1,1) = muA*z(3)/(kSA+z(3))-D;
J(1,3) = muA*kSA*z(1)/(kSA+z(3))^2;
J(2,2) = muA*z(4)/(kSB+z(4))-D;
J(2,4) = muB*kSB*z(2)/(kSB+z(4))^2;
J(3,1) = -kA*muA*z(3)/(kSA+z(3));
J(3,3) = -D-kA*muA*kSA*z(1)/(kSA+z(3))^2;
J(4,1) = kA*muA*z(3)/(kSA + z(3));
J(4,2) = -kB*muB*z(4)/(kSB + z(4));
J(4,3) = kA*muA*kSA*z(1)/(kSA+z(3))^2;
J(4,4) = -D-kB*muB*kSB*z(2)/(kSB+z(4))^2;
J(5,5) = -D;
end