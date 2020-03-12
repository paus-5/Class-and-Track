function [possible_x1, possible_x2, possible_s1, possible_s2, possible_s3] = equilibria(kA,kB,muA,muB,kSA,kSB,D,s_in)
possible_x1 = zeros(3,1);
possible_x2 = zeros(3,1);
possible_s1 = zeros(3,1);
possible_s2 = zeros(3,1);
possible_s3 = zeros(3,1);
%full coexistence
possible_s1(1) = kSA*D/(muA-D);
possible_x1(1) = (s_in-possible_s1(1))/kA;
possible_s2(1) = kSB*D/(muB-D);
possible_x2(1) = (s_in-possible_s1(1)-possible_s2(1))/kB;
possible_s3(1) = s_in - possible_s1(1) - possible_s2(1);
% washout of NOB
possible_s1(2) = kSA*D/(muA-D);
possible_x1(2) = (s_in-possible_s1(2))/kA;
possible_s2(2) = s_in - possible_s1(2);
%washout
possible_s1(3) = s_in;
end