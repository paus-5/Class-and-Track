function [possible_x, possible_s1, possible_s2, possible_s3] = equilibria(A,...
    nA,nB,kA,kB,muA,muB,kSA,kSB,D,s_in)
n = nA+nB;
possible_x = zeros(n, (2^nA-1)*(4*(2^nB-1)+2)+1);
possible_s1 = zeros(1,(2^nA-1)*(4*(2^nB-1)+2)+1);
possible_s2 = zeros(1,(2^nA-1)*(4*(2^nB-1)+2)+1);
possible_s3 = zeros(1,(2^nA-1)*(4*(2^nB-1)+2)+1);
counter = 1;
k = [kA;kB];
vec_aux1 = (D-[muA; muB])./[muA; muB];
vec_aux2 = D*[kSA ; kSB]./[muA; muB];
for i = 1:2^nA-1
    %CN
    for j = 1:2^nB-1
        considered_species = find([de2bi(i,nA) de2bi(j,nB)]);
        A_act = A(considered_species, considered_species);
        B = diag(k(considered_species))/A_act;
        AOB_species = find(considered_species <= nA);
        NOB_species = find(considered_species > nA);
        beta1 = -sum(B(AOB_species,:)*vec_aux1(considered_species));
        beta2 = -sum(B(AOB_species,AOB_species)*vec_aux2(considered_species(AOB_species)));
        beta3 = sum(B(AOB_species,NOB_species)*vec_aux2(considered_species(NOB_species)));
        beta4 = -sum(B(NOB_species,AOB_species)*vec_aux2(considered_species(AOB_species)));
        beta5 = -sum(B(NOB_species,:)*vec_aux1(considered_species));
        beta6 = -sum(B(NOB_species,NOB_species)*vec_aux2(considered_species(NOB_species)));
        c1 = -1/beta3;
        c2 = (s_in+beta1)/beta3;
        c3 = beta2/beta3;
        a4 = (beta3+beta6)*c1^2;
        a3 = (beta3+beta6)*2*c1*c2+(beta5-beta1)*c1;
        a2 = (beta4-beta2)*c1+(beta3+beta6)*(c2^2+2*c1*c3)+(beta5-beta1)*c2-1;
        a1 = (beta4-beta2)*c2+(beta3+beta6)*2*c2*c3+(beta5-beta1)*c3;
        a0 = (beta4-beta2)*c3+(beta3+beta6)*c3^2;
        polySol = roots ([a4 a3 a2 a1 a0])';
        possible_s1(counter:counter+3) = polySol;
        possible_s2(counter:counter+3) = polySol./(c1.*polySol.^2+c2.*polySol+c3);
        possible_s3(counter:counter+3) = s_in - possible_s1(counter:counter+3) - possible_s2(counter:counter+3);
        for count=counter:counter+3
            s1 = possible_s1(count);
            s2 = possible_s2(count);
            mu_inv_aux = [inv_growth_monod(s1,muA,kSA); inv_growth_monod(s2,muB,kSB);];
            possible_x(considered_species,count) = A_act\(mu_inv_aux(considered_species)*D-ones(length(considered_species),1));
        end
        counter = counter +4;
    end
    %PN
    considered_species = find([de2bi(i,nA) de2bi(0,nB)]);
    A_act = A(considered_species, considered_species);
    B = diag(k(considered_species))/A_act;
    AOB_species = find(considered_species <= nA);
    a2 = -1;
    a1 = s_in - sum(B(AOB_species,AOB_species)*vec_aux1(considered_species(AOB_species)));
    a0 = -sum(B(AOB_species,AOB_species)*vec_aux2(considered_species(AOB_species)));
    polySol = roots ([a2 a1 a0])';
    possible_s1(counter:counter+1) = polySol;
    possible_s2(counter:counter+1) =  s_in - possible_s1(counter:counter+1);
    for count=counter:counter+1
        s1 = possible_s1(count);
        s2 = possible_s2(count);
        mu_inv_aux = [inv_growth_monod(s1,muA,kSA); inv_growth_monod(s2,muB,kSB);];
        possible_x(considered_species,count) = A_act\(mu_inv_aux(considered_species)*D-ones(length(considered_species),1));
    end
    counter = counter + 2;
end
possible_s1(end) = s_in;
end