close all
clear
file_name_in = 'Parameters_modified';
load(sprintf('MAT_files\\%s',file_name_in));
file_name_out = 'Parameters_modified';
number_of_points = 200;
d_initial = 0.7/6.5;
d_final = 2.2/6.5;
s_in_initial = 0.5;
s_in_final = 2;
s_in_vector = linspace(s_in_initial,s_in_final,number_of_points);
D_vector =  linspace(d_initial,d_final,number_of_points);
n1 = length(s_in_vector);
n2 = length(D_vector);
%parametersLVarticle;
%Structure to save points.
map_zones = containers.Map;
zones = zeros(n1,n2);
tic
for i = 1:n1
    for j = 1:n2
        [possible_x, possible_s1, possible_s2, possible_s3] = equilibria(A,...
            nA,nB,kA,kB,muA,muB,kSA,kSB,D_vector(j),s_in_vector(i));
        possible_equilibria = [possible_x; possible_s1; possible_s2; possible_s3];
        for k=1:length(possible_equilibria(1,:))
            if all(possible_equilibria(:,k) >=0) && isreal(possible_equilibria(:,k))
                J = dynamic_jacobian(possible_equilibria(:,k),A,nA,nB,kA,kB,muA,muB,kSA,kSB,D_vector(j));
                if all(real(eig(J)) <0)
                    positive_equilibria = (possible_x(:,k)>0)';
                    keyString = num2str(positive_equilibria);
                    if  isKey(map_zones,keyString)
                        zones(i,j) = map_zones(keyString);
                    else
                        count = length(map_zones);
                        map_zones(keyString) = count + 1;
                        zones(i,j) = map_zones(keyString);
                    end
                end
            end
        end
    end
end
toc
save(sprintf('MAT_files\\Operating_Diagram_%s',file_name_out))
