close all
clear
file_name_in = 'parameters_modified_2';
load(sprintf('MAT_files\\%s',file_name_in));
file_name_out = 'parameters_modified_2';
number_of_points = 300;
d_initial = 0.7/6.5;
d_final =  2.2/6.5;
s_in_initial = 0.5;
s_in_final = 2;
s_in_vector = linspace(s_in_initial,s_in_final,number_of_points);
D_vector =  linspace(d_initial,d_final,number_of_points);
n1 = length(s_in_vector);
n2 = length(D_vector);
%Structure to save points.
map_equilibria = containers.Map;
zones = ones(n1,n2) ;
percentage = round(linspace(1,n1,11));
tic
for i = 1:n1
    for j = 1:n2
        [possible_x, possible_s1, possible_s2, possible_s3] = equilibria(A,...
            nA,nB,kA,kB,muA,muB,kSA,kSB,D_vector(j),s_in_vector(i));
        possible_equilibria = [possible_x; possible_s1; possible_s2; possible_s3];
        flag_PN = 0;
        flag_CN = 0;
        flag_washout = 0;
        for k=1:length(possible_equilibria(1,:))
            if all(possible_equilibria(:,k) >=0) && isreal(possible_equilibria(:,k))
                J = dynamic_jacobian(possible_equilibria(:,k),A,nA,nB,kA,kB,muA,muB,kSA,kSB,D_vector(j));
                if all(real(eig(J)) <0)
                    positive_equilibria = (possible_x(:,k)>0)';
                    index_aux = find(positive_equilibria);
                    if isempty(index_aux);
                        flag_washout = 1;
                    elseif max(index_aux) <= nA
                        flag_PN = 1;
                    else
                        flag_CN = 1;
                    end
                end
            end
        end
        if ~flag_PN && flag_CN
            zones(i,j) = 2;
        elseif flag_PN && flag_CN
            zones(i,j) = 3;   
        elseif flag_PN && ~flag_CN
            zones(i,j) = 4;
        elseif flag_washout && flag_PN
            zones(i,j) = 5;           
        elseif flag_washout && ~flag_PN
            zones(i,j) = 6;
        end
    end
    if any( i == percentage)
        fprintf('Progress: %.0f %% ',round(100*i/n1))
        toc
    end
end
save(sprintf('MAT_files\\Operating_Diagram_%s',file_name_out))
