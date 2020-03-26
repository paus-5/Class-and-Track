close all
clear
file_name_in = 'parameters_modified_2';
load(sprintf('MAT_files\\%s',file_name_in));
file_name_out = 'parameters_modified_2';
number_of_points = 400;
d_initial = 0.7/6.5;
d_final =  2.2/6.5;
s_in_initial = 0.5;
s_in_final = 2;
s_in_vector = linspace(s_in_initial,s_in_final,number_of_points);
D_vector =  linspace(d_initial,d_final,number_of_points);
n1 = length(s_in_vector);
n2 = length(D_vector);
%Structure to save points.
map_zones_ED = containers.Map;
zones_OD = ones(n1,n2);
zones_ED = zeros(n1,n2);
percentage = round(linspace(1,n1,11));
tolerance = 1e-9;
tic
for i = 1:n1
    for j = 1:n2
        s_in = s_in_vector(i);
        D = D_vector(j);
%         fun =  dynamic(A,nA,nB,kA,kB,muA,muB,kSA,kSB,@(t) D,@(t) s_in);
        [possible_x, possible_s1, possible_s2, possible_s3] = equilibria(A,...
            nA,nB,kA,kB,muA,muB,kSA,kSB,D,s_in);
        possible_equilibria = [possible_x; possible_s1; possible_s2; possible_s3];
        flag_PN = 0;
        flag_CN = 0;
        flag_washout = 0;
        key_string = '';
        for k=1:length(possible_equilibria(1,:))
            c_positivity = all(possible_equilibria(:,k) >=0);
            c_real_number = isreal(possible_equilibria(:,k));
%             c_evaluation_in_dynamic =all(fun(1,possible_equilibria(:,k))< tolerance);
            if c_positivity && c_real_number % && c_evaluation_in_dynamic
                J = dynamic_jacobian(possible_equilibria(:,k),A,nA,nB,kA,kB,muA,muB,kSA,kSB,D);
                eigenvalues = eig(J);
                if all(real(eigenvalues) <0)
                    positive_equilibria = (possible_x(:,k)>0)';
                    index_aux = find(positive_equilibria);
                    key_string = strcat(key_string,num2str(positive_equilibria),',');
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
            zones_OD(i,j) = 2;
        elseif flag_PN && flag_CN
            zones_OD(i,j) = 3;   
        elseif flag_PN && ~flag_CN
            zones_OD(i,j) = 4;
        elseif flag_washout && flag_PN
            zones_OD(i,j) = 5;           
        elseif flag_washout && ~flag_PN
            zones_OD(i,j) = 6;
        end
        if  isKey(map_zones_ED,key_string)
            zones_ED(i,j) = map_zones_ED(key_string);
        else
            count = length(map_zones_ED);
            map_zones_ED(key_string) = count + 1;
            zones_ED(i,j) = map_zones_ED(key_string);
        end
    end
    if any( i == percentage)
        fprintf('Progress: %.0f %% ',round(100*i/n1))
        toc
    end
end
save(sprintf('MAT_files\\Operating_Diagram_%s',file_name_out))
