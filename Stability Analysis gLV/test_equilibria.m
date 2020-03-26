clear
close all
% generateParameters;
file_name_in = 'parameters_modified_2';
load(sprintf('MAT_files\\%s',file_name_in));
z0 = random('Uniform',0,1,n+3,1);
OTU = 1:n;
options = odeset('NonNegative',1:(nA+nB+3));
s_in = 1.602;
sIn = @(t) s_in;
dilution = 0.1925 ;
D = @(t) dilution;
tF = 20000;
fun =  dynamic(A,nA,nB,kA,kB,muA,muB,kSA,kSB,D,sIn);
[T,Y] = ode15s(fun,[0 tF], z0,options);
end_point_string = sprintf('%1.4f ', Y(end,:)');
eval_end_point =  fun(tF,Y(end,:)');
eval_end_point_string = sprintf('%1.4f ', eval_end_point);
fprintf('Trajectory at T= %d \n Point =%s \n Evaluated in dynamic = %s\n',...
    tF, end_point_string,eval_end_point_string);
J = dynamic_jacobian(Y(end,:)',A,nA,nB,kA,kB,muA,muB,kSA,kSB,dilution);
eigenvalues = eig(J);
if all(real(eigenvalues) <0)
    fprintf('Stable Equilibrium: %s\n',sprintf( '%1.4f', Y(end,:)'));
    disp(eigenvalues)
end
[possible_x, possible_s1, possible_s2, possible_s3]  = equilibria(A,...
    nA,nB,kA,kB,muA,muB,kSA,kSB,dilution,s_in);
possible_equilibria = [possible_x; possible_s1; possible_s2; possible_s3];
tolerance = 1e-6;
for i=1:length(possible_equilibria(1,:))
    c_positivity = all(possible_equilibria(:,i) >=0);
    c_real_number = isreal(possible_equilibria(:,i));
    c_evaluation_in_dynamic =all(fun(1,possible_equilibria(:,i))< tolerance);
    if c_positivity && c_real_number && c_evaluation_in_dynamic
        J = dynamic_jacobian(possible_equilibria(:,i),A,nA,nB,kA,kB,muA,muB,kSA,kSB,dilution);
        eigenvalues = eig(J);
        if all(real(eigenvalues) <0)
            disp(possible_equilibria(:,i)')
            fprintf('Equilibrium %d: Stable Equilibrium! \n Eigenvalues of the Jacobian in the Equilibrium : \n', i);
            disp(eigenvalues)
        end
    end
end