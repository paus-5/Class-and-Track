clear
close all
% generateParameters;
file_name_in = 'parameters_modified';
load(sprintf('MAT_files\\%s',file_name_in));
z0 = random('Uniform',0,1,n+3,1);
OTU = 1:n;
options = odeset('NonNegative',1:(nA+nB+3));
dilution = 0.2906 ;
D = @(t) dilution;
s_in = 1.829;
sIn = @(t) s_in;
tF = 20000;
fun =  dynamic(A,nA,nB,kA,kB,muA,muB,kSA,kSB,D,sIn);
[T,Y] = ode15s(fun,[0 tF], z0,options);
end_point_string = sprintf('%1.4f ', Y(end,:)');
eval_end_point =  fun(tF,Y(end,:)');
eval_end_point_string = sprintf('%1.4f ', eval_end_point);
fprintf('Trajectory at T= %d \n Point =%s \n Evaluated in dynamic = %s\n',...
    tF, end_point_string,eval_end_point_string);
[possible_x, possible_s1, possible_s2, possible_s3]  = equilibria(A,...
    nA,nB,kA,kB,muA,muB,kSA,kSB,dilution,s_in);
possible_equilibria = [possible_x; possible_s1; possible_s2; possible_s3];
num_equilibria = 0;

for i=1:length(possible_equilibria(1,:))
    if all(possible_equilibria(:,i) >=0) && isreal(possible_equilibria(:,i))
        num_equilibria = num_equilibria +1;
%         fprintf('\n Positive Equilibrium %d =', num_equilibria);
%         disp(possible_equilibria(:,i)')
        J = dynamic_jacobian(possible_equilibria(:,i),A,nA,nB,kA,kB,muA,muB,kSA,kSB,dilution);
        eigenvalues = eig(J);
        if all(real(eigenvalues) <0)
            disp(possible_equilibria(:,i)')
            fprintf('Equilibrium %d: Stable Equilibrium! \n Eigenvalues of the Jacobian in the Equilibrium : \n', num_equilibria);
            disp(eigenvalues)
        end
    end
end