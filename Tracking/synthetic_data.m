close all
clear
name_file = '250309_POC';
nA = 1;
nB = 1;
n = nA + nB;
a = ones(nA,1);
b = ones(nB,1);
class_AOB = [ones(nA,1); zeros(nB,1)];
class_NOB = [zeros(nA,1); ones(nB,1)];
time_steps_index = 1;
muARef = 0.77;
muBRef = 1.07;
kSARef = 0.7;
kSBRef = 1.3;
kIRef = 1.217;
x0 = [0.01*random('Unif',0,1,nA,1);0.01*random('Unif',0,1,nB,1); 0.1; 0.1; 0.1];
tF = 300;
t_pert = 150;
s_in_init = 1.25;
s_in_pert = 1.95;
d_interp = @(t) 0.22;
s_in_interp = @(t) s_in_init*(t<t_pert) + s_in_pert*(t >= t_pert);
kA = 1/0.251*ones(nA,1);
kB = 1/0.062*ones(nB,1);
% A = random('Normal',0,1,n,n);
A =[     1    0.1
    -1.5   1];
% A = zeros(n,n);
growth1 = @(X) (1+A(1:nA,:)*X(1:n))*growth_monod(X(n+1),muARef,kSARef);
growth2 = @(X) (1+A(nA+1:n,:)*X(1:n))*growth_monod(X(n+2),muBRef,kSBRef);
dynamic = @(t,X) [diag(growth1(X) - d_interp(t))*X(1:nA);
   diag(growth2(X) - d_interp(t))*X(nA+1:n);
   (s_in_interp(t) - X(n+1))*d_interp(t) - kA'*diag(growth1(X))*X(1:nA);
    -X(n+2)*d_interp(t) + kA'*diag(growth1(X))*X(1:nA)-kB'*diag(growth2(X))*X(nA+1:n);
    -X(n+3)*d_interp(t) + kB'*diag(growth2(X))*X(nA+1:n);];
[t_obs, Y] = ode15s(dynamic,[0 tF],x0);
S1 = Y(:,n+1);
S2 = Y(:,n+2);
S3 = Y(:,n+3);
OTU_interp = Y(:,1:n);
t_OTU = t_obs;
yields_AOB = [kA; zeros(nB,1)];
yields_NOB = [zeros(nA,1); kB];
save(sprintf('MAT_files\\synthetic_data_%s',name_file))