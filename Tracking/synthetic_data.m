close all
clear
name_file = 'no_noise';
nA = 10;
nB = 5;
n = nA + nB;
a = ones(nA,1);
b = ones(nB,1);
indexAOB = 1:nA;
indexNOB = nA+1:n;
time_steps_index = 1;
muARef = 0.77;
muBRef = 1.07;
kSARef = 0.7e-1;
kSBRef = 1.4982e-1;
kIRef = 1.217;
x0 = [0.01*ones(n,1); 0; 0; 0];
d_interp = @(t) 0.3;
s_in_interp = @(t) 0.5;
kA = 0.251*ones(nA,1);
kB = 0.062*ones(nB,1);
A = random('Normal',0,1,n,n);
growth1 = @(X) (1+A(1:nA,:)*X(1:n))*growth_monod(X(n+1),muARef,kSARef);
growth2 = @(X) (1+A(nA+1:n,:)*X(1:n))*growth_monod(X(n+2),muBRef,kSBRef);
dynamic = @(t,X) [diag(growth1(X) - d_interp(t))*X(1:nA);
   diag(growth2(X) - d_interp(t))*X(nA+1:n);
   (s_in_interp(t) - X(n+1))*d_interp(t) - kA'*diag(growth1(X))*X(1:nA);
    -X(n+2)*d_interp(t) + kA'*diag(growth1(X))*X(1:nA)-kB'*diag(growth2(X))*X(nA+1:n);
    -X(n+3)*d_interp(t) + kB'*diag(growth2(X))*X(nA+1:n);];
[t_obs, Y] = ode15s(dynamic,[0 500],x0);

subplot(3,1,1), plot(t_obs, Y(:,1:nA))
subplot(3,1,2), plot(t_obs, Y(:,nA+1:n))
subplot(3,1,3), plot(t_obs, Y(:,n+1:n+3))


save(sprintf('MAT_files\\synthetic_data_%s',name_file))