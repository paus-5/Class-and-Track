function [a, b, yA, yB, errorAOB, errorNOB] = MIP_intlinprog(partialWeights, xG1, xG2,yA_ref,yB_ref,delta)
n = length(partialWeights(:,1));
T = length(partialWeights(1,:));
M1 = 1/(yA_ref*(1-delta)); %Reference Parameter Table 3 Weissmann et al 1994
M2 = 1/(yB_ref*(1-delta)); %big M constraints.
m1 = 1/(yA_ref*(1+delta)); 
m2 = 1/(yB_ref*(1+delta));
A = [-eye(T) zeros(T) partialWeights' zeros(T,3*n);
    -eye(T) zeros(T) -partialWeights' zeros(T,3*n);
    zeros(T) -eye(T) zeros(T,n) partialWeights' zeros(T,2*n);
    zeros(T) -eye(T) zeros(T,n) -partialWeights' zeros(T,2*n);
    zeros(n,2*T) eye(n) zeros(n) -M1*eye(n) zeros(n);
    zeros(n,2*T)  zeros(n) eye(n) zeros(n) -M2*eye(n) ;
    zeros(n,2*T) -eye(n) zeros(n) m1*eye(n) zeros(n);
    zeros(n,2*T)  zeros(n) -eye(n) zeros(n) m2*eye(n) ;
    zeros(n,2*T+2*n) eye(n) eye(n)];
f = [ones(2*T,1); zeros(4*n,1)];
b = [xG1; -xG1; xG2; -xG2; zeros(4*n,1); ones(n,1)];
lb = zeros(2*T+4*n,1);
ub = [inf*ones(2*T,1); M2*ones(2*n,1) ;ones(2*n,1)];
intcon = 2*T+ 2*n+1:2*T+4*n;
opts1=  optimoptions('intlinprog','Display','off');
[sol, ~] =  intlinprog(f,intcon,A,b,[],[],lb,ub,opts1);
errorAOB = sol(1:T);
errorNOB = sol(T+1:2*T);
yA = sol(2*T+1:2*T+n);
yB = sol(2*T+n+1:2*T+2*n);
a = sol(2*T+2*n+1:2*T+3*n);
b = sol(2*T+3*n+1:end);
end
