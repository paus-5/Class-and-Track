function [a, b, y_AOB, y_NOB, error_AOB, error_NOB] = MIP_gurobi(partialWeights, xG1, xG2,yA_ref,yB_ref,delta)
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
b = [xG1; -xG1; xG2; -xG2; zeros(4*n,1); ones(n,1)];
lb = zeros(2*T+4*n,1);
ub = [inf*ones(2*T,1); M2*ones(2*n,1) ;ones(2*n,1)];
intcon = 2*T+ 2*n+1:2*T+4*n;
variables_type = '';
for i=1:2*T+ 2*n
    variables_type = strcat(variables_type,'C');
end
for i=intcon
    variables_type = strcat(variables_type,'B');
end
model.A = sparse(A);
model.Q = sparse(blkdiag(eye((2*T)),zeros(4*n,4*n)));
% model.obj = [ones(2*T,1); zeros(4*n,1)];
model.rhs = b;
model.sense = '<';
model.vtype = variables_type;
model.modelsense = 'min';
model.lb = lb;
model.ub = ub;
gurobi_write(model, 'mip1.lp');
params.outputflag = 0;
result = gurobi(model, params);
sol = result.x;
error_AOB = sol(1:T);
error_NOB = sol(T+1:2*T);
y_AOB = sol(2*T+1:2*T+n);
y_NOB = sol(2*T+n+1:2*T+2*n);
a = sol(2*T+2*n+1:2*T+3*n);
b = sol(2*T+3*n+1:end);
end
