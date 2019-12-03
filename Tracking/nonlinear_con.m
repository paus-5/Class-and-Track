function [c,ceq] = nonlinear_con(interactions,x)
[n, m] = size(x);
c = -reshape(1+interactions*x,n*m,1);     % Compute nonlinear inequalities at x.
ceq = 0;   % Compute nonlinear equalities at x.
end