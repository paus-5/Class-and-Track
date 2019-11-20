function y = inv_growth_monod(s,mu_max,kS)
y = 1./mu_max.*(1+kS./s);
end