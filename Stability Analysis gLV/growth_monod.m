function y = growth_monod(s,mu_max,kS)
y =  mu_max.*s./(s + kS);
end