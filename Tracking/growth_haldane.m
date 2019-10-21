function y = growth_haldane(s,mu_max,kS,kI)
y =  mu_max.*s./(s + kS + s^2./kI);
end