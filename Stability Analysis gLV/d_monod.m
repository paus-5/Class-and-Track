function y=d_monod(s,muMax,kS)
y =  muMax.*kS./(s+kS).^2;
end