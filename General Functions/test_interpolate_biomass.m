clear
close all
addpath('..\\General Functions')
Reactor = 'B';
[t_obs,biomass,s_in,s1,s2,s3,dilutionRate] = load_biomass(Reactor);
[t_OTU,OTU_rel] = load_relative_abundance(Reactor);
OTU_interp = interpolate_biomass(t_obs,biomass,t_OTU,OTU_rel);
plot(t_OTU,OTU_interp,'LineWidth',3)
    xlabel('\fontsize{12}Time [days]'),...
    ylabel('\fontsize{12}Concentration [g/l]'),...
    title(sprintf('OTUs Biomass in Time Reactor %s',Reactor))
    set(gca,'fontsize',15);
