close all
clear
load('MAT_files\Classification_day_183_315_gurobi')
addpath('..\\General Functions')
reactor = 'B';
[t_obs,biomass,s_in,S1,S2,S3,dilution_rate] = load_biomass(reactor);
[t_OTU,OTU_rel] = load_relative_abundance(reactor);
OTU_rel = diag(1./sum(OTU_rel,2))*OTU_rel;
number_OTU = length(OTU_rel(1,:));
minimal_abundance = 0;
t_change = 183; % From dumont et al
delta_t = 132;
time_steps_index = find(t_obs > t_change & t_obs < t_change + delta_t);
% time_steps_index = 1:length(t_obs);
time_steps_index_OTU = find(t_OTU > t_change & t_OTU < t_change + delta_t);
OTU_series = timeseries(OTU_rel(time_steps_index_OTU,:),...
    t_OTU(time_steps_index_OTU));
mean_abundance = mean(OTU_series,'Weighting','time')';
max_abundance =  max(OTU_series)';
considered_OTU = find(max_abundance > minimal_abundance);
OTU_interp = interpolate_biomass(t_obs,biomass,t_OTU,OTU_rel);
OTU_interp_fun = @(t) interp1(t_OTU,OTU_interp(:,considered_OTU),t);
d_interp = @(t) interp1(t_obs,dilution_rate,t);
s_in_interp = @(t) interp1(t_obs,s_in,t,'previous');
s1_interp = @(t) interp1(t_obs,S1,t);
s2_interp = @(t) interp1(t_obs,S2,t);
save('MAT_files\Classification_day183_315_A_for_B')
