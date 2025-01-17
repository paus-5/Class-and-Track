clear
addpath('..\\General Functions')
reactor = 'A';
delta = 0.01; %Between 0 and 1
id = 'day_183_315_gurobi';
file_name = sprintf('Classification_%s',id);
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
dynamic = @(t,z) -1*d_interp(t)*(z - s_in_interp(t));
t0 = t_obs(time_steps_index(1));
tF = t_obs(time_steps_index(end));
[T, Z] = ode15s(dynamic,[t0 tF],s_in(time_steps_index(1)));
invariant = @(t) interp1(T,Z,t);
partial_weights = OTU_interp_fun(t_obs(time_steps_index))';
Z1 = invariant(t_obs(time_steps_index));
% xG1 = max(0,Z1 - S1);
% xG2 = max(0,Z1 - S1 - S2);
xG1 = Z1 - S1(time_steps_index);
xG2 = Z1 - S1(time_steps_index) - S2(time_steps_index);
yA_ref = 0.2229; %Values obtained from y_opt observer_study
yB_ref = 0.0215;
% yA_ref = 0.2039;
% yB_ref = 0.0243;
[a, b, y_AOB, y_NOB, error_AOB, error_NOB] = MIP_gurobi(partial_weights, xG1, xG2,yA_ref,yB_ref,delta);
% [a, b, y_AOB, y_NOB, error_AOB, error_NOB] = MIP_intlinprog(partial_weights, xG1, xG2,yA_ref,yB_ref,delta);
class_AOB = zeros(number_OTU,1);
class_NOB = zeros(number_OTU,1);
yields_AOB = zeros(number_OTU,1);
yields_NOB = zeros(number_OTU,1);
class_AOB(considered_OTU) = round(a);
class_NOB(considered_OTU) = round(b);
yields_AOB(considered_OTU) = y_AOB;
yields_NOB(considered_OTU) = y_NOB;
save(sprintf('MAT_files\\%s',file_name));