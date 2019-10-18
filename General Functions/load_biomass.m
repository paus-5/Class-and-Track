function [t_obs,X,s_in,S1,S2,S3,dilution_rate] = load_biomass(Reactor)
data = xlsread(sprintf('..\\Data\\%s',Reactor));
t_obs = data(:,1)-data(1,1);
X = data(:,2);
s_in = data(:,3);
S1 = data(:,4);
S2 = data(:,5);
S3 = data(:,6);
flow = data(:,7);
dilution_rate = flow/6.5;
end