%Reference: Cimen et al 2002 Nonlinear optimal tracking control
%main: process and filter data to model.
close all
clear
load('MAT_files\Classification_day_183_315_gurobi')
name_file = '191218_Reactor_A';
nA = round(sum(class_AOB));
nB = round(sum(class_NOB)); 
n = nA+nB;
index_t_obs = find(t_obs <315);
index_span = time_steps_index(1):index_t_obs(end);
tS = t_obs(index_span); 
s1_interp = @(t) interp1(tS,S1(index_span),t);
s2_interp = @(t) interp1(tS,S2(index_span),t);
s3_interp = @(t) interp1(tS,S3(index_span),t);
index_AOB = find(class_AOB);
index_NOB = find(class_NOB);
biomass_filtered = OTU_interp(:,[index_AOB; index_NOB]);
min_biomass = min(biomass_filtered(biomass_filtered>0));
biomass_filtered(biomass_filtered == 0) = min_biomass;
z = @(t) interp1(t_OTU,biomass_filtered,t)';
kA = yields_AOB(index_AOB);
kB = yields_NOB(index_NOB);
muARef = 0.77;
muBRef = 1.07;
kSARef = 0.7e-1;
kSBRef = 1.4982e-1;
kIRef = 1.217;
muA = muARef*ones(nA,1);
muB = muBRef*ones(nB,1);
kSA = kSARef*ones(nA,1);
kSB = kSBRef*ones(nB,1);
kI = kIRef*ones(nB,1);
% %%Parameters
% lambda = 0.001;
% lambda = [0.001*ones(nA,1); 0.001*ones(nB,1)];
lambda = ones(n,1)*10^-4;
Q2= [kA; kB];
% Q = Q2/norm(Q2,'inf');
Q = [ones(nA,1); ones(nB,1)];
% vectorWeight = 1./max(biomass);
% vectorWeight = [kA; kB].^2;
% Q = diag(vectorWeight/max(vectorWeight));
% Q = eye(n);
F = zeros(n,1);
colors = mat2cell([autumn(nA); winter(nB)], ones(1,n), 3);
colors2 = mat2cell(hsv(3), ones(1,3),3);
legend_tags = strsplit(strjoin({num2str(1:n),' s_1 s_2 s_3'}));
legend_tags_OTU = strsplit(strjoin({num2str(1:n)}));
tF = tS(end);
t0 = tS(1);
x0 = [biomass_filtered(1,:) s1_interp(t0) s2_interp(t0) s3_interp(t0)]';
dynamic = @(t,X) Ax(X,nA,nB,kA,kB,muA,muB,kSA,kSB,kI,s_in_interp(t),...
    d_interp(t))*X ;
options = odeset('NonNegative',(1:n+3)');
[t_new,x_new] = ode15s(dynamic,[t0 tF], x0,options);
iter = 0;
save(sprintf('MAT_files\\%s_iter_%.0f',name_file,iter))
