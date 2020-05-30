clear
close all
addpath('..\\General Functions')
reactor = 'B';
[t_obs,biomass,s_in,S1,S2,S3,dilution_rate] = load_biomass(reactor);
[t_OTU_rel,OTU_rel] = load_relative_abundance(reactor);
d_inter = @(t) interp1(t_obs,dilution_rate,t,'previous');
s_in_inter = @(t) interp1(t_obs,s_in,t,'previous');
s1_inter = @(t) interp1(t_obs,S1,t);
s2_inter = @(t) interp1(t_obs,S2,t);
s3_inter = @(t) interp1(t_obs,S3,t);
biomas_inter = @(t) interp1(t_obs,biomass,t);
t_initial = 183; %incrase of temperature day = 183
delta_t = 132;
time_steps_index = find(t_obs > t_initial & t_obs < t_initial + delta_t);
yA_ref = 0.147;
yB_ref = 0.042;
Y = [-yA_ref 0; yA_ref -yB_ref;0 yB_ref];
t_used = t_obs(time_steps_index);
z0 = [s1_inter(t_used(1)) s2_inter(t_used(1)) s3_inter(t_used(1))]' ;
dynamic = @(t,z) -1*d_inter(t)*(z - [s_in_inter(t) 0 0]');
% z0 = random('Uniform',0,4,3,1);
t0 = t_obs(time_steps_index(1));
[T, Z] = ode15s(dynamic,[t_used(1) t_used(end)],z0);
figure
hold on
plot(T,Z,'--','LineWidth',2.5);
% [ax,p1,p2] = plotyy(T,arrayfun(s_in_inter,T),T,arrayfun(d_inter,T));
hold on
plot(T,arrayfun(s_in_inter,T),'LineWidth',2,'LineStyle','-.')
xlabel('\fontsize{12}Time [days]'),...
ylabel('\fontsize{12} Concentration [g/l]'),...
% ylabel(ax(2),'\fontsize{12} Dilution Rate [day]^{-1}'),
% set(p1,'LineWidth',2,'LineStyle','-.'),
% set(p2,'LineWidth',2,'LineStyle',':'),
title(sprintf('Invariant in Time Reactor %s',reactor))
set(gca,'fontsize',15);
% set(ax,'fontsize',15);
legend('z_1(t)','z_2(t)','z_3(t)','s_{in}(t)','Location','best');
%     legend(sprintf('Invariant z(0) = %.2f',z0(1)),sprintf('Invariant z(0) = %.2f',z0(2)),...
%         sprintf('Invariant z(0) = %.2f',z0(3)),'s_{in}(t)','D(t)','Location','best');
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print('Images\\Invariant_for_thesis','-dpng','-r0')
% %% Observer Trajectories
z1 = @(t) interp1(T,Z(:,1),t);
z2 = @(t) interp1(T,Z(:,2),t);
z3 = @(t) interp1(T,Z(:,3),t);
%Observer Used for ECC 2019
% xG1 = max(0,arrayfun(invariant,t_used) - arrayfun(s1_inter,t_used));
% xG2 = max(0,arrayfun(invariant,t_used)...
%     - arrayfun(s1_inter,t_used) - arrayfun(s2_inter,t_used));
xG1 = 1/3*max(0,-2*(arrayfun(s1_inter,t_used)-arrayfun(z1,t_used))...
    +arrayfun(s2_inter,t_used)-arrayfun(z2,t_used)...
    +arrayfun(s3_inter,t_used) -arrayfun(z3,t_used));
xG2 =1/3*max(0,-(arrayfun(s1_inter,t_used)-arrayfun(z1,t_used))...
    -(arrayfun(s2_inter,t_used)-arrayfun(z2,t_used))...
    +2*(arrayfun(s3_inter,t_used) -arrayfun(z3,t_used)));
optim_fun = @(x) norm(biomass(time_steps_index) - x(1)*xG1 - x(2)*xG2);
% [y_opt , val] = fmincon(optim_fun,[yA_ref yB_ref],[],[],[],[],[0 0],[],[]);
% [y_opt , val] = simulannealbnd(optim_fun,[yA_ref yB_ref]);
ms = MultiStart;
options = optimset('MaxFunEvals',2000000,'MaxIter',100000)';
problem = createOptimProblem('fmincon','x0',[yA_ref yB_ref],...
    'objective',optim_fun,'options',options,'lb',[0.1 0.03],'ub',[0.45 0.14]);
tic
[y_opt,val,flag,outpt,allmins] = run(ms,problem,300);
toc
% xG1 = Z - arrayfun(s1ObsInterpol,T);
% xG2 = Z - arrayfun(s1ObsInterpol,T) - arrayfun(s2ObsInterpol,T);
figure
hold on
plot(t_used,y_opt(1)*xG1,'--','LineWidth',2.5);
plot(t_used,y_opt(2)*xG2,'--','LineWidth',2.5);
plot(t_used,y_opt(1)*xG1+y_opt(2)*xG2,'-','LineWidth',2.5);
plot(t_used,biomass(time_steps_index),'*','LineWidth',2.5);
xlabel('\fontsize{12}Time [days]'),...
ylabel('\fontsize{12} Concentration [g/l]'),...
title(sprintf('Observers in Time Reactor %s',reactor))
    set(gca,'fontsize',15);
%     set(ax,'fontsize',15);
    l = legend('$$\hat{x}_{G_1}$$','$$\hat{x}_{G_2}$$',...
        '$$\hat{x}_{G_1}+\hat{x}_{G_2}$$ ','Measured Biomass','Location','best');
    set(l,'Interpreter','latex')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print('Images\\Observer_for_thesis','-dpng','-r0')
% k1 = 0.1;
% dynamic_s1 = @(t,s1) k1*(s1ObsInterpol(t)-s1);
% s1_0 = 0.1;
% [T_s1, s1] = ode15s(dynamic_s1,[0 tObs(end)],s1_0);
% figure
% hold on
% plot(T_s1,s1,'-','LineWidth',2.5);
% plot(tObs,S1);
% xlabel('\fontsize{12}Time [days]'),...
% ylabel('\fontsize{12} Concentration [g/l]')
% k2 = 0.2;
% dynamic_s2 = @(t,s2) k1*(s2ObsInterpol(t)-s2);
% s2_0 = 0.;
% [T_s2, s2] = ode15s(dynamic_s2,[0 tObs(end)],s2_0);
% figure
% hold on
% plot(T_s2,s2,'-','LineWidth',2.5);
% plot(tObs,S2);
% xlabel('\fontsize{12}Time [days]'),...
% ylabel('\fontsize{12} Concentration [g/l]')
