clear
close all
addpath('..\\..\\General Functions')
Reactor = 'A';
[tObs,X,Sin,S1,S2,S3,DilutionRate] = loadBiomass(Reactor);
[tEspeces,especes] = loadRelativeAbundance(Reactor);
dInterpol = @(t) interp1(tObs,DilutionRate,t);
sInInterpol = @(t) interp1(tObs,Sin,t);
s1ObsInterpol = @(t) interp1(tObs,S1,t);
s2ObsInterpol = @(t) interp1(tObs,S2,t);
% sInInterpol = @(t) max(Sin)/2*(1+sin(1/10*t));
dynamic = @(t,z) -1*dInterpol(t)*(z - sInInterpol(t));
% z0 = random('Uniform',0,4,3,1);
z0 = [0.3 1 2.9]';
[T, Z] = ode15s(dynamic,[0 tObs(end)],z0(1));
[T2, Z2] = ode15s(dynamic,[0 tObs(end)],z0(2));
[T3, Z3] = ode15s(dynamic,[0 tObs(end)],z0(3));
observer = @(t) interp1(T,Z,t);
tolerance = 0.01;
difference = @(t) abs((observer(t)-sInInterpol(t)))/min(observer(t),sInInterpol(t));
diffEval = arrayfun(difference,tEspeces);
indexOfValidity = find(diffEval< tolerance);
figure
hold on
plot(T,Z,'--','LineWidth',2.5);
plot(T2,Z2,'--','LineWidth',2.5);
plot(T3,Z3,'--','LineWidth',2.5);
[ax,p1,p2] = plotyy(tObs,arrayfun(sInInterpol,tObs),tObs,arrayfun(dInterpol,tObs));
xlabel('\fontsize{12}Time [days]'),...
ylabel('\fontsize{12} Concentration [g/l]'),...
ylabel(ax(2),'\fontsize{12} Dilution Rate [day]^{-1}'),
set(p1,'LineWidth',2,'LineStyle','-.'),
set(p2,'LineWidth',2,'LineStyle',':'),
title(sprintf('Invariant in Time Reactor %s',Reactor))
    set(gca,'fontsize',15);
    set(ax,'fontsize',15);
    legend(sprintf('Invariant z(0) = %.2f',z0(1)),sprintf('Invariant z(0) = %.2f',z0(2)),...
        sprintf('Invariant z(0) = %.2f',z0(3)),'s_{in}(t)','D(t)','Location','best');
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print('..\\..\\Images\\InvariantZ','-dpng','-r0')
% %% Observer Trajectories
xG1 = max(0,Z - arrayfun(s1ObsInterpol,T));
xG2 = max(0,Z - arrayfun(s1ObsInterpol,T) - arrayfun(s2ObsInterpol,T));
% xG1 = Z - arrayfun(s1ObsInterpol,T);
% xG2 = Z - arrayfun(s1ObsInterpol,T) - arrayfun(s2ObsInterpol,T);
figure
hold on
yA_ref = 0.251;
yB_ref = 0.062;
plot(T,yA_ref*xG1,'--','LineWidth',2.5);
plot(T,yB_ref*xG2,'--','LineWidth',2.5);
plot(T,yA_ref*xG1+yB_ref*xG2,'-','LineWidth',2.5);
plot(tObs,X,'*','LineWidth',2.5);
xlabel('\fontsize{12}Time [days]'),...
ylabel('\fontsize{12} Concentration [g/l]'),...
title(sprintf('Observers in Time Reactor %s',Reactor))
    set(gca,'fontsize',15);
    set(ax,'fontsize',15);
    l = legend('$$y^{A}_{ref}\hat{x}_{G_1}$$','$$y^{B}_{ref}\hat{x}_{G_2}$$',...
        '$$y^A_{ref}\hat{x}_{G_1}+y^B_{ref}\hat{x}_{G_2}$$ ','Measured Biomass');
    set(l,'Interpreter','latex')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print('..\\..\\Images\\ObserverXG1XG2','-dpng','-r0')
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

