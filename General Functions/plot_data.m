clear
close all
addpath('..\\General Functions')
Reactor = 'B';
[tObs,X,sIn,s1,s2,s3,dilutionRate] = load_biomass(Reactor);
[tOTU,OTU] = load_relative_abundance(Reactor);
subplot(4,1,1),
plot(tObs,X,'LineWidth',3)
    xlabel('\fontsize{12}Time [days]'),...
    ylabel('\fontsize{12}Concentration [g/l]'),...
    title(sprintf('Total Biomass in Time Reactor %s',Reactor))
    set(gca,'fontsize',15);
subplot(4,1,2),
plot(tObs,s1,'LineWidth',3)
    xlabel('\fontsize{12}Time [days]'),...
    ylabel('\fontsize{12}Concentration [g/l]'),...
    title(sprintf('Total ammonium in Time Reactor %s',Reactor))
    set(gca,'fontsize',15);
subplot(4,1,3),
plot(tObs,s2,'LineWidth',3)
    xlabel('\fontsize{12}Time [days]'),...
    ylabel('\fontsize{12}Concentration [g/l]'),...
    title(sprintf('Total nitrite in Time Reactor %s',Reactor))
    set(gca,'fontsize',15);
subplot(4,1,4),
plot(tObs,s3,'LineWidth',3)
    xlabel('\fontsize{12}Time [days]'),...
    ylabel('\fontsize{12}Concentration [g/l]'),...
    title(sprintf('Total nitrate in Time Reactor %s',Reactor))
    set(gca,'fontsize',15);
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
fileOut = sprintf('Data Reactor %s',Reactor);
% print(sprintf('..\\Images\\%s',fileOut),'-dpng','-r0')
figure
subplot(2,1,1),
plot(tObs,sIn,'LineWidth',3)
    xlabel('\fontsize{12}Time [days]'),...
    ylabel('\fontsize{12}Concentration [g/l]'),...
    title(sprintf('Input ammonium in Time Reactor %s',Reactor))
    set(gca,'fontsize',15);
subplot(2,1,2),
plot(tObs,dilutionRate,'LineWidth',3)
    xlabel('\fontsize{12}Time [days]'),...
    ylabel('\fontsize{12} Rate [day]^{-1}'),...
    title(sprintf('Dilution Rate in Time Reactor %s',Reactor))
    set(gca,'fontsize',15);
figure
numberOfOTU = length(OTU(1,:));
 colors = mat2cell(jet(numberOfOTU), ones(1,numberOfOTU), 3);
 OTUPlot = plot(tOTU,OTU,'LineWidth',1.5);
    xlabel('\fontsize{12}Time [days]'),...
    ylabel('\fontsize{12} % of Abundance'),...
    title(sprintf('Relative abundance of OTU in Time Reactor %s',Reactor))
    set(gca,'fontsize',15),
    set(OTUPlot,{'Color'}, colors)
    legend(strsplit(num2str(1:numberOfOTU)),'fontsize',7);