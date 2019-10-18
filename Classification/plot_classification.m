clear
close all
figure
id = 'all_data_intlinprog';
file_name = sprintf('classification_%s',id);
load(sprintf('MAT_files\\%s',file_name));
file_out = sprintf('Classification OTU %s',id);
position_plot = [2/5 0.5 3/5];
classification = (class_AOB==1)*position_plot(1) + (class_NOB==1)*position_plot(3)+...
    (class_AOB==0 & class_NOB==0)*position_plot(2);
mean_AOB = 100*sum(mean_abundance(class_AOB==1));
mean_NOB = 100*sum(mean_abundance(class_NOB==1));
mean_notClassified = 100*sum(mean_abundance(class_AOB==0 & class_NOB==0));
hold on
spacing = 1:2:number_OTU*2;
add_spacing = rem(1:number_OTU,2);
plot(spacing,classification ,'d','MarkerSize',11,'MarkerEdgeColor','k','MarkerFaceColor','r'),
% text(-1,3/5, sprintf('%.02f %% of mean Total Biomass',mean_AOB));
xlabel('OTU'),...
    title(sprintf('\\fontsize{18} Functional assignment of OTU Reactor %s',reactor)),
legend(sprintf('AOB: %.02f %%, NOB: %.02f %% of Biomass',mean_AOB,mean_NOB))
ax = gca;
axis([0 number_OTU*2+1 1/3 2/3])
ax.XTick = 1:4:number_OTU*2;
ax.XTickLabel = 1:2:number_OTU;
ax.YTick = position_plot;
ax.YTickLabel = {'AOB','Not Classified','NOB'} ;
grid on
ax.LineWidth = 2.5;
set(gca,'fontsize',13)
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',file_out),'-dpng','-r0')
