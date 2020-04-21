clear
close all
name_file_in = '250309_POC';
load(sprintf('MAT_files\\synthetic_data_%s',name_file_in))
name_file_out = '250309_POC';
[num_row, num_col] = size(Y);
Y_cell = mat2cell(Y',num_col,ones(num_row,1));
growth_AOB =  cellfun(growth1,Y_cell)'.*Y(:,1:nA);
growth_NOB =  cellfun(growth2,Y_cell)'.*Y(:,nA+1:n);
%AOB
numbering_AOB = strsplit(num2str(1:nA));
legend_tags_AOB = strsplit(strcat(sprintf('OTU %s,', numbering_AOB{:})),',');
legend_tags_AOB(end) = [];
colors_AOB =  mat2cell(hsv(nA),ones(1,nA), 3);
figure
AOB_plot = plot(t_obs, Y(:,1:nA),'*','LineWidth',0.7);
xlabel('\fontsize{15}Time [days]')
ylabel('\fontsize{15}  Concentration [g/l]')
title(sprintf('Synthetic data simulation AOB'));
set(gca,'fontsize',15),
set(AOB_plot,{'Color'}, colors_AOB)
legend(legend_tags_AOB,'fontsize',10,'Location','bestoutside');
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',sprintf('%s_AOB_plot',name_file_out)),'-dpng','-r0')
%NOB
numbering_NOB = strsplit(num2str(nA+1:n));
legend_tags_NOB = strsplit(strcat(sprintf('OTU %s,', numbering_NOB{:})),',');
legend_tags_NOB(end) = [];
colors_NOB =  mat2cell(hsv(nB),ones(1,nB), 3);
figure
NOB_plot = plot(t_obs, Y(:,nA+1:n),'*','LineWidth',0.7);
xlabel('\fontsize{15}Time [days]')
ylabel('\fontsize{15}  Concentration [g/l]')
title(sprintf('Synthetic data simulation NOB'));
set(gca,'fontsize',15),
set(NOB_plot,{'Color'}, colors_NOB)
legend(legend_tags_NOB,'fontsize',10,'Location','bestoutside');
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',sprintf('%s_NOB_plot',name_file_out)),'-dpng','-r0')
%%Metabolites
figure
metabolite_plot = plot(t_obs,Y(:,(n+1):end),'*','LineWidth',0.7);
xlabel('\fontsize{15}Time [days]'),...
    ylabel('\fontsize{15}  Concentration [g/l]'),...
    title('Synthetic data simulation metabolites');
set(gca,'fontsize',15),
colors_metabolites =  mat2cell(hsv(3),ones(1,3), 3);
set(metabolite_plot,{'Color'},colors_metabolites);
legend('s_1','s_2','s_3','Location','bestoutside')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',sprintf('%s_metabolites',name_file_out)),'-dpng','-r0')
%%Growth Rate AOB
figure
growth_rate_AOB_plot = plot(t_obs,growth_AOB,'*','LineWidth',0.7);
xlabel('\fontsize{15}Time [days]'),...
    ylabel('\fontsize{15}  Rate [gBiomass/day]'),...
    title('Synthetic data growth rate AOB');
set(gca,'fontsize',15),
set(growth_rate_AOB_plot,{'Color'},colors_AOB);
legend('s_1','s_2','s_3','Location','bestoutside')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',sprintf('%s_growth_rate',name_file_out)),'-dpng','-r0')
%%Growth Rate NOB
figure
growth_rate_NOB_plot = plot(t_obs,growth_NOB,'*','LineWidth',0.7);
xlabel('\fontsize{15}Time [days]'),...
    ylabel('\fontsize{15}  Rate [gBiomass/day]'),...
    title('Synthetic data growth rate NOB');
set(gca,'fontsize',15),
set(growth_rate_NOB_plot,{'Color'},colors_NOB);
legend('s_1','s_2','s_3','Location','bestoutside')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',sprintf('%s_growth_rate',name_file_out)),'-dpng','-r0')




