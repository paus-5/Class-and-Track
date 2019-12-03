%plotting
clear
close all
name_file_in = 'MAT_Files\\191130_Reactor_A_iter_10';
load(name_file_in);
numbering_AOB = strsplit(num2str(1:nA)) ;
legend_tags_AOB = strsplit(strcat(sprintf('OTU %s (tracking),', numbering_AOB{:}),...
    sprintf(',OTU %s (data)', numbering_AOB{:})),',');
numbering_NOB = strsplit(num2str(nA+1:n)) ;
legend_tags_NOB = strsplit(strcat(sprintf('OTU %s (tracking),', numbering_NOB{:}),...
    sprintf(',OTU %s (data)', numbering_NOB{:})),',');
legend_tags_biomass{1} = 'Data: Sum of species biomass';
legend_tags_biomass{2} = 'Tracking: Sum of species biomass';
legend_tags_AOB_control = strsplit(sprintf('Control OTU %s,', numbering_AOB{:}),',');
legend_tags_AOB_control(end) = [];
legend_tags_NOB_control = strsplit(sprintf('Control OTU %s,', numbering_NOB{:}),',');
legend_tags_NOB_control(end) = [];
colors_OTU = mat2cell([autumn(nA); winter(nB)], ones(1,n), 3);
colors_AOB = mat2cell(autumn(nA), ones(1,nA), 3);
colors_NOB = mat2cell(winter(nB), ones(1,nB), 3);
colors2 = mat2cell(hsv(3), ones(1,3),3);
% xFun = @(t) interp1(tNew,xNew,t)';    
% control = @(t) trackingControl(lambda,B(t),P(t),xFun(t),sfInterpol(t),n);
% controlEval = arrayfun(control,tNew,'UniformOutput',false);
index_t_OTU = find(t_OTU<max(tS) & t_OTU>min(tS));
%%Total Biomass
figure
hold on
total_biomass_plot = plot(t_OTU(index_t_OTU),sum(OTU_interp(index_t_OTU,:),2),'k*','LineWidth',1.5);
total_biomass_plot2 = plot(t_new,sum(x_new(:,1:n),2),'k.-','LineWidth',1.5);
xlabel('\fontsize{15}Time [days]'),...
    ylabel('\fontsize{15}  Concentration [g/l]'),...
    title(sprintf('Tracking results total biomass iteration:%.0f',iter));
%     title(sprintf('Simulation v(t) = -1'));
set(gca,'fontsize',15),
legend(legend_tags_biomass,'fontsize',12);
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',sprintf('%s_Biomass_Iter_%.0f',...
    name_file,iter)),'-dpng','-r0')
%%AOB
figure
hold on
AOB_plot = plot(t_new,x_new(:,1:nA),'LineWidth',1.5);
AOB_data_plot = plot(t_OTU(index_t_OTU),biomass_filtered(index_t_OTU,1:nA),'*','LineWidth',0.7);
xlabel('\fontsize{15}Time [days]'),...
    ylabel('\fontsize{15}  Concentration [g/l]'),...
    title(sprintf('Tracking results AOB iteration:%.0f',iter));
%     title(sprintf('Simulation v(t) = -1'));
set(gca,'fontsize',15),
set(AOB_plot,{'Color'}, colors_AOB)
set(AOB_data_plot,{'Color'}, colors_AOB) 
legend(legend_tags_AOB,'fontsize',10,'Location','bestoutside');
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',sprintf('%s_AOB_Iter_%.0f',...
    name_file,iter)),'-dpng','-r0')
%%NOB
figure
hold on
NOB_plot = plot(t_new,x_new(:,nA+1:n),'LineWidth',1.5);
NOB_data_plot = plot(t_OTU(index_t_OTU),biomass_filtered(index_t_OTU,nA+1:end),'*','LineWidth',0.7);
xlabel('\fontsize{15}Time [days]'),...
    ylabel('\fontsize{15}  Concentration [g/l]'),...
    title(sprintf('Tracking results NOB iteration:%.0f',iter));
%     title(sprintf('Simulation v(t) = -1'));
set(gca,'fontsize',15),
set(NOB_plot,{'Color'}, colors_NOB)
set(NOB_data_plot,{'Color'}, colors_NOB) 
legend(legend_tags_NOB,'fontsize',10,'Location','bestoutside');
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',sprintf('%s_NOB_Iter_%.0f',...
    name_file,iter)),'-dpng','-r0')
%%Metabolites
figure
hold on
metabolite_plot = plot(t_new,x_new(:,(n+1):end),'LineWidth',1.5);
data_metabolite_plot = plot(tS,[S1(index_span) S2(index_span) S3(index_span)],'*');
xlabel('\fontsize{15}Time [days]'),...
    ylabel('\fontsize{15}  Concentration [g/l]'),...
    title(sprintf('Tracking results iteration:%.0f',iter));
set(gca,'fontsize',15),
set(metabolite_plot,{'Color'},colors2);
set(data_metabolite_plot,{'Color'},colors2);
legend('s_1 tracking','s_2 tracking','s_3 tracking','s_1 data','s_2 data','s_3 data','Location','bestoutside')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',sprintf('%s_metabolites_Iter_%.0f',...
    name_file,iter)),'-dpng','-r0')
%%Control
x_fun = @(t) interp1(t_new,x_new,t)';
B = @(t) Bx(x_fun(t),nA,nB,kA,kB,muA,muB,kSA,kSB,kI);
control = @(t) tracking_control(lambda,B(t),P(t),x_fun(t),sf_interp(t),n);
control_eval = arrayfun(control,t_new,'UniformOutput',false);
control_eval_reshape = reshape(cell2mat(control_eval),n,length(t_new));
%%Control AOB
figure
control_plot_AOB = plot(t_new,control_eval_reshape(1:nA,:)+1,'LineWidth',1.5);
xlabel('\fontsize{12}Time [days]'),...
    ylabel('\fontsize{12}  u(t)'),...
    title(sprintf('Control for AOB iter: %.0f',iter))
set(gca,'fontsize',15),
set(control_plot_AOB,{'Color'}, colors_AOB),
legend(legend_tags_AOB_control,'fontsize',10,'Location','bestoutside');
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',sprintf('%s_Control_AOB_Iter_%.0f',...
    name_file,iter)),'-dpng','-r0')
%%Control NOB
figure
control_plot_NOB = plot(t_new,control_eval_reshape(nA+1:end,:)+1,'LineWidth',1.5);
xlabel('\fontsize{12}Time [days]'),...
    ylabel('\fontsize{12}  u(t)'),...
    title(sprintf('Control for NOB iter: %.0f',iter))
set(gca,'fontsize',15),
set(control_plot_NOB,{'Color'}, colors_NOB),
legend(legend_tags_NOB_control,'fontsize',10,'Location','bestoutside');
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',sprintf('%s_Control_NOB_Iter_%.0f',...
    name_file,iter)),'-dpng','-r0')