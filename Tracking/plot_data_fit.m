%plotting
clear
close all
name_file_in = 'MAT_Files\\200407_iter_20';
load(name_file_in);
legend_tags_biomass{1} = 'Data: Sum of species biomass';
legend_tags_biomass{2} = 'Tracking: Sum of species biomass';
index_t_OTU = find(t_OTU<max(tS) & t_OTU>min(tS));
%%Control
x_fun = @(t) interp1(t_new,x_new,t)';
B = @(t) Bx(x_fun(t),nA,nB,kA,kB,muA,muB,kSA,kSB,kI);
s_vec = x_new(:,(n+1):end);
growth_AOB = @(s) growth_monod(s,muA,kSA)';
growth_NOB = @(s) growth_monod(s,muB,kSB)';
growth_cell_AOB = arrayfun(growth_AOB,s_vec(:,1),'UniformOutput',false);
growth_cell_NOB = arrayfun(growth_NOB,s_vec(:,2),'UniformOutput',false);
growth_mat_AOB = cell2mat(growth_cell_AOB);
growth_mat_NOB = cell2mat(growth_cell_NOB);
growth_mat = [growth_mat_AOB growth_mat_NOB];
control = @(t) tracking_control(lambda,B(t),P(t),x_fun(t),sf_interp(t),n);
tic
control_eval = arrayfun(control,t_new,'UniformOutput',false);
toc
control_eval_reshape = reshape(cell2mat(control_eval),n,length(t_new));
growth_rate = growth_mat'.*x_new(:,1:n)';
growth_times_control = growth_mat'.*(control_eval_reshape+1).*x_new(:,1:n)';
%%Total Biomass
figure
hold on
total_biomass_plot = plot(t_OTU(index_t_OTU),sum(OTU_interp(index_t_OTU,:),2),'k*','LineWidth',1.5);
total_biomass_plot2 = plot(t_new,sum(x_new(:,1:n),2),'k.-','LineWidth',1.5);
xlabel('\fontsize{15}Time [days]'),...
    ylabel('\fontsize{15}  Concentration [g/l]'),...
    title(sprintf('Tracking results total biomass iteration: %.0f',iter));
set(gca,'fontsize',15),
legend(legend_tags_biomass,'fontsize',12);
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',sprintf('%s_iter_%.0f_Biomass',...
    name_file,iter)),'-dpng','-r0')
%%AOB 
number_of_AOB_plots = floor(nA/10) + logical(rem(nA,10));
for k = 1:number_of_AOB_plots
    OTU_plot_index = (10*(k-1)+1):(10*(k-1)+10);
    if k == number_of_AOB_plots && logical(rem(nA,10)) == 1;
        OTU_plot_index = (10*(k-1)+1):(10*(k-1)+rem(nA,10));
    end
    numbering_AOB = strsplit(num2str(OTU_plot_index));
    legend_tags_AOB = strsplit(strcat(sprintf('OTU %s (tracking),', numbering_AOB{:}),...
    sprintf(',OTU %s (data)', numbering_AOB{:})),',');
    colors_AOB =  mat2cell(hsv(length(OTU_plot_index)),...
        ones(1,length(OTU_plot_index)), 3);
    figure
    hold on
    AOB_plot = plot(t_new,x_new(:,OTU_plot_index),'LineWidth',1.5);
    AOB_data_plot = plot(t_OTU(index_t_OTU),biomass_filtered(index_t_OTU,OTU_plot_index),'*','LineWidth',0.7);
    xlabel('\fontsize{15}Time [days]'),...
        ylabel('\fontsize{15}  Concentration [g/l]'),...
        title(sprintf('Tracking results AOB iteration: %.0f',iter));
    set(gca,'fontsize',15),
    set(AOB_plot,{'Color'}, colors_AOB)
    set(AOB_data_plot,{'Color'}, colors_AOB)
    legend(legend_tags_AOB,'fontsize',10,'Location','bestoutside');
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 9 3];
    fig.PaperPositionMode = 'auto';
    print(sprintf('Images\\%s',sprintf('%s_iter_%.0f_AOB_plot_%.0f',...
        name_file,iter,k)),'-dpng','-r0')
%%Control
    legend_tags_AOB_control = strsplit(sprintf('Control OTU %s,', numbering_AOB{:}),',');
    legend_tags_AOB_control(end) = [];
    figure
    control_plot_AOB = plot(t_new,control_eval_reshape(OTU_plot_index,:)+1,'LineWidth',1.5);
    xlabel('\fontsize{12}Time [days]'),...
        ylabel('\fontsize{12}  u(t)'),...
        title(sprintf('Control AOB iter: %.0f',iter))
    set(gca,'fontsize',15),
    set(control_plot_AOB,{'Color'}, colors_AOB),
    legend(legend_tags_AOB_control,'fontsize',10,'Location','bestoutside');
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 9 3];
    fig.PaperPositionMode = 'auto';
    print(sprintf('Images\\%s',sprintf('%s_iter_%.0f_Control_AOB_plot_%.0f',...
        name_file,iter,k)),'-dpng','-r0')
%%f(s)xu
    legend_tags_AOB_growth = strsplit(sprintf('Growth OTU %s,',  numbering_AOB{:}),',');
    legend_tags_AOB_growth(end) = [];
    figure
    growth_plot_control_AOB = plot(t_new,growth_times_control(OTU_plot_index,:),'LineWidth',1.5);
    xlabel('\fontsize{12}Time [days]'),...
        ylabel('\fontsize{12} Growth rate [g Biomass \cdot day^{-1}]'),...
        title(sprintf('Growth rate AOB f_i(s_1(t))u_i(t)x_i(t) iter: %.0f',iter))
    set(gca,'fontsize',15),
    set(growth_plot_control_AOB,{'Color'}, colors_AOB),
    legend(legend_tags_AOB_growth,'fontsize',10,'Location','bestoutside');
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 9 3];
    fig.PaperPositionMode = 'auto';
    print(sprintf('Images\\%s',sprintf('%s_iter_%.0f_growth_control_AOB_plot_%.0f',...
        name_file,iter,k)),'-dpng','-r0')
%%f(s)x
    figure
    growth_plot_AOB = plot(t_new,growth_rate(OTU_plot_index,:),'LineWidth',1.5);
    xlabel('\fontsize{12}Time [days]'),...
        ylabel('\fontsize{12} Growth rate [g Biomass \cdot day^{-1}]'),...
        title(sprintf('Growth rate AOB f_i(s_1(t))x_i(t) iter: %.0f',iter))
    set(gca,'fontsize',15),
    set(growth_plot_AOB,{'Color'}, colors_AOB),
    legend(legend_tags_AOB_growth,'fontsize',10,'Location','bestoutside');
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 9 3];
    fig.PaperPositionMode = 'auto';
    print(sprintf('Images\\%s',sprintf('%s_iter_%.0f_growth_AOB_plot_%.0f',...
        name_file,iter,k)),'-dpng','-r0')
end
%%NOB
number_of_NOB_plots = floor(nB/10) + logical(rem(nB,10));
for k = 1:number_of_NOB_plots
    OTU_plot_index = (10*(k-1)+1):(10*(k-1)+10);
    if k == number_of_NOB_plots && logical(rem(nB,10)) == 1;
        OTU_plot_index = (10*(k-1)+1):(10*(k-1)+rem(nB,10));
    end
    OTU_plot_index = OTU_plot_index + nA;
    numbering_NOB = strsplit(num2str(OTU_plot_index)) ;
    legend_tags_NOB = strsplit(strcat(sprintf('OTU %s (tracking),', numbering_NOB{:}),...
    sprintf(',OTU %s (data)', numbering_NOB{:})),',');
    colors_NOB =  mat2cell(hsv(length(OTU_plot_index)),...
        ones(1,length(OTU_plot_index)), 3);
    figure
    hold on
    NOB_plot = plot(t_new,x_new(:,OTU_plot_index),'LineWidth',1.5);
    NOB_data_plot = plot(t_OTU(index_t_OTU),biomass_filtered(index_t_OTU,OTU_plot_index),'*','LineWidth',0.7);
    xlabel('\fontsize{15}Time [days]'),...
        ylabel('\fontsize{15}  Concentration [g/l]'),...
        title(sprintf('Tracking results NOB iteration: %.0f',iter));
    set(gca,'fontsize',15),
    set(NOB_plot,{'Color'}, colors_NOB)
    set(NOB_data_plot,{'Color'}, colors_NOB)
    legend(legend_tags_NOB,'fontsize',10,'Location','bestoutside');
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 9 3];
    fig.PaperPositionMode = 'auto';
    print(sprintf('Images\\%s',sprintf('%s_iter_%.0f_NOB_plot_%.0f',...
        name_file,iter,k)),'-dpng','-r0')
%%Control
    legend_tags_NOB_control = strsplit(sprintf('Control OTU %s,', numbering_NOB{:}),',');
    legend_tags_NOB_control(end) = [];
    figure
    control_plot_NOB = plot(t_new,control_eval_reshape(OTU_plot_index,:)+1,'LineWidth',1.5);
    xlabel('\fontsize{12}Time [days]'),...
        ylabel('\fontsize{12}  u(t)'),...
        title(sprintf('Control NOB iter: %.0f',iter))
    set(gca,'fontsize',15),
    set(control_plot_NOB,{'Color'}, colors_NOB),
    legend(legend_tags_NOB_control,'fontsize',10,'Location','bestoutside');
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 9 3];
    fig.PaperPositionMode = 'auto';
    print(sprintf('Images\\%s',sprintf('%s_iter_%.0f_Control_NOB_plot_%.0f',...
        name_file,iter,k)),'-dpng','-r0')
%%f(s)ux
    legend_tags_NOB_growth = strsplit(sprintf('Growth OTU %s,', numbering_NOB{:}),',');
    legend_tags_NOB_growth(end) = [];
    figure
    growth_control_plot_NOB = plot(t_new,growth_times_control(OTU_plot_index,:),'LineWidth',1.5);
    xlabel('\fontsize{12}Time [days]'),...
        ylabel('\fontsize{12}  [g Biomass \cdot day^{-1}]'),...
        title(sprintf('Growth rate NOB f_i(s_2(t))u_i(t)x_i(t) iter: %.0f',iter))
    set(gca,'fontsize',15),
    set(growth_control_plot_NOB,{'Color'}, colors_NOB),
    legend(legend_tags_NOB_growth,'fontsize',10,'Location','bestoutside');
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 9 3];
    fig.PaperPositionMode = 'auto';
    print(sprintf('Images\\%s',sprintf('%s_iter_%.0f_growth_control_NOB_plot_%.0f',...
        name_file,iter,k)),'-dpng','-r0')
%%f(s)x
    figure
    growth_plot_NOB = plot(t_new,growth_rate(OTU_plot_index,:),'LineWidth',1.5);
    xlabel('\fontsize{12}Time [days]'),...
        ylabel('\fontsize{12} Growth rate [g Biomass \cdot day^{-1}]'),...
        title(sprintf('Growth NOB f_i(s_2(t))x_i(t) iter: %.0f',iter))
    set(gca,'fontsize',15),
    set(growth_plot_NOB,{'Color'}, colors_NOB),
    legend(legend_tags_NOB_growth,'fontsize',10,'Location','bestoutside');
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 9 3];
    fig.PaperPositionMode = 'auto';
    print(sprintf('Images\\%s',sprintf('%s_iter_%.0f_growth_NOB_plot_%.0f',...
        name_file,iter,k)),'-dpng','-r0')
end
%%Metabolites
figure
hold on
metabolite_plot = plot(t_new,x_new(:,(n+1):end),'LineWidth',1.5);
data_metabolite_plot = plot(tS,[S1(index_span) S2(index_span) S3(index_span)],'*');
xlabel('\fontsize{15}Time [days]'),...
    ylabel('\fontsize{15}  Concentration [g/l]'),...
    title(sprintf('Tracking results iteration: %.0f',iter));
set(gca,'fontsize',15),
set(metabolite_plot,{'Color'},colors2);
set(data_metabolite_plot,{'Color'},colors2);
legend('s_1 tracking','s_2 tracking','s_3 tracking','s_1 data','s_2 data','s_3 data','Location','bestoutside')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',sprintf('%s_iter_%.0f_metabolites',...
    name_file,iter)),'-dpng','-r0')