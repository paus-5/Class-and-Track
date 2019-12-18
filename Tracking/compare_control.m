clear
close all
file_names = {'191210_no_noise_Iter_10', '191210_no_noise_2_Iter_10',...
    '191210_no_noise_3_Iter_10','191209_no_noise_Iter_10',...
    '191209_no_noise_try2_Iter_10','191209_no_noise_try2_Iter_10'};
name_plot_out = 'Comparison_different_control';
figure
index_show_OTU = [6 14];
n_show_OTU = length(index_show_OTU);
number_files = length(file_names);
legend_tags_compare = cell((number_files+1)*n_show_OTU,1);
lw = linspace(0.5,3,length(file_names));
index_tags_control = reshape(1:(number_files+1)*n_show_OTU,n_show_OTU,(number_files+1));
hold on
for k=1:number_files
    file_in = file_names{k};
    load(sprintf('MAT_files/%s',file_in));
    colors = mat2cell(hsv(n_show_OTU),ones(1,n_show_OTU),3);
    index_control = 1:length(t_new);
    t_used = t_new(index_control);
    for k2 = 1:n_show_OTU
        legend_tags_compare{index_tags_control(k2,k)} = sprintf('u_{%.0f} \\lambda = %.1g',index_show_OTU(k2),lambda(1));
    end
    control = @(t) tracking_control(lambda,B(t),P(t),x_fun(t),sf_interp(t),n);
    control_eval = arrayfun(control,t_used,'UniformOutput',false);
    control_matrix = 1+(reshape(cell2mat(control_eval),n,length(index_control)));
    p = plot(t_used,control_matrix(index_show_OTU,:),'--','LineWidth',lw(k));
    set(p,{'Color'},colors)
end
   for k2 = 1:n_show_OTU
        legend_tags_compare{index_tags_control(k2,number_files+1)} = sprintf('I_{%.0f}(x(t))',index_show_OTU(k2));
    end
    A_X_data_trajectories = A*Y(:,1:n)';
    p = plot(t_obs,1+A_X_data_trajectories(index_show_OTU,:),'*');
    legend(legend_tags_compare,'fontsize',12,'Location','bestoutside');
    set(p,{'Color'},colors)
    set(gca,'fontsize',15)
    xlabel('\fontsize{15}Time [days]'),...
    ylabel('\fontsize{15}  Control'),...
    title('Comparison of controls and interaction function');
    fig.PaperUnits = 'inches';
	fig.PaperPosition = [0 0 9 4];
    fig.PaperPositionMode = 'manual';
    print(sprintf('Images\\%s',name_plot_out),'-dpng','-r0')