clear 
close all
file_in = '221119_no_noise_iter_5_interactions';
load(sprintf('MAT_files/%s',file_in));%Simulation and plotting
numbering_AOB = strsplit(num2str(1:nA)) ;
legend_tags_AOB = strsplit(strcat(sprintf('OTU %s (model),', numbering_AOB{:}),...
    sprintf(',OTU %s (data)', numbering_AOB{:})),',');
numbering_NOB = strsplit(num2str(nA+1:n)) ;
legend_tags_NOB = strsplit(strcat(sprintf('OTU %s (model),', numbering_NOB{:}),...
    sprintf(',OTU %s (data)', numbering_NOB{:})),',');
colors_AOB = mat2cell(autumn(nA), ones(1,nA), 3);
colors_NOB = mat2cell(winter(nB), ones(1,nB), 3);
colors_both = mat2cell([autumn(nA); winter(nB)], ones(1,n), 3);
numbering_diff = strsplit(num2str(1:n));
legend_tags_diff = strsplit(sprintf('Difference %s,', numbering_diff{:}),',');
legend_tags_diff(end) = [];
%Quality of the Fit
fit = A*x_new(index_control,1:n)'+1;
control_matrix = reshape(cell2mat(control_eval),n,length(index_control))+1;
diff = abs(control_matrix - fit);
figure
diff_plot = plot(t_used,diff,'LineWidth',1.5);
xlabel('\fontsize{12}Time [days]'),...
    ylabel('\fontsize{12}  Absolute difference between Ax+1 and u(t)'),...
    title('Difference of LV and Control');
%     title(sprintf('Simulation v(t) = -1'));
set(gca,'fontsize',15),
set(diff_plot,{'Color'}, colors_both)
legend(legend_tags_diff,'fontsize',7,'Location','bestoutside');
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',sprintf('%s_difference_fit_%.0f_LV',...
    file_in,iter)),'-dpng','-r0')
%Simulation
muARef = 0.77;
muBRef = 1.07;
kSARef = 0.7e-1;
kSBRef = 1.4982e-1;
growth1 = @(X) (1+A(1:nA,:)*X(1:n))*growth_monod(X(n+1),muARef,kSARef);
growth2 = @(X) (1+A(nA+1:n,:)*X(1:n))*growth_monod(X(n+2),muBRef,kSBRef);
dynamic_inter = @(t,X) [diag(growth1(X) - d_interp(t))*X(1:nA);
   diag(growth2(X) - d_interp(t))*X(nA+1:n);
   (s_in_interp(t) - X(n+1))*d_interp(t) - kA'*diag(growth1(X))*X(1:nA);
    -X(n+2)*d_interp(t) + kA'*diag(growth1(X))*X(1:nA)-kB'*diag(growth2(X))*X(nA+1:n);
    -X(n+3)*d_interp(t) + kB'*diag(growth2(X))*X(nA+1:n);];
t0_interactions = tS(1);
tF_interactions = tS(end);
x0_interactions = x_new(1,:);
ode_options = odeset('Nonnegative',1:(n+3));
[T_inter, Y_inter] = ode45(dynamic_inter,[t0_interactions tF_interactions],x0_interactions,ode_options);
indexTSOTU = find(t_OTU<tF_interactions & t_OTU>t0_interactions);
%AOB
figure
hold on
AOB_plot = plot(T_inter,Y_inter(:,1:nA),'LineWidth',1.5);
dataAOB_plot = plot(t_OTU(indexTSOTU),biomass_filtered(indexTSOTU,1:nA),'*','LineWidth',0.7);
xlabel('\fontsize{12}Time [days]'),...
    ylabel('\fontsize{12}  Concentration [g/l]'),...
    title('AOB abundance simulation of LV model');
set(gca,'fontsize',15),
set(AOB_plot,{'Color'}, colors_AOB)
set(dataAOB_plot,{'Color'}, colors_AOB)
legend(legend_tags_AOB,'fontsize',7,'Location','bestoutside');
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',sprintf('%s_AOB_Iter_%.0f_LV',...
    file_in,iter)),'-dpng','-r0')
%NOB
figure
hold on
NOB_plot = plot(T_inter,Y_inter(:,nA+1:n),'LineWidth',1.5);
dataNOB_plot = plot(t_OTU(indexTSOTU),biomass_filtered(indexTSOTU,nA+1:n),'*','LineWidth',0.7);
xlabel('\fontsize{12}Time [days]'),...
    ylabel('\fontsize{12}  Concentration [g/l]'),...
    title('NOB abundance simulation of LV model');
%     title(sprintf('Simulation v(t) = -1'));
set(gca,'fontsize',15),
set(NOB_plot,{'Color'}, colors_NOB)
set(dataNOB_plot,{'Color'}, colors_NOB)
legend(legend_tags_NOB,'fontsize',7,'Location','bestoutside');
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',sprintf('%s_NOB_Iter_%.0f_LV',...
    file_in,iter)),'-dpng','-r0')
%Metabolites
figure
hold on
metabolitePlot = plot(T_inter,Y_inter(:,(n+1):end),'LineWidth',1.5);
% dataMetabolitePlot = plot(tS,[S1(indexSpan) S2(indexSpan) S3(indexSpan)],'--');
i_0 =index_span(1);
i_end = index_span(end);
dataMetabolitePlot = plot(t_obs(i_0:i_end),[S1(i_0:i_end) S2(i_0:i_end) S3(i_0:i_end)],'*');
xlabel('\fontsize{12}Time [days]'),...
    ylabel('\fontsize{12}  Concentration [g/l]'),...
    title('Metabolites simulation of LV model');
%     title(sprintf('Simulation v(t)=-1'));
set(gca,'fontsize',15),
set(metabolitePlot,{'Color'},colors2);
set(dataMetabolitePlot,{'Color'},colors2);
legend('s_1 predicted','s_2 predicted','s_3 predicted','s_1 data','s_2 data','s_3 data')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 3];
fig.PaperPositionMode = 'auto';
print(sprintf('Images\\%s',sprintf('%s_metabolites_Iter_%.0f_LV',...
    file_in,iter)),'-dpng','-r0')

