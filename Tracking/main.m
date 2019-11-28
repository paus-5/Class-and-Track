%Reference: Cimen et al 2002 Nonlinear optimal tracking control
%main: process and filter data to model.
close all
clear
load('MAT_files\Classification_day183_gurobi')
name_file = '191127_Reactor_A_modified';
nA = round(sum(class_AOB));
nB = round(sum(class_NOB)); 
n = nA+nB;
index_t_obs = find(t_obs <300);
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
muARef = 0.77*1.2;
muBRef = 1.07*1.2;
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
lambda = [0.0001*ones(nA,1); 0.00001*ones(nB,1)];
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
x0 = [0.001*ones(1,n) s1_interp(t0) s2_interp(t0) s3_interp(t0)]';
dynamic = @(t,X) Ax(X,nA,nB,kA,kB,muA,muB,kSA,kSB,kI,s_in_interp(t),...
    d_interp(t))*X ;
options = odeset('NonNegative',(1:n+3)');
[t_new,x_new] = ode15s(dynamic,[t0 tF], x0,options);
x_fun = @(t) interp1(t_new,x_new,t)';
B = @(t) Bx(x_fun(t),nA,nB,kA,kB,muA,muB,kSA,kSB,kI);
PtF = F;
system_ricatti = @(t,P_lin) ricatti_diff(P_lin,x_fun(t),Q,nA,nB,kA,kB,muA,...
    muB,kSA,kSB,kI,s_in_interp(t),d_interp(t),lambda);
tic
[t_Plin, P_lin] = ode15s(system_ricatti,[tF t0], PtF);
toc
P = @(t) interp1(t_Plin,P_lin,t)';
sFF = F.*z(tF);
dynamic_SF = @(t,sF) sF_diff(sF,x_fun(t),z(t),P(t),nA,nB,kA,kB,muA,muB,kSA,...
    kSB,kI,s_in_interp(t),d_interp(t),lambda);
tic
[t_sF, sF] = ode15s(dynamic_SF,[tF t0], sFF);
toc
sf_interp = @(t) interp1(t_sF,sF,t)';
dynamic = @(t,X) Ax(X,nA,nB,kA,kB,muA,muB,kSA,kSB,kI,s_in_interp(t),...
    d_interp(t))*X + Bx(X,nA,nB,kA,kB,muA,muB,kSA,kSB,kI)...
    *tracking_control(lambda,B(t),P(t),X,sf_interp(t),n);
tic
[t_new,x_new] = ode15s(dynamic,[t0 tF], x0,options);
toc
%%SRE
x_cell = arrayfun(x_fun,t_new,'UniformOutput',false);
X = reshape(cell2mat(x_cell),length(t_new),n+3);
x_old = X;
csTolerance = 0.01; %cauchy sequence tolerance
iter = 1;
diff = norm(x_old - x_new);
save(sprintf('MAT_files\\%s_iter_%.0f',name_file,iter))
while diff > csTolerance
    x_fun = @(t) interp1(t_new,x_new,t)';    
%     control = @(t) tracking_control(lambda,B(t),P(t),x_fun(t),sf_interp(t),n);
%     control_eval = arrayfun(control,t_new,'UniformOutput',false);
%     figure
%     hold on
%     OTU_plot = plot(t_new,x_new(:,1:n),'LineWidth',1.5);
%     metabolite_plot = plot(t_new,x_new(:,(n+1):end),'-.','LineWidth',1.5);
%     xlabel('\fontsize{12}Time [days]'),...
%         ylabel('\fontsize{12}  Concentration'),...
%         title(sprintf('Tracking results iteration:%.0f',iter));
%     set(gca,'fontsize',15),
%     set(OTU_plot,{'Color'}, colors)
%     set(metabolite_plot,{'Color'},colors2);
%     legend(legend_tags,'fontsize',7);
%     fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 9 3];
%     fig.PaperPositionMode = 'auto';
%     print(sprintf('Images\\%s',sprintf('%s_Trajectory_Iter_%.0f',...
%         name_file,iter)),'-dpng','-r0')
%     figure 
%     control_plot = plot(t_new,reshape(cell2mat(control_eval),n,...
%         length(t_new)),'LineWidth',1.5,'Visible','on');
%         xlabel('\fontsize{12}Time [days]'),...
%         ylabel('\fontsize{12}  v(t)'),...
%         title(sprintf('Control for the tracking of z iteration: %.0f',iter))
%     set(gca,'fontsize',15),
%     set(control_plot,{'Color'}, colors),
%     legend(legend_tags_OTU,'fontsize',7);
%     fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 9 3];
%     fig.PaperPositionMode = 'auto';
%     print(sprintf('Images\\%s',sprintf('%s_Control_Iter_%.0f',...
%         name_file,iter)),'-dpng','-r0')
    B = @(t) Bx(x_fun(t),nA,nB,kA,kB,muA,muB,kSA,kSB,kI);
    system_ricatti = @(t,P_lin) ricatti_diff(P_lin,x_fun(t),Q,nA,nB,kA,kB,muA,...
        muB,kSA,kSB,kI,s_in_interp(t),d_interp(t),lambda);
    tic
    [t_Plin, P_lin] = ode15s(system_ricatti,[tF t0], PtF);
    toc
    P = @(t) interp1(t_Plin,P_lin,t)';
    sFF = F.*z(tF);
    dynamic_SF = @(t,sF) sF_diff(sF,x_fun(t),z(t),P(t),nA,nB,kA,kB,muA,muB,kSA,...
        kSB,kI,s_in_interp(t),d_interp(t),lambda);
    tic
    [t_sF, sF] = ode15s(dynamic_SF,[tF t0], sFF);
    toc
    sf_interp = @(t) interp1(t_sF,sF,t)';
    dynamic = @(t,X) Ax(X,nA,nB,kA,kB,muA,muB,kSA,kSB,kI,s_in_interp(t),...
        d_interp(t))*X + Bx(X,nA,nB,kA,kB,muA,muB,kSA,kSB,kI)...
        *tracking_control(lambda,B(t),P(t),X,sf_interp(t),n);
    tic
    [t_new,x_new] = ode15s(dynamic,[t0 tF], x0,options);
    toc
    xOldCell = arrayfun(x_fun,t_new,'UniformOutput',false);
    x_old = reshape(cell2mat(xOldCell),length(t_new),n+3);
    iter = iter + 1;
    diff = norm(x_old - x_new);
    save(sprintf('MAT_files\\%s_iter_%.0f',name_file,iter))
end