%main: process and filter data to model.
close all
clear
load('\MAT_files\ClassificationAlgorithm_Tracking_290919')
nameFile = '290919_Try1';
nA = round(sum(a));
nB = round(sum(b)); 
n = nA+nB;
index_t_obs = find(t_obs <300);
index_span = time_steps_index(1):index_t_obs(end);
tS = t_obs(index_span); 
s1_interp = @(t) interp1(tS,S1(index_span),t);
s2_interp = @(t) interp1(tS,S2(index_span),t);
s3_interp = @(t) interp1(tS,S3(index_span),t);
index_AOB = find(class_AOB);
index_NOB = find(class_NOB);
biomass_filtered = x(:,[index_AOB; index_NOB]);
min_biomass = min(biomass_filtered(biomass_filtered>0));
biomass_filtered(biomass_filtered == 0) = min_biomass;
z = @(t) interp1(tOTU,biomass_filtered,t)';
kA = yieldsAOB(index_AOB);
kB = yieldsNOB(index_NOB);
muARef = 0.77;
muBRef = 1.07;
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
legendTags = strsplit(strjoin({num2str(1:n),...
    ' s_1 s_2 s_3'}));
legendTagsOTU = strsplit(strjoin({num2str(1:n)}));
tF = tS(end);
t0 = tS(1);
%x0 = [z(t0)' s1Interpol(t0) s2Interpol(t0) s3Interpol(t0)]';
x0 = [0.01*ones(1,n) s1_interp(t0) s2_interp(t0) s3_interp(t0)]'; 
xFun = @(t) x0;
B = @(t) Bx(xFun(t),nA,nB,kA,kB,muA,muB,kSA,kSB,kI);
PtF = F;
systemRicatti = @(t,Plin) ricattiDiff(Plin,xFun(t),Q,nA,nB,kA,kB,muA,...
    muB,kSA,kSB,kI,sInInterpol(t),dInterpol(t),lambda);
tic
[tPlin, Plin] = ode15s(systemRicatti,[tF t0], PtF);
toc
P = @(t) interp1(tPlin,Plin,t)';
sFF = F.*z(tF);
dynamicSF = @(t,sF) sFdiff(sF,xFun(t),z(t),P(t),nA,nB,kA,kB,muA,muB,kSA,...
    kSB,kI,sInInterpol(t),dInterpol(t),lambda);
tic
[tSF, sF] = ode15s(dynamicSF,[tF t0], sFF);
toc
sfInterpol = @(t) interp1(tSF,sF,t)';
dynamic = @(t,X) Ax(X,nA,nB,kA,kB,muA,muB,kSA,kSB,kI,sInInterpol(t),...
    dInterpol(t))*X + Bx(X,nA,nB,kA,kB,muA,muB,kSA,kSB,kI)...
    *trackingControl(lambda,B(t),P(t),X,sfInterpol(t),n);
tic
options = odeset('NonNegative',(1:n+3)');
[tNew,xNew] = ode15s(dynamic,[t0 tF], x0,options);
toc
%%ASRE
xCell = arrayfun(xFun,tNew,'UniformOutput',false);
X = reshape(cell2mat(xCell),length(tNew),n+3);
xOld = X;
csTolerance = 4; %cauchy sequence tolerance
iter = 1;
diff = norm(xOld - xNew);
save(sprintf('MAT_files\\%s_iter_%.0f',nameFile,iter))
while diff > csTolerance
    xFun = @(t) interp1(tNew,xNew,t)';    
    control = @(t) trackingControl(lambda,B(t),P(t),xFun(t),sfInterpol(t),n);
    controlEval = arrayfun(control,tNew,'UniformOutput',false);
    figure
    hold on
    OTUPlot = plot(tNew,xNew(:,1:n),'LineWidth',1.5);
    metabolitePlot = plot(tNew,xNew(:,(n+1):end),'-.','LineWidth',1.5);
    xlabel('\fontsize{12}Time [days]'),...
        ylabel('\fontsize{12}  Concentration'),...
        title(sprintf('Tracking results iteration:%.0f',iter));
    set(gca,'fontsize',15),
    set(OTUPlot,{'Color'}, colors)
    set(metabolitePlot,{'Color'},colors2);
    legend(legendTags,'fontsize',7);
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 9 3];
    fig.PaperPositionMode = 'auto';
    print(sprintf('..\\..\\Images\\%s',sprintf('%s_Trajectory_Iter_%.0f',...
        nameFile,iter)),'-dpng','-r0')
    figure
    controlPlot = plot(tNew,reshape(cell2mat(controlEval),n,...
        length(tNew)),'LineWidth',1.5);
        xlabel('\fontsize{12}Time [days]'),...
        ylabel('\fontsize{12}  v(t)'),...
        title(sprintf('Control for the tracking of z iteration: %.0f',iter))
    set(gca,'fontsize',15),
    set(controlPlot,{'Color'}, colors),
    legend(legendTagsOTU,'fontsize',7);
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 9 3];
    fig.PaperPositionMode = 'auto';
    print(sprintf('..\\..\\Images\\%s',sprintf('%s_Control_Iter_%.0f',...
        nameFile,iter)),'-dpng','-r0')
    B = @(t) Bx(xFun(t),nA,nB,kA,kB,muA,muB,kSA,kSB,kI);
    systemRicatti = @(t,Plin) ricattiDiff(Plin,xFun(t),Q,nA,nB,kA,kB,muA,...
    muB,kSA,kSB,kI,sInInterpol(t),dInterpol(t),lambda);
    tic
    [tPlin, Plin] = ode15s(systemRicatti,[tF t0], PtF);
    toc
    P = @(t) interp1(tPlin,Plin,t)';
    dynamicSF = @(t,sF) sFdiff(sF,xFun(t),z(t),P(t),nA,nB,kA,kB,muA,muB,kSA,...
    kSB,kI,sInInterpol(t),dInterpol(t),lambda);
    tic
    [tSF, sF] = ode15s(dynamicSF,[tF t0], sFF);
    toc
    sfInterpol = @(t) interp1(tSF,sF,t)';
    dynamic = @(t,X) Ax(X,nA,nB,kA,kB,muA,muB,kSA,kSB,kI,sInInterpol(t),...
    dInterpol(t))*X + Bx(X,nA,nB,kA,kB,muA,muB,kSA,kSB,kI)...
    *trackingControl(lambda,B(t),P(t),X,sfInterpol(t),n);
    tic
    [tNew,xNew] = ode15s(dynamic,[t0 tF], x0,options);
    toc
    xOldCell = arrayfun(xFun,tNew,'UniformOutput',false);
    xOld = reshape(cell2mat(xOldCell),length(tNew),n+3);
    iter = iter + 1;
    diff = norm(xOld - xNew);
    save(sprintf('MAT_files\\%s_iter_%.0f',nameFile,iter))
end