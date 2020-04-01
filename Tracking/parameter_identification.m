%Calculating Interactions
clear
close all
file_in = '300309_POC_try2_iter_3';
load(sprintf('MAT_files/%s',file_in));
% index_control = find(t_new>20 & t_new<80);
index_used = 1:(length(t_new)-30);
t_used = t_new(index_used);
x_fun = @(t) interp1(t_new,x_new,t)';  
x_trajectories = x_new(index_used,1:n)';
B = @(t) Bx(x_fun(t),nA,nB,kA,kB,muA,muB,kSA,kSB,kI);
s_vec = x_new(index_used,(n+1):end);
growth_AOB = @(s) growth_monod(s,muA,kSA)';
growth_NOB = @(s) growth_monod(s,muB,kSB)';
growth_cell_AOB = arrayfun(growth_AOB,s_vec(:,1),'UniformOutput',false);
growth_cell_NOB = arrayfun(growth_NOB,s_vec(:,2),'UniformOutput',false);
growth_mat_AOB = cell2mat(growth_cell_AOB);
growth_mat_NOB = cell2mat(growth_cell_NOB);
growth_mat = [growth_mat_AOB growth_mat_NOB];
control = @(t) tracking_control(lambda,B(t),P(t),x_fun(t),sf_interp(t),n);
tic
control_eval = arrayfun(control,t_used,'UniformOutput',false);
toc
control_eval_reshape = reshape(cell2mat(control_eval),n,length(t_used));
growth_times_control = (growth_mat'.*(control_eval_reshape+1))'; %day^-1
mu_fit = zeros(n,1);
kS_fit = zeros(n,1);
new_A = zeros(n,n);
for i=1:n
    if i <= nA
        optim_fun = @(parameter) norm(growth_times_control(:,i)' -growth_monod(s_vec(:,1),parameter(1),parameter(2))'.*(1+parameter(3:end)*x_trajectories));
        lb = [muARef/2 kSARef/100 -inf*ones(1,n)];
        ub = [2*muARef 100*kSARef inf*ones(1,n)];
        p0 = [muARef kSARef zeros(1,n)];
    else
        optim_fun = @(parameter) norm(growth_times_control(:,i)' -growth_monod(s_vec(:,2),parameter(1),parameter(2))'.*(1+parameter(3:end)*x_trajectories));
        lb = [muBRef/2 kSBRef/100 -inf*ones(1,n)];
        ub = [2*muBRef 100*kSBRef inf*ones(1,n)];
        p0 = [muBRef kSBRef zeros(1,n)];
    end
options = optimset('MaxFunEvals',2000000,'MaxIter',100000)'; 
ms = MultiStart;
problem = createOptimProblem('fmincon','x0',p0,...
    'objective',optim_fun,'lb',lb,'ub',ub,'options',options);
tic
[xmin,fmin,flag,outpt,allmins] = run(ms,problem,100);
toc
mu_fit(i) = xmin(1);
kS_fit(i) = xmin(2);
new_A(i,:) = xmin(3:end);
end
save(sprintf('mat_Files\\%s_iter_%.0f_interactions',name_file,iter));
figure
hold on
p1 = plot(t_used,A*x_trajectories);
set(p1,{'Color'},colors)
p2 = plot(t_used,new_A*x_trajectories,'--');
set(p2,{'Color'},colors)
p3 = plot(t_used,control_eval_reshape,'*');
set(p3,{'Color'},colors)
legend('Original A*Tracked trajectories AOB','Original A*Tracked trajectories NOB',...
    'new A*Tracked trajectories AOB','new A*Tracked trajectories NOB',...
    'Control AOB', 'Control NOB')

