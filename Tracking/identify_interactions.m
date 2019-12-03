%Calculating Interactions
clear
close all
file_in = '221119_no_noise_iter_5';
load(sprintf('MAT_files/%s',file_in));
index_control = find(t_new>20 & t_new<130);
% index_control = 1:length(t_new);
t_used = t_new(index_control);
clear diff
time_intervals = diff(t_used);
time_weights_middle = 0.5*time_intervals(2:end) + 0.5*time_intervals(1:end-1);
time_weight = [time_intervals(1); time_weights_middle; time_intervals(end)];
control = @(t) tracking_control(lambda,B(t),P(t),x_fun(t),sf_interp(t),n);
x_trajectories = x_new(index_control,1:n)';
control_eval = arrayfun(control,t_used,'UniformOutput',false);
optim_fun = @(interactions) norm((reshape(cell2mat(control_eval),n,...
    length(index_control)) - interactions*...
    x_trajectories)*diag(time_weight),'fro');
options = optimset('MaxFunEvals',2000000,'MaxIter',100000)'; 
nonlcon = @(interactions) nonlinear_con(interactions,x_trajectories);
ms = MultiStart;
problem = createOptimProblem('fmincon','x0',zeros(n,n),...
    'objective',optim_fun,'nonlcon',nonlcon,'options',options);
tic
[xmin,fmin,flag,outpt,allmins] = run(ms,problem,1);
toc
% tic
% [xmin, fval] = fmincon(optim_fun,zeros(n,n),[],[],[],[],[],[],nonlcon,options);
% % [xmin, fval] = fminunc(optim_fun,zeros(n,n),options);
% toc
A = xmin;
save(sprintf('mat_Files\\%s_iter_%.0f_interactions',name_file,iter));


