%Calculating Interactions
clear
close all
file_in = '191211_no_noise_iter_10';
load(sprintf('MAT_files/%s',file_in));
index_control = find(t_new>20 & t_new<80);
% index_control = 1:length(t_new);
t_used = t_new(index_control);
clear diff
time_intervals = diff(t_used);
time_weights_middle = 0.5*time_intervals(2:end) + 0.5*time_intervals(1:end-1);
time_weight = [time_intervals(1); time_weights_middle; time_intervals(end)];
control = @(t) tracking_control(lambda,B(t),P(t),x_fun(t),sf_interp(t),n);
x_trajectories = x_new(index_control,1:n)';
control_eval = arrayfun(control,t_used,'UniformOutput',false);
control_matrix = reshape(cell2mat(control_eval),n,length(index_control));
control_matrix_med_filter = medfilt1(control_matrix',200)';
% optim_fun = @(interactions) norm((control_matrix_med_filter - (interactions*...
%     x_trajectories)*diag(time_weight)),'fro');
optim_fun = @(interactions) norm(control_matrix_med_filter - (interactions*...
    x_trajectories),'fro');
nonlcon = @(interactions) nonlinear_con(interactions,x_trajectories);
options = optimset('MaxFunEvals',2000000,'MaxIter',100000)'; 
ms = MultiStart;
problem = createOptimProblem('fmincon','x0',zeros(n,n),...
    'objective',optim_fun,'nonlcon',nonlcon,'options',options);
tic
[xmin,fmin,flag,outpt,allmins] = run(ms,problem,1);
% tic
% [xmin, fval] = fmincon(optim_fun,A,[],[],[],[],[],[],nonlcon,options);
toc
new_A = xmin;
save(sprintf('mat_Files\\%s_iter_%.0f_interactions',name_file,iter));
figure
hold on
p1 = plot(t_used,A*x_trajectories);
set(p1,{'Color'},colors)
p2 = plot(t_used,new_A*x_trajectories,'--');
set(p2,{'Color'},colors)
p3 = plot(t_used,control_matrix,'*');
set(p3,{'Color'},colors)


