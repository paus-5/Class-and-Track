%Calculating Interactions
clear
close all
file_in = '221119_no_noise_iter_5';
load(sprintf('MAT_files/%s',file_in));
control = @(t) tracking_control(lambda,B(t),P(t),x_fun(t),sf_interp(t),n);
control_eval = arrayfun(control,t_new,'UniformOutput',false);
optim_fun = @(interactions) sum(norm(reshape(cell2mat(control_eval),n,...
    length(t_new)) - reshape(interactions,n,n)*x_new(:,1:n)'));
options = optimset('MaxFunEvals',2000000,'MaxIter',100000)';
% ms = MultiStart;
% problem = createOptimProblem('fminunc','x0',zeros(n^2,1),...
%     'objective',optim_fun,'options',options);
% tic
% [xmin,fmin,flag,outpt,allmins] = run(ms,problem,50);
% toc
[xmin, fval] = fminunc(optim_fun,zeros(n^2,1),options);
A = reshape(xmin,n,n);
save(sprintf('mat_Files\\%s_iter_%.0f_interactions_ms',name_file,iter));


