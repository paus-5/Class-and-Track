clear
close all
file_in = '191211_no_noise_iter_10';
load(sprintf('MAT_files/%s',file_in));
index_control = 1:length(t_new);
t_used = t_new(index_control);
clear diff
time_intervals = diff(t_used);
time_weights_middle = 0.5*time_intervals(2:end) + 0.5*time_intervals(1:end-1);
time_weight = [time_intervals(1); time_weights_middle; time_intervals(end)];
colors = mat2cell(hsv(n),ones(1,n),3);
control = @(t) tracking_control(lambda,B(t),P(t),x_fun(t),sf_interp(t),n);
control_eval = arrayfun(control,t_used,'UniformOutput',false);
control_matrix = (reshape(cell2mat(control_eval),n,length(index_control)));
control_matrix_med_filter = medfilt1(control_matrix',300);
feed_f = @(t) feed_forward(lambda,B(t),P(t),x_fun(t),n);
feed_f_eval =  arrayfun(feed_f,t_used,'UniformOutput',false);
feed_f_matrix = (reshape(cell2mat(feed_f_eval),n,length(index_control)));
x_trajectories = x_new(index_control,1:n)';
data_trajectories = Y(:,1:n)';
figure
hold on
p1 = plot(t_used,A*x_trajectories);
set(p1,{'Color'},colors)
p2 = plot(t_used,control_matrix,'--');
set(p2,{'Color'},colors)
% p3 = plot(t_obs,A*data_trajectories,'*');
% set(p3,{'Color'},colors)
p4 = plot(t_used,control_matrix_med_filter,'.');
set(p4,{'Color'},colors)