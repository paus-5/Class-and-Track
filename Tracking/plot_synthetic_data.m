clear
close all
name_file = 'no_noise';
load(sprintf('MAT_files\\synthetic_data_%s',name_file))
subplot(3,1,1), plot(t_obs, Y(:,1:nA))
subplot(3,1,2), plot(t_obs, Y(:,nA+1:n))
subplot(3,1,3), plot(t_obs, Y(:,n+1:n+3))