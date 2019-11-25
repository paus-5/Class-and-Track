clear
close all
name_file = 'no_noise';
name_save_file = 'synthetic_data_noise_2';
load(sprintf('MAT_files\\synthetic_data_%s',name_file))
mean_Y = mean(Y);
noise = zeros(size(Y));
for i=1:length(Y(1,:))
    noise(:,i) = random('unif',0,mean_Y(i)/2,size(noise(:,i)));
end
Y = max(0,Y+noise);
S1 = Y(:,n+1);
S2 = Y(:,n+2);
S3 = Y(:,n+3);
OTU_interp = Y(:,1:n);
save(sprintf('MAT_files\\%s',name_save_file));