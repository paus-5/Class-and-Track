clear
close all
file_name_in = 'parameters_modified';
load(sprintf('MAT_files\\Operating_Diagram_%s',file_name_in))
figure;
present_zones = unique(zones_OD);
number_of_zones = length(present_zones);
cmap1 = hsv(number_of_zones);
[X,Y] = meshgrid(s_in_vector,D_vector);
surfc(X,Y,zones_OD','EdgeColor','none','LineStyle','none','FaceLighting','phong');
%view(2)
colormap(cmap1)
%     tick_labels = {'Partial Nitrification','Complete Nitrification'};
%     tick_position = [0.2, 0.8]+1;
tick_labels = {'No SE','CN','CN & PN','PN','PN & WO','WO'};
% tick_position = (linspace(1,number_of_zones,number_of_zones))*number_of_zones/(number_of_zones+1);
tick_position = [1,2,3,4,5,6]; 
colorbar('Ticks',tick_position(present_zones),...
         'TickLabels',tick_labels(present_zones))
view(2)
xlabel('S_{in} [g/l]');
ylabel('Dilution Rate [day^{-1}]')
title(sprintf('\\fontsize{20} Operating Diagram'))
set(gca,'fontsize',14)
% fig.PaperPosition = [.25 .25 16 12];
fig.PaperPositionMode = 'auto';
print(sprintf('Images/OD_%s',file_name_in),'-dpng','-r0')
