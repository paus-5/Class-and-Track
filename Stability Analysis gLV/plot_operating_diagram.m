clear
close all
file_name_in = 'parameters_Dumont';
load(sprintf('MAT_files\\Operating_Diagram_%s',file_name_in))
figure;
number_of_zones = length(unique(zones));
cmap1 = cool(number_of_zones);
[X,Y] = meshgrid(s_in_vector,D_vector);
surfc(X,Y,zones','EdgeColor','none','LineStyle','none','FaceLighting','phong');%view(2)
colormap(cmap1)
if number_of_zones == 2
    tick_labels = {'Partial Nitrification','Complete Nitrification'};
    tick_position = [0.2, 0.8]+1;
else
    tick_labels = {'Partial Nitrification','Complete Nitrification','Both equilibria present'};
    tick_position = [0.25, 0.5, 0.75]+1;
end
colorbar('Ticks',tick_position,...
         'TickLabels',tick_labels)
view(2)
xlabel('S_{in} [g/l]');
ylabel('Dilution Rate [day^{-1}]')
title(sprintf('\\fontsize{20} Operating Diagram'))
set(gca,'fontsize',18)
% fig.PaperPosition = [.25 .25 16 12];
fig.PaperPositionMode = 'auto';
print(sprintf('Images/OD_%s',file_name_in),'-dpng','-r0')
