clear
close all
file_name_in = 'synthetic_data_no_noise_200121';
load(sprintf('MAT_files\\Operating_Diagram_%s',file_name_in))
figure;
cmap1 = hsv(length(map_zones));
[X,Y] = meshgrid(s_in_vector,D_vector);
surfc(X,Y,zones','EdgeColor','none','LineStyle','none','FaceLighting','phong');%view(2)
colormap(cmap1)
tickLabels = cell(1,length(map_zones));
for i = map_zones.keys
    vec_aux = str2num(char(i));
    index_aux = find(vec_aux);
    str_aux = sprintf('x_{%.0f} ',index_aux);
    if max(index_aux) <= nA;
    label = sprintf('PN: %s',str_aux);
    else
    label = sprintf('CN: %s',str_aux);
    end
    tickLabels(map_zones(char(i))) = cellstr(label);
end
colorbar('Ticks',1:length(map_zones),...
         'TickLabels',tickLabels,...
         'FontSize',9)
view(2)
xlabel('S_{in} [g/l]');
ylabel('Dilution Rate [day^{-1}]')
title(sprintf('\\fontsize{20} Operating Diagram'))
set(gca,'fontsize',18)
% fig.PaperPosition = [.25 .25 16 12];
fig.PaperPositionMode = 'auto';
print(sprintf('Images/OD_%s',file_name_in),'-dpng','-r0')
