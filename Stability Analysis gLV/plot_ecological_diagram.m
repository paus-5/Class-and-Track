clear
close all
file_name_in = 'parameters_modified';
load(sprintf('MAT_files\\Operating_Diagram_%s',file_name_in))
%Only for parameters_modifed
number_of_zones = length(map_zones_ED)-1;
zones_ED(zones_ED == 4) = 2;
%end fix ED
%Normal Zone
% number_of_zones = length(map_zones_ED);
cmap1 = hsv(number_of_zones);
% End_zone
figure;
[X,Y] = meshgrid(s_in_vector,D_vector);
surfc(X,Y,zones_ED','EdgeColor','none','LineStyle','none','FaceLighting','phong');%view(2)
colormap(cmap1)
tickLabels = cell(1,number_of_zones);
for i = map_zones_ED.keys
    vec_aux = str2num(char(i));
    length_vec = length(vec_aux);
    index_aux = find(vec_aux);
    if length_vec <= n
        str_aux = sprintf('x_{%.0f}',index_aux);
        if max(index_aux) == 0;
            label = 'Washout';
        elseif max(index_aux) <= nA;
            label = sprintf('%s',str_aux);
        else
            label = sprintf('%s',str_aux);
        end
    else
        number_of_eq = length_vec/n;
        mat_aux = reshape(vec_aux,n,number_of_eq);
        str_aux = '';
        for j = 1:number_of_eq
            index_aux = find(mat_aux(:,j));
            str_aux = sprintf('%s%0.f) %s\\newline',str_aux,j,sprintf('x_{%.0f}',index_aux));
        end
        str_aux = str_aux(1:end);
%         text(
%         flag_PN = sum(any(mat_aux(nA+1:end,:))) == 0;
%         flag_PN_CN = sum(any(mat_aux(nA+1:end,:))) == 1;
%         flag_CN = all(any(mat_aux(nA+1:end,:)));
%         flag_washout_PN = ~all(any(mat_aux));
        label = sprintf('%s ',str_aux);
    end
    tickLabels(map_zones_ED(char(i))) = cellstr(label);
end
colorbar('Ticks',1:number_of_zones,...
    'TickLabels',tickLabels,...
    'FontSize',18)
view(2)
xlabel('S_{in} [g/l]');
ylabel('Dilution Rate [day^{-1}]')
title(sprintf('\\fontsize{20} Ecological Diagram'))
set(gca,'fontsize',18)
% fig.PaperPosition = [.25 .25 16 12];
fig.PaperPositionMode = 'auto';
print(sprintf('Images/ED_%s',file_name_in),'-dpng','-r0')
