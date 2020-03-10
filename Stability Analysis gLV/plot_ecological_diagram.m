clear
close all
file_name_in = 'default_case';
load(sprintf('MAT_files\\Operating_Diagram_%s',file_name_in))
figure;
number_of_zones = length(map_zones_ED);
cmap1 = hsv(number_of_zones);
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
            label = sprintf('PN: %s',str_aux);
        else
            label = sprintf('CN: %s',str_aux);
        end
    else 
        mat_aux = reshape(vec_aux,n,length_vec/n);
        flag_PN_and_CN = sum(any(mat_aux(nA+1:end,:)))<length_vec/n;
        flag_CN = all(any(mat_aux(nA+1:end,:)));
        flag_washout_PN = ~all(any(mat_aux));
        if flag_CN
            str_aux = 'CN';
        elseif flag_PN_and_CN
            str_aux = 'CN & PN';
        elseif flag_washout_PN
            str_aux = 'Washout & PN';
        end
        label = sprintf('%s: MSE',str_aux);
    end
    tickLabels(map_zones_ED(char(i))) = cellstr(label);
end
colorbar('Ticks',1:number_of_zones,...
    'TickLabels',tickLabels,...
    'FontSize',9)
view(2)
xlabel('S_{in} [g/l]');
ylabel('Dilution Rate [day^{-1}]')
title(sprintf('\\fontsize{20} Ecological Diagram'))
set(gca,'fontsize',18)
% fig.PaperPosition = [.25 .25 16 12];
fig.PaperPositionMode = 'auto';
print(sprintf('Images/ED_%s',file_name_in),'-dpng','-r0')
