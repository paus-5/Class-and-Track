clear
close all
file_in = 'default_case';
load(sprintf('MAT_files/%s',file_in));
labels = cell(n);
for i=1:n
    labels{i} = sprintf('x_{%s}',num2str(i));
end
imagesc(A);
colormap(jet)
text_strings = num2str(A(:), '%1.1f');       % Create strings from the matrix values
text_strings = strtrim(cellstr(text_strings));  % Remove any space padding
[x, y] = meshgrid(1:n);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), text_strings(:), ...  % Plot the strings
                'HorizontalAlignment', 'center','fontsize',12);
mid_value = mean(get(gca, 'CLim'));  % Get the middle value of the color range
text_colors = repmat(A(:) > mid_value, 1, 3);  % Choose white or black for the
                                               %   text color of the strings so
                                               %   they can be easily seen over
                                               %   the background color
set(hStrings, {'Color'}, num2cell(text_colors, 2));  % Change the text colors
colorbar
[cmin, cmax] = caxis;
bound = max(-cmin,cmax);
caxis([-bound bound])
set(gca, 'XTick', 1:n, 'XTickLabel', labels, 'YTick', 1:n, 'YTickLabel', labels)
title(sprintf('\\fontsize{20} Interaction Matrix'))
    set(gca,'fontsize',16)
h.PaperPositionMode = 'manual';
h.PaperPosition = [.25 .25 90 60];
print(sprintf('Images/Interactions_%s',file_in),'-dpng','-r0')
