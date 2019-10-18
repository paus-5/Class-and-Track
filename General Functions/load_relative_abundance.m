function [t_especes,especes] = load_relative_abundance(Reactor)
data = xlsread(sprintf('..\\Data\\especes%s',Reactor));
t_especes = data(:,1)-data(1,1);
especes = data(:,2:end);
end