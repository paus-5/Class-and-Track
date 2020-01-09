close all
clear
load('MAT_files\Classification_day183_gurobi')
class_NOB(28) = 0;
class_AOB(28) = 1;
yields_NOB(28) = 0;
yields_AOB(28) = 4.9539;
save('MAT_files\Classification_day183_gurobi_modified')
