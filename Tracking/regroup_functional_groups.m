close all
clear
load('MAT_files\Classification_day_183_315_gurobi')
index_AOB = find(class_AOB);
index_NOB = find(class_NOB);
total_AOB = sum(OTU_interp(:,index_AOB),2);
total_NOB = sum(OTU_interp(:,index_NOB),2);
OTU_interp = [total_AOB total_NOB];
yields_AOB = [1/yA_ref ; 0];
yields_NOB = [0 ; 1/yB_ref];
class_AOB = [1; 0];
class_NOB = [0; 1];
save('MAT_files\Classification_day183_315_regroup_OTU')