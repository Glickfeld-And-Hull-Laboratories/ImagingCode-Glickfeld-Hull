all_mice = [];
all_depth = [];

%figure1&5-6
clear expt
ds = ['CrossOriRandDirTwoPhase_ExptList'];
eval(ds);
ind_area = find(strcmp([expt.img_loc],'V1'));
ind_sf = find([expt.SF] == 0.05);
mouse_list = [{expt.mouse}];
ind_use = intersect(ind_area,ind_sf);
nexpA = length(ind_use);
miceA = unique(mouse_list(ind_use));
nmiceA = length(miceA);
all_mice = [all_mice miceA];
depth = [expt.z];
all_depth = [all_depth depth(ind_use)];

%figure3
clear expt
ds = ['CrossOriRandPhase_ExptList'];
eval(ds);
expt1 = expt;
clear expt
ds = ['CrossOriRandPhase_lowSF_ExptList'];
eval(ds);
expt2 = expt;
clear expt
expt = combineStructures(expt1,expt2);
mouse_list = [{expt.mouse}];
nexpB = length(mouse_list);
miceB = unique(mouse_list);
nmiceB = length(miceB);
all_mice = [all_mice miceB];
depth = [expt.z];
all_depth = [all_depth depth];

%figure4
clear expt
ds = ['CrossOriRandPhase_ExptList'];
eval(ds);
mouse_list = [{expt.mouse}];
nexpC = length(mouse_list);
miceC = unique(mouse_list);
nmiceC = length(miceC);
all_mice = [all_mice miceC];
depth = [expt.z];
all_depth = [all_depth depth];

%figure7
ds = ['CrossOriSingleStimAdapt_ExptList'];
eval(ds);
mouse_list = [{expt.mouse}];
nexpD = length(mouse_list);
miceD = unique(mouse_list);
nmiceD = length(miceD);
all_mice = [all_mice miceD];
depth = [expt.z];
all_depth = [all_depth depth];

all_mice = unique(all_mice);
nmice_all = length(all_mice);
depths = unique(all_depth);
avg_depth = mean(depths,2);
sem_depth = std(depths,[],2)./sqrt(length(depths));
