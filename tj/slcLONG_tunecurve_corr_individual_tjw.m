%%
%clear everything
clear all
clear all global
clc
close all

%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\slc_long\multi_day_drift'; %folder to save files to
dataset = 'exp_list_slc_long_tjw'; %experiment list to pick files from
eval(dataset); %load dataset


%%

%we start with loading the proper information based on the mouse's info in
%the experiment list

d1 = 2; 
d2 = d1+2; 
d3 = d1+4;
d4 = d1+6;
d5 = d1+8;
mouse = expt(d1).mouse; 
ref_str_d1 = ['runs-',expt(d1).runs];
ref_str_d2 = ['runs-',expt(d2).runs]; 
ref_str_d3 = ['runs-',expt(d3).runs]; 
ref_str_d4 = ['runs-',expt(d4).runs];
ref_str_d5 = ['runs-',expt(d5).runs];
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date_d1 = expt(d1).date; 
date_d2 = expt(d2).date; 
date_d3 = expt(d3).date; 
date_d4 = expt(d4).date; 
date_d5 = expt(d5).date; 
img_folder_d1 = expt(d1).runs; 
img_folder_d2 = expt(d2).runs; 
img_folder_d3 = expt(d3).runs; 
img_folder_d4 = expt(d4).runs; 
img_folder_d5 = expt(d5).runs; 

%load ori info for each day
d1_ori = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningInfo.mat']));
d2_ori = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningInfo.mat']));
d3_ori = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'oriTuningInfo.mat']));
d4_ori = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'oriTuningInfo.mat']));
d5_ori = load(fullfile(fnout, [date_d5 '_' mouse], [date_d5 '_' mouse '_' ref_str_d5], [date_d5 '_' mouse '_' ref_str_d5 '_' 'oriTuningInfo.mat']));

d1_fits = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningAndFits.mat']));
d2_fits = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningAndFits.mat']));
d3_fits = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'oriTuningAndFits.mat']));
d4_fits = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'oriTuningAndFits.mat']));
d5_fits = load(fullfile(fnout, [date_d5 '_' mouse], [date_d5 '_' mouse '_' ref_str_d5], [date_d5 '_' mouse '_' ref_str_d5 '_' 'oriTuningAndFits.mat']));

%load multiday data for day 2
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));
d3_matches = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'multiday_alignment.mat']));
d4_matches = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'multiday_alignment.mat']));
d5_matches = load(fullfile(fnout, [date_d5 '_' mouse], [date_d5 '_' mouse '_' ref_str_d5], [date_d5 '_' mouse '_' ref_str_d5 '_' 'multiday_alignment.mat']));


%% matching indices

tuned_d1 = d1_ori.ind_theta90;
tuned_d2 = d2_ori.ind_theta90;
tuned_d3 = d3_ori.ind_theta90;
tuned_d4 = d4_ori.ind_theta90;
tuned_d5 = d5_ori.ind_theta90;

tuned_all = unique([tuned_d1 tuned_d2 tuned_d3 tuned_d4 tuned_d5]);

match_d2 = find([d2_matches.cellImageAlign.pass]); 
match_d3 = find([d3_matches.cellImageAlign.pass]); 
match_d4 = find([d4_matches.cellImageAlign.pass]); 
match_d5 = find([d5_matches.cellImageAlign.pass]); 

match_all = intersect(intersect(intersect(intersect(match_d2,match_d3),match_d4),match_d5),tuned_all);


%import avg response to each ori for each cell on each day
a = d1_fits.avgResponseEaOri(match_all,:);
b = d2_fits.avgResponseEaOri(match_all,:);
c = d3_fits.avgResponseEaOri(match_all,:);
d = d4_fits.avgResponseEaOri(match_all,:);
e = d5_fits.avgResponseEaOri(match_all,:);

%%
a_cells_corr = corr(a.');
b_cells_corr = corr(b.');
c_cells_corr = corr(c.');
d_cells_corr = corr(d.');
e_cells_corr = corr(e.');

%%
a_idx = eye(size(a_cells_corr));
a_Y = (1-a_idx).*a_cells_corr;
a_Z = a_cells_corr(~a_idx);

b_idx = eye(size(b_cells_corr));
b_Y = (1-b_idx).*b_cells_corr;
b_Z = b_cells_corr(~b_idx);

c_idx = eye(size(c_cells_corr));
c_Y = (1-c_idx).*c_cells_corr;
c_Z = c_cells_corr(~c_idx);

d_idx = eye(size(d_cells_corr));
d_Y = (1-d_idx).*d_cells_corr;
d_Z = d_cells_corr(~d_idx);

e_idx = eye(size(e_cells_corr));
e_Y = (1-e_idx).*e_cells_corr;
e_Z = e_cells_corr(~e_idx);

%%
abcde_mat = [a_Z b_Z c_Z d_Z e_Z];
abcde_corrs = corr(abcde_mat);

%%
figure;
scatter(1,[abcde_corrs(1,2), abcde_corrs(2,3), abcde_corrs(3,4), abcde_corrs(4,5)])
hold on
scatter(2, [abcde_corrs(1,3), abcde_corrs(2,4), abcde_corrs(3,5)])
scatter(3, [abcde_corrs(1,4), abcde_corrs(2,5)])
scatter(4, abcde_corrs(1,5))
xlim([0 5])
ylim([0 1])
xlabel('Days Apart')
ylabel('Correlation')
xticklabels({'','','7','','14','','21','','28'})


%%
%save vars
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_tunecurve_corr.mat']), 'abcde_corrs')


%%

