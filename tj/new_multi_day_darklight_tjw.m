%This code is used to analyze individual subject data from the dark to
%light housing 2p imaging experiment, so that these data are prepared for
%pooled data afterwards.
%%
%clear everything
clear all
clear all global
clc
close all

%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\new_multi_day'; %folder to save files to
dataset = 'exp_list_darklight_actual_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%
%below is the legend for which numbers correspond to which imaging
%timepoints in this study
%img_point: 1 = baseline 1, 2 = baseline 2, 3 = post-dark, 4 = post-dark + 7d

%below are the mouse ID numbers involved in this experiment
%22-25 is 2554, 26-29 is 2555, 30-33 is 2556
%%

%we start with loading the proper information based on the mouse's info in
%the experiment list

%baseline1 to baseline2 comparison
d1 = 22+4; %base1 in expt list
d2 = 23+4; %base2 in expt list
mouse = expt(d1).mouse; 
ref_str_d1 = ['runs-',expt(d1).runs];
ref_str_d2 = ['runs-',expt(d2).runs]; 
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date_d1 = expt(d1).date; 
date_d2 = expt(d2).date; 
img_folder_d1 = expt(d1).runs; 
img_folder_d2 = expt(d2).runs; 


%load ori info for each day
d1_ori = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningInfo.mat']));
d2_ori = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningInfo.mat']));

%load multiday data for day 2
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));

%this function will take in the info loaded above and return the following:
%pref ori of tuned neurons matched across the 2 timepoints, the change in
%pref ori (d score) between days of those cells, the same results for
%broadness of tuning (k) and max dF/F, the cell  IDs for neurons that are
%tuned and matched
[prefori_1_2_match_tune, prefori_2_1_match_tune, d_score_prefori_1_2_match, k_1_2_match, k_2_1_match, dscore_k_1_2_match, max_1_2_match, max_2_1_match, dscore_max_1_2_match, id_match_tune_1_2] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);

%the sections below will do essentially the same thing as herre but for
%other comparison days
%%
%baseline2 to post-dark
d1 = 23+4; %base2 in expt list
d2 = 24+4; %postdark in expt list
mouse = expt(d1).mouse;
ref_str_d1 = ['runs-',expt(d1).runs];
ref_str_d2 = ['runs-',expt(d2).runs]; 
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date_d1 = expt(d1).date; 
date_d2 = expt(d2).date; 
img_folder_d1 = expt(d1).runs; 
img_folder_d2 = expt(d2).runs; 

%load ori info for each day
d1_ori = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningInfo.mat']));
d2_ori = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningInfo.mat']));

%load multiday data for day 2
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));


[prefori_2_3_match_tune, prefori_3_2_match_tune, d_score_prefori_2_3_match, k_2_3_match, k_3_2_match, dscore_k_2_3_match, max_2_3_match, max_3_2_match, dscore_max_2_3_match, id_match_tune_2_3] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);

%%
%post-dark to post-dark 7d
d1 = 24+4; %postdark in expt list
d2 = 25+4; %postdark7d in expt list
mouse = expt(d1).mouse; 
ref_str_d1 = ['runs-',expt(d1).runs];
ref_str_d2 = ['runs-',expt(d2).runs]; 
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date_d1 = expt(d1).date; 
date_d2 = expt(d2).date; 
img_folder_d1 = expt(d1).runs; 
img_folder_d2 = expt(d2).runs; 

%load ori info for each day
d1_ori = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningInfo.mat']));
d2_ori = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningInfo.mat']));

%load multiday data for day 2
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));


[prefori_3_4_match_tune, prefori_4_3_match_tune, d_score_prefori_3_4_match, k_3_4_match, k_4_3_match, dscore_k_3_4_match, max_3_4_match, max_4_3_match, dscore_max_3_4_match, id_match_tune_3_4] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);


%%
%save all of the change scores from above for future use
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_d_scores.mat']), 'd_score_prefori_1_2_match', 'd_score_prefori_2_3_match', 'd_score_prefori_3_4_match', ...
        'dscore_k_1_2_match', 'dscore_k_2_3_match', 'dscore_k_3_4_match', ...
        'dscore_max_1_2_match', 'dscore_max_2_3_match', 'dscore_max_3_4_match')

%%
%save pref ori info
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_ori_changes.mat']), 'prefori_1_2_match_tune', 'prefori_2_3_match_tune', 'prefori_2_1_match_tune', 'prefori_3_2_match_tune', ...
    'prefori_3_4_match_tune', 'prefori_4_3_match_tune')

%%
%save tuned and matched cell IDs
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_id_matches.mat']), 'id_match_tune_1_2', 'id_match_tune_2_3', 'id_match_tune_3_4')

%%
%make a quick plot of pref or changes across sessions
figure;
cdfplot(d_score_prefori_1_2_match);
hold on;
cdfplot(d_score_prefori_2_3_match);
cdfplot(d_score_prefori_3_4_match);
hold off
legend('1-2', '2-3', '3-4')



%%
%now we want to find cells that were matched across 2 sets of comparisons
%that we identified as being key to compare:
%1) 1-2 VS 2-3; base1 to base2 VS base2 to postdark - baseline period vs
%dark housing manipulation
%2) 2-3 VS 3-4; base2 to postdark VS postdark to postdark+7d - dark housing vs
%1 week follow-up


%%
%base1 to base2 VS base2 to postdark

%find cell IDs for cells that are present in both comparison sets; i.e.
%these are cells that are tuned and matched for each respective period AND
%that can be identified between periods
base1_base2_vs_base2_postdark_id = intersect(id_match_tune_1_2, id_match_tune_2_3);

%we will load info here similar to what we did above

d1 = 22; %base1 in expt list
d2 = 23; %base2 in expt list
mouse = expt(d1).mouse; 
ref_str_d1 = ['runs-',expt(d1).runs];
ref_str_d2 = ['runs-',expt(d2).runs]; 
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date_d1 = expt(d1).date; 
date_d2 = expt(d2).date; 
img_folder_d1 = expt(d1).runs; 
img_folder_d2 = expt(d2).runs; 

%load ori info for each day
d1_ori = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningInfo.mat']));
d2_ori = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningInfo.mat']));

%load multiday data for day 2
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));

%do the same as above for the other comparison

d3 = 23; %postdark in expt list
d4 = 24; %postdarkd in expt list
mouse = expt(d3).mouse; 
ref_str_d3 = ['runs-',expt(d3).runs];
ref_str_d4 = ['runs-',expt(d4).runs]; 
img_area = expt(d3).img_loc{1};
img_layer = expt(d3).img_loc{2};
date_d3 = expt(d3).date; 
date_d4 = expt(d4).date; 
img_folder_d3 = expt(d3).runs; 
img_folder_d4 = expt(d4).runs; 

%load ori info for each day
d3_ori = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'oriTuningInfo.mat']));
d4_ori = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'oriTuningInfo.mat']));

%load multiday data for days 2
d4_matches = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'multiday_alignment.mat']));

%load k and max vals
d3_k_max = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'k_and_max_vals.mat']));
d4_k_max = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'k_and_max_vals.mat']));

%find the pref oris for cells on the first comparison that can be found on
%the second comparison as well

comparison_prefori_1_2_match_tune = d1_ori.prefOri(1,base1_base2_vs_base2_postdark_id);
comparison_prefori_2_1_match_tune = d2_ori.prefOri(1,base1_base2_vs_base2_postdark_id);

comparison_prefori_2_3_match_tune = d3_ori.prefOri(1,base1_base2_vs_base2_postdark_id);
comparison_prefori_3_2_match_tune = d4_ori.prefOri(1,base1_base2_vs_base2_postdark_id);

%same as pref ori above but for k here and max below
comparison_k_1_2_match_tune = d1_k_max.k1(1,base1_base2_vs_base2_postdark_id);
comparison_k_2_1_match_tune = d2_k_max.k1(1,base1_base2_vs_base2_postdark_id);

comparison_k_2_3_match_tune = d3_k_max.k1(1,base1_base2_vs_base2_postdark_id);
comparison_k_3_2_match_tune = d4_k_max.k1(1,base1_base2_vs_base2_postdark_id);

comparison_max_1_2_match_tune = d1_k_max.max_dfof(1,base1_base2_vs_base2_postdark_id);
comparison_max_2_1_match_tune = d2_k_max.max_dfof(1,base1_base2_vs_base2_postdark_id);

comparison_max_2_3_match_tune = d3_k_max.max_dfof(1,base1_base2_vs_base2_postdark_id);
comparison_max_3_2_match_tune = d4_k_max.max_dfof(1,base1_base2_vs_base2_postdark_id);

%%
%now we take the info we got above and convert to change in pref ori
%scores, change in k scores, and change in max scores

%d1
comparison_dscore_prefori_d1_d2_match = double(comparison_prefori_1_2_match_tune>90);
comparison_dscore_prefori_d1_d2_match(comparison_dscore_prefori_d1_d2_match>0) = 180;
comparison_dscore_prefori_d1_d2_match = abs(comparison_dscore_prefori_d1_d2_match-comparison_prefori_1_2_match_tune);
%d2
comparison_dscore_prefori_d2_match = double(comparison_prefori_2_1_match_tune>90);
comparison_dscore_prefori_d2_match(comparison_dscore_prefori_d2_match>0) = 180;
comparison_dscore_prefori_d2_match = abs(comparison_dscore_prefori_d2_match-comparison_prefori_2_1_match_tune);

comparison_d_score_prefori_1_2 = abs(comparison_dscore_prefori_d1_d2_match-comparison_dscore_prefori_d2_match);


%d3
comparison_dscore_prefori_d3_d4_match = double(comparison_prefori_2_3_match_tune>90);
comparison_dscore_prefori_d3_d4_match(comparison_dscore_prefori_d3_d4_match>0) = 180;
comparison_dscore_prefori_d3_d4_match = abs(comparison_dscore_prefori_d3_d4_match-comparison_prefori_2_3_match_tune);
%d4
comparison_dscore_prefori_d4_match = double(comparison_prefori_3_2_match_tune>90);
comparison_dscore_prefori_d4_match(comparison_dscore_prefori_d4_match>0) = 180;
comparison_dscore_prefori_d4_match = abs(comparison_dscore_prefori_d4_match-comparison_prefori_3_2_match_tune);

comparison_d_score_prefori_2_3 = abs(comparison_dscore_prefori_d3_d4_match-comparison_dscore_prefori_d4_match);



comparison_d_score_k_1_2 = (comparison_k_1_2_match_tune-comparison_k_2_1_match_tune);


comparison_d_score_k_2_3 = (comparison_k_2_3_match_tune-comparison_k_3_2_match_tune);




comparison_d_score_max_1_2 = (comparison_max_1_2_match_tune-comparison_max_2_1_match_tune);


comparison_d_score_max_2_3 = (comparison_max_2_3_match_tune-comparison_max_3_2_match_tune);






%below we will do the same as above but for the other sets of comparisons
%%
%base2 to postdark (2-3) vs postdark to postdark+7d (3-4)
base2_postdark_vs_postdark_postdark7d_id = intersect(id_match_tune_2_3, id_match_tune_3_4);

d1 = 23; %base1 in expt list
d2 = 24; %base1+4d in expt list
mouse = expt(d1).mouse; 
ref_str_d1 = ['runs-',expt(d1).runs];
ref_str_d2 = ['runs-',expt(d2).runs]; 
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date_d1 = expt(d1).date; 
date_d2 = expt(d2).date; 
img_folder_d1 = expt(d1).runs; 
img_folder_d2 = expt(d2).runs; 

%load ori info for each day
d1_ori = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningInfo.mat']));
d2_ori = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningInfo.mat']));

%load multiday data for day 2
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));

d3 = 24; %base2+4hrs in expt list
d4 = 25; %postdark in expt list
mouse = expt(d3).mouse; 
ref_str_d3 = ['runs-',expt(d3).runs];
ref_str_d4 = ['runs-',expt(d4).runs]; 
img_area = expt(d3).img_loc{1};
img_layer = expt(d3).img_loc{2};
date_d3 = expt(d3).date; 
date_d4 = expt(d4).date; 
img_folder_d3 = expt(d3).runs; 
img_folder_d4 = expt(d4).runs; 

%load ori info for each day
d3_ori = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'oriTuningInfo.mat']));
d4_ori = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'oriTuningInfo.mat']));

%load multiday data for day 2
d4_matches = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'multiday_alignment.mat']));

%load k and max vals
d3_k_max = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'k_and_max_vals.mat']));
d4_k_max = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'k_and_max_vals.mat']));



second_comparison_prefori_2_3_match_tune = d1_ori.prefOri(1,base2_postdark_vs_postdark_postdark7d_id);
second_comparison_prefori_3_2_match_tune = d2_ori.prefOri(1,base2_postdark_vs_postdark_postdark7d_id);

comparison_prefori_3_4_match_tune = d3_ori.prefOri(1,base2_postdark_vs_postdark_postdark7d_id);
comparison_prefori_4_3_match_tune = d4_ori.prefOri(1,base2_postdark_vs_postdark_postdark7d_id);

second_comparison_k_2_3_match_tune = d1_k_max.k1(1,base2_postdark_vs_postdark_postdark7d_id);
second_comparison_k_3_2_match_tune = d2_k_max.k1(1,base2_postdark_vs_postdark_postdark7d_id);

comparison_k_3_4_match_tune = d3_k_max.k1(1,base2_postdark_vs_postdark_postdark7d_id);
comparison_k_4_3_match_tune = d4_k_max.k1(1,base2_postdark_vs_postdark_postdark7d_id);

second_comparison_max_2_3_match_tune = d1_k_max.max_dfof(1,base2_postdark_vs_postdark_postdark7d_id);
second_comparison_max_3_2_match_tune = d2_k_max.max_dfof(1,base2_postdark_vs_postdark_postdark7d_id);

comparison_max_3_4_match_tune = d3_k_max.max_dfof(1,base2_postdark_vs_postdark_postdark7d_id);
comparison_max_4_3_match_tune = d4_k_max.max_dfof(1,base2_postdark_vs_postdark_postdark7d_id);




%%

%d1
comparison_dscore_prefori_d1_d2_match = double(second_comparison_prefori_2_3_match_tune>90);
comparison_dscore_prefori_d1_d2_match(comparison_dscore_prefori_d1_d2_match>0) = 180;
comparison_dscore_prefori_d1_d2_match = abs(comparison_dscore_prefori_d1_d2_match-second_comparison_prefori_2_3_match_tune);
%d2
comparison_dscore_prefori_d2_match = double(second_comparison_prefori_3_2_match_tune>90);
comparison_dscore_prefori_d2_match(comparison_dscore_prefori_d2_match>0) = 180;
comparison_dscore_prefori_d2_match = abs(comparison_dscore_prefori_d2_match-second_comparison_prefori_3_2_match_tune);

second_comparison_d_score_prefori_2_3 = abs(comparison_dscore_prefori_d1_d2_match-comparison_dscore_prefori_d2_match);


%d3
comparison_dscore_prefori_d3_d4_match = double(comparison_prefori_3_4_match_tune>90);
comparison_dscore_prefori_d3_d4_match(comparison_dscore_prefori_d3_d4_match>0) = 180;
comparison_dscore_prefori_d3_d4_match = abs(comparison_dscore_prefori_d3_d4_match-comparison_prefori_3_4_match_tune);
%d4
comparison_dscore_prefori_d4_match = double(comparison_prefori_4_3_match_tune>90);
comparison_dscore_prefori_d4_match(comparison_dscore_prefori_d4_match>0) = 180;
comparison_dscore_prefori_d4_match = abs(comparison_dscore_prefori_d4_match-comparison_prefori_4_3_match_tune);

comparison_d_score_prefori_3_4 = abs(comparison_dscore_prefori_d3_d4_match-comparison_dscore_prefori_d4_match);



second_comparison_d_score_k_2_3 = (second_comparison_k_2_3_match_tune-second_comparison_k_3_2_match_tune);


comparison_d_score_k_3_4 = (comparison_k_3_4_match_tune-comparison_k_4_3_match_tune);


second_comparison_d_score_max_2_3 = (second_comparison_max_2_3_match_tune-second_comparison_max_3_2_match_tune);


comparison_d_score_max_3_4 = (comparison_max_3_4_match_tune-comparison_max_4_3_match_tune);




%%
%now we can make some individual plots for the comparisons we want to make
%for all 3 vars (pref ori, k, max)
fig = figure;
sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(comparison_d_score_prefori_1_2);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(comparison_d_score_prefori_2_3);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',1);
legend(['Base1 - Base2 (n= ', num2str(length(comparison_d_score_prefori_1_2)), ')'], ['Base2 - Post-dark (n= ', num2str(length(comparison_d_score_prefori_2_3)), ')'], 'Location', 'southeast')
xlim([0 90])
xlabel([])
ylabel([])
title([])
hold off

fig = figure;
sgtitle('Pref Ori Changes')
% b=subplot(2,2,2)
h=cdfplot(second_comparison_d_score_prefori_2_3);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(comparison_d_score_prefori_3_4);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',1);
legend(['Base2 - Post-dark (n= ', num2str(length(second_comparison_d_score_prefori_2_3)), ')'], ['Post-dark - Post-dark+7d (n= ', num2str(length(comparison_d_score_prefori_3_4)), ')'], 'Location', 'southeast')
xlabel([])
ylabel([])
title([])
hold off

print(fullfile(newfnout, [mouse,  '_comparison_prefori_change.pdf']), '-dpdf', '-bestfit')

%%
fig = figure;
sgtitle('k Changes')
% a=subplot(2,2,1)
h=cdfplot(comparison_d_score_k_1_2);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(comparison_d_score_k_2_3);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',1);
legend(['Base1 - Base2 (n= ', num2str(length(comparison_d_score_k_1_2)), ')'], ['Base2 - Post-dark (n= ', num2str(length(comparison_d_score_k_2_3)), ')'], 'Location', 'southeast')
xlim([-30 30])
xlabel([])
ylabel([])
title([])
hold off

fig = figure;
sgtitle('k Changes')
% b=subplot(2,2,2)
h=cdfplot(second_comparison_d_score_k_2_3);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(comparison_d_score_k_3_4);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',1);
legend(['Base2 - Post-dark (n= ', num2str(length(second_comparison_d_score_k_2_3)), ')'], ['Post-dark - Post-dark+7d (n= ', num2str(length(comparison_d_score_k_3_4)), ')'], 'Location', 'southeast')
xlabel([])
ylabel([])
title([])
hold off

print(fullfile(newfnout, [mouse,  '_comparison_k_change.pdf']), '-dpdf', '-bestfit')

%%

fig = figure;
sgtitle('Max dF/F Changes')
% a=subplot(2,2,1)
h=cdfplot(comparison_d_score_max_1_2);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(comparison_d_score_max_2_3);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',1);
legend(['Base1 - Base2 (n= ', num2str(length(comparison_d_score_max_1_2)), ')'], ['Base2 - Post-dark (n= ', num2str(length(comparison_d_score_max_2_3)), ')'], 'Location', 'southeast')
xlim([-1 1])
xlabel([])
ylabel([])
title([])
hold off

fig = figure;
sgtitle('Max dF/F Changes')
% b=subplot(2,2,2)
h=cdfplot(second_comparison_d_score_max_2_3);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(comparison_d_score_max_3_4);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',1);
legend(['Base2 - Post-dark (n= ', num2str(length(second_comparison_d_score_max_2_3)), ')'], ['Post-dark - Post-dark+7d (n= ', num2str(length(comparison_d_score_max_3_4)), ')'], 'Location', 'southeast')
xlabel([])
ylabel([])
title([])
hold off

print(fullfile(newfnout, [mouse,  '_comparison_max_change.pdf']), '-dpdf', '-bestfit')





%%
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_same_cell_dscores.mat']), 'comparison_d_score_prefori_1_2', 'comparison_d_score_prefori_2_3', 'second_comparison_d_score_prefori_2_3','comparison_d_score_prefori_3_4', ...
    'comparison_d_score_k_1_2', 'comparison_d_score_k_2_3', 'second_comparison_d_score_k_2_3' ,'comparison_d_score_k_3_4', ...
    'comparison_d_score_max_1_2', 'comparison_d_score_max_2_3', 'second_comparison_d_score_max_2_3', 'comparison_d_score_max_3_4')

%%
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_same_cell_k_vals.mat']), 'comparison_k_1_2_match_tune', 'comparison_k_2_1_match_tune', 'comparison_k_2_3_match_tune', 'comparison_k_3_2_match_tune', ...
    'second_comparison_k_2_3_match_tune', 'second_comparison_k_3_2_match_tune', 'comparison_k_3_4_match_tune', 'comparison_k_4_3_match_tune')

save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_same_cell_max_vals.mat']), 'comparison_max_1_2_match_tune', 'comparison_max_2_1_match_tune', 'comparison_max_2_3_match_tune', 'comparison_max_3_2_match_tune', ...
    'second_comparison_max_2_3_match_tune', 'second_comparison_max_3_2_match_tune', 'comparison_max_3_4_match_tune', 'comparison_max_4_3_match_tune')

save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_same_cell_prefori_vals.mat']), 'comparison_prefori_1_2_match_tune', 'comparison_prefori_2_1_match_tune', 'comparison_prefori_2_3_match_tune', 'comparison_prefori_3_2_match_tune', ...
    'second_comparison_prefori_2_3_match_tune', 'second_comparison_prefori_3_2_match_tune', 'comparison_prefori_3_4_match_tune', 'comparison_prefori_4_3_match_tune')


