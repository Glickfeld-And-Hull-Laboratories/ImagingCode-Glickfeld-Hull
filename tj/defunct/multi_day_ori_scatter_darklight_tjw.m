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
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\multi_day'; %folder to save files to
dataset = 'exp_list_darklight_actual_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%
%below is the legend for which numbers correspond to which imaging
%timepoints in this study
%img_point: 1 = baseline1, 2 = baseline1 + 4d, 3 = baseline2, 4 = baseline2
%+ 4hrs, 5 = post 4d dark, 6 = post dark + 4hrs, 7 = post dark + 7d

%below are the mouse ID numbers involved in this experiment
%1-7 is 2537, 8-14 is 2538, 15-21 is 2543
%%

%we start with loading the proper information based on the mouse's info in
%the experiment list

%baseline1 to baseline1+4d comparison
d1 = 1; %base1 in expt list
d2 = 2; %base1+4d in expt list
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
%baseline1 to baseline2
d1 = 1; %base1 in expt list
d2 = 3; %base2 in expt list
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


[prefori_1_3_match_tune, prefori_3_1_match_tune, d_score_prefori_1_3_match, k_1_3_match, k_3_1_match, dscore_k_1_3_match, max_1_3_match, max_3_1_match, dscore_max_1_3_match, id_match_tune_1_3] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);

%%
%baseline2 to baseline2+4hr
d1 = 3; %base2 in expt list
d2 = 4; %base2+4hr in expt list
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
%baseline2+4hr to postdark
d1 = 4; %base2+4hr in expt list
d2 = 5; %postdark in expt list
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


[prefori_4_5_match_tune, prefori_5_4_match_tune, d_score_prefori_4_5_match, k_4_5_match, k_5_4_match, dscore_k_4_5_match, max_4_5_match, max_5_4_match, dscore_max_4_5_match, id_match_tune_4_5] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);

%%
%postdark to postlight4hr
d1 = 5; %postdark in expt list
d2 = 6; %postlight4hr in expt list
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


[prefori_5_6_match_tune, prefori_6_5_match_tune, d_score_prefori_5_6_match, k_5_6_match, k_6_5_match, dscore_k_5_6_match, max_5_6_match, max_6_5_match, dscore_max_5_6_match, id_match_tune_5_6] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);

%%
%postdark to postlight7d
d1 = 5; %postdark in expt list
d2 = 7; %postlight7d in expt list
mouse = expt(d1).mouse; %mouse
ref_str_d1 = ['runs-',expt(d1).runs];
ref_str_d2 = ['runs-',expt(d2).runs]; 
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date_d1 = expt(d1).date; %day 1 (ref day) date
date_d2 = expt(d2).date; %day 2 date
img_folder_d1 = expt(d1).runs; %img folder of day 1
img_folder_d2 = expt(d2).runs; %img folder of day 2

%load ori info for each day
d1_ori = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningInfo.mat']));
d2_ori = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningInfo.mat']));

%load multiday data for day 2
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));


[prefori_5_7_match_tune, prefori_7_5_match_tune, d_score_prefori_5_7_match, k_5_7_match, k_7_5_match, dscore_k_5_7_match, max_5_7_match, max_7_5_match, dscore_max_5_7_match, id_match_tune_5_7] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);

%%
%save all of the change scores from above for future use
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_d_scores.mat']), 'd_score_prefori_1_2_match', 'd_score_prefori_1_3_match', 'd_score_prefori_3_4_match', 'd_score_prefori_4_5_match', 'd_score_prefori_5_6_match', 'd_score_prefori_5_7_match', ...
        'dscore_k_1_2_match', 'dscore_k_1_3_match', 'dscore_k_3_4_match', 'dscore_k_4_5_match', 'dscore_k_5_6_match', 'dscore_k_5_7_match', ...
        'dscore_max_1_2_match', 'dscore_max_1_3_match', 'dscore_max_3_4_match', 'dscore_max_4_5_match', 'dscore_max_5_6_match', 'dscore_max_5_7_match')

%%
%save pref ori info
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_ori_changes.mat']), 'prefori_1_2_match_tune', 'prefori_1_3_match_tune', 'prefori_2_1_match_tune', 'prefori_3_1_match_tune', ...
    'prefori_3_4_match_tune', 'prefori_4_3_match_tune', 'prefori_4_5_match_tune', 'prefori_5_4_match_tune', 'prefori_5_6_match_tune', 'prefori_5_7_match_tune', 'prefori_6_5_match_tune', ...
    'prefori_7_5_match_tune')

%%
%save tuned and matched cell IDs
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_id_matches.mat']), 'id_match_tune_1_2', 'id_match_tune_1_3', 'id_match_tune_3_4', 'id_match_tune_4_5', 'id_match_tune_5_6', 'id_match_tune_5_7')

%%
%make a quick plot of pref or changes across sessions
figure;
cdfplot(d_score_prefori_1_2_match);
hold on;
cdfplot(d_score_prefori_1_3_match);
cdfplot(d_score_prefori_3_4_match);
cdfplot(d_score_prefori_4_5_match);
cdfplot(d_score_prefori_5_6_match);
cdfplot(d_score_prefori_5_7_match);
hold off
legend('1-2', '1-3', '3-4', '4-5', '5-6', '5-7')


%%
%now we want to find cells that were matched across 2 sets of comparisons
%that we identified as being key to compare:
%1) 1-3 VS 5-7; base1 to base2 VS postdark to postdarkd - light recovery period vs
%baseline period of same length (7d)
%2) 1-2 VS 4-5; base1 to base1+4d VS base2+4hr to postdark - post dark housing vs
%baseline period of 4 d
%3) 3-4 VS 5-6; base2 to base2+4hr vs postdark to postdark+4hr - 4hrs light
%recoery after dark housing vs baseline period of 4hrs


%%
%base1 to base2 (1-3) vs postdark to postdarkd (5-7)

%find cell IDs for cells that are present in both comparison sets; i.e.
%these are cells that are tuned and matched for each respective period AND
%that can be identified between periods
base1_base2_vs_postdark_postdark7d_id = intersect(id_match_tune_1_3, id_match_tune_5_7);

%we will load info here similar to what we did above

d1 = 1; %base1 in expt list
d2 = 3; %base2 in expt list
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

d3 = 5; %postdark in expt list
d4 = 7; %postdarkd in expt list
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
comparison_prefori_1_3_match_tune = d1_ori.prefOri(1,base1_base2_vs_postdark_postdark7d_id);
comparison_prefori_3_1_match_tune = d2_ori.prefOri(1,base1_base2_vs_postdark_postdark7d_id);

comparison_prefori_5_7_match_tune = d3_ori.prefOri(1,base1_base2_vs_postdark_postdark7d_id);
comparison_prefori_7_5_match_tune = d4_ori.prefOri(1,base1_base2_vs_postdark_postdark7d_id);

%same as pref ori above but for k here and max below
comparison_k_1_3_match_tune = d1_k_max.k1(1,base1_base2_vs_postdark_postdark7d_id);
comparison_k_3_1_match_tune = d2_k_max.k1(1,base1_base2_vs_postdark_postdark7d_id);

comparison_k_5_7_match_tune = d3_k_max.k1(1,base1_base2_vs_postdark_postdark7d_id);
comparison_k_7_5_match_tune = d4_k_max.k1(1,base1_base2_vs_postdark_postdark7d_id);

comparison_max_1_3_match_tune = d1_k_max.max_dfof(1,base1_base2_vs_postdark_postdark7d_id);
comparison_max_3_1_match_tune = d2_k_max.max_dfof(1,base1_base2_vs_postdark_postdark7d_id);

comparison_max_5_7_match_tune = d3_k_max.max_dfof(1,base1_base2_vs_postdark_postdark7d_id);
comparison_max_7_5_match_tune = d4_k_max.max_dfof(1,base1_base2_vs_postdark_postdark7d_id);


%%
%now we take the info we got above and convert to change in pref ori
%scores, change in k scores, and change in max scores

%d1
comparison_dscore_prefori_d1_d2_match = double(comparison_prefori_1_3_match_tune>90);
comparison_dscore_prefori_d1_d2_match(comparison_dscore_prefori_d1_d2_match>0) = 180;
comparison_dscore_prefori_d1_d2_match = abs(comparison_dscore_prefori_d1_d2_match-comparison_prefori_1_3_match_tune);
%d2
comparison_dscore_prefori_d2_match = double(comparison_prefori_3_1_match_tune>90);
comparison_dscore_prefori_d2_match(comparison_dscore_prefori_d2_match>0) = 180;
comparison_dscore_prefori_d2_match = abs(comparison_dscore_prefori_d2_match-comparison_prefori_3_1_match_tune);

comparison_d_score_prefori_1_3 = abs(comparison_dscore_prefori_d1_d2_match-comparison_dscore_prefori_d2_match);


%d3
comparison_dscore_prefori_d3_d4_match = double(comparison_prefori_5_7_match_tune>90);
comparison_dscore_prefori_d3_d4_match(comparison_dscore_prefori_d3_d4_match>0) = 180;
comparison_dscore_prefori_d3_d4_match = abs(comparison_dscore_prefori_d3_d4_match-comparison_prefori_5_7_match_tune);
%d4
comparison_dscore_prefori_d4_match = double(comparison_prefori_7_5_match_tune>90);
comparison_dscore_prefori_d4_match(comparison_dscore_prefori_d4_match>0) = 180;
comparison_dscore_prefori_d4_match = abs(comparison_dscore_prefori_d4_match-comparison_prefori_7_5_match_tune);

comparison_d_score_prefori_5_7 = abs(comparison_dscore_prefori_d3_d4_match-comparison_dscore_prefori_d4_match);



comparison_d_score_k_1_3 = (comparison_k_1_3_match_tune-comparison_k_3_1_match_tune);


comparison_d_score_k_5_7 = (comparison_k_5_7_match_tune-comparison_k_7_5_match_tune);




comparison_d_score_max_1_3 = (comparison_max_1_3_match_tune-comparison_max_3_1_match_tune);


comparison_d_score_max_5_7 = (comparison_max_5_7_match_tune-comparison_max_7_5_match_tune);


%below we will do the same as above but for the other sets of comparisons
%%
%base1 to base2 (1-2) vs base2+4hrs to postdark (4-5)
base1_base1plus4d_vs_base2plus4h_postdark_id = intersect(id_match_tune_1_2, id_match_tune_4_5);

d1 = 1; %base1 in expt list
d2 = 2; %base1+4d in expt list
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

d3 = 4; %base2+4hrs in expt list
d4 = 5; %postdark in expt list
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


comparison_prefori_1_2_match_tune = d1_ori.prefOri(1,base1_base1plus4d_vs_base2plus4h_postdark_id);
comparison_prefori_2_1_match_tune = d2_ori.prefOri(1,base1_base1plus4d_vs_base2plus4h_postdark_id);

comparison_prefori_4_5_match_tune = d3_ori.prefOri(1,base1_base1plus4d_vs_base2plus4h_postdark_id);
comparison_prefori_5_4_match_tune = d4_ori.prefOri(1,base1_base1plus4d_vs_base2plus4h_postdark_id);

comparison_k_1_2_match_tune = d1_k_max.k1(1,base1_base1plus4d_vs_base2plus4h_postdark_id);
comparison_k_2_1_match_tune = d2_k_max.k1(1,base1_base1plus4d_vs_base2plus4h_postdark_id);

comparison_k_4_5_match_tune = d3_k_max.k1(1,base1_base1plus4d_vs_base2plus4h_postdark_id);
comparison_k_5_4_match_tune = d4_k_max.k1(1,base1_base1plus4d_vs_base2plus4h_postdark_id);

comparison_max_1_2_match_tune = d1_k_max.max_dfof(1,base1_base1plus4d_vs_base2plus4h_postdark_id);
comparison_max_2_1_match_tune = d2_k_max.max_dfof(1,base1_base1plus4d_vs_base2plus4h_postdark_id);

comparison_max_4_5_match_tune = d3_k_max.max_dfof(1,base1_base1plus4d_vs_base2plus4h_postdark_id);
comparison_max_5_4_match_tune = d4_k_max.max_dfof(1,base1_base1plus4d_vs_base2plus4h_postdark_id);



%%

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
comparison_dscore_prefori_d3_d4_match = double(comparison_prefori_4_5_match_tune>90);
comparison_dscore_prefori_d3_d4_match(comparison_dscore_prefori_d3_d4_match>0) = 180;
comparison_dscore_prefori_d3_d4_match = abs(comparison_dscore_prefori_d3_d4_match-comparison_prefori_4_5_match_tune);
%d4
comparison_dscore_prefori_d4_match = double(comparison_prefori_5_4_match_tune>90);
comparison_dscore_prefori_d4_match(comparison_dscore_prefori_d4_match>0) = 180;
comparison_dscore_prefori_d4_match = abs(comparison_dscore_prefori_d4_match-comparison_prefori_5_4_match_tune);

comparison_d_score_prefori_4_5 = abs(comparison_dscore_prefori_d3_d4_match-comparison_dscore_prefori_d4_match);



comparison_d_score_k_1_2 = (comparison_k_1_2_match_tune-comparison_k_2_1_match_tune);


comparison_d_score_k_4_5 = (comparison_k_4_5_match_tune-comparison_k_5_4_match_tune);


comparison_d_score_max_1_2 = (comparison_max_1_2_match_tune-comparison_max_2_1_match_tune);


comparison_d_score_max_4_5 = (comparison_max_4_5_match_tune-comparison_max_5_4_match_tune);



%%
%base2 to base2+4hr vs postdark to postdark+4hr
base2_base2plus4h_vs_postdark_postdarkplus4h_id = intersect(id_match_tune_3_4, id_match_tune_5_6);

d1 = 3; %base2 in expt list
d2 = 4; %base2+4hr in expt list
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

d3 = 5; %postdark in expt list
d4 = 6; %postdark+4hr in expt list
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


comparison_prefori_3_4_match_tune = d1_ori.prefOri(1,base2_base2plus4h_vs_postdark_postdarkplus4h_id);
comparison_prefori_4_3_match_tune = d2_ori.prefOri(1,base2_base2plus4h_vs_postdark_postdarkplus4h_id);

comparison_prefori_5_6_match_tune = d3_ori.prefOri(1,base2_base2plus4h_vs_postdark_postdarkplus4h_id);
comparison_prefori_6_5_match_tune = d4_ori.prefOri(1,base2_base2plus4h_vs_postdark_postdarkplus4h_id);

comparison_k_3_4_match_tune = d1_k_max.k1(1,base2_base2plus4h_vs_postdark_postdarkplus4h_id);
comparison_k_4_3_match_tune = d2_k_max.k1(1,base2_base2plus4h_vs_postdark_postdarkplus4h_id);

comparison_k_5_6_match_tune = d3_k_max.k1(1,base2_base2plus4h_vs_postdark_postdarkplus4h_id);
comparison_k_6_5_match_tune = d4_k_max.k1(1,base2_base2plus4h_vs_postdark_postdarkplus4h_id);

comparison_max_3_4_match_tune = d1_k_max.max_dfof(1,base2_base2plus4h_vs_postdark_postdarkplus4h_id);
comparison_max_4_3_match_tune = d2_k_max.max_dfof(1,base2_base2plus4h_vs_postdark_postdarkplus4h_id);

comparison_max_5_6_match_tune = d3_k_max.max_dfof(1,base2_base2plus4h_vs_postdark_postdarkplus4h_id);
comparison_max_6_5_match_tune = d4_k_max.max_dfof(1,base2_base2plus4h_vs_postdark_postdarkplus4h_id);



%%

%d1
comparison_dscore_prefori_d1_d2_match = double(comparison_prefori_3_4_match_tune>90);
comparison_dscore_prefori_d1_d2_match(comparison_dscore_prefori_d1_d2_match>0) = 180;
comparison_dscore_prefori_d1_d2_match = abs(comparison_dscore_prefori_d1_d2_match-comparison_prefori_3_4_match_tune);
%d2
comparison_dscore_prefori_d2_match = double(comparison_prefori_4_3_match_tune>90);
comparison_dscore_prefori_d2_match(comparison_dscore_prefori_d2_match>0) = 180;
comparison_dscore_prefori_d2_match = abs(comparison_dscore_prefori_d2_match-comparison_prefori_4_3_match_tune);

comparison_d_score_prefori_3_4 = abs(comparison_dscore_prefori_d1_d2_match-comparison_dscore_prefori_d2_match);


%d3
comparison_dscore_prefori_d3_d4_match = double(comparison_prefori_5_6_match_tune>90);
comparison_dscore_prefori_d3_d4_match(comparison_dscore_prefori_d3_d4_match>0) = 180;
comparison_dscore_prefori_d3_d4_match = abs(comparison_dscore_prefori_d3_d4_match-comparison_prefori_5_6_match_tune);
%d4
comparison_dscore_prefori_d4_match = double(comparison_prefori_6_5_match_tune>90);
comparison_dscore_prefori_d4_match(comparison_dscore_prefori_d4_match>0) = 180;
comparison_dscore_prefori_d4_match = abs(comparison_dscore_prefori_d4_match-comparison_prefori_6_5_match_tune);

comparison_d_score_prefori_5_6 = abs(comparison_dscore_prefori_d3_d4_match-comparison_dscore_prefori_d4_match);



comparison_d_score_k_3_4 = (comparison_k_3_4_match_tune-comparison_k_4_3_match_tune);


comparison_d_score_k_5_6 = (comparison_k_5_6_match_tune-comparison_k_6_5_match_tune);


comparison_d_score_max_3_4 = (comparison_max_3_4_match_tune-comparison_max_4_3_match_tune);


comparison_d_score_max_5_6 = (comparison_max_5_6_match_tune-comparison_max_6_5_match_tune);


%% EXTRA COMPARISON
%base1 to base2 (5-6) vs postdark+4hr to postdarkd (5-7)

%find cell IDs for cells that are present in both comparison sets; i.e.
%these are cells that are tuned and matched for each respective period AND
%that can be identified between periods
base1_base2_vs_postdarkplus4h_postdark7d_id = intersect(id_match_tune_1_3, id_match_tune_5_6);

%we will load info here similar to what we did above

d1 = 1; %base1 in expt list
d2 = 3; %base2 in expt list
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

d3 = 6; %postdark in expt list
d4 = 7; %postdarkd in expt list
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
extra_comparison_prefori_1_3_match_tune = d1_ori.prefOri(1,base1_base2_vs_postdarkplus4h_postdark7d_id);
extra_comparison_prefori_3_1_match_tune = d2_ori.prefOri(1,base1_base2_vs_postdarkplus4h_postdark7d_id);

extra_comparison_prefori_6_7_match_tune = d3_ori.prefOri(1,base1_base2_vs_postdarkplus4h_postdark7d_id);
extra_comparison_prefori_7_6_match_tune = d4_ori.prefOri(1,base1_base2_vs_postdarkplus4h_postdark7d_id);

%same as pref ori above but for k here and max below
extra_comparison_k_1_3_match_tune = d1_k_max.k1(1,base1_base2_vs_postdarkplus4h_postdark7d_id);
extra_comparison_k_3_1_match_tune = d2_k_max.k1(1,base1_base2_vs_postdarkplus4h_postdark7d_id);

extra_comparison_k_6_7_match_tune = d3_k_max.k1(1,base1_base2_vs_postdarkplus4h_postdark7d_id);
extra_comparison_k_7_6_match_tune = d4_k_max.k1(1,base1_base2_vs_postdarkplus4h_postdark7d_id);

extra_comparison_max_1_3_match_tune = d1_k_max.max_dfof(1,base1_base2_vs_postdarkplus4h_postdark7d_id);
extra_comparison_max_3_1_match_tune = d2_k_max.max_dfof(1,base1_base2_vs_postdarkplus4h_postdark7d_id);

extra_comparison_max_6_7_match_tune = d3_k_max.max_dfof(1,base1_base2_vs_postdarkplus4h_postdark7d_id);
extra_comparison_max_7_6_match_tune = d4_k_max.max_dfof(1,base1_base2_vs_postdarkplus4h_postdark7d_id);


%%
%now we take the info we got above and convert to change in pref ori
%scores, change in k scores, and change in max scores

%d1
extra_comparison_dscore_prefori_d1_d2_match = double(extra_comparison_prefori_1_3_match_tune>90);
extra_comparison_dscore_prefori_d1_d2_match(extra_comparison_dscore_prefori_d1_d2_match>0) = 180;
extra_comparison_dscore_prefori_d1_d2_match = abs(extra_comparison_dscore_prefori_d1_d2_match-extra_comparison_prefori_1_3_match_tune);
%d2
extra_comparison_dscore_prefori_d2_match = double(extra_comparison_prefori_3_1_match_tune>90);
extra_comparison_dscore_prefori_d2_match(extra_comparison_dscore_prefori_d2_match>0) = 180;
extra_comparison_dscore_prefori_d2_match = abs(extra_comparison_dscore_prefori_d2_match-extra_comparison_prefori_3_1_match_tune);

extra_comparison_d_score_prefori_1_3 = abs(extra_comparison_dscore_prefori_d1_d2_match-extra_comparison_dscore_prefori_d2_match);


%d3
extra_comparison_dscore_prefori_d3_d4_match = double(extra_comparison_prefori_6_7_match_tune>90);
extra_comparison_dscore_prefori_d3_d4_match(extra_comparison_dscore_prefori_d3_d4_match>0) = 180;
extra_comparison_dscore_prefori_d3_d4_match = abs(extra_comparison_dscore_prefori_d3_d4_match-extra_comparison_prefori_6_7_match_tune);
%d4
extra_comparison_dscore_prefori_d4_match = double(extra_comparison_prefori_7_6_match_tune>90);
extra_comparison_dscore_prefori_d4_match(extra_comparison_dscore_prefori_d4_match>0) = 180;
extra_comparison_dscore_prefori_d4_match = abs(extra_comparison_dscore_prefori_d4_match-extra_comparison_prefori_7_6_match_tune);

extra_comparison_d_score_prefori_6_7 = abs(extra_comparison_dscore_prefori_d3_d4_match-extra_comparison_dscore_prefori_d4_match);



extra_comparison_d_score_k_1_3 = (extra_comparison_k_1_3_match_tune-extra_comparison_k_3_1_match_tune);


extra_comparison_d_score_k_6_7 = (extra_comparison_k_6_7_match_tune-extra_comparison_k_7_6_match_tune);




extra_comparison_d_score_max_1_3 = (extra_comparison_max_1_3_match_tune-extra_comparison_max_3_1_match_tune);


extra_comparison_d_score_max_6_7 = (extra_comparison_max_6_7_match_tune-extra_comparison_max_7_6_match_tune);






%%
%now we can make some individual plots for the comparisons we want to make
%for all 3 vars (pref ori, k, max)
fig = figure;
sgtitle('Pref Ori Changes')
a=subplot(2,2,1)
h=cdfplot(comparison_d_score_prefori_1_2);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(comparison_d_score_prefori_4_5);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',1);
legend(['Base1 - Base2 + 4d (n= ', num2str(length(comparison_d_score_prefori_1_2)), ')'], ['Base2 - Post Dark (n= ', num2str(length(comparison_d_score_prefori_4_5)), ')'], 'Location', 'southeast')
xlim([0 90])
xlabel([])
ylabel([])
title([])
hold off


b=subplot(2,2,2)
h=cdfplot(comparison_d_score_prefori_3_4);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(comparison_d_score_prefori_5_6);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',1);
legend(['Base2 - Base2 + 4hr (n= ', num2str(length(comparison_d_score_prefori_3_4)), ')'], ['Post Dark - Post Light 4hr (n= ', num2str(length(comparison_d_score_prefori_5_6)), ')'], 'Location', 'southeast')
xlabel([])
ylabel([])
title([])
hold off

c=subplot(2,2,3)
h=cdfplot(comparison_d_score_prefori_1_3);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(comparison_d_score_prefori_5_7);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',1);
legend(['Base1 - Base2 (n= ', num2str(length(comparison_d_score_prefori_1_3)), ')'], ['Post Dark - Post Light 7d (n= ', num2str(length(comparison_d_score_prefori_5_7)), ')'], 'Location', 'southeast')
xlabel([])
ylabel([])
title([])
hold off

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
linkaxes([a b c],'xy')
ylabel(han,'% of Cells');
xlabel(han,'Change in Pref Ori');

print(fullfile(newfnout, [mouse,  '_comparison_prefori_change.pdf']), '-dpdf', '-bestfit')

%%
fig = figure;
sgtitle('k Changes')
a=subplot(2,2,1)
h=cdfplot(comparison_d_score_k_1_2);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(comparison_d_score_k_4_5);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',1);
legend(['Base1 - Base2 + 4d (n= ', num2str(length(comparison_d_score_k_1_2)), ')'], ['Base2 - Post Dark (n= ', num2str(length(comparison_d_score_k_4_5)), ')'], 'Location', 'southeast')
xlim([-30 30])
xlabel([])
ylabel([])
title([])
hold off


b=subplot(2,2,2)
h=cdfplot(comparison_d_score_k_3_4);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(comparison_d_score_k_5_6);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',1);
legend(['Base2 - Base2 + 4hr (n= ', num2str(length(comparison_d_score_k_3_4)), ')'], ['Post Dark - Post Light 4hr (n= ', num2str(length(comparison_d_score_k_5_6)), ')'], 'Location', 'southeast')
xlabel([])
ylabel([])
title([])
hold off

c=subplot(2,2,3)
h=cdfplot(comparison_d_score_k_1_3);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(comparison_d_score_k_5_7);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',1);
legend(['Base1 - Base2 (n= ', num2str(length(comparison_d_score_k_1_3)), ')'], ['Post Dark - Post Light 7d (n= ', num2str(length(comparison_d_score_k_5_7)), ')'], 'Location', 'southeast')
xlabel([])
ylabel([])
title([])
hold off

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
linkaxes([a b c],'xy')
ylabel(han,'% of Cells');
xlabel(han,'Change in k Value');

print(fullfile(newfnout, [mouse,  '_comparison_k_change.pdf']), '-dpdf', '-bestfit')

%%

fig = figure;
sgtitle('Max dF/F Changes')
a=subplot(2,2,1)
h=cdfplot(comparison_d_score_max_1_2);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(comparison_d_score_max_4_5);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',1);
legend(['Base1 - Base2 + 4d (n= ', num2str(length(comparison_d_score_max_1_2)), ')'], ['Base2 - Post Dark (n= ', num2str(length(comparison_d_score_max_4_5)), ')'], 'Location', 'southeast')
xlim([-1 1])
xlabel([])
ylabel([])
title([])
hold off


b=subplot(2,2,2)
h=cdfplot(comparison_d_score_max_3_4);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(comparison_d_score_max_5_6);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',1);
legend(['Base2 - Base2 + 4hr (n= ', num2str(length(comparison_d_score_max_3_4)), ')'], ['Post Dark - Post Light 4hr (n= ', num2str(length(comparison_d_score_max_5_6)), ')'], 'Location', 'southeast')
xlabel([])
ylabel([])
title([])
hold off

c=subplot(2,2,3)
h=cdfplot(comparison_d_score_max_1_3);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(comparison_d_score_max_5_7);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',1);
legend(['Base1 - Base2 (n= ', num2str(length(comparison_d_score_max_1_3)), ')'], ['Post Dark - Post Light 7d (n= ', num2str(length(comparison_d_score_max_5_7)), ')'], 'Location', 'southeast')
xlabel([])
ylabel([])
title([])
hold off

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
linkaxes([a b c],'xy')
ylabel(han,'% of Cells');
xlabel(han,'Change in Max dF/F Value');

print(fullfile(newfnout, [mouse,  '_comparison_max_change.pdf']), '-dpdf', '-bestfit')



%%
%make scatter plots of pref ori, k, and max for each pairing 
fig = figure;
a=subplot(3,2,1)
scatter(comparison_k_1_2_match_tune, comparison_k_2_1_match_tune)
xlim([0 30])
ylim([0 30])
refline(1,1)
title('Baseline1 to Baseline1+4d')
b=subplot(3,2,2)
scatter(comparison_k_4_5_match_tune, comparison_k_5_4_match_tune)
refline(1,0)
title('Baseline2 to Postdark')
c=subplot(3,2,3)
scatter(comparison_k_3_4_match_tune, comparison_k_4_3_match_tune)
refline(1,0)
title('Baseline2 to Baseline2+4hr')
d=subplot(3,2,4)
scatter(comparison_k_5_6_match_tune, comparison_k_6_5_match_tune)
refline(1,0)
title('Postdark to Postdark+4hr')
e=subplot(3,2,5)
scatter(comparison_k_1_3_match_tune, comparison_k_3_1_match_tune)
refline(1,0)
title('Baseline1 to Baseline2')
f=subplot(3,2,6)
scatter(comparison_k_5_7_match_tune, comparison_k_7_5_match_tune)
refline(1,0)
title('Postdark to Postdark7d')

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
linkaxes([a b c d e f],'xy')
ylabel(han,'Session 2 Value');
xlabel(han,'Session 1 Value');
title(han,'k Values');

print(fullfile(newfnout, [mouse,  '_comparison_k_scatter.pdf']), '-dpdf', '-bestfit')


%%
%make scatter plots of pref ori, k, and max for each pairing 
fig = figure;
a=subplot(3,2,1)
scatter(comparison_max_1_2_match_tune, comparison_max_2_1_match_tune)
xlim([0 1])
ylim([0 1])
refline(1,0)
title('Baseline1 to Baseline1+4d')
b=subplot(3,2,2)
scatter(comparison_max_4_5_match_tune, comparison_max_5_4_match_tune)
refline(1,0)
title('Baseline2 to Postdark')
c=subplot(3,2,3)
scatter(comparison_max_3_4_match_tune, comparison_max_4_3_match_tune)
refline(1,0)
title('Baseline2 to Baseline2+4hr')
d=subplot(3,2,4)
scatter(comparison_max_5_6_match_tune, comparison_max_6_5_match_tune)
refline(1,0)
title('Postdark to Postdark+4hr')
e=subplot(3,2,5)
scatter(comparison_max_1_3_match_tune, comparison_max_3_1_match_tune)
refline(1,0)
title('Baseline1 to Baseline2')
f=subplot(3,2,6)
scatter(comparison_max_5_7_match_tune, comparison_max_7_5_match_tune)
refline(1,0)
title('Postdark to Postdark7d')

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
linkaxes([a b c d e f],'xy')
ylabel(han,'Session 2 Value');
xlabel(han,'Session 1 Value');
title(han,'Max dF/F Values');

print(fullfile(newfnout, [mouse,  '_comparison_max_scatter.pdf']), '-dpdf', '-bestfit')


%%
%make scatter plots of pref ori, k, and max for each pairing 
fig = figure;
a=subplot(3,2,1)
scatter(comparison_prefori_1_2_match_tune, comparison_prefori_2_1_match_tune)
xlim([0 180])
ylim([0 180])
refline(1,1)
title('Baseline1 to Baseline1+4d')
b=subplot(3,2,2)
scatter(comparison_prefori_4_5_match_tune, comparison_prefori_5_4_match_tune)
refline(1,1)
title('Baseline2 to Postdark')
c=subplot(3,2,3)
scatter(comparison_prefori_3_4_match_tune, comparison_prefori_4_3_match_tune)
refline(1,1)
title('Baseline2 to Baseline2+4hr')
d=subplot(3,2,4)
scatter(comparison_prefori_5_6_match_tune, comparison_prefori_6_5_match_tune)
refline(1,1)
title('Postdark to Postdark+4hr')
e=subplot(3,2,5)
scatter(comparison_prefori_1_3_match_tune, comparison_prefori_3_1_match_tune)
refline(1,1)
title('Baseline1 to Baseline2')
f=subplot(3,2,6)
scatter(comparison_prefori_5_7_match_tune, comparison_prefori_7_5_match_tune)
refline(1,1)
title('Postdark to Postdark7d')

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
linkaxes([a b c d e f],'xy')
ylabel(han,'Session 2 Value');
xlabel(han,'Session 1 Value');
title(han,'Pref Ori Values');

print(fullfile(newfnout, [mouse,  '_comparison_prefori_scatter.pdf']), '-dpdf', '-bestfit')


%%
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_same_cell_dscores.mat']), 'comparison_d_score_prefori_1_2', 'comparison_d_score_prefori_1_3', 'comparison_d_score_prefori_3_4', 'comparison_d_score_prefori_4_5', 'comparison_d_score_prefori_5_6', 'comparison_d_score_prefori_5_7', ...
    'comparison_d_score_k_1_2', 'comparison_d_score_k_1_3', 'comparison_d_score_k_3_4', 'comparison_d_score_k_4_5', 'comparison_d_score_k_5_6', 'comparison_d_score_k_5_7', ...
    'comparison_d_score_max_1_2', 'comparison_d_score_max_1_3', 'comparison_d_score_max_3_4', 'comparison_d_score_max_4_5', 'comparison_d_score_max_5_6', 'comparison_d_score_max_5_7')

%%
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_same_cell_k_vals.mat']), 'comparison_k_1_2_match_tune', 'comparison_k_2_1_match_tune', 'comparison_k_4_5_match_tune', 'comparison_k_5_4_match_tune', ...
    'comparison_k_3_4_match_tune', 'comparison_k_4_3_match_tune', 'comparison_k_5_6_match_tune', 'comparison_k_6_5_match_tune', 'comparison_k_1_3_match_tune', 'comparison_k_3_1_match_tune', ...
    'comparison_k_5_7_match_tune', 'comparison_k_7_5_match_tune')

save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_same_cell_max_vals.mat']), 'comparison_max_1_2_match_tune', 'comparison_max_2_1_match_tune', 'comparison_max_4_5_match_tune', 'comparison_max_5_4_match_tune', ...
    'comparison_max_3_4_match_tune', 'comparison_max_4_3_match_tune', 'comparison_max_5_6_match_tune', 'comparison_max_6_5_match_tune', 'comparison_max_1_3_match_tune', 'comparison_max_3_1_match_tune', ...
    'comparison_max_5_7_match_tune', 'comparison_max_7_5_match_tune')

save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_same_cell_prefori_vals.mat']), 'comparison_prefori_1_2_match_tune', 'comparison_prefori_2_1_match_tune', 'comparison_prefori_4_5_match_tune', 'comparison_prefori_5_4_match_tune', ...
    'comparison_prefori_3_4_match_tune', 'comparison_prefori_4_3_match_tune', 'comparison_prefori_5_6_match_tune', 'comparison_prefori_6_5_match_tune', 'comparison_prefori_1_3_match_tune', 'comparison_prefori_3_1_match_tune', ...
    'comparison_prefori_5_7_match_tune', 'comparison_prefori_7_5_match_tune')

%%
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_extra_comp.mat']),'extra_comparison_d_score_prefori_1_3', 'extra_comparison_d_score_prefori_6_7', 'extra_comparison_d_score_k_1_3', ...
    'extra_comparison_d_score_k_6_7', 'extra_comparison_d_score_max_1_3', 'extra_comparison_d_score_max_6_7', 'extra_comparison_max_1_3_match_tune', 'extra_comparison_max_3_1_match_tune', ...
    'extra_comparison_max_6_7_match_tune', 'extra_comparison_max_7_6_match_tune')

