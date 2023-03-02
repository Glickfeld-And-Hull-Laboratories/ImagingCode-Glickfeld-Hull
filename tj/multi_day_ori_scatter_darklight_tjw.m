%clear everything
clear all
clear all global
clc
close all

%%
%what direction is the k and max change?
%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\multi_day'; %folder to save files to
dataset = 'exp_list_darklight_actual_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%
%img_point: 1 = baseline1, 2 = baseline1 + 4d, 3 = baseline2, 4 = baseline2
%+ 4hrs, 5 = post 4d dark, 6 = post dark + 4hrs, 7 = post dark + 7d

%1-7 is 2537, 8-14 is 2538
%%
%baseline1 to baseline1+4d
d1 = 1+14; %base1 in expt list
d2 = 2+14; %base1 in expt list
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

%load multiday data for days 2 and 3
%d1_matches = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'multiday_alignment.mat']));
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));

[prefori_1_2_match_tune, prefori_2_1_match_tune, d_score_prefori_1_2_match, k_1_2_match, k_2_1_match, dscore_k_1_2_match, max_1_2_match, max_2_1_match, dscore_max_1_2_match, id_match_tune_1_2] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);


%%
%baseline1 to baseline2
d1 = 1+14; %base1 in expt list
d2 = 3+14; %base1 in expt list
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

%load multiday data for days 2 and 3
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));


[prefori_1_3_match_tune, prefori_3_1_match_tune, d_score_prefori_1_3_match, k_1_3_match, k_3_1_match, dscore_k_1_3_match, max_1_3_match, max_3_1_match, dscore_max_1_3_match, id_match_tune_1_3] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);

%%
%baseline2 to baseline2+4hr
d1 = 3+14; %base1 in expt list
d2 = 4+14; %base1 in expt list
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

%load multiday data for days 2 and 3
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));


[prefori_3_4_match_tune, prefori_4_3_match_tune, d_score_prefori_3_4_match, k_3_4_match, k_4_3_match, dscore_k_3_4_match, max_3_4_match, max_4_3_match, dscore_max_3_4_match, id_match_tune_3_4] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);

%%
%baseline2+4hr to postdark
d1 = 4+14; %base1 in expt list
d2 = 5+14; %base1 in expt list
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

%load multiday data for days 2 and 3
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));


[prefori_4_5_match_tune, prefori_5_4_match_tune, d_score_prefori_4_5_match, k_4_5_match, k_5_4_match, dscore_k_4_5_match, max_4_5_match, max_5_4_match, dscore_max_4_5_match, id_match_tune_4_5] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);

%%
%postdark to postlight4hr
d1 = 5+14; %base1 in expt list
d2 = 6+14; %base1 in expt list
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

%load multiday data for days 2 and 3
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));


[prefori_5_6_match_tune, prefori_6_5_match_tune, d_score_prefori_5_6_match, k_5_6_match, k_6_5_match, dscore_k_5_6_match, max_5_6_match, max_6_5_match, dscore_max_5_6_match, id_match_tune_5_6] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);

%%
%postdark to postlight7d
d1 = 5+14; %base1 in expt list
d2 = 7+14; %base1 in expt list
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

%load multiday data for days 2 and 3
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));


[prefori_5_7_match_tune, prefori_7_5_match_tune, d_score_prefori_5_7_match, k_5_7_match, k_7_5_match, dscore_k_5_7_match, max_5_7_match, max_7_5_match, dscore_max_5_7_match, id_match_tune_5_7] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);

%%
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_d_scores.mat']), 'd_score_prefori_1_2_match', 'd_score_prefori_1_3_match', 'd_score_prefori_3_4_match', 'd_score_prefori_4_5_match', 'd_score_prefori_5_6_match', 'd_score_prefori_5_7_match', ...
        'dscore_k_1_2_match', 'dscore_k_1_3_match', 'dscore_k_3_4_match', 'dscore_k_4_5_match', 'dscore_k_5_6_match', 'dscore_k_5_7_match', ...
        'dscore_max_1_2_match', 'dscore_max_1_3_match', 'dscore_max_3_4_match', 'dscore_max_4_5_match', 'dscore_max_5_6_match', 'dscore_max_5_7_match')

%%
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_ori_changes.mat']), 'prefori_1_2_match_tune', 'prefori_1_3_match_tune', 'prefori_2_1_match_tune', 'prefori_3_1_match_tune', ...
    'prefori_3_4_match_tune', 'prefori_4_3_match_tune', 'prefori_4_5_match_tune', 'prefori_5_4_match_tune', 'prefori_5_6_match_tune', 'prefori_5_7_match_tune', 'prefori_6_5_match_tune', ...
    'prefori_7_5_match_tune')

%%
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_id_matches.mat']), 'id_match_tune_1_2', 'id_match_tune_1_3', 'id_match_tune_3_4', 'id_match_tune_4_5', 'id_match_tune_5_6', 'id_match_tune_5_7')

%%
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
%finding only cells that were matched across multiple comparisons

%base1 to base2 vs postdark to postdark+7d
base1_base2_vs_postdark_postdark7d_id = intersect(id_match_tune_1_3, id_match_tune_5_7);

%baseline1 to baseline2
d1 = 1+14; %base1 in expt list
d2 = 3+14; %base1 in expt list
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

%load multiday data for days 2 and 3
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));

d3 = 5+14; %base1 in expt list
d4 = 7+14; %base1 in expt list
mouse = expt(d3).mouse; %mouse
ref_str_d3 = ['runs-',expt(d3).runs];
ref_str_d4 = ['runs-',expt(d4).runs]; 
img_area = expt(d3).img_loc{1};
img_layer = expt(d3).img_loc{2};
date_d3 = expt(d3).date; %day 1 (ref day) date
date_d4 = expt(d4).date; %day 2 date
img_folder_d3 = expt(d3).runs; %img folder of day 1
img_folder_d4 = expt(d4).runs; %img folder of day 2

%load ori info for each day
d3_ori = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'oriTuningInfo.mat']));
d4_ori = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'oriTuningInfo.mat']));

%load multiday data for days 2 and 3
d4_matches = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'multiday_alignment.mat']));

%load k and max vals
d3_k_max = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'k_and_max_vals.mat']));
d4_k_max = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'k_and_max_vals.mat']));


comparison_prefori_1_3_match_tune = d1_ori.prefOri(1,base1_base2_vs_postdark_postdark7d_id);
comparison_prefori_3_1_match_tune = d2_ori.prefOri(1,base1_base2_vs_postdark_postdark7d_id);

comparison_prefori_5_7_match_tune = d3_ori.prefOri(1,base1_base2_vs_postdark_postdark7d_id);
comparison_prefori_7_5_match_tune = d4_ori.prefOri(1,base1_base2_vs_postdark_postdark7d_id);

comparison_k_1_3_match_tune = d1_k_max.k1(1,base1_base2_vs_postdark_postdark7d_id);
comparison_k_3_1_match_tune = d2_k_max.k1(1,base1_base2_vs_postdark_postdark7d_id);

comparison_k_5_7_match_tune = d3_k_max.k1(1,base1_base2_vs_postdark_postdark7d_id);
comparison_k_7_5_match_tune = d4_k_max.k1(1,base1_base2_vs_postdark_postdark7d_id);

comparison_max_1_3_match_tune = d1_k_max.max_dfof(1,base1_base2_vs_postdark_postdark7d_id);
comparison_max_3_1_match_tune = d2_k_max.max_dfof(1,base1_base2_vs_postdark_postdark7d_id);

comparison_max_5_7_match_tune = d3_k_max.max_dfof(1,base1_base2_vs_postdark_postdark7d_id);
comparison_max_7_5_match_tune = d4_k_max.max_dfof(1,base1_base2_vs_postdark_postdark7d_id);






%%

%d1_2
comparison_dscore_prefori_d1_d2_match = double(comparison_prefori_1_3_match_tune>90);
comparison_dscore_prefori_d1_d2_match(comparison_dscore_prefori_d1_d2_match>0) = 180;
comparison_dscore_prefori_d1_d2_match = abs(comparison_dscore_prefori_d1_d2_match-comparison_prefori_1_3_match_tune);
%d2
comparison_dscore_prefori_d2_match = double(comparison_prefori_3_1_match_tune>90);
comparison_dscore_prefori_d2_match(comparison_dscore_prefori_d2_match>0) = 180;
comparison_dscore_prefori_d2_match = abs(comparison_dscore_prefori_d2_match-comparison_prefori_3_1_match_tune);

comparison_d_score_prefori_1_3 = abs(comparison_dscore_prefori_d1_d2_match-comparison_dscore_prefori_d2_match);


%d1_2
comparison_dscore_prefori_d3_d4_match = double(comparison_prefori_5_7_match_tune>90);
comparison_dscore_prefori_d3_d4_match(comparison_dscore_prefori_d3_d4_match>0) = 180;
comparison_dscore_prefori_d3_d4_match = abs(comparison_dscore_prefori_d3_d4_match-comparison_prefori_5_7_match_tune);
%d2
comparison_dscore_prefori_d4_match = double(comparison_prefori_7_5_match_tune>90);
comparison_dscore_prefori_d4_match(comparison_dscore_prefori_d4_match>0) = 180;
comparison_dscore_prefori_d4_match = abs(comparison_dscore_prefori_d4_match-comparison_prefori_7_5_match_tune);

comparison_d_score_prefori_5_7 = abs(comparison_dscore_prefori_d3_d4_match-comparison_dscore_prefori_d4_match);


%d1_2
comparison_dscore_k_d1_d2_match = double(comparison_k_1_3_match_tune>90);
comparison_dscore_k_d1_d2_match(comparison_dscore_k_d1_d2_match>0) = 180;
comparison_dscore_k_d1_d2_match = abs(comparison_dscore_k_d1_d2_match-comparison_k_1_3_match_tune);
%d2
comparison_dscore_k_d2_match = double(comparison_k_3_1_match_tune>90);
comparison_dscore_k_d2_match(comparison_dscore_k_d2_match>0) = 180;
comparison_dscore_k_d2_match = abs(comparison_dscore_k_d2_match-comparison_k_3_1_match_tune);

comparison_d_score_k_1_3 = abs(comparison_dscore_k_d1_d2_match-comparison_dscore_k_d2_match);


%d1_2
comparison_dscore_k_d3_d4_match = double(comparison_k_5_7_match_tune>90);
comparison_dscore_k_d3_d4_match(comparison_dscore_k_d3_d4_match>0) = 180;
comparison_dscore_k_d3_d4_match = abs(comparison_dscore_k_d3_d4_match-comparison_k_5_7_match_tune);
%d2
comparison_dscore_k_d4_match = double(comparison_k_7_5_match_tune>90);
comparison_dscore_k_d4_match(comparison_dscore_k_d4_match>0) = 180;
comparison_dscore_k_d4_match = abs(comparison_dscore_k_d4_match-comparison_k_7_5_match_tune);

comparison_d_score_k_5_7 = abs(comparison_dscore_k_d3_d4_match-comparison_dscore_k_d4_match);



%d1_2
comparison_dscore_max_d1_d2_match = double(comparison_max_1_3_match_tune>90);
comparison_dscore_max_d1_d2_match(comparison_dscore_max_d1_d2_match>0) = 180;
comparison_dscore_max_d1_d2_match = abs(comparison_dscore_max_d1_d2_match-comparison_max_1_3_match_tune);
%d2
comparison_dscore_max_d2_match = double(comparison_max_3_1_match_tune>90);
comparison_dscore_max_d2_match(comparison_dscore_max_d2_match>0) = 180;
comparison_dscore_max_d2_match = abs(comparison_dscore_max_d2_match-comparison_max_3_1_match_tune);

comparison_d_score_max_1_3 = abs(comparison_dscore_max_d1_d2_match-comparison_dscore_max_d2_match);


%d1_2
comparison_dscore_max_d3_d4_match = double(comparison_max_5_7_match_tune>90);
comparison_dscore_max_d3_d4_match(comparison_dscore_max_d3_d4_match>0) = 180;
comparison_dscore_max_d3_d4_match = abs(comparison_dscore_max_d3_d4_match-comparison_max_5_7_match_tune);
%d2
comparison_dscore_max_d4_match = double(comparison_max_7_5_match_tune>90);
comparison_dscore_max_d4_match(comparison_dscore_max_d4_match>0) = 180;
comparison_dscore_max_d4_match = abs(comparison_dscore_max_d4_match-comparison_max_7_5_match_tune);

comparison_d_score_max_5_7 = abs(comparison_dscore_max_d3_d4_match-comparison_dscore_max_d4_match);





figure;
cdfplot(comparison_d_score_prefori_1_3);
hold on
cdfplot(comparison_d_score_prefori_5_7);
hold off
legend('1-3','5-7')
title('Pref Ori Change')

figure;
cdfplot(comparison_d_score_k_1_3);
hold on
cdfplot(comparison_d_score_k_5_7);
hold off
legend('1-3','5-7')
title('k Change')

figure;
cdfplot(comparison_d_score_max_1_3);
hold on
cdfplot(comparison_d_score_max_5_7);
hold off
legend('1-3','5-7')
title('Max dF/F Change')

%%
%base1 to base1+4hr vs base2+4hr to postdark
base1_base1plus4d_vs_base2plus4h_postdark_id = intersect(id_match_tune_1_2, id_match_tune_4_5);

%baseline1 to baseline2
d1 = 1+14; %base1 in expt list
d2 = 2+14; %base1 in expt list
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

%load multiday data for days 2 and 3
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));

d3 = 4+14; %base1 in expt list
d4 = 5+14; %base1 in expt list
mouse = expt(d3).mouse; %mouse
ref_str_d3 = ['runs-',expt(d3).runs];
ref_str_d4 = ['runs-',expt(d4).runs]; 
img_area = expt(d3).img_loc{1};
img_layer = expt(d3).img_loc{2};
date_d3 = expt(d3).date; %day 1 (ref day) date
date_d4 = expt(d4).date; %day 2 date
img_folder_d3 = expt(d3).runs; %img folder of day 1
img_folder_d4 = expt(d4).runs; %img folder of day 2

%load ori info for each day
d3_ori = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'oriTuningInfo.mat']));
d4_ori = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'oriTuningInfo.mat']));

%load multiday data for days 2 and 3
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

%d1_2
comparison_dscore_prefori_d1_d2_match = double(comparison_prefori_1_2_match_tune>90);
comparison_dscore_prefori_d1_d2_match(comparison_dscore_prefori_d1_d2_match>0) = 180;
comparison_dscore_prefori_d1_d2_match = abs(comparison_dscore_prefori_d1_d2_match-comparison_prefori_1_2_match_tune);
%d2
comparison_dscore_prefori_d2_match = double(comparison_prefori_2_1_match_tune>90);
comparison_dscore_prefori_d2_match(comparison_dscore_prefori_d2_match>0) = 180;
comparison_dscore_prefori_d2_match = abs(comparison_dscore_prefori_d2_match-comparison_prefori_2_1_match_tune);

comparison_d_score_prefori_1_2 = abs(comparison_dscore_prefori_d1_d2_match-comparison_dscore_prefori_d2_match);


%d1_2
comparison_dscore_prefori_d3_d4_match = double(comparison_prefori_4_5_match_tune>90);
comparison_dscore_prefori_d3_d4_match(comparison_dscore_prefori_d3_d4_match>0) = 180;
comparison_dscore_prefori_d3_d4_match = abs(comparison_dscore_prefori_d3_d4_match-comparison_prefori_4_5_match_tune);
%d2
comparison_dscore_prefori_d4_match = double(comparison_prefori_5_4_match_tune>90);
comparison_dscore_prefori_d4_match(comparison_dscore_prefori_d4_match>0) = 180;
comparison_dscore_prefori_d4_match = abs(comparison_dscore_prefori_d4_match-comparison_prefori_5_4_match_tune);

comparison_d_score_prefori_4_5 = abs(comparison_dscore_prefori_d3_d4_match-comparison_dscore_prefori_d4_match);


%d1_2
comparison_dscore_k_d1_d2_match = double(comparison_k_1_2_match_tune>90);
comparison_dscore_k_d1_d2_match(comparison_dscore_k_d1_d2_match>0) = 180;
comparison_dscore_k_d1_d2_match = abs(comparison_dscore_k_d1_d2_match-comparison_k_1_2_match_tune);
%d2
comparison_dscore_k_d2_match = double(comparison_k_2_1_match_tune>90);
comparison_dscore_k_d2_match(comparison_dscore_k_d2_match>0) = 180;
comparison_dscore_k_d2_match = abs(comparison_dscore_k_d2_match-comparison_k_2_1_match_tune);

comparison_d_score_k_1_2 = abs(comparison_dscore_k_d1_d2_match-comparison_dscore_k_d2_match);


%d1_2
comparison_dscore_k_d3_d4_match = double(comparison_k_4_5_match_tune>90);
comparison_dscore_k_d3_d4_match(comparison_dscore_k_d3_d4_match>0) = 180;
comparison_dscore_k_d3_d4_match = abs(comparison_dscore_k_d3_d4_match-comparison_k_4_5_match_tune);
%d2
comparison_dscore_k_d4_match = double(comparison_k_5_4_match_tune>90);
comparison_dscore_k_d4_match(comparison_dscore_k_d4_match>0) = 180;
comparison_dscore_k_d4_match = abs(comparison_dscore_k_d4_match-comparison_k_5_4_match_tune);

comparison_d_score_k_4_5 = abs(comparison_dscore_k_d3_d4_match-comparison_dscore_k_d4_match);



%d1_2
comparison_dscore_max_d1_d2_match = double(comparison_max_1_2_match_tune>90);
comparison_dscore_max_d1_d2_match(comparison_dscore_max_d1_d2_match>0) = 180;
comparison_dscore_max_d1_d2_match = abs(comparison_dscore_max_d1_d2_match-comparison_max_1_2_match_tune);
%d2
comparison_dscore_max_d2_match = double(comparison_max_2_1_match_tune>90);
comparison_dscore_max_d2_match(comparison_dscore_max_d2_match>0) = 180;
comparison_dscore_max_d2_match = abs(comparison_dscore_max_d2_match-comparison_max_2_1_match_tune);

comparison_d_score_max_1_2 = abs(comparison_dscore_max_d1_d2_match-comparison_dscore_max_d2_match);


%d1_2
comparison_dscore_max_d3_d4_match = double(comparison_max_4_5_match_tune>90);
comparison_dscore_max_d3_d4_match(comparison_dscore_max_d3_d4_match>0) = 180;
comparison_dscore_max_d3_d4_match = abs(comparison_dscore_max_d3_d4_match-comparison_max_4_5_match_tune);
%d2
comparison_dscore_max_d4_match = double(comparison_max_5_4_match_tune>90);
comparison_dscore_max_d4_match(comparison_dscore_max_d4_match>0) = 180;
comparison_dscore_max_d4_match = abs(comparison_dscore_max_d4_match-comparison_max_5_4_match_tune);

comparison_d_score_max_4_5 = abs(comparison_dscore_max_d3_d4_match-comparison_dscore_max_d4_match);




figure;
cdfplot(comparison_d_score_prefori_1_2);
hold on
cdfplot(comparison_d_score_prefori_4_5);
hold off
title('Pref Ori Change')
legend('1-2', '4-5')

figure;
cdfplot(comparison_d_score_k_1_2);
hold on
cdfplot(comparison_d_score_k_4_5);
hold off
title('k Change')
legend('1-2', '4-5')


figure;
cdfplot(comparison_d_score_max_1_2);
hold on
cdfplot(comparison_d_score_max_4_5);
hold off
title('Max dF/F Change')
legend('1-2', '4-5')

%%
%base2 to base2+4hr vs postdark to postdark+4hr
base2_base2plus4h_vs_postdark_postdarkplus4h_id = intersect(id_match_tune_3_4, id_match_tune_5_6);

%baseline1 to baseline2
d1 = 3+14; %base1 in expt list
d2 = 4+14; %base1 in expt list
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

%load multiday data for days 2 and 3
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));

d3 = 5+14; %base1 in expt list
d4 = 6+14; %base1 in expt list
mouse = expt(d3).mouse; %mouse
ref_str_d3 = ['runs-',expt(d3).runs];
ref_str_d4 = ['runs-',expt(d4).runs]; 
img_area = expt(d3).img_loc{1};
img_layer = expt(d3).img_loc{2};
date_d3 = expt(d3).date; %day 1 (ref day) date
date_d4 = expt(d4).date; %day 2 date
img_folder_d3 = expt(d3).runs; %img folder of day 1
img_folder_d4 = expt(d4).runs; %img folder of day 2

%load ori info for each day
d3_ori = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'oriTuningInfo.mat']));
d4_ori = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'oriTuningInfo.mat']));

%load multiday data for days 2 and 3
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

%d1_2
comparison_dscore_prefori_d1_d2_match = double(comparison_prefori_3_4_match_tune>90);
comparison_dscore_prefori_d1_d2_match(comparison_dscore_prefori_d1_d2_match>0) = 180;
comparison_dscore_prefori_d1_d2_match = abs(comparison_dscore_prefori_d1_d2_match-comparison_prefori_3_4_match_tune);
%d2
comparison_dscore_prefori_d2_match = double(comparison_prefori_4_3_match_tune>90);
comparison_dscore_prefori_d2_match(comparison_dscore_prefori_d2_match>0) = 180;
comparison_dscore_prefori_d2_match = abs(comparison_dscore_prefori_d2_match-comparison_prefori_4_3_match_tune);

comparison_d_score_prefori_3_4 = abs(comparison_dscore_prefori_d1_d2_match-comparison_dscore_prefori_d2_match);


%d1_2
comparison_dscore_prefori_d3_d4_match = double(comparison_prefori_5_6_match_tune>90);
comparison_dscore_prefori_d3_d4_match(comparison_dscore_prefori_d3_d4_match>0) = 180;
comparison_dscore_prefori_d3_d4_match = abs(comparison_dscore_prefori_d3_d4_match-comparison_prefori_5_6_match_tune);
%d2
comparison_dscore_prefori_d4_match = double(comparison_prefori_6_5_match_tune>90);
comparison_dscore_prefori_d4_match(comparison_dscore_prefori_d4_match>0) = 180;
comparison_dscore_prefori_d4_match = abs(comparison_dscore_prefori_d4_match-comparison_prefori_6_5_match_tune);

comparison_d_score_prefori_5_6 = abs(comparison_dscore_prefori_d3_d4_match-comparison_dscore_prefori_d4_match);


%d1_2
comparison_dscore_k_d1_d2_match = double(comparison_k_3_4_match_tune>90);
comparison_dscore_k_d1_d2_match(comparison_dscore_k_d1_d2_match>0) = 180;
comparison_dscore_k_d1_d2_match = abs(comparison_dscore_k_d1_d2_match-comparison_k_3_4_match_tune);
%d2
comparison_dscore_k_d2_match = double(comparison_k_4_3_match_tune>90);
comparison_dscore_k_d2_match(comparison_dscore_k_d2_match>0) = 180;
comparison_dscore_k_d2_match = abs(comparison_dscore_k_d2_match-comparison_k_4_3_match_tune);

comparison_d_score_k_3_4 = abs(comparison_dscore_k_d1_d2_match-comparison_dscore_k_d2_match);


%d1_2
comparison_dscore_k_d3_d4_match = double(comparison_k_5_6_match_tune>90);
comparison_dscore_k_d3_d4_match(comparison_dscore_k_d3_d4_match>0) = 180;
comparison_dscore_k_d3_d4_match = abs(comparison_dscore_k_d3_d4_match-comparison_k_5_6_match_tune);
%d2
comparison_dscore_k_d4_match = double(comparison_k_6_5_match_tune>90);
comparison_dscore_k_d4_match(comparison_dscore_k_d4_match>0) = 180;
comparison_dscore_k_d4_match = abs(comparison_dscore_k_d4_match-comparison_k_6_5_match_tune);

comparison_d_score_k_5_6 = abs(comparison_dscore_k_d3_d4_match-comparison_dscore_k_d4_match);


%d1_2
comparison_dscore_max_d1_d2_match = double(comparison_max_3_4_match_tune>90);
comparison_dscore_max_d1_d2_match(comparison_dscore_max_d1_d2_match>0) = 180;
comparison_dscore_max_d1_d2_match = abs(comparison_dscore_max_d1_d2_match-comparison_max_3_4_match_tune);
%d2
comparison_dscore_max_d2_match = double(comparison_max_4_3_match_tune>90);
comparison_dscore_max_d2_match(comparison_dscore_max_d2_match>0) = 180;
comparison_dscore_max_d2_match = abs(comparison_dscore_max_d2_match-comparison_max_4_3_match_tune);

comparison_d_score_max_3_4 = abs(comparison_dscore_max_d1_d2_match-comparison_dscore_max_d2_match);


%d1_2
comparison_dscore_max_d3_d4_match = double(comparison_max_5_6_match_tune>90);
comparison_dscore_max_d3_d4_match(comparison_dscore_max_d3_d4_match>0) = 180;
comparison_dscore_max_d3_d4_match = abs(comparison_dscore_max_d3_d4_match-comparison_max_5_6_match_tune);
%d2
comparison_dscore_max_d4_match = double(comparison_max_6_5_match_tune>90);
comparison_dscore_max_d4_match(comparison_dscore_max_d4_match>0) = 180;
comparison_dscore_max_d4_match = abs(comparison_dscore_max_d4_match-comparison_max_6_5_match_tune);

comparison_d_score_max_5_6 = abs(comparison_dscore_max_d3_d4_match-comparison_dscore_max_d4_match);



figure;
cdfplot(comparison_d_score_prefori_3_4);
hold on
cdfplot(comparison_d_score_prefori_5_6);
hold off
title('Pref Ori Change')
legend('3-4', '5-6')

figure;
cdfplot(comparison_d_score_k_3_4);
hold on
cdfplot(comparison_d_score_k_5_6);
hold off
title('k Change')
legend('3-4', '5-6')

figure;
cdfplot(comparison_d_score_max_3_4);
hold on
cdfplot(comparison_d_score_max_5_6);
hold off
title('Max dF/F Change')
legend('3-4', '5-6')


%%
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_same_cell_dscores.mat']), 'comparison_d_score_prefori_1_2', 'comparison_d_score_prefori_1_3', 'comparison_d_score_prefori_3_4', 'comparison_d_score_prefori_4_5', 'comparison_d_score_prefori_5_6', 'comparison_d_score_prefori_5_7', ...
    'comparison_d_score_k_1_2', 'comparison_d_score_k_1_3', 'comparison_d_score_k_3_4', 'comparison_d_score_k_4_5', 'comparison_d_score_k_5_6', 'comparison_d_score_k_5_7', ...
    'comparison_d_score_max_1_2', 'comparison_d_score_max_1_3', 'comparison_d_score_max_3_4', 'comparison_d_score_max_4_5', 'comparison_d_score_max_5_6', 'comparison_d_score_max_5_7')
