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
%img_point: 1 = baseline1, 2 = baseline1 + 4d, 3 = baseline2, 4 = baseline2
%+ 4hrs, 5 = post 4d dark, 6 = post dark + 4hrs, 7 = post dark + 7d

%1-7 is 2537, 8-14 is 2538
%%
%baseline1 to baseline1+4d
d1 = 1+7; %base1 in expt list
d2 = 2+7; %base1 in expt list
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


[prefori_1_2_match_tune, prefori_2_1_match_tune, d_score_prefori_1_2_match, k_1_2_match, k_2_1_match, dscore_k_1_2_match, max_1_2_match, max_2_1_match, dscore_max_1_2_match] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);

%%
%baseline1 to baseline2
d1 = 1+7; %base1 in expt list
d2 = 3+7; %base1 in expt list
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


[prefori_1_3_match_tune, prefori_3_1_match_tune, d_score_prefori_1_3_match, k_1_3_match, k_3_1_match, dscore_k_1_3_match, max_1_3_match, max_3_1_match, dscore_max_1_3_match] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);

%%
%baseline2 to baseline2+4hr
d1 = 3+7; %base1 in expt list
d2 = 4+7; %base1 in expt list
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


[prefori_3_4_match_tune, prefori_4_3_match_tune, d_score_prefori_3_4_match, k_3_4_match, k_4_3_match, dscore_k_3_4_match, max_3_4_match, max_4_3_match, dscore_max_3_4_match] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);

%%
%baseline2+4hr to postdark
d1 = 4+7; %base1 in expt list
d2 = 5+7; %base1 in expt list
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


[prefori_4_5_match_tune, prefori_5_4_match_tune, d_score_prefori_4_5_match, k_4_5_match, k_5_4_match, dscore_k_4_5_match, max_4_5_match, max_5_4_match, dscore_max_4_5_match] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);

%%
%postdark to postlight4hr
d1 = 5+7; %base1 in expt list
d2 = 6+7; %base1 in expt list
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


[prefori_5_6_match_tune, prefori_6_5_match_tune, d_score_prefori_5_6_match, k_5_6_match, k_6_5_match, dscore_k_5_6_match, max_5_6_match, max_6_5_match, dscore_max_5_6_match] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);

%%
%postdark to postlight7d
d1 = 5+7; %base1 in expt list
d2 = 7+7; %base1 in expt list
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


[prefori_5_7_match_tune, prefori_7_5_match_tune, d_score_prefori_5_7_match, k_5_7_match, k_7_5_match, dscore_k_5_7_match, max_5_7_match, max_7_5_match, dscore_max_5_7_match] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max);

%%
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_d_scores.mat']), 'd_score_prefori_1_2_match', 'd_score_prefori_1_3_match', 'd_score_prefori_3_4_match', 'd_score_prefori_4_5_match', 'd_score_prefori_5_6_match', 'd_score_prefori_5_7_match', ...
        'dscore_k_1_2_match', 'dscore_k_1_3_match', 'dscore_k_3_4_match', 'dscore_k_4_5_match', 'dscore_k_5_6_match', 'dscore_k_5_7_match', ...
        'dscore_max_1_2_match', 'dscore_max_1_3_match', 'dscore_max_3_4_match', 'dscore_max_4_5_match', 'dscore_max_5_6_match', 'dscore_max_5_7_match')

%%
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_ori_changes.mat']), 'prefori_1_2_match_tune', 'prefori_1_3_match_tune', 'prefori_2_1_match_tune', 'prefori_3_1_match_tune', ...
    'prefori_3_4_match_tune', 'prefori_4_3_match_tune', 'prefori_4_5_match_tune', 'prefori_5_4_match_tune', 'prefori_5_6_match_tune', 'prefori_5_7_match_tune', 'prefori_6_5_match_tune', ...
    'prefori_7_5_match_tune')

%%
% save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_ori_changes.mat']), 'prefori_d1_d2_match_tune', 'prefori_d2_match_tune')
% save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_d_scores.mat']), 'd_scores_all')
% 
% 
% if expt(d1).wheel == 0
%     print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_change_k_scatter.pdf']), '-dpdf', '-bestfit')
%     save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_k_changes.mat']), 'k_d1_d2_match', 'k_d1_d3_match', 'k_d2_match', 'k_d3_match')
% else
%     print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_change_k_scatter.pdf']), '-dpdf', '-bestfit')
%     save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_k_changes.mat']), 'k_d1_d2_match', 'k_d1_d3_match', 'k_d2_match', 'k_d3_match')
% end
% 
% %%
% if expt(d1).wheel == 0
%     print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_change_k_cdf.pdf']), '-dpdf', '-bestfit')
% else
%     print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW_change_k_cdf.pdf']), '-dpdf', '-bestfit')
% end
% 
% %%
% 
% if expt(d1).wheel == 0
%     print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_change_max_scatter.pdf']), '-dpdf', '-bestfit')
%     save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_max_changes.mat']), 'max_d1_d2_match', 'max_d1_d3_match', 'max_d2_match', 'max_d3_match')
% else
%     print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_change_max_scatter.pdf']), '-dpdf', '-bestfit')
%     save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_max_changes.mat']), 'max_d1_d2_match', 'max_d1_d3_match', 'max_d2_match', 'max_d3_match')
% end