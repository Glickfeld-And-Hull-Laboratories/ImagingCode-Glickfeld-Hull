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

d1 = 46; 
d2 = 47; 
d3 = 48;
d4 = 49;
mouse = expt(d1).mouse; 
ref_str_d1 = ['runs-',expt(d1).runs];
ref_str_d2 = ['runs-',expt(d2).runs]; 
ref_str_d3 = ['runs-',expt(d3).runs]; 
ref_str_d4 = ['runs-',expt(d4).runs];
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date_d1 = expt(d1).date; 
date_d2 = expt(d2).date; 
date_d3 = expt(d3).date; 
date_d4 = expt(d4).date; 
img_folder_d1 = expt(d1).runs; 
img_folder_d2 = expt(d2).runs; 
img_folder_d3 = expt(d3).runs; 
img_folder_d4 = expt(d4).runs; 

%load ori info for each day
d1_ori = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningInfo.mat']));
d2_ori = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningInfo.mat']));
d3_ori = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'oriTuningInfo.mat']));
d4_ori = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'oriTuningInfo.mat']));

d1_fits = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningAndFits.mat']));
d2_fits = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningAndFits.mat']));
d3_fits = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'oriTuningAndFits.mat']));
d4_fits = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'oriTuningAndFits.mat']));


%load multiday data for day 2
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));
d3_matches = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'multiday_alignment.mat']));
d4_matches = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));
d3_k_max = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'k_and_max_vals.mat']));
d4_k_max = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'k_and_max_vals.mat']));

%% matching indices

tuned_d1 = d1_ori.ind_theta90;
tuned_d2 = d2_ori.ind_theta90;
tuned_d3 = d3_ori.ind_theta90;
tuned_d4 = d4_ori.ind_theta90;

tuned_all = unique([tuned_d1 tuned_d2 tuned_d3 tuned_d4]);

match_d2 = find([d2_matches.cellImageAlign.pass]); 
match_d3 = find([d3_matches.cellImageAlign.pass]); 
match_d4 = find([d4_matches.cellImageAlign.pass]); 

match_all = intersect(intersect(intersect(match_d2,match_d3),match_d4),tuned_all);


%%
%fit reliability

d1_reliability = d1_fits.fitReliability(match_all);
d2_reliability = d2_fits.fitReliability(match_all);
d3_reliability = d3_fits.fitReliability(match_all);
d4_reliability = d4_fits.fitReliability(match_all);


%%

d1_prefori = d1_ori.prefOri(1,match_all);
d2_prefori = d2_ori.prefOri(1,match_all);
d3_prefori = d3_ori.prefOri(1,match_all);
d4_prefori = d4_ori.prefOri(1,match_all);

dscore_prefori_1 = double(d1_prefori>90);
dscore_prefori_1(dscore_prefori_1>0) = 180;
dscore_prefori_1 = abs(dscore_prefori_1-d1_prefori);

dscore_prefori_2 = double(d2_prefori>90);
dscore_prefori_2(dscore_prefori_2>0) = 180;
dscore_prefori_2 = abs(dscore_prefori_2-d2_prefori);

dscore_prefori_3 = double(d3_prefori>90);
dscore_prefori_3(dscore_prefori_3>0) = 180;
dscore_prefori_3 = abs(dscore_prefori_3-d3_prefori);

dscore_prefori_4 = double(d4_prefori>90);
dscore_prefori_4(dscore_prefori_4>0) = 180;
dscore_prefori_4 = abs(dscore_prefori_4-d4_prefori);


d_score_prefori_d1_d2 = abs(dscore_prefori_1-dscore_prefori_2);
d_score_prefori_d2_d3 = abs(dscore_prefori_2-dscore_prefori_3);
d_score_prefori_d3_d4 = abs(dscore_prefori_3-dscore_prefori_4);
d_score_prefori_d1_d3 = abs(dscore_prefori_1-dscore_prefori_3);
d_score_prefori_d1_d4 = abs(dscore_prefori_1-dscore_prefori_4);


fig = figure;
sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(d_score_prefori_d1_d2);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d_score_prefori_d2_d3);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',1);
l=cdfplot(d_score_prefori_d3_d4);
set(l, 'LineStyle', '-', 'Color', [.8 .8 .8], 'LineWidth',1);
legend(['Base1 - Base2'], ['Base2 - Post-dark'], ['Post-dark - Post-dark+7d'], 'Location', 'southeast')
xlim([0 90])
xlabel([])
ylabel([])
title([])
hold off

print(fullfile(newfnout, [mouse,  '_prefori_change.pdf']), '-dpdf', '-bestfit')

%%
d1_k = d1_k_max.k1(1,match_all);
d2_k = d2_k_max.k1(1,match_all);
d3_k = d3_k_max.k1(1,match_all);
d4_k = d4_k_max.k1(1,match_all);

fig = figure;
sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(d1_k);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d2_k);
set(j, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',1);
l=cdfplot(d3_k);
set(l, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',1);
m=cdfplot(d4_k);
set(m, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',1);
legend(['Base1'], ['Base2'], ['Post-dark'], ['Post-dark+7d'], 'Location', 'southeast')
xlim([0 30])
xlabel([])
ylabel([])
title([])
hold off

print(fullfile(newfnout, [mouse,  '_k_vals.pdf']), '-dpdf', '-bestfit')


%%
d1_max = d1_k_max.max_dfof(1,match_all);
d2_max = d2_k_max.max_dfof(1,match_all);
d3_max = d3_k_max.max_dfof(1,match_all);
d4_max = d4_k_max.max_dfof(1,match_all);

fig = figure;
sgtitle('Max dF/F values')
% a=subplot(2,2,1)
h=cdfplot(d1_max);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d2_max);
set(j, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',1);
l=cdfplot(d3_max);
set(l, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',1);
m=cdfplot(d4_max);
set(m, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',1);
legend(['Base1'], ['Base2'], ['Post-dark'], ['Post-dark+7d'], 'Location', 'southeast')
xlim([0 1])
xlabel([])
ylabel([])
title([])
hold off

print(fullfile(newfnout, [mouse,  '_maxdfof_vals.pdf']), '-dpdf', '-bestfit')



%%
%save all of the change scores from above for future use
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_d_scores.mat']), 'd_score_prefori_d1_d2', 'd_score_prefori_d2_d3', 'd_score_prefori_d3_d4', ...
    'd_score_prefori_d1_d3', 'd_score_prefori_d1_d4')

%save pref ori info
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_ori_changes.mat']), 'd1_prefori', 'd2_prefori', 'd3_prefori', 'd4_prefori')

%save tuned and matched cell IDs
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_id_matches.mat']), 'match_all')

%save k and max vals
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_k_max_changes.mat']), 'd1_k', 'd2_k', 'd3_k', 'd4_k', 'd1_max', 'd2_max', 'd3_max', 'd4_max')

%save fit reliabilities
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_fit_reliability.mat']), 'd1_reliability', 'd2_reliability', 'd3_reliability', 'd4_reliability')

%%
max_all = [d1_max; d2_max; d3_max; d4_max];
grand_max = max(max_all);
grand_med = median(grand_max);

bot_50_d_score_prefori_d1_d2 = d_score_prefori_d1_d2(grand_max < grand_med);
top_50_d_score_prefori_d1_d2 = d_score_prefori_d1_d2(grand_max > grand_med);
bot_50_d_score_prefori_d2_d3 = d_score_prefori_d2_d3(grand_max < grand_med);
top_50_d_score_prefori_d2_d3 = d_score_prefori_d2_d3(grand_max > grand_med);
bot_50_d_score_prefori_d3_d4 = d_score_prefori_d3_d4(grand_max < grand_med);
top_50_d_score_prefori_d3_d4 = d_score_prefori_d3_d4(grand_max > grand_med);
bot_50_d_score_prefori_d1_d3 = d_score_prefori_d1_d3(grand_max < grand_med);
top_50_d_score_prefori_d1_d3 = d_score_prefori_d1_d3(grand_max > grand_med);
bot_50_d_score_prefori_d1_d4 = d_score_prefori_d1_d4(grand_max < grand_med);
top_50_d_score_prefori_d1_d4 = d_score_prefori_d1_d4(grand_max > grand_med);


save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_median_sep_d_scores.mat']), 'bot_50_d_score_prefori_d1_d2', 'top_50_d_score_prefori_d1_d2', 'bot_50_d_score_prefori_d2_d3', ...
    'top_50_d_score_prefori_d2_d3', 'bot_50_d_score_prefori_d3_d4', 'top_50_d_score_prefori_d3_d4', 'bot_50_d_score_prefori_d1_d3', 'top_50_d_score_prefori_d1_d3', ...
    'bot_50_d_score_prefori_d1_d4', 'top_50_d_score_prefori_d1_d4')

%%
k_all = [d1_k; d2_k; d3_k; d4_k];
grand_k = max(k_all);
grand_med_k = median(grand_k);

wide_d_score_prefori_d1_d2 = d_score_prefori_d1_d2(grand_k < grand_med_k);
sharp_d_score_prefori_d1_d2 = d_score_prefori_d1_d2(grand_k > grand_med_k);
wide_d_score_prefori_d2_d3 = d_score_prefori_d2_d3(grand_k < grand_med_k);
sharp_d_score_prefori_d2_d3 = d_score_prefori_d2_d3(grand_k > grand_med_k);
wide_d_score_prefori_d3_d4 = d_score_prefori_d3_d4(grand_k < grand_med_k);
sharp_d_score_prefori_d3_d4 = d_score_prefori_d3_d4(grand_k > grand_med_k);
wide_d_score_prefori_d1_d3 = d_score_prefori_d1_d3(grand_k < grand_med_k);
sharp_d_score_prefori_d1_d3 = d_score_prefori_d1_d3(grand_k > grand_med_k);
wide_d_score_prefori_d1_d4 = d_score_prefori_d1_d4(grand_k < grand_med_k);
sharp_d_score_prefori_d1_d4 = d_score_prefori_d1_d4(grand_k > grand_med_k);

save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_wide_sharp_sep_d_scores.mat']), 'wide_d_score_prefori_d1_d2', 'sharp_d_score_prefori_d1_d2', 'wide_d_score_prefori_d2_d3', ...
    'sharp_d_score_prefori_d2_d3', 'wide_d_score_prefori_d3_d4', 'sharp_d_score_prefori_d3_d4', 'wide_d_score_prefori_d1_d3', 'sharp_d_score_prefori_d1_d3', ...
    'wide_d_score_prefori_d1_d4', 'sharp_d_score_prefori_d1_d4')




