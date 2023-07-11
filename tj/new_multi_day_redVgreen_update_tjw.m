%This code is used to analyze individual subject data from the dark to
%light housing 2p imaging experiment, so that these data are prepared for
%pooled data afterwards.
%%
%clear everything
clear all
clear all global
clc
close all

%% REMEMBER TO CHANGE THE OUT FOLDER BASED ON TUNING CRITERION - NORMAL FOR TUNED ON 1/3 SESSIONS and STRICT FOR TUNED ON ALL 3 ***
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_Arc_greenVred'; %folder to save files to
% newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_Arc_greenVred_STRICT'; %folder to save files to
dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%

%we start with loading the proper information based on the mouse's info in
%the experiment list

d1 = 79; 
d2 = 80; 
d3 = 81;

mouse = expt(d1).mouse; 
group = expt(d1).red_indicator{1};
ref_str_d1 = ['runs-',expt(d1).runs];
ref_str_d2 = ['runs-',expt(d2).runs]; 
ref_str_d3 = ['runs-',expt(d3).runs]; 
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date_d1 = expt(d1).date; 
date_d2 = expt(d2).date; 
date_d3 = expt(d3).date; 
img_folder_d1 = expt(d1).runs; 
img_folder_d2 = expt(d2).runs; 
img_folder_d3 = expt(d3).runs; 

%load ori info for each day
d1_ori = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningInfo.mat']));
d2_ori = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningInfo.mat']));
d3_ori = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'oriTuningInfo.mat']));

d1_fits = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningAndFits.mat']));
d2_fits = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningAndFits.mat']));
d3_fits = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'oriTuningAndFits.mat']));


%load multiday data for day 2
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));
d3_matches = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));
d3_k_max = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'k_and_max_vals.mat']));

%r or gr cell info
d2_matches_tj = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'tj_matches.mat']));
d3_matches_tj = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'tj_matches.mat']));

%significantly responsive cells
d1_sigresp = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'trialData.mat']));
d2_sigresp = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'trialData.mat']));
d3_sigresp = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'trialData.mat']));

%% matching indices *** REMEMBER TO CHANGE THE TUNING CRITERION BASED ON LOOSE (1/3 session) OR STRICT (all 3 session) ***

tuned_d1 = d1_ori.ind_theta90;
tuned_d2 = d2_ori.ind_theta90;
tuned_d3 = d3_ori.ind_theta90;

sig_resp_d1 = d1_sigresp.good_ind;
sig_resp_d2 = d2_sigresp.good_ind;
sig_resp_d3 = d3_sigresp.good_ind;

% tuned_all = unique([tuned_d1 tuned_d2 tuned_d3]);
tuned_all = intersect(intersect(tuned_d1, tuned_d2),tuned_d3);

%find green cells that match from d1 to d2 and d1 to d3 
green_match_d2 = d2_matches_tj.green_match_ind; 
green_match_d3 = d3_matches_tj.green_match_ind; 

%find red cells that match from d1 to d2 and d1 to d3 
red_match_d2 = d2_matches_tj.red_match_ind; 
red_match_d3 = d3_matches_tj.red_match_ind; 

green_match_all = intersect(intersect(green_match_d2,green_match_d3),tuned_all);
red_match_all = intersect(intersect(red_match_d2,red_match_d3),tuned_all);

green_match_sig_d1 = intersect(intersect(green_match_d2,green_match_d3),sig_resp_d1)';
red_match_sig_d1 = intersect(intersect(red_match_d2,red_match_d3),sig_resp_d1)';

%%
%fit reliability

green_d1_reliability = d1_fits.fitReliability(green_match_all);
green_d2_reliability = d2_fits.fitReliability(green_match_all);
green_d3_reliability = d3_fits.fitReliability(green_match_all);

red_d1_reliability = d1_fits.fitReliability(red_match_all);
red_d2_reliability = d2_fits.fitReliability(red_match_all);
red_d3_reliability = d3_fits.fitReliability(red_match_all);

green_d1_matchsig_reliability = d1_fits.fitReliability(green_match_sig_d1);
red_d1_matchsig_reliability = d1_fits.fitReliability(red_match_sig_d1);

if str2double(expt(d1).img_day) == 1
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_matchsig_fits.mat']), 'green_d1_matchsig_reliability', 'red_d1_matchsig_reliability')
else
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_matchsig_fits_2.mat']), 'green_d1_matchsig_reliability', 'red_d1_matchsig_reliability')
end

%%

green_d1_prefori = d1_ori.prefOri(1,green_match_all);
green_d2_prefori = d2_ori.prefOri(1,green_match_all);
green_d3_prefori = d3_ori.prefOri(1,green_match_all);

green_dscore_prefori_1 = double(green_d1_prefori>90);
green_dscore_prefori_1(green_dscore_prefori_1>0) = 180;
green_dscore_prefori_1 = abs(green_dscore_prefori_1-green_d1_prefori);

green_dscore_prefori_2 = double(green_d2_prefori>90);
green_dscore_prefori_2(green_dscore_prefori_2>0) = 180;
green_dscore_prefori_2 = abs(green_dscore_prefori_2-green_d2_prefori);

green_dscore_prefori_3 = double(green_d3_prefori>90);
green_dscore_prefori_3(green_dscore_prefori_3>0) = 180;
green_dscore_prefori_3 = abs(green_dscore_prefori_3-green_d3_prefori);

green_d_score_prefori_d1_d2 = abs(green_dscore_prefori_1-green_dscore_prefori_2);
green_d_score_prefori_d1_d3 = abs(green_dscore_prefori_1-green_dscore_prefori_3);


red_d1_prefori = d1_ori.prefOri(1,red_match_all);
red_d2_prefori = d2_ori.prefOri(1,red_match_all);
red_d3_prefori = d3_ori.prefOri(1,red_match_all);

red_dscore_prefori_1 = double(red_d1_prefori>90);
red_dscore_prefori_1(red_dscore_prefori_1>0) = 180;
red_dscore_prefori_1 = abs(red_dscore_prefori_1-red_d1_prefori);

red_dscore_prefori_2 = double(red_d2_prefori>90);
red_dscore_prefori_2(red_dscore_prefori_2>0) = 180;
red_dscore_prefori_2 = abs(red_dscore_prefori_2-red_d2_prefori);

red_dscore_prefori_3 = double(red_d3_prefori>90);
red_dscore_prefori_3(red_dscore_prefori_3>0) = 180;
red_dscore_prefori_3 = abs(red_dscore_prefori_3-red_d3_prefori);

red_d_score_prefori_d1_d2 = abs(red_dscore_prefori_1-red_dscore_prefori_2);
red_d_score_prefori_d1_d3 = abs(red_dscore_prefori_1-red_dscore_prefori_3);


fig = figure;
sgtitle(['Pref Ori Changes ', mouse, ' ', group])
% a=subplot(2,2,1)
h=cdfplot(green_d_score_prefori_d1_d2);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',1);
hold on
j=cdfplot(green_d_score_prefori_d1_d3);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',1);
l=cdfplot(red_d_score_prefori_d1_d2);
set(l, 'LineStyle', '-', 'Color', 'r', 'LineWidth',1);
m=cdfplot(red_d_score_prefori_d1_d3);
set(m, 'LineStyle', '--', 'Color', 'r', 'LineWidth',1);
legend(['Session 1 v Session 2 (Non-KRAB)'], ['Session 1 v Session 3 (Non-KRAB)'], ['Session 1 v Session 2 (KRAB)'], ['Session 1 v Session 3 (KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title([])
hold off

if str2double(expt(d1).img_day) == 1
    print(fullfile(newfnout, [mouse,  '_prefori_change.pdf']), '-dpdf', '-bestfit')
else
   print(fullfile(newfnout, [mouse,  '_prefori_change_2.pdf']), '-dpdf', '-bestfit') 
end
%%
green_d1_k = d1_k_max.k1(1,green_match_all);
green_d2_k = d2_k_max.k1(1,green_match_all);
green_d3_k = d3_k_max.k1(1,green_match_all);

red_d1_k = d1_k_max.k1(1,red_match_all);
red_d2_k = d2_k_max.k1(1,red_match_all);
red_d3_k = d3_k_max.k1(1,red_match_all);

fig = figure;
sgtitle(['k Values (Non-KRAB) ', mouse, ' ', group])
% a=subplot(2,2,1)
h=cdfplot(green_d1_k);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(green_d2_k);
set(j, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',1);
l=cdfplot(green_d3_k);
set(l, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',1);
legend(['Session 1'], ['Session 2'], ['Session 3'], 'Location', 'southeast')
xlim([0 30])
xlabel(['k Value'])
ylabel(['% of Cells'])
title([])
hold off


if str2double(expt(d1).img_day) == 1
    print(fullfile(newfnout, [mouse,  '_green_k_vals.pdf']), '-dpdf', '-bestfit')
else
   print(fullfile(newfnout, [mouse,  '_green_k_vals_2.pdf']), '-dpdf', '-bestfit') 
end


fig = figure;
sgtitle(['k Values (KRAB) ', mouse, ' ', group])
% a=subplot(2,2,1)
h=cdfplot(red_d1_k);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(red_d2_k);
set(j, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',1);
l=cdfplot(red_d3_k);
set(l, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',1);
legend(['Session 1'], ['Session 2'], ['Session 3'], 'Location', 'southeast')
xlim([0 30])
xlabel(['k Value'])
ylabel(['% of Cells'])
title([])
hold off


if str2double(expt(d1).img_day) == 1
    print(fullfile(newfnout, [mouse,  '_red_k_vals.pdf']), '-dpdf', '-bestfit')
else
   print(fullfile(newfnout, [mouse,  '_red_k_vals_2.pdf']), '-dpdf', '-bestfit') 
end



%%
green_d1_max = d1_k_max.max_dfof(1,green_match_all);
green_d2_max = d2_k_max.max_dfof(1,green_match_all);
green_d3_max = d3_k_max.max_dfof(1,green_match_all);

red_d1_max = d1_k_max.max_dfof(1,red_match_all);
red_d2_max = d2_k_max.max_dfof(1,red_match_all);
red_d3_max = d3_k_max.max_dfof(1,red_match_all);


fig = figure;
sgtitle(['Max dF/F Values (Non-KRAB) ', mouse, ' ', group])
% a=subplot(2,2,1)
h=cdfplot(green_d1_max);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(green_d2_max);
set(j, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',1);
l=cdfplot(green_d3_max);
set(l, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',1);
legend(['Session 1'], ['Session 2'], ['Session 3'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F'])
ylabel(['% of Cells'])
title([])
hold off


if str2double(expt(d1).img_day) == 1
    print(fullfile(newfnout, [mouse,  '_green_maxdfof_vals.pdf']), '-dpdf', '-bestfit')
else
  print(fullfile(newfnout, [mouse,  '_green_maxdfof_vals_2.pdf']), '-dpdf', '-bestfit') 
end


fig = figure;
sgtitle(['Max dF/F Values (Non-KRAB) ', mouse, ' ', group])
% a=subplot(2,2,1)
h=cdfplot(red_d1_max);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(red_d2_max);
set(j, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',1);
l=cdfplot(red_d3_max);
set(l, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',1);
legend(['Session 1'], ['Session 2'], ['Session 3'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F'])
ylabel(['% of Cells'])
title([])
hold off


if str2double(expt(d1).img_day) == 1
    print(fullfile(newfnout, [mouse,  '_red_maxdfof_vals.pdf']), '-dpdf', '-bestfit')
else
    print(fullfile(newfnout, [mouse,  '_red_maxdfof_vals_2.pdf']), '-dpdf', '-bestfit') 
end



%%
if str2double(expt(d1).img_day) == 1
%save all of the change scores from above for future use
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_d_scores.mat']), 'green_d_score_prefori_d1_d2', 'green_d_score_prefori_d1_d3', ...
        'red_d_score_prefori_d1_d2', 'red_d_score_prefori_d1_d3')

%save pref ori info
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_ori_changes.mat']), 'green_d1_prefori', 'green_d2_prefori', 'green_d3_prefori', ...
        'red_d1_prefori', 'red_d2_prefori', 'red_d3_prefori')

%save tuned and matched cell IDs
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_id_matches.mat']), 'green_match_all', 'red_match_all')

%save k and max vals
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_k_max_changes.mat']), 'green_d1_k', 'green_d2_k', 'green_d3_k', ...
        'green_d1_max', 'green_d2_max', 'green_d3_max', 'red_d1_k', 'red_d2_k', 'red_d3_k', ...
        'red_d1_max', 'red_d2_max', 'red_d3_max')

%save fit reliabilities
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_fit_reliability.mat']), 'green_d1_reliability', 'green_d2_reliability', 'green_d3_reliability', ...
        'red_d1_reliability', 'red_d2_reliability', 'red_d3_reliability')

else
%save all of the change scores from above for future use
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_d_scores_2.mat']), 'green_d_score_prefori_d1_d2', 'green_d_score_prefori_d1_d3', ...
        'red_d_score_prefori_d1_d2', 'red_d_score_prefori_d1_d3')

%save pref ori info
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_ori_changes_2.mat']), 'green_d1_prefori', 'green_d2_prefori', 'green_d3_prefori', ...
        'red_d1_prefori', 'red_d2_prefori', 'red_d3_prefori')

%save tuned and matched cell IDs
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_id_matches_2.mat']), 'green_match_all', 'red_match_all')

%save k and max vals
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_k_max_changes_2.mat']), 'green_d1_k', 'green_d2_k', 'green_d3_k', ...
        'green_d1_max', 'green_d2_max', 'green_d3_max', 'red_d1_k', 'red_d2_k', 'red_d3_k', ...
        'red_d1_max', 'red_d2_max', 'red_d3_max')

%save fit reliabilities
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_fit_reliability_2.mat']), 'green_d1_reliability', 'green_d2_reliability', 'green_d3_reliability', ...
        'red_d1_reliability', 'red_d2_reliability', 'red_d3_reliability')

end
%%
green_max_all = [green_d1_max; green_d2_max; green_d3_max];
green_grand_max = max(green_max_all);
green_grand_med = median(green_grand_max);

green_bot_50_d_score_prefori_d1_d2 = green_d_score_prefori_d1_d2(green_grand_max < green_grand_med);
green_top_50_d_score_prefori_d1_d2 = green_d_score_prefori_d1_d2(green_grand_max > green_grand_med);

green_bot_50_d_score_prefori_d1_d3 = green_d_score_prefori_d1_d3(green_grand_max < green_grand_med);
green_top_50_d_score_prefori_d1_d3 = green_d_score_prefori_d1_d3(green_grand_max > green_grand_med);

red_max_all = [red_d1_max; red_d2_max; red_d3_max];
red_grand_max = max(red_max_all);
red_grand_med = median(red_grand_max);

red_bot_50_d_score_prefori_d1_d2 = red_d_score_prefori_d1_d2(red_grand_max < red_grand_med);
red_top_50_d_score_prefori_d1_d2 = red_d_score_prefori_d1_d2(red_grand_max > red_grand_med);

red_bot_50_d_score_prefori_d1_d3 = red_d_score_prefori_d1_d3(red_grand_max < red_grand_med);
red_top_50_d_score_prefori_d1_d3 = red_d_score_prefori_d1_d3(red_grand_max > red_grand_med);


if str2double(expt(d1).img_day) == 1

    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_median_sep_d_scores.mat']), 'green_bot_50_d_score_prefori_d1_d2', 'green_top_50_d_score_prefori_d1_d2', ...
        'green_bot_50_d_score_prefori_d1_d3', 'green_top_50_d_score_prefori_d1_d3', 'red_bot_50_d_score_prefori_d1_d2', 'red_top_50_d_score_prefori_d1_d2', ...
        'red_bot_50_d_score_prefori_d1_d3', 'red_top_50_d_score_prefori_d1_d3')
else
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_median_sep_d_scores_2.mat']), 'green_bot_50_d_score_prefori_d1_d2', 'green_top_50_d_score_prefori_d1_d2', ...
        'green_bot_50_d_score_prefori_d1_d3', 'green_top_50_d_score_prefori_d1_d3', 'red_bot_50_d_score_prefori_d1_d2', 'red_top_50_d_score_prefori_d1_d2', ...
        'red_bot_50_d_score_prefori_d1_d3', 'red_top_50_d_score_prefori_d1_d3')
end

%%
green_k_all = [green_d1_k; green_d2_k; green_d3_k];
green_grand_k = max(green_k_all);
green_grand_med_k = median(green_grand_k);

green_wide_d_score_prefori_d1_d2 = green_d_score_prefori_d1_d2(green_grand_k < green_grand_med_k);
green_sharp_d_score_prefori_d1_d2 = green_d_score_prefori_d1_d2(green_grand_k > green_grand_med_k);

green_wide_d_score_prefori_d1_d3 = green_d_score_prefori_d1_d3(green_grand_k < green_grand_med_k);
green_sharp_d_score_prefori_d1_d3 = green_d_score_prefori_d1_d3(green_grand_k > green_grand_med_k);

red_k_all = [red_d1_k; red_d2_k; red_d3_k];
red_grand_k = max(red_k_all);
red_grand_med_k = median(red_grand_k);

red_wide_d_score_prefori_d1_d2 = red_d_score_prefori_d1_d2(red_grand_k < red_grand_med_k);
red_sharp_d_score_prefori_d1_d2 = red_d_score_prefori_d1_d2(red_grand_k > red_grand_med_k);

red_wide_d_score_prefori_d1_d3 = red_d_score_prefori_d1_d3(red_grand_k < red_grand_med_k);
red_sharp_d_score_prefori_d1_d3 = red_d_score_prefori_d1_d3(red_grand_k > red_grand_med_k);


if str2double(expt(d1).img_day) == 1

    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_wide_sharp_sep_d_scores.mat']), 'green_wide_d_score_prefori_d1_d2', 'green_sharp_d_score_prefori_d1_d2', ...
        'green_wide_d_score_prefori_d1_d3', 'green_sharp_d_score_prefori_d1_d3', 'red_wide_d_score_prefori_d1_d2', 'red_sharp_d_score_prefori_d1_d2', ...
        'red_wide_d_score_prefori_d1_d3', 'red_sharp_d_score_prefori_d1_d3')
else
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_wide_sharp_sep_d_scores_2.mat']), 'green_wide_d_score_prefori_d1_d2', 'green_sharp_d_score_prefori_d1_d2', ...
        'green_wide_d_score_prefori_d1_d3', 'green_sharp_d_score_prefori_d1_d3', 'red_wide_d_score_prefori_d1_d2', 'red_sharp_d_score_prefori_d1_d2', ...
        'red_wide_d_score_prefori_d1_d3', 'red_sharp_d_score_prefori_d1_d3')
end


%% do the median split analyses but with reliability as the dependent variable*** finish this and run it with lacz mice

green_bot_50_reliability_d1 = green_d1_reliability(green_grand_max < green_grand_med);
green_top_50_reliability_d1 = green_d1_reliability(green_grand_max > green_grand_med);
green_bot_50_reliability_d2 = green_d2_reliability(green_grand_max < green_grand_med);
green_top_50_reliability_d2 = green_d2_reliability(green_grand_max > green_grand_med);
green_bot_50_reliability_d3 = green_d3_reliability(green_grand_max < green_grand_med);
green_top_50_reliability_d3 = green_d3_reliability(green_grand_max > green_grand_med);

red_bot_50_reliability_d1 = red_d1_reliability(red_grand_max < red_grand_med);
red_top_50_reliability_d1 = red_d1_reliability(red_grand_max > red_grand_med);
red_bot_50_reliability_d2 = red_d2_reliability(red_grand_max < red_grand_med);
red_top_50_reliability_d2 = red_d2_reliability(red_grand_max > red_grand_med);
red_bot_50_reliability_d3 = red_d3_reliability(red_grand_max < red_grand_med);
red_top_50_reliability_d3 = red_d3_reliability(red_grand_max > red_grand_med);


green_wide_reliability_d1 = green_d1_reliability(green_grand_k < green_grand_med_k);
green_sharp_reliability_d1 = green_d1_reliability(green_grand_k > green_grand_med_k);
green_wide_reliability_d2 = green_d2_reliability(green_grand_k < green_grand_med_k);
green_sharp_reliability_d2 = green_d2_reliability(green_grand_k > green_grand_med_k);
green_wide_reliability_d3 = green_d3_reliability(green_grand_k < green_grand_med_k);
green_sharp_reliability_d3 = green_d3_reliability(green_grand_k > green_grand_med_k);

red_wide_reliability_d1 = red_d1_reliability(red_grand_k < red_grand_med_k);
red_sharp_reliability_d1 = red_d1_reliability(red_grand_k > red_grand_med_k);
red_wide_reliability_d2 = red_d2_reliability(red_grand_k < red_grand_med_k);
red_sharp_reliability_d2 = red_d2_reliability(red_grand_k > red_grand_med_k);
red_wide_reliability_d3 = red_d3_reliability(red_grand_k < red_grand_med_k);
red_sharp_reliability_d3 = red_d3_reliability(red_grand_k > red_grand_med_k);


if str2double(expt(d1).img_day) == 1

    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_max_median_sep_reliability.mat']), 'green_bot_50_reliability_d1', 'green_bot_50_reliability_d2', ...
        'green_bot_50_reliability_d3', 'green_top_50_reliability_d1', 'green_top_50_reliability_d2', 'green_top_50_reliability_d3', ...
        'red_bot_50_reliability_d1', 'red_bot_50_reliability_d2', 'red_bot_50_reliability_d3', 'red_top_50_reliability_d1', 'red_top_50_reliability_d2', 'red_top_50_reliability_d3')
else
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_max_median_sep_reliability_2.mat']), 'green_bot_50_reliability_d1', 'green_bot_50_reliability_d2', ...
        'green_bot_50_reliability_d3', 'green_top_50_reliability_d1', 'green_top_50_reliability_d2', 'green_top_50_reliability_d3', ...
        'red_bot_50_reliability_d1', 'red_bot_50_reliability_d2', 'red_bot_50_reliability_d3', 'red_top_50_reliability_d1', 'red_top_50_reliability_d2', 'red_top_50_reliability_d3')
end


if str2double(expt(d1).img_day) == 1

    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_k_median_sep_reliability.mat']), 'green_wide_reliability_d1', 'green_wide_reliability_d2', ...
        'green_wide_reliability_d3', 'green_sharp_reliability_d1', 'green_sharp_reliability_d2', 'green_sharp_reliability_d3', ...
        'red_wide_reliability_d1', 'red_wide_reliability_d2', 'red_wide_reliability_d3', 'red_sharp_reliability_d1', 'red_sharp_reliability_d2', 'red_sharp_reliability_d3')
else
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_k_median_sep_reliability_2.mat']), 'green_wide_reliability_d1', 'green_wide_reliability_d2', ...
        'green_wide_reliability_d3', 'green_sharp_reliability_d1', 'green_sharp_reliability_d2', 'green_sharp_reliability_d3', ...
        'red_wide_reliability_d1', 'red_wide_reliability_d2', 'red_wide_reliability_d3', 'red_sharp_reliability_d1', 'red_sharp_reliability_d2', 'red_sharp_reliability_d3')
end


%%
%excluding cells pinned to 30 k value - need to do similar to above but exclude BEFORE doing median
%split - still figuring this out

green_not_too_sharp = find(green_grand_k <= 29);
green_too_sharp = find(green_grand_k > 29);
red_not_too_sharp = find(red_grand_k <= 29);
red_too_sharp = find(red_grand_k > 29);

no30k_green_k_all = green_k_all(:,green_not_too_sharp);
no30k_red_k_all = red_k_all(:,red_not_too_sharp);

no30k_green_d_score_prefori_d1_d2 = green_d_score_prefori_d1_d2(green_not_too_sharp);
no30k_green_d_score_prefori_d1_d3 = green_d_score_prefori_d1_d3(green_not_too_sharp);
no30k_red_d_score_prefori_d1_d2 = red_d_score_prefori_d1_d2(red_not_too_sharp);
no30k_red_d_score_prefori_d1_d3 = red_d_score_prefori_d1_d3(red_not_too_sharp);


no30k_green_grand_k = max(no30k_green_k_all);
no30k_green_grand_med_k = median(no30k_green_grand_k);
no30k_red_grand_k = max(no30k_red_k_all);
no30k_red_grand_med_k = median(no30k_red_grand_k);

no30k_green_wide_d_score_prefori_d1_d2 = no30k_green_d_score_prefori_d1_d2(no30k_green_grand_k < no30k_green_grand_med_k);
no30k_green_sharp_d_score_prefori_d1_d2 = no30k_green_d_score_prefori_d1_d2(no30k_green_grand_k > no30k_green_grand_med_k);
no30k_red_wide_d_score_prefori_d1_d2 = no30k_red_d_score_prefori_d1_d2(no30k_red_grand_k < no30k_red_grand_med_k);
no30k_red_sharp_d_score_prefori_d1_d2 = no30k_red_d_score_prefori_d1_d2(no30k_red_grand_k > no30k_red_grand_med_k);
no30k_green_wide_d_score_prefori_d1_d3 = no30k_green_d_score_prefori_d1_d3(no30k_green_grand_k < no30k_green_grand_med_k);
no30k_green_sharp_d_score_prefori_d1_d3 = no30k_green_d_score_prefori_d1_d3(no30k_green_grand_k > no30k_green_grand_med_k);
no30k_red_wide_d_score_prefori_d1_d3 = no30k_red_d_score_prefori_d1_d3(no30k_red_grand_k < no30k_red_grand_med_k);
no30k_red_sharp_d_score_prefori_d1_d3 = no30k_red_d_score_prefori_d1_d3(no30k_red_grand_k > no30k_red_grand_med_k);


if str2double(expt(d1).img_day) == 1

    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_no30k_k_median_sep_prefori.mat']), 'no30k_green_wide_d_score_prefori_d1_d2', 'no30k_green_sharp_d_score_prefori_d1_d2', ...
        'no30k_red_wide_d_score_prefori_d1_d2', 'no30k_red_sharp_d_score_prefori_d1_d2', 'no30k_green_wide_d_score_prefori_d1_d3', 'no30k_green_sharp_d_score_prefori_d1_d3', ...
        'no30k_red_wide_d_score_prefori_d1_d3', 'no30k_red_sharp_d_score_prefori_d1_d3')
else
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_no30k_k_median_sep_prefori_2.mat']), 'no30k_green_wide_d_score_prefori_d1_d2', 'no30k_green_sharp_d_score_prefori_d1_d2', ...
        'no30k_red_wide_d_score_prefori_d1_d2', 'no30k_red_sharp_d_score_prefori_d1_d2', 'no30k_green_wide_d_score_prefori_d1_d3', 'no30k_green_sharp_d_score_prefori_d1_d3', ...
        'no30k_red_wide_d_score_prefori_d1_d3', 'no30k_red_sharp_d_score_prefori_d1_d3')
end


%%
%weakest/sharpest overlap

green_sharpest_ind = find(green_grand_k > green_grand_med_k);
green_weakest_ind = find(green_grand_max < green_grand_med);
green_sharpest_and_weakest_ind = intersect(green_sharpest_ind, green_weakest_ind);

n_green_sharpest_and_weakest = length(green_sharpest_and_weakest_ind);
n_green_total = length(green_grand_k);
overlap_green_sharp_weak = n_green_sharpest_and_weakest/n_green_total;

red_sharpest_ind = find(red_grand_k > red_grand_med_k);
red_weakest_ind = find(red_grand_max < red_grand_med);
red_sharpest_and_weakest_ind = intersect(red_sharpest_ind, red_weakest_ind);

n_red_sharpest_and_weakest = length(red_sharpest_and_weakest_ind);
n_red_total = length(red_grand_k);
overlap_red_sharp_weak = n_red_sharpest_and_weakest/n_red_total;

if str2double(expt(d1).img_day) == 1

    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, 'weak_sharp_overlap.mat']), 'n_green_sharpest_and_weakest', 'n_green_total', ...
        'overlap_green_sharp_weak', 'n_red_sharpest_and_weakest', 'n_red_total', 'overlap_red_sharp_weak')
else
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, 'weak_sharp_overlap_2.mat']), 'n_green_sharpest_and_weakest', 'n_green_total', ...
        'overlap_green_sharp_weak', 'n_red_sharpest_and_weakest', 'n_red_total', 'overlap_red_sharp_weak')
end


%%
