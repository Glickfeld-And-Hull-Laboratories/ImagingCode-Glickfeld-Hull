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
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_virusTG'; %folder to save files to
dataset = 'exp_list_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%

%we start with loading the proper information based on the mouse's info in
%the experiment list

d1 = 51; 
d2 = 53; 
d3 = 55;

mouse = expt(d1).mouse; 
group = expt(d1).indicator{1};
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

%load multiday data for day 2
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));
d3_matches = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));
d3_k_max = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'k_and_max_vals.mat']));




%% matching indices

tuned_d1 = d1_ori.ind_theta90;
tuned_d2 = d2_ori.ind_theta90;
tuned_d3 = d3_ori.ind_theta90;

tuned_all = unique([tuned_d1 tuned_d2 tuned_d3]);

%find green cells that match from d1 to d2 and d1 to d3 
green_match_d2 = find([d2_matches.cellImageAlign.pass]); 
green_match_d3 = find([d3_matches.cellImageAlign.pass]); 

green_match_all = intersect(intersect(green_match_d2,green_match_d3),tuned_all);

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



fig = figure;
sgtitle(['Pref Ori Changes ', mouse, ' ', group])
% a=subplot(2,2,1)
h=cdfplot(green_d_score_prefori_d1_d2);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',1);
hold on
j=cdfplot(green_d_score_prefori_d1_d3);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',1);
legend(['Session 1 v Session 2'], ['Session 1 v Session 3'], 'Location', 'southeast')
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



%%
green_d1_max = d1_k_max.max_dfof(1,green_match_all);
green_d2_max = d2_k_max.max_dfof(1,green_match_all);
green_d3_max = d3_k_max.max_dfof(1,green_match_all);


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




%%
if str2double(expt(d1).img_day) == 1
%save all of the change scores from above for future use
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_d_scores.mat']), 'green_d_score_prefori_d1_d2', 'green_d_score_prefori_d1_d3')

%save pref ori info
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_ori_changes.mat']), 'green_d1_prefori', 'green_d2_prefori', 'green_d3_prefori')

%save tuned and matched cell IDs
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_id_matches.mat']), 'green_match_all')

%save k and max vals
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_k_max_changes.mat']), 'green_d1_k', 'green_d2_k', 'green_d3_k', ...
        'green_d1_max', 'green_d2_max', 'green_d3_max')

else
%save all of the change scores from above for future use
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_d_scores_2.mat']), 'green_d_score_prefori_d1_d2', 'green_d_score_prefori_d1_d3')

%save pref ori info
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_ori_changes_2.mat']), 'green_d1_prefori', 'green_d2_prefori', 'green_d3_prefori')

%save tuned and matched cell IDs
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_id_matches_2.mat']), 'green_match_all')

%save k and max vals
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_k_max_changes_2.mat']), 'green_d1_k', 'green_d2_k', 'green_d3_k', ...
        'green_d1_max', 'green_d2_max', 'green_d3_max')

end
%%
green_max_all = [green_d1_max; green_d2_max; green_d3_max];
green_grand_max = max(green_max_all);
green_grand_med = median(green_grand_max);

green_bot_50_d_score_prefori_d1_d2 = green_d_score_prefori_d1_d2(green_grand_max < green_grand_med);
green_top_50_d_score_prefori_d1_d2 = green_d_score_prefori_d1_d2(green_grand_max > green_grand_med);

green_bot_50_d_score_prefori_d1_d3 = green_d_score_prefori_d1_d3(green_grand_max < green_grand_med);
green_top_50_d_score_prefori_d1_d3 = green_d_score_prefori_d1_d3(green_grand_max > green_grand_med);



if str2double(expt(d1).img_day) == 1

    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_median_sep_d_scores.mat']), 'green_bot_50_d_score_prefori_d1_d2', 'green_top_50_d_score_prefori_d1_d2', ...
        'green_bot_50_d_score_prefori_d1_d3', 'green_top_50_d_score_prefori_d1_d3')
else
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_median_sep_d_scores_2.mat']), 'green_bot_50_d_score_prefori_d1_d2', 'green_top_50_d_score_prefori_d1_d2', ...
        'green_bot_50_d_score_prefori_d1_d3', 'green_top_50_d_score_prefori_d1_d3')
end