%clear everything
clear all
clear all global
clc
%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P\Arc_greenVred'; %folder to save files to
dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
d1 = 46; %day 1 in expt list
d2 = 47; %day 2 in expt list
d3 = 48; %day 3 in expt list
mouse = expt(d1).mouse; %mouse
ref_str = 'runs-003'; %string on file name to load
ref_str_d1 = ['runs-',expt(d1).runs]; %need to fix this part***
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date_d1 = expt(d1).date; %day 1 (ref day) date
date_d2 = expt(d2).date; %day 2 date
date_d3 = expt(d3).date; %day 3 date
img_folder_d1 = expt(d1).runs; %img folder of day 1
img_folder_d2 = expt(d2).runs; %img folder of day 2
img_folder_d3 = expt(d3).runs; %img folder of day 3
%%
%load relevant data (tunings and multi day)

%load ori info for each day
d1_ori = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningInfo.mat']));
d2_ori = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str], [date_d2 '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
d3_ori = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str], [date_d3 '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str], [date_d2 '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
d3_k_max = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str], [date_d3 '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));

%load multiday data for days 2 and 3
d1_matches = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str], [date_d1 '_' mouse '_' ref_str '_' 'multiday_alignment.mat']));
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str], [date_d2 '_' mouse '_' ref_str '_' 'multiday_alignment.mat']));
d3_matches = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str], [date_d3 '_' mouse '_' ref_str '_' 'multiday_alignment.mat']));
d2_tcs = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str], [date_d2 '_' mouse '_' ref_str '_' 'TCs.mat']));
d3_tcs = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str], [date_d3 '_' mouse '_' ref_str '_' 'TCs.mat']));

%green and red cell data
d2_matches_tj = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str], [date_d2 '_' mouse '_' ref_str '_' 'tj_matches.mat']));
d3_matches_tj = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str], [date_d3 '_' mouse '_' ref_str '_' 'tj_matches.mat']));



%%
%separate cells that are tuned and matched from those that are not tuned
%but still matched - plot in separate colors and on separate figures

%find tuned cells for each day and cells that are both tuned and matched
%across days
tuned_d1 = d1_ori.ind_theta90;
tuned_d2 = d2_ori.ind_theta90;
tuned_d3 = d3_ori.ind_theta90;


%find green cells that match from d1 to d2 and d1 to d3 
green_match_d2 = d2_matches_tj.green_match_ind; 
green_match_d3 = d3_matches_tj.green_match_ind; 

%find red cells that match from d1 to d2 and d1 to d3 
red_match_d2 = d2_matches_tj.red_match_ind; 
red_match_d3 = d3_matches_tj.red_match_ind; 


green_lax_tuned_matched_d2 = find(ismember(green_match_d2, tuned_d1) | ismember(green_match_d2, tuned_d2));
green_lax_tuned_matched_d3 = find(ismember(green_match_d3, tuned_d1) | ismember(green_match_d3, tuned_d3));
red_lax_tuned_matched_d2 = find(ismember(red_match_d2, tuned_d1) | ismember(red_match_d2, tuned_d2));
red_lax_tuned_matched_d3 = find(ismember(red_match_d3, tuned_d1) | ismember(red_match_d3, tuned_d3));


green_matched_d1_tuned_d2 = intersect(green_match_d2, tuned_d1);
green_matched_d1_d2_tuned = intersect(green_match_d2, tuned_d2);
green_matched_d1_d2_all_tuned = unique([green_matched_d1_tuned_d2 green_matched_d1_tuned_d2]);
green_matched_d1_tuned_d3 = intersect(green_match_d3, tuned_d1);
green_matched_d1_d3_tuned = intersect(green_match_d3, tuned_d3);
green_matched_d1_d3_all_tuned = unique([green_matched_d1_tuned_d3 green_matched_d1_tuned_d3]);

red_matched_d1_tuned_d2 = intersect(red_match_d2, tuned_d1);
red_matched_d1_d2_tuned = intersect(red_match_d2, tuned_d2);
red_matched_d1_d2_all_tuned = unique([red_matched_d1_tuned_d2 red_matched_d1_tuned_d2]);
red_matched_d1_tuned_d3 = intersect(red_match_d3, tuned_d1);
red_matched_d1_d3_tuned = intersect(red_match_d3, tuned_d3);
red_matched_d1_d3_all_tuned = unique([red_matched_d1_tuned_d3 red_matched_d1_tuned_d3]);

%%
%identify pref oris for each day


%find pref oris of green matched cells for d1 and d2
green_prefori_d1_d2_match = d1_ori.prefOri(1,green_matched_d1_d2_all_tuned );
green_prefori_d2_match = d2_ori.prefOri(1,green_matched_d1_d2_all_tuned);

%find pref oris of red matched cells for d1 and d2
red_prefori_d1_d2_match = d1_ori.prefOri(1,red_matched_d1_d2_all_tuned);
red_prefori_d2_match = d2_ori.prefOri(1,red_matched_d1_d2_all_tuned);

%find pref oris of green matched cells for d1 and d3
green_prefori_d1_d3_match = d1_ori.prefOri(1,green_matched_d1_d3_all_tuned);
green_prefori_d3_match = d3_ori.prefOri(1,green_matched_d1_d3_all_tuned);

%find pref oris of red matched cells for d1 and d3
red_prefori_d1_d3_match = d1_ori.prefOri(1,red_matched_d1_d3_all_tuned);
red_prefori_d3_match = d3_ori.prefOri(1,red_matched_d1_d3_all_tuned);

green_tuned_d1 = intersect(d1_matches.greenCells, tuned_d1);
red_tuned_d1 = intersect(d1_matches.redCells, tuned_d1);



%%

figure('Position', [400 10 650 800]);
sgtitle([mouse, ' ', img_area, ' ', img_layer], 'Interpreter', 'None');
subplot(2,1,1);
scatter(green_prefori_d1_d2_match, green_prefori_d2_match,'g');
hold on
scatter(red_prefori_d1_d2_match, red_prefori_d2_match,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Green (n = ', num2str(length(green_prefori_d1_d2_match)), ')'], ['Red (n = ', num2str(length(red_prefori_d1_d2_match)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(green_prefori_d1_d3_match, green_prefori_d3_match,'g');
hold on
scatter(red_prefori_d1_d3_match, red_prefori_d3_match,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 3 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
legend(['Green (n = ', num2str(length(green_prefori_d1_d3_match)), ')'], ['Red (n = ', num2str(length(red_prefori_d1_d3_match)), ')'], 'Location', 'northwest')

if str2double(expt(d1).img_day) == 1
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_change_ori_scatter.pdf']), '-dpdf', '-bestfit')
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_ori_changes.mat']), 'green_prefori_d1_d2_match', 'green_prefori_d1_d3_match', 'green_prefori_d2_match', 'green_prefori_d3_match', ...
        'red_prefori_d1_d2_match', 'red_prefori_d1_d3_match', 'red_prefori_d2_match', 'red_prefori_d3_match', "green_matched_d1_d2_all_tuned", 'green_matched_d1_d3_all_tuned', ...
        'red_matched_d1_d2_all_tuned', 'red_matched_d1_d3_all_tuned', 'tuned_d1', 'red_tuned_d1', 'green_tuned_d1')
else
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_change_ori_scatter_2.pdf']), '-dpdf', '-bestfit')
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_ori_changes_2.mat']), 'green_prefori_d1_d2_match', 'green_prefori_d1_d3_match', 'green_prefori_d2_match', 'green_prefori_d3_match', ...
        'red_prefori_d1_d2_match', 'red_prefori_d1_d3_match', 'red_prefori_d2_match', 'red_prefori_d3_match', "green_matched_d1_d2_all_tuned", 'green_matched_d1_d3_all_tuned', ...
        'red_matched_d1_d2_all_tuned', 'red_matched_d1_d3_all_tuned', 'tuned_d1', 'red_tuned_d1', 'green_tuned_d1')
end

%%
%diff scores


%d1_2
green_dscore_prefori_d1_d2_match = double(green_prefori_d1_d2_match>90);
green_dscore_prefori_d1_d2_match(green_dscore_prefori_d1_d2_match>0) = 180;
green_dscore_prefori_d1_d2_match = abs(green_dscore_prefori_d1_d2_match-green_prefori_d1_d2_match);

red_dscore_prefori_d1_d2_match = double(red_prefori_d1_d2_match>90);
red_dscore_prefori_d1_d2_match(red_dscore_prefori_d1_d2_match>0) = 180;
red_dscore_prefori_d1_d2_match = abs(red_dscore_prefori_d1_d2_match-red_prefori_d1_d2_match);

%d2
green_dscore_prefori_d2_match = double(green_prefori_d2_match>90);
green_dscore_prefori_d2_match(green_dscore_prefori_d2_match>0) = 180;
green_dscore_prefori_d2_match = abs(green_dscore_prefori_d2_match-green_prefori_d2_match);

red_dscore_prefori_d2_match = double(red_prefori_d2_match>90);
red_dscore_prefori_d2_match(red_dscore_prefori_d2_match>0) = 180;
red_dscore_prefori_d2_match = abs(red_dscore_prefori_d2_match-red_prefori_d2_match);

%d1_d3
green_dscore_prefori_d1_d3_match = double(green_prefori_d1_d3_match>90);
green_dscore_prefori_d1_d3_match(green_dscore_prefori_d1_d3_match>0) = 180;
green_dscore_prefori_d1_d3_match = abs(green_dscore_prefori_d1_d3_match-green_prefori_d1_d3_match);

red_dscore_prefori_d1_d3_match = double(red_prefori_d1_d3_match>90);
red_dscore_prefori_d1_d3_match(red_dscore_prefori_d1_d3_match>0) = 180;
red_dscore_prefori_d1_d3_match = abs(red_dscore_prefori_d1_d3_match-red_prefori_d1_d3_match);

%d3
green_dscore_prefori_d3_match = double(green_prefori_d3_match>90);
green_dscore_prefori_d3_match(green_dscore_prefori_d3_match>0) = 180;
green_dscore_prefori_d3_match = abs(green_dscore_prefori_d3_match-green_prefori_d3_match);

red_dscore_prefori_d3_match = double(red_prefori_d3_match>90);
red_dscore_prefori_d3_match(red_dscore_prefori_d3_match>0) = 180;
red_dscore_prefori_d3_match = abs(red_dscore_prefori_d3_match-red_prefori_d3_match);

green_d_score_prefori_d1_d2 = abs(green_dscore_prefori_d1_d2_match-green_dscore_prefori_d2_match);
green_d_score_prefori_d1_d3 = abs(green_dscore_prefori_d1_d3_match-green_dscore_prefori_d3_match);

red_d_score_prefori_d1_d2 = abs(red_dscore_prefori_d1_d2_match-red_dscore_prefori_d2_match);
red_d_score_prefori_d1_d3 = abs(red_dscore_prefori_d1_d3_match-red_dscore_prefori_d3_match);

green_d_scores_prefori_all = [green_d_score_prefori_d1_d2, green_d_score_prefori_d1_d3];
red_d_scores_prefori_all = [red_d_score_prefori_d1_d2, red_d_score_prefori_d1_d3];

figure; 
h=cdfplot(green_d_score_prefori_d1_d2);
hold on
j=cdfplot(green_d_score_prefori_d1_d3);
k=cdfplot(red_d_score_prefori_d1_d2);
l=cdfplot(red_d_score_prefori_d1_d3);
set(h, 'LineStyle', '-', 'Color', 'g');
set(j, 'LineStyle', '--', 'Color', 'g');
set(k, 'LineStyle', '-', 'Color', 'r');
set(l, 'LineStyle', '--', 'Color', 'r');
hold off
legend(['Day 1 vs. Day 2 Green (n = ', num2str(length(green_d_score_prefori_d1_d2)), ')'], ['Day 1 vs. Day 3 Green (n = ', num2str(length(green_d_score_prefori_d1_d3)), ')'], ...
    ['Day 1 vs. Day 2 Red (n = ', num2str(length(red_d_score_prefori_d1_d2)), ')'], ['Day 1 vs. Day 3 Red (n = ', num2str(length(red_d_score_prefori_d1_d3)), ')'],'Location', 'best')
% legend(['Day 1 v Day 2'; 'Day 1 v Day 3'])
xlabel('Change in Pref Ori')
xlim([0 90]);
title([mouse, ' ', img_area, ' ', img_layer], 'Interpreter', 'None');

if str2double(expt(d1).img_day) == 1
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW_change_ori_cdf.pdf']), '-dpdf', '-bestfit')
else
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW_change_ori_cdf_2.pdf']), '-dpdf', '-bestfit')

end

%% K VALUES

%find pref oris of matched cells for d1 and d2
green_k_d1_d2_match = d1_k_max.k1(1,green_matched_d1_d2_all_tuned);
green_k_d2_match = d2_k_max.k1(1,green_matched_d1_d2_all_tuned);

red_k_d1_d2_match = d1_k_max.k1(1,red_matched_d1_d2_all_tuned);
red_k_d2_match = d2_k_max.k1(1,red_matched_d1_d2_all_tuned);

%find pref oris of matched cells for d1 and d2
green_k_d1_d3_match = d1_k_max.k1(1,green_matched_d1_d3_all_tuned);
green_k_d3_match = d3_k_max.k1(1,green_matched_d1_d3_all_tuned);

red_k_d1_d3_match = d1_k_max.k1(1,red_matched_d1_d3_all_tuned);
red_k_d3_match = d3_k_max.k1(1,red_matched_d1_d3_all_tuned);

figure('Position', [400 10 650 800]);
sgtitle([mouse, ' ', img_area, ' ', img_layer], 'Interpreter', 'None');
subplot(2,1,1);
scatter(green_k_d1_d2_match,green_k_d2_match, 'g');
hold on
scatter(red_k_d1_d2_match,red_k_d2_match,'r');
hold off
xlabel('Day 1 k Value');
ylabel('Day 2 k Value');
xlim([0,30]);
ylim([0,30]);
%xticks(0:20:180);
%yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Green (n = ', num2str(length(green_k_d1_d2_match)), ')'], ['Red (n = ', num2str(length(red_k_d1_d2_match)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(green_k_d1_d3_match,green_k_d3_match, 'g');
hold on
scatter(red_k_d1_d3_match,red_k_d3_match,'r');
hold off
xlabel('Day 1 k Value');
ylabel('Day 3 k Value');
xlim([0,30]);
ylim([0,30]);
%xticks(0:20:180);
%yticks(0:20:180);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
legend(['Green (n = ', num2str(length(green_k_d1_d3_match)), ')'], ['Red (n = ', num2str(length(red_k_d1_d3_match)), ')'], 'Location', 'northwest')

if str2double(expt(d1).img_day) == 1
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW_change_k_scatter.pdf']), '-dpdf', '-bestfit')
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_k_changes.mat']), 'green_k_d1_d2_match', 'green_k_d1_d3_match', 'green_k_d2_match', 'green_k_d3_match', ...
        'red_k_d1_d2_match', 'red_k_d1_d3_match', 'red_k_d2_match', 'red_k_d3_match')
else
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW_change_k_scatter_2.pdf']), '-dpdf', '-bestfit')
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_k_changes_2.mat']), 'green_k_d1_d2_match', 'green_k_d1_d3_match', 'green_k_d2_match', 'green_k_d3_match', ...
        'red_k_d1_d2_match', 'red_k_d1_d3_match', 'red_k_d2_match', 'red_k_d3_match')
end

%d1_2
green_dscore_k_d1_d2_match = abs(green_k_d1_d2_match-green_k_d2_match);
red_dscore_k_d1_d2_match = abs(red_k_d1_d2_match-red_k_d2_match);
%d1_d3
green_dscore_k_d1_d3_match = abs(green_k_d1_d3_match-green_k_d3_match);
red_dscore_k_d1_d3_match = abs(red_k_d1_d3_match-red_k_d3_match);

figure; 
h=cdfplot(green_dscore_k_d1_d2_match);
hold on
j=cdfplot(green_dscore_k_d1_d3_match);
k=cdfplot(red_dscore_k_d1_d2_match);
l=cdfplot(red_dscore_k_d1_d3_match);
set(h, 'LineStyle', '-', 'Color', 'g');
set(j, 'LineStyle', '--', 'Color', 'g');
set(k, 'LineStyle', '-', 'Color', 'r');
set(l, 'LineStyle', '--', 'Color', 'r');
hold off
legend(['Day 1 vs. Day 2 Green (n = ', num2str(length(green_dscore_k_d1_d2_match)), ')'], ['Day 1 vs. Day 3 Green (n = ', num2str(length(green_dscore_k_d1_d3_match)), ')'], ...
    ['Day 1 vs. Day 2 Red (n = ', num2str(length(red_dscore_k_d1_d2_match)), ')'], ['Day 1 vs. Day 3 Red (n = ', num2str(length(red_dscore_k_d1_d3_match)), ')'],'Location', 'best')
xlabel('Change in k Value')
%xlim([0 90]);
title([mouse, ' ', img_area, ' ', img_layer], 'Interpreter', 'None');

if str2double(expt(d1).img_day) == 1
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW_change_k_cdf.pdf']), '-dpdf', '-bestfit')
else
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW_change_k_cdf_2.pdf']), '-dpdf', '-bestfit')
end

%% MAX VALUES

%find pref oris of matched cells for d1 and d2
green_max_d1_d2_match = d1_k_max.max_dfof(1,green_matched_d1_d2_all_tuned);
green_max_d2_match = d2_k_max.max_dfof(1,green_matched_d1_d2_all_tuned);

red_max_d1_d2_match = d1_k_max.max_dfof(1,red_matched_d1_d2_all_tuned);
red_max_d2_match = d2_k_max.max_dfof(1,red_matched_d1_d2_all_tuned);

%find pref oris of matched cells for d1 and d2
green_max_d1_d3_match = d1_k_max.max_dfof(1,green_matched_d1_d3_all_tuned);
green_max_d3_match = d3_k_max.max_dfof(1,green_matched_d1_d3_all_tuned);

red_max_d1_d3_match = d1_k_max.max_dfof(1,red_matched_d1_d3_all_tuned);
red_max_d3_match = d3_k_max.max_dfof(1,red_matched_d1_d3_all_tuned);

figure('Position', [400 10 650 800]);
sgtitle([mouse, ' ', img_area, ' ', img_layer], 'Interpreter', 'None');
subplot(2,1,1);
scatter(green_max_d1_d2_match,green_max_d2_match, 'g');
hold on
scatter(red_max_d1_d2_match,red_max_d2_match,'r');
hold off
xlabel('Day 1 Max dF/F Value');
ylabel('Day 2 Max dF/F Value');
xlim([0,1]);
ylim([0,1]);
%xticks(0:20:180);
%yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Green (n = ', num2str(length(green_max_d1_d2_match)), ')'], ['Red (n = ', num2str(length(red_max_d1_d2_match)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(green_max_d1_d3_match,green_max_d3_match, 'g');
hold on
scatter(red_max_d1_d3_match,red_max_d3_match,'r');
hold off
xlabel('Day 1 Max dF/F Value');
ylabel('Day 3 Max dF/F Value');
xlim([0,1]);
ylim([0,1]);
%xticks(0:20:180);
%yticks(0:20:180);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
legend(['Green (n = ', num2str(length(green_max_d1_d3_match)), ')'], ['Red (n = ', num2str(length(red_max_d1_d3_match)), ')'], 'Location', 'northwest')

if str2double(expt(d1).img_day) == 1
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW_change_max_scatter.pdf']), '-dpdf', '-bestfit')
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_max_changes.mat']), 'green_max_d1_d2_match', 'green_max_d1_d3_match', 'green_max_d2_match', 'green_max_d3_match', ...
        'red_max_d1_d2_match', 'red_max_d1_d3_match', 'red_max_d2_match', 'red_max_d3_match')

else
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW_change_max_scatter_2.pdf']), '-dpdf', '-bestfit')
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_max_changes_2.mat']), 'green_max_d1_d2_match', 'green_max_d1_d3_match', 'green_max_d2_match', 'green_max_d3_match', ...
        'red_max_d1_d2_match', 'red_max_d1_d3_match', 'red_max_d2_match', 'red_max_d3_match')

end

%d1_2
green_dscore_max_d1_d2_match = abs(green_max_d1_d2_match-green_max_d2_match);
red_dscore_max_d1_d2_match = abs(red_max_d1_d2_match-red_max_d2_match);

%d1_d3
green_dscore_max_d1_d3_match = abs(green_max_d1_d3_match-green_max_d3_match);
red_dscore_max_d1_d3_match = abs(red_max_d1_d3_match-red_max_d3_match);


figure; 
h=cdfplot(green_dscore_max_d1_d2_match);
hold on
j=cdfplot(green_dscore_max_d1_d3_match);
k=cdfplot(red_dscore_max_d1_d2_match);
l=cdfplot(red_dscore_max_d1_d3_match);
set(h, 'LineStyle', '-', 'Color', 'g');
set(j, 'LineStyle', '--', 'Color', 'g');
set(k, 'LineStyle', '-', 'Color', 'r');
set(l, 'LineStyle', '--', 'Color', 'r');
hold off
legend(['Day 1 vs. Day 2 Green (n = ', num2str(length(green_dscore_max_d1_d2_match)), ')'], ['Day 1 vs. Day 3 Green (n = ', num2str(length(green_dscore_max_d1_d3_match)), ')'], ...
    ['Day 1 vs. Day 2 Red (n = ', num2str(length(red_dscore_max_d1_d2_match)), ')'], ['Day 1 vs. Day 3 Red (n = ', num2str(length(red_dscore_max_d1_d3_match)), ')'],'Location', 'best')
xlabel('Change in Max dF/F Value')
%xlim([0 90]);
title([mouse, ' ', img_area, ' ', img_layer], 'Interpreter', 'None');



if str2double(expt(d1).img_day) == 1
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW_change_max_cdf.pdf']), '-dpdf', '-bestfit')
else
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW_change_max_cdf_2.pdf']), '-dpdf', '-bestfit')
end

%%
if str2double(expt(d1).img_day) == 1
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_d_scores.mat']), 'green_d_score_prefori_d1_d2', 'green_d_score_prefori_d1_d3', ...
        'green_dscore_k_d1_d2_match', 'green_dscore_k_d1_d3_match', 'green_dscore_max_d1_d2_match', 'green_dscore_max_d1_d3_match', ...
        'red_d_score_prefori_d1_d2', 'red_d_score_prefori_d1_d3', ...
        'red_dscore_k_d1_d2_match', 'red_dscore_k_d1_d3_match', 'red_dscore_max_d1_d2_match', 'red_dscore_max_d1_d3_match')
else
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_d_scores_2.mat']), 'green_d_score_prefori_d1_d2', 'green_d_score_prefori_d1_d3', ...
        'green_dscore_k_d1_d2_match', 'green_dscore_k_d1_d3_match', 'green_dscore_max_d1_d2_match', 'green_dscore_max_d1_d3_match', ...
        'red_d_score_prefori_d1_d2', 'red_d_score_prefori_d1_d3', ...
        'red_dscore_k_d1_d2_match', 'red_dscore_k_d1_d3_match', 'red_dscore_max_d1_d2_match', 'red_dscore_max_d1_d3_match')
end