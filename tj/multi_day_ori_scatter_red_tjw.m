%clear everything
clear all
clear all global
clc
%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P\Arc_Multi_Day_Comparisons'; %folder to save files to
dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
d1 = 64; %day 1 in expt list
d2 = 65; %day 2 in expt list
d3 = 66; %day 3 in expt list
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



%%
%identify matched cells and their pref oris for each day



%find cells that match from d1 to d2 and d1 to d3 
match_d2 = d2_tcs.match_ind; 
match_d3 = d3_tcs.match_ind; 

%find pref oris of matched cells for d1 and d2
prefori_d1_d2_match = d1_ori.prefOri(1,match_d2);
prefori_d2_match = d2_ori.prefOri(1,match_d2);

%find pref oris of matched cells for d1 and d2
prefori_d1_d3_match = d1_ori.prefOri(1,match_d3);
prefori_d3_match = d3_ori.prefOri(1,match_d3);

%%
%FoV avg images of each day

figure;
sgtitle('Multi-Day FoV Avg Images')
subplot(2,2,1);
imagesc(d2_matches.fov_avg{1});
title('Day 1');
subplot(2,2,2);
imagesc(d2_matches.fov_avg{3});
title('Day 2');
subplot(2,2,3);
imagesc(d3_matches.fov_avg{3});
title('Day 3');


%%
%separate cells that are tuned and matched from those that are not tuned
%but still matched - plot in separate colors and on separate figures

%find tuned cells for each day and cells that are both tuned and matched
%across days
tuned_d1 = d1_ori.ind_theta90;
tuned_d2 = d2_ori.ind_theta90;
tuned_d3 = d3_ori.ind_theta90;
lax_tuned_matched_d2 = find(ismember(match_d2, tuned_d1) | ismember(match_d2, tuned_d2));
lax_tuned_matched_d3 = find(ismember(match_d3, tuned_d1) | ismember(match_d3, tuned_d3));

%the below criteria might be too strict
tuned_matched_d2 = find(ismember(match_d2, tuned_d1) & ismember(match_d2, tuned_d2));
tuned_matched_d3 = find(ismember(match_d3, tuned_d1) & ismember(match_d3, tuned_d3));


%%

figure('Position', [400 10 650 800]);
sgtitle([mouse, ' ', img_area, ' ', img_layer], 'Interpreter', 'None');
subplot(2,1,1);
scatter(prefori_d1_d2_match,prefori_d2_match);
hold on
scatter(prefori_d1_d2_match(tuned_matched_d2),prefori_d2_match(tuned_matched_d2),'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Not Well-Fit (n = ', num2str(length(match_d2)-length(tuned_matched_d2)), ')'], ['Well-Fit (n = ', num2str(length(tuned_matched_d2)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(prefori_d1_d3_match,prefori_d3_match);
hold on
scatter(prefori_d1_d3_match(tuned_matched_d3),prefori_d3_match(tuned_matched_d3),'r');
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
legend(['Not Well-Fit (n = ', num2str(length(match_d3)-length(tuned_matched_d3)), ')'], ['Well-Fit (n = ', num2str(length(tuned_matched_d3)), ')'], 'Location', 'northwest')

if str2double(expt(d1).img_day) == 1
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_change_ori_scatter.pdf']), '-dpdf', '-bestfit')
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_ori_changes.mat']), 'prefori_d1_d2_match', 'prefori_d1_d3_match', 'prefori_d2_match', 'prefori_d3_match', 'lax_tuned_matched_d2', 'lax_tuned_matched_d3')
else
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_change_ori_scatter_2.pdf']), '-dpdf', '-bestfit')
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_ori_changes_2.mat']), 'prefori_d1_d2_match', 'prefori_d1_d3_match', 'prefori_d2_match', 'prefori_d3_match', 'lax_tuned_matched_d2', 'lax_tuned_matched_d3')

end

%%
%diff scores


%d1_2
dscore_prefori_d1_d2_match = double(prefori_d1_d2_match>90);
dscore_prefori_d1_d2_match(dscore_prefori_d1_d2_match>0) = 180;
dscore_prefori_d1_d2_match = abs(dscore_prefori_d1_d2_match-prefori_d1_d2_match);
%d2
dscore_prefori_d2_match = double(prefori_d2_match>90);
dscore_prefori_d2_match(dscore_prefori_d2_match>0) = 180;
dscore_prefori_d2_match = abs(dscore_prefori_d2_match-prefori_d2_match);
%d1_d3
dscore_prefori_d1_d3_match = double(prefori_d1_d3_match>90);
dscore_prefori_d1_d3_match(dscore_prefori_d1_d3_match>0) = 180;
dscore_prefori_d1_d3_match = abs(dscore_prefori_d1_d3_match-prefori_d1_d3_match);
%d3
dscore_prefori_d3_match = double(prefori_d3_match>90);
dscore_prefori_d3_match(dscore_prefori_d3_match>0) = 180;
dscore_prefori_d3_match = abs(dscore_prefori_d3_match-prefori_d3_match);



d_score_prefori_d1_d2 = abs(dscore_prefori_d1_d2_match(lax_tuned_matched_d2)-dscore_prefori_d2_match(lax_tuned_matched_d2));
d_score_prefori_d1_d3 = abs(dscore_prefori_d1_d3_match(lax_tuned_matched_d3)-dscore_prefori_d3_match(lax_tuned_matched_d3));

d_scores_prefori_all = [d_score_prefori_d1_d2, d_score_prefori_d1_d3];

figure; 
h=cdfplot(d_score_prefori_d1_d2);
hold on
j=cdfplot(d_score_prefori_d1_d3);
set(h, 'LineStyle', '-', 'Color', 'b');
set(j, 'LineStyle', '--', 'Color', 'b');
hold off
legend(['Day 1 vs. Day 2 (n = ', num2str(length(d_score_prefori_d1_d2)), ')'], ['Day 1 vs. Day 3 (n = ', num2str(length(d_score_prefori_d1_d3)), ')'], 'Location', 'best')
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
k_d1_d2_match = d1_k_max.k1(1,match_d2);
k_d2_match = d2_k_max.k1(1,match_d2);

%find pref oris of matched cells for d1 and d2
k_d1_d3_match = d1_k_max.k1(1,match_d3);
k_d3_match = d3_k_max.k1(1,match_d3);


figure('Position', [400 10 650 800]);
sgtitle([mouse, ' ', img_area, ' ', img_layer], 'Interpreter', 'None');
subplot(2,1,1);
scatter(k_d1_d2_match,k_d2_match);
hold on
scatter(k_d1_d2_match(tuned_matched_d2),k_d2_match(tuned_matched_d2),'r');
hold off
xlabel('Day 1 k Value');
ylabel('Day 2 k Value');
xlim([0,30]);
ylim([0,30]);
%xticks(0:20:180);
%yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Not Well-Fit (n = ', num2str(length(match_d2)-length(tuned_matched_d2)), ')'], ['Well-Fit (n = ', num2str(length(tuned_matched_d2)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(k_d1_d3_match,k_d3_match);
hold on
scatter(k_d1_d3_match(tuned_matched_d3),k_d3_match(tuned_matched_d3),'r');
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
legend(['Not Well-Fit (n = ', num2str(length(match_d3)-length(tuned_matched_d3)), ')'], ['Well-Fit (n = ', num2str(length(tuned_matched_d3)), ')'], 'Location', 'northwest')

if str2double(expt(d1).img_day) == 1
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW_change_k_scatter.pdf']), '-dpdf', '-bestfit')
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_k_changes.mat']), 'k_d1_d2_match', 'k_d1_d3_match', 'k_d2_match', 'k_d3_match')
else
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW_change_k_scatter_2.pdf']), '-dpdf', '-bestfit')
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_k_changes_2.mat']), 'k_d1_d2_match', 'k_d1_d3_match', 'k_d2_match', 'k_d3_match')
end

%d1_2
dscore_k_d1_d2_match = abs(k_d1_d2_match(lax_tuned_matched_d2)-k_d2_match(lax_tuned_matched_d2));
%d1_d3
dscore_k_d1_d3_match = abs(k_d1_d3_match(lax_tuned_matched_d3)-k_d3_match(lax_tuned_matched_d3));

figure; 
h=cdfplot(dscore_k_d1_d2_match);
hold on
j=cdfplot(dscore_k_d1_d3_match);
set(h, 'LineStyle', '-', 'Color', 'b');
set(j, 'LineStyle', '--', 'Color', 'b');
hold off
legend(['Day 1 vs. Day 2 (n = ', num2str(length(dscore_k_d1_d2_match)), ')'], ['Day 1 vs. Day 3 (n = ', num2str(length(dscore_k_d1_d3_match)), ')'], 'Location', 'best')
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
max_d1_d2_match = d1_k_max.max_dfof(1,match_d2);
max_d2_match = d2_k_max.max_dfof(1,match_d2);

%find pref oris of matched cells for d1 and d2
max_d1_d3_match = d1_k_max.max_dfof(1,match_d3);
max_d3_match = d3_k_max.max_dfof(1,match_d3);


figure('Position', [400 10 650 800]);
sgtitle([mouse, ' ', img_area, ' ', img_layer], 'Interpreter', 'None');
subplot(2,1,1);
scatter(max_d1_d2_match,max_d2_match);
hold on
scatter(max_d1_d2_match(tuned_matched_d2),max_d2_match(tuned_matched_d2),'r');
hold off
xlabel('Day 1 Max dF/F Value');
ylabel('Day 2 Max dF/F Value');
xlim([0,1]);
ylim([0,1]);
%xticks(0:20:180);
%yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Not Well-Fit (n = ', num2str(length(match_d2)-length(tuned_matched_d2)), ')'], ['Well-Fit (n = ', num2str(length(tuned_matched_d2)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(max_d1_d3_match,max_d3_match);
hold on
scatter(max_d1_d3_match(tuned_matched_d3),max_d3_match(tuned_matched_d3),'r');
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
legend(['Not Well-Fit (n = ', num2str(length(match_d3)-length(tuned_matched_d3)), ')'], ['Well-Fit (n = ', num2str(length(tuned_matched_d3)), ')'], 'Location', 'northwest')

if str2double(expt(d1).img_day) == 1
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW_change_max_scatter.pdf']), '-dpdf', '-bestfit')
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_max_changes.mat']), 'max_d1_d2_match', 'max_d1_d3_match', 'max_d2_match', 'max_d3_match')

else
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW_change_max_scatter_2.pdf']), '-dpdf', '-bestfit')
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_max_changes_2.mat']), 'max_d1_d2_match', 'max_d1_d3_match', 'max_d2_match', 'max_d3_match')

end

%d1_2
dscore_max_d1_d2_match = abs(max_d1_d2_match(lax_tuned_matched_d2)-max_d2_match(lax_tuned_matched_d2));
%d1_d3
dscore_max_d1_d3_match = abs(max_d1_d3_match(lax_tuned_matched_d3)-max_d3_match(lax_tuned_matched_d3));

figure; 
h=cdfplot(dscore_max_d1_d2_match);
hold on
j=cdfplot(dscore_max_d1_d3_match);
set(h, 'LineStyle', '-', 'Color', 'b');
set(j, 'LineStyle', '--', 'Color', 'b');
hold off
legend(['Day 1 vs. Day 2 (n = ', num2str(length(dscore_max_d1_d2_match)), ')'], ['Day 1 vs. Day 3 (n = ', num2str(length(dscore_max_d1_d3_match)), ')'], 'Location', 'best')
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
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_d_scores.mat']), 'd_score_prefori_d1_d2', 'd_score_prefori_d1_d3', 'dscore_k_d1_d2_match', 'dscore_k_d1_d3_match', 'dscore_max_d1_d2_match', 'dscore_max_d1_d3_match')
else
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_d_scores_2.mat']), 'd_score_prefori_d1_d2', 'd_score_prefori_d1_d3', 'dscore_k_d1_d2_match', 'dscore_k_d1_d3_match', 'dscore_max_d1_d2_match', 'dscore_max_d1_d3_match')
end