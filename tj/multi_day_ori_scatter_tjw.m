%clear everything
clear all
clear all global
clc
%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P\Multi_Day_Comparisons'; %folder to save files to
dataset = 'practice_exp_list_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
d1 = 34; %day 1 in expt list
d2 = 35; %day 2 in expt list
d3 = 36; %day 3 in expt list
mouse = expt(d1).mouse; %mouse
ref_str = 'runs-001'; %string on file name to load
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

%load multiday data for days 2 and 3
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str], [date_d2 '_' mouse '_' ref_str '_' 'multiday_alignment.mat']));
d3_matches = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str], [date_d3 '_' mouse '_' ref_str '_' 'multiday_alignment.mat']));

%%
%identify matched cells and their pref oris for each day

%find cells that match from d1 to d2 and d1 to d3
match_d2 = find([d2_matches.cellImageAlign.pass]); 
match_d3 = find([d3_matches.cellImageAlign.pass]); 

%find pref oris of matched cells for d1 and d2
prefori_d1_d2_match = d1_ori.prefOri(1,match_d2);
prefori_d2_match = d2_ori.prefOri(1,match_d2);

%find pref oris of matched cells for d1 and d2
prefori_d1_d3_match = d1_ori.prefOri(1,match_d3);
prefori_d3_match = d3_ori.prefOri(1,match_d3);

%%
%separate plots of all matched cells (tuned and not tuned)

%plots on different pages
%plot d1 vs d2 pref ori
figure;
scatter(prefori_d1_d2_match,prefori_d2_match);
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
refline(1,0);
title([mouse, ' ', img_area, ' ', img_layer], 'Interpreter', 'None');

%plot d1 vs d3 pref ori
figure;
scatter(prefori_d1_d3_match,prefori_d3_match);
xlabel('Day 1 Pref Ori');
ylabel('Day 3 Pref Ori');
xlim([0,180]);
ylim([0,180]);
refline(1,0);
title([mouse, ' ', img_area, ' ', img_layer], 'Interpreter', 'None');

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
%plots of all cells across days on same page

figure;
sgtitle([mouse, ' ', img_area, ' ', img_layer], 'Interpreter', 'None');
subplot(2,1,1);
scatter(prefori_d1_d2_match,prefori_d2_match);
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
refline(1,0);
subplot(2,1,2);
scatter(prefori_d1_d3_match,prefori_d3_match);
xlabel('Day 1 Pref Ori');
ylabel('Day 3 Pref Ori');
xlim([0,180]);
ylim([0,180]);
refline(1,0);
%%
%separate cells that are tuned and matched from those that are not tuned
%but still matched - plot in separate colors and on separate figures

%find tuned cells for each day and cells that are both tuned and matched
%across days
tuned_d1 = d1_ori.ind_theta90;
tuned_d2 = d2_ori.ind_theta90;
tuned_d3 = d3_ori.ind_theta90;
tuned_matched_d2 = find(ismember(match_d2, tuned_d1) & ismember(match_d2, tuned_d2));
tuned_matched_d3 = find(ismember(match_d3, tuned_d1) & ismember(match_d3, tuned_d3));

%%
%day 1 and 2 matched and tuned vs. not tuned pref oris
figure;
scatter(prefori_d1_d2_match,prefori_d2_match);
hold on
scatter(prefori_d1_d2_match(tuned_matched_d2),prefori_d2_match(tuned_matched_d2),'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
legend('Not tuned', 'Tuned', 'Location', 'Best');
xlim([0,180]);
ylim([0,180]);
line = refline(1,0);
line.DisplayName = 'Stability';
line.Color = 'k';
title([mouse, ' ', img_area, ' ', img_layer], 'Interpreter', 'None');

%day 1 and 3 matched and tuned vs. not tuned pref oris
figure;
scatter(prefori_d1_d3_match,prefori_d3_match);
hold on
scatter(prefori_d1_d3_match(tuned_matched_d3),prefori_d3_match(tuned_matched_d3),'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
legend('Not tuned', 'Tuned', 'Location', 'Best');
xlim([0,180]);
ylim([0,180]);
line = refline(1,0);
line.DisplayName = 'Stability';
line.Color = 'k';
title([mouse, ' ', img_area, ' ', img_layer], 'Interpreter', 'None');

%%


figure('Position', [400 40 650 800]);
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
legend(['Not Tuned (n = ', num2str(length(match_d2)-length(tuned_matched_d2)), ')'], ['Tuned (n = ', num2str(length(tuned_matched_d2)), ')'], 'Location', 'northwest')
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
legend(['Not Tuned (n = ', num2str(length(match_d3)-length(tuned_matched_d3)), ')'], ['Tuned (n = ', num2str(length(tuned_matched_d3)), ')'], 'Location', 'northwest')

if expt(d1).wheel == 0
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_multi_day_matches.pdf']), '-dpdf', '-bestfit')
else
    print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_multi_day_matches.pdf']), '-dpdf', '-bestfit')
end

%%
%diff scores


d_score_d1_d2 = abs(prefori_d1_d2_match(tuned_matched_d2)-prefori_d2_match(tuned_matched_d2));
d_score_d1_d3 = abs(prefori_d1_d3_match(tuned_matched_d3)-prefori_d3_match(tuned_matched_d3));

d_scores_all = [d_score_d1_d2, d_score_d1_d3];

if expt(d1).wheel == 0
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_d_scores.mat']), 'd_scores_all')
else
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW', '_d_scores.mat']), 'd_scores_all')
end


figure; cdfplot(d_scores_all)
xlabel('Change in Pref Ori')
title([mouse, ' ', img_area, ' ', img_layer], 'Interpreter', 'None');

print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_change_ori.pdf']), '-dpdf', '-bestfit')

%%
%ideas
%make subfolders according to driver and tg vs. AAV
%try loops?