%This code is used to analyze individual subject data from the SRP
%experiment
%%
%clear everything
clear all
clear all global
clc
close all

%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\phaserev\multi_day'; %folder to save files to
dataset = 'exp_list_phaserev_tjw'; %experiment list to pick files from
eval(dataset); %load dataset


%%

%we start with loading the proper information based on the mouse's info in
%the experiment list

%baseline1 to baseline1+4d comparison
d1 = 1+7; %base1 in expt list
d2 = 2+7; %base1+4d in expt list
d3 = 3+7;
d4 = 4+7;
d5 = 5+7;
d6 = 6+7;
d7 = 7+7;
mouse = expt(d1).mouse; 
ref_str_d1 = ['runs-',expt(d1).runs];
ref_str_d2 = ['runs-',expt(d2).runs]; 
ref_str_d3 = ['runs-',expt(d3).runs]; 
ref_str_d4 = ['runs-',expt(d4).runs]; 
ref_str_d5 = ['runs-',expt(d5).runs]; 
ref_str_d6 = ['runs-',expt(d6).runs]; 
ref_str_d7 = ['runs-',expt(d7).runs]; 
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date_d1 = expt(d1).date; 
date_d2 = expt(d2).date; 
date_d3 = expt(d3).date; 
date_d4 = expt(d4).date; 
date_d5 = expt(d5).date; 
date_d6 = expt(d6).date; 
date_d7 = expt(d7).date; 


%load multiday data for day 2
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));
d3_matches = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'multiday_alignment.mat']));
d4_matches = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'multiday_alignment.mat']));
d5_matches = load(fullfile(fnout, [date_d5 '_' mouse], [date_d5 '_' mouse '_' ref_str_d5], [date_d5 '_' mouse '_' ref_str_d5 '_' 'multiday_alignment.mat']));
d6_matches = load(fullfile(fnout, [date_d6 '_' mouse], [date_d6 '_' mouse '_' ref_str_d6], [date_d6 '_' mouse '_' ref_str_d6 '_' 'multiday_alignment.mat']));
d7_matches = load(fullfile(fnout, [date_d7 '_' mouse], [date_d7 '_' mouse '_' ref_str_d7], [date_d7 '_' mouse '_' ref_str_d7 '_' 'multiday_alignment.mat']));


%load avg across trial and across cell data
d2_avgs = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'keyvars.mat']));
d3_avgs = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'keyvars.mat']));
d4_avgs = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'keyvars.mat']));
d5_avgs = load(fullfile(fnout, [date_d5 '_' mouse], [date_d5 '_' mouse '_' ref_str_d5], [date_d5 '_' mouse '_' ref_str_d5 '_' 'keyvars.mat']));
d6_avgs = load(fullfile(fnout, [date_d6 '_' mouse], [date_d6 '_' mouse '_' ref_str_d6], [date_d6 '_' mouse '_' ref_str_d6 '_' 'keyvars.mat']));


%load TCs
d1_TCs = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'TCs.mat']));
d2_TCs = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'TCs.mat']));
d3_TCs = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'TCs.mat']));
d4_TCs = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'TCs.mat']));
d5_TCs = load(fullfile(fnout, [date_d5 '_' mouse], [date_d5 '_' mouse '_' ref_str_d5], [date_d5 '_' mouse '_' ref_str_d5 '_' 'TCs.mat']));
d6_TCs = load(fullfile(fnout, [date_d6 '_' mouse], [date_d6 '_' mouse '_' ref_str_d6], [date_d6 '_' mouse '_' ref_str_d6 '_' 'TCs.mat']));
d7_TCs = load(fullfile(fnout, [date_d7 '_' mouse], [date_d7 '_' mouse '_' ref_str_d7], [date_d7 '_' mouse '_' ref_str_d7 '_' 'TCs.mat']));



%load an input file for info
load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'input.mat']));


%% matching indices
match_d2 = find([d2_matches.cellImageAlign.pass]); 
match_d3 = find([d3_matches.cellImageAlign.pass]); 
match_d4 = find([d4_matches.cellImageAlign.pass]); 
match_d5 = find([d5_matches.cellImageAlign.pass]); 
match_d6 = find([d6_matches.cellImageAlign.pass]); 
match_d7 = find([d7_matches.cellImageAlign.pass]); 


match_all =  intersect(intersect(intersect(intersect(match_d2,match_d3),match_d4),match_d5),match_d6);
% match_all_excl = intersect(intersect(intersect(match_d2,match_d4),match_d5),match_d6);
match_all_excl = intersect(intersect(intersect(intersect(match_d2,match_d4),match_d5),match_d6),match_d7); %including session 7 and excluding session 2

%%
nOn = input.nScansOn;
nOff = input.nScansOff;
nTrials = input.stopAfterNTrials;

[d2_pop_avg_per_trial, d2_pop_std_per_trial, d2_pop_se_per_trial, d2_cell_avg_trial, d2_cell_std_trial, d2_cell_se_trial] = tjgetmatchavgs(d2_TCs, match_all_excl, nOn, nOff, nTrials);
%[d3_pop_avg_per_trial, d3_pop_std_per_trial, d3_pop_se_per_trial, d3_cell_avg_trial, d3_cell_std_trial, d3_cell_se_trial] = tjgetmatchavgs(d3_TCs, match_all_excl, nOn, nOff, nTrials);
[d4_pop_avg_per_trial, d4_pop_std_per_trial, d4_pop_se_per_trial, d4_cell_avg_trial, d4_cell_std_trial, d4_cell_se_trial] = tjgetmatchavgs(d4_TCs, match_all_excl, nOn, nOff, nTrials);
[d5_pop_avg_per_trial, d5_pop_std_per_trial, d5_pop_se_per_trial, d5_cell_avg_trial, d5_cell_std_trial, d5_cell_se_trial] = tjgetmatchavgs(d5_TCs, match_all_excl, nOn, nOff, nTrials);
[d6_pop_avg_per_trial, d6_pop_std_per_trial, d6_pop_se_per_trial, d6_cell_avg_trial, d6_cell_std_trial, d6_cell_se_trial] = tjgetmatchavgs(d6_TCs, match_all_excl, nOn, nOff, nTrials);



%%
%these look like crap
% figure; 
% e2 = errorbar(d2_pop_avg_per_trial, d2_pop_se_per_trial);
% set(e2, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
% hold on;
% e4 = errorbar(d4_pop_avg_per_trial, d4_pop_se_per_trial);
% set(e4, 'LineStyle', '-', 'Color', [.25 .25 .25], 'LineWidth',2);
% e5 = errorbar(d5_pop_avg_per_trial, d5_pop_se_per_trial);
% set(e5, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',2);

%%
figure;
a = cdfplot(d2_pop_avg_per_trial);
set(a, 'Color', [0 0 0], 'LineWidth', 2)
hold on;
b = cdfplot(d4_pop_avg_per_trial);
set(b, 'Color', [.25 .25 .25], 'LineWidth', 2)
c = cdfplot(d5_pop_avg_per_trial);
set(c, 'Color', [.5 .5 .5], 'LineWidth', 2)
d = cdfplot(d6_pop_avg_per_trial);
set(d, 'Color', [.75 .75 .75], 'LineWidth', 2)

%%
figure;
subplot(2,2,1)
histogram(d2_pop_avg_per_trial)
subplot(2,2,2)
histogram(d4_pop_avg_per_trial)
subplot(2,2,3)
histogram(d5_pop_avg_per_trial)
subplot(2,2,4)
histogram(d6_pop_avg_per_trial)

%%
figure;
a = cdfplot(d2_cell_avg_trial);
set(a, 'Color', [0 0 0], 'LineWidth', 2)
hold on;
b = cdfplot(d4_cell_avg_trial);
set(b, 'Color', [.25 .25 .25], 'LineWidth', 2)
c = cdfplot(d5_cell_avg_trial);
set(c, 'Color', [.5 .5 .5], 'LineWidth', 2)
d = cdfplot(d6_cell_avg_trial);
set(d, 'Color', [.75 .75 .75], 'LineWidth', 2)

%%
figure;
subplot(2,2,1)
histogram(d2_cell_avg_trial)
subplot(2,2,2)
histogram(d4_cell_avg_trial)
subplot(2,2,3)
histogram(d5_cell_avg_trial)
subplot(2,2,4)
histogram(d6_cell_avg_trial)



%%
%maybe do cells specific to day 1 and day 5 or just do those comparisons

%%
d2_grand_avg = mean(d2_pop_avg_per_trial);
d2_grand_se = std(d2_pop_avg_per_trial)./sqrt(size(d2_pop_avg_per_trial,1));
d4_grand_avg = mean(d4_pop_avg_per_trial);
d4_grand_se = std(d4_pop_avg_per_trial)./sqrt(size(d4_pop_avg_per_trial,1));
d5_grand_avg = mean(d5_pop_avg_per_trial);
d5_grand_se = std(d5_pop_avg_per_trial)./sqrt(size(d5_pop_avg_per_trial,1));
d6_grand_avg = mean(d6_pop_avg_per_trial);
d6_grand_se = std(d6_pop_avg_per_trial)./sqrt(size(d6_pop_avg_per_trial,1));

grand_avgs = [d2_grand_avg d4_grand_avg, d5_grand_avg, d6_grand_avg]
grand_ses = [d2_grand_se d4_grand_se, d5_grand_se, d6_grand_se]


%%
figure;
a = errorbar(grand_avgs, grand_ses);
set(a, 'Color', [0 0 0], 'LineWidth', 2)
xlim([0 5])
ylim([-.02 .1])
% xlabel('Day')
ylabel('Avg dF/F')
xticklabels({'','','Day 1 SRP','','Day 3 SRP','','Day 4 SRP','','Day 5 SRP'})
print(fullfile(newfnout, [mouse,  '_avg_across_days.pdf']), '-dpdf', '-bestfit')


%%
% d2 vs d6 only
match_d2_d6 = intersect(match_d2, match_d6);

[comp_d2_pop_avg_per_trial, comp_d2_pop_std_per_trial, comp_d2_pop_se_per_trial, comp_d2_cell_avg_trial, comp_d2_cell_std_trial, comp_d2_cell_se_trial] = tjgetmatchavgs(d2_TCs, match_d2_d6, nOn, nOff, nTrials);
[comp_d6_pop_avg_per_trial, comp_d6_pop_std_per_trial, comp_d6_pop_se_per_trial, comp_d6_cell_avg_trial, comp_d6_cell_std_trial, comp_d6_cell_se_trial] = tjgetmatchavgs(d6_TCs, match_d2_d6, nOn, nOff, nTrials);


comp_d2_grand_avg = mean(comp_d2_pop_avg_per_trial);
comp_d2_grand_se = std(comp_d2_pop_avg_per_trial)./sqrt(size(comp_d2_pop_avg_per_trial,1));
comp_d6_grand_avg = mean(comp_d6_pop_avg_per_trial);
comp_d6_grand_se = std(comp_d6_pop_avg_per_trial)./sqrt(size(comp_d6_pop_avg_per_trial,1));

comp_grand_avgs = [comp_d2_grand_avg comp_d6_grand_avg];
comp_grand_ses = [comp_d2_grand_se comp_d6_grand_se];

figure;
a = errorbar(comp_grand_avgs, comp_grand_ses);
set(a, 'Color', [0 0 0], 'LineWidth', 2)
xlim([0 3])
xticklabels({'','',1,'',5})
xlabel('Day')
ylabel('Avg dF/F')

%%

[d2_pop_avg_per_trial, d2_pop_std_per_trial, d2_pop_se_per_trial, d2_cell_avg_trial, d2_cell_std_trial, d2_cell_se_trial] = tjgetmatchavgs(d2_TCs, match_d2, nOn, nOff, nTrials);
[d3_pop_avg_per_trial, d3_pop_std_per_trial, d3_pop_se_per_trial, d3_cell_avg_trial, d3_cell_std_trial, d3_cell_se_trial] = tjgetmatchavgs(d3_TCs, match_d3, nOn, nOff, nTrials);
[d4_pop_avg_per_trial, d4_pop_std_per_trial, d4_pop_se_per_trial, d4_cell_avg_trial, d4_cell_std_trial, d4_cell_se_trial] = tjgetmatchavgs(d4_TCs, match_d4, nOn, nOff, nTrials);
[d5_pop_avg_per_trial, d5_pop_std_per_trial, d5_pop_se_per_trial, d5_cell_avg_trial, d5_cell_std_trial, d5_cell_se_trial] = tjgetmatchavgs(d5_TCs, match_d5, nOn, nOff, nTrials);
[d6_pop_avg_per_trial, d6_pop_std_per_trial, d6_pop_se_per_trial, d6_cell_avg_trial, d6_cell_std_trial, d6_cell_se_trial] = tjgetmatchavgs(d6_TCs, match_d6, nOn, nOff, nTrials);

pop_avg_per_trial_all = [d2_pop_avg_per_trial d3_pop_avg_per_trial d4_pop_avg_per_trial d5_pop_avg_per_trial d6_pop_avg_per_trial];
pop_se_per_trial_all = [d2_pop_se_per_trial d3_pop_se_per_trial d4_pop_se_per_trial d5_pop_se_per_trial d6_pop_se_per_trial];
%%
pop_avg_per_trial_all = reshape(pop_avg_per_trial_all,[],1);
pop_se_per_trial_all = reshape(pop_se_per_trial_all,[],1);
%%
figure;
a=errorbar(pop_avg_per_trial_all, pop_se_per_trial_all)
set(a, 'Color', [0 0 0], 'LineWidth', 1)
hold on
vline([50 100 150 200])

%%
figure;
a=plot(pop_avg_per_trial_all)
set(a, 'Color', [0 0 0], 'LineWidth', 1)
hold on
eb = errorbar(pop_avg_per_trial_all,pop_se_per_trial_all,'LineStyle','none'); 
set(eb, 'Color', [.5 .5 .5], 'LineWidth', 1)
hold off
ylim([-.1 .4])
xlabel('Trial')
ylabel('Population dF/F')
vline([50 100 150 200])
print(fullfile(newfnout, [mouse,  '_all_days_avg_across_trial.pdf']), '-dpdf', '-bestfit')


%%
%ALL cells

[d2_pop_avg_per_trial, d2_pop_std_per_trial, d2_pop_se_per_trial, d2_cell_avg_trial, d2_cell_std_trial, d2_cell_se_trial] = tjgetmatchavgs(d2_TCs, match_all, nOn, nOff, nTrials);
[d3_pop_avg_per_trial, d3_pop_std_per_trial, d3_pop_se_per_trial, d3_cell_avg_trial, d3_cell_std_trial, d3_cell_se_trial] = tjgetmatchavgs(d3_TCs, match_all, nOn, nOff, nTrials);
[d4_pop_avg_per_trial, d4_pop_std_per_trial, d4_pop_se_per_trial, d4_cell_avg_trial, d4_cell_std_trial, d4_cell_se_trial] = tjgetmatchavgs(d4_TCs, match_all, nOn, nOff, nTrials);
[d5_pop_avg_per_trial, d5_pop_std_per_trial, d5_pop_se_per_trial, d5_cell_avg_trial, d5_cell_std_trial, d5_cell_se_trial] = tjgetmatchavgs(d5_TCs, match_all, nOn, nOff, nTrials);
[d6_pop_avg_per_trial, d6_pop_std_per_trial, d6_pop_se_per_trial, d6_cell_avg_trial, d6_cell_std_trial, d6_cell_se_trial] = tjgetmatchavgs(d6_TCs, match_all, nOn, nOff, nTrials);

d2_grand_avg = mean(d2_pop_avg_per_trial);
d2_grand_se = std(d2_pop_avg_per_trial)./sqrt(size(d2_pop_avg_per_trial,1));
d3_grand_avg = mean(d3_pop_avg_per_trial);
d3_grand_se = std(d3_pop_avg_per_trial)./sqrt(size(d3_pop_avg_per_trial,1));
d4_grand_avg = mean(d4_pop_avg_per_trial);
d4_grand_se = std(d4_pop_avg_per_trial)./sqrt(size(d4_pop_avg_per_trial,1));
d5_grand_avg = mean(d5_pop_avg_per_trial);
d5_grand_se = std(d5_pop_avg_per_trial)./sqrt(size(d5_pop_avg_per_trial,1));
d6_grand_avg = mean(d6_pop_avg_per_trial);
d6_grand_se = std(d6_pop_avg_per_trial)./sqrt(size(d6_pop_avg_per_trial,1));

grand_avgs = [d2_grand_avg d3_grand_avg d4_grand_avg, d5_grand_avg, d6_grand_avg]
grand_ses = [d2_grand_se d3_grand_se d4_grand_se, d5_grand_se, d6_grand_se]


figure;
a = errorbar(grand_avgs, grand_ses);
set(a, 'Color', [0 0 0], 'LineWidth', 2)
xlim([0 5])
xlabel('Day')
ylabel('Avg dF/F')

%%
%incorporate d1 and d7

d1_pref = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'prefori.mat']));
d7_pref = load(fullfile(fnout, [date_d7 '_' mouse], [date_d7 '_' mouse '_' ref_str_d7], [date_d7 '_' mouse '_' ref_str_d7 '_' 'prefori.mat']));

d1_pref_match = d1_pref.pref_ori(:,match_d7);
d7_pref_match = d7_pref.pref_ori(:,match_d7);

figure;
scatter(d1_pref_match, d7_pref_match)
xlabel('Day 1 Pref Ori')
ylabel('Day 7 Pref Ori')
%%
dscore_prefori_d1_d7_match = double(d1_pref_match>90);
dscore_prefori_d1_d7_match(dscore_prefori_d1_d7_match>0) = 180;
dscore_prefori_d1_d7_match = abs(dscore_prefori_d1_d7_match-d1_pref_match);
%d7
dscore_prefori_d7_match = double(d7_pref_match>90);
dscore_prefori_d7_match(dscore_prefori_d7_match>0) = 180;
dscore_prefori_d7_match = abs(dscore_prefori_d7_match-d7_pref_match);

d_score_prefori_d1_d7 = abs(dscore_prefori_d1_d7_match-dscore_prefori_d7_match);

figure;
cdfplot(d_score_prefori_d1_d7)
xlabel('Change in Pref Ori')
ylabel('% of Cells')
title(mouse)
print(fullfile(newfnout, [mouse,  '_change_pref_d1_d7.pdf']), '-dpdf', '-bestfit')

%%
match_all_excl_plus_oris = intersect(intersect(intersect(intersect(match_d2,match_d4),match_d5),match_d6),match_d7);

input_d1 = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'input.mat']));
input_d7 = load(fullfile(fnout, [date_d7 '_' mouse], [date_d7 '_' mouse '_' ref_str_d7], [date_d7 '_' mouse '_' ref_str_d7 '_' 'input.mat']));

ori_d1 = celleqel2mat_padded(input_d1.input.tGratingDirectionDeg);
ori_d7 = celleqel2mat_padded(input_d7.input.tGratingDirectionDeg);
oris = unique(ori_d1);

d1_45_id = find(ori_d1 == 45);
d7_45_id = find(ori_d7 == 45);
d1_135_id = find(ori_d1 == 135);
d7_135_id = find(ori_d7 == 135);

d1_tcs_match = d1_TCs.npSub_tc(:,match_all_excl_plus_oris); %all frames for matched cells throughout
d7_tcs_match = d7_TCs.cellTCs_match{1,2}(:,match_all_excl_plus_oris);

nCells = size(d1_tcs_match,2);
nTrials = size(input_d1.input.tGratingDirectionDeg,2);
nOff = input_d1.input.nScansOff;
nOn = input_d7.input.nScansOn;

%%

d1_data_trial = permute(reshape(d1_tcs_match,[nOff + nOn nTrials nCells]),[1 3 2]);
d7_data_trial = permute(reshape(d7_tcs_match,[nOff + nOn nTrials nCells]),[1 3 2]);

d1_data_trial_45 = d1_data_trial(:,:,d1_45_id);
d7_data_trial_45 = d7_data_trial(:,:,d7_45_id);

d1_data_f_45 = mean(d1_data_trial_45(nOff-15:nOff,:,:),1);
d1_data_dfof_45 = (d1_data_trial_45-d1_data_f_45)./d1_data_f_45;
d7_data_f_45 = mean(d7_data_trial_45(nOff-15:nOff,:,:),1);
d7_data_dfof_45 = (d7_data_trial_45-d7_data_f_45)./d7_data_f_45;

d1_pop_avg_per_trial_45 = mean(d1_data_dfof_45(nOff+1:nOn,:,:),1, 'omitnan');
d1_pop_std_per_trial_45 = std(d1_pop_avg_per_trial_45, 'omitnan');
d1_pop_se_per_trial_45 = d1_pop_std_per_trial_45./sqrt(size(d1_pop_avg_per_trial_45,2));

d1_pop_avg_per_trial_45 = mean(d1_pop_avg_per_trial_45,2, 'omitnan');

d1_pop_avg_per_trial_45 = squeeze(d1_pop_avg_per_trial_45);
d1_pop_std_per_trial_45 = squeeze(d1_pop_std_per_trial_45);
d1_pop_se_per_trial_45 = squeeze(d1_pop_se_per_trial_45);

d7_pop_avg_per_trial_45 = mean(d7_data_dfof_45(nOff+1:nOn,:,:),1, 'omitnan');
d7_pop_std_per_trial_45 = std(d7_pop_avg_per_trial_45, 'omitnan');
d7_pop_se_per_trial_45 = d7_pop_std_per_trial_45./sqrt(size(d7_pop_avg_per_trial_45,2));

d7_pop_avg_per_trial_45 = mean(d7_pop_avg_per_trial_45,2, 'omitnan');

d7_pop_avg_per_trial_45 = squeeze(d7_pop_avg_per_trial_45);
d7_pop_std_per_trial_45 = squeeze(d7_pop_std_per_trial_45);
d7_pop_se_per_trial_45 = squeeze(d7_pop_se_per_trial_45);

d1_grand_avg_45 = mean(d1_pop_avg_per_trial_45);
d1_grand_se_45 = std(d1_pop_avg_per_trial_45)./sqrt(size(d1_pop_avg_per_trial_45,1));

d7_grand_avg_45 = mean(d7_pop_avg_per_trial_45);
d7_grand_se_45 = std(d7_pop_avg_per_trial_45)./sqrt(size(d7_pop_avg_per_trial_45,1));

figure;
errorbar([d1_grand_avg_45, d7_grand_avg_45], [d1_grand_se_45, d7_grand_se_45])
xlim([0 3])
xticklabels({'','',1,'',5,'',''})


%% now do the above w 135
d1_data_trial_135 = d1_data_trial(:,:,d1_135_id);
d7_data_trial_135 = d7_data_trial(:,:,d7_135_id);

d1_data_f_135 = mean(d1_data_trial_135(nOff-15:nOff,:,:),1);
d1_data_dfof_135 = (d1_data_trial_135-d1_data_f_135)./d1_data_f_135;
d7_data_f_135 = mean(d7_data_trial_135(nOff-15:nOff,:,:),1);
d7_data_dfof_135 = (d7_data_trial_135-d7_data_f_135)./d7_data_f_135;

d1_pop_avg_per_trial_135 = mean(d1_data_dfof_135(nOff+1:nOn,:,:),1, 'omitnan');
d1_pop_std_per_trial_135 = std(d1_pop_avg_per_trial_135, 'omitnan');
d1_pop_se_per_trial_135 = d1_pop_std_per_trial_135./sqrt(size(d1_pop_avg_per_trial_135,2));

d1_pop_avg_per_trial_135 = mean(d1_pop_avg_per_trial_135,2, 'omitnan');

d1_pop_avg_per_trial_135 = squeeze(d1_pop_avg_per_trial_135);
d1_pop_std_per_trial_135 = squeeze(d1_pop_std_per_trial_135);
d1_pop_se_per_trial_135 = squeeze(d1_pop_se_per_trial_135);

d7_pop_avg_per_trial_135 = mean(d7_data_dfof_135(nOff+1:nOn,:,:),1, 'omitnan');
d7_pop_std_per_trial_135 = std(d7_pop_avg_per_trial_135, 'omitnan');
d7_pop_se_per_trial_135 = d7_pop_std_per_trial_135./sqrt(size(d7_pop_avg_per_trial_135,2));

d7_pop_avg_per_trial_135 = mean(d7_pop_avg_per_trial_135,2, 'omitnan');

d7_pop_avg_per_trial_135 = squeeze(d7_pop_avg_per_trial_135);
d7_pop_std_per_trial_135 = squeeze(d7_pop_std_per_trial_135);
d7_pop_se_per_trial_135 = squeeze(d7_pop_se_per_trial_135);

d1_grand_avg_135 = mean(d1_pop_avg_per_trial_135);
d1_grand_se_135 = std(d1_pop_avg_per_trial_135)./sqrt(size(d1_pop_avg_per_trial_135,1));

d7_grand_avg_135 = mean(d7_pop_avg_per_trial_135);
d7_grand_se_135 = std(d7_pop_avg_per_trial_135)./sqrt(size(d7_pop_avg_per_trial_135,1));

figure;
errorbar([d1_grand_avg_135, d7_grand_avg_135], [d1_grand_se_135, d7_grand_se_135])
xlim([0 3])
xticklabels({'','',1,'',5,'',''})

%%
fig = figure;
a=subplot(1,2,1)
c = errorbar([d1_grand_avg_45, d7_grand_avg_45], [d1_grand_se_45, d7_grand_se_45])
set(c, 'Color', [0 0 0], 'LineWidth', 2)
xlim([0 3])
xticklabels({'',1,5,''})
ylim([-.01 .05])
title('45 Degree Stim')
b=subplot(1,2,2)
d = errorbar([d1_grand_avg_135, d7_grand_avg_135], [d1_grand_se_135, d7_grand_se_135])
set(d, 'Color', [0 0 0], 'LineWidth', 2)
xlim([0 3])
xticklabels({'',1,5,''})
ylim([-.01 .05])
title('135 Degree Stim')

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
linkaxes([a b],'xy')
ylabel(han,'Avg dF/F');
xlabel(han,'Day');

%%
figure;
a = errorbar([d1_grand_avg_45, grand_avgs, d7_grand_avg_45], [d1_grand_se_45, grand_ses, d1_grand_se_45]);
set(a, 'Color', [0 0 0], 'LineWidth', 2)
xlim([0 7])
xticklabels({'','Day 1 8ori', 'Day 1 SRP', 'Day 3 SRP', 'Day 4 SRP', 'Day 5 SRP', 'Day 5 8ori'})
ylabel('Avg dF/F')
ylim([-.02 .1])
vline([2 5])
title('Population Response to 45 Degree Stimulus')
print(fullfile(newfnout, [mouse,  '_SRP_timeline.pdf']), '-dpdf', '-bestfit')


%% single cells d1 and d7

% input_d1 = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'input.mat']));
% nCells = size(match_all_excl,2);
% nTrials = size(input_d1.input.tGratingDirectionDeg,2);
% nOff = input_d1.input.nScansOff;
% nOn = input_d1.input.nScansOn;
% 
% d1_tc_match = permute(reshape(d1_TCs.npSub_tc(:,match_all_excl), [nOff + nOn nTrials nCells]),[1 3 2]);
% d7_tc_match = permute(reshape(d7_TCs.cellTCs_match{1,2}(:,match_all_excl), [nOff + nOn nTrials nCells]),[1 3 2]);
% 
% d1_tc_match = mean(d1_tc_match(nOff-15:nOff,:,:),1);
% d1_tc_match = squeeze(d1_tc_match);
% d7_tc_match = mean(d7_tc_match(nOff-15:nOff,:,:),1);
% d7_tc_match = squeeze(d7_tc_match);
% 
% %%
% d1_tc_mean = mean(d1_tc_match(:,d1_45_id),2);
% d1_tc_se = std(d1_tc_match(:,d1_45_id),[],2, 'omitnan')./sqrt(size(d1_tc_match(:,d1_45_id),2))
% 

% 

%%

d1_data_trial = permute(reshape(d1_tcs_match,[nOff + nOn nTrials nCells]),[1 3 2]);
d7_data_trial = permute(reshape(d7_tcs_match,[nOff + nOn nTrials nCells]),[1 3 2]);

d1_data_trial_45 = d1_data_trial(:,:,d1_45_id);
d7_data_trial_45 = d7_data_trial(:,:,d7_45_id);

%%
d1_data_f_45 = mean(d1_data_trial_45(nOff-15:nOff,:,:),1);
d1_data_dfof_45 = (d1_data_trial_45-d1_data_f_45)./d1_data_f_45;
d7_data_f_45 = mean(d7_data_trial_45(nOff-15:nOff,:,:),1);
d7_data_dfof_45 = (d7_data_trial_45-d7_data_f_45)./d7_data_f_45;
%%
d1_data_dfof_45 = mean(d1_data_dfof_45(nOff+1:nOn,:,:),1);
d1_data_dfof_45 = squeeze(d1_data_dfof_45);

d1_mean_data_dfof_45 = mean(d1_data_dfof_45,2);
d1_se_data_dfof_45 = std(d1_data_dfof_45,[],2,'omitnan')./sqrt(size(d1_data_dfof_45,2));

d7_data_dfof_45 = mean(d7_data_dfof_45(nOff+1:nOn,:,:),1);
d7_data_dfof_45 = squeeze(d7_data_dfof_45);

d7_mean_data_dfof_45 = mean(d7_data_dfof_45,2);
d7_se_data_dfof_45 = std(d7_data_dfof_45,[],2,'omitnan')./sqrt(size(d7_data_dfof_45,2));

%% SRP days
d2_tcs_match = d2_TCs.cellTCs_match{1,2}(:,match_all_excl_plus_oris);
d4_tcs_match = d4_TCs.cellTCs_match{1,2}(:,match_all_excl_plus_oris);
d5_tcs_match = d5_TCs.cellTCs_match{1,2}(:,match_all_excl_plus_oris);
d6_tcs_match = d6_TCs.cellTCs_match{1,2}(:,match_all_excl_plus_oris);

nCells = size(d2_tcs_match,2);
nTrials = size(input.tGratingDirectionDeg,2);
nOff = input.nScansOff;
nOn = input.nScansOn;


d2_data_trial = permute(reshape(d2_tcs_match,[nOff + nOn nTrials nCells]),[1 3 2]);
d4_data_trial = permute(reshape(d4_tcs_match,[nOff + nOn nTrials nCells]),[1 3 2]);
d5_data_trial = permute(reshape(d5_tcs_match,[nOff + nOn nTrials nCells]),[1 3 2]);
d6_data_trial = permute(reshape(d6_tcs_match,[nOff + nOn nTrials nCells]),[1 3 2]);


%%
d2_data_f = mean(d2_data_trial(nOff-15:nOff,:,:),1);
d2_data_dfof = (d2_data_trial-d2_data_f)./d2_data_f;
d4_data_f = mean(d4_data_trial(nOff-15:nOff,:,:),1);
d4_data_dfof = (d4_data_trial-d4_data_f)./d4_data_f;
d5_data_f = mean(d5_data_trial(nOff-15:nOff,:,:),1);
d5_data_dfof = (d5_data_trial-d5_data_f)./d5_data_f;
d6_data_f = mean(d6_data_trial(nOff-15:nOff,:,:),1);
d6_data_dfof = (d6_data_trial-d6_data_f)./d6_data_f;


d2_data_dfof = mean(d2_data_dfof(nOff+1:nOn,:,:),1);
d2_data_dfof = squeeze(d2_data_dfof);

d2_mean_data_dfof = mean(d2_data_dfof,2);
d2_se_data_dfof = std(d2_data_dfof,[],2,'omitnan')./sqrt(size(d2_data_dfof,2));


d4_data_dfof = mean(d4_data_dfof(nOff+1:nOn,:,:),1);
d4_data_dfof = squeeze(d4_data_dfof);

d4_mean_data_dfof = mean(d4_data_dfof,2);
d4_se_data_dfof = std(d4_data_dfof,[],2,'omitnan')./sqrt(size(d4_data_dfof,2));

d5_data_dfof = mean(d5_data_dfof(nOff+1:nOn,:,:),1);
d5_data_dfof = squeeze(d5_data_dfof);

d5_mean_data_dfof = mean(d5_data_dfof,2);
d5_se_data_dfof = std(d5_data_dfof,[],2,'omitnan')./sqrt(size(d5_data_dfof,2));

d6_data_dfof = mean(d6_data_dfof(nOff+1:nOn,:,:),1);
d6_data_dfof = squeeze(d6_data_dfof);

d6_mean_data_dfof = mean(d6_data_dfof,2);
d6_se_data_dfof = std(d6_data_dfof,[],2,'omitnan')./sqrt(size(d6_data_dfof,2));


%%
% optional thing for looking at one cell across trials
figure;
subplot(2,2,1)
plot(d2_data_dfof(14,:))
subplot(2,2,2)
plot(d4_data_dfof(14,:))
subplot(2,2,3)
plot(d5_data_dfof(14,:))
subplot(2,2,4)
plot(d6_data_dfof(14,:))


%%
singlecell_means_across_sess = [d1_mean_data_dfof_45 d2_mean_data_dfof d4_mean_data_dfof d5_mean_data_dfof d6_mean_data_dfof d7_mean_data_dfof_45];
singlecell_sems_across_sess = [d1_se_data_dfof_45 d2_se_data_dfof d4_se_data_dfof d5_se_data_dfof d6_se_data_dfof d7_se_data_dfof_45];

%%
[n,n2] = subplotn(nCells);
grand_min = min(min(singlecell_means_across_sess)) - max(max(singlecell_sems_across_sess));
grand_max = max(max(singlecell_means_across_sess)) + max(max(singlecell_sems_across_sess));

figure;
for i = 1:nCells
    subplot(n, n2, i)
    errorbar(singlecell_means_across_sess(i,:), singlecell_sems_across_sess(i,:))
    xlim([0 7])
    ylim([grand_min grand_max])
%     title(['Cell #', match_all_excl(:,i)])
    title(match_all_excl(:,i))
    vline([2 5])

end

print(fullfile(newfnout, [mouse,  '_singlecell_resp.pdf']), '-dpdf', '-fillpage')

