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
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\layer4'; %folder to save files to
dataset = 'exp_list_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%
%below is the legend for which numbers correspond to which imaging
%timepoints in this study
%img_point: 1 = weeek 1/baseline, 2 = week 2, 3 = week 3, 4 = week 4, 5 =
%week 5


%%

%we start with loading the proper information based on the mouse's info in
%the experiment list

d1 = 52; 
d2 = 54; 
d3 = 56;

mouse = expt(d1).mouse; 
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


d1_new_pref = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'new_pref.mat']));
d2_new_pref = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'new_pref.mat']));
d3_new_pref = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'new_pref.mat']));


%% matching indices

tuned_d1 = d1_ori.ind_theta90;
tuned_d2 = d2_ori.ind_theta90;
tuned_d3 = d3_ori.ind_theta90;


tuned_d1_d2 = intersect(tuned_d1, tuned_d2);
tuned_d1_d3 = intersect(tuned_d1, tuned_d3);
tuned_d2_d3 = intersect(tuned_d2, tuned_d3);

tuned_all = unique([tuned_d1 tuned_d2 tuned_d3]);

match_d2 = find([d2_matches.cellImageAlign.pass]); 
match_d3 = find([d3_matches.cellImageAlign.pass]); 

match_d1_d2 = intersect(tuned_d1_d2, match_d2);
match_d1_d3 = intersect(tuned_d1_d3, match_d3);
match_d2_d3 = intersect(intersect(match_d2, match_d3), tuned_d2_d3);

match_all = intersect(intersect(match_d2,match_d3), tuned_all);

%%
%within session pref ori change
d1_within_prefori_1 = changeprefto90TJ(d1_new_pref.prefOri_1);
d1_within_prefori_2 = changeprefto90TJ(d1_new_pref.prefOri_2);
d2_within_prefori_1 = changeprefto90TJ(d2_new_pref.prefOri_1);
d2_within_prefori_2 = changeprefto90TJ(d2_new_pref.prefOri_2);
d3_within_prefori_1 = changeprefto90TJ(d3_new_pref.prefOri_1);
d3_within_prefori_2 = changeprefto90TJ(d3_new_pref.prefOri_2);

d_score_prefori_within_d1 = abs(d1_within_prefori_1-d1_within_prefori_2);
d_score_prefori_within_d1 = d_score_prefori_within_d1(tuned_d1);
d_score_prefori_within_d2 = abs(d2_within_prefori_1-d2_within_prefori_2);
d_score_prefori_within_d2 = d_score_prefori_within_d2(tuned_d2);
d_score_prefori_within_d3 = abs(d3_within_prefori_1-d3_within_prefori_2);
d_score_prefori_within_d3 = d_score_prefori_within_d3(tuned_d3);


fig = figure;
sgtitle('Pref Ori Changes within Session')
% a=subplot(2,2,1)
h=cdfplot(d_score_prefori_within_d1);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d_score_prefori_within_d2);
set(j, 'LineStyle', '-', 'Color', [.4 .4 .4], 'LineWidth',1);
l=cdfplot(d_score_prefori_within_d3);
set(l, 'LineStyle', '-', 'Color', [.8 .8 .8], 'LineWidth',1);
legend(['Session 1'], ['Session 2'], ['Session 3'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori w/i Session'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(newfnout, [mouse,  '_prefori_change_within_sess.pdf']), '-dpdf', '-bestfit')

%%

d1_prefori = d1_ori.prefOri(1,:);
dscore_prefori_1 = double(d1_prefori>90);
dscore_prefori_1(dscore_prefori_1>0) = 180;
dscore_prefori_1 = abs(dscore_prefori_1-d1_prefori);

d2_prefori = d2_ori.prefOri(1,:);
dscore_prefori_2 = double(d2_prefori>90);
dscore_prefori_2(dscore_prefori_2>0) = 180;
dscore_prefori_2 = abs(dscore_prefori_2-d2_prefori);

d3_prefori = d3_ori.prefOri(1,:);
dscore_prefori_3 = double(d3_prefori>90);
dscore_prefori_3(dscore_prefori_3>0) = 180;
dscore_prefori_3 = abs(dscore_prefori_3-d3_prefori);


d_score_prefori_d1_d2 = abs(dscore_prefori_1-dscore_prefori_2);
d_score_prefori_d1_d2 = d_score_prefori_d1_d2(match_d1_d2);

d_score_prefori_d2_d3 = abs(dscore_prefori_2-dscore_prefori_3);
d_score_prefori_d2_d3 = d_score_prefori_d2_d3(match_d2_d3);

d_score_prefori_d1_d3 = abs(dscore_prefori_1-dscore_prefori_3);
d_score_prefori_d1_d3 = d_score_prefori_d1_d3(match_d1_d3);


%%

fig = figure;
sgtitle('Pref Ori Changes from Day 1')
% a=subplot(2,2,1)
h=cdfplot(d_score_prefori_d1_d2);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d_score_prefori_d1_d3);
set(j, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',1);
legend(['Session 1 to 2'], ['Session 1 to 3'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(newfnout, [mouse,  '_prefori_change_fromd1.pdf']), '-dpdf', '-bestfit')

%%
d_score_prefori_1sess = [d_score_prefori_d1_d2 d_score_prefori_d2_d3];
d_score_prefori_2sess = [d_score_prefori_d1_d3];


%%
fig = figure;
sgtitle('Pref Ori Changes All Sessions')
% a=subplot(2,2,1)
h=cdfplot(d_score_prefori_1sess);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d_score_prefori_2sess);
set(j, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',1);
legend(['1 session'], ['2 sessions'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(newfnout, [mouse,  '_prefori_change_all.pdf']), '-dpdf', '-bestfit')

%%
%% NEW matching indices

new_tuned_d1 = d1_new_pref.new_tune;
new_tuned_d2 = d2_new_pref.new_tune;
new_tuned_d3 = d3_new_pref.new_tune;

new_tuned_d1_d2 = intersect(new_tuned_d1, new_tuned_d2);
new_tuned_d1_d3 = intersect(new_tuned_d1, new_tuned_d3);
new_tuned_d2_d3 = intersect(new_tuned_d2, new_tuned_d3);

new_tuned_all = unique([new_tuned_d1 new_tuned_d2 new_tuned_d3]);

new_match_d2 = find([d2_matches.cellImageAlign.pass]); 
new_match_d3 = find([d3_matches.cellImageAlign.pass]); 

new_match_d1_d2 = intersect(new_tuned_d1_d2, new_match_d2);
new_match_d1_d3 = intersect(new_tuned_d1_d3, new_match_d3);
new_match_d2_d3 = intersect(intersect(new_match_d2, new_match_d3), new_tuned_d2_d3);

new_match_all = intersect(intersect(new_match_d2,new_match_d3),new_tuned_all);

%%
%within session theta change
d1_within_new_prefori_1 = changeprefto90TJ(d1_new_pref.theta_1);
d1_within_new_prefori_2 = changeprefto90TJ(d1_new_pref.theta_2);
d2_within_new_prefori_1 = changeprefto90TJ(d2_new_pref.theta_1);
d2_within_new_prefori_2 = changeprefto90TJ(d2_new_pref.theta_2);
d3_within_new_prefori_1 = changeprefto90TJ(d3_new_pref.theta_1);
d3_within_new_prefori_2 = changeprefto90TJ(d3_new_pref.theta_2);

d_score_new_prefori_within_d1 = abs(d1_within_new_prefori_1-d1_within_new_prefori_2);
d_score_new_prefori_within_d1 = d_score_new_prefori_within_d1(tuned_d1);
d_score_new_prefori_within_d2 = abs(d2_within_new_prefori_1-d2_within_new_prefori_2);
d_score_new_prefori_within_d2 = d_score_new_prefori_within_d2(tuned_d2);
d_score_new_prefori_within_d3 = abs(d3_within_new_prefori_1-d3_within_new_prefori_2);
d_score_new_prefori_within_d3 = d_score_new_prefori_within_d3(tuned_d3);

fig = figure;
sgtitle('Theta Changes within Session')
% a=subplot(2,2,1)
h=cdfplot(d_score_new_prefori_within_d1);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d_score_new_prefori_within_d2);
set(j, 'LineStyle', '-', 'Color', [.2 .2 .2], 'LineWidth',1);
l=cdfplot(d_score_new_prefori_within_d3);
set(l, 'LineStyle', '-', 'Color', [.4 .4 .4], 'LineWidth',1);
legend(['Session 1'], ['Session 2'], ['Session 3'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori w/i Session'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(newfnout, [mouse,  '_new_prefori_change_within_sess.pdf']), '-dpdf', '-bestfit')



%%
d1_new_prefori = d1_new_pref.theta;
dscore_new_prefori_1 = double(d1_new_prefori>90);
dscore_new_prefori_1(dscore_new_prefori_1>0) = 180;
dscore_new_prefori_1 = abs(dscore_new_prefori_1-d1_new_prefori);

d2_new_prefori = d2_new_pref.theta;
dscore_new_prefori_2 = double(d2_new_prefori>90);
dscore_new_prefori_2(dscore_new_prefori_2>0) = 180;
dscore_new_prefori_2 = abs(dscore_new_prefori_2-d2_new_prefori);

d3_new_prefori = d3_new_pref.theta;
dscore_new_prefori_3 = double(d3_new_prefori>90);
dscore_new_prefori_3(dscore_new_prefori_3>0) = 180;
dscore_new_prefori_3 = abs(dscore_new_prefori_3-d3_new_prefori);


d_score_new_prefori_d1_d2 = abs(dscore_new_prefori_1-dscore_new_prefori_2);
d_score_new_prefori_d1_d2 = d_score_new_prefori_d1_d2(new_match_d1_d2);

d_score_new_prefori_d2_d3 = abs(dscore_new_prefori_2-dscore_new_prefori_3);
d_score_new_prefori_d2_d3 = d_score_new_prefori_d2_d3(new_match_d2_d3);

d_score_new_prefori_d1_d3 = abs(dscore_new_prefori_1-dscore_new_prefori_3);
d_score_new_prefori_d1_d3 = d_score_new_prefori_d1_d3(new_match_d1_d3);

%%

fig = figure;
sgtitle('Vec Sum Pref Ori Changes from Day 1')
% a=subplot(2,2,1)
h=cdfplot(d_score_new_prefori_d1_d2);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d_score_new_prefori_d1_d3);
set(j, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',1);
legend(['Session 1 to 2'], ['Session 1 to 3'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(newfnout, [mouse,  '_vecsum_prefori_change_fromd1.pdf']), '-dpdf', '-bestfit')

%%

d_score_new_prefori_1sess = [d_score_new_prefori_d1_d2 d_score_new_prefori_d2_d3];
d_score_new_prefori_2sess = [d_score_new_prefori_d1_d3];

%%

fig = figure;
sgtitle('Vec Sum Pref Ori Changes All Sessions')
% a=subplot(2,2,1)
h=cdfplot(d_score_new_prefori_1sess);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d_score_new_prefori_2sess);
set(j, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',1);
legend(['1 session'], ['2 sessions'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(newfnout, [mouse,  '_vecsum_prefori_change_all.pdf']), '-dpdf', '-bestfit')


%%

wk1_overlap = length(intersect(tuned_d1, new_tuned_d1)) / ((length(intersect(tuned_d1, new_tuned_d1))) + (length(setdiff(new_tuned_d1, tuned_d1))) + (length(setdiff(tuned_d1, new_tuned_d1))));
wk1_vecsum =  length(setdiff(new_tuned_d1, tuned_d1)) / ((length(intersect(tuned_d1, new_tuned_d1))) + (length(setdiff(new_tuned_d1, tuned_d1))) + (length(setdiff(tuned_d1, new_tuned_d1))));
wk1_tradition =  length(setdiff(tuned_d1, new_tuned_d1)) / ((length(intersect(tuned_d1, new_tuned_d1))) + (length(setdiff(new_tuned_d1, tuned_d1))) + (length(setdiff(tuned_d1, new_tuned_d1))));

wk2_overlap = length(intersect(tuned_d2, new_tuned_d2)) / ((length(intersect(tuned_d2, new_tuned_d2))) + (length(setdiff(new_tuned_d2, tuned_d2))) + (length(setdiff(tuned_d2, new_tuned_d2))));
wk2_vecsum =  length(setdiff(new_tuned_d2, tuned_d2)) / ((length(intersect(tuned_d2, new_tuned_d2))) + (length(setdiff(new_tuned_d2, tuned_d2))) + (length(setdiff(tuned_d2, new_tuned_d2))));
wk2_tradition =  length(setdiff(tuned_d2, new_tuned_d2)) / ((length(intersect(tuned_d2, new_tuned_d2))) + (length(setdiff(new_tuned_d2, tuned_d2))) + (length(setdiff(tuned_d2, new_tuned_d2))));

wk3_overlap = length(intersect(tuned_d3, new_tuned_d3)) / ((length(intersect(tuned_d3, new_tuned_d3))) + (length(setdiff(new_tuned_d3, tuned_d3))) + (length(setdiff(tuned_d3, new_tuned_d3))));
wk3_vecsum =  length(setdiff(new_tuned_d3, tuned_d3)) / ((length(intersect(tuned_d3, new_tuned_d3))) + (length(setdiff(new_tuned_d3, tuned_d3))) + (length(setdiff(tuned_d3, new_tuned_d3))));
wk3_tradition =  length(setdiff(tuned_d3, new_tuned_d3)) / ((length(intersect(tuned_d3, new_tuned_d3))) + (length(setdiff(new_tuned_d3, tuned_d3))) + (length(setdiff(tuned_d3, new_tuned_d3))));


%%

% figure;
% x_lab = categorical(["Week 1", "Week 2", "Week 3", "Week 4", "Week 5"]);
% bar(x_lab, [wk1_overlap*100, wk2_overlap*100, wk3_overlap*100, wk4_overlap*100, wk5_overlap*100]);
% ylabel("% of Cells");
% title(["Overlap of Tuning for Traditional and Vector Sum"])

figure;
x_lab = categorical(["Week 1", "Week 2", "Week 3"]);
y = [wk1_overlap*100, wk1_vecsum*100, wk1_tradition*100; wk2_overlap*100, wk2_vecsum*100, wk2_tradition*100; wk3_overlap*100, wk3_vecsum*100, ...
    wk3_tradition*100];
b = bar(x_lab, y, 'stacked','FaceColor','flat');
ylabel("% of Cells");
ylim([0 100])
legend(['Overlap'], ['Vector Sum Only'], ['Traditional Only'], 'Location', 'southeast')
title(["Overlap of Tuning for Traditional and Vector Sum"])

print(fullfile(newfnout, [mouse,  '_tuning_overlap.pdf']), '-dpdf', '-bestfit')

% for k = 1:size(y,2)
%     b(k).CData = k;
% end

%%

% figure;
% % subplot(3,2,1)
% h=scatter(dscore_prefori_1(intersect(tuned_d1,new_tuned_d1)), dscore_new_prefori_1(intersect(tuned_d1,new_tuned_d1)));
% set(h, "MarkerEdgeColor", [0 0 0]);
% hold on
% j=scatter(dscore_prefori_2(intersect(tuned_d2,new_tuned_d2)), dscore_new_prefori_2(intersect(tuned_d2,new_tuned_d2)));
% set(j, "MarkerEdgeColor", [.2 .2 .2]);
% l=scatter(dscore_prefori_3(intersect(tuned_d3,new_tuned_d3)), dscore_new_prefori_3(intersect(tuned_d3,new_tuned_d3)));
% set(l, "MarkerEdgeColor", [.4 .4 .4]);
% m=scatter(dscore_prefori_4(intersect(tuned_d4,new_tuned_d4)), dscore_new_prefori_4(intersect(tuned_d4,new_tuned_d4)));
% set(m, "MarkerEdgeColor", [.6 .6 .6]);
% o=scatter(dscore_prefori_5(intersect(tuned_d5,new_tuned_d5)), dscore_new_prefori_5(intersect(tuned_d5,new_tuned_d5)));
% set(o, "MarkerEdgeColor", [.8 .8 .8]);
% xlim([0 90])
% xlabel(['Change in Pref Ori'])
% ylabel(['% of cells'])
% title([])

fig = figure;
subplot(2,2,1)
h=scatter(dscore_prefori_1(intersect(tuned_d1,new_tuned_d1)), dscore_new_prefori_1(intersect(tuned_d1,new_tuned_d1)));
title('Week 1')
hold on
subplot(2,2,2)
j=scatter(dscore_prefori_2(intersect(tuned_d2,new_tuned_d2)), dscore_new_prefori_2(intersect(tuned_d2,new_tuned_d2)));
title('Week 2')
subplot(2,2,3.5)
l=scatter(dscore_prefori_3(intersect(tuned_d3,new_tuned_d3)), dscore_new_prefori_3(intersect(tuned_d3,new_tuned_d3)));
title('Week 3')
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Vec Sum Pref Ori');
xlabel(han,'Traditional Pref Ori');
title(han,'');

print(fullfile(newfnout, [mouse,  '_prefori_overlap.pdf']), '-dpdf', '-bestfit')


%%

%%

a = d1_fits.avgResponseEaOri(match_all,:);
b = d2_fits.avgResponseEaOri(match_all,:);
c = d3_fits.avgResponseEaOri(match_all,:);


%%
a_cells_corr = corr(a.');
b_cells_corr = corr(b.');
c_cells_corr = corr(c.');

%%
a_idx = eye(size(a_cells_corr));
a_Y = (1-a_idx).*a_cells_corr;
a_Z = a_cells_corr(~a_idx);

b_idx = eye(size(b_cells_corr));
b_Y = (1-b_idx).*b_cells_corr;
b_Z = b_cells_corr(~b_idx);

c_idx = eye(size(c_cells_corr));
c_Y = (1-c_idx).*c_cells_corr;
c_Z = c_cells_corr(~c_idx);

%%
abcde_mat = [a_Z b_Z c_Z];
abcde_corrs = corr(abcde_mat);

%%
figure;
scatter(1,[abcde_corrs(1,2), abcde_corrs(2,3)])
hold on
scatter(2, [abcde_corrs(1,3)])
xlim([0 5])
ylim([0 1])
xlabel('Days Apart')
ylabel('Correlation')
xticklabels({'','','7','','14','','21','','28'})

print(fullfile(newfnout, [mouse,  '_tunecurvecorr_bysess.pdf']), '-dpdf', '-bestfit')

%%
figure;
errorbar(1,mean([abcde_corrs(1,2), abcde_corrs(2,3)]), std([abcde_corrs(1,2), abcde_corrs(2,3)]), 'o')
hold on
errorbar(2, mean([abcde_corrs(1,3)]), std([abcde_corrs(1,3)]), 'o')
xlim([0 5])
ylim([0 1])
xlabel('Days Apart')
ylabel('Correlation')
xticklabels({'','','7','','14','','21','','28'})

print(fullfile(newfnout, [mouse,  '_tunecurvecorr_all.pdf']), '-dpdf', '-bestfit')

%%
%save all of the change scores from above for future use
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_d_scores.mat']), 'd_score_prefori_d1_d2', 'd_score_prefori_d2_d3', ...
    'd_score_prefori_d1_d3')

%save pref ori info
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_ori_changes.mat']), 'd1_prefori', 'd2_prefori', 'd3_prefori')

%save tuned and matched cell IDs
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_id_matches.mat']), 'match_all')


