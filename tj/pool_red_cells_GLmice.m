%clear everything
clear all
clear all global
clc
%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P\Arc_Multi_Day_Comparisons'; %folder to save files to
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P\Arc_Pooled_Multi_Day'; %folder to save files to
dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
lacz_list = [49 55 61]; % mice from exp list to use
enh_list = [67 73 76]; % mice from exp list to use

lacz_d_scores = [];
lacz_k_scores = [];
lacz_max_scores = [];
enh_d_scores = [];
enh_k_scores = [];
enh_max_scores = [];
%%
%lacz loop
for iexp = lacz_list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    lacz_dscore_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_d_scores']));
    lacz_dscore_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_d_scores_2']));
    lacz_dscore_ses_all = [lacz_dscore_ses1 lacz_dscore_ses2];
    lacz_d_scores = [lacz_d_scores lacz_dscore_ses_all];
    lacz_k_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_k_changes']));
    lacz_k_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_k_changes_2']));
    lacz_k_ses_all = [lacz_k_ses1 lacz_k_ses2];
    lacz_k_scores = [lacz_k_scores lacz_k_ses_all];
    lacz_max_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_max_changes']));
    lacz_max_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_max_changes_2']));
    lacz_max_ses_all = [lacz_max_ses1 lacz_max_ses2];
    lacz_max_scores = [lacz_max_scores lacz_max_ses_all];
end

%enh loop
for iexp = enh_list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    enh_dscore_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_d_scores']));
    if isfile(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_d_scores_2.mat']))
        enh_dscore_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_d_scores_2']));
    else
        enh_dscore_ses2 = [];
    end
    enh_dscore_ses_all = [enh_dscore_ses1 enh_dscore_ses2];
    enh_d_scores = [enh_d_scores enh_dscore_ses_all];
    enh_k_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_k_changes']));
    if isfile(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_k_changes_2.mat']))
        enh_k_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_k_changes_2']));
    else
        enh_k_ses2 = [];
    end
    enh_k_ses_all = [enh_k_ses1 enh_k_ses2];
    enh_k_scores = [enh_k_scores enh_k_ses_all];
    enh_max_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_max_changes']));
    if isfile(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_max_changes_2.mat']))
        enh_max_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_max_changes_2']));
    else
        enh_max_ses2 = [];
    end
    enh_max_ses_all = [enh_max_ses1 enh_max_ses2];
    enh_max_scores = [enh_max_scores enh_max_ses_all];
end
%%
%lacz
lacz_k_d1_d2_all = [];
lacz_k_d2_all = [];
lacz_k_d1_d3_all = [];
lacz_k_d3_all = [];
lacz_max_d1_d2_all = [];
lacz_max_d2_all = [];
lacz_max_d1_d3_all = [];
lacz_max_d3_all = [];

for idata = 1:length(lacz_k_scores)
    lacz_k_d1_d2 = lacz_k_scores(idata).k_d1_d2_match;
    lacz_k_d2 = lacz_k_scores(idata).k_d2_match;
    lacz_k_d1_d3 = lacz_k_scores(idata).k_d1_d3_match;
    lacz_k_d3 = lacz_k_scores(idata).k_d3_match;
    lacz_k_d1_d2_all = [lacz_k_d1_d2_all lacz_k_d1_d2];
    lacz_k_d2_all = [lacz_k_d2_all lacz_k_d2];
    lacz_k_d1_d3_all = [lacz_k_d1_d3_all lacz_k_d1_d3];
    lacz_k_d3_all = [lacz_k_d3_all lacz_k_d3];
end

for idata = 1:length(lacz_max_scores)
    lacz_max_d1_d2 = lacz_max_scores(idata).max_d1_d2_match;
    lacz_max_d2 = lacz_max_scores(idata).max_d2_match;
    lacz_max_d1_d3 = lacz_max_scores(idata).max_d1_d3_match;
    lacz_max_d3 = lacz_max_scores(idata).max_d3_match;
    lacz_max_d1_d2_all = [lacz_max_d1_d2_all lacz_max_d1_d2];
    lacz_max_d2_all = [lacz_max_d2_all lacz_max_d2];
    lacz_max_d1_d3_all = [lacz_max_d1_d3_all lacz_max_d1_d3];
    lacz_max_d3_all = [lacz_max_d3_all lacz_max_d3];
end

%k vals
figure('Position', [400 10 650 800]);
sgtitle('All LacZ Mice k Values')
subplot(2,1,1);
scatter(lacz_k_d1_d2_all, lacz_k_d2_all, 'filled', 'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
hold on
xlabel('Day 1 k Value');
ylabel('Day 2 k Value');
xlim([0,30]);
ylim([0,30]);
line = refline(1,0);
line.Color = 'k';
subplot(2,1,2);
scatter(lacz_k_d1_d3_all, lacz_k_d3_all, 'filled', 'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
xlabel('Day 1 k Value');
ylabel('Day 3 k Value');
xlim([0,30]);
ylim([0,30]);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
hold off
print(fullfile(realfnout, ['lacz mice', '_k_scatter.pdf']), '-dpdf', '-bestfit')

%k hist
figure('Position', [400 20 650 700]);
subplot(2,2,1);
histogram(lacz_k_d1_d2_all, 20);
hold on
subplot(2,2,2);
histogram(lacz_k_d2_all, 20);
subplot(2,2,3);
histogram(lacz_k_d1_d3_all, 20);
subplot(2,2,4);
histogram(lacz_k_d3_all, 20);
hold off
print(fullfile(realfnout, ['lacz mice', '_k_hist.pdf']), '-dpdf', '-bestfit')

%max vals
figure('Position', [400 10 650 800]);
sgtitle('All LacZ Mice Max dF/F Values')
subplot(2,1,1);
scatter(lacz_max_d1_d2_all, lacz_max_d2_all, 'filled', 'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
hold on
xlabel('Day 1 Max dF/F Value');
ylabel('Day 2 Max dF/F Value');
xlim([0,1]);
ylim([0,1]);
line = refline(1,0);
line.Color = 'k';
subplot(2,1,2);
scatter(lacz_max_d1_d3_all, lacz_max_d3_all, 'filled', 'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
xlabel('Day 1 Max dF/F Value');
ylabel('Day 3 Max dF/F Value');
xlim([0,1]);
ylim([0,1]);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
hold off
print(fullfile(realfnout, ['lacz mice', '_max_scatter.pdf']), '-dpdf', '-bestfit')

%max hist
figure('Position', [400 20 650 700]);
subplot(2,2,1);
histogram(lacz_max_d1_d2_all, 20);
hold on
subplot(2,2,2);
histogram(lacz_max_d2_all, 20);
subplot(2,2,3);
histogram(lacz_max_d1_d3_all, 20);
subplot(2,2,4);
histogram(lacz_max_d3_all, 20);
hold off
print(fullfile(realfnout, ['lacz mice', '_max_hist.pdf']), '-dpdf', '-bestfit')

%%
%enh
enh_k_d1_d2_all = [];
enh_k_d2_all = [];
enh_k_d1_d3_all = [];
enh_k_d3_all = [];
enh_max_d1_d2_all = [];
enh_max_d2_all = [];
enh_max_d1_d3_all = [];
enh_max_d3_all = [];

for idata = 1:length(enh_k_scores)
    enh_k_d1_d2 = enh_k_scores(idata).k_d1_d2_match;
    enh_k_d2 = enh_k_scores(idata).k_d2_match;
    enh_k_d1_d3 = enh_k_scores(idata).k_d1_d3_match;
    enh_k_d3 = enh_k_scores(idata).k_d3_match;
    enh_k_d1_d2_all = [enh_k_d1_d2_all enh_k_d1_d2];
    enh_k_d2_all = [enh_k_d2_all enh_k_d2];
    enh_k_d1_d3_all = [enh_k_d1_d3_all enh_k_d1_d3];
    enh_k_d3_all = [enh_k_d3_all enh_k_d3];
end

for idata = 1:length(enh_max_scores)
    enh_max_d1_d2 = enh_max_scores(idata).max_d1_d2_match;
    enh_max_d2 = enh_max_scores(idata).max_d2_match;
    enh_max_d1_d3 = enh_max_scores(idata).max_d1_d3_match;
    enh_max_d3 = enh_max_scores(idata).max_d3_match;
    enh_max_d1_d2_all = [enh_max_d1_d2_all enh_max_d1_d2];
    enh_max_d2_all = [enh_max_d2_all enh_max_d2];
    enh_max_d1_d3_all = [enh_max_d1_d3_all enh_max_d1_d3];
    enh_max_d3_all = [enh_max_d3_all enh_max_d3];
end

%k
figure('Position', [400 10 650 800]);
sgtitle('All Enhancer Mice k Values')
subplot(2,1,1);
scatter(enh_k_d1_d2_all, enh_k_d2_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'LineWidth', 0.8);
hold on
xlabel('Day 1 k Value');
ylabel('Day 2 k Value');
xlim([0,30]);
ylim([0,30]);
line = refline(1,0);
line.Color = 'k';
subplot(2,1,2);
scatter(enh_k_d1_d3_all, enh_k_d3_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'LineWidth', 0.8);
xlabel('Day 1 k Value');
ylabel('Day 3 k Value');
xlim([0,30]);
ylim([0,30]);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
hold off
print(fullfile(realfnout, ['enh mice', '_k_scatter.pdf']), '-dpdf', '-bestfit')

%k hist
figure('Position', [400 20 650 700]);
subplot(2,2,1);
histogram(enh_k_d1_d2_all, 20);
hold on
subplot(2,2,2);
histogram(enh_k_d2_all, 20);
subplot(2,2,3);
histogram(enh_k_d1_d3_all, 20);
subplot(2,2,4);
histogram(enh_k_d3_all, 20);
hold off
print(fullfile(realfnout, ['enh mice', '_k_hist.pdf']), '-dpdf', '-bestfit')

%max
figure('Position', [400 10 650 800]);
sgtitle('All Enhancer Mice dF/F Values')
subplot(2,1,1);
scatter(enh_max_d1_d2_all, enh_max_d2_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'LineWidth', 0.8);
hold on
xlabel('Day 1 Max dF/F Value');
ylabel('Day 2 Max dF/F Value');
xlim([0,1]);
ylim([0,1]);
line = refline(1,0);
line.Color = 'k';
subplot(2,1,2);
scatter(enh_max_d1_d3_all, enh_max_d3_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'LineWidth', 0.8);
xlabel('Day 1 Max dF/F Value');
ylabel('Day 3 Max dF/F Value');
xlim([0,1]);
ylim([0,1]);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
hold off
print(fullfile(realfnout, ['enh mice', '_max_scatter.pdf']), '-dpdf', '-bestfit')

%max hist
figure('Position', [400 20 650 700]);
subplot(2,2,1);
histogram(enh_max_d1_d2_all, 20);
hold on
subplot(2,2,2);
histogram(enh_max_d2_all, 20);
subplot(2,2,3);
histogram(enh_max_d1_d3_all, 20);
subplot(2,2,4);
histogram(enh_max_d3_all, 20);
hold off
print(fullfile(realfnout, ['enh mice', '_max_hist.pdf']), '-dpdf', '-bestfit')

%%
%prefori changes
lacz_dscore_prefori_d1_d2_all = [];
lacz_dscore_prefori_d1_d3_all = [];
enh_dscore_prefori_d1_d2_all = [];
enh_dscore_prefori_d1_d3_all = [];

for idata = 1:length(lacz_d_scores)
    lacz_dscore_prefori_d1_d2 = lacz_d_scores(idata).d_score_prefori_d1_d2;
    lacz_dscore_prefori_d1_d3 = lacz_d_scores(idata).d_score_prefori_d1_d3;
    lacz_dscore_prefori_d1_d2_all = [lacz_dscore_prefori_d1_d2_all lacz_dscore_prefori_d1_d2];
    lacz_dscore_prefori_d1_d3_all = [lacz_dscore_prefori_d1_d3_all lacz_dscore_prefori_d1_d3];
end

for idata = 1:length(enh_d_scores)
    enh_dscore_prefori_d1_d2 = enh_d_scores(idata).d_score_prefori_d1_d2;
    enh_dscore_prefori_d1_d3 = enh_d_scores(idata).d_score_prefori_d1_d3;
    enh_dscore_prefori_d1_d2_all = [enh_dscore_prefori_d1_d2_all enh_dscore_prefori_d1_d2];
    enh_dscore_prefori_d1_d3_all = [enh_dscore_prefori_d1_d3_all enh_dscore_prefori_d1_d3];
end

figure;
m = cdfplot(enh_dscore_prefori_d1_d2_all);
hold on
n = cdfplot(enh_dscore_prefori_d1_d3_all);
p = cdfplot(lacz_dscore_prefori_d1_d2_all);
q = cdfplot(lacz_dscore_prefori_d1_d3_all);
set(m, 'LineStyle', '-', 'Color', '#00841a', 'LineWidth', 1.0);
set(n, 'LineStyle', '--', 'Color', '#00841a', 'LineWidth', 1.0);
set(p, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 1.0);
set(q, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.0);
legend(['Enhancer D1-D2 (n = ', num2str(length(enh_dscore_prefori_d1_d2_all)), ')'],...
    ['Enhancer D1-D3 (n = ', num2str(length(enh_dscore_prefori_d1_d3_all)), ')'],...
    ['LacZ D1-D2 (n = ', num2str(length(lacz_dscore_prefori_d1_d2_all)), ')'],...
    ['LacZ D1-D3 (n = ', num2str(length(lacz_dscore_prefori_d1_d3_all)), ')'], 'Location', 'Best')
%title('Change in Pref Ori Values')
title('')
xlabel('Change in Pref Ori')
xlim([0 90])
hold off

print(fullfile(realfnout, ['Grace_all_mice', '_change_pref_cdf.pdf']), '-dpdf', '-bestfit')

%%
%individual mice d1 v d2 enhancer and lacz
mouse1_enh_dscores_prefori_d1_d2 = [enh_d_scores(1).d_score_prefori_d1_d2 enh_d_scores(2).d_score_prefori_d1_d2]
mouse2_enh_dscores_prefori_d1_d2 = [enh_d_scores(3).d_score_prefori_d1_d2]
mouse3_enh_dscores_prefori_d1_d2 = [enh_d_scores(4).d_score_prefori_d1_d2 enh_d_scores(5).d_score_prefori_d1_d2]

mouse1_lacz_dscores_prefori_d1_d2 = [lacz_d_scores(1).d_score_prefori_d1_d2 lacz_d_scores(2).d_score_prefori_d1_d2]
mouse2_lacz_dscores_prefori_d1_d2 = [lacz_d_scores(3).d_score_prefori_d1_d2 lacz_d_scores(4).d_score_prefori_d1_d2]
mouse3_lacz_dscores_prefori_d1_d2 = [lacz_d_scores(5).d_score_prefori_d1_d2 lacz_d_scores(6).d_score_prefori_d1_d2]

figure; 
m1_enh_d1_d2 = cdfplot(mouse1_enh_dscores_prefori_d1_d2);
hold on
m2_enh_d1_d2 = cdfplot(mouse2_enh_dscores_prefori_d1_d2);
m3_enh_d1_d2 = cdfplot(mouse3_enh_dscores_prefori_d1_d2);

set(m1_enh_d1_d2, 'Color', '#00841a', 'LineWidth', 1.0);
set(m2_enh_d1_d2, 'Color', '#00841a', 'LineWidth', 1.0);
set(m3_enh_d1_d2, 'Color', '#00841a', 'LineWidth', 1.0);

m1_lacz_d1_d2 = cdfplot(mouse1_lacz_dscores_prefori_d1_d2);
m2_lacz_d1_d2 = cdfplot(mouse2_lacz_dscores_prefori_d1_d2);
m3_lacz_d1_d2 = cdfplot(mouse3_lacz_dscores_prefori_d1_d2);

set(m1_lacz_d1_d2, 'Color', 'r', 'LineWidth', 1.0);
set(m2_lacz_d1_d2, 'Color', 'r', 'LineWidth', 1.0);
set(m3_lacz_d1_d2, 'Color', 'r', 'LineWidth', 1.0);

title('Day 1 vs. Day 2')
xlabel('Change in Pref Ori')
xlim([0 90])
legend([m1_enh_d1_d2(1), m1_lacz_d1_d2(1)], 'Arc Enhancer', ' LacZ')
hold off

print(fullfile(realfnout, ['individual_Grace_mice', '_d1_d2_change_pref_cdf.pdf']), '-dpdf', '-bestfit')


%individual mice d1 v d3 enhancer and lacz
mouse1_enh_dscores_prefori_d1_d3 = [enh_d_scores(1).d_score_prefori_d1_d3 enh_d_scores(2).d_score_prefori_d1_d3]
mouse2_enh_dscores_prefori_d1_d3 = [enh_d_scores(3).d_score_prefori_d1_d3]
mouse3_enh_dscores_prefori_d1_d3 = [enh_d_scores(4).d_score_prefori_d1_d3 enh_d_scores(5).d_score_prefori_d1_d3]

mouse1_lacz_dscores_prefori_d1_d3 = [lacz_d_scores(1).d_score_prefori_d1_d3 lacz_d_scores(2).d_score_prefori_d1_d3]
mouse2_lacz_dscores_prefori_d1_d3 = [lacz_d_scores(3).d_score_prefori_d1_d3 lacz_d_scores(4).d_score_prefori_d1_d3]
mouse3_lacz_dscores_prefori_d1_d3 = [lacz_d_scores(5).d_score_prefori_d1_d3 lacz_d_scores(6).d_score_prefori_d1_d3]

figure; 
m1_enh_d1_d3 = cdfplot(mouse1_enh_dscores_prefori_d1_d3);
hold on
m2_enh_d1_d3 = cdfplot(mouse2_enh_dscores_prefori_d1_d3);
m3_enh_d1_d3 = cdfplot(mouse3_enh_dscores_prefori_d1_d3);

set(m1_enh_d1_d3, 'Color', '#00841a', 'LineWidth', 1.0);
set(m2_enh_d1_d3, 'Color', '#00841a', 'LineWidth', 1.0);
set(m3_enh_d1_d3, 'Color', '#00841a', 'LineWidth', 1.0);

m1_lacz_d1_d3 = cdfplot(mouse1_lacz_dscores_prefori_d1_d3);
m2_lacz_d1_d3 = cdfplot(mouse2_lacz_dscores_prefori_d1_d3);
m3_lacz_d1_d3 = cdfplot(mouse3_lacz_dscores_prefori_d1_d3);

set(m1_lacz_d1_d3, 'Color', 'r', 'LineWidth', 1.0);
set(m2_lacz_d1_d3, 'Color', 'r', 'LineWidth', 1.0);
set(m3_lacz_d1_d3, 'Color', 'r', 'LineWidth', 1.0);

title('Day 1 vs. Day 3')
xlabel('Change in Pref Ori')
xlim([0 90])
legend([m1_enh_d1_d3(1), m1_lacz_d1_d3(1)], 'Arc Enhancer', ' LacZ')
hold off

print(fullfile(realfnout, ['individual_Grace_mice', '_d1_d3_change_pref_cdf.pdf']), '-dpdf', '-bestfit')

%%
%k changes
lacz_dscore_k_d1_d2_all = [];
lacz_dscore_k_d1_d3_all = [];
enh_dscore_k_d1_d2_all = [];
enh_dscore_k_d1_d3_all = [];
for idata = 1:length(lacz_d_scores)
    lacz_dscore_k_d1_d2 = lacz_d_scores(idata).dscore_k_d1_d2_match;
    lacz_dscore_k_d1_d3 = lacz_d_scores(idata).dscore_k_d1_d3_match;
    lacz_dscore_k_d1_d2_all = [lacz_dscore_k_d1_d2_all lacz_dscore_k_d1_d2];
    lacz_dscore_k_d1_d3_all = [lacz_dscore_k_d1_d3_all lacz_dscore_k_d1_d3];
end

for idata = 1:length(enh_d_scores)
    enh_dscore_k_d1_d2 = enh_d_scores(idata).dscore_k_d1_d2_match;
    enh_dscore_k_d1_d3 = enh_d_scores(idata).dscore_k_d1_d3_match;
    enh_dscore_k_d1_d2_all = [enh_dscore_k_d1_d2_all enh_dscore_k_d1_d2];
    enh_dscore_k_d1_d3_all = [enh_dscore_k_d1_d3_all enh_dscore_k_d1_d3];
end

figure; 
m = cdfplot(enh_dscore_k_d1_d2_all);
hold on
n = cdfplot(enh_dscore_k_d1_d3_all);
p = cdfplot(lacz_dscore_k_d1_d2_all);
q = cdfplot(lacz_dscore_k_d1_d3_all);
set(m, 'LineStyle', '-', 'Color', 'b');
set(n, 'LineStyle', '--', 'Color', 'b');
set(p, 'LineStyle', '-', 'Color', 'r');
set(q, 'LineStyle', '--', 'Color', 'r');
legend(['Enhancer D1-D2 (n = ', num2str(length(enh_dscore_prefori_d1_d2_all)), ')'],...
    ['Enhancer D1-D3 (n = ', num2str(length(enh_dscore_prefori_d1_d3_all)), ')'],...
    ['LacZ D1-D2 (n = ', num2str(length(lacz_dscore_prefori_d1_d2_all)), ')'],...
    ['LacZ D1-D3 (n = ', num2str(length(lacz_dscore_prefori_d1_d3_all)), ')'], 'Location', 'Best')
title('Change in Sharpness of Tuning (K) Values')
xlabel('Change in K')
hold off

print(fullfile(realfnout, ['Grace_all_mice', '_change_k_cdf.pdf']), '-dpdf', '-bestfit')

%%
%max df/f changes
lacz_dscore_max_d1_d2_all = [];
lacz_dscore_max_d1_d3_all = [];
enh_dscore_max_d1_d2_all = [];
enh_dscore_max_d1_d3_all = [];

for idata = 1:length(lacz_d_scores)
    lacz_dscore_max_d1_d2 = lacz_d_scores(idata).dscore_max_d1_d2_match;
    lacz_dscore_max_d1_d3 = lacz_d_scores(idata).dscore_max_d1_d3_match;
    lacz_dscore_max_d1_d2_all = [lacz_dscore_max_d1_d2_all lacz_dscore_max_d1_d2];
    lacz_dscore_max_d1_d3_all = [lacz_dscore_max_d1_d3_all lacz_dscore_max_d1_d3];
end

for idata = 1:length(enh_d_scores)
    enh_dscore_max_d1_d2 = enh_d_scores(idata).dscore_max_d1_d2_match;
    enh_dscore_max_d1_d3 = enh_d_scores(idata).dscore_max_d1_d3_match;
    enh_dscore_max_d1_d2_all = [enh_dscore_max_d1_d2_all enh_dscore_max_d1_d2];
    enh_dscore_max_d1_d3_all = [enh_dscore_max_d1_d3_all enh_dscore_max_d1_d3];
end

figure; 
m = cdfplot(enh_dscore_max_d1_d2_all);
hold on
n = cdfplot(enh_dscore_max_d1_d3_all);
p = cdfplot(lacz_dscore_max_d1_d2_all);
q = cdfplot(lacz_dscore_max_d1_d3_all);
set(m, 'LineStyle', '-', 'Color', 'b');
set(n, 'LineStyle', '--', 'Color', 'b');
set(p, 'LineStyle', '-', 'Color', 'r');
set(q, 'LineStyle', '--', 'Color', 'r');
legend(['Enhancer D1-D2 (n = ', num2str(length(enh_dscore_prefori_d1_d2_all)), ')'],...
    ['Enhancer D1-D3 (n = ', num2str(length(enh_dscore_prefori_d1_d3_all)), ')'],...
    ['LacZ D1-D2 (n = ', num2str(length(lacz_dscore_prefori_d1_d2_all)), ')'],...
    ['LacZ D1-D3 (n = ', num2str(length(lacz_dscore_prefori_d1_d3_all)), ')'], 'Location', 'Best')
title('Change in Max dF/F Values')
xlabel('Change in Max dF/F')
hold off

print(fullfile(realfnout, ['Grace_all_mice', '_change_max_cdf.pdf']), '-dpdf', '-bestfit')
