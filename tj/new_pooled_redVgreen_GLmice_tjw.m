%CLEAR EVERYTHING
clear all
clear all global
clc
close all



%%
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_Arc_greenVred'; %folder to save files to
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_pooled_Arc_greenVred_GLmice'; %folder to save files to
dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
eval(dataset); %load dataset


%%
%% ARC PROMOTER MICE
%%
arc_green_d1_pref_all = [];
arc_green_d2_pref_all = [];
arc_green_d3_pref_all = [];
arc_green_d1_d2_pref_d_all = [];
arc_green_d1_d3_pref_d_all = [];
arc_green_d1_k_all = [];
arc_green_d2_k_all = [];
arc_green_d3_k_all = [];
arc_green_d1_max_all = [];
arc_green_d2_max_all = [];
arc_green_d3_max_all = [];

arc_green_d1_reliability_all = [];
arc_green_d2_reliability_all = [];
arc_green_d3_reliability_all = [];

arc_green_top_d1_d2_pref_d_all = [];
arc_green_top_d1_d3_pref_d_all = [];

arc_green_bot_d1_d2_pref_d_all = [];
arc_green_bot_d1_d3_pref_d_all = [];

arc_red_d1_pref_all = [];
arc_red_d2_pref_all = [];
arc_red_d3_pref_all = [];
arc_red_d1_d2_pref_d_all = [];
arc_red_d1_d3_pref_d_all = [];
arc_red_d1_k_all = [];
arc_red_d2_k_all = [];
arc_red_d3_k_all = [];
arc_red_d1_max_all = [];
arc_red_d2_max_all = [];
arc_red_d3_max_all = [];

arc_red_d1_reliability_all = [];
arc_red_d2_reliability_all = [];
arc_red_d3_reliability_all = [];

arc_red_top_d1_d2_pref_d_all = [];
arc_red_top_d1_d3_pref_d_all = [];

arc_red_bot_d1_d2_pref_d_all = [];
arc_red_bot_d1_d3_pref_d_all = [];



arc_prefori_scores = [];
arc_pref_dscores = [];
arc_k_max_scores = [];
arc_median_pref_dscores = [];
arc_reliability_scores = [];
arc_wide_sharp_pref_dscores = [];


list = [67 73 76];
for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};

    arc_prefori_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes_2.mat']));
        arc_prefori_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes_2']));
    else
        arc_prefori_ses2 = [];
    end
    arc_prefori_ses_all = [arc_prefori_ses1 arc_prefori_ses2];
    arc_prefori_scores = [arc_prefori_scores arc_prefori_ses_all];

    arc_pref_dscores_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores_2.mat']));
        arc_pref_dscores_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores_2']));
    else
        arc_pref_dscores_ses_2 = [];
    end
    arc_pref_dscores_ses_all = [arc_pref_dscores_ses_1 arc_pref_dscores_ses_2];
    arc_pref_dscores = [arc_pref_dscores arc_pref_dscores_ses_all];

    arc_k_max_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes_2.mat']));
        arc_k_max_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes_2']));
    else
        arc_k_max_ses_2 = [];
    end
    arc_k_max_ses_all = [arc_k_max_ses_1 arc_k_max_ses_2];
    arc_k_max_scores = [arc_k_max_scores arc_k_max_ses_all];

    arc_reliabilty_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'fit_reliability']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'fit_reliability_2.mat']));
        arc_reliabilty_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'fit_reliability_2']));
    else
        arc_reliabilty_ses_2 = [];
    end
    arc_reliability_ses_all = [arc_reliabilty_ses_1 arc_reliabilty_ses_2];
    arc_reliability_scores = [arc_reliability_scores arc_reliability_ses_all];

    arc_median_sep_pref_dscores_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'median_sep_d_scores']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'median_sep_d_scores_2.mat']));
        arc_median_sep_pref_dscores_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'median_sep_d_scores_2']));
    else
        arc_median_sep_pref_dscores_ses_2 = [];
    end
    arc_median_sep_pref_dscores_ses_all = [arc_median_sep_pref_dscores_ses_1 arc_median_sep_pref_dscores_ses_2];
    arc_median_pref_dscores = [arc_median_pref_dscores arc_median_sep_pref_dscores_ses_all];

    arc_wide_sharp_sep_pref_dscores_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'wide_sharp_sep_d_scores']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'wide_sharp_sep_d_scores_2.mat']));
        arc_wide_sharp_sep_pref_dscores_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'wide_sharp_sep_d_scores_2']));
    else
        arc_wide_sharp_sep_pref_dscores_ses_2 = [];
    end
    arc_wide_sharp_sep_pref_dscores_ses_all = [arc_wide_sharp_sep_pref_dscores_ses_1 arc_wide_sharp_sep_pref_dscores_ses_2];
    arc_wide_sharp_pref_dscores = [arc_wide_sharp_pref_dscores arc_wide_sharp_sep_pref_dscores_ses_all];
end

%%

arc_green_sharp_d1_d2_pref_d_all = [];
arc_green_sharp_d1_d3_pref_d_all = [];
arc_red_sharp_d1_d2_pref_d_all = [];
arc_red_sharp_d1_d3_pref_d_all = [];
arc_green_wide_d1_d2_pref_d_all = [];
arc_green_wide_d1_d3_pref_d_all = [];
arc_red_wide_d1_d2_pref_d_all = [];
arc_red_wide_d1_d3_pref_d_all = [];

for idata = 1:length(arc_prefori_scores)
%     arc_green_d1_pref = arc_prefori_ses_all.green_d1_prefori;
%     arc_green_d1_pref_all = [arc_green_d1_pref_all arc_green_d1_pref];
%     d2_pref = arc_prefori_ses_all.d2_prefori;
%     arc_green_d2_pref_all = [arc_green_d2_pref_all d2_pref];
%     d3_pref = arc_prefori_ses_all.d3_prefori;
%     arc_green_d3_pref_all = [arc_green_d3_pref_all d3_pref];
    

    arc_green_d1_d2_pref_d = arc_pref_dscores(idata).green_d_score_prefori_d1_d2;
    arc_green_d1_d2_pref_d_all = [arc_green_d1_d2_pref_d_all arc_green_d1_d2_pref_d];
    arc_green_d1_d3_pref_d = arc_pref_dscores(idata).green_d_score_prefori_d1_d3;
    arc_green_d1_d3_pref_d_all = [arc_green_d1_d3_pref_d_all arc_green_d1_d3_pref_d];

    arc_red_d1_d2_pref_d = arc_pref_dscores(idata).red_d_score_prefori_d1_d2;
    arc_red_d1_d2_pref_d_all = [arc_red_d1_d2_pref_d_all arc_red_d1_d2_pref_d];
    arc_red_d1_d3_pref_d = arc_pref_dscores(idata).red_d_score_prefori_d1_d3;
    arc_red_d1_d3_pref_d_all = [arc_red_d1_d3_pref_d_all arc_red_d1_d3_pref_d];


    arc_green_top_d1_d2_pref_d = arc_median_pref_dscores(idata).green_top_50_d_score_prefori_d1_d2;
    arc_green_top_d1_d2_pref_d_all = [arc_green_top_d1_d2_pref_d_all arc_green_top_d1_d2_pref_d];
    arc_green_top_d1_d3_pref_d = arc_median_pref_dscores(idata).green_top_50_d_score_prefori_d1_d3;
    arc_green_top_d1_d3_pref_d_all = [arc_green_top_d1_d3_pref_d_all arc_green_top_d1_d3_pref_d];

    arc_red_top_d1_d2_pref_d = arc_median_pref_dscores(idata).red_top_50_d_score_prefori_d1_d2;
    arc_red_top_d1_d2_pref_d_all = [arc_red_top_d1_d2_pref_d_all arc_red_top_d1_d2_pref_d];
    arc_red_top_d1_d3_pref_d = arc_median_pref_dscores(idata).red_top_50_d_score_prefori_d1_d3;
    arc_red_top_d1_d3_pref_d_all = [arc_red_top_d1_d3_pref_d_all arc_red_top_d1_d3_pref_d];

    arc_green_bot_d1_d2_pref_d = arc_median_pref_dscores(idata).green_bot_50_d_score_prefori_d1_d2;
    arc_green_bot_d1_d2_pref_d_all = [arc_green_bot_d1_d2_pref_d_all arc_green_bot_d1_d2_pref_d];
    arc_green_bot_d1_d3_pref_d = arc_median_pref_dscores(idata).green_bot_50_d_score_prefori_d1_d3;
    arc_green_bot_d1_d3_pref_d_all = [arc_green_bot_d1_d3_pref_d_all arc_green_bot_d1_d3_pref_d];

    arc_red_bot_d1_d2_pref_d = arc_median_pref_dscores(idata).red_bot_50_d_score_prefori_d1_d2;
    arc_red_bot_d1_d2_pref_d_all = [arc_red_bot_d1_d2_pref_d_all arc_red_bot_d1_d2_pref_d];
    arc_red_bot_d1_d3_pref_d = arc_median_pref_dscores(idata).red_bot_50_d_score_prefori_d1_d3;
    arc_red_bot_d1_d3_pref_d_all = [arc_red_bot_d1_d3_pref_d_all arc_red_bot_d1_d3_pref_d];

    arc_green_sharp_d1_d2_pref_d = arc_wide_sharp_pref_dscores(idata).green_sharp_d_score_prefori_d1_d2;
    arc_green_sharp_d1_d2_pref_d_all = [arc_green_sharp_d1_d2_pref_d_all arc_green_sharp_d1_d2_pref_d];
    arc_green_sharp_d1_d3_pref_d = arc_wide_sharp_pref_dscores(idata).green_sharp_d_score_prefori_d1_d3;
    arc_green_sharp_d1_d3_pref_d_all = [arc_green_sharp_d1_d3_pref_d_all arc_green_sharp_d1_d3_pref_d];

    arc_red_sharp_d1_d2_pref_d = arc_wide_sharp_pref_dscores(idata).red_sharp_d_score_prefori_d1_d2;
    arc_red_sharp_d1_d2_pref_d_all = [arc_red_sharp_d1_d2_pref_d_all arc_red_sharp_d1_d2_pref_d];
    arc_red_sharp_d1_d3_pref_d = arc_wide_sharp_pref_dscores(idata).red_sharp_d_score_prefori_d1_d3;
    arc_red_sharp_d1_d3_pref_d_all = [arc_red_sharp_d1_d3_pref_d_all arc_red_sharp_d1_d3_pref_d];

    arc_green_wide_d1_d2_pref_d = arc_wide_sharp_pref_dscores(idata).green_wide_d_score_prefori_d1_d2;
    arc_green_wide_d1_d2_pref_d_all = [arc_green_wide_d1_d2_pref_d_all arc_green_wide_d1_d2_pref_d];
    arc_green_wide_d1_d3_pref_d = arc_wide_sharp_pref_dscores(idata).green_wide_d_score_prefori_d1_d3;
    arc_green_wide_d1_d3_pref_d_all = [arc_green_wide_d1_d3_pref_d_all arc_green_wide_d1_d3_pref_d];

    arc_red_wide_d1_d2_pref_d = arc_wide_sharp_pref_dscores(idata).red_wide_d_score_prefori_d1_d2;
    arc_red_wide_d1_d2_pref_d_all = [arc_red_wide_d1_d2_pref_d_all arc_red_wide_d1_d2_pref_d];
    arc_red_wide_d1_d3_pref_d = arc_wide_sharp_pref_dscores(idata).red_wide_d_score_prefori_d1_d3;
    arc_red_wide_d1_d3_pref_d_all = [arc_red_wide_d1_d3_pref_d_all arc_red_wide_d1_d3_pref_d];


    arc_green_d1_k = arc_k_max_scores(idata).green_d1_k;
    arc_green_d1_k_all = [arc_green_d1_k_all arc_green_d1_k];
    arc_green_d2_k = arc_k_max_scores(idata).green_d2_k;
    arc_green_d2_k_all = [arc_green_d2_k_all arc_green_d2_k];
    arc_green_d3_k = arc_k_max_scores(idata).green_d3_k;
    arc_green_d3_k_all = [arc_green_d3_k_all arc_green_d3_k];

    arc_green_d1_max = arc_k_max_scores(idata).green_d1_max;
    arc_green_d1_max_all = [arc_green_d1_max_all arc_green_d1_max];
    arc_green_d2_max = arc_k_max_scores(idata).green_d2_max;
    arc_green_d2_max_all = [arc_green_d2_max_all arc_green_d2_max];
    arc_green_d3_max = arc_k_max_scores(idata).green_d3_max;
    arc_green_d3_max_all = [arc_green_d3_max_all arc_green_d3_max];

    arc_red_d1_k = arc_k_max_scores(idata).red_d1_k;
    arc_red_d1_k_all = [arc_red_d1_k_all arc_red_d1_k];
    arc_red_d2_k = arc_k_max_scores(idata).red_d2_k;
    arc_red_d2_k_all = [arc_red_d2_k_all arc_red_d2_k];
    arc_red_d3_k = arc_k_max_scores(idata).red_d3_k;
    arc_red_d3_k_all = [arc_red_d3_k_all arc_red_d3_k];

    arc_red_d1_max = arc_k_max_scores(idata).red_d1_max;
    arc_red_d1_max_all = [arc_red_d1_max_all arc_red_d1_max];
    arc_red_d2_max = arc_k_max_scores(idata).red_d2_max;
    arc_red_d2_max_all = [arc_red_d2_max_all arc_red_d2_max];
    arc_red_d3_max = arc_k_max_scores(idata).red_d3_max;
    arc_red_d3_max_all = [arc_red_d3_max_all arc_red_d3_max];

    arc_green_d1_reliability = arc_reliability_scores(idata).green_d1_reliability;
    arc_green_d1_reliability_all = [arc_green_d1_reliability_all arc_green_d1_reliability];
    arc_green_d2_reliability = arc_reliability_scores(idata).green_d2_reliability;
    arc_green_d2_reliability_all = [arc_green_d2_reliability_all arc_green_d2_reliability];
    arc_green_d3_reliability = arc_reliability_scores(idata).green_d3_reliability;
    arc_green_d3_reliability_all = [arc_green_d3_reliability_all arc_green_d3_reliability];

    arc_red_d1_reliability = arc_reliability_scores(idata).red_d1_reliability;
    arc_red_d1_reliability_all = [arc_red_d1_reliability_all arc_red_d1_reliability];
    arc_red_d2_reliability = arc_reliability_scores(idata).red_d2_reliability;
    arc_red_d2_reliability_all = [arc_red_d2_reliability_all arc_red_d2_reliability];
    arc_red_d3_reliability = arc_reliability_scores(idata).red_d3_reliability;
    arc_red_d3_reliability_all = [arc_red_d3_reliability_all arc_red_d3_reliability];

  
end

%%
fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(arc_green_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(arc_green_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(arc_red_d1_d2_pref_d_all);
set(l, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
m=cdfplot(arc_red_d1_d3_pref_d_all);
set(m, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 - 2 (Non-KRAB)'], ['Ses1 - 3 (Non-KRAB)'], ['Ses1 - 2 (KRAB)'], ['Ses1 - 3 (KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Arc Enhancer Mice'])
hold off

print(fullfile(realfnout, ['grace_arcprom_changeprefcdf_redVgreen.pdf']), '-dpdf', '-bestfit')


%%
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(arc_green_d1_k_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(arc_green_d2_k_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(arc_green_d3_k_all);
set(l, 'LineStyle', ':', 'Color', 'g', 'LineWidth',2);
m=cdfplot(arc_red_d1_k_all);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
n=cdfplot(arc_red_d2_k_all);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
o=cdfplot(arc_red_d3_k_all);
set(o, 'LineStyle', ':', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 (Non-KRAB)'], ['Ses2 (Non-KRAB)'], ['Ses3 (Non-KRAB)'], ['Ses1 (KRAB)'], ['Ses2 (KRAB)'], ['Ses3 (KRAB)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['k value'])
ylabel(['% of cells'])
title(['Arc Enhancer Mice'])
hold off

print(fullfile(realfnout, ['grace_arcprom_kbysession_redVgreen.pdf']), '-dpdf', '-bestfit')

%%
arc_green_k_all = [arc_green_d1_k_all arc_green_d2_k_all arc_green_d3_k_all];
arc_red_k_all = [arc_red_d1_k_all arc_red_d2_k_all arc_red_d3_k_all];

fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(arc_green_k_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(arc_red_k_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['All Sessions (Non-KRAB)'], ['All Sessions (KRAB)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['k value'])
ylabel(['% of cells'])
title(['Arc Enhancer Mice'])
hold off

print(fullfile(realfnout, ['grace_arcprom_kallsessions_redVgreen.pdf']), '-dpdf', '-bestfit')

%%
%day 1 k
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(arc_green_d1_k_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(arc_red_d1_k_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['All Sessions (Non-KRAB)'], ['All Sessions (KRAB)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['k value'])
ylabel(['% of cells'])
title(['Arc Enhancer Mice'])
hold off
print(fullfile(realfnout, ['grace_arcprom_d1_k_redVgreen.pdf']), '-dpdf', '-bestfit')

%%
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(arc_green_d1_max_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(arc_green_d2_max_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(arc_green_d3_max_all);
set(l, 'LineStyle', ':', 'Color', 'g', 'LineWidth',2);
m=cdfplot(arc_red_d1_max_all);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
n=cdfplot(arc_red_d2_max_all);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
o=cdfplot(arc_red_d3_max_all);
set(o, 'LineStyle', ':', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 (Non-KRAB)'], ['Ses2 (Non-KRAB)'], ['Ses3 (Non-KRAB)'], ['Ses1 (KRAB)'], ['Ses2 (KRAB)'], ['Ses3 (KRAB)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F value'])
ylabel(['% of cells'])
title(['Arc Enhancer Mice'])
hold off

print(fullfile(realfnout, ['grace_arcprom_maxbysession_redVgreen.pdf']), '-dpdf', '-bestfit')

%%
arc_green_max_all = [arc_green_d1_max_all arc_green_d2_max_all arc_green_d3_max_all];
arc_red_max_all = [arc_red_d1_max_all arc_red_d2_max_all arc_red_d3_max_all];

fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(arc_green_max_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(arc_red_max_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['All Sessions (Non-KRAB)'], ['All Sessions (KRAB)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F value'])
ylabel(['% of cells'])
title(['Arc Enhancer Mice'])
hold off

print(fullfile(realfnout, ['grace_arcprom_maxallsessions_redVgreen.pdf']), '-dpdf', '-bestfit')

%%
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(arc_green_d1_max_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(arc_red_d1_max_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['All Sessions (Non-KRAB)'], ['All Sessions (KRAB)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F value'])
ylabel(['% of cells'])
title(['Arc Enhancer Mice'])
hold off
print(fullfile(realfnout, ['grace_arcprom_d1_max_redVgreen.pdf']), '-dpdf', '-bestfit')

%%
%fit reliability
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(arc_green_d1_reliability_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(arc_green_d2_reliability_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(arc_green_d3_reliability_all);
set(l, 'LineStyle', ':', 'Color', 'g', 'LineWidth',2);
m=cdfplot(arc_red_d1_reliability_all);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
n=cdfplot(arc_red_d2_reliability_all);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
o=cdfplot(arc_red_d3_reliability_all);
set(o, 'LineStyle', ':', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 (Non-KRAB)'], ['Ses2 (Non-KRAB)'], ['Ses3 (Non-KRAB)'], ['Ses1 (KRAB)'], ['Ses2 (KRAB)'], ['Ses3 (KRAB)'], 'Location', 'southeast')
% xlim([0 1])
xlabel(['Fit reliability'])
ylabel(['% of cells'])
title(['Arc Promoter Mice'])
hold off
print(fullfile(realfnout, ['grace_arcprom_reliabilitybysession_redVgreen.pdf']), '-dpdf', '-bestfit')

%%
arc_green_reliability_all = [arc_green_d1_reliability_all arc_green_d2_reliability_all arc_green_d3_reliability_all];
arc_red_reliability_all = [arc_red_d1_reliability_all arc_red_d2_reliability_all arc_red_d3_reliability_all];
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(arc_green_reliability_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(arc_red_reliability_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['All Sessions (Non-KRAB)'], ['All Sessions (KRAB)'], 'Location', 'southeast')
% xlim([0 1])
xlabel(['Fit reliability'])
ylabel(['% of cells'])
title(['Arc Promoter Mice'])
hold off
print(fullfile(realfnout, ['grace_arcprom_reliabilityallsessions_redVgreen.pdf']), '-dpdf', '-bestfit')

%%
%median split analysis for 4d dark group
fig = figure;
% sgtitle('Pref Ori Changes')
a=subplot(2,2,1)
h=cdfplot(arc_green_bot_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(arc_green_bot_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
legend(['Ses1 - 2 (Non-KRAB)'], ['Ses1 - 3 (Non-KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Bottom 50% Responsivity (Arc Enh)'])

% sgtitle('Pref Ori Changes')
b=subplot(2,2,2)
h=cdfplot(arc_red_bot_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(arc_red_bot_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 - 2 (KRAB)'], ['Ses1 - 3 (KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Bottom 50% Responsivity (Arc Enh)'])

% sgtitle('Pref Ori Changes')
a=subplot(2,2,3)
h=cdfplot(arc_green_top_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(arc_green_top_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
legend(['Ses1 - 2 (Non-KRAB)'], ['Ses1 - 3 (Non-KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Top 50% Responsivity (Arc Enh)'])

% sgtitle('Pref Ori Changes')
b=subplot(2,2,4)
h=cdfplot(arc_red_top_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(arc_red_top_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 - 2 (KRAB)'], ['Ses1 - 3 (KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Top 50% Responsivity (Arc Enh)'])



print(fullfile(realfnout, ['grace_arcprom_medsplit_changeprefcdf_redVgreen.pdf']), '-dpdf', '-bestfit')

%%
arc_green_bot_pref_d_all = [arc_green_bot_d1_d2_pref_d_all arc_green_bot_d1_d3_pref_d_all];
arc_green_top_pref_d_all = [arc_green_top_d1_d2_pref_d_all arc_green_top_d1_d3_pref_d_all];
arc_red_bot_pref_d_all = [arc_red_bot_d1_d2_pref_d_all arc_red_bot_d1_d3_pref_d_all];
arc_red_top_pref_d_all = [arc_red_top_d1_d2_pref_d_all arc_red_top_d1_d3_pref_d_all];
fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(arc_green_bot_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(arc_green_top_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(arc_red_bot_pref_d_all);
set(l, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
m=cdfplot(arc_red_top_pref_d_all);
set(m, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Bottom 50% (Non-KRAB)'], ['Top 50% (Non-KRAB)'], ['Bottom 50% (KRAB)'], ['Top 50% (KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Arc Enhancer Mice'])
hold off
print(fullfile(realfnout, ['grace_arcprom_medsplit_changeprefcdf_ALLDAYS_redVgreen.pdf']), '-dpdf', '-bestfit')


%%
%% SAME AS ABOVE BUT LacZ
%%

lacz_green_d1_pref_all = [];
lacz_green_d2_pref_all = [];
lacz_green_d3_pref_all = [];
lacz_green_d1_d2_pref_d_all = [];
lacz_green_d1_d3_pref_d_all = [];
lacz_green_d1_k_all = [];
lacz_green_d2_k_all = [];
lacz_green_d3_k_all = [];
lacz_green_d1_max_all = [];
lacz_green_d2_max_all = [];
lacz_green_d3_max_all = [];

lacz_green_d1_reliability_all = [];
lacz_green_d2_reliability_all = [];
lacz_green_d3_reliability_all = [];

lacz_green_top_d1_d2_pref_d_all = [];
lacz_green_top_d1_d3_pref_d_all = [];

lacz_green_bot_d1_d2_pref_d_all = [];
lacz_green_bot_d1_d3_pref_d_all = [];

lacz_red_d1_pref_all = [];
lacz_red_d2_pref_all = [];
lacz_red_d3_pref_all = [];
lacz_red_d1_d2_pref_d_all = [];
lacz_red_d1_d3_pref_d_all = [];
lacz_red_d1_k_all = [];
lacz_red_d2_k_all = [];
lacz_red_d3_k_all = [];
lacz_red_d1_max_all = [];
lacz_red_d2_max_all = [];
lacz_red_d3_max_all = [];

lacz_red_d1_reliability_all = [];
lacz_red_d2_reliability_all = [];
lacz_red_d3_reliability_all = [];

lacz_red_top_d1_d2_pref_d_all = [];
lacz_red_top_d1_d3_pref_d_all = [];

lacz_red_bot_d1_d2_pref_d_all = [];
lacz_red_bot_d1_d3_pref_d_all = [];



lacz_prefori_scores = [];
lacz_pref_dscores = [];
lacz_k_max_scores = [];
lacz_median_pref_dscores = [];
lacz_reliability_scores = [];
lacz_wide_sharp_pref_dscores = [];



list = [49 55 61];
for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};

    lacz_prefori_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes_2.mat']));
        lacz_prefori_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes_2']));
    else
        lacz_prefori_ses2 = [];
    end
    lacz_prefori_ses_all = [lacz_prefori_ses1 lacz_prefori_ses2];
    lacz_prefori_scores = [lacz_prefori_scores lacz_prefori_ses_all];

    lacz_pref_dscores_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores_2.mat']));
        lacz_pref_dscores_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores_2']));
    else
        lacz_pref_dscores_ses_2 = [];
    end
    lacz_pref_dscores_ses_all = [lacz_pref_dscores_ses_1 lacz_pref_dscores_ses_2];
    lacz_pref_dscores = [lacz_pref_dscores lacz_pref_dscores_ses_all];

    lacz_k_max_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes_2.mat']));
        lacz_k_max_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes_2']));
    else
        lacz_k_max_ses_2 = [];
    end
    lacz_k_max_ses_all = [lacz_k_max_ses_1 lacz_k_max_ses_2];
    lacz_k_max_scores = [lacz_k_max_scores lacz_k_max_ses_all];

    lacz_reliabilty_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'fit_reliability']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'fit_reliability_2.mat']));
        lacz_reliabilty_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'fit_reliability_2']));
    else
        lacz_reliabilty_ses_2 = [];
    end
    lacz_reliability_ses_all = [lacz_reliabilty_ses_1 lacz_reliabilty_ses_2];
    lacz_reliability_scores = [lacz_reliability_scores lacz_reliability_ses_all];

    lacz_median_sep_pref_dscores_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'median_sep_d_scores']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'median_sep_d_scores_2.mat']));
        lacz_median_sep_pref_dscores_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'median_sep_d_scores_2']));
    else
        lacz_median_sep_pref_dscores_ses_2 = [];
    end
    lacz_median_sep_pref_dscores_ses_all = [lacz_median_sep_pref_dscores_ses_1 lacz_median_sep_pref_dscores_ses_2];
    lacz_median_pref_dscores = [lacz_median_pref_dscores lacz_median_sep_pref_dscores_ses_all];

    lacz_wide_sharp_sep_pref_dscores_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'wide_sharp_sep_d_scores']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'wide_sharp_sep_d_scores_2.mat']));
        lacz_wide_sharp_sep_pref_dscores_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'wide_sharp_sep_d_scores_2']));
    else
        lacz_wide_sharp_sep_pref_dscores_ses_2 = [];
    end
    lacz_wide_sharp_sep_pref_dscores_ses_all = [lacz_wide_sharp_sep_pref_dscores_ses_1 lacz_wide_sharp_sep_pref_dscores_ses_2];
    lacz_wide_sharp_pref_dscores = [lacz_wide_sharp_pref_dscores lacz_wide_sharp_sep_pref_dscores_ses_all];
end



%%
lacz_green_sharp_d1_d2_pref_d_all = [];
lacz_green_sharp_d1_d3_pref_d_all = [];
lacz_red_sharp_d1_d2_pref_d_all = [];
lacz_red_sharp_d1_d3_pref_d_all = [];
lacz_green_wide_d1_d2_pref_d_all = [];
lacz_green_wide_d1_d3_pref_d_all = [];
lacz_red_wide_d1_d2_pref_d_all = [];
lacz_red_wide_d1_d3_pref_d_all = [];

for idata = 1:length(lacz_prefori_scores)
%     lacz_green_d1_pref = lacz_prefori_ses_all.green_d1_prefori;
%     lacz_green_d1_pref_all = [lacz_green_d1_pref_all lacz_green_d1_pref];
%     d2_pref = lacz_prefori_ses_all.d2_prefori;
%     lacz_green_d2_pref_all = [lacz_green_d2_pref_all d2_pref];
%     d3_pref = lacz_prefori_ses_all.d3_prefori;
%     lacz_green_d3_pref_all = [lacz_green_d3_pref_all d3_pref];
    

    lacz_green_d1_d2_pref_d = lacz_pref_dscores(idata).green_d_score_prefori_d1_d2;
    lacz_green_d1_d2_pref_d_all = [lacz_green_d1_d2_pref_d_all lacz_green_d1_d2_pref_d];
    lacz_green_d1_d3_pref_d = lacz_pref_dscores(idata).green_d_score_prefori_d1_d3;
    lacz_green_d1_d3_pref_d_all = [lacz_green_d1_d3_pref_d_all lacz_green_d1_d3_pref_d];

    lacz_red_d1_d2_pref_d = lacz_pref_dscores(idata).red_d_score_prefori_d1_d2;
    lacz_red_d1_d2_pref_d_all = [lacz_red_d1_d2_pref_d_all lacz_red_d1_d2_pref_d];
    lacz_red_d1_d3_pref_d = lacz_pref_dscores(idata).red_d_score_prefori_d1_d3;
    lacz_red_d1_d3_pref_d_all = [lacz_red_d1_d3_pref_d_all lacz_red_d1_d3_pref_d];


    lacz_green_top_d1_d2_pref_d = lacz_median_pref_dscores(idata).green_top_50_d_score_prefori_d1_d2;
    lacz_green_top_d1_d2_pref_d_all = [lacz_green_top_d1_d2_pref_d_all lacz_green_top_d1_d2_pref_d];
    lacz_green_top_d1_d3_pref_d = lacz_median_pref_dscores(idata).green_top_50_d_score_prefori_d1_d3;
    lacz_green_top_d1_d3_pref_d_all = [lacz_green_top_d1_d3_pref_d_all lacz_green_top_d1_d3_pref_d];

    lacz_red_top_d1_d2_pref_d = lacz_median_pref_dscores(idata).red_top_50_d_score_prefori_d1_d2;
    lacz_red_top_d1_d2_pref_d_all = [lacz_red_top_d1_d2_pref_d_all lacz_red_top_d1_d2_pref_d];
    lacz_red_top_d1_d3_pref_d = lacz_median_pref_dscores(idata).red_top_50_d_score_prefori_d1_d3;
    lacz_red_top_d1_d3_pref_d_all = [lacz_red_top_d1_d3_pref_d_all lacz_red_top_d1_d3_pref_d];

    lacz_green_bot_d1_d2_pref_d = lacz_median_pref_dscores(idata).green_bot_50_d_score_prefori_d1_d2;
    lacz_green_bot_d1_d2_pref_d_all = [lacz_green_bot_d1_d2_pref_d_all lacz_green_bot_d1_d2_pref_d];
    lacz_green_bot_d1_d3_pref_d = lacz_median_pref_dscores(idata).green_bot_50_d_score_prefori_d1_d3;
    lacz_green_bot_d1_d3_pref_d_all = [lacz_green_bot_d1_d3_pref_d_all lacz_green_bot_d1_d3_pref_d];

    lacz_red_bot_d1_d2_pref_d = lacz_median_pref_dscores(idata).red_bot_50_d_score_prefori_d1_d2;
    lacz_red_bot_d1_d2_pref_d_all = [lacz_red_bot_d1_d2_pref_d_all lacz_red_bot_d1_d2_pref_d];
    lacz_red_bot_d1_d3_pref_d = lacz_median_pref_dscores(idata).red_bot_50_d_score_prefori_d1_d3;
    lacz_red_bot_d1_d3_pref_d_all = [lacz_red_bot_d1_d3_pref_d_all lacz_red_bot_d1_d3_pref_d];

    lacz_green_sharp_d1_d2_pref_d = lacz_wide_sharp_pref_dscores(idata).green_sharp_d_score_prefori_d1_d2;
    lacz_green_sharp_d1_d2_pref_d_all = [lacz_green_sharp_d1_d2_pref_d_all lacz_green_sharp_d1_d2_pref_d];
    lacz_green_sharp_d1_d3_pref_d = lacz_wide_sharp_pref_dscores(idata).green_sharp_d_score_prefori_d1_d3;
    lacz_green_sharp_d1_d3_pref_d_all = [lacz_green_sharp_d1_d3_pref_d_all lacz_green_sharp_d1_d3_pref_d];

    lacz_red_sharp_d1_d2_pref_d = lacz_wide_sharp_pref_dscores(idata).red_sharp_d_score_prefori_d1_d2;
    lacz_red_sharp_d1_d2_pref_d_all = [lacz_red_sharp_d1_d2_pref_d_all lacz_red_sharp_d1_d2_pref_d];
    lacz_red_sharp_d1_d3_pref_d = lacz_wide_sharp_pref_dscores(idata).red_sharp_d_score_prefori_d1_d3;
    lacz_red_sharp_d1_d3_pref_d_all = [lacz_red_sharp_d1_d3_pref_d_all lacz_red_sharp_d1_d3_pref_d];

    lacz_green_wide_d1_d2_pref_d = lacz_wide_sharp_pref_dscores(idata).green_wide_d_score_prefori_d1_d2;
    lacz_green_wide_d1_d2_pref_d_all = [lacz_green_wide_d1_d2_pref_d_all lacz_green_wide_d1_d2_pref_d];
    lacz_green_wide_d1_d3_pref_d = lacz_wide_sharp_pref_dscores(idata).green_wide_d_score_prefori_d1_d3;
    lacz_green_wide_d1_d3_pref_d_all = [lacz_green_wide_d1_d3_pref_d_all lacz_green_wide_d1_d3_pref_d];

    lacz_red_wide_d1_d2_pref_d = lacz_wide_sharp_pref_dscores(idata).red_wide_d_score_prefori_d1_d2;
    lacz_red_wide_d1_d2_pref_d_all = [lacz_red_wide_d1_d2_pref_d_all lacz_red_wide_d1_d2_pref_d];
    lacz_red_wide_d1_d3_pref_d = lacz_wide_sharp_pref_dscores(idata).red_wide_d_score_prefori_d1_d3;
    lacz_red_wide_d1_d3_pref_d_all = [lacz_red_wide_d1_d3_pref_d_all lacz_red_wide_d1_d3_pref_d];


    lacz_green_d1_k = lacz_k_max_scores(idata).green_d1_k;
    lacz_green_d1_k_all = [lacz_green_d1_k_all lacz_green_d1_k];
    lacz_green_d2_k = lacz_k_max_scores(idata).green_d2_k;
    lacz_green_d2_k_all = [lacz_green_d2_k_all lacz_green_d2_k];
    lacz_green_d3_k = lacz_k_max_scores(idata).green_d3_k;
    lacz_green_d3_k_all = [lacz_green_d3_k_all lacz_green_d3_k];

    lacz_green_d1_max = lacz_k_max_scores(idata).green_d1_max;
    lacz_green_d1_max_all = [lacz_green_d1_max_all lacz_green_d1_max];
    lacz_green_d2_max = lacz_k_max_scores(idata).green_d2_max;
    lacz_green_d2_max_all = [lacz_green_d2_max_all lacz_green_d2_max];
    lacz_green_d3_max = lacz_k_max_scores(idata).green_d3_max;
    lacz_green_d3_max_all = [lacz_green_d3_max_all lacz_green_d3_max];

    lacz_red_d1_k = lacz_k_max_scores(idata).red_d1_k;
    lacz_red_d1_k_all = [lacz_red_d1_k_all lacz_red_d1_k];
    lacz_red_d2_k = lacz_k_max_scores(idata).red_d2_k;
    lacz_red_d2_k_all = [lacz_red_d2_k_all lacz_red_d2_k];
    lacz_red_d3_k = lacz_k_max_scores(idata).red_d3_k;
    lacz_red_d3_k_all = [lacz_red_d3_k_all lacz_red_d3_k];

    lacz_red_d1_max = lacz_k_max_scores(idata).red_d1_max;
    lacz_red_d1_max_all = [lacz_red_d1_max_all lacz_red_d1_max];
    lacz_red_d2_max = lacz_k_max_scores(idata).red_d2_max;
    lacz_red_d2_max_all = [lacz_red_d2_max_all lacz_red_d2_max];
    lacz_red_d3_max = lacz_k_max_scores(idata).red_d3_max;
    lacz_red_d3_max_all = [lacz_red_d3_max_all lacz_red_d3_max];

    lacz_green_d1_reliability = lacz_reliability_scores(idata).green_d1_reliability;
    lacz_green_d1_reliability_all = [lacz_green_d1_reliability_all lacz_green_d1_reliability];
    lacz_green_d2_reliability = lacz_reliability_scores(idata).green_d2_reliability;
    lacz_green_d2_reliability_all = [lacz_green_d2_reliability_all lacz_green_d2_reliability];
    lacz_green_d3_reliability = lacz_reliability_scores(idata).green_d3_reliability;
    lacz_green_d3_reliability_all = [lacz_green_d3_reliability_all lacz_green_d3_reliability];

    lacz_red_d1_reliability = lacz_reliability_scores(idata).red_d1_reliability;
    lacz_red_d1_reliability_all = [lacz_red_d1_reliability_all lacz_red_d1_reliability];
    lacz_red_d2_reliability = lacz_reliability_scores(idata).red_d2_reliability;
    lacz_red_d2_reliability_all = [lacz_red_d2_reliability_all lacz_red_d2_reliability];
    lacz_red_d3_reliability = lacz_reliability_scores(idata).red_d3_reliability;
    lacz_red_d3_reliability_all = [lacz_red_d3_reliability_all lacz_red_d3_reliability];

  
end



%%
fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(lacz_green_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(lacz_green_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(lacz_red_d1_d2_pref_d_all);
set(l, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
m=cdfplot(lacz_red_d1_d3_pref_d_all);
set(m, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 - 2 (Non-KRAB)'], ['Ses1 - 3 (Non-KRAB)'], ['Ses1 - 2 (KRAB)'], ['Ses1 - 3 (KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Lacz Mice'])
hold off

print(fullfile(realfnout, ['grace_lacz_changeprefcdf_redVgreen.pdf']), '-dpdf', '-bestfit')



%%
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(lacz_green_d1_k_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(lacz_green_d2_k_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(lacz_green_d3_k_all);
set(l, 'LineStyle', ':', 'Color', 'g', 'LineWidth',2);
m=cdfplot(lacz_red_d1_k_all);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
n=cdfplot(lacz_red_d2_k_all);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
o=cdfplot(lacz_red_d3_k_all);
set(o, 'LineStyle', ':', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 (Non-KRAB)'], ['Ses2 (Non-KRAB)'], ['Ses3 (Non-KRAB)'], ['Ses1 (KRAB)'], ['Ses2 (KRAB)'], ['Ses3 (KRAB)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['k value'])
ylabel(['% of cells'])
title(['Lacz Mice'])
hold off

print(fullfile(realfnout, ['grace_lacz_kbysession_redVgreen.pdf']), '-dpdf', '-bestfit')



%%
lacz_green_k_all = [lacz_green_d1_k_all lacz_green_d2_k_all lacz_green_d3_k_all];
lacz_red_k_all = [lacz_red_d1_k_all lacz_red_d2_k_all lacz_red_d3_k_all];

fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(lacz_green_k_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(lacz_red_k_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['All Sessions (Non-KRAB)'], ['All Sessions (KRAB)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['k value'])
ylabel(['% of cells'])
title(['Lacz Mice'])
hold off

print(fullfile(realfnout, ['grace_lacz_kallsessions_redVgreen.pdf']), '-dpdf', '-bestfit')

%%
%day 1 k
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(lacz_green_d1_k_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(lacz_red_d1_k_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['All Sessions (Non-KRAB)'], ['All Sessions (KRAB)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['k value'])
ylabel(['% of cells'])
title(['Lacz Mice'])
hold off
print(fullfile(realfnout, ['grace_lacz_d1_k_redVgreen.pdf']), '-dpdf', '-bestfit')
%%
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(lacz_green_d1_max_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(lacz_green_d2_max_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(lacz_green_d3_max_all);
set(l, 'LineStyle', ':', 'Color', 'g', 'LineWidth',2);
m=cdfplot(lacz_red_d1_max_all);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
n=cdfplot(lacz_red_d2_max_all);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
o=cdfplot(lacz_red_d3_max_all);
set(o, 'LineStyle', ':', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 (Non-KRAB)'], ['Ses2 (Non-KRAB)'], ['Ses3 (Non-KRAB)'], ['Ses1 (KRAB)'], ['Ses2 (KRAB)'], ['Ses3 (KRAB)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F value'])
ylabel(['% of cells'])
title(['Lacz Mice'])
hold off

print(fullfile(realfnout, ['grace_lacz_maxbysession_redVgreen.pdf']), '-dpdf', '-bestfit')



%%
lacz_green_max_all = [lacz_green_d1_max_all lacz_green_d2_max_all lacz_green_d3_max_all];
lacz_red_max_all = [lacz_red_d1_max_all lacz_red_d2_max_all lacz_red_d3_max_all];

fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(lacz_green_max_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(lacz_red_max_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['All Sessions (Non-KRAB)'], ['All Sessions (KRAB)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F value'])
ylabel(['% of cells'])
title(['Lacz Mice'])
hold off

print(fullfile(realfnout, ['grace_lacz_maxallsessions_redVgreen.pdf']), '-dpdf', '-bestfit')

%%
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(lacz_green_d1_max_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(lacz_red_d1_max_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['All Sessions (Non-KRAB)'], ['All Sessions (KRAB)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F value'])
ylabel(['% of cells'])
title(['Lacz Mice'])
hold off
print(fullfile(realfnout, ['grace_lacz_d1_max_redVgreen.pdf']), '-dpdf', '-bestfit')


%%
%fit reliability
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(lacz_green_d1_reliability_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(lacz_green_d2_reliability_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(lacz_green_d3_reliability_all);
set(l, 'LineStyle', ':', 'Color', 'g', 'LineWidth',2);
m=cdfplot(lacz_red_d1_reliability_all);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
n=cdfplot(lacz_red_d2_reliability_all);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
o=cdfplot(lacz_red_d3_reliability_all);
set(o, 'LineStyle', ':', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 (Non-KRAB)'], ['Ses2 (Non-KRAB)'], ['Ses3 (Non-KRAB)'], ['Ses1 (KRAB)'], ['Ses2 (KRAB)'], ['Ses3 (KRAB)'], 'Location', 'southeast')
% xlim([0 1])
xlabel(['Fit reliability'])
ylabel(['% of cells'])
title(['Lacz Mice'])
hold off
print(fullfile(realfnout, ['grace_lacz_reliabilitybysession_redVgreen.pdf']), '-dpdf', '-bestfit')

%%
lacz_green_reliability_all = [lacz_green_d1_reliability_all lacz_green_d2_reliability_all lacz_green_d3_reliability_all];
lacz_red_reliability_all = [lacz_red_d1_reliability_all lacz_red_d2_reliability_all lacz_red_d3_reliability_all];
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(lacz_green_reliability_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(lacz_red_reliability_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['All Sessions (Non-KRAB)'], ['All Sessions (KRAB)'], 'Location', 'southeast')
% xlim([0 1])
xlabel(['Fit reliability'])
ylabel(['% of cells'])
title(['Lacz Mice'])
hold off
print(fullfile(realfnout, ['grace_lacz_reliabilityallsessions_redVgreen.pdf']), '-dpdf', '-bestfit')



%%
%median split analysis for 4d dark group
fig = figure;
% sgtitle('Pref Ori Changes')
a=subplot(2,2,1)
h=cdfplot(lacz_green_bot_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(lacz_green_bot_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
legend(['Ses1 - 2 (Non-KRAB)'], ['Ses1 - 3 (Non-KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Bottom 50% Responsivity (Lacz)'])

% sgtitle('Pref Ori Changes')
b=subplot(2,2,2)
h=cdfplot(lacz_red_bot_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(lacz_red_bot_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 - 2 (KRAB)'], ['Ses1 - 3 (KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Bottom 50% Responsivity (Lacz)'])

% sgtitle('Pref Ori Changes')
a=subplot(2,2,3)
h=cdfplot(lacz_green_top_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(lacz_green_top_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
legend(['Ses1 - 2 (Non-KRAB)'], ['Ses1 - 3 (Non-KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Top 50% Responsivity (Lacz)'])

% sgtitle('Pref Ori Changes')
b=subplot(2,2,4)
h=cdfplot(lacz_red_top_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(lacz_red_top_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 - 2 (KRAB)'], ['Ses1 - 3 (KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Top 50% Responsivity (Lacz)'])

print(fullfile(realfnout, ['grace_lacz_medsplit_changeprefcdf_redVgreen.pdf']), '-dpdf', '-bestfit')

%%
lacz_green_bot_pref_d_all = [lacz_green_bot_d1_d2_pref_d_all lacz_green_bot_d1_d3_pref_d_all];
lacz_green_top_pref_d_all = [lacz_green_top_d1_d2_pref_d_all lacz_green_top_d1_d3_pref_d_all];
lacz_red_bot_pref_d_all = [lacz_red_bot_d1_d2_pref_d_all lacz_red_bot_d1_d3_pref_d_all];
lacz_red_top_pref_d_all = [lacz_red_top_d1_d2_pref_d_all lacz_red_top_d1_d3_pref_d_all];
fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(lacz_green_bot_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(lacz_green_top_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(lacz_red_bot_pref_d_all);
set(l, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
m=cdfplot(lacz_red_top_pref_d_all);
set(m, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Bottom 50% (Non-KRAB)'], ['Top 50% (Non-KRAB)'], ['Bottom 50% (KRAB)'], ['Top 50% (KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['LacZ Mice'])
hold off
print(fullfile(realfnout, ['grace_lacz_medsplit_changeprefcdf_ALLDAYS_redVgreen.pdf']), '-dpdf', '-bestfit')

%% NOW DO SOME LACZ VS ARC COMPARISONS***


%%
arc_green_wide_pref_d_all = [arc_green_wide_d1_d2_pref_d_all arc_green_wide_d1_d3_pref_d_all];
arc_green_sharp_pref_d_all = [arc_green_sharp_d1_d2_pref_d_all arc_green_sharp_d1_d3_pref_d_all];
arc_red_wide_pref_d_all = [arc_red_wide_d1_d2_pref_d_all arc_red_wide_d1_d3_pref_d_all];
arc_red_sharp_pref_d_all = [arc_red_sharp_d1_d2_pref_d_all arc_red_sharp_d1_d3_pref_d_all];

lacz_green_wide_pref_d_all = [lacz_green_wide_d1_d2_pref_d_all lacz_green_wide_d1_d3_pref_d_all];
lacz_green_sharp_pref_d_all = [lacz_green_sharp_d1_d2_pref_d_all lacz_green_sharp_d1_d3_pref_d_all];
lacz_red_wide_pref_d_all = [lacz_red_wide_d1_d2_pref_d_all lacz_red_wide_d1_d3_pref_d_all];
lacz_red_sharp_pref_d_all = [lacz_red_sharp_d1_d2_pref_d_all lacz_red_sharp_d1_d3_pref_d_all];


%%
fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(arc_red_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(arc_red_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
l=cdfplot(lacz_red_d1_d2_pref_d_all);
set(l, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
m=cdfplot(lacz_red_d1_d3_pref_d_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Ses1 - 2 (Arc)'], ['Ses1 - 3 (Arc)'], ['Ses1 - 2 (LacZ)'], ['Ses1 - 3 (LacZ)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['KRAB Cells'])
hold off

print(fullfile(realfnout, ['grace_allmice_changeprefcdf_red_.pdf']), '-dpdf', '-bestfit')

fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(arc_green_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(arc_green_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
l=cdfplot(lacz_green_d1_d2_pref_d_all);
set(l, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
m=cdfplot(lacz_green_d1_d3_pref_d_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Ses1 - 2 (Arc)'], ['Ses1 - 3 (Arc)'], ['Ses1 - 2 (LacZ)'], ['Ses1 - 3 (LacZ)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Non-KRAB Cells'])
hold off

print(fullfile(realfnout, ['grace_allmice_changeprefcdf_green_.pdf']), '-dpdf', '-bestfit')





%%
%median split analysis for 4d dark group
fig = figure;
% sgtitle('Pref Ori Changes')
a=subplot(2,2,1)
h=cdfplot(arc_red_bot_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(arc_red_bot_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
l=cdfplot(lacz_red_bot_d1_d2_pref_d_all);
set(l, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
m=cdfplot(lacz_red_bot_d1_d3_pref_d_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Ses1 - 2 (Arc Enh)'], ['Ses1 - 3 (Arc Enh)'], ['Ses1 - 2 (LacZ)'], ['Ses1 - 3 (LacZ)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Bottom 50% Responsivity (KRAB Cells)'])

a=subplot(2,2,2)
h=cdfplot(arc_red_top_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(arc_red_top_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
l=cdfplot(lacz_red_top_d1_d2_pref_d_all);
set(l, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
m=cdfplot(lacz_red_top_d1_d3_pref_d_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Ses1 - 2 (Arc Enh)'], ['Ses1 - 3 (Arc Enh)'], ['Ses1 - 2 (LacZ)'], ['Ses1 - 3 (LacZ)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Top 50% Responsivity (KRAB Cells)'])


a=subplot(2,2,3)
h=cdfplot(arc_green_bot_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(arc_green_bot_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(lacz_green_bot_d1_d2_pref_d_all);
set(l, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
m=cdfplot(lacz_green_bot_d1_d3_pref_d_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Ses1 - 2 (Arc Enh)'], ['Ses1 - 3 (Arc Enh)'], ['Ses1 - 2 (LacZ)'], ['Ses1 - 3 (LacZ)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Bottom 50% Responsivity (Non-KRAB Cells)'])

a=subplot(2,2,4)
h=cdfplot(arc_green_top_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(arc_green_top_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(lacz_green_top_d1_d2_pref_d_all);
set(l, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
m=cdfplot(lacz_green_top_d1_d3_pref_d_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Ses1 - 2 (Arc Enh)'], ['Ses1 - 3 (Arc Enh)'], ['Ses1 - 2 (LacZ)'], ['Ses1 - 3 (LacZ)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Top 50% Responsivity (Non-KRAB Cells)'])

print(fullfile(realfnout, ['grace_allmice_medsplit_changeprefcdf_green_red_.pdf']), '-dpdf', '-bestfit')

	%%
fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(arc_red_bot_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(arc_red_top_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
l=cdfplot(lacz_red_bot_pref_d_all);
set(l, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
m=cdfplot(lacz_red_top_pref_d_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Bottom 50% (Arc)'], ['Top 50% (Arc)'], ['Bottom 50% (LacZ)'], ['Top 50% (LacZ)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['KRAB Cells'])
hold off
print(fullfile(realfnout, ['grace_allmice_medsplit_changeprefcdf_ALLDAYS_red_.pdf']), '-dpdf', '-bestfit')

%%		
fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(arc_green_bot_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(arc_green_top_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(lacz_green_bot_pref_d_all);
set(l, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
m=cdfplot(lacz_green_top_pref_d_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Bottom 50% (Arc)'], ['Top 50% (Arc)'], ['Bottom 50% (LacZ)'], ['Top 50% (LacZ)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Non-KRAB Cells'])
hold off
print(fullfile(realfnout, ['grace_allmice_medsplit_changeprefcdf_ALLDAYS_green_.pdf']), '-dpdf', '-bestfit')


%%
fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(arc_red_wide_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(arc_red_sharp_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
l=cdfplot(lacz_red_wide_pref_d_all);
set(l, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
m=cdfplot(lacz_red_sharp_pref_d_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Wide (Arc)'], ['Sharp (Arc)'], ['Wide (LacZ)'], ['Sharp (LacZ)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['KRAB Cells'])
hold off
print(fullfile(realfnout, ['grace_allmice_wide_sharp_changeprefcdf_ALLDAYS_red_.pdf']), '-dpdf', '-bestfit')

%%

fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(arc_red_max_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(lacz_red_max_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['All Sessions (Arc Enh)'], ['All Sessions (LacZ)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F value'])
ylabel(['% of cells'])
title(['KRAB Cells'])
hold off
print(fullfile(realfnout, ['grace_allmice_kallsessions_red_.pdf']), '-dpdf', '-bestfit')


fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(arc_green_max_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(lacz_green_max_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
legend(['All Sessions (Arc Enh)'], ['All Sessions (LacZ)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F value'])
ylabel(['% of cells'])
title(['Non-KRAB Cells'])
hold off
print(fullfile(realfnout, ['grace_allmice_kallsessions_green_.pdf']), '-dpdf', '-bestfit')


%%

fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(arc_red_k_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(lacz_red_k_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['All Sessions (Arc Enh)'], ['All Sessions (LacZ)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['k value'])
ylabel(['% of cells'])
title(['KRAB Cells'])
hold off
print(fullfile(realfnout, ['grace_allmice_maxallsessions_red_.pdf']), '-dpdf', '-bestfit')


fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(arc_green_k_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(lacz_green_k_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
legend(['All Sessions (Arc Enh)'], ['All Sessions (LacZ)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['k value'])
ylabel(['% of cells'])
title(['Non-KRAB Cells'])
hold off
print(fullfile(realfnout, ['grace_allmice_maxallsessions_green_.pdf']), '-dpdf', '-bestfit')

%%
fig = figure;
% sgtitle('reliability values')
% a=subplot(2,2,1)
h=cdfplot(arc_red_reliability_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(lacz_red_reliability_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['All Sessions (Arc Prom)'], ['All Sessions (LacZ)'], 'Location', 'southeast')
% xlim([0 30])
xlabel(['Reliability value'])
ylabel(['% of cells'])
title(['KRAB Cells'])
hold off
print(fullfile(realfnout, ['grace_allmice_reliabilityallsessions_red_.pdf']), '-dpdf', '-bestfit')


fig = figure;
% sgtitle('reliability values')
% a=subplot(2,2,1)
h=cdfplot(arc_green_reliability_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(lacz_green_reliability_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
legend(['All Sessions (Arc Prom)'], ['All Sessions (LacZ)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['Reliability value'])
ylabel(['% of cells'])
title(['Non-KRAB Cells'])
hold off
print(fullfile(realfnout, ['grace_allmice_reliabilityallsessions_green_.pdf']), '-dpdf', '-bestfit')

%% all pref ori scores concat

arc_red_all_pref_d = [arc_red_d1_d2_pref_d_all arc_red_d1_d3_pref_d_all];
arc_green_all_pref_d = [arc_green_d1_d2_pref_d_all arc_green_d1_d3_pref_d_all];
lacz_red_all_pref_d = [lacz_red_d1_d2_pref_d_all lacz_red_d1_d3_pref_d_all];
lacz_green_all_pref_d = [lacz_green_d1_d2_pref_d_all lacz_green_d1_d3_pref_d_all];

fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(arc_red_all_pref_d);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(lacz_red_all_pref_d);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Arc Enh (KRAB)'], ['LacZ (KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['KRAB Cells'])
hold off
print(fullfile(realfnout, ['grace_allmice_changepref_red_.pdf']), '-dpdf', '-bestfit')


fig = figure;
l=cdfplot(arc_green_all_pref_d);
set(l, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
m=cdfplot(lacz_green_all_pref_d);
set(m, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
legend(['Arc Enh (Non-KRAB)'], ['LacZ (Non-KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Non-KRAB Cells'])
hold off
print(fullfile(realfnout, ['grace_allmice_changepref_green_.pdf']), '-dpdf', '-bestfit')


%%
%%compare all groups from mine and grace's
%%
%k and max change scores
arc_red_d1_d2_k_all = abs(arc_red_d1_k_all-arc_red_d2_k_all);
arc_red_d1_d3_k_all = abs(arc_red_d1_k_all-arc_red_d3_k_all);
arc_red_d_k_all = [arc_red_d1_d2_k_all arc_red_d1_d3_k_all];

arc_green_d1_d2_k_all = abs(arc_green_d1_k_all-arc_green_d2_k_all);
arc_green_d1_d3_k_all = abs(arc_green_d1_k_all-arc_green_d3_k_all);
arc_green_d_k_all = [arc_green_d1_d2_k_all arc_green_d1_d3_k_all];

arc_red_d1_d2_max_all = abs(arc_red_d1_max_all-arc_red_d2_max_all);
arc_red_d1_d3_max_all = abs(arc_red_d1_max_all-arc_red_d3_max_all);
arc_red_d_max_all = [arc_red_d1_d2_max_all arc_red_d1_d3_max_all];

arc_green_d1_d2_max_all = abs(arc_green_d1_max_all-arc_green_d2_max_all);
arc_green_d1_d3_max_all = abs(arc_green_d1_max_all-arc_green_d3_max_all);
arc_green_d_max_all = [arc_green_d1_d2_max_all arc_green_d1_d3_max_all];


lacz_red_d1_d2_k_all = abs(lacz_red_d1_k_all-lacz_red_d2_k_all);
lacz_red_d1_d3_k_all = abs(lacz_red_d1_k_all-lacz_red_d3_k_all);
lacz_red_d_k_all = [lacz_red_d1_d2_k_all lacz_red_d1_d3_k_all];

lacz_green_d1_d2_k_all = abs(lacz_green_d1_k_all-lacz_green_d2_k_all);
lacz_green_d1_d3_k_all = abs(lacz_green_d1_k_all-lacz_green_d3_k_all);
lacz_green_d_k_all = [lacz_green_d1_d2_k_all lacz_green_d1_d3_k_all];

lacz_red_d1_d2_max_all = abs(lacz_red_d1_max_all-lacz_red_d2_max_all);
lacz_red_d1_d3_max_all = abs(lacz_red_d1_max_all-lacz_red_d3_max_all);
lacz_red_d_max_all = [lacz_red_d1_d2_max_all lacz_red_d1_d3_max_all];

lacz_green_d1_d2_max_all = abs(lacz_green_d1_max_all-lacz_green_d2_max_all);
lacz_green_d1_d3_max_all = abs(lacz_green_d1_max_all-lacz_green_d3_max_all);
lacz_green_d_max_all = [lacz_green_d1_d2_max_all lacz_green_d1_d3_max_all];


%%

fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(arc_red_d_k_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(lacz_red_d_k_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Arc Mice (KRAB)'], ['LacZ Mice (KRAB)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['Change in k Value'])
ylabel(['% of cells'])
title(['KRAB Cells'])
hold off
print(fullfile(realfnout, ['grace_allmice_changek_red.pdf']), '-dpdf', '-bestfit')

fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(arc_red_d_max_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(lacz_red_d_max_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Arc Mice (KRAB)'], ['LacZ Mice (KRAB)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Change in Max dF/F Value'])
ylabel(['% of cells'])
title(['KRAB Cells'])
hold off
print(fullfile(realfnout, ['grace_allmice_changemax_red.pdf']), '-dpdf', '-bestfit')


%%

save(fullfile(realfnout, ['_grace_pooled_pref_data.mat']), 'arc_red_d1_d2_pref_d_all', 'arc_red_d1_d3_pref_d_all', 'arc_green_d1_d2_pref_d_all', 'arc_green_d1_d3_pref_d_all', ...
    'arc_red_all_pref_d', 'arc_green_all_pref_d', 'lacz_red_d1_d2_pref_d_all', 'lacz_red_d1_d3_pref_d_all', 'lacz_green_d1_d2_pref_d_all', 'lacz_green_d1_d3_pref_d_all', ...
    'lacz_red_all_pref_d', 'lacz_green_all_pref_d')

save(fullfile(realfnout, ['_grace_pooled_kmax_data.mat']), 'arc_red_k_all', 'arc_green_k_all', 'lacz_red_k_all', 'lacz_green_k_all', 'arc_red_d_k_all', 'arc_green_d_k_all', ...
    'lacz_red_d_k_all', 'lacz_green_d_k_all', 'arc_red_max_all', 'arc_green_max_all', 'lacz_red_max_all', 'lacz_green_k_all', ...
    'arc_red_d_max_all', 'arc_green_d_max_all', 'lacz_red_d_max_all', 'lacz_green_d_k_all')

save(fullfile(realfnout, ['_grace_pooled_kmax_day1_data.mat']), 'arc_red_d1_k_all', 'arc_green_d1_k_all', 'lacz_red_d1_k_all', 'lacz_green_d1_k_all', ...
    'arc_red_d1_max_all', 'arc_green_d1_max_all', 'lacz_red_d1_max_all', 'lacz_green_d1_max_all')

save(fullfile(realfnout, ['_grace_medsplit_pref_data.mat']), 'arc_red_bot_d1_d2_pref_d_all', 'arc_red_bot_d1_d3_pref_d_all', 'arc_green_bot_d1_d2_pref_d_all', 'arc_green_bot_d1_d3_pref_d_all', ...
    'arc_red_top_d1_d2_pref_d_all', 'arc_red_top_d1_d3_pref_d_all', 'arc_green_top_d1_d2_pref_d_all', 'arc_green_top_d1_d3_pref_d_all', ...
    'arc_red_bot_pref_d_all', 'arc_red_top_pref_d_all', 'arc_green_bot_pref_d_all', 'arc_green_top_pref_d_all', ...
    'lacz_red_bot_d1_d2_pref_d_all', 'lacz_red_bot_d1_d3_pref_d_all', 'lacz_green_bot_d1_d2_pref_d_all', 'lacz_green_bot_d1_d3_pref_d_all', ...
    'lacz_red_top_d1_d2_pref_d_all', 'lacz_red_top_d1_d3_pref_d_all', 'lacz_green_top_d1_d2_pref_d_all', 'lacz_green_top_d1_d3_pref_d_all', ...
    'lacz_red_bot_pref_d_all', 'lacz_red_top_pref_d_all', 'lacz_green_bot_pref_d_all', 'lacz_green_top_pref_d_all')


save(fullfile(realfnout, ['_grace_k_medsplit_pref_data.mat']), 'arc_red_wide_d1_d2_pref_d_all', 'arc_red_wide_d1_d3_pref_d_all', 'arc_red_wide_pref_d_all', 'arc_green_wide_d1_d2_pref_d_all', ...
    'arc_green_wide_d1_d3_pref_d_all', 'arc_green_wide_pref_d_all', 'arc_red_sharp_d1_d2_pref_d_all', 'arc_red_sharp_d1_d3_pref_d_all', ...
    'arc_red_sharp_pref_d_all', 'arc_green_sharp_d1_d2_pref_d_all', 'arc_green_sharp_d1_d3_pref_d_all', 'arc_green_sharp_pref_d_all', ...
    'lacz_red_wide_d1_d2_pref_d_all', 'lacz_red_wide_d1_d3_pref_d_all', 'lacz_red_wide_pref_d_all', 'lacz_green_wide_d1_d2_pref_d_all', ...
    'lacz_green_wide_d1_d3_pref_d_all', 'lacz_green_wide_pref_d_all', 'lacz_red_sharp_d1_d2_pref_d_all', 'lacz_red_sharp_d1_d3_pref_d_all', ...
    'lacz_red_sharp_pref_d_all', 'lacz_green_sharp_d1_d2_pref_d_all', 'lacz_green_sharp_d1_d3_pref_d_all', 'lacz_green_sharp_pref_d_all')

%%
figure;
subplot(2,2,1)
scatter(lacz_green_k_all, lacz_green_max_all, [], 'green')
xlabel(['k Value'])
ylabel(['Max dF/F Value'])
title(['All Green Cells'], ['r = ' num2str(corr(lacz_green_k_all', lacz_green_max_all'))])

% figure;
subplot(2,2,2)
scatter(lacz_red_k_all, lacz_red_max_all, [], 'red')
xlabel(['k Value'])
ylabel(['Max dF/F Value'])
title(['All Red Cells'], ['r = ' num2str(corr(lacz_red_k_all', lacz_red_max_all'))])


lacz_green_notmaxk = find(lacz_green_k_all <= 29);
lacz_green_maxk = find(not(lacz_green_k_all <= 29));
lacz_red_notmaxk = find(lacz_red_k_all <= 29);
lacz_red_maxk = find(not(lacz_red_k_all <= 29));


% figure;
subplot(2,2,3)
scatter(lacz_green_k_all(lacz_green_notmaxk), lacz_green_max_all(lacz_green_notmaxk), [], 'green')
xlabel(['k Value'])
ylabel(['Max dF/F Value'])
title(['Green Cells < 30 k'], ['r = ' num2str(corr(lacz_green_k_all(lacz_green_notmaxk)', lacz_green_max_all(lacz_green_notmaxk)'))])


% figure;
subplot(2,2,4)
scatter(lacz_red_k_all(lacz_red_notmaxk), lacz_red_max_all(lacz_red_notmaxk), [], 'red')
xlabel(['k Value'])
ylabel(['Max dF/F Value'])
title(['Red Cells < 30 k'], ['r = ' num2str(corr(lacz_red_k_all(lacz_red_notmaxk)', lacz_red_max_all(lacz_red_notmaxk)'))])

print(fullfile(realfnout, ['laczmice_k_max_scatter_allcells.pdf']), '-dpdf', '-bestfit')


%continue this - i was trying to explore the values of k that are and are not at max k value***
%maybe look at distribution of maxed out k value cells too*

%%
figure;
subplot(1,2,1)
hist(lacz_green_max_all(lacz_green_maxk))
subplot(1,2,2)
hist(lacz_red_max_all(lacz_red_maxk))

figure;
cdfplot(lacz_green_max_all(lacz_green_maxk))
hold on
cdfplot(lacz_red_max_all(lacz_red_maxk))
hold off
xlim([0 1])

%%
%binning fit reliability
arc_red_grand_reliability = max([arc_red_d1_reliability_all; arc_red_d2_reliability_all; arc_red_d3_reliability_all]);
lacz_red_grand_reliability = max([lacz_red_d1_reliability_all; lacz_red_d2_reliability_all; lacz_red_d3_reliability_all]);

arc_red_bin1_reliability_ind = find(arc_red_grand_reliability >=0 & arc_red_grand_reliability <=30);
arc_red_bin2_reliability_ind = find(arc_red_grand_reliability >=31 & arc_red_grand_reliability <=60);
arc_red_bin3_reliability_ind = find(arc_red_grand_reliability >=61 & arc_red_grand_reliability <=90);

arc_red_d1_d2_pref_bin1_reliability = arc_red_d1_d2_pref_d_all(arc_red_bin1_reliability_ind);
arc_red_d1_d2_pref_bin2_reliability = arc_red_d1_d2_pref_d_all(arc_red_bin2_reliability_ind);
arc_red_d1_d2_pref_bin3_reliability = arc_red_d1_d2_pref_d_all(arc_red_bin3_reliability_ind);

arc_red_d1_d3_pref_bin1_reliability = arc_red_d1_d3_pref_d_all(arc_red_bin1_reliability_ind);
arc_red_d1_d3_pref_bin2_reliability = arc_red_d1_d3_pref_d_all(arc_red_bin2_reliability_ind);
arc_red_d1_d3_pref_bin3_reliability = arc_red_d1_d3_pref_d_all(arc_red_bin3_reliability_ind);

arc_red_pref_bin1_reliability_all = [arc_red_d1_d2_pref_bin1_reliability arc_red_d1_d3_pref_bin1_reliability];
arc_red_pref_bin2_reliability_all = [arc_red_d1_d2_pref_bin2_reliability arc_red_d1_d3_pref_bin2_reliability];
arc_red_pref_bin3_reliability_all = [arc_red_d1_d2_pref_bin3_reliability arc_red_d1_d3_pref_bin3_reliability];

lacz_red_bin1_reliability_ind = find(lacz_red_grand_reliability >=0 & lacz_red_grand_reliability <=30);
lacz_red_bin2_reliability_ind = find(lacz_red_grand_reliability >=31 & lacz_red_grand_reliability <=60);
lacz_red_bin3_reliability_ind = find(lacz_red_grand_reliability >=61 & lacz_red_grand_reliability <=90);

lacz_red_d1_d2_pref_bin1_reliability = lacz_red_d1_d2_pref_d_all(lacz_red_bin1_reliability_ind);
lacz_red_d1_d2_pref_bin2_reliability = lacz_red_d1_d2_pref_d_all(lacz_red_bin2_reliability_ind);
lacz_red_d1_d2_pref_bin3_reliability = lacz_red_d1_d2_pref_d_all(lacz_red_bin3_reliability_ind);

lacz_red_d1_d3_pref_bin1_reliability = lacz_red_d1_d3_pref_d_all(lacz_red_bin1_reliability_ind);
lacz_red_d1_d3_pref_bin2_reliability = lacz_red_d1_d3_pref_d_all(lacz_red_bin2_reliability_ind);
lacz_red_d1_d3_pref_bin3_reliability = lacz_red_d1_d3_pref_d_all(lacz_red_bin3_reliability_ind);

lacz_red_pref_bin1_reliability_all = [lacz_red_d1_d2_pref_bin1_reliability lacz_red_d1_d3_pref_bin1_reliability];
lacz_red_pref_bin2_reliability_all = [lacz_red_d1_d2_pref_bin2_reliability lacz_red_d1_d3_pref_bin2_reliability];
lacz_red_pref_bin3_reliability_all = [lacz_red_d1_d2_pref_bin3_reliability lacz_red_d1_d3_pref_bin3_reliability];

%%
%binning fit reliability
arc_green_grand_reliability = max([arc_green_d1_reliability_all; arc_green_d2_reliability_all; arc_green_d3_reliability_all]);
lacz_green_grand_reliability = max([lacz_green_d1_reliability_all; lacz_green_d2_reliability_all; lacz_green_d3_reliability_all]);

arc_green_bin1_reliability_ind = find(arc_green_grand_reliability >=0 & arc_green_grand_reliability <=30);
arc_green_bin2_reliability_ind = find(arc_green_grand_reliability >=31 & arc_green_grand_reliability <=60);
arc_green_bin3_reliability_ind = find(arc_green_grand_reliability >=61 & arc_green_grand_reliability <=90);

arc_green_d1_d2_pref_bin1_reliability = arc_green_d1_d2_pref_d_all(arc_green_bin1_reliability_ind);
arc_green_d1_d2_pref_bin2_reliability = arc_green_d1_d2_pref_d_all(arc_green_bin2_reliability_ind);
arc_green_d1_d2_pref_bin3_reliability = arc_green_d1_d2_pref_d_all(arc_green_bin3_reliability_ind);

arc_green_d1_d3_pref_bin1_reliability = arc_green_d1_d3_pref_d_all(arc_green_bin1_reliability_ind);
arc_green_d1_d3_pref_bin2_reliability = arc_green_d1_d3_pref_d_all(arc_green_bin2_reliability_ind);
arc_green_d1_d3_pref_bin3_reliability = arc_green_d1_d3_pref_d_all(arc_green_bin3_reliability_ind);

arc_green_pref_bin1_reliability_all = [arc_green_d1_d2_pref_bin1_reliability arc_green_d1_d3_pref_bin1_reliability];
arc_green_pref_bin2_reliability_all = [arc_green_d1_d2_pref_bin2_reliability arc_green_d1_d3_pref_bin2_reliability];
arc_green_pref_bin3_reliability_all = [arc_green_d1_d2_pref_bin3_reliability arc_green_d1_d3_pref_bin3_reliability];

lacz_green_bin1_reliability_ind = find(lacz_green_grand_reliability >=0 & lacz_green_grand_reliability <=30);
lacz_green_bin2_reliability_ind = find(lacz_green_grand_reliability >=31 & lacz_green_grand_reliability <=60);
lacz_green_bin3_reliability_ind = find(lacz_green_grand_reliability >=61 & lacz_green_grand_reliability <=90);

lacz_green_d1_d2_pref_bin1_reliability = lacz_green_d1_d2_pref_d_all(lacz_green_bin1_reliability_ind);
lacz_green_d1_d2_pref_bin2_reliability = lacz_green_d1_d2_pref_d_all(lacz_green_bin2_reliability_ind);
lacz_green_d1_d2_pref_bin3_reliability = lacz_green_d1_d2_pref_d_all(lacz_green_bin3_reliability_ind);

lacz_green_d1_d3_pref_bin1_reliability = lacz_green_d1_d3_pref_d_all(lacz_green_bin1_reliability_ind);
lacz_green_d1_d3_pref_bin2_reliability = lacz_green_d1_d3_pref_d_all(lacz_green_bin2_reliability_ind);
lacz_green_d1_d3_pref_bin3_reliability = lacz_green_d1_d3_pref_d_all(lacz_green_bin3_reliability_ind);

lacz_green_pref_bin1_reliability_all = [lacz_green_d1_d2_pref_bin1_reliability lacz_green_d1_d3_pref_bin1_reliability];
lacz_green_pref_bin2_reliability_all = [lacz_green_d1_d2_pref_bin2_reliability lacz_green_d1_d3_pref_bin2_reliability];
lacz_green_pref_bin3_reliability_all = [lacz_green_d1_d2_pref_bin3_reliability lacz_green_d1_d3_pref_bin3_reliability];



%%

figure;
a = errorbar([mean(arc_red_pref_bin1_reliability_all), mean(arc_red_pref_bin2_reliability_all), mean(arc_red_pref_bin3_reliability_all)], ...
    [std(arc_red_pref_bin1_reliability_all)/sqrt(length(arc_red_pref_bin1_reliability_all)), ... 
    std(arc_red_pref_bin2_reliability_all)/sqrt(length(arc_red_pref_bin2_reliability_all)), ... 
    std(arc_red_pref_bin3_reliability_all)/sqrt(length(arc_red_pref_bin3_reliability_all))]);

hold on

b = errorbar([mean(lacz_red_pref_bin1_reliability_all), mean(lacz_red_pref_bin2_reliability_all), mean(lacz_red_pref_bin3_reliability_all)], ...
    [std(lacz_red_pref_bin1_reliability_all)/sqrt(length(lacz_red_pref_bin1_reliability_all)), ... 
    std(lacz_red_pref_bin2_reliability_all)/sqrt(length(lacz_red_pref_bin2_reliability_all)), ... 
    std(lacz_red_pref_bin3_reliability_all)/sqrt(length(lacz_red_pref_bin3_reliability_all))]);

hold off


a.Marker = 'o';
a.MarkerFaceColor = [1 0 0];
a.Color = [1 0 0];
a.LineStyle = 'none';

b.Marker = 'square';
b.MarkerFaceColor = [.75 0 0];
b.Color = [.75 0 0];
b.LineStyle = 'none';

xlim([0.5 3.5])
ylim([0 45])
xticklabels({'','0-30', '', '31-60', '', '61-90'})
xlabel('Fit Reliability (degrees)')
ylabel('Mean Pref Ori Change (+- s.e.m.)')
legend(['Arc (KRAB)'], ['LacZ (KRAB)'])

%%
figure;
c = errorbar([mean(arc_green_pref_bin1_reliability_all), mean(arc_green_pref_bin2_reliability_all), mean(arc_green_pref_bin3_reliability_all)], ...
    [std(arc_green_pref_bin1_reliability_all)/sqrt(length(arc_green_pref_bin1_reliability_all)), ... 
    std(arc_green_pref_bin2_reliability_all)/sqrt(length(arc_green_pref_bin2_reliability_all)), ... 
    std(arc_green_pref_bin3_reliability_all)/sqrt(length(arc_green_pref_bin3_reliability_all))]);

hold on

d = errorbar([mean(lacz_green_pref_bin1_reliability_all), mean(lacz_green_pref_bin2_reliability_all), mean(lacz_green_pref_bin3_reliability_all)], ...
    [std(lacz_green_pref_bin1_reliability_all)/sqrt(length(lacz_green_pref_bin1_reliability_all)), ... 
    std(lacz_green_pref_bin2_reliability_all)/sqrt(length(lacz_green_pref_bin2_reliability_all)), ... 
    std(lacz_green_pref_bin3_reliability_all)/sqrt(length(lacz_green_pref_bin3_reliability_all))]);
c.Marker = 'o';
c.MarkerFaceColor = [0 1 0];
c.Color = [0 1 0];
c.LineStyle = 'none';

d.Marker = 'square';
d.MarkerFaceColor = [0 .75 0];
d.Color = [0 .75 0];
d.LineStyle = 'none';

xlim([0.5 3.5])
ylim([0 45])
xticklabels({'','0-30', '', '31-60', '', '61-90'})
xlabel('Fit Reliability (degrees)')
ylabel('Mean Pref Ori Change (+- s.e.m.)')
legend(['Arc (Non-KRAB)'], ['LacZ (Non-KRAB)'])

%%
save(fullfile(realfnout, ['_grace_reliabilitysplit_pref_data.mat']), 'arc_red_pref_bin1_reliability_all', 'arc_red_pref_bin2_reliability_all', 'arc_red_pref_bin3_reliability_all', ...
    'lacz_red_pref_bin1_reliability_all', 'lacz_red_pref_bin2_reliability_all', 'lacz_red_pref_bin3_reliability_all', ...
    'arc_green_pref_bin1_reliability_all', 'arc_green_pref_bin2_reliability_all', 'arc_green_pref_bin3_reliability_all', ...
    'lacz_green_pref_bin1_reliability_all', 'lacz_green_pref_bin2_reliability_all', 'lacz_green_pref_bin3_reliability_all')

