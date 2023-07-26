%CLEAR EVERYTHING
clear all
clear all global
clc
close all



%%
%LOAD DATA AND IDENTIFY FOLDERS TO SAVE TO
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_Arc_greenVred'; %folder to save files to
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_pooled_Arc_greenVred'; %folder to save files to
dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%
lacz_maxsep_reliability_scores = [];
lacz_ksep_reliability_scores = [];

% list = [1 4 37 40];
list = [1 4 37 40 49 55 61];
for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};

    lacz_maxsep_reliability_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_max_median_sep_reliability']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_max_median_sep_reliability_2.mat']));
        lacz_maxsep_reliability_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_max_median_sep_reliability_2']));
    else
        lacz_maxsep_reliability_ses2 = [];
    end
    lacz_maxsep_reliability_ses_all = [lacz_maxsep_reliability_ses1 lacz_maxsep_reliability_ses2];
    lacz_maxsep_reliability_scores = [lacz_maxsep_reliability_scores lacz_maxsep_reliability_ses_all];

    lacz_ksep_reliability_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_k_median_sep_reliability']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_k_median_sep_reliability_2.mat']));
        lacz_ksep_reliability_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_k_median_sep_reliability_2']));
    else
        lacz_ksep_reliability_ses2 = [];
    end
    lacz_ksep_reliability_ses_all = [lacz_ksep_reliability_ses1 lacz_ksep_reliability_ses2];
    lacz_ksep_reliability_scores = [lacz_ksep_reliability_scores lacz_ksep_reliability_ses_all];


end

%%

lacz_green_sharp_d1_reliability_all = [];
lacz_green_sharp_d2_reliability_all = [];
lacz_green_sharp_d3_reliability_all = [];

lacz_green_wide_d1_reliability_all = [];
lacz_green_wide_d2_reliability_all = [];
lacz_green_wide_d3_reliability_all = [];

lacz_red_sharp_d1_reliability_all = [];
lacz_red_sharp_d2_reliability_all = [];
lacz_red_sharp_d3_reliability_all = [];

lacz_red_wide_d1_reliability_all = [];
lacz_red_wide_d2_reliability_all = [];
lacz_red_wide_d3_reliability_all = [];

lacz_green_top_50_d1_reliability_all = [];
lacz_green_top_50_d2_reliability_all = [];
lacz_green_top_50_d3_reliability_all = [];

lacz_green_bot_50_d1_reliability_all = [];
lacz_green_bot_50_d2_reliability_all = [];
lacz_green_bot_50_d3_reliability_all = [];

lacz_red_top_50_d1_reliability_all = [];
lacz_red_top_50_d2_reliability_all = [];
lacz_red_top_50_d3_reliability_all = [];

lacz_red_bot_50_d1_reliability_all = [];
lacz_red_bot_50_d2_reliability_all = [];
lacz_red_bot_50_d3_reliability_all = [];





for idata = 1:length(lacz_ksep_reliability_scores)
    lacz_green_sharp_d1_reliability = lacz_ksep_reliability_scores(idata).green_sharp_reliability_d1;
    lacz_green_sharp_d1_reliability_all = [lacz_green_sharp_d1_reliability_all lacz_green_sharp_d1_reliability];
    lacz_green_sharp_d2_reliability = lacz_ksep_reliability_scores(idata).green_sharp_reliability_d2;
    lacz_green_sharp_d2_reliability_all = [lacz_green_sharp_d2_reliability_all lacz_green_sharp_d2_reliability];
    lacz_green_sharp_d3_reliability = lacz_ksep_reliability_scores(idata).green_sharp_reliability_d3;
    lacz_green_sharp_d3_reliability_all = [lacz_green_sharp_d3_reliability_all lacz_green_sharp_d3_reliability];
    
    lacz_green_wide_d1_reliability = lacz_ksep_reliability_scores(idata).green_wide_reliability_d1;
    lacz_green_wide_d1_reliability_all = [lacz_green_wide_d1_reliability_all lacz_green_wide_d1_reliability];
    lacz_green_wide_d2_reliability = lacz_ksep_reliability_scores(idata).green_wide_reliability_d2;
    lacz_green_wide_d2_reliability_all = [lacz_green_wide_d2_reliability_all lacz_green_wide_d2_reliability];
    lacz_green_wide_d3_reliability = lacz_ksep_reliability_scores(idata).green_wide_reliability_d3;
    lacz_green_wide_d3_reliability_all = [lacz_green_wide_d3_reliability_all lacz_green_wide_d3_reliability];
    
    lacz_red_sharp_d1_reliability = lacz_ksep_reliability_scores(idata).red_sharp_reliability_d1;
    lacz_red_sharp_d1_reliability_all = [lacz_red_sharp_d1_reliability_all lacz_red_sharp_d1_reliability];
    lacz_red_sharp_d2_reliability = lacz_ksep_reliability_scores(idata).red_sharp_reliability_d2;
    lacz_red_sharp_d2_reliability_all = [lacz_red_sharp_d2_reliability_all lacz_red_sharp_d2_reliability];
    lacz_red_sharp_d3_reliability = lacz_ksep_reliability_scores(idata).red_sharp_reliability_d3;
    lacz_red_sharp_d3_reliability_all = [lacz_red_sharp_d3_reliability_all lacz_red_sharp_d3_reliability];
    
    lacz_red_wide_d1_reliability = lacz_ksep_reliability_scores(idata).red_wide_reliability_d1;
    lacz_red_wide_d1_reliability_all = [lacz_red_wide_d1_reliability_all lacz_red_wide_d1_reliability];
    lacz_red_wide_d2_reliability = lacz_ksep_reliability_scores(idata).red_wide_reliability_d2;
    lacz_red_wide_d2_reliability_all = [lacz_red_wide_d2_reliability_all lacz_red_wide_d2_reliability];
    lacz_red_wide_d3_reliability = lacz_ksep_reliability_scores(idata).red_wide_reliability_d3;
    lacz_red_wide_d3_reliability_all = [lacz_red_wide_d3_reliability_all lacz_red_wide_d3_reliability];


    lacz_green_top_50_d1_reliability = lacz_maxsep_reliability_scores(idata).green_top_50_reliability_d1;
    lacz_green_top_50_d1_reliability_all = [lacz_green_top_50_d1_reliability_all lacz_green_top_50_d1_reliability];
    lacz_green_top_50_d2_reliability = lacz_maxsep_reliability_scores(idata).green_top_50_reliability_d2;
    lacz_green_top_50_d2_reliability_all = [lacz_green_top_50_d2_reliability_all lacz_green_top_50_d2_reliability];
    lacz_green_top_50_d3_reliability = lacz_maxsep_reliability_scores(idata).green_top_50_reliability_d3;
    lacz_green_top_50_d3_reliability_all = [lacz_green_top_50_d3_reliability_all lacz_green_top_50_d3_reliability];
    
    lacz_green_bot_50_d1_reliability = lacz_maxsep_reliability_scores(idata).green_bot_50_reliability_d1;
    lacz_green_bot_50_d1_reliability_all = [lacz_green_bot_50_d1_reliability_all lacz_green_bot_50_d1_reliability];
    lacz_green_bot_50_d2_reliability = lacz_maxsep_reliability_scores(idata).green_bot_50_reliability_d2;
    lacz_green_bot_50_d2_reliability_all = [lacz_green_bot_50_d2_reliability_all lacz_green_bot_50_d2_reliability];
    lacz_green_bot_50_d3_reliability = lacz_maxsep_reliability_scores(idata).green_bot_50_reliability_d3;
    lacz_green_bot_50_d3_reliability_all = [lacz_green_bot_50_d3_reliability_all lacz_green_bot_50_d3_reliability];
    
    lacz_red_top_50_d1_reliability = lacz_maxsep_reliability_scores(idata).red_top_50_reliability_d1;
    lacz_red_top_50_d1_reliability_all = [lacz_red_top_50_d1_reliability_all lacz_red_top_50_d1_reliability];
    lacz_red_top_50_d2_reliability = lacz_maxsep_reliability_scores(idata).red_top_50_reliability_d2;
    lacz_red_top_50_d2_reliability_all = [lacz_red_top_50_d2_reliability_all lacz_red_top_50_d2_reliability];
    lacz_red_top_50_d3_reliability = lacz_maxsep_reliability_scores(idata).red_top_50_reliability_d3;
    lacz_red_top_50_d3_reliability_all = [lacz_red_top_50_d3_reliability_all lacz_red_top_50_d3_reliability];
    
    lacz_red_bot_50_d1_reliability = lacz_maxsep_reliability_scores(idata).red_bot_50_reliability_d1;
    lacz_red_bot_50_d1_reliability_all = [lacz_red_bot_50_d1_reliability_all lacz_red_bot_50_d1_reliability];
    lacz_red_bot_50_d2_reliability = lacz_maxsep_reliability_scores(idata).red_bot_50_reliability_d2;
    lacz_red_bot_50_d2_reliability_all = [lacz_red_bot_50_d2_reliability_all lacz_red_bot_50_d2_reliability];
    lacz_red_bot_50_d3_reliability = lacz_maxsep_reliability_scores(idata).red_bot_50_reliability_d3;
    lacz_red_bot_50_d3_reliability_all = [lacz_red_bot_50_d3_reliability_all lacz_red_bot_50_d3_reliability];

end

%%

lacz_green_sharp_reliability_all = [lacz_green_sharp_d1_reliability_all lacz_green_sharp_d2_reliability_all lacz_green_sharp_d3_reliability_all];
lacz_green_wide_reliability_all = [lacz_green_wide_d1_reliability_all lacz_green_wide_d2_reliability_all lacz_green_wide_d3_reliability_all];
lacz_red_sharp_reliability_all = [lacz_red_sharp_d1_reliability_all lacz_red_sharp_d2_reliability_all lacz_red_sharp_d3_reliability_all];
lacz_red_wide_reliability_all = [lacz_red_wide_d1_reliability_all lacz_red_wide_d2_reliability_all lacz_red_wide_d3_reliability_all];

lacz_green_top_50_reliability_all = [lacz_green_top_50_d1_reliability_all lacz_green_top_50_d2_reliability_all lacz_green_top_50_d3_reliability_all];
lacz_green_bot_50_reliability_all = [lacz_green_bot_50_d1_reliability_all lacz_green_bot_50_d2_reliability_all lacz_green_bot_50_d3_reliability_all];
lacz_red_top_50_reliability_all = [lacz_red_top_50_d1_reliability_all lacz_red_top_50_d2_reliability_all lacz_red_top_50_d3_reliability_all];
lacz_red_bot_50_reliability_all = [lacz_red_bot_50_d1_reliability_all lacz_red_bot_50_d2_reliability_all lacz_red_bot_50_d3_reliability_all];

%%
fig = figure;
h=cdfplot(lacz_red_wide_reliability_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(lacz_red_sharp_reliability_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
l=cdfplot(lacz_green_wide_reliability_all);
set(l, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
m=cdfplot(lacz_green_sharp_reliability_all);
set(m, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
legend(['Wide (KRAB)'], ['Sharp (KRAB)'], ['Wide (Non-KRAB)'], ['Sharp (Non-KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Fit Reliability'])
ylabel(['% of cells'])
title(['LacZ Mice'])
hold off
print(fullfile(realfnout, ['ALLlaczmice_ksplit_reliability_allcells.pdf']), '-dpdf', '-bestfit')


fig = figure;
h=cdfplot(lacz_red_bot_50_reliability_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(lacz_red_top_50_reliability_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
l=cdfplot(lacz_green_bot_50_reliability_all);
set(l, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
m=cdfplot(lacz_green_top_50_reliability_all);
set(m, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
legend(['Bottom 50% (KRAB)'], ['Top 50% (KRAB)'], ['Bottom 50% (Non-KRAB)'], ['Top 50% (Non-KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Fit Reliability'])
ylabel(['% of cells'])
title(['LacZ Mice'])
hold off
print(fullfile(realfnout, ['ALLlaczmice_maxsplit_reliability_allcells.pdf']), '-dpdf', '-bestfit')







