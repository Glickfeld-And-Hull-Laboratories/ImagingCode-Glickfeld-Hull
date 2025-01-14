%clear everything
clear all
clear all global
clc
close all



%%
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\new_multi_day'; %folder to load files 
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\pooled\four_session_LM'; %folder to save files 
dataset = 'exp_list_darklight_actual_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%

d1_pref_all = [];
d2_pref_all = [];
d3_pref_all = [];
d4_pref_all = [];
d1_d2_pref_d_all = [];
d2_d3_pref_d_all = [];
d3_d4_pref_d_all = [];
d1_d3_pref_d_all = [];
d1_d4_pref_d_all = [];
d1_k_all = [];
d2_k_all = [];
d3_k_all = [];
d4_k_all = [];
d1_max_all = [];
d2_max_all = [];
d3_max_all = [];
d4_max_all = [];
d1_reliability_all = [];
d2_reliability_all = [];
d3_reliability_all = [];
d4_reliability_all = [];

top_d1_d2_pref_d_all = [];
top_d2_d3_pref_d_all = [];
top_d3_d4_pref_d_all = [];
top_d1_d3_pref_d_all = [];
top_d1_d4_pref_d_all = [];

bot_d1_d2_pref_d_all = [];
bot_d2_d3_pref_d_all = [];
bot_d3_d4_pref_d_all = [];
bot_d1_d3_pref_d_all = [];
bot_d1_d4_pref_d_all = [];

sharp_d1_d2_pref_d_all = [];
sharp_d2_d3_pref_d_all = [];
sharp_d3_d4_pref_d_all = [];
sharp_d1_d3_pref_d_all = [];
sharp_d1_d4_pref_d_all = [];

wide_d1_d2_pref_d_all = [];
wide_d2_d3_pref_d_all = [];
wide_d3_d4_pref_d_all = [];
wide_d1_d3_pref_d_all = [];
wide_d1_d4_pref_d_all = [];


list = [38 38+4 38+8];
for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    prefori_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes']));
    pref_dscores_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores']));
    k_max_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes']));
    median_sep_pref_dscores_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'median_sep_d_scores']));
    wide_sharp_sep_pref_dscores_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'wide_sharp_sep_d_scores']));
    reliability_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'fit_reliability']));

    d1_pref = prefori_ses_all.d1_prefori;
    d1_pref_all = [d1_pref_all d1_pref];
    d2_pref = prefori_ses_all.d2_prefori;
    d2_pref_all = [d2_pref_all d2_pref];
    d3_pref = prefori_ses_all.d3_prefori;
    d3_pref_all = [d3_pref_all d3_pref];
    d4_pref = prefori_ses_all.d4_prefori;
    d4_pref_all = [d4_pref_all d4_pref];

    d1_d2_pref_d = pref_dscores_ses_all.d_score_prefori_d1_d2;
    d1_d2_pref_d_all = [d1_d2_pref_d_all d1_d2_pref_d];
    d2_d3_pref_d = pref_dscores_ses_all.d_score_prefori_d2_d3;
    d2_d3_pref_d_all = [d2_d3_pref_d_all d2_d3_pref_d];
    d3_d4_pref_d = pref_dscores_ses_all.d_score_prefori_d3_d4;
    d3_d4_pref_d_all = [d3_d4_pref_d_all d3_d4_pref_d];

    d1_d3_pref_d = pref_dscores_ses_all.d_score_prefori_d1_d3;
    d1_d3_pref_d_all = [d1_d3_pref_d_all d1_d3_pref_d];
    d1_d4_pref_d = pref_dscores_ses_all.d_score_prefori_d1_d4;
    d1_d4_pref_d_all = [d1_d4_pref_d_all d1_d4_pref_d];

    top_d1_d2_pref_d = median_sep_pref_dscores_ses_all.top_50_d_score_prefori_d1_d2;
    top_d1_d2_pref_d_all = [top_d1_d2_pref_d_all top_d1_d2_pref_d];
    top_d2_d3_pref_d = median_sep_pref_dscores_ses_all.top_50_d_score_prefori_d2_d3;
    top_d2_d3_pref_d_all = [top_d2_d3_pref_d_all top_d2_d3_pref_d];
    top_d3_d4_pref_d = median_sep_pref_dscores_ses_all.top_50_d_score_prefori_d3_d4;
    top_d3_d4_pref_d_all = [top_d3_d4_pref_d_all top_d3_d4_pref_d];
    top_d1_d3_pref_d = median_sep_pref_dscores_ses_all.top_50_d_score_prefori_d1_d3;
    top_d1_d3_pref_d_all = [top_d1_d3_pref_d_all top_d1_d3_pref_d];
    top_d1_d4_pref_d = median_sep_pref_dscores_ses_all.top_50_d_score_prefori_d1_d4;
    top_d1_d4_pref_d_all = [top_d1_d4_pref_d_all top_d1_d4_pref_d];

    bot_d1_d2_pref_d = median_sep_pref_dscores_ses_all.bot_50_d_score_prefori_d1_d2;
    bot_d1_d2_pref_d_all = [bot_d1_d2_pref_d_all bot_d1_d2_pref_d];
    bot_d2_d3_pref_d = median_sep_pref_dscores_ses_all.bot_50_d_score_prefori_d2_d3;
    bot_d2_d3_pref_d_all = [bot_d2_d3_pref_d_all bot_d2_d3_pref_d];
    bot_d3_d4_pref_d = median_sep_pref_dscores_ses_all.bot_50_d_score_prefori_d3_d4;
    bot_d3_d4_pref_d_all = [bot_d3_d4_pref_d_all bot_d3_d4_pref_d];
    bot_d1_d3_pref_d = median_sep_pref_dscores_ses_all.bot_50_d_score_prefori_d1_d3;
    bot_d1_d3_pref_d_all = [bot_d1_d3_pref_d_all bot_d1_d3_pref_d];
    bot_d1_d4_pref_d = median_sep_pref_dscores_ses_all.bot_50_d_score_prefori_d1_d4;
    bot_d1_d4_pref_d_all = [bot_d1_d4_pref_d_all bot_d1_d4_pref_d];

    sharp_d1_d2_pref_d = wide_sharp_sep_pref_dscores_ses_all.sharp_d_score_prefori_d1_d2;
    sharp_d1_d2_pref_d_all = [sharp_d1_d2_pref_d_all sharp_d1_d2_pref_d];
    sharp_d2_d3_pref_d = wide_sharp_sep_pref_dscores_ses_all.sharp_d_score_prefori_d2_d3;
    sharp_d2_d3_pref_d_all = [sharp_d2_d3_pref_d_all sharp_d2_d3_pref_d];
    sharp_d3_d4_pref_d = wide_sharp_sep_pref_dscores_ses_all.sharp_d_score_prefori_d3_d4;
    sharp_d3_d4_pref_d_all = [sharp_d3_d4_pref_d_all sharp_d3_d4_pref_d];
    sharp_d1_d3_pref_d = wide_sharp_sep_pref_dscores_ses_all.sharp_d_score_prefori_d1_d3;
    sharp_d1_d3_pref_d_all = [sharp_d1_d3_pref_d_all sharp_d1_d3_pref_d];
    sharp_d1_d4_pref_d = wide_sharp_sep_pref_dscores_ses_all.sharp_d_score_prefori_d1_d4;
    sharp_d1_d4_pref_d_all = [sharp_d1_d4_pref_d_all sharp_d1_d4_pref_d];

    wide_d1_d2_pref_d = wide_sharp_sep_pref_dscores_ses_all.wide_d_score_prefori_d1_d2;
    wide_d1_d2_pref_d_all = [wide_d1_d2_pref_d_all wide_d1_d2_pref_d];
    wide_d2_d3_pref_d = wide_sharp_sep_pref_dscores_ses_all.wide_d_score_prefori_d2_d3;
    wide_d2_d3_pref_d_all = [wide_d2_d3_pref_d_all wide_d2_d3_pref_d];
    wide_d3_d4_pref_d = wide_sharp_sep_pref_dscores_ses_all.wide_d_score_prefori_d3_d4;
    wide_d3_d4_pref_d_all = [wide_d3_d4_pref_d_all wide_d3_d4_pref_d];
    wide_d1_d3_pref_d = wide_sharp_sep_pref_dscores_ses_all.wide_d_score_prefori_d1_d3;
    wide_d1_d3_pref_d_all = [wide_d1_d3_pref_d_all wide_d1_d3_pref_d];
    wide_d1_d4_pref_d = wide_sharp_sep_pref_dscores_ses_all.wide_d_score_prefori_d1_d4;
    wide_d1_d4_pref_d_all = [wide_d1_d4_pref_d_all wide_d1_d4_pref_d];

    d1_k = k_max_ses_all.d1_k;
    d1_k_all = [d1_k_all d1_k];
    d2_k = k_max_ses_all.d2_k;
    d2_k_all = [d2_k_all d2_k];
    d3_k = k_max_ses_all.d3_k;
    d3_k_all = [d3_k_all d3_k];
    d4_k = k_max_ses_all.d4_k;
    d4_k_all = [d4_k_all d4_k];

    d1_max = k_max_ses_all.d1_max;
    d1_max_all = [d1_max_all d1_max];
    d2_max = k_max_ses_all.d2_max;
    d2_max_all = [d2_max_all d2_max];
    d3_max = k_max_ses_all.d3_max;
    d3_max_all = [d3_max_all d3_max];
    d4_max = k_max_ses_all.d4_max;
    d4_max_all = [d4_max_all d4_max];

    d1_reliability = reliability_ses_all.d1_reliability;
    d1_reliability_all = [d1_reliability_all d1_reliability];
    d2_reliability = reliability_ses_all.d2_reliability;
    d2_reliability_all = [d2_reliability_all d2_reliability];
    d3_reliability = reliability_ses_all.d3_reliability;
    d3_reliability_all = [d3_reliability_all d3_reliability];
    d4_reliability = reliability_ses_all.d4_reliability;
    d4_reliability_all = [d4_reliability_all d4_reliability];

  
end

%%
fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(d2_d3_pref_d_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',2);
l=cdfplot(d3_d4_pref_d_all);
set(l, 'LineStyle', '-', 'Color', [.8 .8 .8], 'LineWidth',2);
legend(['Base1 - Base2'], ['Base2 - Post-dark'], ['Post-dark - Post-dark+7d'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['7d_dark_pool_prefori_change.pdf']), '-dpdf', '-bestfit')


%%
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(d1_k_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(d2_k_all);
set(j, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(d3_k_all);
set(l, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
m=cdfplot(d4_k_all);
set(m, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
legend(['Base1'], ['Base2'], ['Post-dark'], ['Post-dark+7d'], 'Location', 'southeast')
xlim([0 30])
xlabel(['k value'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['7d_dark_pool_k_vals.pdf']), '-dpdf', '-bestfit')

%%
fig = figure;
% sgtitle('Max dF/F values')
% a=subplot(2,2,1)
h=cdfplot(d1_max_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(d2_max_all);
set(j, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(d3_max_all);
set(l, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
m=cdfplot(d4_max_all);
set(m, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
legend(['Base1'], ['Base2'], ['Post-dark'], ['Post-dark+7d'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['7d_dark_pool_max_vals.pdf']), '-dpdf', '-bestfit')

%%
fig = figure;
% sgtitle('Reliability dF/F values')
% a=subplot(2,2,1)
h=cdfplot(d1_reliability_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(d2_reliability_all);
set(j, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(d3_reliability_all);
set(l, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
m=cdfplot(d4_reliability_all);
set(m, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
legend(['Base1'], ['Base2'], ['Post-dark'], ['Post-dark+7d'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Fit Reliability'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['7d_dark_pool_reliability_vals.pdf']), '-dpdf', '-bestfit')



%%
fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(d1_d3_pref_d_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',2);
l=cdfplot(d1_d4_pref_d_all);
set(l, 'LineStyle', '-', 'Color', [.8 .8 .8], 'LineWidth',2);
legend(['Base1 - Base2'], ['Base1 - Post-dark'], ['Base1 - Post-dark+7d'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['d1_to_others_7d_dark_pool_prefori_change.pdf']), '-dpdf', '-bestfit')


%%
%median split analysis for 4d dark group
fig = figure;
% sgtitle('Pref Ori Changes')
a=subplot(2,2,1)
h=cdfplot(bot_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(bot_d2_d3_pref_d_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',2);
l=cdfplot(bot_d3_d4_pref_d_all);
set(l, 'LineStyle', '-', 'Color', [.8 .8 .8], 'LineWidth',2);
legend(['Base1 - Base2'], ['Base2 - Post-dark'], ['Post-dark - Post-dark+7d'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Bottom 50% Responsivity'])

% sgtitle('Pref Ori Changes')
b=subplot(2,2,2)
h=cdfplot(bot_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(bot_d1_d3_pref_d_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',2);
l=cdfplot(bot_d1_d4_pref_d_all);
set(l, 'LineStyle', '-', 'Color', [.8 .8 .8], 'LineWidth',2);
legend(['Base1 - Base2'], ['Base1 - Post-dark'], ['Base1 - Post-dark+7d'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Bottom 50% Responsivity'])

c=subplot(2,2,3)
h=cdfplot(top_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(top_d2_d3_pref_d_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',2);
l=cdfplot(top_d3_d4_pref_d_all);
set(l, 'LineStyle', '-', 'Color', [.8 .8 .8], 'LineWidth',2);
legend(['Base1 - Base2'], ['Base2 - Post-dark'], ['Post-dark - Post-dark+7d'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Top 50% Responsivity'])

d=subplot(2,2,4)
h=cdfplot(top_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(top_d1_d3_pref_d_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',2);
l=cdfplot(top_d1_d4_pref_d_all);
set(l, 'LineStyle', '-', 'Color', [.8 .8 .8], 'LineWidth',2);
legend(['Base1 - Base2'], ['Base1 - Post-dark'], ['Base1 - Post-dark+7d'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Top 50% Responsivity'])
hold off

print(fullfile(realfnout, ['median_split_7d_dark_pool_prefori_change.pdf']), '-dpdf', '-bestfit')


