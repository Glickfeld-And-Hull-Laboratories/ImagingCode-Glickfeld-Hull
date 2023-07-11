%CLEAR EVERYTHINg
clear all
clear all global
clc
close all



%%
%LOAD DATA AND IDENTIFY FOLDERS TO SAVE TO
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_virusTG'; %folder to save files to
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_pooled_RW'; %folder to save files to
dataset = 'exp_list_tjw'; %experiment list to pick files from
eval(dataset); %load dataset


%%
%% ARC PROMOTER MICE
%%
preRW_d1_pref_all = [];
preRW_d2_pref_all = [];
preRW_d3_pref_all = [];
preRW_d1_d2_pref_d_all = [];
preRW_d1_d3_pref_d_all = [];
preRW_d1_k_all = [];
preRW_d2_k_all = [];
preRW_d3_k_all = [];
preRW_d1_max_all = [];
preRW_d2_max_all = [];
preRW_d3_max_all = [];

preRW_top_d1_d2_pref_d_all = [];
preRW_top_d1_d3_pref_d_all = [];

preRW_bot_d1_d2_pref_d_all = [];
preRW_bot_d1_d3_pref_d_all = [];

postRW_d1_pref_all = [];
postRW_d2_pref_all = [];
postRW_d3_pref_all = [];
postRW_d1_d2_pref_d_all = [];
postRW_d1_d3_pref_d_all = [];
postRW_d1_k_all = [];
postRW_d2_k_all = [];
postRW_d3_k_all = [];
postRW_d1_max_all = [];
postRW_d2_max_all = [];
postRW_d3_max_all = [];

postRW_top_d1_d2_pref_d_all = [];
postRW_top_d1_d3_pref_d_all = [];

postRW_bot_d1_d2_pref_d_all = [];
postRW_bot_d1_d3_pref_d_all = [];

preRW_prefori_ses_all = [];
preRW_prefori_scores = [];
preRW_pref_dscores_ses_all = [];
preRW_pref_dscores = [];
preRW_k_max_ses_all = [];
preRW_k_max_scores = [];
preRW_median_sep_pref_dscores_ses_all = [];
preRW_median_pref_dscores = [];

postRW_prefori_ses_all = [];
postRW_prefori_scores = [];
postRW_pref_dscores_ses_all = [];
postRW_pref_dscores = [];
postRW_k_max_ses_all = [];
postRW_k_max_scores = [];
postRW_median_sep_pref_dscores_ses_all = [];
postRW_median_pref_dscores = [];



%%

preRW_list = [19 22 40 43];
for iexp = preRW_list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};

    preRW_prefori_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes']));
    
    preRW_prefori_ses2 = [];
    
    preRW_prefori_ses_all = [preRW_prefori_ses1 preRW_prefori_ses2];
    preRW_prefori_scores = [preRW_prefori_scores preRW_prefori_ses_all];

    preRW_pref_dscores_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores']));
   
    preRW_pref_dscores_ses_2 = [];
    
    preRW_pref_dscores_ses_all = [preRW_pref_dscores_ses_1 preRW_pref_dscores_ses_2];
    preRW_pref_dscores = [preRW_pref_dscores preRW_pref_dscores_ses_all];

    preRW_k_max_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes']));
    
    preRW_k_max_ses_2 = [];
    
    preRW_k_max_ses_all = [preRW_k_max_ses_1 preRW_k_max_ses_2];
    preRW_k_max_scores = [preRW_k_max_scores preRW_k_max_ses_all];

    preRW_median_sep_pref_dscores_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'median_sep_d_scores']));
   
    preRW_median_sep_pref_dscores_ses_2 = [];
    
    preRW_median_sep_pref_dscores_ses_all = [preRW_median_sep_pref_dscores_ses_1 preRW_median_sep_pref_dscores_ses_2];
    preRW_median_pref_dscores = [preRW_median_pref_dscores preRW_median_sep_pref_dscores_ses_all];
end

%%
for idata = 1:length(preRW_prefori_scores)

    preRW_d1_d2_pref_d = preRW_pref_dscores(idata).green_d_score_prefori_d1_d2;
    preRW_d1_d2_pref_d_all = [preRW_d1_d2_pref_d_all preRW_d1_d2_pref_d];
    preRW_d1_d3_pref_d = preRW_pref_dscores(idata).green_d_score_prefori_d1_d3;
    preRW_d1_d3_pref_d_all = [preRW_d1_d3_pref_d_all preRW_d1_d3_pref_d];

    preRW_top_d1_d2_pref_d = preRW_median_pref_dscores(idata).green_top_50_d_score_prefori_d1_d2;
    preRW_top_d1_d2_pref_d_all = [preRW_top_d1_d2_pref_d_all preRW_top_d1_d2_pref_d];
    preRW_top_d1_d3_pref_d = preRW_median_pref_dscores(idata).green_top_50_d_score_prefori_d1_d3;
    preRW_top_d1_d3_pref_d_all = [preRW_top_d1_d3_pref_d_all preRW_top_d1_d3_pref_d];

    preRW_bot_d1_d2_pref_d = preRW_median_pref_dscores(idata).green_bot_50_d_score_prefori_d1_d2;
    preRW_bot_d1_d2_pref_d_all = [preRW_bot_d1_d2_pref_d_all preRW_bot_d1_d2_pref_d];
    preRW_bot_d1_d3_pref_d = preRW_median_pref_dscores(idata).green_bot_50_d_score_prefori_d1_d3;
    preRW_bot_d1_d3_pref_d_all = [preRW_bot_d1_d3_pref_d_all preRW_bot_d1_d3_pref_d];


    preRW_d1_k = preRW_k_max_scores(idata).green_d1_k;
    preRW_d1_k_all = [preRW_d1_k_all preRW_d1_k];
    preRW_d2_k = preRW_k_max_scores(idata).green_d2_k;
    preRW_d2_k_all = [preRW_d2_k_all preRW_d2_k];
    preRW_d3_k = preRW_k_max_scores(idata).green_d3_k;
    preRW_d3_k_all = [preRW_d3_k_all preRW_d3_k];

    preRW_d1_max = preRW_k_max_scores(idata).green_d1_max;
    preRW_d1_max_all = [preRW_d1_max_all preRW_d1_max];
    preRW_d2_max = preRW_k_max_scores(idata).green_d2_max;
    preRW_d2_max_all = [preRW_d2_max_all preRW_d2_max];
    preRW_d3_max = preRW_k_max_scores(idata).green_d3_max;
    preRW_d3_max_all = [preRW_d3_max_all preRW_d3_max];

  
end


%%
postRW_list = [34 37 46 49];
for iexp = postRW_list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};

    postRW_prefori_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes_2']));
   
    postRW_prefori_ses2 = [];
    
    postRW_prefori_ses_all = [postRW_prefori_ses1 postRW_prefori_ses2];
    postRW_prefori_scores = [postRW_prefori_scores postRW_prefori_ses_all];

    postRW_pref_dscores_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores_2']));
    
    postRW_pref_dscores_ses_2 = [];
    
    postRW_pref_dscores_ses_all = [postRW_pref_dscores_ses_1 postRW_pref_dscores_ses_2];
    postRW_pref_dscores = [postRW_pref_dscores postRW_pref_dscores_ses_all];

    postRW_k_max_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes_2']));
    
    postRW_k_max_ses_2 = [];
    
    postRW_k_max_ses_all = [postRW_k_max_ses_1 postRW_k_max_ses_2];
    postRW_k_max_scores = [postRW_k_max_scores postRW_k_max_ses_all];

    postRW_median_sep_pref_dscores_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'median_sep_d_scores_2']));
    
    postRW_median_sep_pref_dscores_ses_2 = [];
    
    postRW_median_sep_pref_dscores_ses_all = [postRW_median_sep_pref_dscores_ses_1 postRW_median_sep_pref_dscores_ses_2];
    postRW_median_pref_dscores = [postRW_median_pref_dscores postRW_median_sep_pref_dscores_ses_all];
end

%%
for idata = 1:length(postRW_prefori_scores)

    postRW_d1_d2_pref_d = postRW_pref_dscores(idata).green_d_score_prefori_d1_d2;
    postRW_d1_d2_pref_d_all = [postRW_d1_d2_pref_d_all postRW_d1_d2_pref_d];
    postRW_d1_d3_pref_d = postRW_pref_dscores(idata).green_d_score_prefori_d1_d3;
    postRW_d1_d3_pref_d_all = [postRW_d1_d3_pref_d_all postRW_d1_d3_pref_d];

    postRW_top_d1_d2_pref_d = postRW_median_pref_dscores(idata).green_top_50_d_score_prefori_d1_d2;
    postRW_top_d1_d2_pref_d_all = [postRW_top_d1_d2_pref_d_all postRW_top_d1_d2_pref_d];
    postRW_top_d1_d3_pref_d = postRW_median_pref_dscores(idata).green_top_50_d_score_prefori_d1_d3;
    postRW_top_d1_d3_pref_d_all = [postRW_top_d1_d3_pref_d_all postRW_top_d1_d3_pref_d];

    postRW_bot_d1_d2_pref_d = postRW_median_pref_dscores(idata).green_bot_50_d_score_prefori_d1_d2;
    postRW_bot_d1_d2_pref_d_all = [postRW_bot_d1_d2_pref_d_all postRW_bot_d1_d2_pref_d];
    postRW_bot_d1_d3_pref_d = postRW_median_pref_dscores(idata).green_bot_50_d_score_prefori_d1_d3;
    postRW_bot_d1_d3_pref_d_all = [postRW_bot_d1_d3_pref_d_all postRW_bot_d1_d3_pref_d];


    postRW_d1_k = postRW_k_max_scores(idata).green_d1_k;
    postRW_d1_k_all = [postRW_d1_k_all postRW_d1_k];
    postRW_d2_k = postRW_k_max_scores(idata).green_d2_k;
    postRW_d2_k_all = [postRW_d2_k_all postRW_d2_k];
    postRW_d3_k = postRW_k_max_scores(idata).green_d3_k;
    postRW_d3_k_all = [postRW_d3_k_all postRW_d3_k];

    postRW_d1_max = postRW_k_max_scores(idata).green_d1_max;
    postRW_d1_max_all = [postRW_d1_max_all postRW_d1_max];
    postRW_d2_max = postRW_k_max_scores(idata).green_d2_max;
    postRW_d2_max_all = [postRW_d2_max_all postRW_d2_max];
    postRW_d3_max = postRW_k_max_scores(idata).green_d3_max;
    postRW_d3_max_all = [postRW_d3_max_all postRW_d3_max];

  
end



%%
fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(preRW_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(preRW_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(postRW_d1_d2_pref_d_all);
set(l, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
m=cdfplot(postRW_d1_d3_pref_d_all);
set(m, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 - 2 (preRW)'], ['Ses1 - 3 (preRW)'], ['Ses1 - 2 (postRW)'], ['Ses1 - 3 (postRW)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['preVpostRW_changeprefcdf.pdf']), '-dpdf', '-bestfit')

%%
preRW_pref_d_all = [preRW_d1_d2_pref_d_all preRW_d1_d3_pref_d_all];
postRW_pref_d_all = [postRW_d1_d2_pref_d_all postRW_d1_d3_pref_d_all];


fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(preRW_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(postRW_pref_d_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['All Days (preRW)'], ['All days (postRW)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['preVpostRW_ALLDAYS_changeprefcdf.pdf']), '-dpdf', '-bestfit')


%%
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(preRW_d1_k_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(preRW_d2_k_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(preRW_d3_k_all);
set(l, 'LineStyle', ':', 'Color', 'g', 'LineWidth',2);
m=cdfplot(postRW_d1_k_all);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
n=cdfplot(postRW_d2_k_all);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
o=cdfplot(postRW_d3_k_all);
set(o, 'LineStyle', ':', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 (preRW)'], ['Ses2 (preRW)'], ['Ses3 (preRW)'], ['Ses1 (postRW)'], ['Ses2 (postRW)'], ['Ses3 (postRW)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['k value'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['preVpostRW_kbysession.pdf']), '-dpdf', '-bestfit')

%%
preRW_k_all = [preRW_d1_k_all preRW_d2_k_all preRW_d3_k_all];
postRW_k_all = [postRW_d1_k_all postRW_d2_k_all postRW_d3_k_all];

fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(preRW_k_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(postRW_k_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['All Sessions (preRW)'], ['All Sessions (postRW)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['k value'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['preVpostRW_kallsessions.pdf']), '-dpdf', '-bestfit')

%%
%day 1 k
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(preRW_d1_k_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(postRW_d1_k_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['Day 1 (preRW)'], ['Day 1 (postRW)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['k value'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['preVpostRW_d1_k.pdf']), '-dpdf', '-bestfit')


%%
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(preRW_d1_max_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(preRW_d2_max_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(preRW_d3_max_all);
set(l, 'LineStyle', ':', 'Color', 'g', 'LineWidth',2);
m=cdfplot(postRW_d1_max_all);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
n=cdfplot(postRW_d2_max_all);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
o=cdfplot(postRW_d3_max_all);
set(o, 'LineStyle', ':', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 (preRW)'], ['Ses2 (preRW)'], ['Ses3 (preRW)'], ['Ses1 (postRW)'], ['Ses2 (postRW)'], ['Ses3 (postRW)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F value'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['preVpostRW_maxbysession.pdf']), '-dpdf', '-bestfit')

%%
preRW_max_all = [preRW_d1_max_all preRW_d2_max_all preRW_d3_max_all];
postRW_max_all = [postRW_d1_max_all postRW_d2_max_all postRW_d3_max_all];

fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(preRW_max_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(postRW_max_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['All Sessions (preRW)'], ['All Sessions (postRW)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F value'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['preVpostRW_maxallsessions.pdf']), '-dpdf', '-bestfit')

%%
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(preRW_d1_max_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(postRW_d1_max_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['Day 1 (preRW)'], ['Day 1 (postRW)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F value'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['preVpostRW_d1_max.pdf']), '-dpdf', '-bestfit')


%%
%median split analysis for 4d dark group
fig = figure;
% sgtitle('Pref Ori Changes')
a=subplot(2,2,1)
h=cdfplot(preRW_bot_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(preRW_bot_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
legend(['Ses1 - 2 (preRW)'], ['Ses1 - 3 (preRW)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Bottom 50% Responsivity'])

% sgtitle('Pref Ori Changes')
b=subplot(2,2,2)
h=cdfplot(postRW_bot_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(postRW_bot_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 - 2 (postRW)'], ['Ses1 - 3 (postRW)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Bottom 50% Responsivity'])

% sgtitle('Pref Ori Changes')
a=subplot(2,2,3)
h=cdfplot(preRW_top_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(preRW_top_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
legend(['Ses1 - 2 (preRW)'], ['Ses1 - 3 (preRW)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Top 50% Responsivity'])

% sgtitle('Pref Ori Changes')
b=subplot(2,2,4)
h=cdfplot(postRW_top_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(postRW_top_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 - 2 (postRW)'], ['Ses1 - 3 (postRW)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Top 50% Responsivity'])



print(fullfile(realfnout, ['preVpostRW_medsplit_changeprefcdf.pdf']), '-dpdf', '-bestfit')

%%
preRW_bot_pref_d_all = [preRW_bot_d1_d2_pref_d_all preRW_bot_d1_d3_pref_d_all];
preRW_top_pref_d_all = [preRW_top_d1_d2_pref_d_all preRW_top_d1_d3_pref_d_all];
postRW_bot_pref_d_all = [postRW_bot_d1_d2_pref_d_all postRW_bot_d1_d3_pref_d_all];
postRW_top_pref_d_all = [postRW_top_d1_d2_pref_d_all postRW_top_d1_d3_pref_d_all];

fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(preRW_bot_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(preRW_top_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(postRW_bot_pref_d_all);
set(l, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
m=cdfplot(postRW_top_pref_d_all);
set(m, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Bottom 50% (preRW)'], ['Top 50% (preRW)'], ['Bottom 50% (postRW)'], ['Top 50% (postRW)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['preVpostRW_medsplit_changeprefcdf_ALLDAYS.pdf']), '-dpdf', '-bestfit')




%%
%k and max change scores
postRW_d1_d2_k_all = abs(postRW_d1_k_all-postRW_d2_k_all);
postRW_d1_d3_k_all = abs(postRW_d1_k_all-postRW_d3_k_all);
postRW_d_k_all = [postRW_d1_d2_k_all postRW_d1_d3_k_all];

preRW_d1_d2_k_all = abs(preRW_d1_k_all-preRW_d2_k_all);
preRW_d1_d3_k_all = abs(preRW_d1_k_all-preRW_d3_k_all);
preRW_d_k_all = [preRW_d1_d2_k_all preRW_d1_d3_k_all];

postRW_d1_d2_max_all = abs(postRW_d1_max_all-postRW_d2_max_all);
postRW_d1_d3_max_all = abs(postRW_d1_max_all-postRW_d3_max_all);
postRW_d_max_all = [postRW_d1_d2_max_all postRW_d1_d3_max_all];

preRW_d1_d2_max_all = abs(preRW_d1_max_all-preRW_d2_max_all);
preRW_d1_d3_max_all = abs(preRW_d1_max_all-preRW_d3_max_all);
preRW_d_max_all = [preRW_d1_d2_max_all preRW_d1_d3_max_all];

%%
fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(preRW_d_k_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(postRW_d_k_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['preRW'], ['postRW'], 'Location', 'southeast')
xlim([0 30])
xlabel(['Change in k Value'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['preVpostRW_changek.pdf']), '-dpdf', '-bestfit')

%%

fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(preRW_d_max_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(postRW_d_max_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['preRW'], ['postRW'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Change in Max dF/F Value'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['preVpostRW_changemax.pdf']), '-dpdf', '-bestfit')



%%
%%
%%


