%CLEAR EVERYTHINg
clear all
clear all global
clc
close all



%%
%LOAD DATA AND IDENTIFY FOLDERS TO SAVE TO
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_virusTG'; %folder to save files to
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_pooled_virusTG'; %folder to save files to
dataset = 'exp_list_tjw'; %experiment list to pick files from
eval(dataset); %load dataset


%%
%% ARC PROMOTER MICE
%%
virus_d1_pref_all = [];
virus_d2_pref_all = [];
virus_d3_pref_all = [];
virus_d1_d2_pref_d_all = [];
virus_d1_d3_pref_d_all = [];
virus_d1_k_all = [];
virus_d2_k_all = [];
virus_d3_k_all = [];
virus_d1_max_all = [];
virus_d2_max_all = [];
virus_d3_max_all = [];

virus_top_d1_d2_pref_d_all = [];
virus_top_d1_d3_pref_d_all = [];

virus_bot_d1_d2_pref_d_all = [];
virus_bot_d1_d3_pref_d_all = [];

TG_d1_pref_all = [];
TG_d2_pref_all = [];
TG_d3_pref_all = [];
TG_d1_d2_pref_d_all = [];
TG_d1_d3_pref_d_all = [];
TG_d1_k_all = [];
TG_d2_k_all = [];
TG_d3_k_all = [];
TG_d1_max_all = [];
TG_d2_max_all = [];
TG_d3_max_all = [];

TG_top_d1_d2_pref_d_all = [];
TG_top_d1_d3_pref_d_all = [];

TG_bot_d1_d2_pref_d_all = [];
TG_bot_d1_d3_pref_d_all = [];

virus_prefori_ses_all = [];
virus_prefori_scores = [];
virus_pref_dscores_ses_all = [];
virus_pref_dscores = [];
virus_k_max_ses_all = [];
virus_k_max_scores = [];
virus_median_sep_pref_dscores_ses_all = [];
virus_median_pref_dscores = [];

TG_prefori_ses_all = [];
TG_prefori_scores = [];
TG_pref_dscores_ses_all = [];
TG_pref_dscores = [];
TG_k_max_ses_all = [];
TG_k_max_scores = [];
TG_median_sep_pref_dscores_ses_all = [];
TG_median_pref_dscores = [];



%%

virus_list = [1 7 13 16 19 22 25 34 37];
for iexp = virus_list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};

    virus_prefori_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes_2.mat']));
        virus_prefori_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes_2']));
    else
        virus_prefori_ses2 = [];
    end
    virus_prefori_ses_all = [virus_prefori_ses1 virus_prefori_ses2];
    virus_prefori_scores = [virus_prefori_scores virus_prefori_ses_all];

    virus_pref_dscores_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores_2.mat']));
        virus_pref_dscores_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores_2']));
    else
        virus_pref_dscores_ses_2 = [];
    end
    virus_pref_dscores_ses_all = [virus_pref_dscores_ses_1 virus_pref_dscores_ses_2];
    virus_pref_dscores = [virus_pref_dscores virus_pref_dscores_ses_all];

    virus_k_max_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes_2.mat']));
        virus_k_max_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes_2']));
    else
        virus_k_max_ses_2 = [];
    end
    virus_k_max_ses_all = [virus_k_max_ses_1 virus_k_max_ses_2];
    virus_k_max_scores = [virus_k_max_scores virus_k_max_ses_all];

    virus_median_sep_pref_dscores_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'median_sep_d_scores']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'median_sep_d_scores_2.mat']));
        virus_median_sep_pref_dscores_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'median_sep_d_scores_2']));
    else
        virus_median_sep_pref_dscores_ses_2 = [];
    end
    virus_median_sep_pref_dscores_ses_all = [virus_median_sep_pref_dscores_ses_1 virus_median_sep_pref_dscores_ses_2];
    virus_median_pref_dscores = [virus_median_pref_dscores virus_median_sep_pref_dscores_ses_all];
end

%%
for idata = 1:length(virus_prefori_scores)

    virus_d1_d2_pref_d = virus_pref_dscores(idata).green_d_score_prefori_d1_d2;
    virus_d1_d2_pref_d_all = [virus_d1_d2_pref_d_all virus_d1_d2_pref_d];
    virus_d1_d3_pref_d = virus_pref_dscores(idata).green_d_score_prefori_d1_d3;
    virus_d1_d3_pref_d_all = [virus_d1_d3_pref_d_all virus_d1_d3_pref_d];

    virus_top_d1_d2_pref_d = virus_median_pref_dscores(idata).green_top_50_d_score_prefori_d1_d2;
    virus_top_d1_d2_pref_d_all = [virus_top_d1_d2_pref_d_all virus_top_d1_d2_pref_d];
    virus_top_d1_d3_pref_d = virus_median_pref_dscores(idata).green_top_50_d_score_prefori_d1_d3;
    virus_top_d1_d3_pref_d_all = [virus_top_d1_d3_pref_d_all virus_top_d1_d3_pref_d];

    virus_bot_d1_d2_pref_d = virus_median_pref_dscores(idata).green_bot_50_d_score_prefori_d1_d2;
    virus_bot_d1_d2_pref_d_all = [virus_bot_d1_d2_pref_d_all virus_bot_d1_d2_pref_d];
    virus_bot_d1_d3_pref_d = virus_median_pref_dscores(idata).green_bot_50_d_score_prefori_d1_d3;
    virus_bot_d1_d3_pref_d_all = [virus_bot_d1_d3_pref_d_all virus_bot_d1_d3_pref_d];


    virus_d1_k = virus_k_max_scores(idata).green_d1_k;
    virus_d1_k_all = [virus_d1_k_all virus_d1_k];
    virus_d2_k = virus_k_max_scores(idata).green_d2_k;
    virus_d2_k_all = [virus_d2_k_all virus_d2_k];
    virus_d3_k = virus_k_max_scores(idata).green_d3_k;
    virus_d3_k_all = [virus_d3_k_all virus_d3_k];

    virus_d1_max = virus_k_max_scores(idata).green_d1_max;
    virus_d1_max_all = [virus_d1_max_all virus_d1_max];
    virus_d2_max = virus_k_max_scores(idata).green_d2_max;
    virus_d2_max_all = [virus_d2_max_all virus_d2_max];
    virus_d3_max = virus_k_max_scores(idata).green_d3_max;
    virus_d3_max_all = [virus_d3_max_all virus_d3_max];

  
end


%%
TG_list = [28 31 40 43 46 49];
for iexp = TG_list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};

    TG_prefori_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes_2.mat']));
        TG_prefori_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes_2']));
    else
        TG_prefori_ses2 = [];
    end
    TG_prefori_ses_all = [TG_prefori_ses1 TG_prefori_ses2];
    TG_prefori_scores = [TG_prefori_scores TG_prefori_ses_all];

    TG_pref_dscores_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores_2.mat']));
        TG_pref_dscores_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores_2']));
    else
        TG_pref_dscores_ses_2 = [];
    end
    TG_pref_dscores_ses_all = [TG_pref_dscores_ses_1 TG_pref_dscores_ses_2];
    TG_pref_dscores = [TG_pref_dscores TG_pref_dscores_ses_all];

    TG_k_max_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes_2.mat']));
        TG_k_max_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes_2']));
    else
        TG_k_max_ses_2 = [];
    end
    TG_k_max_ses_all = [TG_k_max_ses_1 TG_k_max_ses_2];
    TG_k_max_scores = [TG_k_max_scores TG_k_max_ses_all];

    TG_median_sep_pref_dscores_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'median_sep_d_scores']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'median_sep_d_scores_2.mat']));
        TG_median_sep_pref_dscores_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'median_sep_d_scores_2']));
    else
        TG_median_sep_pref_dscores_ses_2 = [];
    end
    TG_median_sep_pref_dscores_ses_all = [TG_median_sep_pref_dscores_ses_1 TG_median_sep_pref_dscores_ses_2];
    TG_median_pref_dscores = [TG_median_pref_dscores TG_median_sep_pref_dscores_ses_all];
end

%%
for idata = 1:length(TG_prefori_scores)

    TG_d1_d2_pref_d = TG_pref_dscores(idata).green_d_score_prefori_d1_d2;
    TG_d1_d2_pref_d_all = [TG_d1_d2_pref_d_all TG_d1_d2_pref_d];
    TG_d1_d3_pref_d = TG_pref_dscores(idata).green_d_score_prefori_d1_d3;
    TG_d1_d3_pref_d_all = [TG_d1_d3_pref_d_all TG_d1_d3_pref_d];

    TG_top_d1_d2_pref_d = TG_median_pref_dscores(idata).green_top_50_d_score_prefori_d1_d2;
    TG_top_d1_d2_pref_d_all = [TG_top_d1_d2_pref_d_all TG_top_d1_d2_pref_d];
    TG_top_d1_d3_pref_d = TG_median_pref_dscores(idata).green_top_50_d_score_prefori_d1_d3;
    TG_top_d1_d3_pref_d_all = [TG_top_d1_d3_pref_d_all TG_top_d1_d3_pref_d];

    TG_bot_d1_d2_pref_d = TG_median_pref_dscores(idata).green_bot_50_d_score_prefori_d1_d2;
    TG_bot_d1_d2_pref_d_all = [TG_bot_d1_d2_pref_d_all TG_bot_d1_d2_pref_d];
    TG_bot_d1_d3_pref_d = TG_median_pref_dscores(idata).green_bot_50_d_score_prefori_d1_d3;
    TG_bot_d1_d3_pref_d_all = [TG_bot_d1_d3_pref_d_all TG_bot_d1_d3_pref_d];


    TG_d1_k = TG_k_max_scores(idata).green_d1_k;
    TG_d1_k_all = [TG_d1_k_all TG_d1_k];
    TG_d2_k = TG_k_max_scores(idata).green_d2_k;
    TG_d2_k_all = [TG_d2_k_all TG_d2_k];
    TG_d3_k = TG_k_max_scores(idata).green_d3_k;
    TG_d3_k_all = [TG_d3_k_all TG_d3_k];

    TG_d1_max = TG_k_max_scores(idata).green_d1_max;
    TG_d1_max_all = [TG_d1_max_all TG_d1_max];
    TG_d2_max = TG_k_max_scores(idata).green_d2_max;
    TG_d2_max_all = [TG_d2_max_all TG_d2_max];
    TG_d3_max = TG_k_max_scores(idata).green_d3_max;
    TG_d3_max_all = [TG_d3_max_all TG_d3_max];

  
end



%%
fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(virus_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(virus_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(TG_d1_d2_pref_d_all);
set(l, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
m=cdfplot(TG_d1_d3_pref_d_all);
set(m, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 - 2 (Virus)'], ['Ses1 - 3 (Virus)'], ['Ses1 - 2 (TG)'], ['Ses1 - 3 (TG)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['TGvirus_changeprefcdf.pdf']), '-dpdf', '-bestfit')

%%
virus_pref_d_all = [virus_d1_d2_pref_d_all virus_d1_d3_pref_d_all];
TG_pref_d_all = [TG_d1_d2_pref_d_all TG_d1_d3_pref_d_all];


fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(virus_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(TG_pref_d_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['All Days (Virus)'], ['All days (TG)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['TGvirus_ALLDAYS_changeprefcdf.pdf']), '-dpdf', '-bestfit')


%%
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(virus_d1_k_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(virus_d2_k_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(virus_d3_k_all);
set(l, 'LineStyle', ':', 'Color', 'g', 'LineWidth',2);
m=cdfplot(TG_d1_k_all);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
n=cdfplot(TG_d2_k_all);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
o=cdfplot(TG_d3_k_all);
set(o, 'LineStyle', ':', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 (Virus)'], ['Ses2 (Virus)'], ['Ses3 (Virus)'], ['Ses1 (TG)'], ['Ses2 (TG)'], ['Ses3 (TG)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['k value'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['TGvirus_kbysession.pdf']), '-dpdf', '-bestfit')

%%
virus_k_all = [virus_d1_k_all virus_d2_k_all virus_d3_k_all];
TG_k_all = [TG_d1_k_all TG_d2_k_all TG_d3_k_all];

fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(virus_k_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(TG_k_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['All Sessions (Virus)'], ['All Sessions (TG)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['k value'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['TGvirus_kallsessions.pdf']), '-dpdf', '-bestfit')

%%
%day 1 k
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(virus_d1_k_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(TG_d1_k_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['Day 1 (Virus)'], ['Day 1 (TG)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['k value'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['TGvirus_d1_k.pdf']), '-dpdf', '-bestfit')


%%
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(virus_d1_max_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(virus_d2_max_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(virus_d3_max_all);
set(l, 'LineStyle', ':', 'Color', 'g', 'LineWidth',2);
m=cdfplot(TG_d1_max_all);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
n=cdfplot(TG_d2_max_all);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
o=cdfplot(TG_d3_max_all);
set(o, 'LineStyle', ':', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 (Virus)'], ['Ses2 (Virus)'], ['Ses3 (Virus)'], ['Ses1 (TG)'], ['Ses2 (TG)'], ['Ses3 (TG)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F value'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['TGvirus_maxbysession.pdf']), '-dpdf', '-bestfit')

%%
virus_max_all = [virus_d1_max_all virus_d2_max_all virus_d3_max_all];
TG_max_all = [TG_d1_max_all TG_d2_max_all TG_d3_max_all];

fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(virus_max_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(TG_max_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['All Sessions (Virus)'], ['All Sessions (TG)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F value'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['TGvirus_maxallsessions.pdf']), '-dpdf', '-bestfit')

%%
fig = figure;
% sgtitle('k values')
% a=subplot(2,2,1)
h=cdfplot(virus_d1_max_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(TG_d1_max_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['Day 1 (Virus)'], ['Day 1 (TG)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F value'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['TGvirus_d1_max.pdf']), '-dpdf', '-bestfit')


%%
%median split analysis for 4d dark group
fig = figure;
% sgtitle('Pref Ori Changes')
a=subplot(2,2,1)
h=cdfplot(virus_bot_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(virus_bot_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
legend(['Ses1 - 2 (Virus)'], ['Ses1 - 3 (Virus)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Bottom 50% Responsivity'])

% sgtitle('Pref Ori Changes')
b=subplot(2,2,2)
h=cdfplot(TG_bot_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(TG_bot_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 - 2 (TG)'], ['Ses1 - 3 (TG)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Bottom 50% Responsivity'])

% sgtitle('Pref Ori Changes')
a=subplot(2,2,3)
h=cdfplot(virus_top_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(virus_top_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
legend(['Ses1 - 2 (Virus)'], ['Ses1 - 3 (Virus)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Top 50% Responsivity'])

% sgtitle('Pref Ori Changes')
b=subplot(2,2,4)
h=cdfplot(TG_top_d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(TG_top_d1_d3_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Ses1 - 2 (TG)'], ['Ses1 - 3 (TG)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title(['Top 50% Responsivity'])



print(fullfile(realfnout, ['TGvirus_medsplit_changeprefcdf.pdf']), '-dpdf', '-bestfit')

%%
virus_bot_pref_d_all = [virus_bot_d1_d2_pref_d_all virus_bot_d1_d3_pref_d_all];
virus_top_pref_d_all = [virus_top_d1_d2_pref_d_all virus_top_d1_d3_pref_d_all];
TG_bot_pref_d_all = [TG_bot_d1_d2_pref_d_all TG_bot_d1_d3_pref_d_all];
TG_top_pref_d_all = [TG_top_d1_d2_pref_d_all TG_top_d1_d3_pref_d_all];

fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(virus_bot_pref_d_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(virus_top_pref_d_all);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
l=cdfplot(TG_bot_pref_d_all);
set(l, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
m=cdfplot(TG_top_pref_d_all);
set(m, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Bottom 50% (Virus)'], ['Top 50% (Virus)'], ['Bottom 50% (TG)'], ['Top 50% (TG)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['TGvirus_medsplit_changeprefcdf_ALLDAYS.pdf']), '-dpdf', '-bestfit')




%%
%k and max change scores
TG_d1_d2_k_all = abs(TG_d1_k_all-TG_d2_k_all);
TG_d1_d3_k_all = abs(TG_d1_k_all-TG_d3_k_all);
TG_d_k_all = [TG_d1_d2_k_all TG_d1_d3_k_all];

virus_d1_d2_k_all = abs(virus_d1_k_all-virus_d2_k_all);
virus_d1_d3_k_all = abs(virus_d1_k_all-virus_d3_k_all);
virus_d_k_all = [virus_d1_d2_k_all virus_d1_d3_k_all];

TG_d1_d2_max_all = abs(TG_d1_max_all-TG_d2_max_all);
TG_d1_d3_max_all = abs(TG_d1_max_all-TG_d3_max_all);
TG_d_max_all = [TG_d1_d2_max_all TG_d1_d3_max_all];

virus_d1_d2_max_all = abs(virus_d1_max_all-virus_d2_max_all);
virus_d1_d3_max_all = abs(virus_d1_max_all-virus_d3_max_all);
virus_d_max_all = [virus_d1_d2_max_all virus_d1_d3_max_all];

%%
fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(virus_d_k_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(TG_d_k_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['Virus'], ['TG'], 'Location', 'southeast')
xlim([0 30])
xlabel(['Change in k Value'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['TGvirus_changek.pdf']), '-dpdf', '-bestfit')

%%

fig = figure;
% sgtitle('Pref Ori Changes')
% a=subplot(2,2,1)
h=cdfplot(virus_d_max_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold on
j=cdfplot(TG_d_max_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
legend(['Virus'], ['TG'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Change in Max dF/F Value'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(realfnout, ['TGvirus_changemax.pdf']), '-dpdf', '-bestfit')



%%
%%
%%


