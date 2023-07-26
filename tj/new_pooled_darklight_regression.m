%CLEAR EVERYTHINg
clear all
clear all global
clc
close all


%*** SAVE ALL MODELS AND GRAPHS DESIRED***%


%%
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\new_multi_day'; %folder to load files 
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\pooled\regress_four_session'; %folder to save files 
dataset = 'exp_list_darklight_actual_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%
%%
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


list = [1 1+7 1+14 22 22+4 22+8 22+12];
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
%% SAME AS ABOVE BUT LM
%%

d1_pref_LM_all = [];
d2_pref_LM_all = [];
d3_pref_LM_all = [];
d4_pref_LM_all = [];
d1_d2_pref_d_LM_all = [];
d2_d3_pref_d_LM_all = [];
d3_d4_pref_d_LM_all = [];
d1_d3_pref_d_LM_all = [];
d1_d4_pref_d_LM_all = [];
d1_k_LM_all = [];
d2_k_LM_all = [];
d3_k_LM_all = [];
d4_k_LM_all = [];
d1_max_LM_all = [];
d2_max_LM_all = [];
d3_max_LM_all = [];
d4_max_LM_all = [];
d1_reliability_LM_all = [];
d2_reliability_LM_all = [];
d3_reliability_LM_all = [];
d4_reliability_LM_all = [];


list = [38 38+4 38+8];
for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    prefori_ses_LM_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes']));
    pref_dscores_ses_LM_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores']));
    k_max_ses_LM_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes']));
    median_sep_pref_dscores_ses_LM_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'median_sep_d_scores']));
    wide_sharp_sep_pref_dscores_ses_LM_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'wide_sharp_sep_d_scores']));
    reliability_ses_LM_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'fit_reliability']));

    d1_pref = prefori_ses_LM_all.d1_prefori;
    d1_pref_LM_all = [d1_pref_LM_all d1_pref];
    d2_pref = prefori_ses_LM_all.d2_prefori;
    d2_pref_LM_all = [d2_pref_LM_all d2_pref];
    d3_pref = prefori_ses_LM_all.d3_prefori;
    d3_pref_LM_all = [d3_pref_LM_all d3_pref];
    d4_pref = prefori_ses_LM_all.d4_prefori;
    d4_pref_LM_all = [d4_pref_LM_all d4_pref];

    d1_d2_pref_d = pref_dscores_ses_LM_all.d_score_prefori_d1_d2;
    d1_d2_pref_d_LM_all = [d1_d2_pref_d_LM_all d1_d2_pref_d];
    d2_d3_pref_d = pref_dscores_ses_LM_all.d_score_prefori_d2_d3;
    d2_d3_pref_d_LM_all = [d2_d3_pref_d_LM_all d2_d3_pref_d];
    d3_d4_pref_d = pref_dscores_ses_LM_all.d_score_prefori_d3_d4;
    d3_d4_pref_d_LM_all = [d3_d4_pref_d_LM_all d3_d4_pref_d];

    d1_d3_pref_d = pref_dscores_ses_LM_all.d_score_prefori_d1_d3;
    d1_d3_pref_d_LM_all = [d1_d3_pref_d_LM_all d1_d3_pref_d];
    d1_d4_pref_d = pref_dscores_ses_LM_all.d_score_prefori_d1_d4;
    d1_d4_pref_d_LM_all = [d1_d4_pref_d_LM_all d1_d4_pref_d];

    d1_k = k_max_ses_LM_all.d1_k;
    d1_k_LM_all = [d1_k_LM_all d1_k];
    d2_k = k_max_ses_LM_all.d2_k;
    d2_k_LM_all = [d2_k_LM_all d2_k];
    d3_k = k_max_ses_LM_all.d3_k;
    d3_k_LM_all = [d3_k_LM_all d3_k];
    d4_k = k_max_ses_LM_all.d4_k;
    d4_k_LM_all = [d4_k_LM_all d4_k];

    d1_max = k_max_ses_LM_all.d1_max;
    d1_max_LM_all = [d1_max_LM_all d1_max];
    d2_max = k_max_ses_LM_all.d2_max;
    d2_max_LM_all = [d2_max_LM_all d2_max];
    d3_max = k_max_ses_LM_all.d3_max;
    d3_max_LM_all = [d3_max_LM_all d3_max];
    d4_max = k_max_ses_LM_all.d4_max;
    d4_max_LM_all = [d4_max_LM_all d4_max];

    d1_reliability = reliability_ses_LM_all.d1_reliability;
    d1_reliability_LM_all = [d1_reliability_LM_all d1_reliability];
    d2_reliability = reliability_ses_LM_all.d2_reliability;
    d2_reliability_LM_all = [d2_reliability_LM_all d2_reliability];
    d3_reliability = reliability_ses_LM_all.d3_reliability;
    d3_reliability_LM_all = [d3_reliability_LM_all d3_reliability];
    d4_reliability = reliability_ses_LM_all.d4_reliability;
    d4_reliability_LM_all = [d4_reliability_LM_all d4_reliability];
  
end

%%
fig = figure;
a=subplot(2,2,1)
h=cdfplot(d1_reliability_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d1_reliability_LM_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',1);
legend(['Baseline 1 (V1)'], ['Baseline 1 (LM)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Fit Reliability'])
ylabel(['% of cells'])
title(['Baseline 1 Fit Reliability'])

a=subplot(2,2,2)
h=cdfplot(d2_reliability_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d2_reliability_LM_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',1);
legend(['Baseline 2 (V1)'], ['Baseline 2 (LM)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Fit Reliability'])
ylabel(['% of cells'])
title(['Baseline 2 Fit Reliability'])

a=subplot(2,2,3)
h=cdfplot(d3_reliability_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d3_reliability_LM_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',1);
legend(['Post Dark 15mins (V1)'], ['Post Dark 15mins (LM)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Fit Reliability'])
ylabel(['% of cells'])
title(['Post Dark 15mins Fit Reliability'])

a=subplot(2,2,4)
h=cdfplot(d4_reliability_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d4_reliability_LM_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',1);
legend(['Post Dark 7days (V1)'], ['Post Dark 7days (LM)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Fit Reliability'])
ylabel(['% of cells'])
title(['Post Dark 7days Fit Reliability'])
hold off

print(fullfile(realfnout, ['reliability_by_session_v1_vs_lm.pdf']), '-dpdf', '-bestfit')

%%
fig = figure;
a=subplot(2,2,1)
h=cdfplot(d1_max_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d1_max_LM_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',1);
legend(['Baseline 1 (V1)'], ['Baseline 1 (LM)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F'])
ylabel(['% of cells'])
title(['Baseline 1 Max dF/F'])

a=subplot(2,2,2)
h=cdfplot(d2_max_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d2_max_LM_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',1);
legend(['Baseline 2 (V1)'], ['Baseline 2 (LM)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F'])
ylabel(['% of cells'])
title(['Baseline 2 Max dF/F'])

a=subplot(2,2,3)
h=cdfplot(d3_max_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d3_max_LM_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',1);
legend(['Post Dark 15mins (V1)'], ['Post Dark 15mins (LM)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F'])
ylabel(['% of cells'])
title(['Post Dark 15mins Max dF/F'])

a=subplot(2,2,4)
h=cdfplot(d4_max_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d4_max_LM_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',1);
legend(['Post Dark 7days (V1)'], ['Post Dark 7days (LM)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F'])
ylabel(['% of cells'])
title(['Post Dark 7days Max dF/F'])
hold off

print(fullfile(realfnout, ['max_by_session_v1_vs_lm.pdf']), '-dpdf', '-bestfit')

%%
fig = figure;
a=subplot(2,2,1)
h=cdfplot(d1_k_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d1_k_LM_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',1);
legend(['Baseline 1 (V1)'], ['Baseline 1 (LM)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['K value'])
ylabel(['% of cells'])
title(['Baseline 1 K value'])

a=subplot(2,2,2)
h=cdfplot(d2_k_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d2_k_LM_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',1);
legend(['Baseline 2 (V1)'], ['Baseline 2 (LM)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['K value'])
ylabel(['% of cells'])
title(['Baseline 2 K value'])

a=subplot(2,2,3)
h=cdfplot(d3_k_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d3_k_LM_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',1);
legend(['Post Dark 15mins (V1)'], ['Post Dark 15mins (LM)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['K value'])
ylabel(['% of cells'])
title(['Post Dark 15mins K value'])

a=subplot(2,2,4)
h=cdfplot(d4_k_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d4_k_LM_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',1);
legend(['Post Dark 7days (V1)'], ['Post Dark 7days (LM)'], 'Location', 'southeast')
xlim([0 30])
xlabel(['K value'])
ylabel(['% of cells'])
title(['Post Dark 7days K value'])
hold off

print(fullfile(realfnout, ['k_by_session_v1_vs_lm.pdf']), '-dpdf', '-bestfit')

%%
fig = figure;
a=subplot(2,2,1)
h=cdfplot(d1_d2_pref_d_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d1_d2_pref_d_LM_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',1);
legend(['Base1 - Base2 (V1)'], ['Base1 - Base2 (LM)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Pref Ori Change'])
ylabel(['% of cells'])
title(['Base1 - Base2 Pref Ori Change'])

a=subplot(2,2,2)
h=cdfplot(d1_d3_pref_d_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d1_d3_pref_d_LM_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',1);
legend(['Base1 - Post Dark (15min) (V1)'], ['Base1 - Post Dark (15min) (LM)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Pref Ori Change'])
ylabel(['% of cells'])
title(['Base2 - Post Dark (15min) Pref Ori Change'])

a=subplot(2,2,3)
h=cdfplot(d1_d4_pref_d_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d1_d4_pref_d_LM_all);
set(j, 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth',1);
legend(['Base1 - Post Dark (7days) (V1)'], ['Base1 - Post Dark (7days) (LM)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Pref Ori Change'])
ylabel(['% of cells'])
title(['Post Dark (15min) - Post Dark (7days) Pref Ori Change'])

print(fullfile(realfnout, ['ori_change_by_session_v1_vs_lm.pdf']), '-dpdf', '-bestfit')



%%
%% regression
%%

zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

%***STUCK HERE TRYING TO MAKE CELL IDS***%
% cell_ids = cell(length(d1_d2_pref_d_all),1);
% for i = 1:length(cell_ids)
%     cell_ids{i} = string(i);
% end
% cell_ids = cell_ids';
% cell_ids_all = [cell_ids cell_ids cell_ids];

d1_d2_labels = cell(length(d1_d2_pref_d_all),1);
d1_d2_labels(:) = {'base1_base2'};
d1_d2_labels = d1_d2_labels';

d1_d3_labels = cell(length(d1_d3_pref_d_all),1);
d1_d3_labels(:) = {'base1_postdark15mins'};
d1_d3_labels = d1_d3_labels';

d1_d4_labels = cell(length(d1_d4_pref_d_all),1);
d1_d4_labels(:) = {'base1_postdark1week'};
d1_d4_labels = d1_d4_labels';

all_d_labels = categorical([d1_d2_labels d1_d3_labels d1_d4_labels]');
% all_cell_labels = categorical(cell_ids_all)

alldays_reliability = [d1_reliability_all; d2_reliability_all; d3_reliability_all; d4_reliability_all];
alldays_avg_reliability = nanmean(alldays_reliability);
alldays_max_reliability = max(alldays_reliability);

std_avg_predictors = [zscore(alldays_avg_reliability)]';
std_d1_predictors = [zscor_xnan(d1_reliability_all)]';
std_max_predictors = [zscore(alldays_max_reliability)]';

d1_d2_avg_reliability = nanmean([d1_reliability_all;d2_reliability_all]);
d1_d3_avg_reliability = nanmean([d1_reliability_all;d3_reliability_all]);
d1_d4_avg_reliability = nanmean([d1_reliability_all;d4_reliability_all]);

d1_d2_avg_max = nanmean([d1_max_all;d2_max_all]);
d1_d3_avg_max = nanmean([d1_max_all;d3_max_all]);
d1_d4_avg_max = nanmean([d1_max_all;d4_max_all]);

d1_d2_avg_k = nanmean([d1_k_all;d2_k_all]);
d1_d3_avg_k = nanmean([d1_k_all;d3_k_all]);
d1_d4_avg_k = nanmean([d1_k_all;d4_k_all]);


% v1_std_table = table(all_d_labels, [std_avg_predictors; std_avg_predictors; std_avg_predictors], [std_max_predictors; std_max_predictors; std_max_predictors], [std_d1_predictors; std_d1_predictors; std_d1_predictors], ...
%     zscore([d1_d2_pref_d_all'; d1_d3_pref_d_all'; d1_d4_pref_d_all']), [d1_d2_pref_d_all'; d1_d3_pref_d_all'; d1_d4_pref_d_all'], 'VariableNames', {'timepoint', 'avg_fit_reliability', 'max_fit_reliability', 'd1_fit_reliability', 'prefori_change', 'raw_pref_change'})

v1_std_table = table(all_d_labels, [std_avg_predictors; std_avg_predictors; std_avg_predictors], [std_max_predictors; std_max_predictors; std_max_predictors], [std_d1_predictors; std_d1_predictors; std_d1_predictors], ...
    zscore([d1_d2_pref_d_all'; d1_d3_pref_d_all'; d1_d4_pref_d_all']), 'VariableNames', {'timepoint', 'avg_fit_reliability', 'max_fit_reliability', 'd1_fit_reliability', 'prefori_change'})

%%
%testing which reliability measures gives the best model

v1_avg_std_model = fitlm(v1_std_table, 'prefori_change~avg_fit_reliability') %using average fit gives the best model -> we will use this measure then
v1_avg_std_model = fitlm(v1_std_table, 'prefori_change~max_fit_reliability')
v1_d1_std_model = fitlm(v1_std_table, 'prefori_change~d1_fit_reliability')

v1_avg_std_model = fitlm(v1_std_table, 'prefori_change~timepoint+avg_fit_reliability') 


%%
v1_simple_table = table(all_d_labels, [std_avg_predictors; std_avg_predictors; std_avg_predictors], zscore([d1_d2_pref_d_all'; d1_d3_pref_d_all'; d1_d4_pref_d_all']), 'VariableNames', {'timepoint', 'avg_fit_reliability', 'prefori_change'})
v1_stepwise_model = stepwiselm(v1_simple_table,'interactions')


plotSlice(v1_stepwise_model)
plotEffects(v1_stepwise_model)
plotInteraction(v1_stepwise_model,'timepoint','avg_fit_reliability', 'predictions')


%%

d1_d2_labels_LM = cell(length(d1_d2_pref_d_LM_all),1)
d1_d2_labels_LM(:) = {'base1_base2'}
d1_d2_labels_LM = d1_d2_labels_LM'

d1_d3_labels_LM = cell(length(d1_d3_pref_d_LM_all),1)
d1_d3_labels_LM(:) = {'base1_postdark15mins'}
d1_d3_labels_LM = d1_d3_labels_LM'

d1_d4_labels_LM = cell(length(d1_d4_pref_d_LM_all),1)
d1_d4_labels_LM(:) = {'base1_postdark1week'}
d1_d4_labels_LM = d1_d4_labels_LM'

all_d_labels_LM = categorical([d1_d2_labels_LM d1_d3_labels_LM d1_d4_labels_LM]')

alldays_reliability_LM = [d1_reliability_LM_all; d2_reliability_LM_all; d3_reliability_LM_all; d4_reliability_LM_all];
alldays_avg_reliability_LM = nanmean(alldays_reliability_LM)
alldays_max_reliability_LM = max(alldays_reliability_LM)

std_avg_predictors_LM = [zscore(alldays_avg_reliability_LM)]'
std_d1_predictors_LM = [zscor_xnan(d1_reliability_LM_all)]'
std_max_predictors_LM = [zscore(alldays_max_reliability_LM)]'

d1_d2_avg_reliability_LM = nanmean([d1_reliability_LM_all;d2_reliability_LM_all]);
d1_d3_avg_reliability_LM = nanmean([d1_reliability_LM_all;d3_reliability_LM_all]);
d1_d4_avg_reliability_LM = nanmean([d1_reliability_LM_all;d4_reliability_LM_all]);

d1_d2_avg_max_LM = nanmean([d1_max_LM_all;d2_max_LM_all]);
d1_d3_avg_max_LM = nanmean([d1_max_LM_all;d3_max_LM_all]);
d1_d4_avg_max_LM = nanmean([d1_max_LM_all;d4_max_LM_all]);

d1_d2_avg_k_LM = nanmean([d1_k_LM_all;d2_k_LM_all]);
d1_d3_avg_k_LM = nanmean([d1_k_LM_all;d3_k_LM_all]);
d1_d4_avg_k_LM = nanmean([d1_k_LM_all;d4_k_LM_all]);


% LM_std_table = table(all_d_labels_LM, [alldays_avg_reliability_LM'; alldays_avg_reliability_LM'; alldays_avg_reliability_LM'], [alldays_max_reliability_LM'; alldays_max_reliability_LM'; alldays_max_reliability_LM'], [d1_reliability_LM_all'; d1_reliability_LM_all'; d1_reliability_LM_all'], ...
%     zscore([d1_d2_pref_d_LM_all'; d1_d3_pref_d_LM_all'; d1_d4_pref_d_LM_all']), [d1_d2_pref_d_LM_all'; d1_d3_pref_d_LM_all'; d1_d4_pref_d_LM_all'], 'VariableNames', {'timepoint', 'avg_fit_reliability', 'max_fit_reliability', 'd1_fit_reliability', 'prefori_change', 'raw_pref_change'})

LM_std_table = table(all_d_labels_LM, [std_avg_predictors_LM; std_avg_predictors_LM; std_avg_predictors_LM], [std_max_predictors_LM; std_max_predictors_LM; std_max_predictors_LM], [std_d1_predictors_LM; std_d1_predictors_LM; std_d1_predictors_LM], ...
    zscore([d1_d2_pref_d_LM_all'; d1_d3_pref_d_LM_all'; d1_d4_pref_d_LM_all']), 'VariableNames', {'timepoint', 'avg_fit_reliability', 'max_fit_reliability', 'd1_fit_reliability', 'prefori_change'})

%%

LM_avg_std_model = fitlm(LM_std_table, 'prefori_change~avg_fit_reliability')
LM_max_std_model = fitlm(LM_std_table, 'prefori_change~max_fit_reliability')
LM_d1_std_model = fitlm(LM_std_table, 'prefori_change~d1_fit_reliability')

LM_avg_std_model = fitlm(LM_std_table, 'prefori_change~timepoint+avg_fit_reliability') 


%%
LM_simple_table = table(all_d_labels_LM, [std_avg_predictors_LM; std_avg_predictors_LM; std_avg_predictors_LM], zscore([d1_d2_pref_d_LM_all'; d1_d3_pref_d_LM_all'; d1_d4_pref_d_LM_all']), 'VariableNames', {'timepoint', 'avg_fit_reliability', 'prefori_change'})
LM_stepwise_model = stepwiselm(LM_simple_table,'interactions')

LM_full_std_model = fitlm(LM_simple_table,'interactions','ResponseVar','prefori_change',...
    'PredictorVars',{'timepoint','avg_fit_reliability'},...
    'CategoricalVar',{'timepoint'})

plotSlice(LM_stepwise_model)
plotEffects(LM_stepwise_model)
plotInteraction(LM_stepwise_model,'timepoint','avg_fit_reliability', 'predictions')

plotInteraction(LM_full_std_model,'timepoint','avg_fit_reliability', 'predictions')


%%
%model with group included (imaging region)
v1_group_labels = cell(length(v1_std_table.prefori_change),1);
v1_group_labels(:) = {'v1'};

LM_group_labels = cell(length(LM_std_table.prefori_change),1);
LM_group_labels(:) = {'LM'};


merge_v1_table = table(v1_group_labels, all_d_labels, [alldays_avg_reliability'; alldays_avg_reliability'; alldays_avg_reliability'], ...
    [d1_d2_pref_d_all'; d1_d3_pref_d_all'; d1_d4_pref_d_all'], 'VariableNames', {'group', 'timepoint', 'avg_fit_reliability', 'prefori_change'});

merge_lm_table = table(LM_group_labels, all_d_labels_LM, [alldays_avg_reliability_LM'; alldays_avg_reliability_LM'; alldays_avg_reliability_LM'], ...
    [d1_d2_pref_d_LM_all'; d1_d3_pref_d_LM_all'; d1_d4_pref_d_LM_all'], 'VariableNames', {'group', 'timepoint', 'avg_fit_reliability', 'prefori_change'});

merge_full_table = [merge_v1_table; merge_lm_table];

merge_full_table.avg_fit_reliability = zscor_xnan(merge_full_table.avg_fit_reliability);
merge_full_table.prefori_change = zscor_xnan(merge_full_table.prefori_change);

merge_full_table

%%

merge_full_std_model = fitlm(merge_full_table,'interactions','ResponseVar','prefori_change',...
    'PredictorVars',{'group','timepoint','avg_fit_reliability'},...
    'CategoricalVar',{'group','timepoint'})

merge_full_std_step_model = stepwiselm(merge_full_table, 'interactions')

plotEffects(merge_full_std_model)
plotInteraction(merge_full_std_model,'group','timepoint','predictions')
plotInteraction(merge_full_std_model,'group','avg_fit_reliability','predictions')



%%
%residuals

figure;
histogram(v1_stepwise_model.Residuals.Raw)

figure;
histogram(LM_stepwise_model.Residuals.Raw)

figure;
histogram(merge_full_std_model.Residuals.Raw)


%%
lme = fitlme(merge_full_table,'prefori_change ~ 1 + group + avg_fit_reliability + (1|timepoint)')

lme2 = fitlme(merge_full_table,'prefori_change~group+avg_fit_reliability+(avg_fit_reliability|timepoint)');



%%
%% *UPDATED/IMPROVED MODEL WITH AVGS OF 2 TIMEPOINTS
%%

new_v1_table = table(all_d_labels, zscor_xnan([d1_d2_avg_reliability'; d1_d3_avg_reliability'; d1_d4_avg_reliability']), ...
    zscor_xnan([d1_d2_avg_max'; d1_d3_avg_max'; d1_d4_avg_max']), zscor_xnan([d1_d2_avg_k'; d1_d3_avg_k'; d1_d4_avg_k']), ...
    zscor_xnan([d1_d2_pref_d_all'; d1_d3_pref_d_all'; d1_d4_pref_d_all']), ...
    'VariableNames', {'timepoint', 'avg_reliability', 'avg_max_dfof', 'avg_k_val', 'prefori_change'})

new_v1_avg_std_model = fitlm(new_v1_table, 'prefori_change~timepoint+avg_reliability+avg_max_dfof+avg_k_val') 

new_v1_interaction_avg_std_model = fitlm(new_v1_table, 'prefori_change~timepoint+avg_reliability+avg_max_dfof+avg_k_val+timepoint*avg_reliability') 
plotInteraction(new_v1_interaction_avg_std_model,'timepoint','avg_reliability','predictions')

new_v1_step_model = stepwiselm(new_v1_table, 'interactions')

%%
new_LM_table = table(all_d_labels_LM, zscor_xnan([d1_d2_avg_reliability_LM'; d1_d3_avg_reliability_LM'; d1_d4_avg_reliability_LM']), ...
    zscor_xnan([d1_d2_avg_max_LM'; d1_d3_avg_max_LM'; d1_d4_avg_max_LM']), zscor_xnan([d1_d2_avg_k_LM'; d1_d3_avg_k_LM'; d1_d4_avg_k_LM']), ...
    zscor_xnan([d1_d2_pref_d_LM_all'; d1_d3_pref_d_LM_all'; d1_d4_pref_d_LM_all']), ...
    'VariableNames', {'timepoint', 'avg_reliability', 'avg_max_dfof', 'avg_k_val', 'prefori_change'})

new_LM_avg_std_model = fitlm(new_LM_table, 'prefori_change~timepoint+avg_reliability+avg_max_dfof+avg_k_val') 

new_LM_interaction_avg_std_model = fitlm(new_LM_table, 'prefori_change~timepoint+avg_reliability+avg_max_dfof+avg_k_val+timepoint*avg_reliability') 
plotInteraction(new_LM_interaction_avg_std_model,'timepoint','avg_reliability','predictions')

new_LM_step_model = stepwiselm(new_LM_table, 'interactions')


LM_interaction_only_model = fitlm(new_LM_table, 'prefori_change~timepoint:avg_reliability') 

%%
%group and timepoint effect on change ori

new_v1_merge_table = table(v1_group_labels, all_d_labels, [d1_d2_avg_reliability'; d1_d3_avg_reliability'; d1_d4_avg_reliability'], ...
    [d1_d2_avg_max'; d1_d3_avg_max'; d1_d4_avg_max'], [d1_d2_avg_k'; d1_d3_avg_k'; d1_d4_avg_k'], ...
    [d1_d2_pref_d_all'; d1_d3_pref_d_all'; d1_d4_pref_d_all'], ...
    'VariableNames', {'group', 'timepoint', 'avg_reliability', 'avg_max_dfof', 'avg_k_val', 'prefori_change'})

new_LM_merge_table = table(LM_group_labels, all_d_labels_LM, [d1_d2_avg_reliability_LM'; d1_d3_avg_reliability_LM'; d1_d4_avg_reliability_LM'], ...
    [d1_d2_avg_max_LM'; d1_d3_avg_max_LM'; d1_d4_avg_max_LM'], [d1_d2_avg_k_LM'; d1_d3_avg_k_LM'; d1_d4_avg_k_LM'], ...
    [d1_d2_pref_d_LM_all'; d1_d3_pref_d_LM_all'; d1_d4_pref_d_LM_all'], ...
    'VariableNames', {'group', 'timepoint', 'avg_reliability', 'avg_max_dfof', 'avg_k_val', 'prefori_change'})

new_full_merge_table = [new_v1_merge_table; new_LM_merge_table]

new_full_merge_table.avg_reliability = zscor_xnan(new_full_merge_table.avg_reliability);
new_full_merge_table.avg_max_dfof = zscor_xnan(new_full_merge_table.avg_max_dfof);
new_full_merge_table.avg_k_val = zscor_xnan(new_full_merge_table.avg_k_val);
new_full_merge_table.prefori_change = zscor_xnan(new_full_merge_table.prefori_change);

new_full_merge_table


new_merge_ori_model = fitlm(new_full_merge_table, 'prefori_change~group*timepoint')
plotInteraction(new_merge_ori_model,'timepoint','group','predictions')
plotInteraction(new_merge_ori_model,'group','timepoint','predictions')

%% setting up new tables for d1-d4 values
d1_labels = cell(length(d1_reliability_all),1);
d1_labels(:) = {'base1'};
d1_labels = d1_labels';

d2_labels = cell(length(d2_reliability_all),1);
d2_labels(:) = {'base2'};
d2_labels = d2_labels';

d3_labels = cell(length(d3_reliability_all),1);
d3_labels(:) = {'postdark_15mins'};
d3_labels = d3_labels';

d4_labels = cell(length(d4_reliability_all),1);
d4_labels(:) = {'postdark_1week'};
d4_labels = d4_labels';

eachday_d_label = categorical([d1_labels d2_labels d3_labels d4_labels]');

eachday_group_label = cell(length(eachday_d_label),1);
eachday_group_label(:) = {'v1'};
eachday_group_label = eachday_group_label';

d1_labels_LM = cell(length(d1_reliability_LM_all),1);
d1_labels_LM(:) = {'base1'};
d1_labels_LM = d1_labels_LM';

d2_labels_LM = cell(length(d2_reliability_LM_all),1);
d2_labels_LM(:) = {'base2'};
d2_labels_LM = d2_labels_LM';

d3_labels_LM = cell(length(d3_reliability_LM_all),1);
d3_labels_LM(:) = {'postdark_15mins'};
d3_labels_LM = d3_labels_LM';

d4_labels_LM = cell(length(d4_reliability_LM_all),1);
d4_labels_LM(:) = {'postdark_1week'};
d4_labels_LM = d4_labels_LM';

eachday_d_label_LM = categorical([d1_labels_LM d2_labels_LM d3_labels_LM d4_labels_LM]');

eachday_group_label_LM = cell(length(eachday_d_label_LM),1);
eachday_group_label_LM(:) = {'LM'};
eachday_group_label_LM = eachday_group_label_LM';

new_eachday_merge_table = table([eachday_group_label'; eachday_group_label_LM'], [eachday_d_label; eachday_d_label_LM], ...
    [d1_reliability_all'; d2_reliability_all'; d3_reliability_all'; d4_reliability_all'; d1_reliability_LM_all'; d2_reliability_LM_all'; d3_reliability_LM_all'; d4_reliability_LM_all';], ...
    [d1_max_all'; d2_max_all'; d3_max_all'; d4_max_all'; d1_max_LM_all'; d2_max_LM_all'; d3_max_LM_all'; d4_max_LM_all';], ...
    [d1_k_all'; d2_k_all'; d3_k_all'; d4_k_all'; d1_k_LM_all'; d2_k_LM_all'; d3_k_LM_all'; d4_k_LM_all';], ...
    'VariableNames', {'group', 'timepoint', 'fit_reliability', 'max_dfof', 'k_val'})

new_eachday_merge_table.fit_reliability = zscor_xnan(new_eachday_merge_table.fit_reliability);
new_eachday_merge_table.max_dfof = zscor_xnan(new_eachday_merge_table.max_dfof);
new_eachday_merge_table.k_val = zscor_xnan(new_eachday_merge_table.k_val);

new_eachday_merge_table

%group and timepoint on fit reliability
new_merge_reliability_model = fitlm(new_eachday_merge_table, 'fit_reliability~group*timepoint')
plotInteraction(new_merge_reliability_model,'group','timepoint','predictions')

%group and timepoint on max df/f
new_merge_max_model = fitlm(new_eachday_merge_table, 'max_dfof~group*timepoint')
plotInteraction(new_merge_max_model,'group','timepoint','predictions')

%group and timepoint on k_val
new_merge_k_model = fitlm(new_eachday_merge_table, 'k_val~group*timepoint')
plotInteraction(new_merge_k_model,'group','timepoint','predictions')

%%
%% huge model

new_full_merge_table

new_merge_full_std_model = fitlm(new_full_merge_table,'interactions','ResponseVar','prefori_change',...
    'PredictorVars',{'group','timepoint','avg_reliability', 'avg_max_dfof', 'avg_k_val'},...
    'CategoricalVar',{'group','timepoint'})

new_merge_full_std_step_model = stepwiselm(new_full_merge_table, 'interactions')

new_merge_full_std_model_simpler = fitlm(new_full_merge_table, 'prefori_change~group+timepoint+avg_reliability+avg_max_dfof+avg_k_val')