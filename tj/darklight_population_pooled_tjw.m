%%
%clear everything
clear all
clear all global
clc
close all

%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\multi_day_drift'; %folder to save files to
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\pooled\drift'; %folder to save files to
dataset = 'exp_list_darklight_actual_tjw'; %experiment list to pick files from
eval(dataset); %load dataset


%%

mean_diag_mat_all = [];
mean_off_diag_mat_all = [];
mean_diag_mat_prev_ses_all = [];
mean_off_diag_22point5diff_all = [];
mean_off_diag_45diff_all = [];
mean_off_diag_67point5diff_all = [];
mean_off_diag_90diff_all = [];

list = [1 8 15 22 26 30 34];
for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    diag_matrices = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'diag_matrices']));

    mean_diag_mat = mean(diag_matrices.diag_mat);
    mean_diag_mat_all = [mean_diag_mat_all; mean_diag_mat];
    mean_off_diag_mat = mean(diag_matrices.off_diag_mat);
    mean_off_diag_mat_all = [mean_off_diag_mat_all; mean_off_diag_mat];
    mean_diag_mat_prev_ses = mean(diag_matrices.diag_mat_prev_ses);
    mean_diag_mat_prev_ses_all = [mean_diag_mat_prev_ses_all; mean_diag_mat_prev_ses];
    mean_off_diag_22point5diff = mean(diag_matrices.off_diag_22point5diff);
    mean_off_diag_22point5diff_all = [mean_off_diag_22point5diff_all; mean_off_diag_22point5diff];
    mean_off_diag_45diff = mean(diag_matrices.off_diag_45diff);
    mean_off_diag_45diff_all = [mean_off_diag_45diff_all; mean_off_diag_45diff];
    mean_off_diag_67point5diff = mean(diag_matrices.off_diag_67point5diff);
    mean_off_diag_67point5diff_all = [mean_off_diag_67point5diff_all; mean_off_diag_67point5diff];
    mean_off_diag_90diff = mean(diag_matrices.off_diag_90diff);
    mean_off_diag_90diff_all = [mean_off_diag_90diff_all; mean_off_diag_90diff];

  
end





%%
figure;
errorbar([2,4,6], mean(mean_diag_mat_all), (std(mean_diag_mat_all)/sqrt(length(mean_diag_mat_all))), '-o', 'Color', 'red', "MarkerSize", 9)
hold on;
errorbar([1,3,5,7], mean(mean_off_diag_22point5diff_all), (std(mean_off_diag_22point5diff_all)/sqrt(length(mean_off_diag_22point5diff_all))), '-o', 'Color', [0 0 0], "MarkerSize", 9)
errorbar([1,3,5,7], mean(mean_off_diag_45diff_all), (std(mean_off_diag_45diff_all)/sqrt(length(mean_off_diag_45diff_all))), '-diamond', 'Color', [.3 .3 .3], "MarkerSize", 9)
errorbar([1,3,5,7], mean(mean_off_diag_67point5diff_all), (std(mean_off_diag_67point5diff_all)/sqrt(length(mean_off_diag_67point5diff_all))), '-^', 'Color', [.6 .6 .6], "MarkerSize", 9)
errorbar([1,3,5,7], mean(mean_off_diag_90diff_all), (std(mean_off_diag_90diff_all)/sqrt(length(mean_off_diag_90diff_all))), '-square', 'Color', [.9 .9 .9], "MarkerSize", 10)
xlim([0 8])
xticks([0:8])
xticklabels({'','Base1 w/self','Base1 to Base2','Base2 w/self','Base1 to Post-dark','Post-dark w/self','Base1 to Post-dark+7d','Post-dark+7d w/self',''})
legend(['Same ori (b/w session)'],['22.5 degrees away (w/i session)'], ['45 degrees away (w/i session)'], ['67.5 degrees away (w/i session)'],  ['90 degrees away (w/i session)'], 'FontSize', 12, 'Location', 'Best')
ylabel('Population Correlation')
print(fullfile(realfnout, ['drift_scores_bw_wi.pdf']), '-dpdf', '-bestfit')


%%
figure;
errorbar(mean(mean_diag_mat_all), (std(mean_diag_mat_all)/sqrt(length(mean_diag_mat_all))))
xlim([0 4])
ylim([0 1])
hold on
errorbar(mean(mean_diag_mat_prev_ses_all), (std(mean_diag_mat_prev_ses_all)/sqrt(length(mean_diag_mat_prev_ses_all))))
legend(['Same ori compared to base1'],['Same ori compared to prev session'], 'FontSize', 12, 'Location', 'Best')
print(fullfile(realfnout, ['drift_scores_ses1_vs_prevses.pdf']), '-dpdf', '-bestfit')
