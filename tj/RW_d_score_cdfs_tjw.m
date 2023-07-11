%this script was made to compare pre and post running wheel (RW) data for the same mice on an individual mouse basis
%
%% 
%clear everything
clear all
clear all global
clc
%%
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\Multi_Day_Comparisons'; %folder to save files to
dataset = 'exp_list_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
pre_rw = 22;
post_rw = 37; %day 1 in expt list
mouse = expt(pre_rw).mouse; %mouse
img_area = expt(pre_rw).img_loc{1};
img_layer = expt(pre_rw).img_loc{2};

d_scores_pre_rw = load(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_', 'd_scores.mat']));
d_scores_post_rw = load(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_RW_', 'd_scores.mat']));
d_scores_pre_rw = d_scores_pre_rw.d_scores_all;
d_scores_post_rw = d_scores_post_rw.d_scores_all;


figure;
cdfplot(d_scores_pre_rw)
hold on
cdfplot(d_scores_post_rw)
hold off
xlabel('Change in Pref Ori')
legend('Pre-RW', 'Post-RW', 'Location', 'southeast');
title([mouse, ' ', img_area, ' ', img_layer], 'Interpreter', 'None');

print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_ori_cdf.pdf']), '-dpdf', '-bestfit')
