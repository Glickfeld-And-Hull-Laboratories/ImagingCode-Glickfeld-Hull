%I am honestly not sure why this code was written, so it is going in the defunct folder

%%
%clear everything
clear all
clear all global
clc
close all
%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight'; %folder to save files to
dataset = 'exp_list_darklight_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
d1 = 2; %day 1 in expt list
% d2 = 32; %day 2 in expt list
% d3 = 33; %day 3 in expt list
mouse = expt(d1).mouse; %mouse
ref_str = 'runs-001'; %string on file name to load
ref_str_d1 = ['runs-',expt(d1).runs]; %need to fix this part***
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date_d1 = expt(d1).date; %day 1 (ref day) date
% date_d2 = expt(d2).date; %day 2 date
% date_d3 = expt(d3).date; %day 3 date
img_folder_d1 = expt(d1).runs; %img folder of day 1
% img_folder_d2 = expt(d2).runs; %img folder of day 2
% img_folder_d3 = expt(d3).runs; %img folder of day 3
%%
%load relevant data (tunings and multi day)

%load ori info for each day
d1_ori = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningInfo.mat']));
tuned_d1 = d1_ori.ind_theta90;
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d1_k_tuned_green = d1_k_max.k1(tuned_d1);
d1_max_tuned_green = d1_k_max.max_dfof(tuned_d1);

figure;
h = cdfplot(d1_k_tuned_green);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
xlabel('k Values')
xlim([0 30]);
title('Baseline K Values', mouse);
print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_baseline_k.pdf']), '-dpdf', '-bestfit')

figure;
j = cdfplot(d1_max_tuned_green);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
xlabel('Max dF/F Values')
xlim([0 1]);
title('Baseline Max dF/F Values', mouse);
print(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_baseline_max.pdf']), '-dpdf', '-bestfit')

save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, 'baseline_k_max.mat']), 'd1_k_tuned_green', 'd1_max_tuned_green')

%%
%pooled

clear all
clear all global
clc
close all
%%
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight'; %folder to load files 
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\pooled'; %folder to save files 
dataset = 'exp_list_darklight_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
%%
ref_str = 'runs-001';
d1_ori_all = [];
d1_k_max_all = [];
d1_k_tuned_green_all = [];
d1_max_tuned_green_all = [];

d1 = [1 2];

for isess = d1
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    d1_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    d1_ori_all = [d1_ori_all d1_ori];
    d1_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    d1_tuned = d1_ori.ind_theta90; %only taking tuned cells
    d1_k_max.k1 = d1_k_max.k1(d1_tuned);
    d1_k_max.max_dfof = d1_k_max.max_dfof(d1_tuned);
    d1_k_max_all = [d1_k_max_all d1_k_max];
    d1_k_tuned_green = d1_k_max.k1;
    d1_k_tuned_green_all = [d1_k_tuned_green_all d1_k_tuned_green];
    d1_max_tuned_green = d1_k_max.max_dfof;
    d1_max_tuned_green_all = [d1_max_tuned_green_all d1_max_tuned_green];
%     d1_matches_all = [d1_matches_all d1_matches];
%     arc_d1_green_tuned_index = intersect(arc_d1_matches.greenCells, arc_d1_ori.ind_theta90);
%     arc_d1_red_tuned_index = intersect(arc_d1_matches.redCells, arc_d1_ori.ind_theta90);
%     arc_d1_green_tuned_index_all = [arc_d1_green_tuned_index_all arc_d1_green_tuned_index];
%     arc_d1_red_tuned_index_all = [arc_d1_red_tuned_index_all arc_d1_red_tuned_index];

end

%%
figure; 
h=cdfplot(d1_k_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
xlabel('k Values')
xlim([0 30]);
title('Baseline K Values');
print(fullfile(realfnout, ['pooldarklight_baseline_k.pdf']), '-dpdf', '-bestfit')


figure;
j = cdfplot(d1_max_tuned_green_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
xlabel('Max dF/F Values')
xlim([0 1]);
title('Baseline Max dF/F Values');
print(fullfile(realfnout, ['pooldarklight_baseline_max.pdf']), '-dpdf', '-bestfit')

%%
%ind plots on one graph
len1 = length(d1_k_max_all(1).k1);
len2 = length(d1_k_max_all(2).k1);

figure;
h = cdfplot(d1_k_tuned_green_all(1:len1))
hold on
j = cdfplot(d1_k_tuned_green_all((len1+1):(len1+len2)))
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
set(j, 'LineStyle', '-.', 'Color', 'r', 'LineWidth',2);
xlabel('k Values')
xlim([0 30]);
legend(['i2537'], ['i2538'])
title('Baseline K Values');
hold off

figure;
h = cdfplot(d1_max_tuned_green_all(1:len1))
hold on
j = cdfplot(d1_max_tuned_green_all((len1+1):(len1+len2)))
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
set(j, 'LineStyle', '-.', 'Color', 'r', 'LineWidth',2);
xlabel('Max dF/F Values')
xlim([0 1]);
legend(['i2537'], ['i2538'])
title('Baseline Max dF/F Values');
hold off
