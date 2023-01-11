%clear everything
clear all
clear all global
clc
close all


%need to do pref ori and individuals
%%
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\multi_day'; %folder to load files 
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\pooled\update'; %folder to save files 
dataset = 'exp_list_darklight_actual_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%
%rad2deg(acos((1./k) .* log((exp(k)./2) + (exp(-k)./2)))) this is the formula
%to convert k value into full width at half max (FWHM)
%%
%baseline1
base1_ori_all = [];
base1_k_max_all = [];
base1_k_tuned_green_all = [];
base1_max_tuned_green_all = [];
base1_mouse_id_all = [];
base1_avgtuningcurve_all = [];


base1 = [1 8];

for isess = base1
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    base1_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    base1_ori_all = [base1_ori_all base1_ori];
    base1_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    base1_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    base1_tuned = base1_ori.ind_theta90; %only taking tuned cells
    base1_k_max.k1 = base1_k_max.k1(base1_tuned);
    base1_k_max.max_dfof = base1_k_max.max_dfof(base1_tuned);
    base1_k_max_all = [base1_k_max_all base1_k_max];
    base1_k_tuned_green = base1_k_max.k1;
    base1_k_tuned_green_all = [base1_k_tuned_green_all base1_k_tuned_green];
    base1_max_tuned_green = base1_k_max.max_dfof;
    base1_max_tuned_green_all = [base1_max_tuned_green_all base1_max_tuned_green];
    base1_ncells = length(base1_k_tuned_green);
    base1_mouse_id = string(repmat(mouse,base1_ncells,1))';
    base1_mouse_id_all = [base1_mouse_id_all base1_mouse_id];
    base1_avgtuningcurve = base1_newavg.newAvg(base1_tuned,:)';
    base1_avgtuningcurve_all = [base1_avgtuningcurve_all base1_avgtuningcurve];
end

base1_avgtuningcurve_all = [base1_avgtuningcurve_all]';

base1_tuning_width_all = rad2deg(acos((1./base1_k_tuned_green_all) .* log((exp(base1_k_tuned_green_all)./2) + (exp(-base1_k_tuned_green_all)./2))));
%%
%baseline1 + 4d
base1_4d_ori_all = [];
base1_4d_k_max_all = [];
base1_4d_k_tuned_green_all = [];
base1_4d_max_tuned_green_all = [];
base1_4d_mouse_id_all = [];
base1_4d_avgtuningcurve_all = [];


base1_4d = [2 9];

for isess = base1_4d
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    base1_4d_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    base1_4d_ori_all = [base1_4d_ori_all base1_4d_ori];
    base1_4d_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    base1_4d_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    base1_4d_tuned = base1_4d_ori.ind_theta90; %only taking tuned cells
    base1_4d_k_max.k1 = base1_4d_k_max.k1(base1_4d_tuned);
    base1_4d_k_max.max_dfof = base1_4d_k_max.max_dfof(base1_4d_tuned);
    base1_4d_k_max_all = [base1_4d_k_max_all base1_4d_k_max];
    base1_4d_k_tuned_green = base1_4d_k_max.k1;
    base1_4d_k_tuned_green_all = [base1_4d_k_tuned_green_all base1_4d_k_tuned_green];
    base1_4d_max_tuned_green = base1_4d_k_max.max_dfof;
    base1_4d_max_tuned_green_all = [base1_4d_max_tuned_green_all base1_4d_max_tuned_green];
    base1_4d_ncells = length(base1_4d_k_tuned_green);
    base1_4d_mouse_id = string(repmat(mouse,base1_4d_ncells,1))';
    base1_4d_mouse_id_all = [base1_4d_mouse_id_all base1_4d_mouse_id];
    base1_4d_avgtuningcurve = base1_4d_newavg.newAvg(base1_4d_tuned,:)';
    base1_4d_avgtuningcurve_all = [base1_4d_avgtuningcurve_all base1_4d_avgtuningcurve];
end

base1_4d_avgtuningcurve_all = [base1_4d_avgtuningcurve_all]';

base1_4d_tuning_width_all = rad2deg(acos((1./base1_4d_k_tuned_green_all) .* log((exp(base1_4d_k_tuned_green_all)./2) + (exp(-base1_4d_k_tuned_green_all)./2))));

%%
%baseline2
base2_ori_all = [];
base2_k_max_all = [];
base2_k_tuned_green_all = [];
base2_max_tuned_green_all = [];
base2_mouse_id_all = [];
base2_avgtuningcurve_all = [];


base2 = [3 10];

for isess = base2
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    base2_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    base2_ori_all = [base2_ori_all base2_ori];
    base2_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    base2_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    base2_tuned = base2_ori.ind_theta90; %only taking tuned cells
    base2_k_max.k1 = base2_k_max.k1(base2_tuned);
    base2_k_max.max_dfof = base2_k_max.max_dfof(base2_tuned);
    base2_k_max_all = [base2_k_max_all base2_k_max];
    base2_k_tuned_green = base2_k_max.k1;
    base2_k_tuned_green_all = [base2_k_tuned_green_all base2_k_tuned_green];
    base2_max_tuned_green = base2_k_max.max_dfof;
    base2_max_tuned_green_all = [base2_max_tuned_green_all base2_max_tuned_green];
    base2_ncells = length(base2_k_tuned_green);
    base2_mouse_id = string(repmat(mouse,base2_ncells,1))';
    base2_mouse_id_all = [base2_mouse_id_all base2_mouse_id];
    base2_avgtuningcurve = base2_newavg.newAvg(base2_tuned,:)';
    base2_avgtuningcurve_all = [base2_avgtuningcurve_all base2_avgtuningcurve];
end

base2_avgtuningcurve_all = [base2_avgtuningcurve_all]';

base2_tuning_width_all = rad2deg(acos((1./base2_k_tuned_green_all) .* log((exp(base2_k_tuned_green_all)./2) + (exp(-base2_k_tuned_green_all)./2))));

%%
%baseline2 + 4hr
base2_4hr_ori_all = [];
base2_4hr_k_max_all = [];
base2_4hr_k_tuned_green_all = [];
base2_4hr_max_tuned_green_all = [];
base2_4hr_mouse_id_all = [];
base2_4hr_avgtuningcurve_all = [];


base2_4hr = [4 11];

for isess = base2_4hr
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    base2_4hr_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    base2_4hr_ori_all = [base2_4hr_ori_all base2_4hr_ori];
    base2_4hr_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    base2_4hr_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    base2_4hr_tuned = base2_4hr_ori.ind_theta90; %only taking tuned cells
    base2_4hr_k_max.k1 = base2_4hr_k_max.k1(base2_4hr_tuned);
    base2_4hr_k_max.max_dfof = base2_4hr_k_max.max_dfof(base2_4hr_tuned);
    base2_4hr_k_max_all = [base2_4hr_k_max_all base2_4hr_k_max];
    base2_4hr_k_tuned_green = base2_4hr_k_max.k1;
    base2_4hr_k_tuned_green_all = [base2_4hr_k_tuned_green_all base2_4hr_k_tuned_green];
    base2_4hr_max_tuned_green = base2_4hr_k_max.max_dfof;
    base2_4hr_max_tuned_green_all = [base2_4hr_max_tuned_green_all base2_4hr_max_tuned_green];
    base2_4hr_ncells = length(base2_4hr_k_tuned_green);
    base2_4hr_mouse_id = string(repmat(mouse,base2_4hr_ncells,1))';
    base2_4hr_mouse_id_all = [base2_4hr_mouse_id_all base2_4hr_mouse_id];
    base2_4hr_avgtuningcurve = base2_4hr_newavg.newAvg(base2_4hr_tuned,:)';
    base2_4hr_avgtuningcurve_all = [base2_4hr_avgtuningcurve_all base2_4hr_avgtuningcurve];
end

base2_4hr_avgtuningcurve_all = [base2_4hr_avgtuningcurve_all]';

base2_4hr_tuning_width_all = rad2deg(acos((1./base2_4hr_k_tuned_green_all) .* log((exp(base2_4hr_k_tuned_green_all)./2) + (exp(-base2_4hr_k_tuned_green_all)./2))));

%%
%postdark + 15min
post_dark_ori_all = [];
post_dark_k_max_all = [];
post_dark_k_tuned_green_all = [];
post_dark_max_tuned_green_all = [];
post_dark_mouse_id_all = [];
post_dark_avgtuningcurve_all = [];


post_dark = [5 12];

for isess = post_dark
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    post_dark_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    post_dark_ori_all = [post_dark_ori_all post_dark_ori];
    post_dark_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    post_dark_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    post_dark_tuned = post_dark_ori.ind_theta90; %only taking tuned cells
    post_dark_k_max.k1 = post_dark_k_max.k1(post_dark_tuned);
    post_dark_k_max.max_dfof = post_dark_k_max.max_dfof(post_dark_tuned);
    post_dark_k_max_all = [post_dark_k_max_all post_dark_k_max];
    post_dark_k_tuned_green = post_dark_k_max.k1;
    post_dark_k_tuned_green_all = [post_dark_k_tuned_green_all post_dark_k_tuned_green];
    post_dark_max_tuned_green = post_dark_k_max.max_dfof;
    post_dark_max_tuned_green_all = [post_dark_max_tuned_green_all post_dark_max_tuned_green];
    post_dark_ncells = length(post_dark_k_tuned_green);
    post_dark_mouse_id = string(repmat(mouse,post_dark_ncells,1))';
    post_dark_mouse_id_all = [post_dark_mouse_id_all post_dark_mouse_id];
    post_dark_avgtuningcurve = post_dark_newavg.newAvg(post_dark_tuned,:)';
    post_dark_avgtuningcurve_all = [post_dark_avgtuningcurve_all post_dark_avgtuningcurve];
end

post_dark_avgtuningcurve_all = [post_dark_avgtuningcurve_all]';

post_dark_tuning_width_all = rad2deg(acos((1./post_dark_k_tuned_green_all) .* log((exp(post_dark_k_tuned_green_all)./2) + (exp(-post_dark_k_tuned_green_all)./2))));

%%
%post 4hr light
post_light_ori_all = [];
post_light_k_max_all = [];
post_light_k_tuned_green_all = [];
post_light_max_tuned_green_all = [];
post_light_mouse_id_all = [];
post_light_avgtuningcurve_all = [];


post_light = [6 13];

for isess = post_light
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    post_light_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    post_light_ori_all = [post_light_ori_all post_light_ori];
    post_light_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    post_light_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    post_light_tuned = post_light_ori.ind_theta90; %only taking tuned cells
    post_light_k_max.k1 = post_light_k_max.k1(post_light_tuned);
    post_light_k_max.max_dfof = post_light_k_max.max_dfof(post_light_tuned);
    post_light_k_max_all = [post_light_k_max_all post_light_k_max];
    post_light_k_tuned_green = post_light_k_max.k1;
    post_light_k_tuned_green_all = [post_light_k_tuned_green_all post_light_k_tuned_green];
    post_light_max_tuned_green = post_light_k_max.max_dfof;
    post_light_max_tuned_green_all = [post_light_max_tuned_green_all post_light_max_tuned_green];
    post_light_ncells = length(post_light_k_tuned_green);
    post_light_mouse_id = string(repmat(mouse,post_light_ncells,1))';
    post_light_mouse_id_all = [post_light_mouse_id_all post_light_mouse_id];
    post_light_avgtuningcurve = post_light_newavg.newAvg(post_light_tuned,:)';
    post_light_avgtuningcurve_all = [post_light_avgtuningcurve_all post_light_avgtuningcurve];
end

post_light_avgtuningcurve_all = [post_light_avgtuningcurve_all]';

post_light_tuning_width_all = rad2deg(acos((1./post_light_k_tuned_green_all) .* log((exp(post_light_k_tuned_green_all)./2) + (exp(-post_light_k_tuned_green_all)./2))));

%%
%post 4hr light + 7d
post_light_7d_ori_all = [];
post_light_7d_k_max_all = [];
post_light_7d_k_tuned_green_all = [];
post_light_7d_max_tuned_green_all = [];
post_light_7d_mouse_id_all = [];
post_light_7d_avgtuningcurve_all = [];


post_light_7d = [7 14];

for isess = post_light_7d
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    post_light_7d_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    post_light_7d_ori_all = [post_light_7d_ori_all post_light_7d_ori];
    post_light_7d_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    post_light_7d_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    post_light_7d_tuned = post_light_7d_ori.ind_theta90; %only taking tuned cells
    post_light_7d_k_max.k1 = post_light_7d_k_max.k1(post_light_7d_tuned);
    post_light_7d_k_max.max_dfof = post_light_7d_k_max.max_dfof(post_light_7d_tuned);
    post_light_7d_k_max_all = [post_light_7d_k_max_all post_light_7d_k_max];
    post_light_7d_k_tuned_green = post_light_7d_k_max.k1;
    post_light_7d_k_tuned_green_all = [post_light_7d_k_tuned_green_all post_light_7d_k_tuned_green];
    post_light_7d_max_tuned_green = post_light_7d_k_max.max_dfof;
    post_light_7d_max_tuned_green_all = [post_light_7d_max_tuned_green_all post_light_7d_max_tuned_green];
    post_light_7d_ncells = length(post_light_7d_k_tuned_green);
    post_light_7d_mouse_id = string(repmat(mouse,post_light_7d_ncells,1))';
    post_light_7d_mouse_id_all = [post_light_7d_mouse_id_all post_light_7d_mouse_id];
    post_light_7d_avgtuningcurve = post_light_7d_newavg.newAvg(post_light_7d_tuned,:)';
    post_light_7d_avgtuningcurve_all = [post_light_7d_avgtuningcurve_all post_light_7d_avgtuningcurve];
end

post_light_7d_avgtuningcurve_all = [post_light_7d_avgtuningcurve_all]';

post_light_7d_tuning_width_all = rad2deg(acos((1./post_light_7d_k_tuned_green_all) .* log((exp(post_light_7d_k_tuned_green_all)./2) + (exp(-post_light_7d_k_tuned_green_all)./2))));

%%
figure; 
h=cdfplot(base1_k_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(base1_4d_k_tuned_green_all);
set(j, 'LineStyle', '-', 'Color', [.15 .15 .15], 'LineWidth',2);
k=cdfplot(base2_k_tuned_green_all);
set(k, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(base2_4hr_k_tuned_green_all);
set(l, 'LineStyle', '-', 'Color', [.45 .45 .45], 'LineWidth',2);
m=cdfplot(post_dark_k_tuned_green_all);
set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
n=cdfplot(post_light_k_tuned_green_all);
set(n, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
o=cdfplot(post_light_7d_k_tuned_green_all);
set(o, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('k Values')
xlim([0 30]);
title('K Values');
legend(['Baseline1'], ['Baseline1 + 4d'], ['Baseline2'], ['Baseline2 + 4hr'], ['Post Dark 15min'], ['Post Light 4hr'], ['Post Light + 7d'])
hold off
print(fullfile(realfnout, ['pooldarklight_allpoints_k.pdf']), '-dpdf', '-bestfit')

figure; 
h=cdfplot(base1_max_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(base1_4d_max_tuned_green_all);
set(j, 'LineStyle', '-', 'Color', [.15 .15 .15], 'LineWidth',2);
k=cdfplot(base2_max_tuned_green_all);
set(k, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(base2_4hr_max_tuned_green_all);
set(l, 'LineStyle', '-', 'Color', [.45 .45 .45], 'LineWidth',2);
m=cdfplot(post_dark_max_tuned_green_all);
set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
n=cdfplot(post_light_max_tuned_green_all);
set(n, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
o=cdfplot(post_light_7d_max_tuned_green_all);
set(o, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('Max dF/F Values')
xlim([0 1]);
title('Max dF/F Values');
legend(['Baseline1'], ['Baseline1 + 4d'], ['Baseline2'], ['Baseline2 + 4hr'], ['Post Dark 15min'], ['Post Light 4hr'], ['Post Light + 7d'])
hold off
print(fullfile(realfnout, ['pooldarklight_allpoints_max.pdf']), '-dpdf', '-bestfit')

figure; 
h=cdfplot(base1_tuning_width_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(base1_4d_tuning_width_all);
set(j, 'LineStyle', '-', 'Color', [.15 .15 .15], 'LineWidth',2);
k=cdfplot(base2_tuning_width_all);
set(k, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(base2_4hr_tuning_width_all);
set(l, 'LineStyle', '-', 'Color', [.45 .45 .45], 'LineWidth',2);
m=cdfplot(post_dark_tuning_width_all);
set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
n=cdfplot(post_light_tuning_width_all);
set(n, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
o=cdfplot(post_light_7d_tuning_width_all);
set(o, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('Tuning Width Values')
title('Tuning Width Values');
legend(['Baseline1'], ['Baseline1 + 4d'], ['Baseline2'], ['Baseline2 + 4hr'], ['Post Dark 15min'], ['Post Light 4hr'], ['Post Light + 7d'])
hold off
print(fullfile(realfnout, ['pooldarklight_allpoints_tuning_width.pdf']), '-dpdf', '-bestfit')


%%
%avg tuning curves
figure;
fast_errbar(1:8,base1_avgtuningcurve_all,1,'color',[0 0 0]);hold on;fast_errbar(1:8,base1_4d_avgtuningcurve_all,1,'color',[.15 .15 .15]);fast_errbar(1:8,base2_avgtuningcurve_all,1,'color',[.3 .3 .3]); fast_errbar(1:8,base2_4hr_avgtuningcurve_all,1,'color',[.45 .45 .45]); fast_errbar(1:8,post_dark_avgtuningcurve_all,1,'color',[.6 .6 .6]); fast_errbar(1:8,post_light_avgtuningcurve_all,1,'color',[.75 .75 .75]); fast_errbar(1:8,post_light_7d_avgtuningcurve_all,1,'color',[.9 .9 .9]);
xlabel('Orientation')
ylabel('dF/F')
legend(['Baseline1'], ['Baseline1 + 4d'], ['Baseline2'], ['Baseline2 + 4hr'], ['Post Dark 15min'], ['Post Light 4hr'], ['Post Light + 7d'])
hold off;
print(fullfile(realfnout, ['pooldarklight_allpoints_avgtuning.pdf']), '-dpdf', '-bestfit')

%%
prefori_scores = [];
pref_dscores_all = [];
k_scores = [];
max_scores = [];

list = [1 8];
for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    prefori_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes']));
    prefori_scores = [prefori_scores prefori_ses_all];
    pref_dscores_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores']));
    pref_dscores_all = [pref_dscores_all pref_dscores_ses_all];
%     k_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_k_changes']));
%     k_scores = [k_scores k_ses_all];
%     max_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_max_changes']));
%     max_scores = [max_scores max_ses_all];
end

%%
pref_d_1_2_all = [];
pref_d_1_3_all = [];
pref_d_3_4_all = [];
pref_d_4_5_all = [];
pref_d_5_6_all = [];
pref_d_5_7_all = [];


for idata = 1:length(pref_dscores_all)
    pref_d_1_2 = pref_dscores_all(idata).d_score_prefori_1_2_match;
    pref_d_1_2_all = [pref_d_1_2_all pref_d_1_2];
    pref_d_1_3 = pref_dscores_all(idata).d_score_prefori_1_3_match;
    pref_d_1_3_all = [pref_d_1_3_all pref_d_1_3];
    pref_d_3_4 = pref_dscores_all(idata).d_score_prefori_3_4_match;
    pref_d_3_4_all = [pref_d_3_4_all pref_d_3_4];
    pref_d_4_5 = pref_dscores_all(idata).d_score_prefori_4_5_match;
    pref_d_4_5_all = [pref_d_4_5_all pref_d_4_5];
    pref_d_5_6 = pref_dscores_all(idata).d_score_prefori_5_6_match;
    pref_d_5_6_all = [pref_d_5_6_all pref_d_5_6];
    pref_d_5_7 = pref_dscores_all(idata).d_score_prefori_5_7_match;
    pref_d_5_7_all = [pref_d_5_7_all pref_d_5_7];
end

%%
figure; 
h=cdfplot(pref_d_1_2_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(pref_d_1_3_all);
set(j, 'LineStyle', '-', 'Color', [.15 .15 .15], 'LineWidth',2);
k=cdfplot(pref_d_3_4_all);
set(k, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(pref_d_4_5_all);
set(l, 'LineStyle', '-', 'Color', [.45 .45 .45], 'LineWidth',2);
m=cdfplot(pref_d_5_6_all);
set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
n=cdfplot(pref_d_5_7_all);
set(n, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
xlabel('Change in Pref Ori')
title('Change in Pref Ori');
legend(['Baseline1 to Baseline1 + 4d'], ['Baseline1 to Baseline2'], ['Baseline2 to Baseline2 + 4hr'], ['Baseline2 + 4hr to Post Dark'], ['Post Dark to Post Light 4hr'], ['Post Dark to Post Light 7d'])
hold off
print(fullfile(realfnout, ['pooldarklight_allpoints_change_pref_cdf.pdf']), '-dpdf', '-bestfit')

%%
%indexing - use find
i2537_index_base1 = find(base1_mouse_id_all == 'i2537');
i2538_index_base1 = find(base1_mouse_id_all == 'i2538');
i2537_index_base1_4d = find(base1_4d_mouse_id_all == 'i2537');
i2538_index_base1_4d = find(base1_4d_mouse_id_all == 'i2538');
i2537_index_base2 = find(base2_mouse_id_all == 'i2537');
i2538_index_base2 = find(base2_mouse_id_all == 'i2538');
i2537_index_base2_4hr= find(base2_4hr_mouse_id_all == 'i2537');
i2538_index_base2_4hr = find(base2_4hr_mouse_id_all == 'i2538');
i2537_index_post_dark= find(post_dark_mouse_id_all == 'i2537');
i2538_index_post_dark = find(post_dark_mouse_id_all == 'i2538');
i2537_index_post_light= find(post_light_mouse_id_all == 'i2537');
i2538_index_post_light = find(post_light_mouse_id_all == 'i2538');
i2537_index_post_light_7d= find(post_light_7d_mouse_id_all == 'i2537');
i2538_index_post_light_7d = find(post_light_7d_mouse_id_all == 'i2538');

%%
%individual plots

figure; 
subplot(1,2,1)
h=cdfplot(base1_tuning_width_all(i2537_index_base1));
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(base1_4d_tuning_width_all(i2537_index_base1_4d));
set(j, 'LineStyle', '-', 'Color', [.15 .15 .15], 'LineWidth',2);
k=cdfplot(base2_tuning_width_all(i2537_index_base2));
set(k, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(base2_4hr_tuning_width_all(i2537_index_base2_4hr));
set(l, 'LineStyle', '-', 'Color', [.45 .45 .45], 'LineWidth',2);
m=cdfplot(post_dark_tuning_width_all(i2537_index_post_dark));
set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
n=cdfplot(post_light_tuning_width_all(i2537_index_post_light));
set(n, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
o=cdfplot(post_light_7d_tuning_width_all(i2537_index_post_light_7d));
set(o, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('Tuning Width Values')
title('Tuning Width Values');
legend(['Baseline1'], ['Baseline1 + 4d'], ['Baseline2'], ['Baseline2 + 4hr'], ['Post Dark 15min'], ['Post Light 4hr'], ['Post Light + 7d'])

 
subplot(1,2,2)
h=cdfplot(base1_tuning_width_all(i2538_index_base1));
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(base1_4d_tuning_width_all(i2538_index_base1_4d));
set(j, 'LineStyle', '-', 'Color', [.15 .15 .15], 'LineWidth',2);
k=cdfplot(base2_tuning_width_all(i2538_index_base2));
set(k, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(base2_4hr_tuning_width_all(i2538_index_base2_4hr));
set(l, 'LineStyle', '-', 'Color', [.45 .45 .45], 'LineWidth',2);
m=cdfplot(post_dark_tuning_width_all(i2538_index_post_dark));
set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
n=cdfplot(post_light_tuning_width_all(i2538_index_post_light));
set(n, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
o=cdfplot(post_light_7d_tuning_width_all(i2538_index_post_light_7d));
set(o, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('Tuning Width Values')
title('Tuning Width Values');
legend(['Baseline1'], ['Baseline1 + 4d'], ['Baseline2'], ['Baseline2 + 4hr'], ['Post Dark 15min'], ['Post Light 4hr'], ['Post Light + 7d'])
hold off

