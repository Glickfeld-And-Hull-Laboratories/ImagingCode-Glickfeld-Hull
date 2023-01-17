%clear everything
clear all
clear all global
clc
close all

%need to do avg tuning curve for postlight1 and 2
%%
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight'; %folder to load files 
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\pooled'; %folder to save files 
dataset = 'exp_list_darklight_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
base1_ori_all = [];
base1_k_max_all = [];
base1_k_tuned_green_all = [];
base1_max_tuned_green_all = [];
base1_mouse_id_all = [];
base1_avgtuningcurve_all = [];


base1 = [1 2];

for isess = base1
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    base1_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    base1_ori_all = [base1_ori_all base1_ori];
    base1_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    base1_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str 'newAvg.mat']));
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



%%
figure; 
h=cdfplot(base1_k_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
xlabel('k Values')
xlim([0 30]);
title('Baseline K Values');
print(fullfile(realfnout, ['pooldarklight_baseline_k.pdf']), '-dpdf', '-bestfit')


figure;
j = cdfplot(base1_max_tuned_green_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
xlabel('Max dF/F Values')
xlim([0 1]);
title('Baseline Max dF/F Values');
print(fullfile(realfnout, ['pooldarklight_baseline_max.pdf']), '-dpdf', '-bestfit')

%%
%ind plots on one graph
len1 = length(base1_k_max_all(1).k1);
len2 = length(base1_k_max_all(2).k1);

figure;
h = cdfplot(base1_k_tuned_green_all(1:len1));
hold on
j = cdfplot(base1_k_tuned_green_all((len1+1):(len1+len2)));
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
set(j, 'LineStyle', '-.', 'Color', 'r', 'LineWidth',2);
xlabel('k Values')
xlim([0 30]);
legend(['i2537'], ['i2538'])
title('Baseline K Values');
hold off

figure;
h = cdfplot(base1_max_tuned_green_all(1:len1));
hold on
j = cdfplot(base1_max_tuned_green_all((len1+1):(len1+len2)));
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
set(j, 'LineStyle', '-.', 'Color', 'r', 'LineWidth',2);
xlabel('Max dF/F Values')
xlim([0 1]);
legend(['i2537'], ['i2538'])
title('Baseline Max dF/F Values');
hold off

%%
postdark_ori_all = [];
postdark_k_max_all = [];
postdark_k_tuned_green_all = [];
postdark_max_tuned_green_all = [];
postdark_mouse_id_all = [];
postdark_avgtuningcurve_all = [];

postdark = [3 4];

for isess = postdark
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    postdark_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    postdark_ori_all = [postdark_ori_all postdark_ori];
    postdark_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    postdark_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str 'newAvg.mat']));
    postdark_tuned = postdark_ori.ind_theta90; %only taking tuned cells
    postdark_k_max.k1 = postdark_k_max.k1(postdark_tuned);
    postdark_k_max.max_dfof = postdark_k_max.max_dfof(postdark_tuned);
    postdark_k_max_all = [postdark_k_max_all postdark_k_max];
    postdark_k_tuned_green = postdark_k_max.k1;
    postdark_k_tuned_green_all = [postdark_k_tuned_green_all postdark_k_tuned_green];
    postdark_max_tuned_green = postdark_k_max.max_dfof;
    postdark_max_tuned_green_all = [postdark_max_tuned_green_all postdark_max_tuned_green];
    postdark_ncells = length(postdark_k_tuned_green);
    postdark_mouse_id = string(repmat(mouse,postdark_ncells,1))';
    postdark_mouse_id_all = [postdark_mouse_id_all postdark_mouse_id];
    postdark_avgtuningcurve = postdark_newavg.newAvg(postdark_tuned,:)';
    postdark_avgtuningcurve_all = [postdark_avgtuningcurve_all postdark_avgtuningcurve];
end

postdark_avgtuningcurve_all = [postdark_avgtuningcurve_all]';

%%
figure; 
h=cdfplot(postdark_k_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
xlabel('k Values')
xlim([0 30]);
title('Post Dark K Values');
print(fullfile(realfnout, ['pooldarklight_postdark_k.pdf']), '-dpdf', '-bestfit')


figure;
j = cdfplot(postdark_max_tuned_green_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
xlabel('Max dF/F Values')
xlim([0 1]);
title('Post Dark Max dF/F Values');
print(fullfile(realfnout, ['pooldarklight_postdark_max.pdf']), '-dpdf', '-bestfit')

%%
postlight1_ori_all = [];
postlight1_k_max_all = [];
postlight1_k_tuned_green_all = [];
postlight1_max_tuned_green_all = [];
postlight1_mouse_id_all = [];
postlight1_avgtuningcurve_all = [];

postlight1 = [5 6];

for isess = postlight1
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    postlight1_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    postlight1_ori_all = [postlight1_ori_all postlight1_ori];
    postlight1_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    postlight1_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str 'newAvg.mat']));
    postlight1_tuned = postlight1_ori.ind_theta90; %only taking tuned cells
    postlight1_k_max.k1 = postlight1_k_max.k1(postlight1_tuned);
    postlight1_k_max.max_dfof = postlight1_k_max.max_dfof(postlight1_tuned);
    postlight1_k_max_all = [postlight1_k_max_all postlight1_k_max];
    postlight1_k_tuned_green = postlight1_k_max.k1;
    postlight1_k_tuned_green_all = [postlight1_k_tuned_green_all postlight1_k_tuned_green];
    postlight1_max_tuned_green = postlight1_k_max.max_dfof;
    postlight1_max_tuned_green_all = [postlight1_max_tuned_green_all postlight1_max_tuned_green];
    postlight1_ncells = length(postlight1_k_tuned_green);
    postlight1_mouse_id = string(repmat(mouse,postlight1_ncells,1))';
    postlight1_mouse_id_all = [postlight1_mouse_id_all postlight1_mouse_id];
    postlight1_ncells = length(postlight1_k_tuned_green);
    postlight1_mouse_id = string(repmat(mouse,postlight1_ncells,1))';
    postlight1_avgtuningcurve = postlight1_newavg.newAvg(postlight1_tuned,:)';
    postlight1_avgtuningcurve_all = [postlight1_avgtuningcurve_all postlight1_avgtuningcurve];
end

postlight1_avgtuningcurve_all = [postlight1_avgtuningcurve_all]';

%%
figure; 
h=cdfplot(postlight1_k_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
xlabel('k Values')
xlim([0 30]);
title('Post Light K Values');
print(fullfile(realfnout, ['pooldarklight_postlight_k.pdf']), '-dpdf', '-bestfit')


figure;
j = cdfplot(postlight1_max_tuned_green_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
xlabel('Max dF/F Values')
xlim([0 1]);
title('Post Light Max dF/F Values');
print(fullfile(realfnout, ['pooldarklight_postlight_max.pdf']), '-dpdf', '-bestfit')

%%
postlight2_ori_all = [];
postlight2_k_max_all = [];
postlight2_k_tuned_green_all = [];
postlight2_max_tuned_green_all = [];
postlight2_mouse_id_all = [];
postlight2_avgtuningcurve_all = [];

postlight2 = [7 8];

for isess = postlight2
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    postlight2_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    postlight2_ori_all = [postlight2_ori_all postlight2_ori];
    postlight2_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    postlight2_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str 'newAvg.mat']));
    postlight2_tuned = postlight2_ori.ind_theta90; %only taking tuned cells
    postlight2_k_max.k1 = postlight2_k_max.k1(postlight2_tuned);
    postlight2_k_max.max_dfof = postlight2_k_max.max_dfof(postlight2_tuned);
    postlight2_k_max_all = [postlight2_k_max_all postlight2_k_max];
    postlight2_k_tuned_green = postlight2_k_max.k1;
    postlight2_k_tuned_green_all = [postlight2_k_tuned_green_all postlight2_k_tuned_green];
    postlight2_max_tuned_green = postlight2_k_max.max_dfof;
    postlight2_max_tuned_green_all = [postlight2_max_tuned_green_all postlight2_max_tuned_green];
    postlight2_ncells = length(postlight2_k_tuned_green);
    postlight2_mouse_id = string(repmat(mouse,postlight2_ncells,1))';
    postlight2_mouse_id_all = [postlight2_mouse_id_all postlight2_mouse_id];
    postlight2_ncells = length(postlight2_k_tuned_green);
    postlight2_mouse_id = string(repmat(mouse,postlight2_ncells,1))';
    postlight2_avgtuningcurve = postlight2_newavg.newAvg(postlight2_tuned,:)';
    postlight2_avgtuningcurve_all = [postlight2_avgtuningcurve_all postlight2_avgtuningcurve];
end

postlight2_avgtuningcurve_all = [postlight2_avgtuningcurve_all]';


%%
figure; 
h=cdfplot(postlight2_k_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
xlabel('k Values')
xlim([0 30]);
title('Post Light 1 Week K Values');
print(fullfile(realfnout, ['pooldarklight_postlight1week_k.pdf']), '-dpdf', '-bestfit')


figure;
j = cdfplot(postlight2_max_tuned_green_all);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
xlabel('Max dF/F Values')
xlim([0 1]);
title('Post Light 1 Week Max dF/F Values');
print(fullfile(realfnout, ['pooldarklight_postlight1week_max.pdf']), '-dpdf', '-bestfit')

%%
figure; 
h=cdfplot(base1_k_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(postdark_k_tuned_green_all);
set(j, 'LineStyle', '-', 'Color', [.35 .35 .35], 'LineWidth',2);
k=cdfplot(postlight1_k_tuned_green_all);
set(k, 'LineStyle', '-', 'Color', [.65 .65 .65], 'LineWidth',2);
l=cdfplot(postlight2_k_tuned_green_all);
set(l, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('k Values')
xlim([0 30]);
title('K Values');
legend(['Baseline'], ['Post Dark'], ['Post Light 4 Hours'], ['Post Light 1 Week'])
hold off
print(fullfile(realfnout, ['pooldarklight_allpoints_k.pdf']), '-dpdf', '-bestfit')

figure; 
h=cdfplot(base1_max_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(postdark_max_tuned_green_all);
set(j, 'LineStyle', '-', 'Color', [.35 .35 .35], 'LineWidth',2);
k=cdfplot(postlight1_max_tuned_green_all);
set(k, 'LineStyle', '-', 'Color', [.65 .65 .65], 'LineWidth',2);
l=cdfplot(postlight2_max_tuned_green_all);
set(l, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('Max dF/F Values')
xlim([0 1]);
title('Max dF/F Values');
legend(['Baseline'], ['Post Dark'], ['Post Light 4 Hours'], ['Post Light 1 Week'])
hold off
print(fullfile(realfnout, ['pooldarklight_allpoints_max.pdf']), '-dpdf', '-bestfit')


%%
%indexing - use find
i2537_index_base1 = find(base1_mouse_id_all == 'i2537');
i2538_index_base1 = find(base1_mouse_id_all == 'i2538');
i2537_index_postdark = find(postdark_mouse_id_all == 'i2537');
i2538_index_postdark = find(postdark_mouse_id_all == 'i2538');
i2537_index_postlight1 = find(postlight1_mouse_id_all == 'i2537');
i2538_index_postlight1 = find(postlight1_mouse_id_all == 'i2538');
i2537_index_postlight2 = find(postlight2_mouse_id_all == 'i2537');
i2538_index_postlight2 = find(postlight2_mouse_id_all == 'i2538');

%%
%individual plots
figure; 
subplot(1,2,1)
h=cdfplot(base1_max_tuned_green_all(i2537_index_base1));
set(h, 'LineStyle', '--', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(postdark_max_tuned_green_all(i2537_index_postdark));
set(j, 'LineStyle', '-', 'Color', [.35 .35 .35], 'LineWidth',2);
k=cdfplot(postlight1_max_tuned_green_all(i2537_index_postlight1));
set(k, 'LineStyle', '-', 'Color', [.65 .65 .65], 'LineWidth',2);
l=cdfplot(postlight2_max_tuned_green_all(i2537_index_postlight2));
set(l, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('Max dF/F Values')
xlim([0 1]);
title('Max dF/F Values');
legend(['Baseline'], ['Post Dark'], ['Post Light 4 Hours'], ['Post Light 1 Week'])
hold off

subplot(1,2,2)
h=cdfplot(base1_max_tuned_green_all(i2538_index_base1));
set(h, 'LineStyle', '--', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(postdark_max_tuned_green_all(i2538_index_postdark));
set(j, 'LineStyle', '-', 'Color', [.35 .35 .35], 'LineWidth',2);
k=cdfplot(postlight1_max_tuned_green_all(i2538_index_postlight1));
set(k, 'LineStyle', '-', 'Color', [.65 .65 .65], 'LineWidth',2);
l=cdfplot(postlight2_max_tuned_green_all(i2538_index_postlight2));
set(l, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('Max dF/F Values')
xlim([0 1]);
title('Max dF/F Values');
legend(['Baseline'], ['Post Dark'], ['Post Light 4 Hours'], ['Post Light 1 Week'])
hold off

figure; 
subplot(1,2,1)
h=cdfplot(base1_k_tuned_green_all(i2537_index_base1));
set(h, 'LineStyle', '--', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(postdark_k_tuned_green_all(i2537_index_postdark));
set(j, 'LineStyle', '-', 'Color', [.35 .35 .35], 'LineWidth',2);
k=cdfplot(postlight1_k_tuned_green_all(i2537_index_postlight1));
set(k, 'LineStyle', '-', 'Color', [.65 .65 .65], 'LineWidth',2);
l=cdfplot(postlight2_k_tuned_green_all(i2537_index_postlight2));
set(l, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('K Values')
xlim([0 30]);
title('K Values');
legend(['Baseline'], ['Post Dark'], ['Post Light 4 Hours'], ['Post Light 1 Week'])
hold off

subplot(1,2,2)
h=cdfplot(base1_k_tuned_green_all(i2538_index_base1));
set(h, 'LineStyle', '--', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(postdark_k_tuned_green_all(i2538_index_postdark));
set(j, 'LineStyle', '-', 'Color', [.35 .35 .35], 'LineWidth',2);
k=cdfplot(postlight1_k_tuned_green_all(i2538_index_postlight1));
set(k, 'LineStyle', '-', 'Color', [.65 .65 .65], 'LineWidth',2);
l=cdfplot(postlight2_k_tuned_green_all(i2538_index_postlight2));
set(l, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('K Values')
xlim([0 30]);
title('K Values');
legend(['Baseline'], ['Post Dark'], ['Post Light 4 Hours'], ['Post Light 1 Week'])
hold off




%%
%avg tuning curves
figure;
fast_errbar(1:8,base1_avgtuningcurve_all,1,'color',[0 0 0]);hold on;fast_errbar(1:8,postdark_avgtuningcurve_all,1,'color',[.35 .35 .35]);fast_errbar(1:8,postlight1_avgtuningcurve_all,1,'color',[.65 .65 .65]); fast_errbar(1:8,postlight2_avgtuningcurve_all,1,'color',[.9 .9 .9]);
xlabel('Orientation')
ylabel('dF/F')
legend(['Baseline'], ['Post Dark'], ['Post Light 4 Hours'], ['Post Light 1 Week']);
hold off;
%fix_axes(gcf,16,'Orientation','dF/F'); axis square

%%
%individual avg tuning curves
figure;
subplot(1,2,1)
fast_errbar(1:8,base1_avgtuningcurve_all(i2537_index_base1,:),1,'color',[0 0 0]);hold on;fast_errbar(1:8,postdark_avgtuningcurve_all(i2537_index_postdark,:),1,'color',[.35 .35 .35]);fast_errbar(1:8,postlight1_avgtuningcurve_all(i2537_index_postlight1,:),1,'color',[.65 .65 .65]); fast_errbar(1:8,postlight2_avgtuningcurve_all(i2537_index_postlight2,:),1,'color',[.9 .9 .9]);
subplot(1,2,2)
fast_errbar(1:8,base1_avgtuningcurve_all(i2538_index_base1,:),1,'color',[0 0 0]);hold on;fast_errbar(1:8,postdark_avgtuningcurve_all(i2538_index_postdark,:),1,'color',[.35 .35 .35]);fast_errbar(1:8,postlight1_avgtuningcurve_all(i2538_index_postlight1,:),1,'color',[.65 .65 .65]); fast_errbar(1:8,postlight2_avgtuningcurve_all(i2538_index_postlight2,:),1,'color',[.9 .9 .9]);
legend(['Baseline'], ['Post Dark'], ['Post Light 4 Hours'], ['Post Light 1 Week']);


%%
figure;
fast_errbar(1:8,base1_avgtuningcurve_all(i2538_index_base1,:),1,'color',[0 0 0]);hold on;fast_errbar(1:8,postdark_avgtuningcurve_all(i2538_index_postdark,:),1,'color',[.25 .25 .25]);fast_errbar(1:8,postlight1_avgtuningcurve_all(i2538_index_postlight1,:),1,'color',[.5 .5 .5]); fast_errbar(1:8,postlight2_avgtuningcurve_all(i2538_index_postlight2,:),1,'color',[.75 .75 .75]);

%%
%ori pref dark to light

postdark_ori_all = [];
postlight1_ori_all = [];
postdark_matches_all = [];
postlight1_matches_all = [];
postdark_tc_all = [];
postlight1_tc_all = [];
