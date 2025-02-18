%clear everything - THE BOTTOM OF THIS NEEDS FIXED (SCATTERS)
clear all
clear all global
clc
%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\Arc_greenVred'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\Arc_greenVred_Pooled'; %folder to save files to
dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
%%
ref_str = 'runs-003';
arc_d1_ori_all = [];
arc_d2_ori_all = [];
arc_d3_ori_all = [];
arc_d1_matches_all = [];
arc_d2_matches_all = [];
arc_d3_matches_all = [];
arc_d1_tc_all = [];
arc_d2_tc_all = [];
arc_d3_tc_all = [];
arc_d1_k_all = [];
arc_d1_max_all = [];
arc_d1_k_red_all = [];
arc_d1_k_tuned_red_all = [];
arc_d1_max_red_all = [];
arc_d1_max_tuned_red_all = [];
arc_d1_k_max_all = [];
arc_d2_k_all = [];
arc_d2_max_all = [];
arc_d2_k_max_all = [];
arc_d3_k_all = [];
arc_d3_max_all = [];
arc_d3_k_max_all = [];
arc_d2_k_red_all = [];
arc_d2_max_red_all = [];

lacz_d1_ori_all = [];
lacz_d2_ori_all = [];
lacz_d3_ori_all = [];
lacz_d1_matches_all = [];
lacz_d2_matches_all = [];
lacz_d3_matches_all = [];
lacz_d1_tc_all = [];
lacz_d2_tc_all = [];
lacz_d3_tc_all = [];
lacz_d1_k_all = [];
lacz_d1_max_all = [];
lacz_d1_k_red_all = [];
lacz_d1_k_tuned_red_all = [];
lacz_d1_max_red_all = [];
lacz_d1_max_tuned_red_all = [];
lacz_d1_k_max_all = [];
lacz_d2_k_all = [];
lacz_d2_max_all = [];
lacz_d2_k_max_all = [];
lacz_d3_k_all = [];
lacz_d3_max_all = [];
lacz_d3_k_max_all = [];
lacz_d2_k_red_all = [];
lacz_d2_max_red_all = [];


arc_d1_k_green_all = [];
arc_d1_k_tuned_green_all = [];
arc_d1_max_tuned_green_all = [];
arc_d2_k_green_all = [];
arc_d2_max_green_all = [];

arc_d1_green_tuned_index_all = [];
arc_d1_red_tuned_index_all = [];


lacz_d1_k_green_all = [];
lacz_d1_k_tuned_green_all = [];
lacz_d1_max_tuned_green_all = [];
lacz_d2_k_green_all = [];
lacz_d2_max_green_all = [];

lacz_d1_green_tuned_index_all = [];
lacz_d1_red_tuned_index_all = [];


arc_d1 = [67 70 73 76 79];
arc_d2 = arc_d1+1;
arc_d3 = arc_d2+1;

lacz_d1 = [49 52 55 58 61 64];
lacz_d2 = lacz_d1+1;
lacz_d3 = lacz_d2+1;

arc_prefori_ses_all = [];
arc_prefori_scores = [];
lacz_prefori_ses_all = [];
lacz_prefori_scores = [];
%%

%arc d1
for isess = arc_d1
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    arc_d1_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    arc_d1_ori_all = [arc_d1_ori_all arc_d1_ori];
    arc_d1_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    arc_d1_matches = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'multiday_alignment.mat']));
    arc_d1_k_max_all = [arc_d1_k_max_all arc_d1_k_max];
    arc_d1_k_red = arc_d1_k_max.k1(arc_d1_matches.redCells); %only red cells
    arc_d1_k_red_all = [arc_d1_k_red_all arc_d1_k_red];
    arc_d1_max_red = arc_d1_k_max.max_dfof(arc_d1_matches.redCells); %only red cells
    arc_d1_max_red_all = [arc_d1_max_red_all arc_d1_max_red];
    arc_d1_k_green = arc_d1_k_max.k1(arc_d1_matches.greenCells); %only red cells
    arc_d1_k_green_all = [arc_d1_k_green_all arc_d1_k_green];
    arc_d1_max_red = arc_d1_k_max.max_dfof(arc_d1_matches.redCells); %only red cells
    arc_d1_k_red_tuned_index = intersect(arc_d1_matches.redCells, arc_d1_ori.ind_theta90);
    arc_d1_k_tuned_red = arc_d1_k_max.k1(arc_d1_k_red_tuned_index);
    arc_d1_k_tuned_red_all = [arc_d1_k_tuned_red_all arc_d1_k_tuned_red];
    arc_d1_k_green_tuned_index = intersect(arc_d1_matches.greenCells, arc_d1_ori.ind_theta90);
    arc_d1_k_tuned_green = arc_d1_k_max.k1(arc_d1_k_green_tuned_index);
    arc_d1_k_tuned_green_all = [arc_d1_k_tuned_green_all arc_d1_k_tuned_green];
    arc_d1_max_tuned_red = arc_d1_k_max.max_dfof(arc_d1_k_red_tuned_index);
    arc_d1_max_tuned_red_all = [arc_d1_max_tuned_red_all arc_d1_max_tuned_red];
    arc_d1_max_tuned_green = arc_d1_k_max.max_dfof(arc_d1_k_green_tuned_index);
    arc_d1_max_tuned_green_all = [arc_d1_max_tuned_green_all arc_d1_max_tuned_green];
    arc_d1_matches_all = [arc_d1_matches_all arc_d1_matches];
    arc_d1_green_tuned_index = intersect(arc_d1_matches.greenCells, arc_d1_ori.ind_theta90);
    arc_d1_red_tuned_index = intersect(arc_d1_matches.redCells, arc_d1_ori.ind_theta90);
    arc_d1_green_tuned_index_all = [arc_d1_green_tuned_index_all arc_d1_green_tuned_index];
    arc_d1_red_tuned_index_all = [arc_d1_red_tuned_index_all arc_d1_red_tuned_index];

end
%%
%lacz
for isess = lacz_d1
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    lacz_d1_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    lacz_d1_ori_all = [lacz_d1_ori_all lacz_d1_ori];
    lacz_d1_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    lacz_d1_matches = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'multiday_alignment.mat']));
    lacz_d1_k_max_all = [lacz_d1_k_max_all lacz_d1_k_max];
    lacz_d1_k_red = lacz_d1_k_max.k1(lacz_d1_matches.redCells); %only red cells
    lacz_d1_k_red_all = [lacz_d1_k_red_all lacz_d1_k_red];
    lacz_d1_max_red = lacz_d1_k_max.max_dfof(lacz_d1_matches.redCells); %only red cells
    lacz_d1_max_red_all = [lacz_d1_max_red_all lacz_d1_max_red];
    lacz_d1_k_green = lacz_d1_k_max.k1(lacz_d1_matches.greenCells); %only red cells
    lacz_d1_k_green_all = [lacz_d1_k_green_all lacz_d1_k_green];
    lacz_d1_max_red = lacz_d1_k_max.max_dfof(lacz_d1_matches.redCells); %only red cells
    lacz_d1_k_red_tuned_index = intersect(lacz_d1_matches.redCells, lacz_d1_ori.ind_theta90);
    lacz_d1_k_tuned_red = lacz_d1_k_max.k1(lacz_d1_k_red_tuned_index);
    lacz_d1_k_tuned_red_all = [lacz_d1_k_tuned_red_all lacz_d1_k_tuned_red];
    lacz_d1_k_green_tuned_index = intersect(lacz_d1_matches.greenCells, lacz_d1_ori.ind_theta90);
    lacz_d1_k_tuned_green = lacz_d1_k_max.k1(lacz_d1_k_green_tuned_index);
    lacz_d1_k_tuned_green_all = [lacz_d1_k_tuned_green_all lacz_d1_k_tuned_green];
    lacz_d1_max_tuned_red = lacz_d1_k_max.max_dfof(lacz_d1_k_red_tuned_index);
    lacz_d1_max_tuned_red_all = [lacz_d1_max_tuned_red_all lacz_d1_max_tuned_red];
    lacz_d1_max_tuned_green = lacz_d1_k_max.max_dfof(lacz_d1_k_green_tuned_index);
    lacz_d1_max_tuned_green_all = [lacz_d1_max_tuned_green_all lacz_d1_max_tuned_green];
    lacz_d1_matches_all = [lacz_d1_matches_all lacz_d1_matches];
    lacz_d1_green_tuned_index = intersect(lacz_d1_matches.greenCells, lacz_d1_ori.ind_theta90);
    lacz_d1_red_tuned_index = intersect(lacz_d1_matches.redCells, lacz_d1_ori.ind_theta90);
    lacz_d1_green_tuned_index_all = [lacz_d1_green_tuned_index_all lacz_d1_green_tuned_index];
    lacz_d1_red_tuned_index_all = [lacz_d1_red_tuned_index_all lacz_d1_red_tuned_index];

end
%%
%MAX DF/F
%all mice all (red and green cells)
figure; 
h=cdfplot(arc_d1_max_tuned_red_all);
hold on
j=cdfplot(arc_d1_max_tuned_green_all);
k=cdfplot(lacz_d1_max_tuned_red_all);
l=cdfplot(lacz_d1_max_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
set(k, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
set(l, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
hold off
%legend(['Arc Red (n = ', num2str(length(arc_d1_max_tuned_red_all)), ')'], ['Arc Green (n = ', num2str(length(arc_d1_max_tuned_green_all)), ')'], ...
%    ['LacZ Red (n = ', num2str(length(lacz_d1_max_tuned_red_all)), ')'], ['LacZ Green (n = ', num2str(length(lacz_d1_max_tuned_green_all)), ')'],'Location', 'best')
% legend(['Day 1 v Day 2'; 'Day 1 v Day 3'])
legend('Arc Enhancer KRAB Cells', 'Arc Enhancer non-KRAB Cells', 'LacZ KRAB Cells', 'LacZ non-KRAB Cells', 'Location', 'best')
xlabel('Max dF/F Values')
xlim([0 1]);
title('Day 1 Max dF/F Values');
print(fullfile(newfnout, ['Grace mice all (green and red) cells', '_maxdfof_cdf.pdf']), '-dpdf', '-bestfit')

%all mice red cells
figure; 
h=cdfplot(arc_d1_max_tuned_red_all);
hold on
j=cdfplot(lacz_d1_max_tuned_red_all);
set(h, 'LineStyle', '-', 'Color', 'b', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold off
legend(['Arc Enhancer KRAB Cells (n = ', num2str(length(arc_d1_max_tuned_red_all)), ')'], ...
    ['LacZ KRAB Cells (n = ', num2str(length(lacz_d1_max_tuned_red_all)), ')'],'Location', 'best')
xlabel('Max dF/F Values')
xlim([0 1]);
title('Day 1 Max dF/F Values');
print(fullfile(newfnout, ['Grace mice red cells', '_maxdfof_cdf.pdf']), '-dpdf', '-bestfit')

%all mice green cells 
figure; 
h=cdfplot(arc_d1_max_tuned_green_all);
hold on
j=cdfplot(lacz_d1_max_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', 'b');
set(j, 'LineStyle', '-', 'Color', 'r');
hold off
legend(['Arc Green (n = ', num2str(length(arc_d1_max_tuned_green_all)), ')'], ...
    ['LacZ Green (n = ', num2str(length(lacz_d1_max_tuned_green_all)), ')'],'Location', 'best')
xlabel('Max dF/F Values')
xlim([0 1]);
title('Day 1 Max dF/F Values');
print(fullfile(newfnout, ['Grace mice green cells', '_maxdfof_cdf.pdf']), '-dpdf', '-bestfit')

%arc mice all cells
figure; 
h=cdfplot(arc_d1_max_tuned_red_all);
hold on
j=cdfplot(arc_d1_max_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold off
legend(['Arc Enhancer KRAB Cells (n = ', num2str(length(arc_d1_max_tuned_red_all)), ')'], ...
    ['Arc Enhancer non-KRAB (n = ', num2str(length(arc_d1_max_tuned_green_all)), ')'],'Location', 'best')
xlabel('Max dF/F Values')
xlim([0 1]);
title('Day 1 Max dF/F Values');
print(fullfile(newfnout, ['Grace arc mice all (red and green) cells', '_maxdfof_cdf.pdf']), '-dpdf', '-bestfit')

%lacz mice all cells
figure; 
h=cdfplot(lacz_d1_max_tuned_red_all);
hold on
j=cdfplot(lacz_d1_max_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', 'r');
set(j, 'LineStyle', '-', 'Color', 'g');
hold off
legend(['LacZ Red (n = ', num2str(length(lacz_d1_max_tuned_red_all)), ')'], ...
    ['LacZ Green (n = ', num2str(length(lacz_d1_max_tuned_green_all)), ')'],'Location', 'best')
xlabel('Max dF/F Values')
xlim([0 1]);
title('Day 1 Max dF/F Values');
print(fullfile(newfnout, ['Grace lacz mice all (red and green) cells', '_maxdfof_cdf.pdf']), '-dpdf', '-bestfit')

%%
%K VALS
%all mice all (red and green) cells
figure; 
h=cdfplot(arc_d1_k_tuned_red_all);
hold on
j=cdfplot(arc_d1_k_tuned_green_all);
k=cdfplot(lacz_d1_k_tuned_red_all);
l=cdfplot(lacz_d1_k_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
set(k, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
set(l, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
hold off
legend(['Arc Red (n = ', num2str(length(arc_d1_k_tuned_red_all)), ')'], ['Arc Green (n = ', num2str(length(arc_d1_k_tuned_green_all)), ')'], ...
    ['LacZ Red (n = ', num2str(length(lacz_d1_k_tuned_red_all)), ')'], ['LacZ Green (n = ', num2str(length(lacz_d1_k_tuned_green_all)), ')'],'Location', 'best')
% legend(['Day 1 v Day 2'; 'Day 1 v Day 3'])
legend('Arc Enhancer KRAB Cells', 'Arc Enhancer non-KRAB Cells', 'LacZ KRAB Cells', 'LacZ non-KRAB Cells', 'Location', 'best')
xlabel('k Values')
xlim([0 30]);
title('Day 1 K Values');
print(fullfile(newfnout, ['Grace mice all (green and red) cells', '_k_cdf.pdf']), '-dpdf', '-bestfit')

%all mice red cells
figure; 
h=cdfplot(arc_d1_k_tuned_red_all);
hold on
j=cdfplot(lacz_d1_k_tuned_red_all);
set(h, 'LineStyle', '-', 'Color', 'b', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold off
% legend(['Arc Enhancer KRAB Cells (n = ', num2str(length(arc_d1_k_tuned_red_all)), ')'], ...
%     ['LacZ KRAB Cells (n = ', num2str(length(lacz_d1_k_tuned_red_all)), ')'],'Location', 'best')
xlabel('k Values')
xlim([0 30]);
title('Day 1 k Values');
print(fullfile(newfnout, ['Grace mice red cells', '_k_cdf.pdf']), '-dpdf', '-bestfit')

%all mice green cells 
figure; 
h=cdfplot(arc_d1_k_tuned_green_all);
hold on
j=cdfplot(lacz_d1_k_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', 'b');
set(j, 'LineStyle', '-', 'Color', 'r');
hold off
legend(['Arc Green (n = ', num2str(length(arc_d1_k_tuned_green_all)), ')'], ...
    ['LacZ Green (n = ', num2str(length(lacz_d1_k_tuned_green_all)), ')'],'Location', 'best')
xlabel('k Values')
xlim([0 30]);
title('Day 1 k Values');
print(fullfile(newfnout, ['Grace mice green cells', '_k_cdf.pdf']), '-dpdf', '-bestfit')

%arc mice all cells
figure; 
h=cdfplot(arc_d1_k_tuned_red_all);
hold on
j=cdfplot(arc_d1_k_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold off
% legend(['Arc Red (n = ', num2str(length(arc_d1_k_tuned_red_all)), ')'], ...
%     ['Arc Green (n = ', num2str(length(arc_d1_k_tuned_green_all)), ')'],'Location', 'best')
xlabel('k Values')
xlim([0 30]);
title('Day 1 k Values');
print(fullfile(newfnout, ['Grace arc mice all (red and green) cells', '_k_cdf.pdf']), '-dpdf', '-bestfit')

%lacz mice all cells
figure; 
h=cdfplot(lacz_d1_k_tuned_red_all);
hold on
j=cdfplot(lacz_d1_k_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', 'r');
set(j, 'LineStyle', '-', 'Color', 'g');
hold off
legend(['LacZ Red (n = ', num2str(length(lacz_d1_k_tuned_red_all)), ')'], ...
    ['LacZ Green (n = ', num2str(length(lacz_d1_k_tuned_green_all)), ')'],'Location', 'best')
xlabel('k Values')
xlim([0 30]);
title('Day 1 k Values');
print(fullfile(newfnout, ['Grace lacz mice all (red and green) cells', '_k_cdf.pdf']), '-dpdf', '-bestfit')

%% POSTER

%all mice red cells max
figure; 
h=cdfplot(arc_d1_max_tuned_red_all);
hold on
j=cdfplot(lacz_d1_max_tuned_red_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
hold off
% legend(['Arc Enhancer KRAB Cells (n = ', num2str(length(arc_d1_max_tuned_red_all)), ')'], ...
%     ['LacZ KRAB Cells (n = ', num2str(length(lacz_d1_max_tuned_red_all)), ')'], 'Location', 'best')
% legend(['Arc Red (n = ', num2str(length(arc_d1_max_tuned_red_all)), ')'], ...
%     ['LacZ Red (n = ', num2str(length(lacz_d1_max_tuned_red_all)), ')'],'Location', 'best')
xlabel('Response Amplitude')
ylabel('Proportion of Cells')
xlim([0 1]);
%title('Day 1 Max dF/F Values');
title('')
print(fullfile(newfnout, ['poster grace all mice red cells', '_maxdfof_cdf.pdf']), '-dpdf', '-bestfit')

%all mice red cells k
figure; 
h=cdfplot(arc_d1_k_tuned_red_all);
hold on
j=cdfplot(lacz_d1_k_tuned_red_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
hold off
% legend(['Arc Red (n = ', num2str(length(arc_d1_k_tuned_red_all)), ')'], ...
%     ['LacZ Red (n = ', num2str(length(lacz_d1_k_tuned_red_all)), ')'],'Location', 'best')
xlabel('Tuning Sharpness')
ylabel('Proportion of Cells')
xlim([0 30]);
%title('Day 1 k Values');
title('')
print(fullfile(newfnout, ['poster grace all mice red cells', '_k_cdf.pdf']), '-dpdf', '-bestfit')

%arc mice all cells k
figure; 
h=cdfplot(arc_d1_k_tuned_red_all);
hold on
j=cdfplot(arc_d1_k_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold off
% legend(['Arc Red (n = ', num2str(length(arc_d1_k_tuned_red_all)), ')'], ...
%     ['Arc Green (n = ', num2str(length(arc_d1_k_tuned_green_all)), ')'],'Location', 'best')
xlabel('Tuning Sharpness')
ylabel('Proportion of Cells')
xlim([0 30]);
title('')
%title('Day 1 k Values');
print(fullfile(newfnout, ['poster grace arc mice all cells', '_k_cdf.pdf']), '-dpdf', '-bestfit')

%arc mice all cells max
figure; 
h=cdfplot(arc_d1_max_tuned_red_all);
hold on
j=cdfplot(arc_d1_max_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
hold off
% legend(['Arc Enhancer KRAB Cells (n = ', num2str(length(arc_d1_max_tuned_red_all)), ')'], ...
%     ['Arc Enhancer non-KRAB Cells (n = ', num2str(length(arc_d1_max_tuned_green_all)), ')'],'Location', 'best')
xlabel('Response Amplitude')
ylabel('Proportion of Cells')
xlim([0 1]);
title('')
%title('Day 1 Max dF/F Values');
print(fullfile(newfnout, ['poster grace arc mice all cells', '_maxdfof_cdf.pdf']), '-dpdf', '-bestfit')


%%
%K-S tests
k_s_max_arcredVlaczred = kstest2(arc_d1_max_tuned_red_all, lacz_d1_max_tuned_red_all);
k_s_max_arcredVarcgreen = kstest2(arc_d1_max_tuned_red_all, arc_d1_max_tuned_green_all);
k_s_k_arcredVlaczred = kstest2(arc_d1_k_tuned_red_all, lacz_d1_k_tuned_red_all);
k_s_k_arcredVarcgreen = kstest2(arc_d1_k_tuned_red_all, arc_d1_k_tuned_green_all);


%% 

arc_green_prefori_d1_d2_match_all = [];
arc_green_prefori_d2_match_all = [];
arc_green_prefori_d1_d3_match_all = [];
arc_green_prefori_d3_match_all = [];
arc_red_prefori_d1_d2_match_all = [];
arc_red_prefori_d2_match_all = [];
arc_red_prefori_d1_d3_match_all = [];
arc_red_prefori_d3_match_all = [];

arc_pref_dscores_all = [];
arc_pref_dscores_ses_all = [];

arc_green_pref_d_d1_d2_all = [];
arc_green_pref_d_d1_d3_all = [];
arc_red_pref_d_d1_d2_all = [];
arc_red_pref_d_d1_d3_all = [];


arc_list = [67 73 76];
for iexp = arc_list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    arc_prefori_ses1 = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_ori_changes']));
    if isfile(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_ori_changes_2.mat']))
        arc_prefori_ses2 = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_ori_changes_2']));
    else
        arc_prefori_ses2 = [];
    end
    arc_prefori_ses_all = [arc_prefori_ses1 arc_prefori_ses2];
    arc_prefori_scores = [arc_prefori_scores arc_prefori_ses_all];
    arc_pref_dscores_ses1 = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_d_scores']));
    if isfile(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_d_scores_2.mat']))
        arc_pref_dscores_ses2 = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_d_scores_2']));
    else
        arc_pref_dscores_ses2 = [];
    end
    arc_pref_dscores_ses_all = [arc_pref_dscores_ses1 arc_pref_dscores_ses2];
    arc_pref_dscores_all = [arc_pref_dscores_all arc_pref_dscores_ses_all];
end

for idata = 1:length(arc_prefori_scores)
    arc_green_prefori_d1_d2_match = arc_prefori_scores(idata).green_prefori_d1_d2_match;
    arc_green_prefori_d2_match = arc_prefori_scores(idata).green_prefori_d2_match;
    arc_green_prefori_d1_d3_match = arc_prefori_scores(idata).green_prefori_d1_d3_match;
    arc_green_prefori_d3_match = arc_prefori_scores(idata).green_prefori_d3_match;
    arc_green_prefori_d1_d2_match_all = [arc_green_prefori_d1_d2_match_all arc_green_prefori_d1_d2_match];
    arc_green_prefori_d2_match_all = [arc_green_prefori_d2_match_all arc_green_prefori_d2_match];
    arc_green_prefori_d1_d3_match_all = [arc_green_prefori_d1_d3_match_all arc_green_prefori_d1_d3_match];
    arc_green_prefori_d3_match_all = [arc_green_prefori_d3_match_all arc_green_prefori_d3_match];
    arc_red_prefori_d1_d2_match = arc_prefori_scores(idata).red_prefori_d1_d2_match;
    arc_red_prefori_d2_match = arc_prefori_scores(idata).red_prefori_d2_match;
    arc_red_prefori_d1_d3_match = arc_prefori_scores(idata).red_prefori_d1_d3_match;
    arc_red_prefori_d3_match = arc_prefori_scores(idata).red_prefori_d3_match;
    arc_red_prefori_d1_d2_match_all = [arc_red_prefori_d1_d2_match_all arc_red_prefori_d1_d2_match];
    arc_red_prefori_d2_match_all = [arc_red_prefori_d2_match_all arc_red_prefori_d2_match];
    arc_red_prefori_d1_d3_match_all = [arc_red_prefori_d1_d3_match_all arc_red_prefori_d1_d3_match];
    arc_red_prefori_d3_match_all = [arc_red_prefori_d3_match_all arc_red_prefori_d3_match];
end

for idata = 1:length(arc_pref_dscores_all)
    arc_green_pref_d_d1_d2 = arc_pref_dscores_all(idata).green_d_score_prefori_d1_d2;
    arc_green_pref_d_d1_d2_all = [arc_green_pref_d_d1_d2_all arc_green_pref_d_d1_d2];
    arc_green_pref_d_d1_d3 = arc_pref_dscores_all(idata).green_d_score_prefori_d1_d3;
    arc_green_pref_d_d1_d3_all = [arc_green_pref_d_d1_d3_all arc_green_pref_d_d1_d3];
    arc_red_pref_d_d1_d2 = arc_pref_dscores_all(idata).red_d_score_prefori_d1_d2;
    arc_red_pref_d_d1_d2_all = [arc_red_pref_d_d1_d2_all arc_red_pref_d_d1_d2];
    arc_red_pref_d_d1_d3 = arc_pref_dscores_all(idata).red_d_score_prefori_d1_d3;
    arc_red_pref_d_d1_d3_all = [arc_red_pref_d_d1_d3_all arc_red_pref_d_d1_d3];
end
%% 

lacz_green_prefori_d1_d2_match_all = [];
lacz_green_prefori_d2_match_all = [];
lacz_green_prefori_d1_d3_match_all = [];
lacz_green_prefori_d3_match_all = [];
lacz_red_prefori_d1_d2_match_all = [];
lacz_red_prefori_d2_match_all = [];
lacz_red_prefori_d1_d3_match_all = [];
lacz_red_prefori_d3_match_all = [];

lacz_pref_dscores_all = [];
lacz_pref_dscores_ses_all = [];

lacz_green_pref_d_d1_d2_all = [];
lacz_green_pref_d_d1_d3_all = [];
lacz_red_pref_d_d1_d2_all = [];
lacz_red_pref_d_d1_d3_all = [];


lacz_list = [49 55 61];
for iexp = lacz_list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    lacz_prefori_ses1 = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_ori_changes']));
    if isfile(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_ori_changes_2.mat']))
        lacz_prefori_ses2 = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_ori_changes_2']));
    else
        lacz_prefori_ses2 = [];
    end
    lacz_prefori_ses_all = [lacz_prefori_ses1 lacz_prefori_ses2];
    lacz_prefori_scores = [lacz_prefori_scores lacz_prefori_ses_all];
    lacz_pref_dscores_ses1 = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_d_scores']));
    if isfile(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_d_scores_2.mat']))
        lacz_pref_dscores_ses2 = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_d_scores_2']));
    else
        lacz_pref_dscores_ses2 = [];
    end
    lacz_pref_dscores_ses_all = [lacz_pref_dscores_ses1 lacz_pref_dscores_ses2];
    lacz_pref_dscores_all = [lacz_pref_dscores_all lacz_pref_dscores_ses_all];
end

for idata = 1:length(lacz_prefori_scores)
    lacz_green_prefori_d1_d2_match = lacz_prefori_scores(idata).green_prefori_d1_d2_match;
    lacz_green_prefori_d2_match = lacz_prefori_scores(idata).green_prefori_d2_match;
    lacz_green_prefori_d1_d3_match = lacz_prefori_scores(idata).green_prefori_d1_d3_match;
    lacz_green_prefori_d3_match = lacz_prefori_scores(idata).green_prefori_d3_match;
    lacz_green_prefori_d1_d2_match_all = [lacz_green_prefori_d1_d2_match_all lacz_green_prefori_d1_d2_match];
    lacz_green_prefori_d2_match_all = [lacz_green_prefori_d2_match_all lacz_green_prefori_d2_match];
    lacz_green_prefori_d1_d3_match_all = [lacz_green_prefori_d1_d3_match_all lacz_green_prefori_d1_d3_match];
    lacz_green_prefori_d3_match_all = [lacz_green_prefori_d3_match_all lacz_green_prefori_d3_match];
    lacz_red_prefori_d1_d2_match = lacz_prefori_scores(idata).red_prefori_d1_d2_match;
    lacz_red_prefori_d2_match = lacz_prefori_scores(idata).red_prefori_d2_match;
    lacz_red_prefori_d1_d3_match = lacz_prefori_scores(idata).red_prefori_d1_d3_match;
    lacz_red_prefori_d3_match = lacz_prefori_scores(idata).red_prefori_d3_match;
    lacz_red_prefori_d1_d2_match_all = [lacz_red_prefori_d1_d2_match_all lacz_red_prefori_d1_d2_match];
    lacz_red_prefori_d2_match_all = [lacz_red_prefori_d2_match_all lacz_red_prefori_d2_match];
    lacz_red_prefori_d1_d3_match_all = [lacz_red_prefori_d1_d3_match_all lacz_red_prefori_d1_d3_match];
    lacz_red_prefori_d3_match_all = [lacz_red_prefori_d3_match_all lacz_red_prefori_d3_match];
end

for idata = 1:length(lacz_pref_dscores_all)
    lacz_green_pref_d_d1_d2 = lacz_pref_dscores_all(idata).green_d_score_prefori_d1_d2;
    lacz_green_pref_d_d1_d2_all = [lacz_green_pref_d_d1_d2_all lacz_green_pref_d_d1_d2];
    lacz_green_pref_d_d1_d3 = lacz_pref_dscores_all(idata).green_d_score_prefori_d1_d3;
    lacz_green_pref_d_d1_d3_all = [lacz_green_pref_d_d1_d3_all lacz_green_pref_d_d1_d3];
    lacz_red_pref_d_d1_d2 = lacz_pref_dscores_all(idata).red_d_score_prefori_d1_d2;
    lacz_red_pref_d_d1_d2_all = [lacz_red_pref_d_d1_d2_all lacz_red_pref_d_d1_d2];
    lacz_red_pref_d_d1_d3 = lacz_pref_dscores_all(idata).red_d_score_prefori_d1_d3;
    lacz_red_pref_d_d1_d3_all = [lacz_red_pref_d_d1_d3_all lacz_red_pref_d_d1_d3];
end

%% THIS IS NOT RIGHT - THE CELL NUMBERS ARE WRONG
%READY TO START SCATTERS AND DSCORES
%all arc both days all cells
figure('Position', [400 20 650 700]);
sgtitle(['All Arc Enhancer Mice'], 'Interpreter', 'None');
subplot(2,1,1);
scatter(arc_green_prefori_d1_d2_match_all, arc_green_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 0.8);
hold on
scatter(arc_red_prefori_d1_d2_match_all, arc_red_prefori_d2_match_all,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Green cells (n = ', num2str(length(arc_green_prefori_d1_d2_match_all)), ')'], ['Red cells (n = ', num2str(length(arc_red_prefori_d1_d2_match_all)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(arc_green_prefori_d1_d3_match_all, arc_green_prefori_d3_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 0.8);
hold on
scatter(arc_red_prefori_d1_d3_match_all, arc_red_prefori_d3_match_all,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 3 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
legend(['Green cells (n = ', num2str(length(arc_green_prefori_d1_d3_match_all)), ')'], ['Red cells (n = ', num2str(length(arc_red_prefori_d1_d3_match_all)), ')'], 'Location', 'northwest')
print(fullfile(newfnout, ['Grace arc mice all (green and red) cells', '_prefori_scatter.pdf']), '-dpdf', '-bestfit')

figure();
scatter(arc_green_prefori_d1_d2_match_all, arc_green_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 0.8);
hold on
scatter(arc_red_prefori_d1_d2_match_all, arc_red_prefori_d2_match_all,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Green cells (n = ', num2str(length(arc_green_prefori_d1_d2_match_all)), ')'], ['Red cells (n = ', num2str(length(arc_red_prefori_d1_d2_match_all)), ')'], 'Location', 'northwest')
print(fullfile(newfnout, ['Grace arc mice all (green and red) D2 cells', '_prefori_scatter.pdf']), '-dpdf', '-bestfit')

figure();
scatter(arc_green_prefori_d1_d3_match_all, arc_green_prefori_d3_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 0.8);
hold on
scatter(arc_red_prefori_d1_d3_match_all, arc_red_prefori_d3_match_all,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 3 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Green cells (n = ', num2str(length(arc_green_prefori_d1_d3_match_all)), ')'], ['Red cells (n = ', num2str(length(arc_red_prefori_d1_d3_match_all)), ')'], 'Location', 'northwest')
print(fullfile(newfnout, ['Grace arc mice all (green and red) D3 cells', '_prefori_scatter.pdf']), '-dpdf', '-bestfit')

figure();
scatter(arc_green_prefori_d1_d2_match_all, arc_green_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g', 'LineWidth', 0.8);
hold on
scatter(arc_red_prefori_d1_d2_match_all, arc_red_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'LineWidth', 0.8, 'MarkerFaceAlpha', 0.5);
scatter(arc_green_prefori_d1_d3_match_all, arc_green_prefori_d3_match_all, 'g', 'LineWidth', 0.8);
scatter(arc_red_prefori_d1_d3_match_all, arc_red_prefori_d3_match_all, 'r', 'LineWidth', 0.8);
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day X Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Green Cells Day 1-Day 2 (n = ', num2str(length(arc_green_prefori_d1_d2_match_all)), ')'], ['Red Cells Day 1-Day 2 (n = ', num2str(length(arc_red_prefori_d1_d2_match_all)), ')'], ... 
    ['Green Cells Day 1-Day 3 (n = ', num2str(length(arc_green_prefori_d1_d3_match_all)), ')'], ['Red Cells Day 1-Day 2 (n = ', num2str(length(arc_red_prefori_d1_d2_match_all)), ')'], 'Location', 'best')

%% POSTER
figure();
scatter(arc_red_prefori_d1_d2_match_all, arc_red_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'LineWidth', 0.8);
hold on
scatter(arc_red_prefori_d1_d3_match_all, arc_red_prefori_d3_match_all, 'filled', 'r', 'LineWidth', 0.8, 'MarkerEdgeAlpha', 0.5, 'MarkerFaceAlpha', 0.5);
scatter(arc_green_prefori_d1_d2_match_all, arc_green_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 0.8);
scatter(arc_green_prefori_d1_d3_match_all, arc_green_prefori_d3_match_all, 'filled' ,'g', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day X Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
% legend(['Arc Enhancer KRAB Cells D1-D2 (n = ', num2str(length(arc_red_prefori_d1_d2_match_all)), ')'], ['Arc Enhancer KRAB Cells D1-D3 (n = ', num2str(length(arc_red_prefori_d1_d3_match_all)), ')'], ... 
%     ['Arc Enhancer non-KRAB Cells D1-D2 (n = ', num2str(length(arc_green_prefori_d1_d2_match_all)), ')'], ['Arc Enhancer non-KRAB Cells D1-D3 (n = ', num2str(length(arc_green_prefori_d1_d3_match_all)), ')'], 'Location', 'best')
print(fullfile(newfnout, ['poster grace arc mice both days all cells', '_prefori_scatter.pdf']), '-dpdf', '-bestfit')


%% 
%lacz scatters
figure('Position', [400 20 650 700]);
sgtitle(['All Lacz Mice'], 'Interpreter', 'None');
subplot(2,1,1);
scatter(lacz_green_prefori_d1_d2_match_all, lacz_green_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 0.8);
hold on
scatter(lacz_red_prefori_d1_d2_match_all, lacz_red_prefori_d2_match_all,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Green cells (n = ', num2str(length(lacz_green_prefori_d1_d2_match_all)), ')'], ['Red cells (n = ', num2str(length(lacz_red_prefori_d1_d2_match_all)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(lacz_green_prefori_d1_d3_match_all, lacz_green_prefori_d3_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 0.8);
hold on
scatter(lacz_red_prefori_d1_d3_match_all, lacz_red_prefori_d3_match_all,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 3 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
legend(['Green cells (n = ', num2str(length(lacz_green_prefori_d1_d3_match_all)), ')'], ['Red cells (n = ', num2str(length(lacz_red_prefori_d1_d3_match_all)), ')'], 'Location', 'northwest')
print(fullfile(newfnout, ['Grace lacz mice all (green and red) cells', '_prefori_scatter.pdf']), '-dpdf', '-bestfit')

figure();
scatter(lacz_green_prefori_d1_d2_match_all, lacz_green_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 0.8);
hold on
scatter(lacz_red_prefori_d1_d2_match_all, lacz_red_prefori_d2_match_all,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Green cells (n = ', num2str(length(lacz_green_prefori_d1_d2_match_all)), ')'], ['Red cells (n = ', num2str(length(lacz_red_prefori_d1_d2_match_all)), ')'], 'Location', 'northwest')
print(fullfile(newfnout, ['Grace lacz mice all (green and red) D2 cells', '_prefori_scatter.pdf']), '-dpdf', '-bestfit')


figure();
scatter(lacz_green_prefori_d1_d3_match_all, lacz_green_prefori_d3_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 0.8);
hold on
scatter(lacz_red_prefori_d1_d3_match_all, lacz_red_prefori_d3_match_all,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 3 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Green cells (n = ', num2str(length(lacz_green_prefori_d1_d3_match_all)), ')'], ['Red cells (n = ', num2str(length(lacz_red_prefori_d1_d3_match_all)), ')'], 'Location', 'northwest')
print(fullfile(newfnout, ['Grace lacz mice all (green and red) D3 cells', '_prefori_scatter.pdf']), '-dpdf', '-bestfit')

figure();
scatter(lacz_green_prefori_d1_d2_match_all, lacz_green_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g', 'LineWidth', 0.8);
hold on
scatter(lacz_red_prefori_d1_d2_match_all, lacz_red_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'LineWidth', 0.8, 'MarkerFaceAlpha', 0.5);
scatter(lacz_green_prefori_d1_d3_match_all, lacz_green_prefori_d3_match_all, 'g', 'LineWidth', 0.8);
scatter(lacz_red_prefori_d1_d3_match_all, lacz_red_prefori_d3_match_all, 'r', 'LineWidth', 0.8);
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day X Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Green Cells Day 1-Day 2 (n = ', num2str(length(lacz_green_prefori_d1_d2_match_all)), ')'], ['Red Cells Day 1-Day 2 (n = ', num2str(length(lacz_red_prefori_d1_d2_match_all)), ')'], ... 
    ['Green Cells Day 1-Day 3 (n = ', num2str(length(lacz_green_prefori_d1_d3_match_all)), ')'], ['Red Cells Day 1-Day 2 (n = ', num2str(length(lacz_red_prefori_d1_d2_match_all)), ')'], 'Location', 'best')


%%
%arc vs lacz red
figure('Position', [400 20 650 700]);
sgtitle(['Arc Enhancer v LacZ Mice'], 'Interpreter', 'None');
subplot(2,1,1);
scatter(arc_red_prefori_d1_d2_match_all, arc_red_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 0.8);
hold on
scatter(lacz_red_prefori_d1_d2_match_all, lacz_red_prefori_d2_match_all,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
% legend(['Arc Enhancer (n = ', num2str(length(arc_red_prefori_d1_d2_match_all)), ')'], ['LacZ (n = ', num2str(length(lacz_red_prefori_d1_d2_match_all)), ')'], 'Location', 'northwest')
legend(['Arc Enhancer KRAB Cells'], ['LacZ LacZ KRAB Cells'], 'Location', 'northwest')
subplot(2,1,2);
scatter(arc_red_prefori_d1_d3_match_all, arc_red_prefori_d3_match_all, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 0.8);
hold on
scatter(lacz_red_prefori_d1_d3_match_all, lacz_red_prefori_d3_match_all,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 3 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
%legend(['Arc Enhancer (n = ', num2str(length(arc_red_prefori_d1_d3_match_all)), ')'], ['LacZ (n = ', num2str(length(lacz_red_prefori_d1_d3_match_all)), ')'], 'Location', 'northwest')
legend(['Arc Enhancer KRAB Cells'], ['LacZ LacZ KRAB Cells'], 'Location', 'northwest')
print(fullfile(newfnout, ['Grace mice red cells', '_prefori_scatter.pdf']), '-dpdf', '-bestfit')

%arc v lacz red d1-d2
figure();
scatter(arc_red_prefori_d1_d2_match_all, arc_red_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 0.8);
hold on
scatter(lacz_red_prefori_d1_d2_match_all, lacz_red_prefori_d2_match_all,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
%legend(['Arc Enhancer (n = ', num2str(length(arc_red_prefori_d1_d2_match_all)), ')'], ['LacZ (n = ', num2str(length(lacz_red_prefori_d1_d2_match_all)), ')'], 'Location', 'northwest')
legend(['Arc Enhancer KRAB Cells'], ['LacZ LacZ KRAB Cells'], 'Location', 'northwest')
print(fullfile(newfnout, ['Grace mice red D2 cells', '_prefori_scatter.pdf']), '-dpdf', '-bestfit')

%arc v lacz red d1-d3
figure();
scatter(arc_red_prefori_d1_d3_match_all, arc_red_prefori_d3_match_all, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 0.8);
hold on
scatter(lacz_red_prefori_d1_d3_match_all, lacz_red_prefori_d3_match_all,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 3 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
%legend(['Arc Enhancer (n = ', num2str(length(arc_red_prefori_d1_d3_match_all)), ')'], ['LacZ (n = ', num2str(length(lacz_red_prefori_d1_d3_match_all)), ')'], 'Location', 'northwest')
legend(['Arc Enhancer KRAB Cells'], ['LacZ LacZ KRAB Cells'], 'Location', 'northwest')
print(fullfile(newfnout, ['Grace mice red D3 cells', '_prefori_scatter.pdf']), '-dpdf', '-bestfit')

%%
figure();
scatter(arc_red_prefori_d1_d2_match_all, arc_red_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'LineWidth', 0.8);
hold on
scatter(arc_red_prefori_d1_d3_match_all, arc_red_prefori_d3_match_all, 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'LineWidth', 0.8, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
scatter(lacz_red_prefori_d1_d2_match_all, lacz_red_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
scatter(lacz_red_prefori_d1_d3_match_all, lacz_red_prefori_d3_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 0.8, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day x Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
% legend(['Arc Enhancer KRAB Cells D1-D2 (n = ', num2str(length(arc_red_prefori_d1_d2_match_all)), ')'], ['Arc Enhancer KRAB Cells D1-D3 (n = ', num2str(length(arc_red_prefori_d1_d3_match_all)), ')'], ...
%    ['LacZ KRAB Cells D1-D2 (n = ', num2str(length(lacz_red_prefori_d1_d2_match_all)), ')'], ['LacZ KRAB Cells D1-D3 (n = ', num2str(length(lacz_red_prefori_d1_d3_match_all)), ')'], 'Location', 'best')
%legend(['Arc Promoter KRAB Cells D1-D2'], ['Arc Promoter KRAB Cells D1-D3'], ['LacZ KRAB Cells D1-D2'], ['LacZ KRAB Cells D1-D3'], 'Location', 'best')
print(fullfile(newfnout, ['poster grace all mice red both days cells', '_prefori_scatter.pdf']), '-dpdf', '-bestfit')

%%
%arc vs lacz green
figure('Position', [400 20 650 700]);
sgtitle(['Arc Enhancer v LacZ Mice'], 'Interpreter', 'None');
subplot(2,1,1);
scatter(arc_green_prefori_d1_d2_match_all, arc_green_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 0.8);
hold on
scatter(lacz_green_prefori_d1_d2_match_all, lacz_green_prefori_d2_match_all,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
%legend(['Arc Enhancer (n = ', num2str(length(arc_green_prefori_d1_d2_match_all)), ')'], ['LacZ (n = ', num2str(length(lacz_green_prefori_d1_d2_match_all)), ')'], 'Location', 'northwest')
legend(['Arc Enhancer non-KRAB Cells'], ['LacZ LacZ non-KRAB Cells'], 'Location', 'northwest')
subplot(2,1,2);
scatter(arc_green_prefori_d1_d3_match_all, arc_green_prefori_d3_match_all, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 0.8);
hold on
scatter(lacz_green_prefori_d1_d3_match_all, lacz_green_prefori_d3_match_all,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 3 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
%legend(['Arc Enhancer (n = ', num2str(length(arc_green_prefori_d1_d3_match_all)), ')'], ['LacZ (n = ', num2str(length(lacz_green_prefori_d1_d3_match_all)), ')'], 'Location', 'northwest')
legend(['Arc Enhancer non-KRAB Cells'], ['LacZ LacZ non-KRAB Cells'], 'Location', 'northwest')
print(fullfile(newfnout, ['Grace mice green cells', '_prefori_scatter.pdf']), '-dpdf', '-bestfit')

%arc v lacz green d1-d2
figure();
scatter(arc_green_prefori_d1_d2_match_all, arc_green_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 0.8);
hold on
scatter(lacz_green_prefori_d1_d2_match_all, lacz_green_prefori_d2_match_all,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
%legend(['Arc Enhancer (n = ', num2str(length(arc_green_prefori_d1_d2_match_all)), ')'], ['LacZ (n = ', num2str(length(lacz_green_prefori_d1_d2_match_all)), ')'], 'Location', 'northwest')
legend(['Arc Enhancer non-KRAB Cells'], ['LacZ LacZ non-KRAB Cells'], 'Location', 'northwest')
print(fullfile(newfnout, ['Grace mice green D2 cells', '_prefori_scatter.pdf']), '-dpdf', '-bestfit')

%arc v lacz green d1-d3
figure();
scatter(arc_green_prefori_d1_d3_match_all, arc_green_prefori_d3_match_all, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 0.8);
hold on
scatter(lacz_green_prefori_d1_d3_match_all, lacz_green_prefori_d3_match_all,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 3 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
%legend(['Arc Enhancer (n = ', num2str(length(arc_green_prefori_d1_d3_match_all)), ')'], ['LacZ (n = ', num2str(length(lacz_green_prefori_d1_d3_match_all)), ')'], 'Location', 'northwest')
legend(['Arc Enhancer non-KRAB Cells'], ['LacZ LacZ non-KRAB Cells'], 'Location', 'northwest')
print(fullfile(newfnout, ['Grace mice green D3 cells', '_prefori_scatter.pdf']), '-dpdf', '-bestfit')

%%
%dscores
figure; 
h = cdfplot(arc_green_pref_d_d1_d2_all);
hold on
j = cdfplot(arc_green_pref_d_d1_d3_all);
m = cdfplot(arc_red_pref_d_d1_d2_all);
n = cdfplot(arc_red_pref_d_d1_d3_all);
% p = cdfplot(enh_dscore_prefori_d1_d2_all);
% q = cdfplot(enh_dscore_prefori_d1_d3_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth', 1.0);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth', 1.0);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 1.0);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.0);
% set(p, 'LineStyle', '-', 'Color', 'g');
% set(q, 'LineStyle', '--', 'Color', 'g');
%legend('Arc D1-D2', 'Arc D1-D3', 'LacZ D1-D2', 'LacZ D1-D3', 'Arc-Enh D1-D2', 'Arc-Enh D1-D3', 'Location', 'Best')
legend(['Arc Green Cell D1-D2 (n = ', num2str(length(arc_green_pref_d_d1_d2_all)), ')'],...
    ['Arc Green Cell D1-D3 (n = ', num2str(length(arc_green_pref_d_d1_d3_all)), ')'],...
    ['Arc Red Cell D1-D2 (n = ', num2str(length(arc_red_pref_d_d1_d2_all)), ')'],...
    ['Arc Red Cell D1-D3 (n = ', num2str(length(arc_red_pref_d_d1_d3_all)), ')'], 'Location', 'Best')
%title('Change in Pref Ori Values')
title('')
xlabel('Change in Pref Ori')
xlim([0 90])
hold off
print(fullfile(newfnout, ['Grace arc mice green and red cells', '_change_pref_cdf.pdf']), '-dpdf', '-bestfit')


%%
%lacz
figure; 
h = cdfplot(lacz_green_pref_d_d1_d2_all);
hold on
j = cdfplot(lacz_green_pref_d_d1_d3_all);
m = cdfplot(lacz_red_pref_d_d1_d2_all);
n = cdfplot(lacz_red_pref_d_d1_d3_all);
% p = cdfplot(enh_dscore_prefori_d1_d2_all);
% q = cdfplot(enh_dscore_prefori_d1_d3_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth', 1.0);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth', 1.0);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 1.0);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.0);
% set(p, 'LineStyle', '-', 'Color', 'g');
% set(q, 'LineStyle', '--', 'Color', 'g');
%legend('Lacz D1-D2', 'Lacz D1-D3', 'LacZ D1-D2', 'LacZ D1-D3', 'Lacz-Enh D1-D2', 'Lacz-Enh D1-D3', 'Location', 'Best')
legend(['Lacz Green Cell D1-D2 (n = ', num2str(length(lacz_green_pref_d_d1_d2_all)), ')'],...
    ['Lacz Green Cell D1-D3 (n = ', num2str(length(lacz_green_pref_d_d1_d3_all)), ')'],...
    ['LacZ Red Cell D1-D2 (n = ', num2str(length(lacz_red_pref_d_d1_d2_all)), ')'],...
    ['LacZ Red Cell D1-D3 (n = ', num2str(length(lacz_red_pref_d_d1_d3_all)), ')'], 'Location', 'Best')
%title('Change in Pref Ori Values')
title('')
xlabel('Change in Pref Ori')
xlim([0 90])
hold off
print(fullfile(newfnout, ['Grace lacz mice green and red cells', '_change_pref_cdf.pdf']), '-dpdf', '-bestfit')

%%
%arc v lacz red
figure; 
h = cdfplot(arc_red_pref_d_d1_d2_all);
hold on
j = cdfplot(arc_red_pref_d_d1_d3_all);
m = cdfplot(lacz_red_pref_d_d1_d2_all);
n = cdfplot(lacz_red_pref_d_d1_d3_all);
% p = cdfplot(enh_dscore_prefori_d1_d2_all);
% q = cdfplot(enh_dscore_prefori_d1_d3_all);
set(h, 'LineStyle', '-', 'Color', 'b', 'LineWidth', 2.0);
set(j, 'LineStyle', '--', 'Color', 'b', 'LineWidth', 2.0);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 2.0);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 2.0);
% set(p, 'LineStyle', '-', 'Color', 'g');
% set(q, 'LineStyle', '--', 'Color', 'g');
%legend('Arc D1-D2', 'Arc D1-D3', 'LacZ D1-D2', 'LacZ D1-D3', 'Arc-Enh D1-D2', 'Arc-Enh D1-D3', 'Location', 'Best')
% legend(['Arc Enh D1-D2 (n = ', num2str(length(arc_red_pref_d_d1_d2_all)), ')'],...
%     ['Arc Enh D1-D3 (n = ', num2str(length(arc_red_pref_d_d1_d3_all)), ')'],...
%     ['LacZ D1-D2 (n = ', num2str(length(lacz_red_pref_d_d1_d2_all)), ')'],...
%     ['LacZ D1-D3 (n = ', num2str(length(lacz_red_pref_d_d1_d3_all)), ')'], 'Location', 'Best')
legend(['Arc Enh KRAB Day 1-Day 2'],...
    ['Arc Enh KRAB Day 1-Day 3'],...
    ['LacZ KRAB Day 1-Day 2'],...
    ['LacZ KRAB Day 1-Day 3'], 'Location', 'Best')
%title('Change in Pref Ori Values')
title('')
xlabel('Change in Pref Ori')
xlim([0 90])
hold off
print(fullfile(newfnout, ['Grace all mice red cells', '_change_pref_cdf.pdf']), '-dpdf', '-bestfit')


%%
%arc v lacz green
figure; 
h = cdfplot(arc_green_pref_d_d1_d2_all);
hold on
j = cdfplot(arc_green_pref_d_d1_d3_all);
m = cdfplot(lacz_green_pref_d_d1_d2_all);
n = cdfplot(lacz_green_pref_d_d1_d3_all);
% p = cdfplot(enh_dscore_prefori_d1_d2_all);
% q = cdfplot(enh_dscore_prefori_d1_d3_all);
set(h, 'LineStyle', '-', 'Color', 'b', 'LineWidth', 2.0);
set(j, 'LineStyle', '--', 'Color', 'b', 'LineWidth', 2.0);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 2.0);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 2.0);
% set(p, 'LineStyle', '-', 'Color', 'g');
% set(q, 'LineStyle', '--', 'Color', 'g');
%legend('Arc D1-D2', 'Arc D1-D3', 'LacZ D1-D2', 'LacZ D1-D3', 'Arc-Enh D1-D2', 'Arc-Enh D1-D3', 'Location', 'Best')
% legend(['Arc Enh D1-D2 (n = ', num2str(length(arc_green_pref_d_d1_d2_all)), ')'],...
%     ['Arc Enh D1-D3 (n = ', num2str(length(arc_green_pref_d_d1_d3_all)), ')'],...
%     ['LacZ D1-D2 (n = ', num2str(length(lacz_green_pref_d_d1_d2_all)), ')'],...
%     ['LacZ D1-D3 (n = ', num2str(length(lacz_green_pref_d_d1_d3_all)), ')'], 'Location', 'Best')
legend(['Arc Enh non-KRAB Day 1-Day 2'],...
    ['Arc Enh non-KRAB Day 1-Day 3'],...
    ['LacZ non-KRAB Day 1-Day 2'],...
    ['LacZ non-KRAB Day 1-Day 3'], 'Location', 'Best')
%title('Change in Pref Ori Values')
title('')
xlabel('Change in Pref Ori')
xlim([0 90])
hold off
print(fullfile(newfnout, ['Grace all mice green cells', '_change_pref_cdf.pdf']), '-dpdf', '-bestfit')

%%
concat_dscores_arc_red = [arc_red_pref_d_d1_d2_all arc_red_pref_d_d1_d3_all];
concat_dscores_lacz_red = [lacz_red_pref_d_d1_d2_all lacz_red_pref_d_d1_d3_all];
concat_dscores_arc_green = [arc_green_pref_d_d1_d2_all arc_green_pref_d_d1_d3_all];

figure;
y = cdfplot(concat_dscores_arc_red)
hold on
z = cdfplot(concat_dscores_lacz_red)
%legend('Arc Enhancer KRAB Cells', 'LacZ KRAB Cells', 'Location', 'best')
set(y, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 2.0);
set(z, 'LineStyle', '-', 'Color', 'k', 'LineWidth', 2.0);
title('')
xlabel('Change in Pref Ori')
ylabel('Proportion of Cells')
xlim([0 90])
hold off
print(fullfile(newfnout, ['poster grace all mice red cells', '_change_pref_cdf.pdf']), '-dpdf', '-bestfit')

figure;
yy = cdfplot(concat_dscores_arc_red)
hold on
zz = cdfplot(concat_dscores_arc_green)
%legend('Arc Enhancer KRAB Cells', 'Arc Enhancer non-KRAB Cells', 'Location', 'best')
set(yy, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 2.0);
set(zz, 'LineStyle', '-', 'Color', 'g', 'LineWidth', 2.0);
title('')
xlabel('Change in Pref Ori')
ylabel('Proportion of Cells')
xlim([0 90])
hold off
print(fullfile(newfnout, ['poster grace arc mice all cells', '_change_pref_cdf.pdf']), '-dpdf', '-bestfit')

%%
%prefori dscores k-s tests

k_s_prefd_arcredVlaczred = kstest2(concat_dscores_arc_red, concat_dscores_lacz_red);
k_s_prefd_arcredVarcgreen = kstest2(concat_dscores_arc_red, concat_dscores_arc_green);

