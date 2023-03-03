%clear everything
clear all
clear all global
clc
close all
%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\Multi_Day_Comparisons'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\RW_Pooled_Multi_Day'; %folder to save files to
dataset = 'exp_list_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%
ref_str = 'runs-001';
d1_ori_all = [];
d2_ori_all = [];
d3_ori_all = [];
d1_matches_all = [];
d2_matches_all = [];
d3_matches_all = [];
d1_tc_all = [];
d2_tc_all = [];
d3_tc_all = [];
d1_k_all = [];
d1_max_all = [];
d1_kall = [];
d1_k_tuned_all = [];
d1_maxall = [];
d1_max_tuned_all = [];
d1_k_max_all = [];
d2_k_all = [];
d2_k_tuned_all = [];
d2_max_tuned_all = [];
d2_max_all = [];
d2_k_max_all = [];
d3_k_all = [];
d3_max_all = [];
d3_k_max_all = [];
d3_k_tuned_all = [];
d3_max_tuned_all = [];


RW_d1_ori_all = [];
RW_d2_ori_all = [];
RW_d3_ori_all = [];
RW_d1_matches_all = [];
RW_d2_matches_all = [];
RW_d3_matches_all = [];
RW_d1_tc_all = [];
RW_d2_tc_all = [];
RW_d3_tc_all = [];
RW_d1_k_all = [];
RW_d1_max_all = [];
RW_d1_kall = [];
RW_d1_k_tuned_all = [];
RW_d1_ma_xall = [];
RW_d1_max_tuned_all = [];
RW_d1_k_max_all = [];
RW_d2_k_all = [];
RW_d2_k_tuned_all = [];
RW_d2_max_tuned_all = [];
RW_d2_max_all = [];
RW_d2_k_max_all = [];
RW_d3_k_all = [];
RW_d3_max_all = [];
RW_d3_k_max_all = [];
RW_d3_k_tuned_all = [];
RW_d3_max_tuned_all = [];

d1_avgtuningcurve_all = [];
d2_avgtuningcurve_all = [];
d3_avgtuningcurve_all = [];

RW_d1_avgtuningcurve_all = [];
RW_d2_avgtuningcurve_all = [];
RW_d3_avgtuningcurve_all = [];


d1 = [19 22 40 43];
d2 = d1+1;
d3 = d2+1;

RW_d1 = [34 37 46 49];
RW_d2 = RW_d1+1;
RW_d3 = RW_d2+1;

prefori_ses_all = [];
prefori_scores = [];
RW_prefori_ses_all = [];
RW_prefori_scores = [];

%%

d1_k_all = [];
d1_k_tuned_all = [];
d1_max_all = [];
d1_max_tuned_all = [];

%no RW d1
for isess = d1
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    d1_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    d1_ori_all = [d1_ori_all d1_ori];
    d1_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    d1_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    d1_k_max_all = [d1_k_max_all d1_k_max];
    d1_k = d1_k_max.k1;
    d1_k_all = [d1_k_all d1_k];
    d1_k_tuned_index = d1_ori.ind_theta90;
    d1_k_tuned = d1_k_max.k1(d1_k_tuned_index);
    d1_k_tuned_all = [d1_k_tuned_all d1_k_tuned];
    d1_max = d1_k_max.max_dfof;
    d1_max_all = [d1_max_all d1_max];
    d1_max_tuned_index = d1_ori.ind_theta90;
    d1_max_tuned = d1_k_max.max_dfof(d1_max_tuned_index);
    d1_max_tuned_all = [d1_max_tuned_all d1_max_tuned];
    d1_avgtuningcurve = d1_newavg.newAvg(d1_k_tuned_index,:)';
    d1_avgtuningcurve_all = [d1_avgtuningcurve_all d1_avgtuningcurve];
end

d1_avgtuningcurve_all = [d1_avgtuningcurve_all]';

d1_tuning_width_all = rad2deg(acos((1./d1_k_tuned_all) .* log((exp(d1_k_tuned_all)./2) + (exp(-d1_k_tuned_all)./2))));


%%
for isess = d2
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    d2_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    d2_ori_all = [d2_ori_all d2_ori];
    d2_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    d2_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    d2_k_max_all = [d2_k_max_all d2_k_max];
    d2_k = d2_k_max.k1;
    d2_k_all = [d2_k_all d2_k];
    d2_k_tuned_index = d2_ori.ind_theta90;
    d2_k_tuned = d2_k_max.k1(d2_k_tuned_index);
    d2_k_tuned_all = [d2_k_tuned_all d2_k_tuned];
    d2_max = d2_k_max.max_dfof;
    d2_max_all = [d2_max_all d2_max];
    d2_max_tuned_index = d2_ori.ind_theta90;
    d2_max_tuned = d2_k_max.max_dfof(d2_max_tuned_index);
    d2_max_tuned_all = [d2_max_tuned_all d2_max_tuned];
    d2_avgtuningcurve = d2_newavg.newAvg(d2_k_tuned_index,:)';
    d2_avgtuningcurve_all = [d2_avgtuningcurve_all d2_avgtuningcurve];
end


d2_avgtuningcurve_all = [d2_avgtuningcurve_all]';

d2_tuning_width_all = rad2deg(acos((1./d2_k_tuned_all) .* log((exp(d2_k_tuned_all)./2) + (exp(-d2_k_tuned_all)./2))));


%%
for isess = d3
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    d3_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    d3_ori_all = [d3_ori_all d3_ori];
    d3_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    d3_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    d3_k_max_all = [d3_k_max_all d3_k_max];
    d3_k = d3_k_max.k1;
    d3_k_all = [d3_k_all d3_k];
    d3_k_tuned_index = d3_ori.ind_theta90;
    d3_k_tuned = d3_k_max.k1(d3_k_tuned_index);
    d3_k_tuned_all = [d3_k_tuned_all d3_k_tuned];
    d3_max = d3_k_max.max_dfof;
    d3_max_all = [d3_max_all d3_max];
    d3_max_tuned_index = d3_ori.ind_theta90;
    d3_max_tuned = d3_k_max.max_dfof(d3_max_tuned_index);
    d3_max_tuned_all = [d3_max_tuned_all d3_max_tuned];
    d3_avgtuningcurve = d3_newavg.newAvg(d3_k_tuned_index,:)';
    d3_avgtuningcurve_all = [d3_avgtuningcurve_all d3_avgtuningcurve];
end


d3_avgtuningcurve_all = [d3_avgtuningcurve_all]';

d3_tuning_width_all = rad2deg(acos((1./d3_k_tuned_all) .* log((exp(d3_k_tuned_all)./2) + (exp(-d3_k_tuned_all)./2))));


%%
%RW
for isess = RW_d1
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    RW_d1_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    RW_d1_ori_all = [RW_d1_ori_all RW_d1_ori];
    RW_d1_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    RW_d1_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    RW_d1_k_max_all = [RW_d1_k_max_all RW_d1_k_max];
    RW_d1_k = RW_d1_k_max.k1;
    RW_d1_k_all = [RW_d1_k_all RW_d1_k];
    RW_d1_k_tuned_index = RW_d1_ori.ind_theta90;
    RW_d1_k_tuned = RW_d1_k_max.k1(RW_d1_k_tuned_index);
    RW_d1_k_tuned_all = [RW_d1_k_tuned_all RW_d1_k_tuned];
    RW_d1_max = RW_d1_k_max.max_dfof;
    RW_d1_max_all = [RW_d1_max_all RW_d1_max];
    RW_d1_max_tuned_index = RW_d1_ori.ind_theta90;
    RW_d1_max_tuned = RW_d1_k_max.max_dfof(RW_d1_max_tuned_index);
    RW_d1_max_tuned_all = [RW_d1_max_tuned_all RW_d1_max_tuned];
    RW_d1_avgtuningcurve = RW_d1_newavg.newAvg(RW_d1_k_tuned_index,:)';
    RW_d1_avgtuningcurve_all = [RW_d1_avgtuningcurve_all RW_d1_avgtuningcurve];

end

RW_d1_avgtuningcurve_all = [RW_d1_avgtuningcurve_all]';

RW_d1_tuning_width_all = rad2deg(acos((1./RW_d1_k_tuned_all) .* log((exp(RW_d1_k_tuned_all)./2) + (exp(-RW_d1_k_tuned_all)./2))));

%%
for isess = RW_d2
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    RW_d2_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    RW_d2_ori_all = [RW_d2_ori_all RW_d2_ori];
    RW_d2_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    RW_d2_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    RW_d2_k_max_all = [RW_d2_k_max_all RW_d2_k_max];
    RW_d2_k = RW_d2_k_max.k1;
    RW_d2_k_all = [RW_d2_k_all RW_d2_k];
    RW_d2_k_tuned_index = RW_d2_ori.ind_theta90;
    RW_d2_k_tuned = RW_d2_k_max.k1(RW_d2_k_tuned_index);
    RW_d2_k_tuned_all = [RW_d2_k_tuned_all RW_d2_k_tuned];
    RW_d2_max = RW_d2_k_max.max_dfof;
    RW_d2_max_all = [RW_d2_max_all RW_d2_max];
    RW_d2_max_tuned_index = RW_d2_ori.ind_theta90;
    RW_d2_max_tuned = RW_d2_k_max.max_dfof(RW_d2_max_tuned_index);
    RW_d2_max_tuned_all = [RW_d2_max_tuned_all RW_d2_max_tuned];
    RW_d2_avgtuningcurve = RW_d2_newavg.newAvg(RW_d2_k_tuned_index,:)';
    RW_d2_avgtuningcurve_all = [RW_d2_avgtuningcurve_all RW_d2_avgtuningcurve];
end


RW_d2_avgtuningcurve_all = [RW_d2_avgtuningcurve_all]';

RW_d2_tuning_width_all = rad2deg(acos((1./RW_d2_k_tuned_all) .* log((exp(RW_d2_k_tuned_all)./2) + (exp(-RW_d2_k_tuned_all)./2))));


%%
for isess = RW_d3
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    RW_d3_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    RW_d3_ori_all = [RW_d3_ori_all RW_d3_ori];
    RW_d3_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    RW_d3_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    RW_d3_k_max_all = [RW_d3_k_max_all RW_d3_k_max];
    RW_d3_k = RW_d3_k_max.k1;
    RW_d3_k_all = [RW_d3_k_all RW_d3_k];
    RW_d3_k_tuned_index = RW_d3_ori.ind_theta90;
    RW_d3_k_tuned = RW_d3_k_max.k1(RW_d3_k_tuned_index);
    RW_d3_k_tuned_all = [RW_d3_k_tuned_all RW_d3_k_tuned];
    RW_d3_max = RW_d3_k_max.max_dfof;
    RW_d3_max_all = [RW_d3_max_all RW_d3_max];
    RW_d3_max_tuned_index = RW_d3_ori.ind_theta90;
    RW_d3_max_tuned = RW_d3_k_max.max_dfof(RW_d3_max_tuned_index);
    RW_d3_max_tuned_all = [RW_d3_max_tuned_all RW_d3_max_tuned];
    RW_d3_avgtuningcurve = RW_d3_newavg.newAvg(RW_d3_k_tuned_index,:)';
    RW_d3_avgtuningcurve_all = [RW_d3_avgtuningcurve_all RW_d3_avgtuningcurve];
end


RW_d3_avgtuningcurve_all = [RW_d3_avgtuningcurve_all]';

RW_d3_tuning_width_all = rad2deg(acos((1./RW_d3_k_tuned_all) .* log((exp(RW_d3_k_tuned_all)./2) + (exp(-RW_d3_k_tuned_all)./2))));



%%
%MAX DF/F
%all mice all (red and green cells)
figure; 
h=cdfplot(d1_max_tuned_all);
hold on
j=cdfplot(RW_d1_max_tuned_all);
% k=cdfplot(RW_d1_max_tuned_all);
% l=cdfplot(RW_d1_max_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
% set(k, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
% set(l, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
hold off
%legend(['Arc Red (n = ', num2str(length(d1_max_tunedall)), ')'], ['Arc Green (n = ', num2str(length(d1_max_tuned_all)), ')'], ...
%    ['RW Red (n = ', num2str(length(RW_d1_max_tunedall)), ')'], ['RW Green (n = ', num2str(length(RW_d1_max_tuned_all)), ')'],'Location', 'best')
% legend(['Day 1 v Day 2'; 'Day 1 v Day 3'])
legend(['Pre RW (n = ', num2str(length(d1_max_tuned_all)), ')'], ['Post RW (n = ', num2str(length(RW_d1_max_tuned_all)), ')'], 'Location', 'best')
xlabel('Max dF/F Values')
xlim([0 1]);
title('Day 1 Max dF/F Values');
print(fullfile(newfnout, ['tj all mice all cells', '_d1_maxdfof_cdf.pdf']), '-dpdf', '-bestfit')

%MAX DF/F
%all mice all (red and green cells)
figure; 
h=cdfplot(d2_max_tuned_all);
hold on
j=cdfplot(RW_d2_max_tuned_all);
% k=cdfplot(RW_d1_max_tuned_all);
% l=cdfplot(RW_d1_max_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
% set(k, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
% set(l, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
hold off
%legend(['Arc Red (n = ', num2str(length(d1_max_tunedall)), ')'], ['Arc Green (n = ', num2str(length(d1_max_tuned_all)), ')'], ...
%    ['RW Red (n = ', num2str(length(RW_d1_max_tunedall)), ')'], ['RW Green (n = ', num2str(length(RW_d1_max_tuned_all)), ')'],'Location', 'best')
% legend(['Day 1 v Day 2'; 'Day 1 v Day 3'])
legend(['Pre RW (n = ', num2str(length(d2_max_tuned_all)), ')'], ['Post RW (n = ', num2str(length(RW_d2_max_tuned_all)), ')'], 'Location', 'best')
xlabel('Max dF/F Values')
xlim([0 1]);
title('Day 2 Max dF/F Values');
print(fullfile(newfnout, ['tj all mice all cells', '_d2_maxdfof_cdf.pdf']), '-dpdf', '-bestfit')

%MAX DF/F
%all mice all (red and green cells)
figure; 
h=cdfplot(d3_max_tuned_all);
hold on
j=cdfplot(RW_d3_max_tuned_all);
% k=cdfplot(RW_d1_max_tuned_all);
% l=cdfplot(RW_d1_max_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
% set(k, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
% set(l, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
hold off
%legend(['Arc Red (n = ', num2str(length(d1_max_tunedall)), ')'], ['Arc Green (n = ', num2str(length(d1_max_tuned_all)), ')'], ...
%    ['RW Red (n = ', num2str(length(RW_d1_max_tunedall)), ')'], ['RW Green (n = ', num2str(length(RW_d1_max_tuned_all)), ')'],'Location', 'best')
% legend(['Day 1 v Day 2'; 'Day 1 v Day 3'])
legend(['Pre RW (n = ', num2str(length(d3_max_tuned_all)), ')'], ['Post RW (n = ', num2str(length(RW_d3_max_tuned_all)), ')'], 'Location', 'best')
xlabel('Max dF/F Values')
xlim([0 1]);
title('Day 3 Max dF/F Values');
print(fullfile(newfnout, ['tj all mice all cells', '_d3_maxdfof_cdf.pdf']), '-dpdf', '-bestfit')



%%
%K VALS
%all mice all (red and green) cells
figure; 
h=cdfplot(d1_k_tuned_all);
hold on
j=cdfplot(RW_d1_k_tuned_all);
% k=cdfplot(RW_d1_k_tuned_all);
% l=cdfplot(RW_d1_k_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
% set(k, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
% set(l, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
hold off
%legend(['Arc Red (n = ', num2str(length(d1_k_tunedall)), ')'], ['Arc Green (n = ', num2str(length(d1_k_tuned_all)), ')'], ...
%    ['RW Red (n = ', num2str(length(RW_d1_k_tunedall)), ')'], ['RW Green (n = ', num2str(length(RW_d1_k_tuned_all)), ')'],'Location', 'best')
% legend(['Day 1 v Day 2'; 'Day 1 v Day 3'])
legend('Pre RW', 'Post RW', 'Location', 'best')
xlabel('k Values')
xlim([0 30]);
title('Day 1 K Values');
print(fullfile(newfnout, ['tj all mice all cells', '_d1_k_cdf.pdf']), '-dpdf', '-bestfit')

%K VALS
%all mice all (red and green) cells
figure; 
h=cdfplot(d2_k_tuned_all);
hold on
j=cdfplot(RW_d2_k_tuned_all);
% k=cdfplot(RW_d2_k_tuned_all);
% l=cdfplot(RW_d2_k_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
% set(k, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
% set(l, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
hold off
%legend(['Arc Red (n = ', num2str(length(d2_k_tunedall)), ')'], ['Arc Green (n = ', num2str(length(d2_k_tuned_all)), ')'], ...
%    ['RW Red (n = ', num2str(length(RW_d2_k_tunedall)), ')'], ['RW Green (n = ', num2str(length(RW_d2_k_tuned_all)), ')'],'Location', 'best')
% legend(['Day 1 v Day 2'; 'Day 1 v Day 3'])
legend('Pre RW', 'Post RW', 'Location', 'best')
xlabel('k Values')
xlim([0 30]);
title('Day 2 K Values');
print(fullfile(newfnout, ['tj all mice all cells', '_d2_k_cdf.pdf']), '-dpdf', '-bestfit')

%K VALS
%all mice all (red and green) cells
figure; 
h=cdfplot(d3_k_tuned_all);
hold on
j=cdfplot(RW_d3_k_tuned_all);
% k=cdfplot(RW_d3_k_tuned_all);
% l=cdfplot(RW_d3_k_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
% set(k, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
% set(l, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
hold off
%legend(['Arc Red (n = ', num2str(length(d3_k_tunedall)), ')'], ['Arc Green (n = ', num2str(length(d3_k_tuned_all)), ')'], ...
%    ['RW Red (n = ', num2str(length(RW_d3_k_tunedall)), ')'], ['RW Green (n = ', num2str(length(RW_d3_k_tuned_all)), ')'],'Location', 'best')
% legend(['Day 1 v Day 2'; 'Day 1 v Day 3'])
legend('Pre RW', 'Post RW', 'Location', 'best')
xlabel('k Values')
xlim([0 30]);
title('Day 3 K Values');
print(fullfile(newfnout, ['tj all mice all cells', '_d3_k_cdf.pdf']), '-dpdf', '-bestfit')




%%
%tuning width 
%all mice all (red and green) cells
figure; 
h=cdfplot(d1_tuning_width_all);
hold on
j=cdfplot(RW_d1_tuning_width_all);
% k=cdfplot(RW_d1_k_tuned_all);
% l=cdfplot(RW_d1_k_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
% set(k, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
% set(l, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
hold off
%legend(['Arc Red (n = ', num2str(length(d1_k_tunedall)), ')'], ['Arc Green (n = ', num2str(length(d1_k_tuned_all)), ')'], ...
%    ['RW Red (n = ', num2str(length(RW_d1_k_tunedall)), ')'], ['RW Green (n = ', num2str(length(RW_d1_k_tuned_all)), ')'],'Location', 'best')
% legend(['Day 1 v Day 2'; 'Day 1 v Day 3'])
legend('Pre RW', 'Post RW', 'Location', 'best')
xlabel('Tuning Width Values')
%xlim([0 30]);
title('Day 1 Tuning Width Values');
print(fullfile(newfnout, ['tj all mice all cells', '_d1_tuning_width_cdf.pdf']), '-dpdf', '-bestfit')

%tuning width 
%all mice all (red and green) cells
figure; 
h=cdfplot(d2_tuning_width_all);
hold on
j=cdfplot(RW_d2_tuning_width_all);
% k=cdfplot(RW_d2_k_tuned_all);
% l=cdfplot(RW_d2_k_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
% set(k, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
% set(l, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
hold off
%legend(['Arc Red (n = ', num2str(length(d2_k_tunedall)), ')'], ['Arc Green (n = ', num2str(length(d2_k_tuned_all)), ')'], ...
%    ['RW Red (n = ', num2str(length(RW_d2_k_tunedall)), ')'], ['RW Green (n = ', num2str(length(RW_d2_k_tuned_all)), ')'],'Location', 'best')
% legend(['Day 1 v Day 2'; 'Day 1 v Day 3'])
legend('Pre RW', 'Post RW', 'Location', 'best')
xlabel('Tuning Width Values')
%xlim([0 30]);
title('Day 2 Tuning Width Values');
print(fullfile(newfnout, ['tj all mice all cells', '_d2_tuning_width_cdf.pdf']), '-dpdf', '-bestfit')

%tuning width 
%all mice all (red and green) cells
figure; 
h=cdfplot(d3_tuning_width_all);
hold on
j=cdfplot(RW_d3_tuning_width_all);
% k=cdfplot(RW_d3_k_tuned_all);
% l=cdfplot(RW_d3_k_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
set(j, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
% set(k, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
% set(l, 'LineStyle', '--', 'Color', 'g', 'LineWidth',2);
hold off
%legend(['Arc Red (n = ', num2str(length(d3_k_tunedall)), ')'], ['Arc Green (n = ', num2str(length(d3_k_tuned_all)), ')'], ...
%    ['RW Red (n = ', num2str(length(RW_d3_k_tunedall)), ')'], ['RW Green (n = ', num2str(length(RW_d3_k_tuned_all)), ')'],'Location', 'best')
% legend(['Day 1 v Day 2'; 'Day 1 v Day 3'])
legend('Pre RW', 'Post RW', 'Location', 'best')
xlabel('Tuning Width Values')
%xlim([0 30]);
title('Day 3 Tuning Width Values');
print(fullfile(newfnout, ['tj all mice all cells', '_d3_tuning_width_cdf.pdf']), '-dpdf', '-bestfit')




%%
%K-S tests
k_s_max_prevspostRW = kstest2(d1_max_tuned_all, RW_d1_max_tuned_all)
k_s_k_prevspostRW = kstest2(d1_k_tuned_all, RW_d1_k_tuned_all)
k_s_width_prevspostRW = kstest2(d1_tuning_width_all, RW_d1_tuning_width_all)


%%
%avg tuning curves
figure;
fast_errbar(1:8,d1_avgtuningcurve_all,1,'color',[0 1 0]);hold on;fast_errbar(1:8,RW_d1_avgtuningcurve_all,1,'color',[1 0 0]);
xlabel('Orientation')
ylabel('dF/F')
legend(['Pre RW'], ['Post RW'])
hold off;
print(fullfile(newfnout, ['tj all mice all cells', '_d1_avgtuning.pdf']), '-dpdf', '-bestfit')


figure;
fast_errbar(1:8,d2_avgtuningcurve_all,1,'color',[0 1 0]);hold on;fast_errbar(1:8,RW_d2_avgtuningcurve_all,1,'color',[1 0 0]);
xlabel('Orientation')
ylabel('dF/F')
legend(['Pre RW'], ['Post RW'])
hold off;
print(fullfile(newfnout, ['tj all mice all cells', '_d2_avgtuning.pdf']), '-dpdf', '-bestfit')


figure;
fast_errbar(1:8,d3_avgtuningcurve_all,1,'color',[0 1 0]);hold on;fast_errbar(1:8,RW_d3_avgtuningcurve_all,1,'color',[1 0 0]);
xlabel('Orientation')
ylabel('dF/F')
legend(['Pre RW'], ['Post RW'])
hold off;
print(fullfile(newfnout, ['tj all mice all cells', '_d3_avgtuning.pdf']), '-dpdf', '-bestfit')




%% 

prefori_d1_d2_match_all = [];
prefori_d2_match_all = [];
prefori_d1_d3_match_all = [];
prefori_d3_match_all = [];


k_d1_d2_match_all = [];
k_d2_match_all = [];
k_d1_d3_match_all = [];
k_d3_match_all = [];


max_d1_d2_match_all = [];
max_d2_match_all = [];
max_d1_d3_match_all = [];
max_d3_match_all = [];

pref_dscores_all = [];
pref_dscores_ses_all = [];

pref_d_d1_d2_all = [];
pref_d_d1_d3_all = [];


k_d_d1_d2_all = [];
k_d_d1_d3_all = [];


max_d_d1_d2_all = [];
max_d_d1_d3_all = [];



k_scores = [];
max_scores = [];


list = [40 43];
for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    prefori_ses_all = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes']));
    prefori_scores = [prefori_scores prefori_ses_all];
    pref_dscores_ses_all = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores']));
    pref_dscores_all = [pref_dscores_all pref_dscores_ses_all];
    k_ses_all = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_k_changes']));
    k_scores = [k_scores k_ses_all];
    max_ses_all = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_max_changes']));
    max_scores = [max_scores max_ses_all];
end

for idata = 1:length(prefori_scores)
    prefori_d1_d2_match = prefori_scores(idata).prefori_d1_d2_match_tune;
    prefori_d2_match = prefori_scores(idata).prefori_d2_match_tune;
    prefori_d1_d3_match = prefori_scores(idata).prefori_d1_d3_match_tune;
    prefori_d3_match = prefori_scores(idata).prefori_d3_match_tune;
    prefori_d1_d2_match_all = [prefori_d1_d2_match_all prefori_d1_d2_match];
    prefori_d2_match_all = [prefori_d2_match_all prefori_d2_match];
    prefori_d1_d3_match_all = [prefori_d1_d3_match_all prefori_d1_d3_match];
    prefori_d3_match_all = [prefori_d3_match_all prefori_d3_match];

    %these are not tuned
%     k_d1_d2_match = k_scores(idata).k_d1_d2_match;
%     k_d2_match = k_scores(idata).k_d2_match;
%     k_d1_d3_match = k_scores(idata).k_d1_d3_match;
%     k_d3_match = k_scores(idata).k_d3_match;
%     k_d1_d2_match_all = [k_d1_d2_match_all k_d1_d2_match];
%     k_d2_match_all = [k_d2_match_all k_d2_match];
%     k_d1_d3_match_all = [k_d1_d3_match_all k_d1_d3_match];
%     k_d3_match_all = [k_d3_match_all k_d3_match];


%     max_d1_d2_match = max_scores(idata).max_d1_d2_match;
%     max_d2_match = max_scores(idata).max_d2_match;
%     max_d1_d3_match = max_scores(idata).max_d1_d3_match;
%     max_d3_match = max_scores(idata).max_d3_match;
%     max_d1_d2_match_all = [max_d1_d2_match_all max_d1_d2_match];
%     max_d2_match_all = [max_d2_match_all max_d2_match];
%     max_d1_d3_match_all = [max_d1_d3_match_all max_d1_d3_match];
%     max_d3_match_all = [max_d3_match_all max_d3_match];


end

for idata = 1:length(pref_dscores_all)
    pref_d_d1_d2 = pref_dscores_all(idata).d_score_prefori_d1_d2;
    pref_d_d1_d2_all = [pref_d_d1_d2_all pref_d_d1_d2];
    pref_d_d1_d3 = pref_dscores_all(idata).d_score_prefori_d1_d3;
    pref_d_d1_d3_all = [pref_d_d1_d3_all pref_d_d1_d3];


%     k_d_d1_d2 = pref_dscores_all(idata).dscore_k_d1_d2_match;
%     k_d_d1_d2_all = [k_d_d1_d2_all k_d_d1_d2];
%     k_d_d1_d3 = pref_dscores_all(idata).dscore_k_d1_d3_match;
%     k_d_d1_d3_all = [k_d_d1_d3_all k_d_d1_d3];


    
%     max_d_d1_d2 = pref_dscores_all(idata).dscore_max_d1_d2_match;
%     max_d_d1_d2_all = [max_d_d1_d2_all max_d_d1_d2];
%     max_d_d1_d3 = pref_dscores_all(idata).dscore_max_d1_d3_match;
%     max_d_d1_d3_all = [max_d_d1_d3_all max_d_d1_d3];


    
end
%%
%READY TO START SCATTERS AND DSCORES
%all arc both days all cells
figure('Position', [400 20 650 700]);
sgtitle(['All Pre-RW Mice'], 'Interpreter', 'None');
subplot(2,1,1);
scatter(prefori_d1_d2_match_all, prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 0.8);
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Green cells (n = ', num2str(length(prefori_d1_d2_match_all)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(prefori_d1_d3_match_all, prefori_d3_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 0.8);
xlabel('Day 1 Pref Ori');
ylabel('Day 3 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
legend(['Green cells (n = ', num2str(length(prefori_d1_d3_match_all)), ')'], 'Location', 'northwest')
print(fullfile(newfnout, ['tj mice pre RW', '_prefori_scatter.pdf']), '-dpdf', '-bestfit')




%% 
%KEEP DOING HERE WHAT WAS DONE ABOVE NOW WITH POST RW MICE
RW_prefori_d1_d2_match_all = [];
RW_prefori_d2_match_all = [];
RW_prefori_d1_d3_match_all = [];
RW_prefori_d3_match_all = [];


RW_pref_dscores_all = [];
RW_pref_dscores_ses_all = [];

RW_pref_d_d1_d2_all = [];
RW_pref_d_d1_d3_all = [];



RW_list = [46 49];
for iexp = RW_list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    RW_prefori_ses_all = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_ori_changes']));
    RW_prefori_scores = [RW_prefori_scores RW_prefori_ses_all];
    RW_pref_dscores_ses_all = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_d_scores']));
    RW_pref_dscores_all = [RW_pref_dscores_all RW_pref_dscores_ses_all];


end

for idata = 1:length(RW_prefori_scores)
    RW_prefori_d1_d2_match = RW_prefori_scores(idata).prefori_d1_d2_match_tune;
    RW_prefori_d2_match = RW_prefori_scores(idata).prefori_d2_match_tune;
    RW_prefori_d1_d3_match = RW_prefori_scores(idata).prefori_d1_d3_match_tune;
    RW_prefori_d3_match = RW_prefori_scores(idata).prefori_d3_match_tune;
    RW_prefori_d1_d2_match_all = [RW_prefori_d1_d2_match_all RW_prefori_d1_d2_match];
    RW_prefori_d2_match_all = [RW_prefori_d2_match_all RW_prefori_d2_match];
    RW_prefori_d1_d3_match_all = [RW_prefori_d1_d3_match_all RW_prefori_d1_d3_match];
    RW_prefori_d3_match_all = [RW_prefori_d3_match_all RW_prefori_d3_match];
end

for idata = 1:length(RW_pref_dscores_all)
    RW_pref_d_d1_d2 = RW_pref_dscores_all(idata).d_score_prefori_d1_d2;
    RW_pref_d_d1_d2_all = [RW_pref_d_d1_d2_all RW_pref_d_d1_d2];
    RW_pref_d_d1_d3 = RW_pref_dscores_all(idata).d_score_prefori_d1_d3;
    RW_pref_d_d1_d3_all = [RW_pref_d_d1_d3_all RW_pref_d_d1_d3];
end


%%
%RW scatters
figure('Position', [400 20 650 700]);
sgtitle(['All Post RW Mice'], 'Interpreter', 'None');
subplot(2,1,1);
scatter(RW_prefori_d1_d2_match_all, RW_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 0.8);
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Green cells (n = ', num2str(length(RW_prefori_d1_d2_match_all)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(RW_prefori_d1_d3_match_all, RW_prefori_d3_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 0.8);
xlabel('Day 1 Pref Ori');
ylabel('Day 3 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
legend(['Green cells (n = ', num2str(length(RW_prefori_d1_d3_match_all)), ')'], 'Location', 'northwest')
print(fullfile(newfnout, ['tj post RW mice all cells', '_prefori_scatter.pdf']), '-dpdf', '-bestfit')


%%
%pre vs post RW 
figure('Position', [400 20 650 700]);
sgtitle(['Arc Promoter v RW Mice'], 'Interpreter', 'None');
subplot(2,1,1);
scatter(prefori_d1_d2_match_all, prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g', 'LineWidth', 0.8);
hold on
scatter(RW_prefori_d1_d2_match_all, RW_prefori_d2_match_all,'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Pre RW (n = ', num2str(length(prefori_d1_d2_match_all)), ')'], ['Post RW (n = ', num2str(length(RW_prefori_d1_d2_match_all)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(prefori_d1_d3_match_all, prefori_d3_match_all, 'filled', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g', 'LineWidth', 0.8);
hold on
scatter(RW_prefori_d1_d3_match_all, RW_prefori_d3_match_all,'r');
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
%legend(['Arc Promoter (n = ', num2str(length(red_prefori_d1_d3_match_all)), ')'], ['RW (n = ', num2str(length(RWprefori_d1_d3_match_all)), ')'], 'Location', 'northwest')
legend(['Pre RW (n = ', num2str(length(prefori_d1_d3_match_all)), ')'], ['Post RW (n = ', num2str(length(RW_prefori_d1_d3_match_all)), ')'], 'Location', 'northwest')
print(fullfile(newfnout, ['tj all mice pre vs post', '_prefori_scatter.pdf']), '-dpdf', '-bestfit')


%%
%dscores
figure; 
h = cdfplot(pref_d_d1_d2_all);
hold on
j = cdfplot(pref_d_d1_d3_all);
m = cdfplot(RW_pref_d_d1_d2_all);
n = cdfplot(RW_pref_d_d1_d3_all);
% p = cdfplot(enh_dscore_prefori_d1_d2_all);
% q = cdfplot(enh_dscore_prefori_d1_d3_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth', 2.0);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth', 2.0);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 2.0);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 2.0);
% set(p, 'LineStyle', '-', 'Color', 'g');
% set(q, 'LineStyle', '--', 'Color', 'g');
%legend('Arc D1-D2', 'Arc D1-D3', 'RW D1-D2', 'RW D1-D3', 'Arc-Enh D1-D2', 'Arc-Enh D1-D3', 'Location', 'Best')
legend(['Pre RW D1-D2 (n = ', num2str(length(pref_d_d1_d2_all)), ')'],...
    ['Pre RW D1-D3 (n = ', num2str(length(pref_d_d1_d3_all)), ')'],...
    ['Post RW D1-D2 (n = ', num2str(length(RW_pref_d_d1_d2_all)), ')'],...
    ['Post RW D1-D3 (n = ', num2str(length(RW_pref_d_d1_d3_all)), ')'], 'Location', 'Best')
%title('Change in Pref Ori Values')
title('')
xlabel('Change in Pref Ori')
xlim([0 90])
hold off
print(fullfile(newfnout, ['tj pre vs post', '_change_pref_cdf.pdf']), '-dpdf', '-bestfit')




%%
concat_dscores_RW = [RW_pref_d_d1_d2_all RW_pref_d_d1_d3_all];
concat_dscores = [pref_d_d1_d2_all pref_d_d1_d3_all];

figure;
y = cdfplot(concat_dscores)
hold on
z = cdfplot(concat_dscores_RW)
%legend('Arc Promoter KRAB Cells', 'RW KRAB Cells', 'Location', 'best')
legend(['Pre RW (n = ', num2str(length(concat_dscores)), ')'],...
    ['Post RW (n = ', num2str(length(concat_dscores_RW)), ')'],...
    'Location', 'Best')
set(y, 'LineStyle', '-', 'Color', 'g', 'LineWidth', 2.0);
set(z, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 2.0);
title('')
xlabel('Change in Pref Ori')
ylabel('Proportion of Cells')
xlim([0 90])
hold off
print(fullfile(newfnout, ['tj concat dscores pre vs post RW', '_change_pref_cdf.pdf']), '-dpdf', '-bestfit')


%%
%prefori dscores k-s tests

k_s_prefd_preVSpostRW = kstest2(concat_dscores, concat_dscores_RW);
[ks, ks1, ks2] = kstest2(concat_dscores, concat_dscores_RW)

%%
%STOPPED HERE FOR NOW
%
%%

RW_k_d1_d2_match_all = [];
RW_k_d2_match_all = [];
RW_k_d1_d3_match_all = [];
RW_k_d3_match_all = [];
RWk_d1_d2_match_all = [];
RWk_d2_match_all = [];
RWk_d1_d3_match_all = [];
RWk_d3_match_all = [];

RW_max_d1_d2_match_all = [];
RW_max_d2_match_all = [];
RW_max_d1_d3_match_all = [];
RW_max_d3_match_all = [];
RWmax_d1_d2_match_all = [];
RWmax_d2_match_all = [];
RWmax_d1_d3_match_all = [];
RWmax_d3_match_all = [];

RW_k_d_d1_d2_all = [];
RW_k_d_d1_d3_all = [];
RWk_d_d1_d2_all = [];
RWk_d_d1_d3_all = [];

RW_max_d_d1_d2_all = [];
RW_max_d_d1_d3_all = [];
RWmax_d_d1_d2_all = [];
RWmax_d_d1_d3_all = [];


RW_k_scores = [];
RW_max_scores = [];


RW_list = [1 4 37 40];
for iexp = RW_list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    RW_k_ses1 = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_k_changes']));
    RW_k_ses2 = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_k_changes_2']));
    RW_k_ses_all = [RW_k_ses1 RW_k_ses2];
    RW_k_scores = [RW_k_scores RW_k_ses_all];
    RW_max_ses1 = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_max_changes']));
    RW_max_ses2 = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'RW_max_changes_2']));
    RW_max_ses_all = [RW_max_ses1 RW_max_ses2];
    RW_max_scores = [RW_max_scores RW_max_ses_all];
end

for idata = 1:length(RW_prefori_scores)
    
    RW_k_d1_d2_match = RW_k_scores(idata).k_d1_d2_match;
    RW_k_d2_match = RW_k_scores(idata).k_d2_match;
    RW_k_d1_d3_match = RW_k_scores(idata).k_d1_d3_match;
    RW_k_d3_match = RW_k_scores(idata).k_d3_match;
    RW_k_d1_d2_match_all = [RW_k_d1_d2_match_all RW_k_d1_d2_match];
    RW_k_d2_match_all = [RW_k_d2_match_all RW_k_d2_match];
    RW_k_d1_d3_match_all = [RW_k_d1_d3_match_all RW_k_d1_d3_match];
    RW_k_d3_match_all = [RW_k_d3_match_all RW_k_d3_match];
    RWk_d1_d2_match = RW_k_scores(idata).red_k_d1_d2_match;
    RWk_d2_match = RW_k_scores(idata).red_k_d2_match;
    RWk_d1_d3_match = RW_k_scores(idata).red_k_d1_d3_match;
    RWk_d3_match = RW_k_scores(idata).red_k_d3_match;
    RWk_d1_d2_match_all = [RWk_d1_d2_match_all RWk_d1_d2_match];
    RWk_d2_match_all = [RWk_d2_match_all RWk_d2_match];
    RWk_d1_d3_match_all = [RWk_d1_d3_match_all RWk_d1_d3_match];
    RWk_d3_match_all = [RWk_d3_match_all RWk_d3_match];

    RW_max_d1_d2_match = RW_max_scores(idata).max_d1_d2_match;
    RW_max_d2_match = RW_max_scores(idata).max_d2_match;
    RW_max_d1_d3_match = RW_max_scores(idata).max_d1_d3_match;
    RW_max_d3_match = RW_max_scores(idata).max_d3_match;
    RW_max_d1_d2_match_all = [RW_max_d1_d2_match_all RW_max_d1_d2_match];
    RW_max_d2_match_all = [RW_max_d2_match_all RW_max_d2_match];
    RW_max_d1_d3_match_all = [RW_max_d1_d3_match_all RW_max_d1_d3_match];
    RW_max_d3_match_all = [RW_max_d3_match_all RW_max_d3_match];
    RWmax_d1_d2_match = RW_max_scores(idata).red_max_d1_d2_match;
    RWmax_d2_match = RW_max_scores(idata).red_max_d2_match;
    RWmax_d1_d3_match = RW_max_scores(idata).red_max_d1_d3_match;
    RWmax_d3_match = RW_max_scores(idata).red_max_d3_match;
    RWmax_d1_d2_match_all = [RWmax_d1_d2_match_all RWmax_d1_d2_match];
    RWmax_d2_match_all = [RWmax_d2_match_all RWmax_d2_match];
    RWmax_d1_d3_match_all = [RWmax_d1_d3_match_all RWmax_d1_d3_match];
    RWmax_d3_match_all = [RWmax_d3_match_all RWmax_d3_match];

end

for idata = 1:length(RW_pref_dscores_all)
   
    RW_k_d_d1_d2 = RW_pref_dscores_all(idata).dscore_k_d1_d2_match;
    RW_k_d_d1_d2_all = [RW_k_d_d1_d2_all RW_k_d_d1_d2];
    RW_k_d_d1_d3 = RW_pref_dscores_all(idata).dscore_k_d1_d3_match;
    RW_k_d_d1_d3_all = [RW_k_d_d1_d3_all RW_k_d_d1_d3];
    RWk_d_d1_d2 = RW_pref_dscores_all(idata).red_dscore_k_d1_d2_match;
    RWk_d_d1_d2_all = [RWk_d_d1_d2_all RWk_d_d1_d2];
    RWk_d_d1_d3 = RW_pref_dscores_all(idata).red_dscore_k_d1_d3_match;
    RWk_d_d1_d3_all = [RWk_d_d1_d3_all RWk_d_d1_d3];

    
    RW_max_d_d1_d2 = RW_pref_dscores_all(idata).dscore_max_d1_d2_match;
    RW_max_d_d1_d2_all = [RW_max_d_d1_d2_all RW_max_d_d1_d2];
    RW_max_d_d1_d3 = RW_pref_dscores_all(idata).dscore_max_d1_d3_match;
    RW_max_d_d1_d3_all = [RW_max_d_d1_d3_all RW_max_d_d1_d3];
    RWmax_d_d1_d2 = RW_pref_dscores_all(idata).red_dscore_max_d1_d2_match;
    RWmax_d_d1_d2_all = [RWmax_d_d1_d2_all RWmax_d_d1_d2];
    RWmax_d_d1_d3 = RW_pref_dscores_all(idata).red_dscore_max_d1_d3_match;
    RWmax_d_d1_d3_all = [RWmax_d_d1_d3_all RWmax_d_d1_d3];

    
end



%%
%arc vs RW red k scatters
figure('Position', [400 20 650 700]);
sgtitle(['Arc Promoter v RW Mice (Red Cells)'], 'Interpreter', 'None');
subplot(2,1,1);
scatter(red_k_d1_d2_match_all, red_k_d2_match_all, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 0.8);
hold on
scatter(RWk_d1_d2_match_all, RWk_d2_match_all,'r');
hold off
xlabel('Day 1 k');
ylabel('Day 2 k');
xlim([0,30]);
ylim([0,30]);
%xticks(0:20:180);
%yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Arc Promoter (n = ', num2str(length(red_k_d1_d2_match_all)), ')'], ['RW (n = ', num2str(length(RWk_d1_d2_match_all)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(red_k_d1_d3_match_all, red_k_d3_match_all, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 0.8);
hold on
scatter(RWk_d1_d3_match_all, RWk_d3_match_all,'r');
hold off
xlabel('Day 1 k');
ylabel('Day 3 k');
xlim([0,30]);
ylim([0,30]);
%xticks(0:20:180);
%yticks(0:20:180);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
legend(['Arc Promoter (n = ', num2str(length(red_k_d1_d3_match_all)), ')'], ['RW (n = ', num2str(length(RWk_d1_d3_match_all)), ')'], 'Location', 'northwest')
print(fullfile(newfnout, ['tj all mice red cells', '_k_scatter.pdf']), '-dpdf', '-bestfit')

%arc vs RW green k scatters
figure('Position', [400 20 650 700]);
sgtitle(['Arc Promoter v RW Mice (Green Cells)'], 'Interpreter', 'None');
subplot(2,1,1);
scatter(k_d1_d2_match_all, k_d2_match_all, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 0.8);
hold on
scatter(RW_k_d1_d2_match_all, RW_k_d2_match_all,'r');
hold off
xlabel('Day 1 k');
ylabel('Day 2 k');
xlim([0,30]);
ylim([0,30]);
%xticks(0:20:180);
%yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Arc Promoter (n = ', num2str(length(k_d1_d2_match_all)), ')'], ['RW (n = ', num2str(length(RW_k_d1_d2_match_all)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(k_d1_d3_match_all, k_d3_match_all, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 0.8);
hold on
scatter(RW_k_d1_d3_match_all, RW_k_d3_match_all,'r');
hold off
xlabel('Day 1 k');
ylabel('Day 3 k');
xlim([0,30]);
ylim([0,30]);
%xticks(0:20:180);
%yticks(0:20:180);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
legend(['Arc Promoter (n = ', num2str(length(k_d1_d3_match_all)), ')'], ['RW (n = ', num2str(length(RW_k_d1_d3_match_all)), ')'], 'Location', 'northwest')
print(fullfile(newfnout, ['tj all mice green cells', '_k_scatter.pdf']), '-dpdf', '-bestfit')


%%
%arc vs RW red max scatters
figure('Position', [400 20 650 700]);
sgtitle(['Arc Promoter v RW Mice (Red Cells)'], 'Interpreter', 'None');
subplot(2,1,1);
scatter(red_max_d1_d2_match_all, red_max_d2_match_all, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 0.8);
hold on
scatter(RWmax_d1_d2_match_all, RWmax_d2_match_all,'r');
hold off
xlabel('Day 1 Max');
ylabel('Day 2 Max');
%xlim([0,30]);
%ylim([0,30]);
%xticks(0:20:180);
%yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Arc Promoter (n = ', num2str(length(red_max_d1_d2_match_all)), ')'], ['RW (n = ', num2str(length(RWmax_d1_d2_match_all)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(red_max_d1_d3_match_all, red_max_d3_match_all, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 0.8);
hold on
scatter(RWmax_d1_d3_match_all, RWmax_d3_match_all,'r');
hold off
xlabel('Day 1 Max');
ylabel('Day 3 Max');
%xlim([0,30]);
%ylim([0,30]);
%xticks(0:20:180);
%yticks(0:20:180);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
legend(['Arc Promoter (n = ', num2str(length(red_max_d1_d3_match_all)), ')'], ['RW (n = ', num2str(length(RWmax_d1_d3_match_all)), ')'], 'Location', 'northwest')
print(fullfile(newfnout, ['tj all mice red cells', '_max_scatter.pdf']), '-dpdf', '-bestfit')

%arc vs RW green max scatters
figure('Position', [400 20 650 700]);
sgtitle(['Arc Promoter v RW Mice (Green Cells)'], 'Interpreter', 'None');
subplot(2,1,1);
scatter(max_d1_d2_match_all, max_d2_match_all, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 0.8);
hold on
scatter(RW_max_d1_d2_match_all, RW_max_d2_match_all,'r');
hold off
xlabel('Day 1 Max');
ylabel('Day 2 Max');
%xlim([0,30]);
%ylim([0,30]);
%xticks(0:20:180);
%yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Arc Promoter (n = ', num2str(length(max_d1_d2_match_all)), ')'], ['RW (n = ', num2str(length(RW_max_d1_d2_match_all)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(max_d1_d3_match_all, max_d3_match_all, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 0.8);
hold on
scatter(RW_max_d1_d3_match_all, RW_max_d3_match_all,'r');
hold off
xlabel('Day 1 Max');
ylabel('Day 3 Max');
%xlim([0,30]);
%ylim([0,30]);
%xticks(0:20:180);
%yticks(0:20:180);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
legend(['Arc Promoter (n = ', num2str(length(max_d1_d3_match_all)), ')'], ['RW (n = ', num2str(length(RW_max_d1_d3_match_all)), ')'], 'Location', 'northwest')
print(fullfile(newfnout, ['tj all mice green cells', '_max_scatter.pdf']), '-dpdf', '-bestfit')


%%
%DO THE MAX AND K D SCORE CDFS
%DO THIS ALL FOR ENHANCER TOO
%TRY THE SPATIAL PREF ORI EXPLORATION
%WORK ON INDIVIDUAL MICE


%%
%dscores for k
figure; 
h = cdfplot(k_d_d1_d2_all);
hold on
j = cdfplot(k_d_d1_d3_all);
m = cdfplot(red_k_d_d1_d2_all);
n = cdfplot(red_k_d_d1_d3_all);
% p = cdfplot(enh_dscore_kori_d1_d2_all);
% q = cdfplot(enh_dscore_kori_d1_d3_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth', 1.0);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth', 1.0);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 1.0);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.0);
% set(p, 'LineStyle', '-', 'Color', 'g');
% set(q, 'LineStyle', '--', 'Color', 'g');
%legend('Arc D1-D2', 'Arc D1-D3', 'RW D1-D2', 'RW D1-D3', 'Arc-Enh D1-D2', 'Arc-Enh D1-D3', 'Location', 'Best')
legend(['Arc Green Cell D1-D2 (n = ', num2str(length(k_d_d1_d2_all)), ')'],...
    ['Arc Green Cell D1-D3 (n = ', num2str(length(k_d_d1_d3_all)), ')'],...
    ['Arc Red Cell D1-D2 (n = ', num2str(length(red_k_d_d1_d2_all)), ')'],...
    ['Arc Red Cell D1-D3 (n = ', num2str(length(red_k_d_d1_d3_all)), ')'], 'Location', 'Best')
%title('Change in K Ori Values')
title('')
xlabel('Change in k')
xlim([0 30])
hold off
print(fullfile(newfnout, ['tj arc mice all cells', '_change_k_cdf.pdf']), '-dpdf', '-bestfit')

%%
figure; 
h = cdfplot(RW_k_d_d1_d2_all);
hold on
j = cdfplot(RW_k_d_d1_d3_all);
m = cdfplot(RWk_d_d1_d2_all);
n = cdfplot(RWk_d_d1_d3_all);
% p = cdfplot(enh_dscore_kori_d1_d2_all);
% q = cdfplot(enh_dscore_kori_d1_d3_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth', 1.0);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth', 1.0);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 1.0);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.0);
% set(p, 'LineStyle', '-', 'Color', 'g');
% set(q, 'LineStyle', '--', 'Color', 'g');
%legend('RW D1-D2', 'RW D1-D3', 'RW D1-D2', 'RW D1-D3', 'RW-Enh D1-D2', 'RW-Enh D1-D3', 'Location', 'Best')
legend(['RW Green Cell D1-D2 (n = ', num2str(length(RW_k_d_d1_d2_all)), ')'],...
    ['RW Green Cell D1-D3 (n = ', num2str(length(RW_k_d_d1_d3_all)), ')'],...
    ['RW Red Cell D1-D2 (n = ', num2str(length(RWk_d_d1_d2_all)), ')'],...
    ['RW Red Cell D1-D3 (n = ', num2str(length(RWk_d_d1_d3_all)), ')'], 'Location', 'Best')
%title('Change in K Ori Values')
title('')
xlabel('Change in k')
xlim([0 30])
hold off
print(fullfile(newfnout, ['tj RW mice all cells', '_change_k_cdf.pdf']), '-dpdf', '-bestfit')

%%
%arc v RW red
figure; 
h = cdfplot(red_k_d_d1_d2_all);
hold on
j = cdfplot(red_k_d_d1_d3_all);
m = cdfplot(RWk_d_d1_d2_all);
n = cdfplot(RWk_d_d1_d3_all);
% p = cdfplot(enh_dscore_kori_d1_d2_all);
% q = cdfplot(enh_dscore_kori_d1_d3_all);
set(h, 'LineStyle', '-', 'Color', 'b', 'LineWidth', 1.0);
set(j, 'LineStyle', '--', 'Color', 'b', 'LineWidth', 1.0);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 1.0);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.0);
% set(p, 'LineStyle', '-', 'Color', 'g');
% set(q, 'LineStyle', '--', 'Color', 'g');
%legend('Arc D1-D2', 'Arc D1-D3', 'RW D1-D2', 'RW D1-D3', 'Arc-Enh D1-D2', 'Arc-Enh D1-D3', 'Location', 'Best')
legend(['Arc Prom D1-D2 (n = ', num2str(length(red_k_d_d1_d2_all)), ')'],...
    ['Arc Prom D1-D3 (n = ', num2str(length(red_k_d_d1_d3_all)), ')'],...
    ['RW D1-D2 (n = ', num2str(length(RWk_d_d1_d2_all)), ')'],...
    ['RW D1-D3 (n = ', num2str(length(RWk_d_d1_d3_all)), ')'], 'Location', 'Best')
%title('Change in K Ori Values')
title('')
xlabel('Change in K Ori')
xlim([0 30])
hold off
print(fullfile(newfnout, ['tj all mice red cells', '_change_k_cdf.pdf']), '-dpdf', '-bestfit')

%started exploring percentiles to quantify differences
[prctile(red_k_d_d1_d2_all, 90), prctile(red_k_d_d1_d3_all, 90)];
[prctile(RWk_d_d1_d2_all, 90), prctile(RWk_d_d1_d3_all, 90)];
%%
%arc v RW green
figure; 
h = cdfplot(k_d_d1_d2_all);
hold on
j = cdfplot(k_d_d1_d3_all);
m = cdfplot(RW_k_d_d1_d2_all);
n = cdfplot(RW_k_d_d1_d3_all);
% p = cdfplot(enh_dscore_kori_d1_d2_all);
% q = cdfplot(enh_dscore_kori_d1_d3_all);
set(h, 'LineStyle', '-', 'Color', 'b', 'LineWidth', 1.0);
set(j, 'LineStyle', '--', 'Color', 'b', 'LineWidth', 1.0);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 1.0);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.0);
% set(p, 'LineStyle', '-', 'Color', 'g');
% set(q, 'LineStyle', '--', 'Color', 'g');
%legend('Arc D1-D2', 'Arc D1-D3', 'RW D1-D2', 'RW D1-D3', 'Arc-Enh D1-D2', 'Arc-Enh D1-D3', 'Location', 'Best')
legend(['Arc Prom D1-D2 (n = ', num2str(length(k_d_d1_d2_all)), ')'],...
    ['Arc Prom D1-D3 (n = ', num2str(length(k_d_d1_d3_all)), ')'],...
    ['RW D1-D2 (n = ', num2str(length(RW_k_d_d1_d2_all)), ')'],...
    ['RW D1-D3 (n = ', num2str(length(RW_k_d_d1_d3_all)), ')'], 'Location', 'Best')
%title('Change in K Ori Values')
title('')
xlabel('Change in K Ori')
xlim([0 30])
hold off
print(fullfile(newfnout, ['tj all mice green cells', '_change_k_cdf.pdf']), '-dpdf', '-bestfit')


%%
%dscores for max
figure; 
h = cdfplot(max_d_d1_d2_all);
hold on
j = cdfplot(max_d_d1_d3_all);
m = cdfplot(red_max_d_d1_d2_all);
n = cdfplot(red_max_d_d1_d3_all);
% p = cdfplot(enh_dscore_kori_d1_d2_all);
% q = cdfplot(enh_dscore_kori_d1_d3_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth', 1.0);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth', 1.0);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 1.0);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.0);
% set(p, 'LineStyle', '-', 'Color', 'g');
% set(q, 'LineStyle', '--', 'Color', 'g');
%legend('Arc D1-D2', 'Arc D1-D3', 'RW D1-D2', 'RW D1-D3', 'Arc-Enh D1-D2', 'Arc-Enh D1-D3', 'Location', 'Best')
legend(['Arc Green Cell D1-D2 (n = ', num2str(length(max_d_d1_d2_all)), ')'],...
    ['Arc Green Cell D1-D3 (n = ', num2str(length(max_d_d1_d3_all)), ')'],...
    ['Arc Red Cell D1-D2 (n = ', num2str(length(red_max_d_d1_d2_all)), ')'],...
    ['Arc Red Cell D1-D3 (n = ', num2str(length(red_max_d_d1_d3_all)), ')'], 'Location', 'Best')
%title('Change in K Ori Values')
title('')
xlabel('Change in Max dF/F')
%xlim([0 30])
hold off
print(fullfile(newfnout, ['tj arc mice all cells', '_change_max_cdf.pdf']), '-dpdf', '-bestfit')

%%
figure; 
h = cdfplot(RW_max_d_d1_d2_all);
hold on
j = cdfplot(RW_max_d_d1_d3_all);
m = cdfplot(RWmax_d_d1_d2_all);
n = cdfplot(RWmax_d_d1_d3_all);
% p = cdfplot(enh_dscore_kori_d1_d2_all);
% q = cdfplot(enh_dscore_kori_d1_d3_all);
set(h, 'LineStyle', '-', 'Color', 'g', 'LineWidth', 1.0);
set(j, 'LineStyle', '--', 'Color', 'g', 'LineWidth', 1.0);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 1.0);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.0);
% set(p, 'LineStyle', '-', 'Color', 'g');
% set(q, 'LineStyle', '--', 'Color', 'g');
%legend('RW D1-D2', 'RW D1-D3', 'RW D1-D2', 'RW D1-D3', 'RW-Enh D1-D2', 'RW-Enh D1-D3', 'Location', 'Best')
legend(['RW Green Cell D1-D2 (n = ', num2str(length(RW_max_d_d1_d2_all)), ')'],...
    ['RW Green Cell D1-D3 (n = ', num2str(length(RW_max_d_d1_d3_all)), ')'],...
    ['RW Red Cell D1-D2 (n = ', num2str(length(RWmax_d_d1_d2_all)), ')'],...
    ['RW Red Cell D1-D3 (n = ', num2str(length(RWmax_d_d1_d3_all)), ')'], 'Location', 'Best')
%title('Change in K Ori Values')
title('')
xlabel('Change in Max dF/F')
%xlim([0 30])
hold off
print(fullfile(newfnout, ['tj RW mice all cells', '_change_max_cdf.pdf']), '-dpdf', '-bestfit')

%%
%arc v RW red
figure; 
h = cdfplot(red_max_d_d1_d2_all);
hold on
j = cdfplot(red_max_d_d1_d3_all);
m = cdfplot(RWmax_d_d1_d2_all);
n = cdfplot(RWmax_d_d1_d3_all);
% p = cdfplot(enh_dscore_kori_d1_d2_all);
% q = cdfplot(enh_dscore_kori_d1_d3_all);
set(h, 'LineStyle', '-', 'Color', 'b', 'LineWidth', 1.0);
set(j, 'LineStyle', '--', 'Color', 'b', 'LineWidth', 1.0);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 1.0);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.0);
% set(p, 'LineStyle', '-', 'Color', 'g');
% set(q, 'LineStyle', '--', 'Color', 'g');
%legend('Arc D1-D2', 'Arc D1-D3', 'RW D1-D2', 'RW D1-D3', 'Arc-Enh D1-D2', 'Arc-Enh D1-D3', 'Location', 'Best')
legend(['Arc Prom D1-D2 (n = ', num2str(length(red_max_d_d1_d2_all)), ')'],...
    ['Arc Prom D1-D3 (n = ', num2str(length(red_max_d_d1_d3_all)), ')'],...
    ['RW D1-D2 (n = ', num2str(length(RWmax_d_d1_d2_all)), ')'],...
    ['RW D1-D3 (n = ', num2str(length(RWmax_d_d1_d3_all)), ')'], 'Location', 'Best')
%title('Change in K Ori Values')
title('')
xlabel('Change in Max dF/F')
%xlim([0 30])
hold off
print(fullfile(newfnout, ['tj all mice red cells', '_change_max_cdf.pdf']), '-dpdf', '-bestfit')

%%
%arc v RW green
figure; 
h = cdfplot(max_d_d1_d2_all);
hold on
j = cdfplot(max_d_d1_d3_all);
m = cdfplot(RW_max_d_d1_d2_all);
n = cdfplot(RW_max_d_d1_d3_all);
% p = cdfplot(enh_dscore_kori_d1_d2_all);
% q = cdfplot(enh_dscore_kori_d1_d3_all);
set(h, 'LineStyle', '-', 'Color', 'b', 'LineWidth', 1.0);
set(j, 'LineStyle', '--', 'Color', 'b', 'LineWidth', 1.0);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 1.0);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.0);
% set(p, 'LineStyle', '-', 'Color', 'g');
% set(q, 'LineStyle', '--', 'Color', 'g');
%legend('Arc D1-D2', 'Arc D1-D3', 'RW D1-D2', 'RW D1-D3', 'Arc-Enh D1-D2', 'Arc-Enh D1-D3', 'Location', 'Best')
legend(['Arc Prom D1-D2 (n = ', num2str(length(max_d_d1_d2_all)), ')'],...
    ['Arc Prom D1-D3 (n = ', num2str(length(max_d_d1_d3_all)), ')'],...
    ['RW D1-D2 (n = ', num2str(length(RW_max_d_d1_d2_all)), ')'],...
    ['RW D1-D3 (n = ', num2str(length(RW_max_d_d1_d3_all)), ')'], 'Location', 'Best')
%title('Change in K Ori Values')
title('')
xlabel('Change in Max dF/F')
%xlim([0 30])
hold off
print(fullfile(newfnout, ['tj all mice green cells', '_change_max_cdf.pdf']), '-dpdf', '-bestfit')


%%
%try to optimize individual mice


%%
%individual mice d1 v d2 arc and RW
mouse1_d1_ktuned_index_ses1 = intersect(d1_matches_all(1).redCells, d1_ori_all(1).ind_theta90);
mouse1_k_d1_ses1 = d1_k_max_all(1).k1(mouse1_d1_ktuned_index_ses1);
mouse1_max_d1_ses1 = d1_k_max_all(1).max_dfof(mouse1_d1_ktuned_index_ses1);
mouse1_d1_ktuned_index_ses2 = intersect(d1_matches_all(2).redCells, d1_ori_all(2).ind_theta90);
mouse1_k_d1_ses2 = d1_k_max_all(2).k1(mouse1_d1_ktuned_index_ses2);
mouse1_max_d1_ses2 = d1_k_max_all(2).max_dfof(mouse1_d1_ktuned_index_ses2);
mouse1_k_d1_all = [mouse1_k_d1_ses1 mouse1_k_d1_ses2];
mouse1_max_d1_all = [mouse1_max_d1_ses1 mouse1_max_d1_ses2];

mouse2_d1_ktuned_index_ses1 = intersect(d1_matches_all(3).redCells, d1_ori_all(3).ind_theta90);
mouse2_k_d1_ses1 = d1_k_max_all(3).k1(mouse2_d1_ktuned_index_ses1);
mouse2_max_d1_ses1 = d1_k_max_all(3).max_dfof(mouse2_d1_ktuned_index_ses1);
mouse2_d1_ktuned_index_ses2 = intersect(d1_matches_all(4).redCells, d1_ori_all(4).ind_theta90);
mouse2_k_d1_ses2 = d1_k_max_all(4).k1(mouse2_d1_ktuned_index_ses2);
mouse2_max_d1_ses2 = d1_k_max_all(4).max_dfof(mouse2_d1_ktuned_index_ses2);
mouse2_k_d1_all = [mouse2_k_d1_ses1 mouse2_k_d1_ses2];
mouse2_max_d1_all = [mouse2_max_d1_ses1 mouse2_max_d1_ses2];

mouse3_d1_ktuned_index_ses1 = intersect(d1_matches_all(5).redCells, d1_ori_all(5).ind_theta90);
mouse3_k_d1_ses1 = d1_k_max_all(5).k1(mouse3_d1_ktuned_index_ses1);
mouse3_max_d1_ses1 = d1_k_max_all(5).max_dfof(mouse3_d1_ktuned_index_ses1);
mouse3_d1_ktuned_index_ses2 = intersect(d1_matches_all(6).redCells, d1_ori_all(6).ind_theta90);
mouse3_max_d1_ses2 = d1_k_max_all(6).max_dfof(mouse3_d1_ktuned_index_ses2);
mouse3_k_d1_ses2 = d1_k_max_all(6).k1(mouse3_d1_ktuned_index_ses2);
mouse3_k_d1_all = [mouse3_k_d1_ses1 mouse3_k_d1_ses2];
mouse3_max_d1_all = [mouse3_max_d1_ses1 mouse3_max_d1_ses2];

mouse4_d1_ktuned_index_ses1 = intersect(d1_matches_all(7).redCells, d1_ori_all(7).ind_theta90);
mouse4_k_d1_ses1 = d1_k_max_all(7).k1(mouse4_d1_ktuned_index_ses1);
mouse4_max_d1_ses1 = d1_k_max_all(7).max_dfof(mouse4_d1_ktuned_index_ses1);
mouse4_d1_ktuned_index_ses2 = intersect(d1_matches_all(8).redCells, d1_ori_all(8).ind_theta90);
mouse4_k_d1_ses2 = d1_k_max_all(8).k1(mouse4_d1_ktuned_index_ses2);
mouse4_max_d1_ses2 = d1_k_max_all(8).max_dfof(mouse4_d1_ktuned_index_ses2);
mouse4_k_d1_all = [mouse4_k_d1_ses1 mouse4_k_d1_ses2];
mouse4_max_d1_all = [mouse4_max_d1_ses1 mouse4_max_d1_ses2];

mouse1_RW_d1_ktuned_index_ses1 = intersect(RW_d1_matches_all(1).redCells, RW_d1_ori_all(1).ind_theta90);
mouse1_RW_k_d1_ses1 = RW_d1_k_max_all(1).k1(mouse1_RW_d1_ktuned_index_ses1);
mouse1_RW_max_d1_ses1 = RW_d1_k_max_all(1).max_dfof(mouse1_RW_d1_ktuned_index_ses1);
mouse1_RW_d1_ktuned_index_ses2 = intersect(RW_d1_matches_all(2).redCells, RW_d1_ori_all(2).ind_theta90);
mouse1_RW_k_d1_ses2 = RW_d1_k_max_all(2).k1(mouse1_RW_d1_ktuned_index_ses2);
mouse1_RW_max_d1_ses2 = RW_d1_k_max_all(2).max_dfof(mouse1_RW_d1_ktuned_index_ses2);
mouse1_RW_k_d1_all = [mouse1_RW_k_d1_ses1 mouse1_RW_k_d1_ses2];
mouse1_RW_max_d1_all = [mouse1_RW_max_d1_ses1 mouse1_RW_max_d1_ses2];

mouse2_RW_d1_ktuned_index_ses1 = intersect(RW_d1_matches_all(3).redCells, RW_d1_ori_all(3).ind_theta90);
mouse2_RW_k_d1_ses1 = RW_d1_k_max_all(3).k1(mouse2_RW_d1_ktuned_index_ses1);
mouse2_RW_max_d1_ses1 = RW_d1_k_max_all(3).max_dfof(mouse2_RW_d1_ktuned_index_ses1);
mouse2_RW_d1_ktuned_index_ses2 = intersect(RW_d1_matches_all(4).redCells, RW_d1_ori_all(4).ind_theta90);
mouse2_RW_k_d1_ses2 = RW_d1_k_max_all(4).k1(mouse2_RW_d1_ktuned_index_ses2);
mouse2_RW_max_d1_ses2 = RW_d1_k_max_all(4).max_dfof(mouse2_RW_d1_ktuned_index_ses2);
mouse2_RW_k_d1_all = [mouse2_RW_k_d1_ses1 mouse2_RW_k_d1_ses2];
mouse2_RW_max_d1_all = [mouse2_RW_max_d1_ses1 mouse2_RW_max_d1_ses2];

mouse3_RW_d1_ktuned_index_ses1 = intersect(RW_d1_matches_all(5).redCells, RW_d1_ori_all(5).ind_theta90);
mouse3_RW_k_d1_ses1 = RW_d1_k_max_all(5).k1(mouse3_RW_d1_ktuned_index_ses1);
mouse3_RW_max_d1_ses1 = RW_d1_k_max_all(5).max_dfof(mouse3_RW_d1_ktuned_index_ses1);
mouse3_RW_d1_ktuned_index_ses2 = intersect(RW_d1_matches_all(6).redCells, RW_d1_ori_all(6).ind_theta90);
mouse3_RW_k_d1_ses2 = RW_d1_k_max_all(6).k1(mouse3_RW_d1_ktuned_index_ses2);
mouse3_RW_max_d1_ses2 = RW_d1_k_max_all(6).max_dfof(mouse3_RW_d1_ktuned_index_ses2);
mouse3_RW_k_d1_all = [mouse3_RW_k_d1_ses1 mouse3_RW_k_d1_ses2];
mouse3_RW_max_d1_all = [mouse3_RW_max_d1_ses1 mouse3_RW_max_d1_ses2];

mouse4_RW_d1_ktuned_index_ses1 = intersect(RW_d1_matches_all(7).redCells, RW_d1_ori_all(7).ind_theta90);
mouse4_RW_k_d1_ses1 = RW_d1_k_max_all(7).k1(mouse4_RW_d1_ktuned_index_ses1);
mouse4_RW_max_d1_ses1 = RW_d1_k_max_all(7).max_dfof(mouse4_RW_d1_ktuned_index_ses1);
mouse4_RW_d1_ktuned_index_ses2 = intersect(RW_d1_matches_all(8).redCells, RW_d1_ori_all(8).ind_theta90);
mouse4_RW_k_d1_ses2 = RW_d1_k_max_all(8).k1(mouse4_RW_d1_ktuned_index_ses2);
mouse4_RW_max_d1_ses2 = RW_d1_k_max_all(8).max_dfof(mouse4_RW_d1_ktuned_index_ses2);
mouse4_RW_k_d1_all = [mouse4_RW_k_d1_ses1 mouse4_RW_k_d1_ses2];
mouse4_RW_max_d1_all = [mouse4_RW_max_d1_ses1 mouse4_RW_max_d1_ses2];


figure; 
m1_k_d1 = cdfplot(mouse1_k_d1_all);
hold on
m2_k_d1 = cdfplot(mouse2_k_d1_all);
m3_k_d1 = cdfplot(mouse3_k_d1_all);
m4_k_d1 = cdfplot(mouse4_k_d1_all);
set(m1_k_d1, 'Color', 'b', 'LineWidth', 1.0);
set(m2_k_d1, 'Color', 'b', 'LineWidth', 1.0);
set(m3_k_d1, 'Color', 'b', 'LineWidth', 1.0);
set(m4_k_d1, 'Color', 'b', 'LineWidth', 1.0);

m1_RW_k_d1 = cdfplot(mouse1_RW_k_d1_all);
m2_RW_k_d1 = cdfplot(mouse2_RW_k_d1_all);
m3_RW_k_d1 = cdfplot(mouse3_RW_k_d1_all);
m4_RW_k_d1 = cdfplot(mouse4_RW_k_d1_all);
set(m1_RW_k_d1, 'Color', 'r', 'LineWidth', 1.0);
set(m2_RW_k_d1, 'Color', 'r', 'LineWidth', 1.0);
set(m3_RW_k_d1, 'Color', 'r', 'LineWidth', 1.0);
set(m4_RW_k_d1, 'Color', 'r', 'LineWidth', 1.0);

title('')
xlabel('k Value')
xlim([0 30])
legend([m1_k_d1(1), m1_RW_k_d1(1)], 'Arc Promoter', ' RW')
hold off

%print(fullfile(newfnout, ['individual_all_mice', '_d1_k_cdf.pdf']), '-dpdf', '-bestfit')

figure; 
m1_max_d1 = cdfplot(mouse1_max_d1_all);
hold on
m2_max_d1 = cdfplot(mouse2_max_d1_all);
m3_max_d1 = cdfplot(mouse3_max_d1_all);
m4_max_d1 = cdfplot(mouse4_max_d1_all);
set(m1_max_d1, 'Color', 'b', 'LineWidth', 1.0);
set(m2_max_d1, 'Color', 'b', 'LineWidth', 1.0);
set(m3_max_d1, 'Color', 'b', 'LineWidth', 1.0);
set(m4_max_d1, 'Color', 'b', 'LineWidth', 1.0);

m1_RW_max_d1 = cdfplot(mouse1_RW_max_d1_all);
m2_RW_max_d1 = cdfplot(mouse2_RW_max_d1_all);
m3_RW_max_d1 = cdfplot(mouse3_RW_max_d1_all);
m4_RW_max_d1 = cdfplot(mouse4_RW_max_d1_all);
set(m1_RW_max_d1, 'Color', 'r', 'LineWidth', 1.0);
set(m2_RW_max_d1, 'Color', 'r', 'LineWidth', 1.0);
set(m3_RW_max_d1, 'Color', 'r', 'LineWidth', 1.0);
set(m4_RW_max_d1, 'Color', 'r', 'LineWidth', 1.0);

title('')
xlabel('Max dF/F Value')
%xlim([0 30])
legend([m1_max_d1(1), m1_RW_max_d1(1)], 'Arc Promoter', ' RW')
hold off

%print(fullfile(newfnout, ['individual_all_mice', '_d1_max_cdf.pdf']), '-dpdf', '-bestfit')

%now do this for enhancer vs. RW! - also try and separate red and green
%cells for individual mice!