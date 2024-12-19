%clear everything
clear all
clear all global
close all
clc
%%
%***NEED TO DO POOLED W/I SESSION DRIFT***
%not sure the population analysis is the best way to go about this, as
%there are so few cells which don't really represent a full population
%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_KRAB\multi_day'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_KRAB\red_pooled'; %folder to save files to
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

lacz_d1_ori_all = [];
lacz_d2_ori_all = [];
lacz_d3_ori_all = [];
lacz_d1_matches_all = [];
lacz_d2_matches_all = [];
lacz_d3_matches_all = [];
lacz_d1_tc_all = [];
lacz_d2_tc_all = [];
lacz_d3_tc_all = [];


arc_d1_green_tuned_index_all = [];
arc_d1_red_tuned_index_all = [];

lacz_d1_green_tuned_index_all = [];
lacz_d1_red_tuned_index_all = [];


% arc_d1 = [7 10 19 22 25 28 31 34];
arc_d1 = [7 10];
arc_d2 = arc_d1+1;
arc_d3 = arc_d2+1;

% lacz_d1 = [1 4 13 16 37 40 43 46];
lacz_d1 = [1 4];
lacz_d2 = lacz_d1+1;
lacz_d3 = lacz_d2+1;

arc_prefori_ses_all = [];
arc_prefori_scores = [];
lacz_prefori_ses_all = [];
lacz_prefori_scores = [];


%% 
arc_red_prefori_d1_d2_match_all = [];
arc_red_prefori_d2_match_all = [];
arc_red_prefori_d1_d3_match_all = [];
arc_red_prefori_d3_match_all = [];
arc_pref_dscores_all = [];
arc_pref_dscores_ses_all = [];

arc_red_pref_d_d1_d2_all = [];
arc_red_pref_d_d1_d3_all = [];

arc_matches_all = [];
arc_pref_dscores_all = [];
arc_pref_dscores_newpref_all = [];
arc_tc_corrs_all = [];


arc_list = [7 10 19 25 28 34];
% arc_list = [67 70 73 76];
for iexp = arc_list
    mouse = expt(iexp).mouse;
    if str2double(expt(iexp).img_day) == 1
        mouse = mouse;
    else
        mouse = [mouse '_2'];
    end
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    arc_matches = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'id_matches']));
    arc_pref_dscores = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores']));
    arc_pref_dscores_newpref = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores_newpref']));
    arc_tc_corrs = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'tc_corrs']));
    arc_matches_all = [arc_matches_all arc_matches];
    arc_pref_dscores_all = [arc_pref_dscores_all arc_pref_dscores];
    arc_pref_dscores_newpref_all = [arc_pref_dscores_newpref_all arc_pref_dscores_newpref];
    arc_tc_corrs_all = [arc_tc_corrs_all arc_tc_corrs];
end


arc_red_prefori_d2_d3_match_all = [];
arc_red_prefori_within_d1_all = [];
arc_red_prefori_within_d2_all = [];
arc_red_prefori_within_d3_all = [];

for idata = 1:length(arc_pref_dscores_all)
    arc_red_prefori_d1_d2_match = arc_pref_dscores_all(idata).d_score_prefori_d1_d2;
    arc_red_prefori_d1_d3_match = arc_pref_dscores_all(idata).d_score_prefori_d1_d3;
    arc_red_prefori_d2_d3_match = arc_pref_dscores_all(idata).d_score_prefori_d2_d3;
    arc_red_prefori_d1_d2_match_all = [arc_red_prefori_d1_d2_match_all arc_red_prefori_d1_d2_match];
    arc_red_prefori_d1_d3_match_all = [arc_red_prefori_d1_d3_match_all arc_red_prefori_d1_d3_match];
    arc_red_prefori_d2_d3_match_all = [arc_red_prefori_d2_d3_match_all arc_red_prefori_d2_d3_match];

    arc_red_prefori_within_d1 = arc_pref_dscores_all(idata).d_score_prefori_within_d1;
    arc_red_prefori_within_d2 = arc_pref_dscores_all(idata).d_score_prefori_within_d2;
    arc_red_prefori_within_d3 = arc_pref_dscores_all(idata).d_score_prefori_within_d3;
    arc_red_prefori_within_d1_all = [arc_red_prefori_within_d1_all arc_red_prefori_within_d1];
    arc_red_prefori_within_d2_all = [arc_red_prefori_within_d2_all arc_red_prefori_within_d2];
    arc_red_prefori_within_d3_all = [arc_red_prefori_within_d3_all arc_red_prefori_within_d3];
end

arc_red_newprefori_d1_d2_match_all = [];
arc_red_newprefori_d1_d3_match_all = [];
arc_red_newprefori_d2_d3_match_all =[];
arc_red_newprefori_within_d1_all = [];
arc_red_newprefori_within_d2_all = [];
arc_red_newprefori_within_d3_all = [];

for idata = 1:length(arc_pref_dscores_newpref_all)
    arc_red_newprefori_d1_d2_match = arc_pref_dscores_newpref_all(idata).d_score_new_prefori_d1_d2;
    arc_red_newprefori_d1_d3_match = arc_pref_dscores_newpref_all(idata).d_score_new_prefori_d1_d3;
    arc_red_newprefori_d2_d3_match = arc_pref_dscores_newpref_all(idata).d_score_new_prefori_d2_d3;
    arc_red_newprefori_d1_d2_match_all = [arc_red_newprefori_d1_d2_match_all arc_red_newprefori_d1_d2_match];
    arc_red_newprefori_d1_d3_match_all = [arc_red_newprefori_d1_d3_match_all arc_red_newprefori_d1_d3_match];
    arc_red_newprefori_d2_d3_match_all = [arc_red_newprefori_d2_d3_match_all arc_red_newprefori_d2_d3_match];
    
    arc_red_newprefori_within_d1 = arc_pref_dscores_newpref_all(idata).d_score_new_prefori_within_d1;
    arc_red_newprefori_within_d2 = arc_pref_dscores_newpref_all(idata).d_score_new_prefori_within_d2;
    arc_red_newprefori_within_d3 = arc_pref_dscores_newpref_all(idata).d_score_new_prefori_within_d3;
    arc_red_newprefori_within_d1_all = [arc_red_newprefori_within_d1_all arc_red_newprefori_within_d1];
    arc_red_newprefori_within_d2_all = [arc_red_newprefori_within_d2_all arc_red_newprefori_within_d2];
    arc_red_newprefori_within_d3_all = [arc_red_newprefori_within_d3_all arc_red_newprefori_within_d3];

end

arc_red_corrs_all = [];
arc_n_days_all = [];
arc_n_cells_all = [];

for idata = 1:length(arc_tc_corrs_all)
    arc_red_corr_vals = [arc_tc_corrs_all(idata).abc_corrs(1,2) arc_tc_corrs_all(idata).abc_corrs(2,3) arc_tc_corrs_all(idata).abc_corrs(1,3)];
    arc_red_corrs_all = [arc_red_corrs_all; arc_red_corr_vals];
    arc_n_days = [arc_tc_corrs_all(idata).n_days_d1_d2 arc_tc_corrs_all(idata).n_days_d2_d3 arc_tc_corrs_all(idata).n_days_d1_d3];
    arc_n_days_all = [arc_n_days_all; arc_n_days];
    arc_n_cells = length(arc_matches_all(idata).match_all);
    arc_n_cells_all = [arc_n_cells_all arc_n_cells];
end



%% 
lacz_red_prefori_d1_d2_match_all = [];
lacz_red_prefori_d2_match_all = [];
lacz_red_prefori_d1_d3_match_all = [];
lacz_red_prefori_d3_match_all = [];
lacz_pref_dscores_all = [];
lacz_pref_dscores_ses_all = [];

lacz_red_pref_d_d1_d2_all = [];
lacz_red_pref_d_d1_d3_all = [];

lacz_matches_all = [];
lacz_pref_dscores_all = [];
lacz_pref_dscores_newpref_all = [];
lacz_tc_corrs_all = [];


lacz_list = [1 4 13 16 37 43 46];
% lacz_list = [49 52 55 58 61 64];

for iexp = lacz_list
    mouse = expt(iexp).mouse;
    if str2double(expt(iexp).img_day) == 1
        mouse = mouse;
    else
        mouse = [mouse '_2'];
    end
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    lacz_matches = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'id_matches']));
    lacz_pref_dscores = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores']));
    lacz_pref_dscores_newpref = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores_newpref']));
    lacz_tc_corrs = load(fullfile(realfnout, [mouse '_' img_area '_' img_layer '_' 'tc_corrs']));
    lacz_matches_all = [lacz_matches_all lacz_matches];
    lacz_pref_dscores_all = [lacz_pref_dscores_all lacz_pref_dscores];
    lacz_pref_dscores_newpref_all = [lacz_pref_dscores_newpref_all lacz_pref_dscores_newpref];
    lacz_tc_corrs_all = [lacz_tc_corrs_all lacz_tc_corrs];
end


lacz_red_prefori_d2_d3_match_all = [];
lacz_red_prefori_within_d1_all = [];
lacz_red_prefori_within_d2_all = [];
lacz_red_prefori_within_d3_all = [];

for idata = 1:length(lacz_pref_dscores_all)
    lacz_red_prefori_d1_d2_match = lacz_pref_dscores_all(idata).d_score_prefori_d1_d2;
    lacz_red_prefori_d1_d3_match = lacz_pref_dscores_all(idata).d_score_prefori_d1_d3;
    lacz_red_prefori_d2_d3_match = lacz_pref_dscores_all(idata).d_score_prefori_d2_d3;
    lacz_red_prefori_d1_d2_match_all = [lacz_red_prefori_d1_d2_match_all lacz_red_prefori_d1_d2_match];
    lacz_red_prefori_d1_d3_match_all = [lacz_red_prefori_d1_d3_match_all lacz_red_prefori_d1_d3_match];
    lacz_red_prefori_d2_d3_match_all = [lacz_red_prefori_d2_d3_match_all lacz_red_prefori_d2_d3_match];

    lacz_red_prefori_within_d1 = lacz_pref_dscores_all(idata).d_score_prefori_within_d1;
    lacz_red_prefori_within_d2 = lacz_pref_dscores_all(idata).d_score_prefori_within_d2;
    lacz_red_prefori_within_d3 = lacz_pref_dscores_all(idata).d_score_prefori_within_d3;
    lacz_red_prefori_within_d1_all = [lacz_red_prefori_within_d1_all lacz_red_prefori_within_d1];
    lacz_red_prefori_within_d2_all = [lacz_red_prefori_within_d2_all lacz_red_prefori_within_d2];
    lacz_red_prefori_within_d3_all = [lacz_red_prefori_within_d3_all lacz_red_prefori_within_d3];
end

lacz_red_newprefori_d1_d2_match_all = [];
lacz_red_newprefori_d1_d3_match_all = [];
lacz_red_newprefori_d2_d3_match_all =[];
lacz_red_newprefori_within_d1_all = [];
lacz_red_newprefori_within_d2_all = [];
lacz_red_newprefori_within_d3_all = [];

for idata = 1:length(lacz_pref_dscores_newpref_all)
    lacz_red_newprefori_d1_d2_match = lacz_pref_dscores_newpref_all(idata).d_score_new_prefori_d1_d2;
    lacz_red_newprefori_d1_d3_match = lacz_pref_dscores_newpref_all(idata).d_score_new_prefori_d1_d3;
    lacz_red_newprefori_d2_d3_match = lacz_pref_dscores_newpref_all(idata).d_score_new_prefori_d2_d3;
    lacz_red_newprefori_d1_d2_match_all = [lacz_red_newprefori_d1_d2_match_all lacz_red_newprefori_d1_d2_match];
    lacz_red_newprefori_d1_d3_match_all = [lacz_red_newprefori_d1_d3_match_all lacz_red_newprefori_d1_d3_match];
    lacz_red_newprefori_d2_d3_match_all = [lacz_red_newprefori_d2_d3_match_all lacz_red_newprefori_d2_d3_match];
    
    lacz_red_newprefori_within_d1 = lacz_pref_dscores_newpref_all(idata).d_score_new_prefori_within_d1;
    lacz_red_newprefori_within_d2 = lacz_pref_dscores_newpref_all(idata).d_score_new_prefori_within_d2;
    lacz_red_newprefori_within_d3 = lacz_pref_dscores_newpref_all(idata).d_score_new_prefori_within_d3;
    lacz_red_newprefori_within_d1_all = [lacz_red_newprefori_within_d1_all lacz_red_newprefori_within_d1];
    lacz_red_newprefori_within_d2_all = [lacz_red_newprefori_within_d2_all lacz_red_newprefori_within_d2];
    lacz_red_newprefori_within_d3_all = [lacz_red_newprefori_within_d3_all lacz_red_newprefori_within_d3];

end


lacz_red_corrs_all = [];
lacz_n_days_all = [];
lacz_n_cells_all = [];

for idata = 1:length(lacz_tc_corrs_all)
    lacz_red_corr_vals = [lacz_tc_corrs_all(idata).abc_corrs(1,2) lacz_tc_corrs_all(idata).abc_corrs(2,3) lacz_tc_corrs_all(idata).abc_corrs(1,3)];
    lacz_red_corrs_all = [lacz_red_corrs_all; lacz_red_corr_vals];
    lacz_n_days = [lacz_tc_corrs_all(idata).n_days_d1_d2 lacz_tc_corrs_all(idata).n_days_d2_d3 lacz_tc_corrs_all(idata).n_days_d1_d3];
    lacz_n_days_all = [lacz_n_days_all; lacz_n_days];
    lacz_n_cells = length(lacz_matches_all(idata).match_all);
    lacz_n_cells_all = [lacz_n_cells_all lacz_n_cells];
end

%%
%READY TO START SCATTERS AND DSCORES

%%

%%
%arc v lacz red
figure; 
h = cdfplot(arc_red_prefori_d1_d2_match_all);
hold on
j = cdfplot(arc_red_prefori_d1_d3_match_all);
m = cdfplot(lacz_red_prefori_d1_d2_match_all);
n = cdfplot(lacz_red_prefori_d1_d3_match_all);
% p = cdfplot(enh_dscore_prefori_d1_d2_all);
% q = cdfplot(enh_dscore_prefori_d1_d3_all);
set(h, 'LineStyle', '-', 'Color', 'b', 'LineWidth', 2.0);
set(j, 'LineStyle', '--', 'Color', 'b', 'LineWidth', 2.0);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 2.0);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 2.0);
% set(p, 'LineStyle', '-', 'Color', 'g');
% set(q, 'LineStyle', '--', 'Color', 'g');
%legend('Arc D1-D2', 'Arc D1-D3', 'LacZ D1-D2', 'LacZ D1-D3', 'Arc-Enh D1-D2', 'Arc-Enh D1-D3', 'Location', 'Best')
legend(['Arc Prom D1-D2 (n = ', num2str(length(arc_red_prefori_d1_d2_match_all)), ')'],...
    ['Arc Prom D1-D3 (n = ', num2str(length(arc_red_prefori_d1_d3_match_all)), ')'],...
    ['LacZ D1-D2 (n = ', num2str(length(lacz_red_prefori_d1_d2_match_all)), ')'],...
    ['LacZ D1-D3 (n = ', num2str(length(lacz_red_prefori_d1_d3_match_all)), ')'], 'Location', 'Best')
% legend(['Arc Prom KRAB Day 1-Day 2'],...
%     ['Arc Prom KRAB Day 1-Day 3'],...
%     ['LacZ KRAB Day 1-Day 2'],...
%     ['LacZ KRAB Day 1-Day 3'], 'Location', 'Best')
%title('Change in Pref Ori Values')
title('')
xlabel('Change in Pref Ori')
xlim([0 90])
hold off
% print(fullfile(newfnout, ['tj all mice red cells', '_change_pref_cdf.pdf']), '-dpdf', '-bestfit')
print(fullfile(newfnout, ['grace all mice red cells', '_change_pref_cdf.pdf']), '-dpdf', '-bestfit')

%%
%arc v lacz red new pref
figure; 
h = cdfplot(arc_red_newprefori_d1_d2_match_all);
hold on
j = cdfplot(arc_red_newprefori_d1_d3_match_all);
m = cdfplot(lacz_red_newprefori_d1_d2_match_all);
n = cdfplot(lacz_red_newprefori_d1_d3_match_all);
% p = cdfplot(enh_dscore_prefori_d1_d2_all);
% q = cdfplot(enh_dscore_prefori_d1_d3_all);
set(h, 'LineStyle', '-', 'Color', 'b', 'LineWidth', 2.0);
set(j, 'LineStyle', '--', 'Color', 'b', 'LineWidth', 2.0);
set(m, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 2.0);
set(n, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 2.0);
% set(p, 'LineStyle', '-', 'Color', 'g');
% set(q, 'LineStyle', '--', 'Color', 'g');
%legend('Arc D1-D2', 'Arc D1-D3', 'LacZ D1-D2', 'LacZ D1-D3', 'Arc-Enh D1-D2', 'Arc-Enh D1-D3', 'Location', 'Best')
legend(['Arc Prom D1-D2 (n = ', num2str(length(arc_red_newprefori_d1_d2_match_all)), ')'],...
    ['Arc Prom D1-D3 (n = ', num2str(length(arc_red_newprefori_d1_d3_match_all)), ')'],...
    ['LacZ D1-D2 (n = ', num2str(length(lacz_red_newprefori_d1_d2_match_all)), ')'],...
    ['LacZ D1-D3 (n = ', num2str(length(lacz_red_newprefori_d1_d3_match_all)), ')'], 'Location', 'Best')
% legend(['Arc Prom KRAB Day 1-Day 2'],...
%     ['Arc Prom KRAB Day 1-Day 3'],...
%     ['LacZ KRAB Day 1-Day 2'],...
%     ['LacZ KRAB Day 1-Day 3'], 'Location', 'Best')
%title('Change in Pref Ori Values')
title('')
xlabel('Change in Pref Ori')
xlim([0 90])
hold off
% print(fullfile(newfnout, ['tj all mice red cells', '_change_newpref_cdf.pdf']), '-dpdf', '-bestfit')
print(fullfile(newfnout, ['grace all mice red cells', '_change_newpref_cdf.pdf']), '-dpdf', '-bestfit')

%%

%arc v lacz red
figure; 
h = cdfplot(arc_red_prefori_within_d1_all);
hold on
j = cdfplot(arc_red_prefori_within_d2_all);
m = cdfplot(arc_red_prefori_within_d3_all);
n = cdfplot(lacz_red_prefori_within_d1_all);
p = cdfplot(lacz_red_prefori_within_d2_all);
q = cdfplot(lacz_red_prefori_within_d3_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 2.0);
set(j, 'LineStyle', '--', 'Color', [.4 .4 .4], 'LineWidth', 2.0);
set(m, 'LineStyle', ':', 'Color', [.8 .8 .8], 'LineWidth', 2.0);
set(n, 'LineStyle', '-', 'Color', [1 0 0], 'LineWidth', 2.0);
set(p, 'LineStyle', '--', 'Color', [.8 0 0]);
set(q, 'LineStyle', ':', 'Color', [.4 0 0]);
%legend('Arc D1-D2', 'Arc D1-D3', 'LacZ D1-D2', 'LacZ D1-D3', 'Arc-Enh D1-D2', 'Arc-Enh D1-D3', 'Location', 'Best')
legend(['Arc Prom D1 (n = ', num2str(length(arc_red_prefori_within_d1_all)), ')'],...
    ['Arc Prom D2 (n = ', num2str(length(arc_red_prefori_within_d2_all)), ')'],...
    ['Arc Prom D3 (n = ', num2str(length(arc_red_prefori_within_d3_all)), ')'],...
    ['LacZ D1 (n = ', num2str(length(lacz_red_prefori_within_d1_all)), ')'],...
    ['LacZ D2 (n = ', num2str(length(lacz_red_prefori_within_d2_all)), ')'],...
    ['LacZ D3 (n = ', num2str(length(lacz_red_prefori_within_d3_all)), ')'], 'Location', 'Best')
% legend(['Arc Prom KRAB Day 1-Day 2'],...
%     ['Arc Prom KRAB Day 1-Day 3'],...
%     ['LacZ KRAB Day 1-Day 2'],...
%     ['LacZ KRAB Day 1-Day 3'], 'Location', 'Best')
%title('Change in Pref Ori Values')
title('')
xlabel('Change in Pref Ori w/i Session')
xlim([0 90])
hold off
% print(fullfile(newfnout, ['tj all mice red cells', '_change_withinpref_cdf.pdf']), '-dpdf', '-bestfit')
print(fullfile(newfnout, ['grace all mice red cells', '_change_withinpref_cdf.pdf']), '-dpdf', '-bestfit')

%%

%arc v lacz red
figure; 
h = cdfplot(arc_red_newprefori_within_d1_all);
hold on
j = cdfplot(arc_red_newprefori_within_d2_all);
m = cdfplot(arc_red_newprefori_within_d3_all);
n = cdfplot(lacz_red_newprefori_within_d1_all);
p = cdfplot(lacz_red_newprefori_within_d2_all);
q = cdfplot(lacz_red_newprefori_within_d3_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 2.0);
set(j, 'LineStyle', '--', 'Color', [.4 .4 .4], 'LineWidth', 2.0);
set(m, 'LineStyle', ':', 'Color', [.8 .8 .8], 'LineWidth', 2.0);
set(n, 'LineStyle', '-', 'Color', [1 0 0], 'LineWidth', 2.0);
set(p, 'LineStyle', '--', 'Color', [.8 0 0]);
set(q, 'LineStyle', ':', 'Color', [.4 0 0]);
%legend('Arc D1-D2', 'Arc D1-D3', 'LacZ D1-D2', 'LacZ D1-D3', 'Arc-Enh D1-D2', 'Arc-Enh D1-D3', 'Location', 'Best')
legend(['Arc Prom D1 (n = ', num2str(length(arc_red_newprefori_within_d1_all)), ')'],...
    ['Arc Prom D2 (n = ', num2str(length(arc_red_newprefori_within_d2_all)), ')'],...
    ['Arc Prom D3 (n = ', num2str(length(arc_red_newprefori_within_d3_all)), ')'],...
    ['LacZ D1 (n = ', num2str(length(lacz_red_newprefori_within_d1_all)), ')'],...
    ['LacZ D2 (n = ', num2str(length(lacz_red_newprefori_within_d2_all)), ')'],...
    ['LacZ D3 (n = ', num2str(length(lacz_red_newprefori_within_d3_all)), ')'], 'Location', 'Best')
% legend(['Arc Prom KRAB Day 1-Day 2'],...
%     ['Arc Prom KRAB Day 1-Day 3'],...
%     ['LacZ KRAB Day 1-Day 2'],...
%     ['LacZ KRAB Day 1-Day 3'], 'Location', 'Best')
%title('Change in Pref Ori Values')
title('')
xlabel('Change in Pref Ori w/i Session')
xlim([0 90])
hold off
% print(fullfile(newfnout, ['tj all mice red cells', '_change_newwithinpref_cdf.pdf']), '-dpdf', '-bestfit')
print(fullfile(newfnout, ['grace all mice red cells', '_change_newwithinpref_cdf.pdf']), '-dpdf', '-bestfit')

%%

% figure;
% scatter(arc_n_days_all, arc_red_corrs_all, 'red')
% hold on
% scatter(lacz_n_days_all, lacz_red_corrs_all, 'black')
% xlim([0 14])
% ylim([-0.2 1])
% xlabel('Days Apart')
% ylabel('Correlation')


%%
arc_days_table = [arc_n_days_all(:,1); arc_n_days_all(:,2); arc_n_days_all(:,3)]; 
arc_corrs_table = [arc_red_corrs_all(:,1); arc_red_corrs_all(:,2); arc_red_corrs_all(:,3)];
arc_cells_table = [repmat(arc_n_cells_all',3,1)];
arc_master_table = table(arc_days_table, arc_cells_table, arc_corrs_table);
arc_master_table = arc_master_table(arc_master_table.arc_cells_table>10,:)
arc_mean_table = groupsummary(arc_master_table, 'arc_days_table', [0 2 4 6 7], {'mean', 'median'});

lacz_days_table = [lacz_n_days_all(:,1); lacz_n_days_all(:,2); lacz_n_days_all(:,3)]; 
lacz_corrs_table = [lacz_red_corrs_all(:,1); lacz_red_corrs_all(:,2); lacz_red_corrs_all(:,3)];
lacz_cells_table = [repmat(lacz_n_cells_all',3,1)];
lacz_master_table = table(lacz_days_table, lacz_cells_table, lacz_corrs_table);
lacz_master_table = lacz_master_table(lacz_master_table.lacz_cells_table>10,:)
lacz_mean_table = groupsummary(lacz_master_table, 'lacz_days_table', [0 2 4 6 7], {'mean', 'median'});

%%
figure;
scatter(arc_master_table, "arc_days_table", "arc_corrs_table", 'MarkerEdgeColor', 'red', 'XJitter','density', 'Marker','square')
hold on
scatter(lacz_master_table, "lacz_days_table", "lacz_corrs_table", 'MarkerEdgeColor', 'black', 'XJitter','density')
xlim([1 8])
ylim([0 1])
xlabel('Days Apart')
ylabel('Correlation')
%%
figure; 
scatter(arc_mean_table, "disc_arc_days_table", "median_arc_corrs_table", 'MarkerEdgeColor', 'red')
hold on 
scatter(lacz_mean_table, "disc_lacz_days_table", "median_lacz_corrs_table", 'MarkerEdgeColor', 'black')
xlabel('Session')
xticklabels({'', '1-2', '2-3', '1-3'})
ylabel('Median Pooled TC Corr')
legend(['Arc Prom'], ['LacZ'])

%%
%do within session prefori too***
%do enhancer data too