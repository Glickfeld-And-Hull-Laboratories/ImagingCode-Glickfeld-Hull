%CLEAR EVERYTHINg
clear all
clear all global
clc
close all


%% *** REMEMBER TO CHANGE NEW AND REALFNOUT BASED ON LOOSE OR STRICT TUNING CRITERION **AND TJ OR GRACE MICE
%LOAD DATA AND IDENTIFY FOLDERS TO SAVE TO
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_Arc_greenVred'; %folder to save files to
% newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_Arc_greenVred_STRICT'; %folder to save files to
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_pooled_Arc_greenVred_GLmice'; %folder to save files to
% realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_pooled_Arc_greenVred'; %folder to save files to
% realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_pooled_Arc_greenVred_STRICT'; %folder to save files to
dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
eval(dataset); %load dataset


%%
arc_sigmatch_reliability_scores = [];

% list = [7 10 25 28]; %arc promoter list
list = [67 73 76]; %arc enhancer list

for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};

    arc_sigmatch_reliability_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'matchsig_fits']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'matchsig_fits_2.mat']));
        arc_sigmatch_reliability_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'matchsig_fits_2']));
    else
        arc_sigmatch_reliability_ses2 = [];
    end
    arc_sigmatch_reliability_ses_all = [arc_sigmatch_reliability_ses1 arc_sigmatch_reliability_ses2];
    arc_sigmatch_reliability_scores = [arc_sigmatch_reliability_scores arc_sigmatch_reliability_ses_all];
end


lacz_sigmatch_reliability_scores = [];

% list = [1 4 37 40]; %tj lacz list
list = [49 55 61]; %tj lacz list

for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};

    lacz_sigmatch_reliability_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'matchsig_fits']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'matchsig_fits_2.mat']));
        lacz_sigmatch_reliability_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'matchsig_fits_2']));
    else
        lacz_sigmatch_reliability_ses2 = [];
    end
    lacz_sigmatch_reliability_ses_all = [lacz_sigmatch_reliability_ses1 lacz_sigmatch_reliability_ses2];
    lacz_sigmatch_reliability_scores = [lacz_sigmatch_reliability_scores lacz_sigmatch_reliability_ses_all];
end


%%

arc_red_d1_matchsig_reliability_d1_all = [];

for idata = 1:length(arc_sigmatch_reliability_scores)
    arc_red_d1_matchsig_reliability_d1 = arc_sigmatch_reliability_scores(idata).red_d1_matchsig_reliability;
    arc_red_d1_matchsig_reliability_d1_all = [arc_red_d1_matchsig_reliability_d1_all arc_red_d1_matchsig_reliability_d1 ];
end

lacz_red_d1_matchsig_reliability_d1_all = [];

for idata = 1:length(lacz_sigmatch_reliability_scores)
    lacz_red_d1_matchsig_reliability_d1 = lacz_sigmatch_reliability_scores(idata).red_d1_matchsig_reliability;
    lacz_red_d1_matchsig_reliability_d1_all = [lacz_red_d1_matchsig_reliability_d1_all lacz_red_d1_matchsig_reliability_d1 ];
end

%%
figure;
l=cdfplot(arc_red_d1_matchsig_reliability_d1_all);
set(l, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
m=cdfplot(lacz_red_d1_matchsig_reliability_d1_all);
set(m, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
legend(['Arc Enhancer (KRAB)'], ['LacZ (KRAB)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Fit Reliability (degrees)'])
ylabel(['% of cells'])
title(['Matched and Responsive Cells'])
hold off

print(fullfile(realfnout, ['tj_allmice_matchsig_reliability.pdf']), '-dpdf', '-bestfit')
