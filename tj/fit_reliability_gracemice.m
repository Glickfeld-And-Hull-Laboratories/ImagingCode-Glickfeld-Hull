%CLEAR EVERYTHINg
clear all
clear all global
clc
close all



%%
%LOAD DATA AND IDENTIFY FOLDERS TO SAVE TO
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_Arc_greenVred'; %folder to save files to
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_pooled_Arc_greenVred_GLmice'; %folder to save files to
dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%
arc_maxsep_reliability_scores = [];
arc_ksep_reliability_scores = [];
arc_no30k_ksep_prefori_scores = [];
arc_weak_sharp_overlap_scores = [];

list = [67 73 76];
for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};

    arc_maxsep_reliability_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_max_median_sep_reliability']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_max_median_sep_reliability_2.mat']));
        arc_maxsep_reliability_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_max_median_sep_reliability_2']));
    else
        arc_maxsep_reliability_ses2 = [];
    end
    arc_maxsep_reliability_ses_all = [arc_maxsep_reliability_ses1 arc_maxsep_reliability_ses2];
    arc_maxsep_reliability_scores = [arc_maxsep_reliability_scores arc_maxsep_reliability_ses_all];

    arc_ksep_reliability_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_k_median_sep_reliability']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_k_median_sep_reliability_2.mat']));
        arc_ksep_reliability_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_k_median_sep_reliability_2']));
    else
        arc_ksep_reliability_ses2 = [];
    end
    arc_ksep_reliability_ses_all = [arc_ksep_reliability_ses1 arc_ksep_reliability_ses2];
    arc_ksep_reliability_scores = [arc_ksep_reliability_scores arc_ksep_reliability_ses_all];


    arc_no30k_ksep_prefori_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_no30k_k_median_sep_prefori']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_no30k_k_median_sep_prefori_2.mat']));
        arc_no30k_ksep_prefori_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_no30k_k_median_sep_prefori_2']));
    else
        arc_no30k_ksep_prefori_ses2 = [];
    end
    arc_no30k_ksep_prefori_ses_all = [arc_no30k_ksep_prefori_ses1 arc_no30k_ksep_prefori_ses2];
    arc_no30k_ksep_prefori_scores = [arc_no30k_ksep_prefori_scores arc_no30k_ksep_prefori_ses_all];


    arc_weak_sharp_overlap_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer 'weak_sharp_overlap']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer 'weak_sharp_overlap_2.mat']));
        arc_weak_sharp_overlap_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer 'weak_sharp_overlap_2']));
    else
        arc_weak_sharp_overlap_ses2 = [];
    end
    arc_weak_sharp_overlap_ses_all = [arc_weak_sharp_overlap_ses1 arc_weak_sharp_overlap_ses2];
    arc_weak_sharp_overlap_scores = [arc_weak_sharp_overlap_scores arc_weak_sharp_overlap_ses_all];

end


%%

arc_green_sharp_d1_reliability_all = [];
arc_green_sharp_d2_reliability_all = [];
arc_green_sharp_d3_reliability_all = [];

arc_green_wide_d1_reliability_all = [];
arc_green_wide_d2_reliability_all = [];
arc_green_wide_d3_reliability_all = [];

arc_red_sharp_d1_reliability_all = [];
arc_red_sharp_d2_reliability_all = [];
arc_red_sharp_d3_reliability_all = [];

arc_red_wide_d1_reliability_all = [];
arc_red_wide_d2_reliability_all = [];
arc_red_wide_d3_reliability_all = [];

arc_green_top_50_d1_reliability_all = [];
arc_green_top_50_d2_reliability_all = [];
arc_green_top_50_d3_reliability_all = [];

arc_green_bot_50_d1_reliability_all = [];
arc_green_bot_50_d2_reliability_all = [];
arc_green_bot_50_d3_reliability_all = [];

arc_red_top_50_d1_reliability_all = [];
arc_red_top_50_d2_reliability_all = [];
arc_red_top_50_d3_reliability_all = [];

arc_red_bot_50_d1_reliability_all = [];
arc_red_bot_50_d2_reliability_all = [];
arc_red_bot_50_d3_reliability_all = [];





for idata = 1:length(arc_ksep_reliability_scores)
    arc_green_sharp_d1_reliability = arc_ksep_reliability_scores(idata).green_sharp_reliability_d1;
    arc_green_sharp_d1_reliability_all = [arc_green_sharp_d1_reliability_all arc_green_sharp_d1_reliability];
    arc_green_sharp_d2_reliability = arc_ksep_reliability_scores(idata).green_sharp_reliability_d2;
    arc_green_sharp_d2_reliability_all = [arc_green_sharp_d2_reliability_all arc_green_sharp_d2_reliability];
    arc_green_sharp_d3_reliability = arc_ksep_reliability_scores(idata).green_sharp_reliability_d3;
    arc_green_sharp_d3_reliability_all = [arc_green_sharp_d3_reliability_all arc_green_sharp_d3_reliability];
    
    arc_green_wide_d1_reliability = arc_ksep_reliability_scores(idata).green_wide_reliability_d1;
    arc_green_wide_d1_reliability_all = [arc_green_wide_d1_reliability_all arc_green_wide_d1_reliability];
    arc_green_wide_d2_reliability = arc_ksep_reliability_scores(idata).green_wide_reliability_d2;
    arc_green_wide_d2_reliability_all = [arc_green_wide_d2_reliability_all arc_green_wide_d2_reliability];
    arc_green_wide_d3_reliability = arc_ksep_reliability_scores(idata).green_wide_reliability_d3;
    arc_green_wide_d3_reliability_all = [arc_green_wide_d3_reliability_all arc_green_wide_d3_reliability];
    
    arc_red_sharp_d1_reliability = arc_ksep_reliability_scores(idata).red_sharp_reliability_d1;
    arc_red_sharp_d1_reliability_all = [arc_red_sharp_d1_reliability_all arc_red_sharp_d1_reliability];
    arc_red_sharp_d2_reliability = arc_ksep_reliability_scores(idata).red_sharp_reliability_d2;
    arc_red_sharp_d2_reliability_all = [arc_red_sharp_d2_reliability_all arc_red_sharp_d2_reliability];
    arc_red_sharp_d3_reliability = arc_ksep_reliability_scores(idata).red_sharp_reliability_d3;
    arc_red_sharp_d3_reliability_all = [arc_red_sharp_d3_reliability_all arc_red_sharp_d3_reliability];
    
    arc_red_wide_d1_reliability = arc_ksep_reliability_scores(idata).red_wide_reliability_d1;
    arc_red_wide_d1_reliability_all = [arc_red_wide_d1_reliability_all arc_red_wide_d1_reliability];
    arc_red_wide_d2_reliability = arc_ksep_reliability_scores(idata).red_wide_reliability_d2;
    arc_red_wide_d2_reliability_all = [arc_red_wide_d2_reliability_all arc_red_wide_d2_reliability];
    arc_red_wide_d3_reliability = arc_ksep_reliability_scores(idata).red_wide_reliability_d3;
    arc_red_wide_d3_reliability_all = [arc_red_wide_d3_reliability_all arc_red_wide_d3_reliability];


    arc_green_top_50_d1_reliability = arc_maxsep_reliability_scores(idata).green_top_50_reliability_d1;
    arc_green_top_50_d1_reliability_all = [arc_green_top_50_d1_reliability_all arc_green_top_50_d1_reliability];
    arc_green_top_50_d2_reliability = arc_maxsep_reliability_scores(idata).green_top_50_reliability_d2;
    arc_green_top_50_d2_reliability_all = [arc_green_top_50_d2_reliability_all arc_green_top_50_d2_reliability];
    arc_green_top_50_d3_reliability = arc_maxsep_reliability_scores(idata).green_top_50_reliability_d3;
    arc_green_top_50_d3_reliability_all = [arc_green_top_50_d3_reliability_all arc_green_top_50_d3_reliability];
    
    arc_green_bot_50_d1_reliability = arc_maxsep_reliability_scores(idata).green_bot_50_reliability_d1;
    arc_green_bot_50_d1_reliability_all = [arc_green_bot_50_d1_reliability_all arc_green_bot_50_d1_reliability];
    arc_green_bot_50_d2_reliability = arc_maxsep_reliability_scores(idata).green_bot_50_reliability_d2;
    arc_green_bot_50_d2_reliability_all = [arc_green_bot_50_d2_reliability_all arc_green_bot_50_d2_reliability];
    arc_green_bot_50_d3_reliability = arc_maxsep_reliability_scores(idata).green_bot_50_reliability_d3;
    arc_green_bot_50_d3_reliability_all = [arc_green_bot_50_d3_reliability_all arc_green_bot_50_d3_reliability];
    
    arc_red_top_50_d1_reliability = arc_maxsep_reliability_scores(idata).red_top_50_reliability_d1;
    arc_red_top_50_d1_reliability_all = [arc_red_top_50_d1_reliability_all arc_red_top_50_d1_reliability];
    arc_red_top_50_d2_reliability = arc_maxsep_reliability_scores(idata).red_top_50_reliability_d2;
    arc_red_top_50_d2_reliability_all = [arc_red_top_50_d2_reliability_all arc_red_top_50_d2_reliability];
    arc_red_top_50_d3_reliability = arc_maxsep_reliability_scores(idata).red_top_50_reliability_d3;
    arc_red_top_50_d3_reliability_all = [arc_red_top_50_d3_reliability_all arc_red_top_50_d3_reliability];
    
    arc_red_bot_50_d1_reliability = arc_maxsep_reliability_scores(idata).red_bot_50_reliability_d1;
    arc_red_bot_50_d1_reliability_all = [arc_red_bot_50_d1_reliability_all arc_red_bot_50_d1_reliability];
    arc_red_bot_50_d2_reliability = arc_maxsep_reliability_scores(idata).red_bot_50_reliability_d2;
    arc_red_bot_50_d2_reliability_all = [arc_red_bot_50_d2_reliability_all arc_red_bot_50_d2_reliability];
    arc_red_bot_50_d3_reliability = arc_maxsep_reliability_scores(idata).red_bot_50_reliability_d3;
    arc_red_bot_50_d3_reliability_all = [arc_red_bot_50_d3_reliability_all arc_red_bot_50_d3_reliability];

end


%%

arc_green_sharp_reliability_all = [arc_green_sharp_d1_reliability_all arc_green_sharp_d2_reliability_all arc_green_sharp_d3_reliability_all];
arc_green_wide_reliability_all = [arc_green_wide_d1_reliability_all arc_green_wide_d2_reliability_all arc_green_wide_d3_reliability_all];
arc_red_sharp_reliability_all = [arc_red_sharp_d1_reliability_all arc_red_sharp_d2_reliability_all arc_red_sharp_d3_reliability_all];
arc_red_wide_reliability_all = [arc_red_wide_d1_reliability_all arc_red_wide_d2_reliability_all arc_red_wide_d3_reliability_all];

arc_green_top_50_reliability_all = [arc_green_top_50_d1_reliability_all arc_green_top_50_d2_reliability_all arc_green_top_50_d3_reliability_all];
arc_green_bot_50_reliability_all = [arc_green_bot_50_d1_reliability_all arc_green_bot_50_d2_reliability_all arc_green_bot_50_d3_reliability_all];
arc_red_top_50_reliability_all = [arc_red_top_50_d1_reliability_all arc_red_top_50_d2_reliability_all arc_red_top_50_d3_reliability_all];
arc_red_bot_50_reliability_all = [arc_red_bot_50_d1_reliability_all arc_red_bot_50_d2_reliability_all arc_red_bot_50_d3_reliability_all];


%%
lacz_maxsep_reliability_scores = [];
lacz_ksep_reliability_scores = [];
lacz_no30k_ksep_prefori_scores = [];
lacz_weak_sharp_overlap_scores = [];

list = [49 55 61];
for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};

    lacz_maxsep_reliability_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_max_median_sep_reliability']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_max_median_sep_reliability_2.mat']));
        lacz_maxsep_reliability_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_max_median_sep_reliability_2']));
    else
        lacz_maxsep_reliability_ses2 = [];
    end
    lacz_maxsep_reliability_ses_all = [lacz_maxsep_reliability_ses1 lacz_maxsep_reliability_ses2];
    lacz_maxsep_reliability_scores = [lacz_maxsep_reliability_scores lacz_maxsep_reliability_ses_all];

    lacz_ksep_reliability_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_k_median_sep_reliability']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_k_median_sep_reliability_2.mat']));
        lacz_ksep_reliability_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_k_median_sep_reliability_2']));
    else
        lacz_ksep_reliability_ses2 = [];
    end
    lacz_ksep_reliability_ses_all = [lacz_ksep_reliability_ses1 lacz_ksep_reliability_ses2];
    lacz_ksep_reliability_scores = [lacz_ksep_reliability_scores lacz_ksep_reliability_ses_all];


    lacz_no30k_ksep_prefori_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_no30k_k_median_sep_prefori']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_no30k_k_median_sep_prefori_2.mat']));
        lacz_no30k_ksep_prefori_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_no30k_k_median_sep_prefori_2']));
    else
        lacz_no30k_ksep_prefori_ses2 = [];
    end
    lacz_no30k_ksep_prefori_ses_all = [lacz_no30k_ksep_prefori_ses1 lacz_no30k_ksep_prefori_ses2];
    lacz_no30k_ksep_prefori_scores = [lacz_no30k_ksep_prefori_scores lacz_no30k_ksep_prefori_ses_all];


    lacz_weak_sharp_overlap_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer 'weak_sharp_overlap']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer 'weak_sharp_overlap_2.mat']));
        lacz_weak_sharp_overlap_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer 'weak_sharp_overlap_2']));
    else
        lacz_weak_sharp_overlap_ses2 = [];
    end
    lacz_weak_sharp_overlap_ses_all = [lacz_weak_sharp_overlap_ses1 lacz_weak_sharp_overlap_ses2];
    lacz_weak_sharp_overlap_scores = [lacz_weak_sharp_overlap_scores lacz_weak_sharp_overlap_ses_all];

end


%%

lacz_green_sharp_d1_reliability_all = [];
lacz_green_sharp_d2_reliability_all = [];
lacz_green_sharp_d3_reliability_all = [];

lacz_green_wide_d1_reliability_all = [];
lacz_green_wide_d2_reliability_all = [];
lacz_green_wide_d3_reliability_all = [];

lacz_red_sharp_d1_reliability_all = [];
lacz_red_sharp_d2_reliability_all = [];
lacz_red_sharp_d3_reliability_all = [];

lacz_red_wide_d1_reliability_all = [];
lacz_red_wide_d2_reliability_all = [];
lacz_red_wide_d3_reliability_all = [];

lacz_green_top_50_d1_reliability_all = [];
lacz_green_top_50_d2_reliability_all = [];
lacz_green_top_50_d3_reliability_all = [];

lacz_green_bot_50_d1_reliability_all = [];
lacz_green_bot_50_d2_reliability_all = [];
lacz_green_bot_50_d3_reliability_all = [];

lacz_red_top_50_d1_reliability_all = [];
lacz_red_top_50_d2_reliability_all = [];
lacz_red_top_50_d3_reliability_all = [];

lacz_red_bot_50_d1_reliability_all = [];
lacz_red_bot_50_d2_reliability_all = [];
lacz_red_bot_50_d3_reliability_all = [];





for idata = 1:length(lacz_ksep_reliability_scores)
    lacz_green_sharp_d1_reliability = lacz_ksep_reliability_scores(idata).green_sharp_reliability_d1;
    lacz_green_sharp_d1_reliability_all = [lacz_green_sharp_d1_reliability_all lacz_green_sharp_d1_reliability];
    lacz_green_sharp_d2_reliability = lacz_ksep_reliability_scores(idata).green_sharp_reliability_d2;
    lacz_green_sharp_d2_reliability_all = [lacz_green_sharp_d2_reliability_all lacz_green_sharp_d2_reliability];
    lacz_green_sharp_d3_reliability = lacz_ksep_reliability_scores(idata).green_sharp_reliability_d3;
    lacz_green_sharp_d3_reliability_all = [lacz_green_sharp_d3_reliability_all lacz_green_sharp_d3_reliability];
    
    lacz_green_wide_d1_reliability = lacz_ksep_reliability_scores(idata).green_wide_reliability_d1;
    lacz_green_wide_d1_reliability_all = [lacz_green_wide_d1_reliability_all lacz_green_wide_d1_reliability];
    lacz_green_wide_d2_reliability = lacz_ksep_reliability_scores(idata).green_wide_reliability_d2;
    lacz_green_wide_d2_reliability_all = [lacz_green_wide_d2_reliability_all lacz_green_wide_d2_reliability];
    lacz_green_wide_d3_reliability = lacz_ksep_reliability_scores(idata).green_wide_reliability_d3;
    lacz_green_wide_d3_reliability_all = [lacz_green_wide_d3_reliability_all lacz_green_wide_d3_reliability];
    
    lacz_red_sharp_d1_reliability = lacz_ksep_reliability_scores(idata).red_sharp_reliability_d1;
    lacz_red_sharp_d1_reliability_all = [lacz_red_sharp_d1_reliability_all lacz_red_sharp_d1_reliability];
    lacz_red_sharp_d2_reliability = lacz_ksep_reliability_scores(idata).red_sharp_reliability_d2;
    lacz_red_sharp_d2_reliability_all = [lacz_red_sharp_d2_reliability_all lacz_red_sharp_d2_reliability];
    lacz_red_sharp_d3_reliability = lacz_ksep_reliability_scores(idata).red_sharp_reliability_d3;
    lacz_red_sharp_d3_reliability_all = [lacz_red_sharp_d3_reliability_all lacz_red_sharp_d3_reliability];
    
    lacz_red_wide_d1_reliability = lacz_ksep_reliability_scores(idata).red_wide_reliability_d1;
    lacz_red_wide_d1_reliability_all = [lacz_red_wide_d1_reliability_all lacz_red_wide_d1_reliability];
    lacz_red_wide_d2_reliability = lacz_ksep_reliability_scores(idata).red_wide_reliability_d2;
    lacz_red_wide_d2_reliability_all = [lacz_red_wide_d2_reliability_all lacz_red_wide_d2_reliability];
    lacz_red_wide_d3_reliability = lacz_ksep_reliability_scores(idata).red_wide_reliability_d3;
    lacz_red_wide_d3_reliability_all = [lacz_red_wide_d3_reliability_all lacz_red_wide_d3_reliability];


    lacz_green_top_50_d1_reliability = lacz_maxsep_reliability_scores(idata).green_top_50_reliability_d1;
    lacz_green_top_50_d1_reliability_all = [lacz_green_top_50_d1_reliability_all lacz_green_top_50_d1_reliability];
    lacz_green_top_50_d2_reliability = lacz_maxsep_reliability_scores(idata).green_top_50_reliability_d2;
    lacz_green_top_50_d2_reliability_all = [lacz_green_top_50_d2_reliability_all lacz_green_top_50_d2_reliability];
    lacz_green_top_50_d3_reliability = lacz_maxsep_reliability_scores(idata).green_top_50_reliability_d3;
    lacz_green_top_50_d3_reliability_all = [lacz_green_top_50_d3_reliability_all lacz_green_top_50_d3_reliability];
    
    lacz_green_bot_50_d1_reliability = lacz_maxsep_reliability_scores(idata).green_bot_50_reliability_d1;
    lacz_green_bot_50_d1_reliability_all = [lacz_green_bot_50_d1_reliability_all lacz_green_bot_50_d1_reliability];
    lacz_green_bot_50_d2_reliability = lacz_maxsep_reliability_scores(idata).green_bot_50_reliability_d2;
    lacz_green_bot_50_d2_reliability_all = [lacz_green_bot_50_d2_reliability_all lacz_green_bot_50_d2_reliability];
    lacz_green_bot_50_d3_reliability = lacz_maxsep_reliability_scores(idata).green_bot_50_reliability_d3;
    lacz_green_bot_50_d3_reliability_all = [lacz_green_bot_50_d3_reliability_all lacz_green_bot_50_d3_reliability];
    
    lacz_red_top_50_d1_reliability = lacz_maxsep_reliability_scores(idata).red_top_50_reliability_d1;
    lacz_red_top_50_d1_reliability_all = [lacz_red_top_50_d1_reliability_all lacz_red_top_50_d1_reliability];
    lacz_red_top_50_d2_reliability = lacz_maxsep_reliability_scores(idata).red_top_50_reliability_d2;
    lacz_red_top_50_d2_reliability_all = [lacz_red_top_50_d2_reliability_all lacz_red_top_50_d2_reliability];
    lacz_red_top_50_d3_reliability = lacz_maxsep_reliability_scores(idata).red_top_50_reliability_d3;
    lacz_red_top_50_d3_reliability_all = [lacz_red_top_50_d3_reliability_all lacz_red_top_50_d3_reliability];
    
    lacz_red_bot_50_d1_reliability = lacz_maxsep_reliability_scores(idata).red_bot_50_reliability_d1;
    lacz_red_bot_50_d1_reliability_all = [lacz_red_bot_50_d1_reliability_all lacz_red_bot_50_d1_reliability];
    lacz_red_bot_50_d2_reliability = lacz_maxsep_reliability_scores(idata).red_bot_50_reliability_d2;
    lacz_red_bot_50_d2_reliability_all = [lacz_red_bot_50_d2_reliability_all lacz_red_bot_50_d2_reliability];
    lacz_red_bot_50_d3_reliability = lacz_maxsep_reliability_scores(idata).red_bot_50_reliability_d3;
    lacz_red_bot_50_d3_reliability_all = [lacz_red_bot_50_d3_reliability_all lacz_red_bot_50_d3_reliability];

end


%%

lacz_green_sharp_reliability_all = [lacz_green_sharp_d1_reliability_all lacz_green_sharp_d2_reliability_all lacz_green_sharp_d3_reliability_all];
lacz_green_wide_reliability_all = [lacz_green_wide_d1_reliability_all lacz_green_wide_d2_reliability_all lacz_green_wide_d3_reliability_all];
lacz_red_sharp_reliability_all = [lacz_red_sharp_d1_reliability_all lacz_red_sharp_d2_reliability_all lacz_red_sharp_d3_reliability_all];
lacz_red_wide_reliability_all = [lacz_red_wide_d1_reliability_all lacz_red_wide_d2_reliability_all lacz_red_wide_d3_reliability_all];

lacz_green_top_50_reliability_all = [lacz_green_top_50_d1_reliability_all lacz_green_top_50_d2_reliability_all lacz_green_top_50_d3_reliability_all];
lacz_green_bot_50_reliability_all = [lacz_green_bot_50_d1_reliability_all lacz_green_bot_50_d2_reliability_all lacz_green_bot_50_d3_reliability_all];
lacz_red_top_50_reliability_all = [lacz_red_top_50_d1_reliability_all lacz_red_top_50_d2_reliability_all lacz_red_top_50_d3_reliability_all];
lacz_red_bot_50_reliability_all = [lacz_red_bot_50_d1_reliability_all lacz_red_bot_50_d2_reliability_all lacz_red_bot_50_d3_reliability_all];




%%
fig = figure;
h=cdfplot(arc_red_wide_reliability_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(arc_red_sharp_reliability_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
l=cdfplot(lacz_red_wide_reliability_all);
set(l, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
m=cdfplot(lacz_red_sharp_reliability_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Wide (Arc Promoter)'], ['Sharp (Arc Promoter)'], ['Wide (LacZ)'], ['Sharp (LacZ)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Fit Reliability'])
ylabel(['% of cells'])
title(['KRAB Cells'])
hold off
print(fullfile(realfnout, ['gracemice_ksplit_reliability_redcells.pdf']), '-dpdf', '-bestfit')

fig = figure;
h=cdfplot(arc_red_bot_50_reliability_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(arc_red_top_50_reliability_all);
set(j, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
l=cdfplot(lacz_red_bot_50_reliability_all);
set(l, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
m=cdfplot(lacz_red_top_50_reliability_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Bottom 50% (Arc Promoter)'], ['Top 50% (Arc Promoter)'], ['Bottom 50% (LacZ)'], ['Top 50% (LacZ)'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Fit Reliability'])
ylabel(['% of cells'])
title(['KRAB Cells'])
hold off
print(fullfile(realfnout, ['gracemice_maxsplit_reliability_redcells.pdf']), '-dpdf', '-bestfit')


%%
arc_n_green_sharp_weak=0;
arc_n_green_total = 0;
arc_n_red_sharp_weak=0;
arc_n_red_total = 0;

for i=1:length(arc_weak_sharp_overlap_scores)
    arc_n_green_sharp_weak=arc_n_green_sharp_weak + arc_weak_sharp_overlap_scores(i).n_green_sharpest_and_weakest;
    arc_n_green_total = arc_n_green_total + arc_weak_sharp_overlap_scores(i).n_green_total;
    arc_n_red_sharp_weak=arc_n_red_sharp_weak + arc_weak_sharp_overlap_scores(i).n_red_sharpest_and_weakest;
    arc_n_red_total = arc_n_red_total + arc_weak_sharp_overlap_scores(i).n_red_total;
end

lacz_n_green_sharp_weak=0;
lacz_n_green_total = 0;
lacz_n_red_sharp_weak=0;
lacz_n_red_total = 0;

for i=1:length(lacz_weak_sharp_overlap_scores)
    lacz_n_green_sharp_weak=lacz_n_green_sharp_weak + lacz_weak_sharp_overlap_scores(i).n_green_sharpest_and_weakest;
    lacz_n_green_total = lacz_n_green_total + lacz_weak_sharp_overlap_scores(i).n_green_total;
    lacz_n_red_sharp_weak=lacz_n_red_sharp_weak + lacz_weak_sharp_overlap_scores(i).n_red_sharpest_and_weakest;
    lacz_n_red_total = lacz_n_red_total + lacz_weak_sharp_overlap_scores(i).n_red_total;
end


arc_green_percent_overlap = arc_n_green_sharp_weak/arc_n_green_total
arc_red_percent_overlap = arc_n_red_sharp_weak/arc_n_red_total
lacz_green_percent_overlap = lacz_n_green_sharp_weak/lacz_n_green_total
lacz_red_percent_overlap = lacz_n_red_sharp_weak/lacz_n_red_total

figure;
bar([arc_green_percent_overlap arc_red_percent_overlap lacz_green_percent_overlap lacz_red_percent_overlap])
xticklabels({'Arc (Non-Krab)', 'Arc (KRAB)', 'LacZ (Non-Krab)', 'LacZ (KRAB)'})
ylim([0 1])
ylabel('% of Weak/Sharp Overlap')
print(fullfile(realfnout, ['gracemice_weaksharp_overlap_allcells_bar.pdf']), '-dpdf', '-bestfit')

%%

save(fullfile(realfnout, ['_grace_pooled_reliability_overlap_data.mat']), 'arc_red_wide_reliability_all', 'arc_red_sharp_reliability_all', 'lacz_red_wide_reliability_all', 'lacz_red_sharp_reliability_all', ...
    'arc_red_bot_50_reliability_all', 'arc_red_top_50_reliability_all', 'lacz_red_bot_50_reliability_all', 'lacz_red_top_50_reliability_all', 'arc_green_percent_overlap', 'arc_red_percent_overlap', ...
    'lacz_green_percent_overlap', 'lacz_red_percent_overlap')