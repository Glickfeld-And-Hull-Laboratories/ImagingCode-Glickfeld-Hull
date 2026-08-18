%CLEAR EVERYTHINg
clear all
clear all global
clc
close all



%% *** REMEMBER TO CHANGE NEW AND REALFNOUT BASED ON LOOSE OR STRICT TUNING CRITERION
%LOAD DATA AND IDENTIFY FOLDERS TO SAVE TO
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_Arc_allVred_LG'; %folder to save files to
% newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_Arc_allVred_STRICT'; %folder to save files to
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_pooled_Arc_allVred_GLmice_LG'; %folder to save files to
% realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_pooled_Arc_allVred_STRICT'; %folder to save files to
dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
if ~exist(realfnout)
    mkdir(realfnout)
end

%% ARC Enhancer MICE

arc_all_d1_d2_pref_d_all = [];
arc_all_d1_d3_pref_d_all = [];
arc_all_d1_d2_pref_d_tuned_all = [];
arc_all_d1_d3_pref_d_tuned_all = [];

arc_all_d1_k_all = [];
arc_all_d1vd2_k_all = [];
arc_all_d1vd3_k_all = [];
arc_all_d2_k_all = [];
arc_all_d3_k_all = [];

arc_all_d1_k_tuned_all = [];
arc_all_d1vd2_k_tuned_all = [];
arc_all_d1vd3_k_tuned_all = [];
arc_all_d2_k_tuned_all = [];
arc_all_d3_k_tuned_all = [];

arc_all_d1_max_all = [];
arc_all_d1_max_tuned_all = [];
arc_all_d1vd2_max_all = [];
arc_all_d1vd3_max_all = [];
arc_all_d2_max_all = [];
arc_all_d3_max_all = [];

arc_all_d1_reliability_all = [];
arc_all_d1vd2_reliability_all = [];
arc_all_d1vd3_reliability_all = [];
arc_all_d2_reliability_all = [];
arc_all_d3_reliability_all = [];

arc_red_d1_d2_pref_d_all = [];
arc_red_d1_d3_pref_d_all = [];
arc_red_d1_d2_pref_d_tuned_all = [];
arc_red_d1_d3_pref_d_tuned_all = [];

arc_red_d1_k_all = [];
arc_red_d1vd2_k_all = [];
arc_red_d1vd3_k_all = [];
arc_red_d2_k_all = [];
arc_red_d3_k_all = [];

arc_red_d1_k_tuned_all = [];
arc_red_d1vd2_k_tuned_all = [];
arc_red_d1vd3_k_tuned_all = [];
arc_red_d2_k_tuned_all = [];
arc_red_d3_k_tuned_all = [];

arc_all_d1_respEaOri_all = [];
arc_all_d1_respEaOri_tuned_all = [];
arc_red_d1_respEaOri_all = [];
arc_red_d1_respEaOri_tuned_all = [];

arc_all_d1_vonMises_all = [];
arc_all_d1_vonMises_tuned_all = [];
arc_red_d1_vonMises_all = [];
arc_red_d1_vonMises_tuned_all = [];

arc_red_d1_max_all = [];
arc_red_d1_max_tuned_all = [];
arc_red_d1vd2_max_all = [];
arc_red_d1vd3_max_all = [];
arc_red_d2_max_all = [];
arc_red_d3_max_all = [];

arc_red_d1_reliability_all = [];
arc_red_d1vd2_reliability_all = [];
arc_red_d1vd3_reliability_all = [];
arc_red_d2_reliability_all = [];
arc_red_d3_reliability_all = [];

arc_pref_dscores = [];
arc_k_max_scores = [];
arc_reliability_scores= [];

 
list = [67 73 76];
for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};

    arc_pref_dscores_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores_2.mat']));
        arc_pref_dscores_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores_2']));
    else
        arc_pref_dscores_ses_2 = [];
    end
    arc_pref_dscores_ses_all = [arc_pref_dscores_ses_1 arc_pref_dscores_ses_2];
    arc_pref_dscores = [arc_pref_dscores arc_pref_dscores_ses_all];

    arc_k_max_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes_2.mat']));
        arc_k_max_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes_2']));
    else
        arc_k_max_ses_2 = [];
    end
    arc_k_max_ses_all = [arc_k_max_ses_1 arc_k_max_ses_2];
    arc_k_max_scores = [arc_k_max_scores arc_k_max_ses_all];

    arc_reliabilty_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'fit_reliability']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'fit_reliability_2.mat']));
        arc_reliabilty_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'fit_reliability_2']));
    else
        arc_reliabilty_ses_2 = [];
    end
    arc_reliability_ses_all = [arc_reliabilty_ses_1 arc_reliabilty_ses_2];
    arc_reliability_scores = [arc_reliability_scores arc_reliability_ses_all];
end

%%

for idata = 1:length(arc_pref_dscores)

    arc_all_d1_d2_pref_d = arc_pref_dscores(idata).all_d_score_prefori_d1_d2;
    arc_all_d1_d2_pref_d_all = [arc_all_d1_d2_pref_d_all arc_all_d1_d2_pref_d];
    arc_all_d1_d3_pref_d = arc_pref_dscores(idata).all_d_score_prefori_d1_d3;
    arc_all_d1_d3_pref_d_all = [arc_all_d1_d3_pref_d_all arc_all_d1_d3_pref_d];

    arc_red_d1_d2_pref_d = arc_pref_dscores(idata).red_d_score_prefori_d1_d2;
    arc_red_d1_d2_pref_d_all = [arc_red_d1_d2_pref_d_all arc_red_d1_d2_pref_d];
    arc_red_d1_d3_pref_d = arc_pref_dscores(idata).red_d_score_prefori_d1_d3;
    arc_red_d1_d3_pref_d_all = [arc_red_d1_d3_pref_d_all arc_red_d1_d3_pref_d];

    arc_all_d1_d2_pref_d_tuned = arc_pref_dscores(idata).all_d_score_prefori_d1_d2_tuned;
    arc_all_d1_d2_pref_d_tuned_all = [arc_all_d1_d2_pref_d_tuned_all arc_all_d1_d2_pref_d_tuned];
    arc_all_d1_d3_pref_d_tuned = arc_pref_dscores(idata).all_d_score_prefori_d1_d3_tuned;
    arc_all_d1_d3_pref_d_tuned_all = [arc_all_d1_d3_pref_d_tuned_all arc_all_d1_d3_pref_d_tuned];

    arc_red_d1_d2_pref_d_tuned = arc_pref_dscores(idata).red_d_score_prefori_d1_d2_tuned;
    arc_red_d1_d2_pref_d_tuned_all = [arc_red_d1_d2_pref_d_tuned_all arc_red_d1_d2_pref_d_tuned];
    arc_red_d1_d3_pref_d_tuned = arc_pref_dscores(idata).red_d_score_prefori_d1_d3_tuned;
    arc_red_d1_d3_pref_d_tuned_all = [arc_red_d1_d3_pref_d_tuned_all arc_red_d1_d3_pref_d_tuned];

    arc_all_d1_k = arc_k_max_scores(idata).all_d1_k;
    arc_all_d1_k_all = [arc_all_d1_k_all arc_all_d1_k];
    arc_all_d1vd2_k = arc_k_max_scores(idata).all_d1vd2_k;
    arc_all_d1vd2_k_all = [arc_all_d1vd2_k_all arc_all_d1vd2_k];
    arc_all_d1vd3_k = arc_k_max_scores(idata).all_d1vd3_k;
    arc_all_d1vd3_k_all = [arc_all_d1vd3_k_all arc_all_d1vd3_k];
    arc_all_d2_k = arc_k_max_scores(idata).all_d2_k;
    arc_all_d2_k_all = [arc_all_d2_k_all arc_all_d2_k];
    arc_all_d3_k = arc_k_max_scores(idata).all_d3_k;
    arc_all_d3_k_all = [arc_all_d3_k_all arc_all_d3_k];

    arc_all_d1_k_tuned = arc_k_max_scores(idata).all_d1_k_tuned;
    arc_all_d1_k_tuned_all = [arc_all_d1_k_tuned_all arc_all_d1_k_tuned];
    arc_all_d1vd2_k_tuned = arc_k_max_scores(idata).all_d1vd2_k_tuned;
    arc_all_d1vd2_k_tuned_all = [arc_all_d1vd2_k_tuned_all arc_all_d1vd2_k_tuned];
    arc_all_d1vd3_k_tuned = arc_k_max_scores(idata).all_d1vd3_k_tuned;
    arc_all_d1vd3_k_tuned_all = [arc_all_d1vd3_k_tuned_all arc_all_d1vd3_k_tuned];
    arc_all_d2_k_tuned = arc_k_max_scores(idata).all_d2_k_tuned;
    arc_all_d2_k_tuned_all = [arc_all_d2_k_tuned_all arc_all_d2_k_tuned];
    arc_all_d3_k_tuned = arc_k_max_scores(idata).all_d3_k_tuned;
    arc_all_d3_k_tuned_all = [arc_all_d3_k_tuned_all arc_all_d3_k_tuned];

    arc_all_d1_respEaOri = arc_k_max_scores(idata).all_d1_respEaOri;
    arc_all_d1_respEaOri_all = [arc_all_d1_respEaOri_all; arc_all_d1_respEaOri];
    arc_all_d1_respEaOri_tuned = arc_k_max_scores(idata).all_d1_respEaOri_tuned;
    arc_all_d1_respEaOri_tuned_all = [arc_all_d1_respEaOri_tuned_all; arc_all_d1_respEaOri_tuned];

    arc_red_d1_respEaOri = arc_k_max_scores(idata).red_d1_respEaOri;
    arc_red_d1_respEaOri_all = [arc_red_d1_respEaOri_all; arc_red_d1_respEaOri];
    arc_red_d1_respEaOri_tuned = arc_k_max_scores(idata).red_d1_respEaOri_tuned;
    arc_red_d1_respEaOri_tuned_all = [arc_red_d1_respEaOri_tuned_all; arc_red_d1_respEaOri_tuned];

    arc_all_d1_vonMises = arc_k_max_scores(idata).all_d1_vonMises';
    arc_all_d1_vonMises_all = [arc_all_d1_vonMises_all; arc_all_d1_vonMises];
    arc_all_d1_vonMises_tuned = arc_k_max_scores(idata).all_d1_vonMises_tuned';
    arc_all_d1_vonMises_tuned_all = [arc_all_d1_vonMises_tuned_all; arc_all_d1_vonMises_tuned];

    arc_red_d1_vonMises = arc_k_max_scores(idata).red_d1_vonMises';
    arc_red_d1_vonMises_all = [arc_red_d1_vonMises_all; arc_red_d1_vonMises];
    arc_red_d1_vonMises_tuned = arc_k_max_scores(idata).red_d1_vonMises_tuned';
    arc_red_d1_vonMises_tuned_all = [arc_red_d1_vonMises_tuned_all; arc_red_d1_vonMises_tuned];

    arc_red_d1_k = arc_k_max_scores(idata).red_d1_k;
    arc_red_d1_k_all = [arc_red_d1_k_all arc_red_d1_k];
    arc_red_d1vd2_k = arc_k_max_scores(idata).red_d1vd2_k;
    arc_red_d1vd2_k_all = [arc_red_d1vd2_k_all arc_red_d1vd2_k];
    arc_red_d1vd3_k = arc_k_max_scores(idata).red_d1vd3_k;
    arc_red_d1vd3_k_all = [arc_red_d1vd3_k_all arc_red_d1vd3_k];
    arc_red_d2_k = arc_k_max_scores(idata).red_d2_k;
    arc_red_d2_k_all = [arc_red_d2_k_all arc_red_d2_k];
    arc_red_d3_k = arc_k_max_scores(idata).red_d3_k;
    arc_red_d3_k_all = [arc_red_d3_k_all arc_red_d3_k];

    arc_red_d1_k_tuned = arc_k_max_scores(idata).red_d1_k_tuned;
    arc_red_d1_k_tuned_all = [arc_red_d1_k_tuned_all arc_red_d1_k_tuned];
    arc_red_d1vd2_k_tuned = arc_k_max_scores(idata).red_d1vd2_k_tuned;
    arc_red_d1vd2_k_tuned_all = [arc_red_d1vd2_k_tuned_all arc_red_d1vd2_k_tuned];
    arc_red_d1vd3_k_tuned = arc_k_max_scores(idata).red_d1vd3_k_tuned;
    arc_red_d1vd3_k_tuned_all = [arc_red_d1vd3_k_tuned_all arc_red_d1vd3_k_tuned];
    arc_red_d2_k_tuned = arc_k_max_scores(idata).red_d2_k_tuned;
    arc_red_d2_k_tuned_all = [arc_red_d2_k_tuned_all arc_red_d2_k_tuned];
    arc_red_d3_k_tuned = arc_k_max_scores(idata).red_d3_k_tuned;
    arc_red_d3_k_tuned_all = [arc_red_d3_k_tuned_all arc_red_d3_k_tuned];

    arc_all_d1_max = arc_k_max_scores(idata).all_d1_max;
    arc_all_d1_max_all = [arc_all_d1_max_all arc_all_d1_max];
    arc_all_d1_max_tuned = arc_k_max_scores(idata).all_d1_max_tuned;
    arc_all_d1_max_tuned_all = [arc_all_d1_max_tuned_all arc_all_d1_max_tuned];
    arc_all_d1vd2_max = arc_k_max_scores(idata).all_d1vd2_max;
    arc_all_d1vd2_max_all = [arc_all_d1vd2_max_all arc_all_d1vd2_max];
    arc_all_d1vd3_max = arc_k_max_scores(idata).all_d1vd3_max;
    arc_all_d1vd3_max_all = [arc_all_d1vd3_max_all arc_all_d1vd3_max];
    arc_all_d2_max = arc_k_max_scores(idata).all_d2_max;
    arc_all_d2_max_all = [arc_all_d2_max_all arc_all_d2_max];
    arc_all_d3_max = arc_k_max_scores(idata).all_d3_max;
    arc_all_d3_max_all = [arc_all_d3_max_all arc_all_d3_max];

    arc_red_d1_max = arc_k_max_scores(idata).red_d1_max;
    arc_red_d1_max_all = [arc_red_d1_max_all arc_red_d1_max];
    arc_red_d1_max_tuned = arc_k_max_scores(idata).red_d1_max_tuned;
    arc_red_d1_max_tuned_all = [arc_red_d1_max_tuned_all arc_red_d1_max_tuned];
    arc_red_d1vd2_max = arc_k_max_scores(idata).red_d1vd2_max;
    arc_red_d1vd2_max_all = [arc_red_d1vd2_max_all arc_red_d1vd2_max];
    arc_red_d1vd3_max = arc_k_max_scores(idata).red_d1vd3_max;
    arc_red_d1vd3_max_all = [arc_red_d1vd3_max_all arc_red_d1vd3_max];
    arc_red_d2_max = arc_k_max_scores(idata).red_d2_max;
    arc_red_d2_max_all = [arc_red_d2_max_all arc_red_d2_max];
    arc_red_d3_max = arc_k_max_scores(idata).red_d3_max;
    arc_red_d3_max_all = [arc_red_d3_max_all arc_red_d3_max];
    
    arc_all_d1_reliability = arc_reliability_scores(idata).all_d1_reliability;
    arc_all_d1_reliability_all = [arc_all_d1_reliability_all arc_all_d1_reliability];
    arc_all_d1vd2_reliability = arc_reliability_scores(idata).all_d1vd2_reliability;
    arc_all_d1vd2_reliability_all = [arc_all_d1vd2_reliability_all arc_all_d1vd2_reliability];
    arc_all_d1vd3_reliability = arc_reliability_scores(idata).all_d1vd3_reliability;
    arc_all_d1vd3_reliability_all = [arc_all_d1vd3_reliability_all arc_all_d1vd3_reliability];
    arc_all_d2_reliability = arc_reliability_scores(idata).all_d2_reliability;
    arc_all_d2_reliability_all = [arc_all_d2_reliability_all arc_all_d2_reliability];
    arc_all_d3_reliability = arc_reliability_scores(idata).all_d3_reliability;
    arc_all_d3_reliability_all = [arc_all_d3_reliability_all arc_all_d3_reliability];

    arc_red_d1_reliability = arc_reliability_scores(idata).red_d1_reliability;
    arc_red_d1_reliability_all = [arc_red_d1_reliability_all arc_red_d1_reliability];
    arc_red_d1vd2_reliability = arc_reliability_scores(idata).red_d1vd2_reliability;
    arc_red_d1vd2_reliability_all = [arc_red_d1vd2_reliability_all arc_red_d1vd2_reliability];
    arc_red_d1vd3_reliability = arc_reliability_scores(idata).red_d1vd3_reliability;
    arc_red_d1vd3_reliability_all = [arc_red_d1vd3_reliability_all arc_red_d1vd3_reliability];
    arc_red_d2_reliability = arc_reliability_scores(idata).red_d2_reliability;
    arc_red_d2_reliability_all = [arc_red_d2_reliability_all arc_red_d2_reliability];
    arc_red_d3_reliability = arc_reliability_scores(idata).red_d3_reliability;
    arc_red_d3_reliability_all = [arc_red_d3_reliability_all arc_red_d3_reliability];
end

%% SAME AS ABOVE BUT LacZ

lacz_all_d1_d2_pref_d_all = [];
lacz_all_d1_d3_pref_d_all = [];
lacz_all_d1_d2_pref_d_tuned_all = [];
lacz_all_d1_d3_pref_d_tuned_all = [];

lacz_all_d1_k_all = [];
lacz_all_d1vd2_k_all = [];
lacz_all_d1vd3_k_all = [];
lacz_all_d2_k_all = [];
lacz_all_d3_k_all = [];

lacz_all_d1_k_tuned_all = [];
lacz_all_d1vd2_k_tuned_all = [];
lacz_all_d1vd3_k_tuned_all = [];
lacz_all_d2_k_tuned_all = [];
lacz_all_d3_k_tuned_all = [];

lacz_all_d1_max_all = [];
lacz_all_d1_max_tuned_all = [];
lacz_all_d1vd2_max_all = [];
lacz_all_d1vd3_max_all = [];
lacz_all_d2_max_all = [];
lacz_all_d3_max_all = [];

lacz_all_d1_reliability_all = [];
lacz_all_d1vd2_reliability_all = [];
lacz_all_d1vd3_reliability_all = [];
lacz_all_d2_reliability_all = [];
lacz_all_d3_reliability_all = [];

lacz_red_d1_d2_pref_d_all = [];
lacz_red_d1_d3_pref_d_all = [];
lacz_red_d1_d2_pref_d_tuned_all = [];
lacz_red_d1_d3_pref_d_tuned_all = [];

lacz_red_d1_k_all = [];
lacz_red_d1vd2_k_all = [];
lacz_red_d1vd3_k_all = [];
lacz_red_d2_k_all = [];
lacz_red_d3_k_all = [];

lacz_red_d1_k_tuned_all = [];
lacz_red_d1vd2_k_tuned_all = [];
lacz_red_d1vd3_k_tuned_all = [];
lacz_red_d2_k_tuned_all = [];
lacz_red_d3_k_tuned_all = [];

lacz_all_d1_respEaOri_all = [];
lacz_all_d1_respEaOri_tuned_all = [];
lacz_red_d1_respEaOri_all = [];
lacz_red_d1_respEaOri_tuned_all = [];

lacz_all_d1_vonMises_all = [];
lacz_all_d1_vonMises_tuned_all = [];
lacz_red_d1_vonMises_all = [];
lacz_red_d1_vonMises_tuned_all = [];

lacz_red_d1_max_all = [];
lacz_red_d1_max_tuned_all = [];
lacz_red_d1vd2_max_all = [];
lacz_red_d1vd3_max_all = [];
lacz_red_d2_max_all = [];
lacz_red_d3_max_all = [];

lacz_red_d1_reliability_all = [];
lacz_red_d1vd2_reliability_all = [];
lacz_red_d1vd3_reliability_all = [];
lacz_red_d2_reliability_all = [];
lacz_red_d3_reliability_all = [];

lacz_pref_dscores = [];
lacz_k_max_scores = [];
lacz_reliability_scores= [];

list = [49 55 31 79];
for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};

    lacz_pref_dscores_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores_2.mat']));
        lacz_pref_dscores_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores_2']));
    else
        lacz_pref_dscores_ses_2 = [];
    end
    lacz_pref_dscores_ses_all = [lacz_pref_dscores_ses_1 lacz_pref_dscores_ses_2];
    lacz_pref_dscores = [lacz_pref_dscores lacz_pref_dscores_ses_all];

    lacz_k_max_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes_2.mat']));
        lacz_k_max_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'k_max_changes_2']));
    else
        lacz_k_max_ses_2 = [];
    end
    lacz_k_max_ses_all = [lacz_k_max_ses_1 lacz_k_max_ses_2];
    lacz_k_max_scores = [lacz_k_max_scores lacz_k_max_ses_all];

    lacz_reliabilty_ses_1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'fit_reliability']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'fit_reliability_2.mat']));
        lacz_reliabilty_ses_2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'fit_reliability_2']));
    else
        lacz_reliabilty_ses_2 = [];
    end
    lacz_reliability_ses_all = [lacz_reliabilty_ses_1 lacz_reliabilty_ses_2];
    lacz_reliability_scores = [lacz_reliability_scores lacz_reliability_ses_all];
end



%%

for idata = 1:length(lacz_pref_dscores)
    lacz_all_d1_d2_pref_d = lacz_pref_dscores(idata).all_d_score_prefori_d1_d2;
    lacz_all_d1_d2_pref_d_all = [lacz_all_d1_d2_pref_d_all lacz_all_d1_d2_pref_d];
    lacz_all_d1_d3_pref_d = lacz_pref_dscores(idata).all_d_score_prefori_d1_d3;
    lacz_all_d1_d3_pref_d_all = [lacz_all_d1_d3_pref_d_all lacz_all_d1_d3_pref_d];

    lacz_red_d1_d2_pref_d = lacz_pref_dscores(idata).red_d_score_prefori_d1_d2;
    lacz_red_d1_d2_pref_d_all = [lacz_red_d1_d2_pref_d_all lacz_red_d1_d2_pref_d];
    lacz_red_d1_d3_pref_d = lacz_pref_dscores(idata).red_d_score_prefori_d1_d3;
    lacz_red_d1_d3_pref_d_all = [lacz_red_d1_d3_pref_d_all lacz_red_d1_d3_pref_d];

    lacz_all_d1_d2_pref_d_tuned = lacz_pref_dscores(idata).all_d_score_prefori_d1_d2_tuned;
    lacz_all_d1_d2_pref_d_tuned_all = [lacz_all_d1_d2_pref_d_tuned_all lacz_all_d1_d2_pref_d_tuned];
    lacz_all_d1_d3_pref_d_tuned = lacz_pref_dscores(idata).all_d_score_prefori_d1_d3_tuned;
    lacz_all_d1_d3_pref_d_tuned_all = [lacz_all_d1_d3_pref_d_tuned_all lacz_all_d1_d3_pref_d_tuned];

    lacz_red_d1_d2_pref_d_tuned = lacz_pref_dscores(idata).red_d_score_prefori_d1_d2_tuned;
    lacz_red_d1_d2_pref_d_tuned_all = [lacz_red_d1_d2_pref_d_tuned_all lacz_red_d1_d2_pref_d_tuned];
    lacz_red_d1_d3_pref_d_tuned = lacz_pref_dscores(idata).red_d_score_prefori_d1_d3_tuned;
    lacz_red_d1_d3_pref_d_tuned_all = [lacz_red_d1_d3_pref_d_tuned_all lacz_red_d1_d3_pref_d_tuned];

    lacz_all_d1_k = lacz_k_max_scores(idata).all_d1_k;
    lacz_all_d1_k_all = [lacz_all_d1_k_all lacz_all_d1_k];
    lacz_all_d1vd2_k = lacz_k_max_scores(idata).all_d1vd2_k;
    lacz_all_d1vd2_k_all = [lacz_all_d1vd2_k_all lacz_all_d1vd2_k];
    lacz_all_d1vd3_k = lacz_k_max_scores(idata).all_d1vd3_k;
    lacz_all_d1vd3_k_all = [lacz_all_d1vd3_k_all lacz_all_d1vd3_k];
    lacz_all_d2_k = lacz_k_max_scores(idata).all_d2_k;
    lacz_all_d2_k_all = [lacz_all_d2_k_all lacz_all_d2_k];
    lacz_all_d3_k = lacz_k_max_scores(idata).all_d3_k;
    lacz_all_d3_k_all = [lacz_all_d3_k_all lacz_all_d3_k];

    lacz_all_d1_k_tuned = lacz_k_max_scores(idata).all_d1_k_tuned;
    lacz_all_d1_k_tuned_all = [lacz_all_d1_k_tuned_all lacz_all_d1_k_tuned];
    lacz_all_d1vd2_k_tuned = lacz_k_max_scores(idata).all_d1vd2_k_tuned;
    lacz_all_d1vd2_k_tuned_all = [lacz_all_d1vd2_k_tuned_all lacz_all_d1vd2_k_tuned];
    lacz_all_d1vd3_k_tuned = lacz_k_max_scores(idata).all_d1vd3_k_tuned;
    lacz_all_d1vd3_k_tuned_all = [lacz_all_d1vd3_k_tuned_all lacz_all_d1vd3_k_tuned];
    lacz_all_d2_k_tuned = lacz_k_max_scores(idata).all_d2_k_tuned;
    lacz_all_d2_k_tuned_all = [lacz_all_d2_k_tuned_all lacz_all_d2_k_tuned];
    lacz_all_d3_k_tuned = lacz_k_max_scores(idata).all_d3_k_tuned;
    lacz_all_d3_k_tuned_all = [lacz_all_d3_k_tuned_all lacz_all_d3_k_tuned];

    lacz_all_d1_respEaOri = lacz_k_max_scores(idata).all_d1_respEaOri;
    lacz_all_d1_respEaOri_all = [lacz_all_d1_respEaOri_all; lacz_all_d1_respEaOri];
    lacz_all_d1_respEaOri_tuned = lacz_k_max_scores(idata).all_d1_respEaOri_tuned;
    lacz_all_d1_respEaOri_tuned_all = [lacz_all_d1_respEaOri_tuned_all; lacz_all_d1_respEaOri_tuned];

    lacz_red_d1_respEaOri = lacz_k_max_scores(idata).red_d1_respEaOri;
    lacz_red_d1_respEaOri_all = [lacz_red_d1_respEaOri_all; lacz_red_d1_respEaOri];
    lacz_red_d1_respEaOri_tuned = lacz_k_max_scores(idata).red_d1_respEaOri_tuned;
    lacz_red_d1_respEaOri_tuned_all = [lacz_red_d1_respEaOri_tuned_all; lacz_red_d1_respEaOri_tuned];

    lacz_all_d1_vonMises = lacz_k_max_scores(idata).all_d1_vonMises';
    lacz_all_d1_vonMises_all = [lacz_all_d1_vonMises_all; lacz_all_d1_vonMises];
    lacz_all_d1_vonMises_tuned = lacz_k_max_scores(idata).all_d1_vonMises_tuned';
    lacz_all_d1_vonMises_tuned_all = [lacz_all_d1_vonMises_tuned_all; lacz_all_d1_vonMises_tuned];

    lacz_red_d1_vonMises = lacz_k_max_scores(idata).red_d1_vonMises';
    lacz_red_d1_vonMises_all = [lacz_red_d1_vonMises_all; lacz_red_d1_vonMises];
    lacz_red_d1_vonMises_tuned = lacz_k_max_scores(idata).red_d1_vonMises_tuned';
    lacz_red_d1_vonMises_tuned_all = [lacz_red_d1_vonMises_tuned_all; lacz_red_d1_vonMises_tuned];

    lacz_red_d1_k = lacz_k_max_scores(idata).red_d1_k;
    lacz_red_d1_k_all = [lacz_red_d1_k_all lacz_red_d1_k];
    lacz_red_d1vd2_k = lacz_k_max_scores(idata).red_d1vd2_k;
    lacz_red_d1vd2_k_all = [lacz_red_d1vd2_k_all lacz_red_d1vd2_k];
    lacz_red_d1vd3_k = lacz_k_max_scores(idata).red_d1vd3_k;
    lacz_red_d1vd3_k_all = [lacz_red_d1vd3_k_all lacz_red_d1vd3_k];
    lacz_red_d2_k = lacz_k_max_scores(idata).red_d2_k;
    lacz_red_d2_k_all = [lacz_red_d2_k_all lacz_red_d2_k];
    lacz_red_d3_k = lacz_k_max_scores(idata).red_d3_k;
    lacz_red_d3_k_all = [lacz_red_d3_k_all lacz_red_d3_k];

    lacz_red_d1_k_tuned = lacz_k_max_scores(idata).red_d1_k_tuned;
    lacz_red_d1_k_tuned_all = [lacz_red_d1_k_tuned_all lacz_red_d1_k_tuned];
    lacz_red_d1vd2_k_tuned = lacz_k_max_scores(idata).red_d1vd2_k_tuned;
    lacz_red_d1vd2_k_tuned_all = [lacz_red_d1vd2_k_tuned_all lacz_red_d1vd2_k_tuned];
    lacz_red_d1vd3_k_tuned = lacz_k_max_scores(idata).red_d1vd3_k_tuned;
    lacz_red_d1vd3_k_tuned_all = [lacz_red_d1vd3_k_tuned_all lacz_red_d1vd3_k_tuned];
    lacz_red_d2_k_tuned = lacz_k_max_scores(idata).red_d2_k_tuned;
    lacz_red_d2_k_tuned_all = [lacz_red_d2_k_tuned_all lacz_red_d2_k_tuned];
    lacz_red_d3_k_tuned = lacz_k_max_scores(idata).red_d3_k_tuned;
    lacz_red_d3_k_tuned_all = [lacz_red_d3_k_tuned_all lacz_red_d3_k_tuned];

    lacz_all_d1_max = lacz_k_max_scores(idata).all_d1_max;
    lacz_all_d1_max_all = [lacz_all_d1_max_all lacz_all_d1_max];
    lacz_all_d1_max_tuned = lacz_k_max_scores(idata).all_d1_max_tuned;
    lacz_all_d1_max_tuned_all = [lacz_all_d1_max_tuned_all lacz_all_d1_max_tuned];
    lacz_all_d1vd2_max = lacz_k_max_scores(idata).all_d1vd2_max;
    lacz_all_d1vd2_max_all = [lacz_all_d1vd2_max_all lacz_all_d1vd2_max];
    lacz_all_d1vd3_max = lacz_k_max_scores(idata).all_d1vd3_max;
    lacz_all_d1vd3_max_all = [lacz_all_d1vd3_max_all lacz_all_d1vd3_max];
    lacz_all_d2_max = lacz_k_max_scores(idata).all_d2_max;
    lacz_all_d2_max_all = [lacz_all_d2_max_all lacz_all_d2_max];
    lacz_all_d3_max = lacz_k_max_scores(idata).all_d3_max;
    lacz_all_d3_max_all = [lacz_all_d3_max_all lacz_all_d3_max];

    lacz_red_d1_max = lacz_k_max_scores(idata).red_d1_max;
    lacz_red_d1_max_all = [lacz_red_d1_max_all lacz_red_d1_max];
    lacz_red_d1_max_tuned = lacz_k_max_scores(idata).red_d1_max_tuned;
    lacz_red_d1_max_tuned_all = [lacz_red_d1_max_tuned_all lacz_red_d1_max_tuned];
    lacz_red_d1vd2_max = lacz_k_max_scores(idata).red_d1vd2_max;
    lacz_red_d1vd2_max_all = [lacz_red_d1vd2_max_all lacz_red_d1vd2_max];
    lacz_red_d1vd3_max = lacz_k_max_scores(idata).red_d1vd3_max;
    lacz_red_d1vd3_max_all = [lacz_red_d1vd3_max_all lacz_red_d1vd3_max];
    lacz_red_d2_max = lacz_k_max_scores(idata).red_d2_max;
    lacz_red_d2_max_all = [lacz_red_d2_max_all lacz_red_d2_max];
    lacz_red_d3_max = lacz_k_max_scores(idata).red_d3_max;
    lacz_red_d3_max_all = [lacz_red_d3_max_all lacz_red_d3_max];
    
    lacz_all_d1_reliability = lacz_reliability_scores(idata).all_d1_reliability;
    lacz_all_d1_reliability_all = [lacz_all_d1_reliability_all lacz_all_d1_reliability];
    lacz_all_d1vd2_reliability = lacz_reliability_scores(idata).all_d1vd2_reliability;
    lacz_all_d1vd2_reliability_all = [lacz_all_d1vd2_reliability_all lacz_all_d1vd2_reliability];
    lacz_all_d1vd3_reliability = lacz_reliability_scores(idata).all_d1vd3_reliability;
    lacz_all_d1vd3_reliability_all = [lacz_all_d1vd3_reliability_all lacz_all_d1vd3_reliability];
    lacz_all_d2_reliability = lacz_reliability_scores(idata).all_d2_reliability;
    lacz_all_d2_reliability_all = [lacz_all_d2_reliability_all lacz_all_d2_reliability];
    lacz_all_d3_reliability = lacz_reliability_scores(idata).all_d3_reliability;
    lacz_all_d3_reliability_all = [lacz_all_d3_reliability_all lacz_all_d3_reliability];

    lacz_red_d1_reliability = lacz_reliability_scores(idata).red_d1_reliability;
    lacz_red_d1_reliability_all = [lacz_red_d1_reliability_all lacz_red_d1_reliability];
    lacz_red_d1vd2_reliability = lacz_reliability_scores(idata).red_d1vd2_reliability;
    lacz_red_d1vd2_reliability_all = [lacz_red_d1vd2_reliability_all lacz_red_d1vd2_reliability];
    lacz_red_d1vd3_reliability = lacz_reliability_scores(idata).red_d1vd3_reliability;
    lacz_red_d1vd3_reliability_all = [lacz_red_d1vd3_reliability_all lacz_red_d1vd3_reliability];
    lacz_red_d2_reliability = lacz_reliability_scores(idata).red_d2_reliability;
    lacz_red_d2_reliability_all = [lacz_red_d2_reliability_all lacz_red_d2_reliability];
    lacz_red_d3_reliability = lacz_reliability_scores(idata).red_d3_reliability;
    lacz_red_d3_reliability_all = [lacz_red_d3_reliability_all lacz_red_d3_reliability];
end

%% compare all groups from mine and grace's
%k and max change scores

arc_red_d1_d2_k_all = abs(arc_red_d1vd2_k_all-arc_red_d2_k_all);
arc_red_d1_d3_k_all = abs(arc_red_d1vd3_k_all-arc_red_d3_k_all);
arc_red_d1_d2_k_tuned_all = abs(arc_red_d1vd2_k_tuned_all-arc_red_d2_k_tuned_all);
arc_red_d1_d3_k_tuned_all = abs(arc_red_d1vd3_k_tuned_all-arc_red_d3_k_tuned_all);

lacz_red_d1_d2_k_all = abs(lacz_red_d1vd2_k_all-lacz_red_d2_k_all);
lacz_red_d1_d3_k_all = abs(lacz_red_d1vd3_k_all-lacz_red_d3_k_all);
lacz_red_d1_d2_k_tuned_all = abs(lacz_red_d1vd2_k_tuned_all-lacz_red_d2_k_tuned_all);
lacz_red_d1_d3_k_tuned_all = abs(lacz_red_d1vd3_k_tuned_all-lacz_red_d3_k_tuned_all);

arc_all_d1_d2_k_all = abs(arc_all_d1vd2_k_all-arc_all_d2_k_all);
arc_all_d1_d3_k_all = abs(arc_all_d1vd3_k_all-arc_all_d3_k_all);
arc_all_d1_d2_k_tuned_all = abs(arc_all_d1vd2_k_tuned_all-arc_all_d2_k_tuned_all);
arc_all_d1_d3_k_tuned_all = abs(arc_all_d1vd3_k_tuned_all-arc_all_d3_k_tuned_all);

lacz_all_d1_d2_k_all = abs(lacz_all_d1vd2_k_all-lacz_all_d2_k_all);
lacz_all_d1_d3_k_all = abs(lacz_all_d1vd3_k_all-lacz_all_d3_k_all);
lacz_all_d1_d2_k_tuned_all = abs(lacz_all_d1vd2_k_tuned_all-lacz_all_d2_k_tuned_all);
lacz_all_d1_d3_k_tuned_all = abs(lacz_all_d1vd3_k_tuned_all-lacz_all_d3_k_tuned_all);

arc_red_d1_d2_max_all = abs(arc_red_d1vd2_max_all-arc_red_d2_max_all);
arc_red_d1_d3_max_all = abs(arc_red_d1vd3_max_all-arc_red_d3_max_all);
lacz_red_d1_d2_max_all = abs(lacz_red_d1vd2_max_all-lacz_red_d2_max_all);
lacz_red_d1_d3_max_all = abs(lacz_red_d1vd3_max_all-lacz_red_d3_max_all);

arc_all_d1_d2_max_all = abs(arc_all_d1vd2_max_all-arc_all_d2_max_all);
arc_all_d1_d3_max_all = abs(arc_all_d1vd3_max_all-arc_all_d3_max_all);
lacz_all_d1_d2_max_all = abs(lacz_all_d1vd2_max_all-lacz_all_d2_max_all);
lacz_all_d1_d3_max_all = abs(lacz_all_d1vd3_max_all-lacz_all_d3_max_all);


%%

save(fullfile(realfnout, ['_grace_pooled_pref_data.mat']), 'arc_red_d1_d2_pref_d_all', 'arc_red_d1_d3_pref_d_all', 'arc_all_d1_d2_pref_d_all', 'arc_all_d1_d3_pref_d_all', ...
    'lacz_red_d1_d2_pref_d_all', 'lacz_red_d1_d3_pref_d_all', 'lacz_all_d1_d2_pref_d_all', 'lacz_all_d1_d3_pref_d_all', ...
   'arc_red_d1_d2_pref_d_tuned_all', 'arc_red_d1_d3_pref_d_tuned_all', 'arc_all_d1_d2_pref_d_tuned_all', 'arc_all_d1_d3_pref_d_tuned_all', ...
    'lacz_red_d1_d2_pref_d_tuned_all', 'lacz_red_d1_d3_pref_d_tuned_all', 'lacz_all_d1_d2_pref_d_tuned_all', 'lacz_all_d1_d3_pref_d_tuned_all')

save(fullfile(realfnout, ['_grace_pooled_kmax_data.mat']), 'arc_red_d1_d2_k_all', 'arc_all_d1_d2_k_all','lacz_red_d1_d2_k_all', 'lacz_all_d1_d2_k_all',...
    'arc_red_d1_d2_max_all', 'arc_all_d1_d2_max_all', 'lacz_red_d1_d2_max_all', 'lacz_all_d1_d2_max_all',...
    'arc_red_d1_d3_k_all', 'arc_all_d1_d3_k_all','lacz_red_d1_d3_k_all', 'lacz_all_d1_d3_k_all',...
    'arc_red_d1_d3_max_all', 'arc_all_d1_d3_max_all', 'lacz_red_d1_d3_max_all', 'lacz_all_d1_d3_max_all',...
    'arc_red_d1_d2_k_tuned_all', 'arc_all_d1_d2_k_tuned_all','lacz_red_d1_d2_k_tuned_all', 'lacz_all_d1_d2_k_tuned_all',...
    'arc_red_d1_d3_k_tuned_all', 'arc_all_d1_d3_k_tuned_all','lacz_red_d1_d3_k_tuned_all', 'lacz_all_d1_d3_k_tuned_all')

save(fullfile(realfnout, ['_grace_pooled_kmax_day1_data.mat']), 'arc_red_d1_k_all', 'arc_all_d1_k_all','lacz_red_d1_k_all', 'lacz_all_d1_k_all',...
    'arc_red_d1_k_tuned_all', 'arc_all_d1_k_tuned_all','lacz_red_d1_k_tuned_all', 'lacz_all_d1_k_tuned_all',...
    'arc_red_d1_max_all', 'arc_all_d1_max_all','lacz_red_d1_max_all', 'lacz_all_d1_max_all',...
    'arc_red_d1_max_tuned_all', 'arc_all_d1_max_tuned_all','lacz_red_d1_max_tuned_all', 'lacz_all_d1_max_tuned_all',...
    'arc_red_d1_reliability_all', 'lacz_red_d1_reliability_all',...
    'arc_all_d1_respEaOri_all', 'arc_red_d1_respEaOri_all','lacz_all_d1_respEaOri_all', 'lacz_red_d1_respEaOri_all', ...
    'arc_all_d1_respEaOri_tuned_all', 'arc_red_d1_respEaOri_tuned_all','lacz_all_d1_respEaOri_tuned_all', 'lacz_red_d1_respEaOri_tuned_all',...
    'arc_all_d1_vonMises_all', 'arc_red_d1_vonMises_all','lacz_all_d1_vonMises_all', 'lacz_red_d1_vonMises_all', ...
    'arc_all_d1_vonMises_tuned_all', 'arc_red_d1_vonMises_tuned_all','lacz_all_d1_vonMises_tuned_all', 'lacz_red_d1_vonMises_tuned_all')

%%
%binning fit reliability

arc_red_bin1_d1vd2_reliability_ind = find(arc_red_d1vd2_reliability_all >=0 & arc_red_d1vd2_reliability_all <=10);
arc_red_bin2_d1vd2_reliability_ind = find(arc_red_d1vd2_reliability_all >=11 & arc_red_d1vd2_reliability_all <=30);
arc_red_bin3_d1vd2_reliability_ind = find(arc_red_d1vd2_reliability_all >=31 & arc_red_d1vd2_reliability_all <=90);

arc_red_d1_d2_pref_bin1_reliability = arc_red_d1_d2_pref_d_all(arc_red_bin1_d1vd2_reliability_ind);
arc_red_d1_d2_pref_bin2_reliability = arc_red_d1_d2_pref_d_all(arc_red_bin2_d1vd2_reliability_ind);
arc_red_d1_d2_pref_bin3_reliability = arc_red_d1_d2_pref_d_all(arc_red_bin3_d1vd2_reliability_ind);

arc_red_bin1_d1d3_reliability_ind = find(arc_red_d1vd3_reliability_all >=0 & arc_red_d1vd3_reliability_all <=10);
arc_red_bin2_d1d3_reliability_ind = find(arc_red_d1vd3_reliability_all >=11 & arc_red_d1vd3_reliability_all <=30);
arc_red_bin3_d1d3_reliability_ind = find(arc_red_d1vd3_reliability_all >=31 & arc_red_d1vd3_reliability_all <=90);

arc_red_d1_d3_pref_bin1_reliability = arc_red_d1_d3_pref_d_all(arc_red_bin1_d1d3_reliability_ind);
arc_red_d1_d3_pref_bin2_reliability = arc_red_d1_d3_pref_d_all(arc_red_bin2_d1d3_reliability_ind);
arc_red_d1_d3_pref_bin3_reliability = arc_red_d1_d3_pref_d_all(arc_red_bin3_d1d3_reliability_ind);

lacz_red_bin1_d1vd2_reliability_ind = find(lacz_red_d1vd2_reliability_all >=0 & lacz_red_d1vd2_reliability_all <=10);
lacz_red_bin2_d1vd2_reliability_ind = find(lacz_red_d1vd2_reliability_all >=11 & lacz_red_d1vd2_reliability_all <=30);
lacz_red_bin3_d1vd2_reliability_ind = find(lacz_red_d1vd2_reliability_all >=31 & lacz_red_d1vd2_reliability_all <=90);

lacz_red_d1_d2_pref_bin1_reliability = lacz_red_d1_d2_pref_d_all(lacz_red_bin1_d1vd2_reliability_ind);
lacz_red_d1_d2_pref_bin2_reliability = lacz_red_d1_d2_pref_d_all(lacz_red_bin2_d1vd2_reliability_ind);
lacz_red_d1_d2_pref_bin3_reliability = lacz_red_d1_d2_pref_d_all(lacz_red_bin3_d1vd2_reliability_ind);

lacz_red_bin1_d1d3_reliability_ind = find(lacz_red_d1vd3_reliability_all >=0 & lacz_red_d1vd3_reliability_all <=10);
lacz_red_bin2_d1d3_reliability_ind = find(lacz_red_d1vd3_reliability_all >=11 & lacz_red_d1vd3_reliability_all <=30);
lacz_red_bin3_d1d3_reliability_ind = find(lacz_red_d1vd3_reliability_all >=31 & lacz_red_d1vd3_reliability_all <=90);

lacz_red_d1_d3_pref_bin1_reliability = lacz_red_d1_d3_pref_d_all(lacz_red_bin1_d1d3_reliability_ind);
lacz_red_d1_d3_pref_bin2_reliability = lacz_red_d1_d3_pref_d_all(lacz_red_bin2_d1d3_reliability_ind);
lacz_red_d1_d3_pref_bin3_reliability = lacz_red_d1_d3_pref_d_all(lacz_red_bin3_d1d3_reliability_ind);


arc_all_bin1_d1vd2_reliability_ind = find(arc_all_d1vd2_reliability_all >=0 & arc_all_d1vd2_reliability_all <=10);
arc_all_bin2_d1vd2_reliability_ind = find(arc_all_d1vd2_reliability_all >=11 & arc_all_d1vd2_reliability_all <=30);
arc_all_bin3_d1vd2_reliability_ind = find(arc_all_d1vd2_reliability_all >=31 & arc_all_d1vd2_reliability_all <=90);

arc_all_d1_d2_pref_bin1_reliability = arc_all_d1_d2_pref_d_all(arc_all_bin1_d1vd2_reliability_ind);
arc_all_d1_d2_pref_bin2_reliability = arc_all_d1_d2_pref_d_all(arc_all_bin2_d1vd2_reliability_ind);
arc_all_d1_d2_pref_bin3_reliability = arc_all_d1_d2_pref_d_all(arc_all_bin3_d1vd2_reliability_ind);

arc_all_bin1_d1d3_reliability_ind = find(arc_all_d1vd3_reliability_all >=0 & arc_all_d1vd3_reliability_all <=10);
arc_all_bin2_d1d3_reliability_ind = find(arc_all_d1vd3_reliability_all >=11 & arc_all_d1vd3_reliability_all <=30);
arc_all_bin3_d1d3_reliability_ind = find(arc_all_d1vd3_reliability_all >=31 & arc_all_d1vd3_reliability_all <=90);

arc_all_d1_d3_pref_bin1_reliability = arc_all_d1_d3_pref_d_all(arc_all_bin1_d1d3_reliability_ind);
arc_all_d1_d3_pref_bin2_reliability = arc_all_d1_d3_pref_d_all(arc_all_bin2_d1d3_reliability_ind);
arc_all_d1_d3_pref_bin3_reliability = arc_all_d1_d3_pref_d_all(arc_all_bin3_d1d3_reliability_ind);

lacz_all_bin1_d1vd2_reliability_ind = find(lacz_all_d1vd2_reliability_all >=0 & lacz_all_d1vd2_reliability_all <=10);
lacz_all_bin2_d1vd2_reliability_ind = find(lacz_all_d1vd2_reliability_all >=11 & lacz_all_d1vd2_reliability_all <=30);
lacz_all_bin3_d1vd2_reliability_ind = find(lacz_all_d1vd2_reliability_all >=31 & lacz_all_d1vd2_reliability_all <=90);

lacz_all_d1_d2_pref_bin1_reliability = lacz_all_d1_d2_pref_d_all(lacz_all_bin1_d1vd2_reliability_ind);
lacz_all_d1_d2_pref_bin2_reliability = lacz_all_d1_d2_pref_d_all(lacz_all_bin2_d1vd2_reliability_ind);
lacz_all_d1_d2_pref_bin3_reliability = lacz_all_d1_d2_pref_d_all(lacz_all_bin3_d1vd2_reliability_ind);

lacz_all_bin1_d1d3_reliability_ind = find(lacz_all_d1vd3_reliability_all >=0 & lacz_all_d1vd3_reliability_all <=10);
lacz_all_bin2_d1d3_reliability_ind = find(lacz_all_d1vd3_reliability_all >=11 & lacz_all_d1vd3_reliability_all <=30);
lacz_all_bin3_d1d3_reliability_ind = find(lacz_all_d1vd3_reliability_all >=31 & lacz_all_d1vd3_reliability_all <=90);

lacz_all_d1_d3_pref_bin1_reliability = lacz_all_d1_d3_pref_d_all(lacz_all_bin1_d1d3_reliability_ind);
lacz_all_d1_d3_pref_bin2_reliability = lacz_all_d1_d3_pref_d_all(lacz_all_bin2_d1d3_reliability_ind);
lacz_all_d1_d3_pref_bin3_reliability = lacz_all_d1_d3_pref_d_all(lacz_all_bin3_d1d3_reliability_ind);


%%
save(fullfile(realfnout, ['_grace_reliabilitysplit_pref_data.mat']), 'arc_red_d1_d2_pref_bin1_reliability', 'arc_red_d1_d2_pref_bin2_reliability', 'arc_red_d1_d2_pref_bin3_reliability', ...
    'arc_red_d1_d3_pref_bin1_reliability', 'arc_red_d1_d3_pref_bin2_reliability', 'arc_red_d1_d3_pref_bin3_reliability', ...
    'arc_all_d1_d2_pref_bin1_reliability', 'arc_all_d1_d2_pref_bin2_reliability', 'arc_all_d1_d2_pref_bin3_reliability', ...
    'arc_all_d1_d3_pref_bin1_reliability', 'arc_all_d1_d3_pref_bin2_reliability', 'arc_all_d1_d3_pref_bin3_reliability', ...
    'lacz_red_d1_d2_pref_bin1_reliability', 'lacz_red_d1_d2_pref_bin2_reliability', 'lacz_red_d1_d2_pref_bin3_reliability', ...
    'lacz_red_d1_d3_pref_bin1_reliability', 'lacz_red_d1_d3_pref_bin2_reliability', 'lacz_red_d1_d3_pref_bin3_reliability', ...
    'lacz_all_d1_d2_pref_bin1_reliability', 'lacz_all_d1_d2_pref_bin2_reliability', 'lacz_all_d1_d2_pref_bin3_reliability', ...
    'lacz_all_d1_d3_pref_bin1_reliability', 'lacz_all_d1_d3_pref_bin2_reliability', 'lacz_all_d1_d3_pref_bin3_reliability')
