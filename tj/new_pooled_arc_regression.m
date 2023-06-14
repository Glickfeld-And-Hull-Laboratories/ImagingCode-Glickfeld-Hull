%CLEAR EVERYTHINg
clear all
clear all global
clc
close all



%% *** REMEMBER TO CHANGE NEW AND REALFNOUT BASED ON LOOSE OR STRICT TUNING CRITERION
%LOAD DATA AND IDENTIFY FOLDERS TO SAVE TO
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_Arc_greenVred'; %folder to save files to
% newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_Arc_greenVred_STRICT'; %folder to save files to
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\regress_new_pooled_Arc_greenVred'; %folder to save files to
dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%
%% ARC PROMOTER MICE
%%
arc_green_d1_pref_all = [];
arc_green_d2_pref_all = [];
arc_green_d3_pref_all = [];
arc_green_d1_d2_pref_d_all = [];
arc_green_d1_d3_pref_d_all = [];
arc_green_d1_k_all = [];
arc_green_d2_k_all = [];
arc_green_d3_k_all = [];
arc_green_d1_max_all = [];
arc_green_d2_max_all = [];
arc_green_d3_max_all = [];
arc_green_d1_reliability_all = [];
arc_green_d2_reliability_all = [];
arc_green_d3_reliability_all = [];

arc_green_top_d1_d2_pref_d_all = [];
arc_green_top_d1_d3_pref_d_all = [];

arc_green_bot_d1_d2_pref_d_all = [];
arc_green_bot_d1_d3_pref_d_all = [];

arc_red_d1_pref_all = [];
arc_red_d2_pref_all = [];
arc_red_d3_pref_all = [];
arc_red_d1_d2_pref_d_all = [];
arc_red_d1_d3_pref_d_all = [];
arc_red_d1_k_all = [];
arc_red_d2_k_all = [];
arc_red_d3_k_all = [];
arc_red_d1_max_all = [];
arc_red_d2_max_all = [];
arc_red_d3_max_all = [];
arc_red_d1_reliability_all = [];
arc_red_d2_reliability_all = [];
arc_red_d3_reliability_all = [];


arc_red_top_d1_d2_pref_d_all = [];
arc_red_top_d1_d3_pref_d_all = [];

arc_red_bot_d1_d2_pref_d_all = [];
arc_red_bot_d1_d3_pref_d_all = [];


arc_prefori_scores = [];
arc_pref_dscores = [];
arc_k_max_scores = [];
arc_median_pref_dscores = [];
arc_reliability_scores = [];
arc_wide_sharp_pref_dscores = [];


list = [7 10 25 28];
for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};

    arc_prefori_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes_2.mat']));
        arc_prefori_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes_2']));
    else
        arc_prefori_ses2 = [];
    end
    arc_prefori_ses_all = [arc_prefori_ses1 arc_prefori_ses2];
    arc_prefori_scores = [arc_prefori_scores arc_prefori_ses_all];

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

for idata = 1:length(arc_prefori_scores)
%     arc_green_d1_pref = arc_prefori_ses_all.green_d1_prefori;
%     arc_green_d1_pref_all = [arc_green_d1_pref_all arc_green_d1_pref];
%     d2_pref = arc_prefori_ses_all.d2_prefori;
%     arc_green_d2_pref_all = [arc_green_d2_pref_all d2_pref];
%     d3_pref = arc_prefori_ses_all.d3_prefori;
%     arc_green_d3_pref_all = [arc_green_d3_pref_all d3_pref];
    

    arc_green_d1_d2_pref_d = arc_pref_dscores(idata).green_d_score_prefori_d1_d2;
    arc_green_d1_d2_pref_d_all = [arc_green_d1_d2_pref_d_all arc_green_d1_d2_pref_d];
    arc_green_d1_d3_pref_d = arc_pref_dscores(idata).green_d_score_prefori_d1_d3;
    arc_green_d1_d3_pref_d_all = [arc_green_d1_d3_pref_d_all arc_green_d1_d3_pref_d];

    arc_red_d1_d2_pref_d = arc_pref_dscores(idata).red_d_score_prefori_d1_d2;
    arc_red_d1_d2_pref_d_all = [arc_red_d1_d2_pref_d_all arc_red_d1_d2_pref_d];
    arc_red_d1_d3_pref_d = arc_pref_dscores(idata).red_d_score_prefori_d1_d3;
    arc_red_d1_d3_pref_d_all = [arc_red_d1_d3_pref_d_all arc_red_d1_d3_pref_d];

    arc_green_d1_k = arc_k_max_scores(idata).green_d1_k;
    arc_green_d1_k_all = [arc_green_d1_k_all arc_green_d1_k];
    arc_green_d2_k = arc_k_max_scores(idata).green_d2_k;
    arc_green_d2_k_all = [arc_green_d2_k_all arc_green_d2_k];
    arc_green_d3_k = arc_k_max_scores(idata).green_d3_k;
    arc_green_d3_k_all = [arc_green_d3_k_all arc_green_d3_k];

    arc_green_d1_max = arc_k_max_scores(idata).green_d1_max;
    arc_green_d1_max_all = [arc_green_d1_max_all arc_green_d1_max];
    arc_green_d2_max = arc_k_max_scores(idata).green_d2_max;
    arc_green_d2_max_all = [arc_green_d2_max_all arc_green_d2_max];
    arc_green_d3_max = arc_k_max_scores(idata).green_d3_max;
    arc_green_d3_max_all = [arc_green_d3_max_all arc_green_d3_max];

    arc_red_d1_k = arc_k_max_scores(idata).red_d1_k;
    arc_red_d1_k_all = [arc_red_d1_k_all arc_red_d1_k];
    arc_red_d2_k = arc_k_max_scores(idata).red_d2_k;
    arc_red_d2_k_all = [arc_red_d2_k_all arc_red_d2_k];
    arc_red_d3_k = arc_k_max_scores(idata).red_d3_k;
    arc_red_d3_k_all = [arc_red_d3_k_all arc_red_d3_k];

    arc_red_d1_max = arc_k_max_scores(idata).red_d1_max;
    arc_red_d1_max_all = [arc_red_d1_max_all arc_red_d1_max];
    arc_red_d2_max = arc_k_max_scores(idata).red_d2_max;
    arc_red_d2_max_all = [arc_red_d2_max_all arc_red_d2_max];
    arc_red_d3_max = arc_k_max_scores(idata).red_d3_max;
    arc_red_d3_max_all = [arc_red_d3_max_all arc_red_d3_max];
    

    arc_green_d1_reliability = arc_reliability_scores(idata).green_d1_reliability;
    arc_green_d1_reliability_all = [arc_green_d1_reliability_all arc_green_d1_reliability];
    arc_green_d2_reliability = arc_reliability_scores(idata).green_d2_reliability;
    arc_green_d2_reliability_all = [arc_green_d2_reliability_all arc_green_d2_reliability];
    arc_green_d3_reliability = arc_reliability_scores(idata).green_d3_reliability;
    arc_green_d3_reliability_all = [arc_green_d3_reliability_all arc_green_d3_reliability];

    arc_red_d1_reliability = arc_reliability_scores(idata).red_d1_reliability;
    arc_red_d1_reliability_all = [arc_red_d1_reliability_all arc_red_d1_reliability];
    arc_red_d2_reliability = arc_reliability_scores(idata).red_d2_reliability;
    arc_red_d2_reliability_all = [arc_red_d2_reliability_all arc_red_d2_reliability];
    arc_red_d3_reliability = arc_reliability_scores(idata).red_d3_reliability;
    arc_red_d3_reliability_all = [arc_red_d3_reliability_all arc_red_d3_reliability];    

  
end


%%
%% SAME AS ABOVE BUT LacZ
%%

lacz_green_d1_pref_all = [];
lacz_green_d2_pref_all = [];
lacz_green_d3_pref_all = [];
lacz_green_d1_d2_pref_d_all = [];
lacz_green_d1_d3_pref_d_all = [];
lacz_green_d1_k_all = [];
lacz_green_d2_k_all = [];
lacz_green_d3_k_all = [];
lacz_green_d1_max_all = [];
lacz_green_d2_max_all = [];
lacz_green_d3_max_all = [];
lacz_green_d1_reliability_all = [];
lacz_green_d2_reliability_all = [];
lacz_green_d3_reliability_all = [];

lacz_red_d1_pref_all = [];
lacz_red_d2_pref_all = [];
lacz_red_d3_pref_all = [];
lacz_red_d1_d2_pref_d_all = [];
lacz_red_d1_d3_pref_d_all = [];
lacz_red_d1_k_all = [];
lacz_red_d2_k_all = [];
lacz_red_d3_k_all = [];
lacz_red_d1_max_all = [];
lacz_red_d2_max_all = [];
lacz_red_d3_max_all = [];
lacz_red_d1_reliability_all = [];
lacz_red_d2_reliability_all = [];
lacz_red_d3_reliability_all = [];


lacz_prefori_scores = [];
lacz_pref_dscores = [];
lacz_k_max_scores = [];
lacz_median_pref_dscores = [];
lacz_reliability_scores = [];
lacz_wide_sharp_pref_dscores = [];


list = [1 4 37 40];
for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};

    lacz_prefori_ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes']));
    if exist(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes_2.mat']));
        lacz_prefori_ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes_2']));
    else
        lacz_prefori_ses2 = [];
    end
    lacz_prefori_ses_all = [lacz_prefori_ses1 lacz_prefori_ses2];
    lacz_prefori_scores = [lacz_prefori_scores lacz_prefori_ses_all];

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
%%


for idata = 1:length(lacz_prefori_scores)
%     lacz_green_d1_pref = lacz_prefori_ses_all.green_d1_prefori;
%     lacz_green_d1_pref_all = [lacz_green_d1_pref_all lacz_green_d1_pref];
%     d2_pref = lacz_prefori_ses_all.d2_prefori;
%     lacz_green_d2_pref_all = [lacz_green_d2_pref_all d2_pref];
%     d3_pref = lacz_prefori_ses_all.d3_prefori;
%     lacz_green_d3_pref_all = [lacz_green_d3_pref_all d3_pref];
    

    lacz_green_d1_d2_pref_d = lacz_pref_dscores(idata).green_d_score_prefori_d1_d2;
    lacz_green_d1_d2_pref_d_all = [lacz_green_d1_d2_pref_d_all lacz_green_d1_d2_pref_d];
    lacz_green_d1_d3_pref_d = lacz_pref_dscores(idata).green_d_score_prefori_d1_d3;
    lacz_green_d1_d3_pref_d_all = [lacz_green_d1_d3_pref_d_all lacz_green_d1_d3_pref_d];

    lacz_red_d1_d2_pref_d = lacz_pref_dscores(idata).red_d_score_prefori_d1_d2;
    lacz_red_d1_d2_pref_d_all = [lacz_red_d1_d2_pref_d_all lacz_red_d1_d2_pref_d];
    lacz_red_d1_d3_pref_d = lacz_pref_dscores(idata).red_d_score_prefori_d1_d3;
    lacz_red_d1_d3_pref_d_all = [lacz_red_d1_d3_pref_d_all lacz_red_d1_d3_pref_d];

    lacz_green_d1_k = lacz_k_max_scores(idata).green_d1_k;
    lacz_green_d1_k_all = [lacz_green_d1_k_all lacz_green_d1_k];
    lacz_green_d2_k = lacz_k_max_scores(idata).green_d2_k;
    lacz_green_d2_k_all = [lacz_green_d2_k_all lacz_green_d2_k];
    lacz_green_d3_k = lacz_k_max_scores(idata).green_d3_k;
    lacz_green_d3_k_all = [lacz_green_d3_k_all lacz_green_d3_k];

    lacz_green_d1_max = lacz_k_max_scores(idata).green_d1_max;
    lacz_green_d1_max_all = [lacz_green_d1_max_all lacz_green_d1_max];
    lacz_green_d2_max = lacz_k_max_scores(idata).green_d2_max;
    lacz_green_d2_max_all = [lacz_green_d2_max_all lacz_green_d2_max];
    lacz_green_d3_max = lacz_k_max_scores(idata).green_d3_max;
    lacz_green_d3_max_all = [lacz_green_d3_max_all lacz_green_d3_max];

    lacz_red_d1_k = lacz_k_max_scores(idata).red_d1_k;
    lacz_red_d1_k_all = [lacz_red_d1_k_all lacz_red_d1_k];
    lacz_red_d2_k = lacz_k_max_scores(idata).red_d2_k;
    lacz_red_d2_k_all = [lacz_red_d2_k_all lacz_red_d2_k];
    lacz_red_d3_k = lacz_k_max_scores(idata).red_d3_k;
    lacz_red_d3_k_all = [lacz_red_d3_k_all lacz_red_d3_k];

    lacz_red_d1_max = lacz_k_max_scores(idata).red_d1_max;
    lacz_red_d1_max_all = [lacz_red_d1_max_all lacz_red_d1_max];
    lacz_red_d2_max = lacz_k_max_scores(idata).red_d2_max;
    lacz_red_d2_max_all = [lacz_red_d2_max_all lacz_red_d2_max];
    lacz_red_d3_max = lacz_k_max_scores(idata).red_d3_max;
    lacz_red_d3_max_all = [lacz_red_d3_max_all lacz_red_d3_max];
    

    lacz_green_d1_reliability = lacz_reliability_scores(idata).green_d1_reliability;
    lacz_green_d1_reliability_all = [lacz_green_d1_reliability_all lacz_green_d1_reliability];
    lacz_green_d2_reliability = lacz_reliability_scores(idata).green_d2_reliability;
    lacz_green_d2_reliability_all = [lacz_green_d2_reliability_all lacz_green_d2_reliability];
    lacz_green_d3_reliability = lacz_reliability_scores(idata).green_d3_reliability;
    lacz_green_d3_reliability_all = [lacz_green_d3_reliability_all lacz_green_d3_reliability];

    lacz_red_d1_reliability = lacz_reliability_scores(idata).red_d1_reliability;
    lacz_red_d1_reliability_all = [lacz_red_d1_reliability_all lacz_red_d1_reliability];
    lacz_red_d2_reliability = lacz_reliability_scores(idata).red_d2_reliability;
    lacz_red_d2_reliability_all = [lacz_red_d2_reliability_all lacz_red_d2_reliability];
    lacz_red_d3_reliability = lacz_reliability_scores(idata).red_d3_reliability;
    lacz_red_d3_reliability_all = [lacz_red_d3_reliability_all lacz_red_d3_reliability];    

  
end


%%
%regression

%arc day 1 values full model
X_arc_red_full = [arc_red_d1_reliability_all; arc_red_d1_max_all; arc_red_d1_k_all]';
lm_arc_red_full = fitlm(X_arc_red_full,arc_red_d1_d3_pref_d_all')

%lacz day 1 values full model
X_lac_red_full = [lacz_red_d1_reliability_all; lacz_red_d1_max_all; lacz_red_d1_k_all]';
lm_lacz_red_full = fitlm(X_lac_red_full,lacz_red_d1_d3_pref_d_all')


%residuals
figure;
histogram(lm_arc_red_full.Residuals.Raw)
hold on
histogram(lm_lacz_red_full.Residuals.Raw)

%predictions
figure;
scatter(arc_red_d1_d3_pref_d_all, lm_arc_red_full.Fitted)
hold on
scatter(lacz_red_d1_d3_pref_d_all, lm_lacz_red_full.Fitted)


%%
%% no 30k
%%

%indexing 30k values out
arc_red_no30k_index = find(arc_red_d1_k_all<29);
lacz_red_no30k_index = find(lacz_red_d1_k_all<29);


X_arc_red_no30 = [arc_red_d1_reliability_all(arc_red_no30k_index); arc_red_d1_max_all(arc_red_no30k_index); arc_red_d1_k_all(arc_red_no30k_index)]';
lm_arc_red_no30k = fitlm(X_arc_red_no30, arc_red_d1_d3_pref_d_all(arc_red_no30k_index));

X_lacz_red_no30 = [lacz_red_d1_reliability_all(lacz_red_no30k_index); lacz_red_d1_max_all(lacz_red_no30k_index); lacz_red_d1_k_all(lacz_red_no30k_index)]';
lm_lacz_red_no30k = fitlm(X_lacz_red_no30, lacz_red_d1_d3_pref_d_all(lacz_red_no30k_index));


%residuals
figure;
histogram(lm_arc_red_no30k.Residuals.Raw)
hold on
histogram(lm_lacz_red_no30k.Residuals.Raw)

%predictions
figure;
scatter(arc_red_d1_d3_pref_d_all(arc_red_no30k_index), lm_arc_red_no30k.Fitted)
hold on
scatter(lacz_red_d1_d3_pref_d_all(lacz_red_no30k_index), lm_lacz_red_no30k.Fitted)

%%
figure;
m=cdfplot(arc_red_d1_d3_pref_d_all);
hold on
n = cdfplot(lacz_red_d1_d3_pref_d_all);

%%
figure;
m=cdfplot(arc_red_d1_d3_pref_d_all(arc_red_no30k_index));
hold on
n = cdfplot(lacz_red_d1_d3_pref_d_all(lacz_red_no30k_index));


%%
%using predictor averages across 3 sessions
avg_arc_red_reliability = mean([arc_red_d1_reliability_all; arc_red_d2_reliability_all; arc_red_d3_reliability_all]);
avg_arc_red_max = mean([arc_red_d1_max_all; arc_red_d2_max_all; arc_red_d3_max_all]);
avg_arc_red_k = mean([arc_red_d1_k_all; arc_red_d2_k_all; arc_red_d3_k_all]);

X_avg_arc_red = [avg_arc_red_reliability; avg_arc_red_max; avg_arc_red_k]';
lm_avg_arc_red = fitlm(X_avg_arc_red,arc_red_d1_d3_pref_d_all')


avg_lacz_red_reliability = mean([lacz_red_d1_reliability_all; lacz_red_d2_reliability_all; lacz_red_d3_reliability_all]);
avg_lacz_red_max = mean([lacz_red_d1_max_all; lacz_red_d2_max_all; lacz_red_d3_max_all]);
avg_lacz_red_k = mean([lacz_red_d1_k_all; lacz_red_d2_k_all; lacz_red_d3_k_all]);

X_avg_lacz_red = [avg_lacz_red_reliability; avg_lacz_red_max; avg_lacz_red_k]';
lm_avg_lacz_red = fitlm(X_avg_lacz_red,lacz_red_d1_d3_pref_d_all')

%%
%same as above but no 30k values
no30k_avg_arc_red_reliability = mean([arc_red_d1_reliability_all(arc_red_no30k_index); arc_red_d2_reliability_all(arc_red_no30k_index); arc_red_d3_reliability_all(arc_red_no30k_index)]);
no30k_avg_arc_red_max = mean([arc_red_d1_max_all(arc_red_no30k_index); arc_red_d2_max_all(arc_red_no30k_index); arc_red_d3_max_all(arc_red_no30k_index)]);
no30k_avg_arc_red_k = mean([arc_red_d1_k_all(arc_red_no30k_index); arc_red_d2_k_all(arc_red_no30k_index); arc_red_d3_k_all(arc_red_no30k_index)]);

no30k_X_avg_arc_red = [no30k_avg_arc_red_reliability; no30k_avg_arc_red_max; no30k_avg_arc_red_k]';
fitlm(no30k_X_avg_arc_red, arc_red_d1_d3_pref_d_all(arc_red_no30k_index)')

no30k_avg_lacz_red_reliability = mean([lacz_red_d1_reliability_all(lacz_red_no30k_index); lacz_red_d2_reliability_all(lacz_red_no30k_index); lacz_red_d3_reliability_all(lacz_red_no30k_index)]);
no30k_avg_lacz_red_max = mean([lacz_red_d1_max_all(lacz_red_no30k_index); lacz_red_d2_max_all(lacz_red_no30k_index); lacz_red_d3_max_all(lacz_red_no30k_index)]);
no30k_avg_lacz_red_k = mean([lacz_red_d1_k_all(lacz_red_no30k_index); lacz_red_d2_k_all(lacz_red_no30k_index); lacz_red_d3_k_all(lacz_red_no30k_index)]);

no30k_X_avg_lacz_red = [no30k_avg_lacz_red_reliability; no30k_avg_lacz_red_max; no30k_avg_lacz_red_k]';
fitlm(no30k_X_avg_lacz_red, lacz_red_d1_d3_pref_d_all(lacz_red_no30k_index)')

%%
%residuals

a = lm_arc_red_no30k.Residuals.Raw
b = lm_lacz_red_no30k.Residuals.Raw
aa = arc_red_d1_d3_pref_d_all(arc_red_no30k_index)'
ab = lm_arc_red_no30k.Fitted

figure;
scatter(aa,ab)
%%
%standardized models

%arc all
X_std_arc_all = [zscore(arc_red_d1_reliability_all); zscore(arc_red_d1_max_all); zscore(arc_red_d1_k_all)]';
lm_std_arc_all = fitlm(X_std_arc_all, zscore(arc_red_d1_d3_pref_d_all))

zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

%lacz all
X_std_lacz_all = [zscor_xnan(lacz_red_d1_reliability_all); zscore(lacz_red_d1_max_all); zscore(lacz_red_d1_k_all)]';
lm_std_lacz_all = fitlm(X_std_lacz_all, zscore(lacz_red_d1_d3_pref_d_all))

%%
arc_labels_red = cell(length(arc_red_d1_d2_pref_d_all),1)
arc_labels_red(:) = {'arc'}
arc_labels_red = arc_labels_red'

lacz_labels_red = cell(length(lacz_red_d1_d2_pref_d_all),1)
lacz_labels_red(:) = {'lacz'}
lacz_labels_red = lacz_labels_red'

all_labels_red = categorical([arc_labels_red lacz_labels_red]')
std_avg_predictors_arc_red = [zscor_xnan(avg_arc_red_reliability); zscor_xnan(avg_arc_red_max); zscor_xnan(avg_arc_red_k)]'
std_avg_predictors_lacz_red = [zscor_xnan(avg_lacz_red_reliability); zscor_xnan(avg_lacz_red_max); zscor_xnan(avg_lacz_red_k)]'

full_std_table_red = table(all_labels_red,zscor_xnan([avg_arc_red_reliability';avg_lacz_red_reliability']), zscor_xnan([avg_arc_red_max'; avg_lacz_red_max']), zscor_xnan([avg_arc_red_k'; avg_lacz_red_k']), ...
    zscor_xnan([arc_red_d1_d3_pref_d_all';lacz_red_d1_d3_pref_d_all']), 'VariableNames', {'group', 'fit_reliability', 'max_dff', 'k_val', 'prefori_change_d1_d3'})

%%

full_std_lm_red = fitlm(full_std_table_red, 'prefori_change_d1_d3~group+fit_reliability+max_dff+k_val')

full_std_lm_int_red = fitlm(full_std_table_red, 'prefori_change_d1_d3~group+fit_reliability+max_dff+k_val+group:fit_reliability')

merge_full_std_step_model = stepwiselm(full_std_table_red, 'interactions')

merge_full_std_constant_model = stepwiselm(full_std_table_red, 'constant', 'ResponseVar', 'prefori_change_d1_d3')

figure;
plotInteraction(full_std_lm_int_red,'group','fit_reliability','predictions')


%%
figure;
histogram(full_std_lm_red.Residuals.Raw)

figure;
scatter(full_std_lm_red.Variables.prefori_change_d1_d3, full_std_lm_red.Fitted)


%%
%% GREEN CELLS

%%
%using predictor averages across 3 sessions
avg_arc_green_reliability = mean([arc_green_d1_reliability_all; arc_green_d2_reliability_all; arc_green_d3_reliability_all]);
avg_arc_green_max = mean([arc_green_d1_max_all; arc_green_d2_max_all; arc_green_d3_max_all]);
avg_arc_green_k = mean([arc_green_d1_k_all; arc_green_d2_k_all; arc_green_d3_k_all]);

X_avg_arc_green = [avg_arc_green_reliability; avg_arc_green_max; avg_arc_green_k]';
lm_avg_arc_green = fitlm(X_avg_arc_green,arc_green_d1_d3_pref_d_all')


avg_lacz_green_reliability = mean([lacz_green_d1_reliability_all; lacz_green_d2_reliability_all; lacz_green_d3_reliability_all]);
avg_lacz_green_max = mean([lacz_green_d1_max_all; lacz_green_d2_max_all; lacz_green_d3_max_all]);
avg_lacz_green_k = mean([lacz_green_d1_k_all; lacz_green_d2_k_all; lacz_green_d3_k_all]);

X_avg_lacz_green = [avg_lacz_green_reliability; avg_lacz_green_max; avg_lacz_green_k]';
lm_avg_lacz_green = fitlm(X_avg_lacz_green,lacz_green_d1_d3_pref_d_all')

%%
%same as above but no 30k values

%indexing 30k values out
arc_green_no30k_index = find(arc_green_d1_k_all<29);
lacz_green_no30k_index = find(lacz_green_d1_k_all<29);

no30k_avg_arc_green_reliability = mean([arc_green_d1_reliability_all(arc_green_no30k_index); arc_green_d2_reliability_all(arc_green_no30k_index); arc_green_d3_reliability_all(arc_green_no30k_index)]);
no30k_avg_arc_green_max = mean([arc_green_d1_max_all(arc_green_no30k_index); arc_green_d2_max_all(arc_green_no30k_index); arc_green_d3_max_all(arc_green_no30k_index)]);
no30k_avg_arc_green_k = mean([arc_green_d1_k_all(arc_green_no30k_index); arc_green_d2_k_all(arc_green_no30k_index); arc_green_d3_k_all(arc_green_no30k_index)]);

no30k_X_avg_arc_green = [no30k_avg_arc_green_reliability; no30k_avg_arc_green_max; no30k_avg_arc_green_k]';
fitlm(no30k_X_avg_arc_green, arc_green_d1_d3_pref_d_all(arc_green_no30k_index)')

no30k_avg_lacz_green_reliability = mean([lacz_green_d1_reliability_all(lacz_green_no30k_index); lacz_green_d2_reliability_all(lacz_green_no30k_index); lacz_green_d3_reliability_all(lacz_green_no30k_index)]);
no30k_avg_lacz_green_max = mean([lacz_green_d1_max_all(lacz_green_no30k_index); lacz_green_d2_max_all(lacz_green_no30k_index); lacz_green_d3_max_all(lacz_green_no30k_index)]);
no30k_avg_lacz_green_k = mean([lacz_green_d1_k_all(lacz_green_no30k_index); lacz_green_d2_k_all(lacz_green_no30k_index); lacz_green_d3_k_all(lacz_green_no30k_index)]);

no30k_X_avg_lacz_green = [no30k_avg_lacz_green_reliability; no30k_avg_lacz_green_max; no30k_avg_lacz_green_k]';
fitlm(no30k_X_avg_lacz_green, lacz_green_d1_d3_pref_d_all(lacz_green_no30k_index)')


%%
arc_labels_green = cell(length(arc_green_d1_d2_pref_d_all),1)
arc_labels_green(:) = {'arc_green'}
arc_labels_green = arc_labels_green'

lacz_labels_green = cell(length(lacz_green_d1_d2_pref_d_all),1)
lacz_labels_green(:) = {'lacz_green'}
lacz_labels_green = lacz_labels_green'

all_labels_green = categorical([arc_labels_green lacz_labels_green]')
std_avg_predictors_arc_green = [zscor_xnan(avg_arc_green_reliability); zscor_xnan(avg_arc_green_max); zscor_xnan(avg_arc_green_k)]'
std_avg_predictors_lacz_green = [zscor_xnan(avg_lacz_green_reliability); zscor_xnan(avg_lacz_green_max); zscor_xnan(avg_lacz_green_k)]'

full_std_table_green = table(all_labels_green,zscor_xnan([avg_arc_green_reliability';avg_lacz_green_reliability']), zscor_xnan([avg_arc_green_max'; avg_lacz_green_max']), zscor_xnan([avg_arc_green_k'; avg_lacz_green_k']), ...
    zscor_xnan([arc_green_d1_d3_pref_d_all';lacz_green_d1_d3_pref_d_all']), 'VariableNames', {'group', 'fit_reliability', 'max_dff', 'k_val', 'prefori_change_d1_d3'})

full_std_lm_green = fitlm(full_std_table_green, 'prefori_change_d1_d3~group+fit_reliability+max_dff+k_val')

%%

full_std_lm_int_green = fitlm(full_std_table_green, 'prefori_change_d1_d3~group+fit_reliability+max_dff+k_val+group:fit_reliability')

merge_full_std_step_model_green = stepwiselm(full_std_table_green, 'interactions')

merge_full_std_constant_model_green = stepwiselm(full_std_table_green, 'constant', 'ResponseVar', 'prefori_change_d1_d3')

figure;
plotInteraction(full_std_lm_int_green,'group','fit_reliability','predictions')
%%
%red v green
%%

all_labels_arc_red_green = categorical([arc_labels_red arc_labels_green]')

full_std_table_arc_red_green = table(all_labels_arc_red_green,zscor_xnan([avg_arc_red_reliability';avg_arc_green_reliability']), zscor_xnan([avg_arc_red_max'; avg_arc_green_max']), zscor_xnan([avg_arc_red_k'; avg_arc_green_k']), ...
    zscor_xnan([arc_red_d1_d3_pref_d_all';arc_green_d1_d3_pref_d_all']), 'VariableNames', {'group', 'fit_reliability', 'max_dff', 'k_val', 'prefori_change_d1_d3'})

full_std_lm_arc_red_green = fitlm(full_std_table_arc_red_green, 'prefori_change_d1_d3~group+fit_reliability+max_dff+k_val')

full_std_lm_int_arc_red_green = fitlm(full_std_table_arc_red_green, 'prefori_change_d1_d3~group+fit_reliability+max_dff+k_val+group:fit_reliability')

merge_full_std_step_model_red_green = stepwiselm(full_std_table_arc_red_green, 'interactions')

figure;
plotInteraction(full_std_lm_int_arc_red_green,'group','fit_reliability','predictions')

figure;
plotInteraction(full_std_lm_int_arc_red_green,'fit_reliability','k_val','predictions')

%%
all_labels_lacz_red_green = categorical([lacz_labels_red lacz_labels_green]')


arc_group_labels = cell(length(all_labels_arc_red_green),1)
arc_group_labels(:) = {'arc'}
arc_group_labels = arc_group_labels'

lacz_group_labels = cell(length(all_labels_lacz_red_green),1)
lacz_group_labels(:) = {'lacz'}
lacz_group_labels = lacz_group_labels'

arc_red_labels = cell(length(arc_labels_red),1)
arc_red_labels(:) = {'red'}
arc_red_labels = arc_red_labels'

arc_green_labels = cell(length(arc_labels_green),1)
arc_green_labels(:) = {'green'}
arc_green_labels = arc_green_labels'

lacz_red_labels = cell(length(lacz_labels_red),1)
lacz_red_labels(:) = {'red'}
lacz_red_labels = lacz_red_labels'

lacz_green_labels = cell(length(lacz_labels_green),1)
lacz_green_labels(:) = {'green'}
lacz_green_labels = lacz_green_labels'


total_table = table([arc_group_labels'; lacz_group_labels'], [arc_red_labels'; arc_green_labels'; lacz_red_labels'; lacz_green_labels'], ...
    zscor_xnan([avg_arc_red_reliability';avg_arc_green_reliability';avg_lacz_red_reliability'; avg_lacz_green_reliability']), ...
    zscor_xnan([arc_red_d1_d3_pref_d_all';arc_green_d1_d3_pref_d_all';lacz_red_d1_d3_pref_d_all';lacz_green_d1_d3_pref_d_all']), ...
    'VariableNames', {'gRNA', 'cell_color', 'avg_fit_reliability', 'prefori_change_d1_d3'})


total_full_std_model = fitlm(total_table,'interactions','ResponseVar','prefori_change_d1_d3',...
    'PredictorVars',{'gRNA','cell_color','avg_fit_reliability'},...
    'CategoricalVar',{'gRNA','cell_color'})


total_std_step_model_red_green = stepwiselm(total_table, 'interactions')

figure;
plotInteraction(total_full_std_model,'gRNA','cell_color','predictions')


%%
