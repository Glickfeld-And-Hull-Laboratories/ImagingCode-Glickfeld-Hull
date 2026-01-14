%This code is used to analyze individual subject data from the dark to
%light housing 2p imaging experiment, so that these data are prepared for
%pooled data afterwards.
%%
%clear everything
clear all
clear all global
clc
close all

%% REMEMBER TO CHANGE THE OUT FOLDER BASED ON TUNING CRITERION - NORMAL FOR TUNED ON 1/3 SESSIONS and STRICT FOR TUNED ON ALL 3 ***
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_Arc_allVred_LG'; %folder to save files to
% newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_Arc_greenVred_STRICT'; %folder to save files to
dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
if ~exist(newfnout)
    mkdir(newfnout)
end

%%

%we start with loading the proper information based on the mouse's info in
%the experiment list

d1_list = [1:3:79]; 
d2_list = [2:3:80]; 
d3_list = [3:3:81];

for i = 1:length(d1_list)
    d1 = d1_list(i);
    d2 = d2_list(i);
    d3 = d3_list(i);

mouse = expt(d1).mouse; 
group = expt(d1).red_indicator{1};
ref_str_d1 = ['runs-',expt(d1).runs];
ref_str_d2 = ['runs-',expt(d2).runs]; 
ref_str_d3 = ['runs-',expt(d3).runs]; 
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date_d1 = expt(d1).date; 
date_d2 = expt(d2).date; 
date_d3 = expt(d3).date; 
img_folder_d1 = expt(d1).runs; 
img_folder_d2 = expt(d2).runs; 
img_folder_d3 = expt(d3).runs; 

fprintf([mouse ' ' date_d1 '\n'])

%load ori info for each day
d1_ori = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningInfo.mat']));
d2_ori = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningInfo.mat']));
d3_ori = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'oriTuningInfo.mat']));

d1_fits = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningAndFits.mat']));
d2_fits = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningAndFits.mat']));
d3_fits = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'oriTuningAndFits.mat']));


%load multiday data for day 2
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));
d3_matches = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));
d3_k_max = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'k_and_max_vals.mat']));

%r or gr cell info
d2_matches_tj = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'tj_matches.mat']));
d3_matches_tj = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'tj_matches.mat']));

%significantly responsive cells
d1_sigresp = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'trialData.mat']));
d2_sigresp = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'trialData.mat']));
d3_sigresp = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'trialData.mat']));


tuned_d1 = d1_ori.ind_theta90;
tuned_d2 = d2_ori.ind_theta90;
tuned_d3 = d3_ori.ind_theta90;

sig_resp_d1 = d1_sigresp.good_ind;
sig_resp_d2 = d2_sigresp.good_ind;
sig_resp_d3 = d3_sigresp.good_ind;

%find green cells that match from d1 to d2 and d1 to d3 
green_match_d2 = d2_matches_tj.green_match_ind; 
green_match_d3 = d3_matches_tj.green_match_ind; 

%find red cells that match from d1 to d2 and d1 to d3 
red_match_d2 = d2_matches_tj.red_match_ind; 
red_match_d3 = d3_matches_tj.red_match_ind; 

red_tuned_d1 = intersect(intersect(unique([red_match_d2 red_match_d2]),sig_resp_d1),tuned_d1);
red_match_tuned_d2 = intersect(intersect(red_match_d2,sig_resp_d1),tuned_d1);
red_match_tuned_d3 = intersect(intersect(red_match_d3,sig_resp_d1),tuned_d1);
red_match_tuned_alld = intersect(intersect(intersect(red_match_d2,red_match_d3),sig_resp_d1),tuned_d1);

all_tuned_d1 = intersect(sig_resp_d1,tuned_d1);
all_match_tuned_d2 = intersect(intersect([green_match_d2 red_match_d2],sig_resp_d1),tuned_d1);
all_match_tuned_d3 = intersect(intersect([green_match_d3 red_match_d3],sig_resp_d1),tuned_d1);

red_sig_d1 = intersect(unique([red_match_d2 red_match_d3]),sig_resp_d1);
red_match_sig_d2 = intersect(red_match_d2,sig_resp_d1);
red_match_sig_d3 = intersect(red_match_d3,sig_resp_d1);

all_sig_d1 = sig_resp_d1;
all_match_sig_d2 = intersect([green_match_d2 red_match_d2],sig_resp_d1);
all_match_sig_d3 = intersect([green_match_d3 red_match_d3],sig_resp_d1);

%fit reliability

all_d1_reliability = d1_fits.fitReliability(all_sig_d1);
all_d1vd2_reliability = d1_fits.fitReliability(all_match_sig_d2);
all_d1vd3_reliability = d1_fits.fitReliability(all_match_sig_d3);
all_d2_reliability = d2_fits.fitReliability(all_match_sig_d2);
all_d3_reliability = d3_fits.fitReliability(all_match_sig_d3);

red_d1_reliability = d1_fits.fitReliability(red_sig_d1);
red_d1vd2_reliability = d1_fits.fitReliability(red_match_sig_d2);
red_d1vd3_reliability = d1_fits.fitReliability(red_match_sig_d3);
red_d2_reliability = d2_fits.fitReliability(red_match_sig_d2);
red_d3_reliability = d3_fits.fitReliability(red_match_sig_d3);

%pref ori

d1_prefori = d1_ori.prefOri(1,:);
d2_prefori = d2_ori.prefOri(1,:);
d3_prefori = d3_ori.prefOri(1,:);

d_prefori_d1_d2 = abs(d1_prefori-d2_prefori);
d_prefori_d1_d3 = abs(d1_prefori-d3_prefori);
d_prefori_d1_d2(find(d_prefori_d1_d2>90)) = 180-d_prefori_d1_d2(find(d_prefori_d1_d2>90));
d_prefori_d1_d3(find(d_prefori_d1_d3>90)) = 180-d_prefori_d1_d3(find(d_prefori_d1_d3>90));

all_d_score_prefori_d1_d2 = d_prefori_d1_d2(all_match_sig_d2);
all_d_score_prefori_d1_d3 = d_prefori_d1_d3(all_match_sig_d3);
red_d_score_prefori_d1_d2 = d_prefori_d1_d2(red_match_sig_d2);
red_d_score_prefori_d1_d3 = d_prefori_d1_d3(red_match_sig_d3);

all_d_score_prefori_d1_d2_tuned = d_prefori_d1_d2(all_match_tuned_d2);
all_d_score_prefori_d1_d3_tuned = d_prefori_d1_d3(all_match_tuned_d3);
red_d_score_prefori_d1_d2_tuned = d_prefori_d1_d2(red_match_tuned_d2);
red_d_score_prefori_d1_d3_tuned = d_prefori_d1_d3(red_match_tuned_d3);

% d1 tuning curve
all_d1_respEaOri = d1_fits.avgResponseEaOri(all_sig_d1,:);
all_d1_respEaOri_tuned = d1_fits.avgResponseEaOri(all_tuned_d1,:);
red_d1_respEaOri = d1_fits.avgResponseEaOri(red_sig_d1,:);
red_d1_respEaOri_tuned = d1_fits.avgResponseEaOri(red_tuned_d1,:);
red_alld_respEaOri_tuned = [d1_fits.avgResponseEaOri(red_match_tuned_alld,:); d2_fits.avgResponseEaOri(red_match_tuned_alld,:); d3_fits.avgResponseEaOri(red_match_tuned_alld,:)];

all_d1_vonMises = squeeze(d1_fits.vonMisesFitAllCellsAllBoots(:,1,all_sig_d1));
all_d1_vonMises_tuned = squeeze(d1_fits.vonMisesFitAllCellsAllBoots(:,1,all_tuned_d1));
red_d1_vonMises = squeeze(d1_fits.vonMisesFitAllCellsAllBoots(:,1,red_sig_d1));
red_d1_vonMises_tuned = squeeze(d1_fits.vonMisesFitAllCellsAllBoots(:,1,red_tuned_d1));

% k and max
all_d1_k_tuned = d1_k_max.k1(1,all_tuned_d1);
all_d1vd2_k_tuned = d1_k_max.k1(1,all_match_tuned_d2);
all_d1vd3_k_tuned = d1_k_max.k1(1,all_match_tuned_d3);
all_d2_k_tuned = d2_k_max.k1(1,all_match_tuned_d2);
all_d3_k_tuned = d3_k_max.k1(1,all_match_tuned_d3);

red_d1_k_tuned = d1_k_max.k1(1,red_tuned_d1);
red_d1vd2_k_tuned = d1_k_max.k1(1,red_match_tuned_d2);
red_d1vd3_k_tuned = d1_k_max.k1(1,red_match_tuned_d3);
red_d2_k_tuned = d2_k_max.k1(1,red_match_tuned_d2);
red_d3_k_tuned = d3_k_max.k1(1,red_match_tuned_d3);

all_d1_k = d1_k_max.k1(1,all_sig_d1);
all_d1vd2_k = d1_k_max.k1(1,all_match_sig_d2);
all_d1vd3_k = d1_k_max.k1(1,all_match_sig_d3);
all_d2_k = d2_k_max.k1(1,all_match_sig_d2);
all_d3_k = d3_k_max.k1(1,all_match_sig_d3);

red_d1_k = d1_k_max.k1(1,red_sig_d1);
red_d1vd2_k = d1_k_max.k1(1,red_match_sig_d2);
red_d1vd3_k = d1_k_max.k1(1,red_match_sig_d3);
red_d2_k = d2_k_max.k1(1,red_match_sig_d2);
red_d3_k = d3_k_max.k1(1,red_match_sig_d3);

all_d1_max = d1_k_max.max_dfof(1,all_sig_d1);
all_d1_max_tuned = d1_k_max.max_dfof(1,all_tuned_d1);
all_d1vd2_max = d1_k_max.max_dfof(1,all_match_sig_d2);
all_d1vd3_max = d1_k_max.max_dfof(1,all_match_sig_d3);
all_d2_max = d2_k_max.max_dfof(1,all_match_sig_d2);
all_d3_max = d3_k_max.max_dfof(1,all_match_sig_d3);

red_d1_max = d1_k_max.max_dfof(1,red_sig_d1);
red_d1_max_tuned = d1_k_max.max_dfof(1,red_tuned_d1);
red_d1vd2_max = d1_k_max.max_dfof(1,red_match_sig_d2);
red_d1vd3_max = d1_k_max.max_dfof(1,red_match_sig_d3);
red_d2_max = d2_k_max.max_dfof(1,red_match_sig_d2);
red_d3_max = d3_k_max.max_dfof(1,red_match_sig_d3);


% saving
if str2double(expt(d1).img_day) == 1
%save all of the change scores from above for future use
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_d_scores.mat']), 'all_d_score_prefori_d1_d2', 'all_d_score_prefori_d1_d3', ...
        'red_d_score_prefori_d1_d2', 'red_d_score_prefori_d1_d3', 'all_d_score_prefori_d1_d2_tuned', 'all_d_score_prefori_d1_d3_tuned', ...
        'red_d_score_prefori_d1_d2_tuned', 'red_d_score_prefori_d1_d3_tuned')

%save pref ori info
    % save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_ori_changes.mat']), 'green_d1_prefori', 'green_d2_prefori', 'green_d3_prefori', ...
    %     'red_d1_prefori', 'red_d2_prefori', 'red_d3_prefori')

%save tuned and matched cell IDs
    % save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_id_matches.mat']), 'green_match_all', 'red_match_all')

%save k and max vals
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_k_max_changes.mat']), 'all_d1_k_tuned', 'all_d2_k_tuned', 'all_d3_k_tuned', ...
        'all_d1vd2_k_tuned', 'all_d1vd3_k_tuned', 'red_d1_k_tuned', 'red_d2_k_tuned', 'red_d3_k_tuned', ...
        'red_d1vd2_k_tuned', 'red_d1vd3_k_tuned','all_d1_k', 'all_d2_k', 'all_d3_k', ...
        'all_d1vd2_k', 'all_d1vd3_k', 'red_d1_k', 'red_d2_k', 'red_d3_k', ...
        'red_d1vd2_k', 'red_d1vd3_k','all_d1_max_tuned', 'all_d1_max', 'all_d2_max', 'all_d3_max', ...
        'all_d1vd2_max', 'all_d1vd3_max', 'red_d1_max_tuned', 'red_d1_max', 'red_d2_max', 'red_d3_max', ...
        'red_d1vd2_max', 'red_d1vd3_max',...
        'all_d1_respEaOri','all_d1_respEaOri_tuned','red_d1_respEaOri','red_d1_respEaOri_tuned','red_alld_respEaOri_tuned',...
        'all_d1_vonMises','all_d1_vonMises_tuned','red_d1_vonMises','red_d1_vonMises_tuned','red_match_tuned_alld')

%save fit reliabilities
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_fit_reliability.mat']), 'all_d1_reliability', 'all_d1vd2_reliability', 'all_d1vd3_reliability', ...
        'all_d2_reliability', 'all_d3_reliability', 'red_d1vd2_reliability', 'red_d1vd3_reliability',...
        'red_d1_reliability', 'red_d2_reliability', 'red_d3_reliability')

else
%save all of the change scores from above for future use
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_d_scores_2.mat']), 'all_d_score_prefori_d1_d2', 'all_d_score_prefori_d1_d3', ...
        'red_d_score_prefori_d1_d2', 'red_d_score_prefori_d1_d3', 'all_d_score_prefori_d1_d2_tuned', 'all_d_score_prefori_d1_d3_tuned', ...
        'red_d_score_prefori_d1_d2_tuned', 'red_d_score_prefori_d1_d3_tuned')

%save pref ori info
    % save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_ori_changes_2.mat']), 'green_d1_prefori', 'green_d2_prefori', 'green_d3_prefori', ...
    %     'red_d1_prefori', 'red_d2_prefori', 'red_d3_prefori')

%save tuned and matched cell IDs
    % save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_id_matches_2.mat']), 'green_match_all', 'red_match_all')

%save k and max vals
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_k_max_changes_2.mat']), 'all_d1_k_tuned', 'all_d2_k_tuned', 'all_d3_k_tuned', ...
        'all_d1vd2_k_tuned', 'all_d1vd3_k_tuned', 'red_d1_k_tuned', 'red_d2_k_tuned', 'red_d3_k_tuned', ...
        'red_d1vd2_k_tuned', 'red_d1vd3_k_tuned','all_d1_k', 'all_d2_k', 'all_d3_k', ...
        'all_d1vd2_k', 'all_d1vd3_k', 'red_d1_k', 'red_d2_k', 'red_d3_k', ...
        'red_d1vd2_k', 'red_d1vd3_k','all_d1_max_tuned', 'all_d1_max', 'all_d2_max', 'all_d3_max', ...
        'all_d1vd2_max', 'all_d1vd3_max', 'red_d1_max_tuned', 'red_d1_max', 'red_d2_max', 'red_d3_max', ...
        'red_d1vd2_max', 'red_d1vd3_max',...
        'all_d1_respEaOri','all_d1_respEaOri_tuned','red_d1_respEaOri','red_d1_respEaOri_tuned','red_alld_respEaOri_tuned',...
        'all_d1_vonMises','all_d1_vonMises_tuned','red_d1_vonMises','red_d1_vonMises_tuned','red_match_tuned_alld')

%save fit reliabilities
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_fit_reliability_2.mat']), 'all_d1_reliability', 'all_d1vd2_reliability', 'all_d1vd3_reliability', ...
        'all_d2_reliability', 'all_d3_reliability', 'red_d1vd2_reliability', 'red_d1vd3_reliability',...
        'red_d1_reliability', 'red_d2_reliability', 'red_d3_reliability')

    end
end
