%This script was designed to be able to directly compare the data from TJ's and Grace's portions of
%the Arc imaging project after the data have been individually processed per mouse. 
%These comparisons are done with an emphasis on between-group - that is - comparing the KRAB
%expressing cells of each separate group of mice. No green cell analyses are done here, as Grace's
%cohort had very low expression of green-only cells.

%Prior processing was done with 'new_multi_day_redVgreen_update_tjw' to identify key variables
%per mouse for green and red cells that were matched across all 3 imaging sessions. Cells
%were also required to be well-fit with tuning curve on any one of the days.

% TJ's data are from KRAB mice injected w/ GCaMP7f in green and either Arc Prom
%gRNA or LacZ gRNA in red, such that green cells are Non-KRAB-expressing and red are
%KRAB-expressing. Grace's data are the same but with gRNA to the Arc Enh, and a separate cohort
%of LacZ injected-mice.


% *** DO THE RELIABILITY SPLIT FIGURES ***
%%
%clear everything
clear all
clear all global
clc
close all


%%
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_pooled_Arc_allVred_LG'; %folder to load files from
realfnoutGL = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_pooled_Arc_allVred_GLmice_LG'; %folder to load files from
compfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_pooled_Arc_Enh_Prom_comps\LGcomp'; %folder to save to
if ~exist(compfnout)
    mkdir(compfnout)
end
% dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
% eval(dataset); %load dataset

%%
%load data that were already pooled and analyzed by 'new_pooled_redVgreen_GLmice_LG'
% or new_pooled_redVgreen_LG, which consist of key vars for matched, tuned cells
%key vars are changes in pref ori, changes in pref ori after splitting cells based on median max
%responsivity, k and max values, and day1 k (d1) and max vals

tj_pref = load(fullfile(realfnout, ['_tj_pooled_pref_data']));
tj_kmax = load(fullfile(realfnout, ['_tj_pooled_kmax_data']));
tj_d1_kmax = load(fullfile(realfnout, ['_tj_pooled_kmax_day1_data']));
tj_reliability_split_pref = load(fullfile(realfnout, ['_tj_reliabilitysplit_pref_data.mat']));
grace_pref = load(fullfile(realfnoutGL, ['_grace_pooled_pref_data']));
grace_kmax = load(fullfile(realfnoutGL, ['_grace_pooled_kmax_data']));
grace_d1_kmax = load(fullfile(realfnoutGL, ['_grace_pooled_kmax_day1_data']));
grace_reliability_split_pref = load(fullfile(realfnoutGL, ['_grace_reliabilitysplit_pref_data.mat']));

%%
%this part is just extracting certain variables for plotting from the loaded data structures
%arc_prom is for arc Prom while arc_enh is for arc Enh. note that tj and grace have
%separate lacz cohorts that will be compared.
%k is sharpness
%max is max dF/F
%tuned is cells that are both responsive and tuned on D1
%respEaOri is one value for each orienation measured ncells x 8
%vonMises is fit of respEaOri necells x 181
%d1 comparisons

arc_prom_red_d1_k_all = tj_d1_kmax.arc_red_d1_k_all;
arc_prom_red_d1_k_tuned_all = tj_d1_kmax.arc_red_d1_k_tuned_all;
arc_prom_red_d1_max_all = tj_d1_kmax.arc_red_d1_max_all;
arc_prom_red_d1_max_tuned_all = tj_d1_kmax.arc_red_d1_max_tuned_all;
arc_prom_red_d1_reliability_all = tj_d1_kmax.arc_red_d1_reliability_all;
tj_lacz_red_d1_k_all = tj_d1_kmax.lacz_red_d1_k_all;
tj_lacz_red_d1_k_tuned_all = tj_d1_kmax.lacz_red_d1_k_tuned_all;
tj_lacz_red_d1_max_all = tj_d1_kmax.lacz_red_d1_max_all;
tj_lacz_red_d1_max_tuned_all = tj_d1_kmax.lacz_red_d1_max_tuned_all;
tj_lacz_red_d1_reliability_all = tj_d1_kmax.lacz_red_d1_reliability_all;

arc_enh_red_d1_k_all = grace_d1_kmax.arc_red_d1_k_all;
arc_enh_red_d1_k_tuned_all = grace_d1_kmax.arc_red_d1_k_tuned_all;
arc_enh_red_d1_max_all = grace_d1_kmax.arc_red_d1_max_all;
arc_enh_red_d1_max_tuned_all = grace_d1_kmax.arc_red_d1_max_tuned_all;
arc_enh_red_d1_reliability_all = grace_d1_kmax.arc_red_d1_reliability_all;
grace_lacz_red_d1_k_all = grace_d1_kmax.lacz_red_d1_k_all;
grace_lacz_red_d1_k_tuned_all = grace_d1_kmax.lacz_red_d1_k_tuned_all;
grace_lacz_red_d1_max_all = grace_d1_kmax.lacz_red_d1_max_all;
grace_lacz_red_d1_max_tuned_all = grace_d1_kmax.lacz_red_d1_max_tuned_all;
grace_lacz_red_d1_reliability_all = tj_d1_kmax.lacz_red_d1_reliability_all;

arc_enh_red_d1_respEaOri_all = grace_d1_kmax.arc_red_d1_respEaOri_all;
arc_enh_red_d1_respEaOri_tuned_all = grace_d1_kmax.arc_red_d1_respEaOri_tuned_all;
grace_lacz_red_d1_respEaOri_all = grace_d1_kmax.lacz_red_d1_respEaOri_all;
grace_lacz_red_d1_respEaOri_tuned_all = grace_d1_kmax.lacz_red_d1_respEaOri_tuned_all;
arc_prom_red_d1_respEaOri_all = tj_d1_kmax.arc_red_d1_respEaOri_all;
arc_prom_red_d1_respEaOri_tuned_all = tj_d1_kmax.arc_red_d1_respEaOri_tuned_all;
tj_lacz_red_d1_respEaOri_all = tj_d1_kmax.lacz_red_d1_respEaOri_all;
tj_lacz_red_d1_respEaOri_tuned_all = tj_d1_kmax.lacz_red_d1_respEaOri_tuned_all;

arc_enh_all_d1_respEaOri_all = grace_d1_kmax.arc_all_d1_respEaOri_all;
arc_enh_all_d1_respEaOri_tuned_all = grace_d1_kmax.arc_all_d1_respEaOri_tuned_all;
grace_lacz_all_d1_respEaOri_all = grace_d1_kmax.lacz_all_d1_respEaOri_all;
grace_lacz_all_d1_respEaOri_tuned_all = grace_d1_kmax.lacz_all_d1_respEaOri_tuned_all;
arc_prom_all_d1_respEaOri_all = tj_d1_kmax.arc_all_d1_respEaOri_all;
arc_prom_all_d1_respEaOri_tuned_all = tj_d1_kmax.arc_all_d1_respEaOri_tuned_all;
tj_lacz_all_d1_respEaOri_all = tj_d1_kmax.lacz_all_d1_respEaOri_all;
tj_lacz_all_d1_respEaOri_tuned_all = tj_d1_kmax.lacz_all_d1_respEaOri_tuned_all;

arc_enh_red_d1_vonMises_all = grace_d1_kmax.arc_red_d1_vonMises_all;
arc_enh_red_d1_vonMises_tuned_all = grace_d1_kmax.arc_red_d1_vonMises_tuned_all;
grace_lacz_red_d1_vonMises_all = grace_d1_kmax.lacz_red_d1_vonMises_all;
grace_lacz_red_d1_vonMises_tuned_all = grace_d1_kmax.lacz_red_d1_vonMises_tuned_all;
arc_prom_red_d1_vonMises_all = tj_d1_kmax.arc_red_d1_vonMises_all;
arc_prom_red_d1_vonMises_tuned_all = tj_d1_kmax.arc_red_d1_vonMises_tuned_all;
tj_lacz_red_d1_vonMises_all = tj_d1_kmax.lacz_red_d1_vonMises_all;
tj_lacz_red_d1_vonMises_tuned_all = tj_d1_kmax.lacz_red_d1_vonMises_tuned_all;

arc_enh_all_d1_vonMises_all = grace_d1_kmax.arc_all_d1_vonMises_all;
arc_enh_all_d1_vonMises_tuned_all = grace_d1_kmax.arc_all_d1_vonMises_tuned_all;
grace_lacz_all_d1_vonMises_all = grace_d1_kmax.lacz_all_d1_vonMises_all;
grace_lacz_all_d1_vonMises_tuned_all = grace_d1_kmax.lacz_all_d1_vonMises_tuned_all;
arc_prom_all_d1_vonMises_all = tj_d1_kmax.arc_all_d1_vonMises_all;
arc_prom_all_d1_vonMises_tuned_all = tj_d1_kmax.arc_all_d1_vonMises_tuned_all;
tj_lacz_all_d1_vonMises_all = tj_d1_kmax.lacz_all_d1_vonMises_all;
tj_lacz_all_d1_vonMises_tuned_all = tj_d1_kmax.lacz_all_d1_vonMises_tuned_all;

arc_prom_all_d1_k_all = tj_d1_kmax.arc_all_d1_k_all;
arc_prom_all_d1_k_tuned_all = tj_d1_kmax.arc_all_d1_k_tuned_all;
arc_prom_all_d1_max_all = tj_d1_kmax.arc_all_d1_max_all;
tj_lacz_all_d1_k_all = tj_d1_kmax.lacz_all_d1_k_all;
tj_lacz_all_d1_k_tuned_all = tj_d1_kmax.lacz_all_d1_k_tuned_all;
tj_lacz_all_d1_max_all = tj_d1_kmax.lacz_all_d1_max_all;

arc_enh_all_d1_k_all = grace_d1_kmax.arc_all_d1_k_all;
arc_enh_all_d1_k_tuned_all = grace_d1_kmax.arc_all_d1_k_tuned_all;
arc_enh_all_d1_max_all = grace_d1_kmax.arc_all_d1_max_all;
grace_lacz_all_d1_k_all = grace_d1_kmax.lacz_all_d1_k_all;
grace_lacz_all_d1_k_tuned_all = grace_d1_kmax.lacz_all_d1_k_tuned_all;
grace_lacz_all_d1_max_all = grace_d1_kmax.lacz_all_d1_max_all;

% d1 vs d2
% ori change
arc_prom_red_d1_d2_pref_d_all = tj_pref.arc_red_d1_d2_pref_d_all;
tj_lacz_red_d1_d2_pref_d_all = tj_pref.lacz_red_d1_d2_pref_d_all;
arc_enh_red_d1_d2_pref_d_all = grace_pref.arc_red_d1_d2_pref_d_all;
grace_lacz_red_d1_d2_pref_d_all = grace_pref.lacz_red_d1_d2_pref_d_all;

arc_prom_red_d1_d2_pref_d_tuned_all = tj_pref.arc_red_d1_d2_pref_d_tuned_all;
tj_lacz_red_d1_d2_pref_d_tuned_all = tj_pref.lacz_red_d1_d2_pref_d_tuned_all;
arc_enh_red_d1_d2_pref_d_tuned_all = grace_pref.arc_red_d1_d2_pref_d_tuned_all;
grace_lacz_red_d1_d2_pref_d_tuned_all = grace_pref.lacz_red_d1_d2_pref_d_tuned_all;

arc_prom_all_d1_d2_pref_d_all = tj_pref.arc_all_d1_d2_pref_d_all;
tj_lacz_all_d1_d2_pref_d_all = tj_pref.lacz_all_d1_d2_pref_d_all;
arc_enh_all_d1_d2_pref_d_all = grace_pref.arc_all_d1_d2_pref_d_all;
grace_lacz_all_d1_d2_pref_d_all = grace_pref.lacz_all_d1_d2_pref_d_all;

arc_prom_all_d1_d2_pref_d_tuned_all = tj_pref.arc_all_d1_d2_pref_d_tuned_all;
tj_lacz_all_d1_d2_pref_d_tuned_all = tj_pref.lacz_all_d1_d2_pref_d_tuned_all;
arc_enh_all_d1_d2_pref_d_tuned_all = grace_pref.arc_all_d1_d2_pref_d_tuned_all;
grace_lacz_all_d1_d2_pref_d_tuned_all = grace_pref.lacz_all_d1_d2_pref_d_tuned_all;

%max and k change
arc_prom_red_d1_d2_d_k_all = tj_kmax.arc_red_d1_d2_k_all;
arc_prom_red_d1_d2_d_k_tuned_all = tj_kmax.arc_red_d1_d2_k_tuned_all;
arc_prom_red_d1_d2_d_max_all = tj_kmax.arc_red_d1_d2_max_all;
tj_lacz_red_d1_d2_d_k_all = tj_kmax.lacz_red_d1_d2_k_all;
tj_lacz_red_d1_d2_d_k_tuned_all = tj_kmax.lacz_red_d1_d2_k_tuned_all;
tj_lacz_red_d1_d2_d_max_all = tj_kmax.lacz_red_d1_d2_max_all;

arc_enh_red_d1_d2_d_k_all = grace_kmax.arc_red_d1_d2_k_all;
arc_enh_red_d1_d2_d_k_tuned_all = grace_kmax.arc_red_d1_d2_k_tuned_all;
arc_enh_red_d1_d2_d_max_all = grace_kmax.arc_red_d1_d2_max_all;
grace_lacz_red_d1_d2_d_k_all = grace_kmax.lacz_red_d1_d2_k_all;
grace_lacz_red_d1_d2_d_k_tuned_all = grace_kmax.lacz_red_d1_d2_k_tuned_all;
grace_lacz_red_d1_d2_d_max_all = grace_kmax.lacz_red_d1_d2_max_all;

arc_prom_all_d1_d2_d_k_all = tj_kmax.arc_all_d1_d2_k_all;
arc_prom_all_d1_d2_d_k_tuned_all = tj_kmax.arc_all_d1_d2_k_tuned_all;
arc_prom_all_d1_d2_d_max_all = tj_kmax.arc_all_d1_d2_max_all;
tj_lacz_all_d1_d2_d_k_all = tj_kmax.lacz_all_d1_d2_k_all;
tj_lacz_all_d1_d2_d_k_tuned_all = tj_kmax.lacz_all_d1_d2_k_tuned_all;
tj_lacz_all_d1_d2_d_max_all = tj_kmax.lacz_all_d1_d2_max_all;

arc_enh_all_d1_d2_d_k_all = grace_kmax.arc_all_d1_d2_k_all;
arc_enh_all_d1_d2_d_k_tuned_all = grace_kmax.arc_all_d1_d2_k_tuned_all;
arc_enh_all_d1_d2_d_max_all = grace_kmax.arc_all_d1_d2_max_all;
grace_lacz_all_d1_d2_d_k_all = grace_kmax.lacz_all_d1_d2_k_all;
grace_lacz_all_d1_d2_d_k_tuned_all = grace_kmax.lacz_all_d1_d2_k_tuned_all;
grace_lacz_all_d1_d2_d_max_all = grace_kmax.lacz_all_d1_d2_max_all;

%reliability bins
arc_prom_red_d1_d2_pref_bin1 = tj_reliability_split_pref.arc_red_d1_d2_pref_bin1_reliability;
arc_prom_red_d1_d2_pref_bin2 = tj_reliability_split_pref.arc_red_d1_d2_pref_bin2_reliability;
arc_prom_red_d1_d2_pref_bin3 = tj_reliability_split_pref.arc_red_d1_d2_pref_bin3_reliability;

tj_lacz_red_d1_d2_pref_bin1 = tj_reliability_split_pref.lacz_red_d1_d2_pref_bin1_reliability;
tj_lacz_red_d1_d2_pref_bin2 = tj_reliability_split_pref.lacz_red_d1_d2_pref_bin2_reliability;
tj_lacz_red_d1_d2_pref_bin3 = tj_reliability_split_pref.lacz_red_d1_d2_pref_bin3_reliability;

arc_enh_red_d1_d2_pref_bin1 = grace_reliability_split_pref.arc_red_d1_d2_pref_bin1_reliability;
arc_enh_red_d1_d2_pref_bin2 = grace_reliability_split_pref.arc_red_d1_d2_pref_bin2_reliability;
arc_enh_red_d1_d2_pref_bin3 = grace_reliability_split_pref.arc_red_d1_d2_pref_bin3_reliability;

grace_lacz_red_d1_d2_pref_bin1 = grace_reliability_split_pref.lacz_red_d1_d2_pref_bin1_reliability;
grace_lacz_red_d1_d2_pref_bin2 = grace_reliability_split_pref.lacz_red_d1_d2_pref_bin2_reliability;
grace_lacz_red_d1_d2_pref_bin3 = grace_reliability_split_pref.lacz_red_d1_d2_pref_bin3_reliability;

arc_prom_all_d1_d2_pref_bin1 = tj_reliability_split_pref.arc_all_d1_d2_pref_bin1_reliability;
arc_prom_all_d1_d2_pref_bin2 = tj_reliability_split_pref.arc_all_d1_d2_pref_bin2_reliability;
arc_prom_all_d1_d2_pref_bin3 = tj_reliability_split_pref.arc_all_d1_d2_pref_bin3_reliability;

tj_lacz_all_d1_d2_pref_bin1 = tj_reliability_split_pref.lacz_all_d1_d2_pref_bin1_reliability;
tj_lacz_all_d1_d2_pref_bin2 = tj_reliability_split_pref.lacz_all_d1_d2_pref_bin2_reliability;
tj_lacz_all_d1_d2_pref_bin3 = tj_reliability_split_pref.lacz_all_d1_d2_pref_bin3_reliability;

arc_enh_all_d1_d2_pref_bin1 = grace_reliability_split_pref.arc_all_d1_d2_pref_bin1_reliability;
arc_enh_all_d1_d2_pref_bin2 = grace_reliability_split_pref.arc_all_d1_d2_pref_bin2_reliability;
arc_enh_all_d1_d2_pref_bin3 = grace_reliability_split_pref.arc_all_d1_d2_pref_bin3_reliability;

grace_lacz_all_d1_d2_pref_bin1 = grace_reliability_split_pref.lacz_all_d1_d2_pref_bin1_reliability;
grace_lacz_all_d1_d2_pref_bin2 = grace_reliability_split_pref.lacz_all_d1_d2_pref_bin2_reliability;
grace_lacz_all_d1_d2_pref_bin3 = grace_reliability_split_pref.lacz_all_d1_d2_pref_bin3_reliability;


% d1 vs d3
% ori change
arc_prom_red_d1_d3_pref_d_all = tj_pref.arc_red_d1_d3_pref_d_all;
tj_lacz_red_d1_d3_pref_d_all = tj_pref.lacz_red_d1_d3_pref_d_all;
arc_enh_red_d1_d3_pref_d_all = grace_pref.arc_red_d1_d3_pref_d_all;
grace_lacz_red_d1_d3_pref_d_all = grace_pref.lacz_red_d1_d3_pref_d_all;

arc_prom_red_d1_d3_pref_d_tuned_all = tj_pref.arc_red_d1_d3_pref_d_tuned_all;
tj_lacz_red_d1_d3_pref_d_tuned_all = tj_pref.lacz_red_d1_d3_pref_d_tuned_all;
arc_enh_red_d1_d3_pref_d_tuned_all = grace_pref.arc_red_d1_d3_pref_d_tuned_all;
grace_lacz_red_d1_d3_pref_d_tuned_all = grace_pref.lacz_red_d1_d3_pref_d_tuned_all;

arc_prom_all_d1_d3_pref_d_all = tj_pref.arc_all_d1_d3_pref_d_all;
tj_lacz_all_d1_d3_pref_d_all = tj_pref.lacz_all_d1_d3_pref_d_all;
arc_enh_all_d1_d3_pref_d_all = grace_pref.arc_all_d1_d3_pref_d_all;
grace_lacz_all_d1_d3_pref_d_all = grace_pref.lacz_all_d1_d3_pref_d_all;

arc_prom_all_d1_d3_pref_d_tuned_all = tj_pref.arc_all_d1_d3_pref_d_tuned_all;
tj_lacz_all_d1_d3_pref_d_tuned_all = tj_pref.lacz_all_d1_d3_pref_d_tuned_all;
arc_enh_all_d1_d3_pref_d_tuned_all = grace_pref.arc_all_d1_d3_pref_d_tuned_all;
grace_lacz_all_d1_d3_pref_d_tuned_all = grace_pref.lacz_all_d1_d3_pref_d_tuned_all;

%max and k change
arc_prom_red_d1_d3_d_k_all = tj_kmax.arc_red_d1_d3_k_all;
arc_prom_red_d1_d3_d_k_tuned_all = tj_kmax.arc_red_d1_d3_k_tuned_all;
arc_prom_red_d1_d3_d_max_all = tj_kmax.arc_red_d1_d3_max_all;
tj_lacz_red_d1_d3_d_k_all = tj_kmax.lacz_red_d1_d3_k_all;
tj_lacz_red_d1_d3_d_k_tuned_all = tj_kmax.lacz_red_d1_d3_k_tuned_all;
tj_lacz_red_d1_d3_d_max_all = tj_kmax.lacz_red_d1_d3_max_all;

arc_enh_red_d1_d3_d_k_all = grace_kmax.arc_red_d1_d3_k_all;
arc_enh_red_d1_d3_d_k_tuned_all = grace_kmax.arc_red_d1_d3_k_tuned_all;
arc_enh_red_d1_d3_d_max_all = grace_kmax.arc_red_d1_d3_max_all;
grace_lacz_red_d1_d3_d_k_all = grace_kmax.lacz_red_d1_d3_k_all;
grace_lacz_red_d1_d3_d_k_tuned_all = grace_kmax.lacz_red_d1_d3_k_tuned_all;
grace_lacz_red_d1_d3_d_max_all = grace_kmax.lacz_red_d1_d3_max_all;

arc_prom_all_d1_d3_d_k_all = tj_kmax.arc_all_d1_d3_k_all;
arc_prom_all_d1_d3_d_k_tuned_all = tj_kmax.arc_all_d1_d3_k_tuned_all;
arc_prom_all_d1_d3_d_max_all = tj_kmax.arc_all_d1_d3_max_all;
tj_lacz_all_d1_d3_d_k_all = tj_kmax.lacz_all_d1_d3_k_all;
tj_lacz_all_d1_d3_d_k_tuned_all = tj_kmax.lacz_all_d1_d3_k_tuned_all;
tj_lacz_all_d1_d3_d_max_all = tj_kmax.lacz_all_d1_d3_max_all;

arc_enh_all_d1_d3_d_k_all = grace_kmax.arc_all_d1_d3_k_all;
arc_enh_all_d1_d3_d_k_tuned_all = grace_kmax.arc_all_d1_d3_k_tuned_all;
arc_enh_all_d1_d3_d_max_all = grace_kmax.arc_all_d1_d3_max_all;
grace_lacz_all_d1_d3_d_k_all = grace_kmax.lacz_all_d1_d3_k_all;
grace_lacz_all_d1_d3_d_k_tuned_all = grace_kmax.lacz_all_d1_d3_k_tuned_all;
grace_lacz_all_d1_d3_d_max_all = grace_kmax.lacz_all_d1_d3_max_all;

%reliability bins
arc_prom_red_d1_d3_pref_bin1 = tj_reliability_split_pref.arc_red_d1_d3_pref_bin1_reliability;
arc_prom_red_d1_d3_pref_bin2 = tj_reliability_split_pref.arc_red_d1_d3_pref_bin2_reliability;
arc_prom_red_d1_d3_pref_bin3 = tj_reliability_split_pref.arc_red_d1_d3_pref_bin3_reliability;

tj_lacz_red_d1_d3_pref_bin1 = tj_reliability_split_pref.lacz_red_d1_d3_pref_bin1_reliability;
tj_lacz_red_d1_d3_pref_bin2 = tj_reliability_split_pref.lacz_red_d1_d3_pref_bin2_reliability;
tj_lacz_red_d1_d3_pref_bin3 = tj_reliability_split_pref.lacz_red_d1_d3_pref_bin3_reliability;

arc_enh_red_d1_d3_pref_bin1 = grace_reliability_split_pref.arc_red_d1_d3_pref_bin1_reliability;
arc_enh_red_d1_d3_pref_bin2 = grace_reliability_split_pref.arc_red_d1_d3_pref_bin2_reliability;
arc_enh_red_d1_d3_pref_bin3 = grace_reliability_split_pref.arc_red_d1_d3_pref_bin3_reliability;

grace_lacz_red_d1_d3_pref_bin1 = grace_reliability_split_pref.lacz_red_d1_d3_pref_bin1_reliability;
grace_lacz_red_d1_d3_pref_bin2 = grace_reliability_split_pref.lacz_red_d1_d3_pref_bin2_reliability;
grace_lacz_red_d1_d3_pref_bin3 = grace_reliability_split_pref.lacz_red_d1_d3_pref_bin3_reliability;

arc_prom_all_d1_d3_pref_bin1 = tj_reliability_split_pref.arc_all_d1_d3_pref_bin1_reliability;
arc_prom_all_d1_d3_pref_bin2 = tj_reliability_split_pref.arc_all_d1_d3_pref_bin2_reliability;
arc_prom_all_d1_d3_pref_bin3 = tj_reliability_split_pref.arc_all_d1_d3_pref_bin3_reliability;

tj_lacz_all_d1_d3_pref_bin1 = tj_reliability_split_pref.lacz_all_d1_d3_pref_bin1_reliability;
tj_lacz_all_d1_d3_pref_bin2 = tj_reliability_split_pref.lacz_all_d1_d3_pref_bin2_reliability;
tj_lacz_all_d1_d3_pref_bin3 = tj_reliability_split_pref.lacz_all_d1_d3_pref_bin3_reliability;

arc_enh_all_d1_d3_pref_bin1 = grace_reliability_split_pref.arc_all_d1_d3_pref_bin1_reliability;
arc_enh_all_d1_d3_pref_bin2 = grace_reliability_split_pref.arc_all_d1_d3_pref_bin2_reliability;
arc_enh_all_d1_d3_pref_bin3 = grace_reliability_split_pref.arc_all_d1_d3_pref_bin3_reliability;

grace_lacz_all_d1_d3_pref_bin1 = grace_reliability_split_pref.lacz_all_d1_d3_pref_bin1_reliability;
grace_lacz_all_d1_d3_pref_bin2 = grace_reliability_split_pref.lacz_all_d1_d3_pref_bin2_reliability;
grace_lacz_all_d1_d3_pref_bin3 = grace_reliability_split_pref.lacz_all_d1_d3_pref_bin3_reliability;


%%
%here we plot a comparison of how the same cells change their pref ori across both sessions1 and 2 as well
%as between 1 and 3. the comparison between groups is done on KRAB-expressing (red) cells only

fig = figure;
subplot(2,2,1)
h=cdfplot(arc_prom_red_d1_d2_pref_d_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(tj_lacz_red_d1_d2_pref_d_tuned_all);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_d2_pref_d_tuned_all);
set(l, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
m=cdfplot(grace_lacz_red_d1_d2_pref_d_tuned_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Arc Prom- n = ' num2str(size(arc_prom_red_d1_d2_pref_d_tuned_all,2))], ['LacZ- n = ' num2str(size(tj_lacz_red_d1_d2_pref_d_tuned_all,2))], ['Arc Enh- n = ' num2str(size(arc_enh_red_d1_d2_pref_d_tuned_all,2))], ['LacZ- n = ' num2str(size(grace_lacz_red_d1_d2_pref_d_tuned_all,2))], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_d2_pref_d_tuned_all,tj_lacz_red_d1_d2_pref_d_tuned_all);
[h_enh p_enh] = ttest2(arc_enh_red_d1_d2_pref_d_tuned_all,grace_lacz_red_d1_d2_pref_d_tuned_all);
title([{'Red- D1 v D2- Prom'} {['p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}])
hold off

subplot(2,2,2)
h=cdfplot(arc_prom_red_d1_d3_pref_d_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(tj_lacz_red_d1_d3_pref_d_tuned_all);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_d3_pref_d_tuned_all);
set(l, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
m=cdfplot(grace_lacz_red_d1_d3_pref_d_tuned_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Arc Prom- n = ' num2str(size(arc_prom_red_d1_d3_pref_d_tuned_all,2))], ['LacZ- n = ' num2str(size(tj_lacz_red_d1_d3_pref_d_tuned_all,2))], ['Arc Enh- n = ' num2str(size(arc_enh_red_d1_d3_pref_d_tuned_all,2))], ['LacZ- n = ' num2str(size(grace_lacz_red_d1_d3_pref_d_tuned_all,2))], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_d3_pref_d_tuned_all,tj_lacz_red_d1_d3_pref_d_tuned_all);
[h_enh p_enh] = ttest2(arc_enh_red_d1_d3_pref_d_tuned_all,grace_lacz_red_d1_d3_pref_d_tuned_all);
title([{'Red- D1 v D3- Prom'} {['p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}])
hold off

subplot(2,2,3)
h=cdfplot(arc_prom_all_d1_d2_pref_d_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(tj_lacz_all_d1_d2_pref_d_tuned_all);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
l=cdfplot(arc_enh_all_d1_d2_pref_d_tuned_all);
set(l, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
m=cdfplot(grace_lacz_all_d1_d2_pref_d_tuned_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Arc Prom- n = ' num2str(size(arc_prom_all_d1_d2_pref_d_tuned_all,2))], ['LacZ- n = ' num2str(size(tj_lacz_all_d1_d2_pref_d_tuned_all,2))], ['Arc Enh- n = ' num2str(size(arc_enh_all_d1_d2_pref_d_tuned_all,2))], ['LacZ- n = ' num2str(size(grace_lacz_all_d1_d2_pref_d_tuned_all,2))], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_all_d1_d2_pref_d_tuned_all,tj_lacz_all_d1_d2_pref_d_tuned_all);
[h_enh p_enh] = ttest2(arc_enh_all_d1_d2_pref_d_tuned_all,grace_lacz_all_d1_d2_pref_d_tuned_all);
title([{'All- D1 v D2- Prom'} {['p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}])
hold off

subplot(2,2,4)
h=cdfplot(arc_prom_all_d1_d3_pref_d_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(tj_lacz_all_d1_d3_pref_d_tuned_all);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
l=cdfplot(arc_enh_all_d1_d3_pref_d_tuned_all);
set(l, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
m=cdfplot(grace_lacz_all_d1_d3_pref_d_tuned_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Arc Prom- n = ' num2str(size(arc_prom_all_d1_d3_pref_d_tuned_all,2))], ['LacZ- n = ' num2str(size(tj_lacz_all_d1_d3_pref_d_tuned_all,2))], ['Arc Enh- n = ' num2str(size(arc_enh_all_d1_d3_pref_d_tuned_all,2))], ['LacZ- n = ' num2str(size(grace_lacz_all_d1_d3_pref_d_tuned_all,2))], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_all_d1_d3_pref_d_tuned_all,tj_lacz_all_d1_d3_pref_d_tuned_all);
[h_enh p_enh] = ttest2(arc_enh_all_d1_d3_pref_d_tuned_all,grace_lacz_all_d1_d3_pref_d_tuned_all);
title([{'All- D1 v D3- Prom'} {['p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}])
hold off

print(fullfile(compfnout, ['allmice_changepref.pdf']), '-dpdf', '-bestfit')

%% LacZ D1D2 vs D1D3 change in pref ori
fig = figure;
subplot(2,2,1)
h=cdfplot([tj_lacz_red_d1_d2_pref_d_tuned_all grace_lacz_red_d1_d2_pref_d_tuned_all]);
set(h, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
hold on
j=cdfplot([tj_lacz_red_d1_d3_pref_d_tuned_all grace_lacz_red_d1_d3_pref_d_tuned_all]);
set(j, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['D1 vs D2- n = ' num2str(length([tj_lacz_red_d1_d2_pref_d_tuned_all grace_lacz_red_d1_d2_pref_d_tuned_all]))], ...
    ['D1 vs D3- n = ' num2str(length([tj_lacz_red_d1_d3_pref_d_tuned_all grace_lacz_red_d1_d3_pref_d_tuned_all]))], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
[h_d2d3 p_d2d3] = ttest2([tj_lacz_red_d1_d2_pref_d_tuned_all grace_lacz_red_d1_d2_pref_d_tuned_all],[tj_lacz_red_d1_d3_pref_d_tuned_all grace_lacz_red_d1_d3_pref_d_tuned_all]);
title([{'Red Lacz'} {['p = ' num2str(chop(p_d2d3,2))]}])
hold off

subplot(2,2,2)
h=cdfplot([arc_prom_red_d1_d2_pref_d_tuned_all]);
set(h, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
hold on
j=cdfplot([arc_prom_red_d1_d3_pref_d_tuned_all]);
set(j, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
%legend(['D1 vs D2'], ['D1 vs D3'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
legend(['D1 vs D2- n = ' num2str(length(arc_prom_red_d1_d2_pref_d_tuned_all))], ...
    ['D1 vs D3- n = ' num2str(length(arc_prom_red_d1_d3_pref_d_tuned_all))], 'Location', 'southeast')
[h_d2d3 p_d2d3] = ttest2(arc_prom_red_d1_d2_pref_d_tuned_all,arc_prom_red_d1_d3_pref_d_tuned_all);
title([{'Red Prom'} {['p = ' num2str(chop(p_d2d3,2))]}])
hold off

subplot(2,2,3)
h=cdfplot([arc_enh_red_d1_d2_pref_d_tuned_all]);
set(h, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
hold on
j=cdfplot([arc_enh_red_d1_d3_pref_d_tuned_all]);
set(j, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['D1 vs D2- n = ' num2str(length(arc_enh_red_d1_d2_pref_d_tuned_all))], ...
    ['D1 vs D3- n = ' num2str(length(arc_enh_red_d1_d3_pref_d_tuned_all))], 'Location', 'southeast');
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
[h_d2d3 p_d2d3] = ttest2(arc_enh_red_d1_d2_pref_d_tuned_all,arc_enh_red_d1_d3_pref_d_tuned_all);
title([{'Red Enh'} {['p = ' num2str(chop(p_d2d3,2))]}])
hold off
sgtitle('Change in pref ori- D1D2 v D1D3')

print(fullfile(compfnout, ['allmice_changepref_d2vd3.pdf']), '-dpdf', '-bestfit')
%%
%this plot shows how much the value of k and max changed across all sessions

subplot(2,2,1)
h=cdfplot(arc_prom_red_d1_d2_d_k_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(tj_lacz_red_d1_d2_d_k_all);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_d2_d_k_all);
set(l, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
m=cdfplot(grace_lacz_red_d1_d2_d_k_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Arc Prom- n = ' num2str(size(arc_prom_red_d1_d2_d_k_all,2))], ['LacZ- n = ' num2str(size(tj_lacz_red_d1_d2_d_k_all,2))], ['Arc Enh- n = ' num2str(size(arc_enh_red_d1_d2_d_k_all,2))], ['LacZ- n = ' num2str(size(grace_lacz_red_d1_d2_d_k_all,2))], 'Location', 'southeast')
xlim([0 30])
xlabel(['Change in k Value'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_d2_d_k_all,tj_lacz_red_d1_d2_d_k_all);
[h_enh p_enh] = ttest2(arc_enh_red_d1_d2_d_k_all,grace_lacz_red_d1_d2_d_k_all);
title([{'Red- D1 v D2'} {['Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}])
hold off

subplot(2,2,2)
h=cdfplot(arc_prom_red_d1_d2_d_max_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(tj_lacz_red_d1_d2_d_max_all);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_d2_d_max_all);
set(l, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
m=cdfplot(grace_lacz_red_d1_d2_d_max_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
%legend(['Arc Prom (KRAB)'], ['LacZ (KRAB-TJ)'], ['Arc Enh (KRAB)'], ['LacZ (KRAB-Grace)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Change in Max dF/F Value'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_d2_d_max_all,tj_lacz_red_d1_d2_d_max_all);
[h_enh p_enh] = ttest2(arc_enh_red_d1_d2_d_max_all,grace_lacz_red_d1_d2_d_max_all);
title([{'Red- D1 v D2'} {['Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}])
hold off

subplot(2,2,3)
h=cdfplot(arc_prom_red_d1_d3_d_k_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(tj_lacz_red_d1_d3_d_k_all);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_d3_d_k_all);
set(l, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
m=cdfplot(grace_lacz_red_d1_d3_d_k_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Arc Prom- n = ' num2str(size(arc_prom_red_d1_d3_d_k_all,2))], ['LacZ- n = ' num2str(size(tj_lacz_red_d1_d3_d_k_all,2))], ['Arc Enh- n = ' num2str(size(arc_enh_red_d1_d3_d_k_all,2))], ['LacZ- n = ' num2str(size(grace_lacz_red_d1_d3_d_k_all,2))], 'Location', 'southeast')
xlim([0 30])
xlabel(['Change in k Value'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_d3_d_k_all,tj_lacz_red_d1_d3_d_k_all);
[h_enh p_enh] = ttest2(arc_enh_red_d1_d3_d_k_all,grace_lacz_red_d1_d3_d_k_all);
title([{'Red- D1 v D3'} {['Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}])
hold off

subplot(2,2,4)
h=cdfplot(arc_prom_red_d1_d3_d_max_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(tj_lacz_red_d1_d3_d_max_all);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_d3_d_max_all);
set(l, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
m=cdfplot(grace_lacz_red_d1_d3_d_max_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
%legend(['Arc Prom (KRAB)'], ['LacZ (KRAB-TJ)'], ['Arc Enh (KRAB)'], ['LacZ (KRAB-Grace)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Change in Max dF/F Value'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_d3_d_max_all,tj_lacz_red_d1_d3_d_max_all);
[h_enh p_enh] = ttest2(arc_enh_red_d1_d3_d_max_all,grace_lacz_red_d1_d3_d_max_all);
title([{'Red- D1 v D3'} {['Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}])
hold off

%print(fullfile(compfnout, ['allmice_changemaxall_red.pdf']), '-dpdf', '-bestfit')
sgtitle('D1 vs DN Change')
exportgraphics(gcf,fullfile(compfnout, ['allmice_change_red.pdf']))

all_prefori_diff = [tj_lacz_red_d1_d2_pref_d_tuned_all grace_lacz_red_d1_d2_pref_d_tuned_all tj_lacz_red_d1_d3_pref_d_tuned_all grace_lacz_red_d1_d3_pref_d_tuned_all arc_prom_red_d1_d2_pref_d_tuned_all arc_enh_red_d1_d2_pref_d_tuned_all arc_prom_red_d1_d3_pref_d_tuned_all arc_enh_red_d1_d3_pref_d_tuned_all];
all_prefori_diff_grna = [ones(size(tj_lacz_red_d1_d2_pref_d_tuned_all)) ones(size(grace_lacz_red_d1_d2_pref_d_tuned_all)) ones(size(tj_lacz_red_d1_d3_pref_d_tuned_all)) ones(size(grace_lacz_red_d1_d3_pref_d_tuned_all)) 2*ones(size(arc_prom_red_d1_d2_pref_d_tuned_all)) 3*ones(size(arc_enh_red_d1_d2_pref_d_tuned_all)) 2*ones(size(arc_prom_red_d1_d3_pref_d_tuned_all)) 3*ones(size(arc_enh_red_d1_d3_pref_d_tuned_all))];
all_prefori_diff_day = [ones(size(tj_lacz_red_d1_d2_pref_d_tuned_all)) ones(size(grace_lacz_red_d1_d2_pref_d_tuned_all)) 2*ones(size(tj_lacz_red_d1_d3_pref_d_tuned_all)) 2*ones(size(grace_lacz_red_d1_d3_pref_d_tuned_all)) ones(size(arc_prom_red_d1_d2_pref_d_tuned_all)) ones(size(arc_enh_red_d1_d2_pref_d_tuned_all)) 2*ones(size(arc_prom_red_d1_d3_pref_d_tuned_all)) 2*ones(size(arc_enh_red_d1_d3_pref_d_tuned_all))];
[a b c] = anovan(all_prefori_diff',[all_prefori_diff_grna;all_prefori_diff_day]');
[results, ~,~,gnames] = multcompare(c,"Dimension",[1 2]);

%%
%same as the other k value plot, but only for the distribution of values on day 1 (d1)

fig = figure;
subplot(2,2,1)
h=cdfplot(arc_prom_red_d1_k_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(tj_lacz_red_d1_k_all);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_k_all);
set(l, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
m=cdfplot(grace_lacz_red_d1_k_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Arc Prom- n = ' num2str(size(arc_prom_red_d1_k_all,2))], ['LacZ- n = ' num2str(size(tj_lacz_red_d1_k_all,2))], ['Arc Enh- n = ' num2str(size(arc_enh_red_d1_k_all,2))], ['LacZ- n = ' num2str(size(grace_lacz_red_d1_k_all,2))], 'Location', 'southeast')
xlim([0 30])
xlabel(['Day 1 k Value'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_k_all,tj_lacz_red_d1_k_all);
[h_enh p_enh] = ttest2(arc_enh_red_d1_k_all,grace_lacz_red_d1_k_all);
title([{'Red Responsive'},{['Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}]);
hold off

subplot(2,2,2)
h=cdfplot(arc_prom_red_d1_k_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(tj_lacz_red_d1_k_tuned_all);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_k_tuned_all);
set(l, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
m=cdfplot(grace_lacz_red_d1_k_tuned_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
legend(['Arc Prom- n = ' num2str(size(arc_prom_red_d1_k_tuned_all,2))], ['LacZ- n = ' num2str(size(tj_lacz_red_d1_k_tuned_all,2))], ['Arc Enh- n = ' num2str(size(arc_enh_red_d1_k_tuned_all,2))], ['LacZ- n = ' num2str(size(grace_lacz_red_d1_k_tuned_all,2))], 'Location', 'southeast')
xlim([0 30])
xlabel(['Day 1 k Value'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_k_tuned_all,tj_lacz_red_d1_k_tuned_all);
[h_enh p_enh] = ttest2(arc_enh_red_d1_k_tuned_all,grace_lacz_red_d1_k_tuned_all);
title([{'Red Tuned'},{['Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}]);
hold off


subplot(2,2,3)
h=cdfplot(arc_prom_red_d1_max_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(tj_lacz_red_d1_max_all);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_max_all);
set(l, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
m=cdfplot(grace_lacz_red_d1_max_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
%legend(['Arc Prom (KRAB)'], ['LacZ (KRAB-TJ)'], ['Arc Enh (KRAB)'], ['LacZ (KRAB-Grace)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Day 1 Max dF/F Value'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_max_all,tj_lacz_red_d1_max_all);
[h_enh p_enh] = ttest2(arc_enh_red_d1_max_all,grace_lacz_red_d1_max_all);
title([{'Red Responsive'},{['Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}]);
hold off

subplot(2,2,4)
h=cdfplot(arc_prom_red_d1_max_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold on
j=cdfplot(tj_lacz_red_d1_max_tuned_all);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_max_tuned_all);
set(l, 'LineStyle', '--', 'Color', 'r', 'LineWidth',2);
m=cdfplot(grace_lacz_red_d1_max_tuned_all);
set(m, 'LineStyle', '--', 'Color', 'k', 'LineWidth',2);
%legend(['Arc Prom (KRAB)'], ['LacZ (KRAB-TJ)'], ['Arc Enh (KRAB)'], ['LacZ (KRAB-Grace)'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Day 1 Max dF/F Value'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_max_tuned_all,tj_lacz_red_d1_max_tuned_all);
[h_enh p_enh] = ttest2(arc_enh_red_d1_max_tuned_all,grace_lacz_red_d1_max_tuned_all);
title([{'Red Tuned'},{['Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}]);
hold off

sgtitle('D1 props')
exportgraphics(gcf,fullfile(compfnout, ['allmice_day1.pdf']))

%% D1 tuning curves
% respEaOri
arc_prom_red_d1_respEaOri_all_sort = arc_prom_red_d1_respEaOri_all;
[v,m] = max(arc_prom_red_d1_respEaOri_all,[],2);
for i = 1:size(m,1)
    arc_prom_red_d1_respEaOri_all_sort(i,:) = circshift(arc_prom_red_d1_respEaOri_all(i,:),5-m(i))./v(i); 
end
arc_enh_red_d1_respEaOri_all_sort = arc_enh_red_d1_respEaOri_all;
[v,m] = max(arc_enh_red_d1_respEaOri_all,[],2);
for i = 1:size(m,1)
    arc_enh_red_d1_respEaOri_all_sort(i,:) = circshift(arc_enh_red_d1_respEaOri_all(i,:),5-m(i))./v(i); 
end
grace_lacz_red_d1_respEaOri_all_sort = grace_lacz_red_d1_respEaOri_all;
[v,m] = max(grace_lacz_red_d1_respEaOri_all,[],2);
for i = 1:size(m,1)
    grace_lacz_red_d1_respEaOri_all_sort(i,:) = circshift(grace_lacz_red_d1_respEaOri_all(i,:),5-m(i))./v(i); 
end
tj_lacz_red_d1_respEaOri_all_sort = tj_lacz_red_d1_respEaOri_all;
[v,m] = max(tj_lacz_red_d1_respEaOri_all,[],2);
for i = 1:size(m,1)
    tj_lacz_red_d1_respEaOri_all_sort(i,:) = circshift(tj_lacz_red_d1_respEaOri_all(i,:),5-m(i))./v(i); 
end
tj_lacz_red_d1_respEaOri_all_sort(find(tj_lacz_red_d1_respEaOri_all_sort(:,1)<-10),:) = [];
arc_prom_red_d1_respEaOri_tuned_all_sort = arc_prom_red_d1_respEaOri_tuned_all;
[v,m] = max(arc_prom_red_d1_respEaOri_tuned_all,[],2);
for i = 1:size(m,1)
    arc_prom_red_d1_respEaOri_tuned_all_sort(i,:) = circshift(arc_prom_red_d1_respEaOri_tuned_all(i,:),5-m(i))./v(i); 
end
arc_enh_red_d1_respEaOri_tuned_all_sort = arc_enh_red_d1_respEaOri_tuned_all;
[v,m] = max(arc_enh_red_d1_respEaOri_tuned_all,[],2);
for i = 1:size(m,1)
    arc_enh_red_d1_respEaOri_tuned_all_sort(i,:) = circshift(arc_enh_red_d1_respEaOri_tuned_all(i,:),5-m(i))./v(i); 
end
grace_lacz_red_d1_respEaOri_tuned_all_sort = grace_lacz_red_d1_respEaOri_tuned_all;
[v,m] = max(grace_lacz_red_d1_respEaOri_tuned_all,[],2);
for i = 1:size(m,1)
    grace_lacz_red_d1_respEaOri_tuned_all_sort(i,:) = circshift(grace_lacz_red_d1_respEaOri_tuned_all(i,:),5-m(i))./v(i); 
end
tj_lacz_red_d1_respEaOri_tuned_all_sort = tj_lacz_red_d1_respEaOri_tuned_all;
[v,m] = max(tj_lacz_red_d1_respEaOri_tuned_all,[],2);
for i = 1:size(m,1)
    tj_lacz_red_d1_respEaOri_tuned_all_sort(i,:) = circshift(tj_lacz_red_d1_respEaOri_tuned_all(i,:),5-m(i))./v(i); 
end


figure;
subplot(2,2,1)
errorbar(0:22.5:157.5, mean(arc_prom_red_d1_respEaOri_all_sort,1),std(arc_prom_red_d1_respEaOri_all_sort,[],1)./sqrt(size(arc_prom_red_d1_respEaOri_all_sort,1)),'-or')
hold on
errorbar(0:22.5:157.5, mean(tj_lacz_red_d1_respEaOri_all_sort,1),std(tj_lacz_red_d1_respEaOri_all_sort,[],1)./sqrt(size(tj_lacz_red_d1_respEaOri_all_sort,1)),'-ok')
ylim([-0.1 1.25])
ylabel('Normalized dF/F')
xlabel('Orientation (deg)')
title('Responsive')
subplot(2,2,2)
errorbar(0:22.5:157.5, mean(arc_enh_red_d1_respEaOri_all_sort,1),std(arc_enh_red_d1_respEaOri_all_sort,[],1)./sqrt(size(arc_enh_red_d1_respEaOri_all_sort,1)),'-or')
hold on
errorbar(0:22.5:157.5, mean(grace_lacz_red_d1_respEaOri_all_sort,1),std(grace_lacz_red_d1_respEaOri_all_sort,[],1)./sqrt(size(grace_lacz_red_d1_respEaOri_all_sort,1)),'-ok')
ylim([-0.1 1.25])
ylabel('Normalized dF/F')
xlabel('Orientation (deg)')
title('Responsive')
subplot(2,2,3)
errorbar(0:22.5:157.5, mean(arc_prom_red_d1_respEaOri_tuned_all_sort,1),std(arc_prom_red_d1_respEaOri_tuned_all_sort,[],1)./sqrt(size(arc_prom_red_d1_respEaOri_tuned_all_sort,1)),'-or')
hold on
errorbar(0:22.5:157.5, mean(tj_lacz_red_d1_respEaOri_tuned_all_sort,1),std(tj_lacz_red_d1_respEaOri_tuned_all_sort,[],1)./sqrt(size(tj_lacz_red_d1_respEaOri_tuned_all_sort,1)),'-ok')
ylim([-0.1 1.25])
ylabel('Normalized dF/F')
xlabel('Orientation (deg)')
title('Tuned')
subplot(2,2,4)
errorbar(0:22.5:157.5, mean(arc_enh_red_d1_respEaOri_tuned_all_sort,1),std(arc_enh_red_d1_respEaOri_tuned_all_sort,[],1)./sqrt(size(arc_enh_red_d1_respEaOri_tuned_all_sort,1)),'-or')
hold on
errorbar(0:22.5:157.5, mean(grace_lacz_red_d1_respEaOri_tuned_all_sort,1),std(grace_lacz_red_d1_respEaOri_tuned_all_sort,[],1)./sqrt(size(grace_lacz_red_d1_respEaOri_tuned_all_sort,1)),'-ok')
ylim([-0.1 1.25])
ylabel('Normalized dF/F')
xlabel('Orientation (deg)')
title('Tuned')

sgtitle('D1 props')
exportgraphics(gcf,fullfile(compfnout, ['allmice_day1_tuningcurve.pdf']))

% vonMises
arc_prom_red_d1_vonMises_all_sort = arc_prom_red_d1_vonMises_all;
[v,m] = max(arc_prom_red_d1_vonMises_all,[],2);
for i = 1:size(m,1)
    arc_prom_red_d1_vonMises_all_sort(i,:) = circshift(arc_prom_red_d1_vonMises_all(i,:),90-m(i))./v(i); 
end
arc_enh_red_d1_vonMises_all_sort = arc_enh_red_d1_vonMises_all;
[v,m] = max(arc_enh_red_d1_vonMises_all,[],2);
for i = 1:size(m,1)
    arc_enh_red_d1_vonMises_all_sort(i,:) = circshift(arc_enh_red_d1_vonMises_all(i,:),90-m(i))./v(i); 
end
grace_lacz_red_d1_vonMises_all_sort = grace_lacz_red_d1_vonMises_all;
[v,m] = max(grace_lacz_red_d1_vonMises_all,[],2);
for i = 1:size(m,1)
    grace_lacz_red_d1_vonMises_all_sort(i,:) = circshift(grace_lacz_red_d1_vonMises_all(i,:),90-m(i))./v(i); 
end
tj_lacz_red_d1_vonMises_all_sort = tj_lacz_red_d1_vonMises_all;
[v,m] = max(tj_lacz_red_d1_vonMises_all,[],2);
for i = 1:size(m,1)
    tj_lacz_red_d1_vonMises_all_sort(i,:) = circshift(tj_lacz_red_d1_vonMises_all(i,:),90-m(i))./v(i); 
end
tj_lacz_red_d1_vonMises_all_sort(find(tj_lacz_red_d1_vonMises_all_sort(:,1)<-10),:) = [];
arc_prom_red_d1_vonMises_tuned_all_sort = arc_prom_red_d1_vonMises_tuned_all;
[v,m] = max(arc_prom_red_d1_vonMises_tuned_all,[],2);
for i = 1:size(m,1)
    arc_prom_red_d1_vonMises_tuned_all_sort(i,:) = circshift(arc_prom_red_d1_vonMises_tuned_all(i,:),90-m(i))./v(i); 
end
arc_enh_red_d1_vonMises_tuned_all_sort = arc_enh_red_d1_vonMises_tuned_all;
[v,m] = max(arc_enh_red_d1_vonMises_tuned_all,[],2);
for i = 1:size(m,1)
    arc_enh_red_d1_vonMises_tuned_all_sort(i,:) = circshift(arc_enh_red_d1_vonMises_tuned_all(i,:),90-m(i))./v(i); 
end
grace_lacz_red_d1_vonMises_tuned_all_sort = grace_lacz_red_d1_vonMises_tuned_all;
[v,m] = max(grace_lacz_red_d1_vonMises_tuned_all,[],2);
for i = 1:size(m,1)
    grace_lacz_red_d1_vonMises_tuned_all_sort(i,:) = circshift(grace_lacz_red_d1_vonMises_tuned_all(i,:),90-m(i))./v(i); 
end
tj_lacz_red_d1_vonMises_tuned_all_sort = tj_lacz_red_d1_vonMises_tuned_all;
[v,m] = max(tj_lacz_red_d1_vonMises_tuned_all,[],2);
for i = 1:size(m,1)
    tj_lacz_red_d1_vonMises_tuned_all_sort(i,:) = circshift(tj_lacz_red_d1_vonMises_tuned_all(i,:),90-m(i))./v(i); 
end

figure;
subplot(2,2,1)
shadedErrorBar(0:180, mean(arc_prom_red_d1_vonMises_all_sort,1),std(arc_prom_red_d1_vonMises_all_sort,[],1)./sqrt(size(arc_prom_red_d1_vonMises_all_sort,1)),'-r');
hold on
shadedErrorBar(0:180, mean(tj_lacz_red_d1_vonMises_all_sort,1),std(tj_lacz_red_d1_vonMises_all_sort,[],1)./sqrt(size(tj_lacz_red_d1_vonMises_all_sort,1)),'-k');
ylim([-0.1 1.25])
ylabel('Normalized dF/F')
xlabel('Orientation (deg)')
title('Responsive')
subplot(2,2,2)
shadedErrorBar(0:180, mean(arc_enh_red_d1_vonMises_all_sort,1),std(arc_enh_red_d1_vonMises_all_sort,[],1)./sqrt(size(arc_enh_red_d1_vonMises_all_sort,1)),'-r');
hold on
shadedErrorBar(0:180, mean(grace_lacz_red_d1_vonMises_all_sort,1),std(grace_lacz_red_d1_vonMises_all_sort,[],1)./sqrt(size(grace_lacz_red_d1_vonMises_all_sort,1)),'-k');
ylim([-0.1 1.25])
ylabel('Normalized dF/F')
xlabel('Orientation (deg)')
title('Responsive')
subplot(2,2,3)
shadedErrorBar(0:180, mean(arc_prom_red_d1_vonMises_tuned_all_sort,1),std(arc_prom_red_d1_vonMises_tuned_all_sort,[],1)./sqrt(size(arc_prom_red_d1_vonMises_tuned_all_sort,1)),'-r');
hold on
shadedErrorBar(0:180, mean(tj_lacz_red_d1_vonMises_tuned_all_sort,1),std(tj_lacz_red_d1_vonMises_tuned_all_sort,[],1)./sqrt(size(tj_lacz_red_d1_vonMises_tuned_all_sort,1)),'-k');
ylim([-0.1 1.25])
ylabel('Normalized dF/F')
xlabel('Orientation (deg)')
title('Tuned')
subplot(2,2,4)
shadedErrorBar(0:180, mean(arc_enh_red_d1_vonMises_tuned_all_sort,1),std(arc_enh_red_d1_vonMises_tuned_all_sort,[],1)./sqrt(size(arc_enh_red_d1_vonMises_tuned_all_sort,1)),'-r');
hold on
shadedErrorBar(0:180, mean(grace_lacz_red_d1_vonMises_tuned_all_sort,1),std(grace_lacz_red_d1_vonMises_tuned_all_sort,[],1)./sqrt(size(grace_lacz_red_d1_vonMises_tuned_all_sort,1)),'-k');
ylim([-0.1 1.25])
ylabel('Normalized dF/F')
xlabel('Orientation (deg)')
title('Tuned')

sgtitle('D1 props')
exportgraphics(gcf,fullfile(compfnout, ['allmice_day1_tuningfit.pdf']))
%%
bin1_prom = [arc_prom_red_d1_d2_pref_bin1' ones(size(arc_prom_red_d1_d2_pref_bin1')) ones(size(arc_prom_red_d1_d2_pref_bin1'))];
bin2_prom = [arc_prom_red_d1_d2_pref_bin2' 2.*ones(size(arc_prom_red_d1_d2_pref_bin2')) ones(size(arc_prom_red_d1_d2_pref_bin2'))];
bin3_prom = [arc_prom_red_d1_d2_pref_bin3' 3.*ones(size(arc_prom_red_d1_d2_pref_bin3')) ones(size(arc_prom_red_d1_d2_pref_bin3'))];

bin1_lacz = [tj_lacz_red_d1_d2_pref_bin1' ones(size(tj_lacz_red_d1_d2_pref_bin1')) 2.*ones(size(tj_lacz_red_d1_d2_pref_bin1'))];
bin2_lacz = [tj_lacz_red_d1_d2_pref_bin2' 2.*ones(size(tj_lacz_red_d1_d2_pref_bin2')) 2.*ones(size(tj_lacz_red_d1_d2_pref_bin2'))];
bin3_lacz = [tj_lacz_red_d1_d2_pref_bin3' 3.*ones(size(tj_lacz_red_d1_d2_pref_bin3')) 2.*ones(size(tj_lacz_red_d1_d2_pref_bin3'))];

promVlacz = [bin1_prom;bin2_prom;bin3_prom;bin1_lacz;bin2_lacz;bin3_lacz];
[p,tbl,stats] = anovan(promVlacz(:,1),promVlacz(:,2:3),'model','interaction');

figure;
subplot(2,2,1)
a = errorbar([mean(arc_prom_red_d1_d2_pref_bin1), mean(arc_prom_red_d1_d2_pref_bin2), mean(arc_prom_red_d1_d2_pref_bin3)], ...
    [std(arc_prom_red_d1_d2_pref_bin1)/sqrt(length(arc_prom_red_d1_d2_pref_bin1)), ... 
    std(arc_prom_red_d1_d2_pref_bin2)/sqrt(length(arc_prom_red_d1_d2_pref_bin2)), ... 
    std(arc_prom_red_d1_d2_pref_bin3)/sqrt(length(arc_prom_red_d1_d2_pref_bin3))]);

hold on

b = errorbar([mean(tj_lacz_red_d1_d2_pref_bin1), mean(tj_lacz_red_d1_d2_pref_bin2), mean(tj_lacz_red_d1_d2_pref_bin3)], ...
    [std(tj_lacz_red_d1_d2_pref_bin1)/sqrt(length(tj_lacz_red_d1_d2_pref_bin1)), ... 
    std(tj_lacz_red_d1_d2_pref_bin2)/sqrt(length(tj_lacz_red_d1_d2_pref_bin2)), ... 
    std(tj_lacz_red_d1_d2_pref_bin3)/sqrt(length(tj_lacz_red_d1_d2_pref_bin3))]);

hold off

a.Marker = 'o';
a.MarkerFaceColor = [1 0 0];
a.Color = [1 0 0];
a.LineStyle = 'none';
a.MarkerSize = 5;

b.Marker = 'o';
b.MarkerFaceColor = [0 0 0];
b.Color = [0 0 0];
b.LineStyle = 'none';
b.MarkerSize = 5;

xlim([0.5 3.5])
ylim([0 90])
xticklabels({'0-10', '11-30','31-90'})
xlabel('Fit Reliability (degrees)')
ylabel('Mean Pref Ori Change (+- s.e.m.)')
%legend(['Arc Prom (KRAB)'], ['TJ-LacZ (KRAB)'])
title([{['Reliability- p = ' num2str(chop(p(1),2))]} {['Prom- p = ' num2str(chop(p(2),2))]} {['Interaction- p =' num2str(chop(p(3),2))]}])

%print(fullfile(compfnout, ['tjmice_reliabilitysplit_prefori_change_red.pdf']), '-dpdf', '-bestfit')

%%
bin1_enh = [arc_enh_red_d1_d2_pref_bin1' ones(size(arc_enh_red_d1_d2_pref_bin1')) ones(size(arc_enh_red_d1_d2_pref_bin1'))];
bin2_enh = [arc_enh_red_d1_d2_pref_bin2' 2.*ones(size(arc_enh_red_d1_d2_pref_bin2')) ones(size(arc_enh_red_d1_d2_pref_bin2'))];
bin3_enh = [arc_enh_red_d1_d2_pref_bin3' 3.*ones(size(arc_enh_red_d1_d2_pref_bin3')) ones(size(arc_enh_red_d1_d2_pref_bin3'))];

bin1_lacz = [grace_lacz_red_d1_d2_pref_bin1' ones(size(grace_lacz_red_d1_d2_pref_bin1')) 2.*ones(size(grace_lacz_red_d1_d2_pref_bin1'))];
bin2_lacz = [grace_lacz_red_d1_d2_pref_bin2' 2.*ones(size(grace_lacz_red_d1_d2_pref_bin2')) 2.*ones(size(grace_lacz_red_d1_d2_pref_bin2'))];
bin3_lacz = [grace_lacz_red_d1_d2_pref_bin3' 3.*ones(size(grace_lacz_red_d1_d2_pref_bin3')) 2.*ones(size(grace_lacz_red_d1_d2_pref_bin3'))];

enhVlacz = [bin1_enh;bin2_enh;bin3_enh;bin1_lacz;bin2_lacz;bin3_lacz];
[p,tbl,stats] = anovan(enhVlacz(:,1),enhVlacz(:,2:3),'model','interaction');

subplot(2,2,2)
c = errorbar([mean(arc_enh_red_d1_d2_pref_bin1), mean(arc_enh_red_d1_d2_pref_bin2), mean(arc_enh_red_d1_d2_pref_bin3)], ...
    [std(arc_enh_red_d1_d2_pref_bin1)/sqrt(length(arc_enh_red_d1_d2_pref_bin1)), ... 
    std(arc_enh_red_d1_d2_pref_bin2)/sqrt(length(arc_enh_red_d1_d2_pref_bin2)), ... 
    std(arc_enh_red_d1_d2_pref_bin3)/sqrt(length(arc_enh_red_d1_d2_pref_bin3))]);

hold on

d = errorbar([mean(grace_lacz_red_d1_d2_pref_bin1), mean(grace_lacz_red_d1_d2_pref_bin2), mean(grace_lacz_red_d1_d2_pref_bin3)], ...
    [std(grace_lacz_red_d1_d2_pref_bin1)/sqrt(length(grace_lacz_red_d1_d2_pref_bin1)), ... 
    std(grace_lacz_red_d1_d2_pref_bin2)/sqrt(length(grace_lacz_red_d1_d2_pref_bin2)), ... 
    std(grace_lacz_red_d1_d2_pref_bin3)/sqrt(length(grace_lacz_red_d1_d2_pref_bin3))]);

hold off

c.Marker = 'o';
c.MarkerFaceColor = [.75 0 0];
c.Color = [.75 0 0];
c.LineStyle = 'none';
c.MarkerSize = 5;

d.Marker = 'o';
d.MarkerFaceColor = [.25 .25 .25];
d.Color = [.25 .25 .25];
d.LineStyle = 'none';
d.MarkerSize = 5;

xlim([0.5 3.5])
ylim([0 90])
xticklabels({'0-10','11-30','31-90'})
xlabel('Fit Reliability (degrees)')
ylabel('Mean Pref Ori Change (+- s.e.m.)')
%legend(['Arc Enh (KRAB)'], ['Grace-LacZ (KRAB)'])

title([{['Reliability- p = ' num2str(chop(p(1),2))]} {['Enh- p = ' num2str(chop(p(2),2))]} {['Interaction- p =' num2str(chop(p(3),2))]}])
%print(fullfile(compfnout, ['gracemice_reliabilitysplit_prefori_change_red.pdf']), '-dpdf', '-bestfit')
exportgraphics(gcf,fullfile(compfnout, ['reliabilitysplit_prefori_change_red.pdf']))

%% combine lacz into single dataset
%change pref ori
fig = figure;
subplot(2,2,1)
j=cdfplot([tj_lacz_red_d1_d2_pref_d_tuned_all grace_lacz_red_d1_d2_pref_d_tuned_all]);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
hold on
h=cdfplot(arc_prom_red_d1_d2_pref_d_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_d2_pref_d_tuned_all);
set(l, 'LineStyle', '-', 'Color', 'm', 'LineWidth',2);
legend(['LacZ- n = ' num2str(length([tj_lacz_red_d1_d2_pref_d_tuned_all grace_lacz_red_d1_d2_pref_d_tuned_all]))],...
    ['Prom- n = ' num2str(length(arc_prom_red_d1_d2_pref_d_tuned_all))],...
    ['Enh- n = ' num2str(length(arc_enh_red_d1_d2_pref_d_tuned_all))], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_d2_pref_d_tuned_all,[grace_lacz_red_d1_d2_pref_d_tuned_all tj_lacz_red_d1_d2_pref_d_tuned_all]);
[h_enh p_enh] = ttest2(arc_enh_red_d1_d2_pref_d_tuned_all,[grace_lacz_red_d1_d2_pref_d_tuned_all tj_lacz_red_d1_d2_pref_d_tuned_all]);
title([{'Red- D1 v D2'} {[' Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}])
hold off

subplot(2,2,2)
j=cdfplot([tj_lacz_red_d1_d3_pref_d_tuned_all grace_lacz_red_d1_d3_pref_d_tuned_all]);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
hold on
h=cdfplot(arc_prom_red_d1_d3_pref_d_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_d3_pref_d_tuned_all);
set(l, 'LineStyle', '-', 'Color', 'm', 'LineWidth',2);
legend(['LacZ'],['Prom'],['Enh'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
legend(['LacZ- n = ' num2str(length([tj_lacz_red_d1_d3_pref_d_tuned_all grace_lacz_red_d1_d3_pref_d_tuned_all]))],...
    ['Prom- n = ' num2str(length(arc_prom_red_d1_d3_pref_d_tuned_all))],...
    ['Enh- n = ' num2str(length(arc_enh_red_d1_d3_pref_d_tuned_all))], 'Location', 'southeast')
[h_prom p_prom] = ttest2(arc_prom_red_d1_d3_pref_d_tuned_all,[grace_lacz_red_d1_d3_pref_d_tuned_all tj_lacz_red_d1_d3_pref_d_tuned_all]);
[h_enh p_enh] = ttest2(arc_enh_red_d1_d3_pref_d_tuned_all,[grace_lacz_red_d1_d3_pref_d_tuned_all tj_lacz_red_d1_d3_pref_d_tuned_all]);
title([{'Red- D1 v D3'} {[' Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}])
hold off

print(fullfile(compfnout, ['allmice_changepref_combLacZ.pdf']), '-dpdf', '-bestfit')

%% change k and max
figure;
subplot(2,2,1)
j=cdfplot([tj_lacz_red_d1_d2_d_k_all grace_lacz_red_d1_d2_d_k_all]);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
hold on
h=cdfplot(arc_prom_red_d1_d2_d_k_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_d2_d_k_all);
set(l, 'LineStyle', '-', 'Color', 'm', 'LineWidth',2);
legend(['LacZ- n = ' num2str(length([tj_lacz_red_d1_d2_d_k_all grace_lacz_red_d1_d2_d_k_all]))],...
    ['Prom- n = ' num2str(length(arc_prom_red_d1_d2_d_k_all))],...
    ['Enh- n = ' num2str(length(arc_enh_red_d1_d2_d_k_all))], 'Location', 'southeast')
xlim([0 30])
xlabel(['Change in k Value'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_d2_d_k_all,[tj_lacz_red_d1_d2_d_k_all grace_lacz_red_d1_d2_d_k_all]);
[h_enh p_enh] = ttest2(arc_enh_red_d1_d2_d_k_all,[tj_lacz_red_d1_d2_d_k_all grace_lacz_red_d1_d2_d_k_all]);
title([{'Red- D1 v D2'} {['Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}])
hold off

subplot(2,2,2)
j=cdfplot([tj_lacz_red_d1_d2_d_max_all grace_lacz_red_d1_d2_d_max_all]);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
hold on
h=cdfplot(arc_prom_red_d1_d2_d_max_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_d2_d_max_all);
set(l, 'LineStyle', '-', 'Color', 'm', 'LineWidth',2);
xlim([0 1])
xlabel(['Change in max Value'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_d2_d_max_all,[tj_lacz_red_d1_d2_d_max_all grace_lacz_red_d1_d2_d_max_all]);
[h_enh p_enh] = ttest2(arc_enh_red_d1_d2_d_max_all,[tj_lacz_red_d1_d2_d_max_all grace_lacz_red_d1_d2_d_max_all]);
title([{'Red- D1 v D2'} {['Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}])
hold off

subplot(2,2,3)
j=cdfplot([tj_lacz_red_d1_d3_d_k_all grace_lacz_red_d1_d3_d_k_all]);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
hold on
h=cdfplot(arc_prom_red_d1_d3_d_k_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_d3_d_k_all);
set(l, 'LineStyle', '-', 'Color', 'm', 'LineWidth',2);
legend(['LacZ- n = ' num2str(length([tj_lacz_red_d1_d3_d_k_all grace_lacz_red_d1_d3_d_k_all]))],...
    ['Prom- n = ' num2str(length(arc_prom_red_d1_d3_d_k_all))],...
    ['Enh- n = ' num2str(length(arc_enh_red_d1_d3_d_k_all))], 'Location', 'southeast')
xlim([0 30])
xlabel(['Change in k Value'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_d3_d_k_all,[tj_lacz_red_d1_d3_d_k_all grace_lacz_red_d1_d3_d_k_all]);
[h_enh p_enh] = ttest2(arc_enh_red_d1_d3_d_k_all,[tj_lacz_red_d1_d3_d_k_all grace_lacz_red_d1_d3_d_k_all]);
title([{'Red- D1 v D3'} {['Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}])
hold off

subplot(2,2,4)
j=cdfplot([tj_lacz_red_d1_d3_d_max_all grace_lacz_red_d1_d3_d_max_all]);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
hold on
h=cdfplot(arc_prom_red_d1_d3_d_max_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_d3_d_max_all);
set(l, 'LineStyle', '-', 'Color', 'm', 'LineWidth',2);
xlim([0 1])
xlabel(['Change in max Value'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_d3_d_max_all,[tj_lacz_red_d1_d3_d_max_all grace_lacz_red_d1_d3_d_max_all]);
[h_enh p_enh] = ttest2(arc_enh_red_d1_d3_d_max_all,[tj_lacz_red_d1_d3_d_max_all grace_lacz_red_d1_d3_d_max_all]);
title([{'Red- D1 v D3'} {['Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}])
hold off
sgtitle('D1 vs DN Change')
exportgraphics(gcf,fullfile(compfnout, ['allmice_change_kandmax_combLacZ.pdf']))

%% D1 comps
fig = figure;
subplot(2,2,1)
j=cdfplot([tj_lacz_red_d1_k_all grace_lacz_red_d1_k_all]);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
hold on
h=cdfplot(arc_prom_red_d1_k_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_k_all);
set(l, 'LineStyle', '-', 'Color', 'm', 'LineWidth',2);
legend(['LacZ- n = ' num2str(length([grace_lacz_red_d1_k_all tj_lacz_red_d1_k_all]))],...
    ['Prom- n = ' num2str(length(arc_prom_red_d1_k_all))],...
    ['Enh- n = ' num2str(length(arc_enh_red_d1_k_all))], 'Location', 'southeast')
xlim([0 30])
xlabel(['Day 1 k Value'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_k_all,[grace_lacz_red_d1_k_all tj_lacz_red_d1_k_all]);
[h_enh p_enh] = ttest2(arc_enh_red_d1_k_all,[tj_lacz_red_d1_k_all grace_lacz_red_d1_k_all]);
title([{'Red Responsive'},{['Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}]);
hold off

subplot(2,2,2)
j=cdfplot([tj_lacz_red_d1_max_all grace_lacz_red_d1_max_all]);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
hold on
h=cdfplot(arc_prom_red_d1_max_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_max_all);
set(l, 'LineStyle', '-', 'Color', 'm', 'LineWidth',2);
%legend(['LacZ'],['Prom'],['Enh'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Day 1 max Value'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_max_all,[grace_lacz_red_d1_max_all tj_lacz_red_d1_max_all]);
[h_enh p_enh] = ttest2(arc_enh_red_d1_max_all,[tj_lacz_red_d1_max_all grace_lacz_red_d1_max_all]);
title([{'Red Responsive'},{['Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}]);
hold off

subplot(2,2,3)
j=cdfplot([tj_lacz_red_d1_reliability_all grace_lacz_red_d1_reliability_all]);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
hold on
h=cdfplot(arc_prom_red_d1_reliability_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_reliability_all);
set(l, 'LineStyle', '-', 'Color', 'm', 'LineWidth',2);
%legend(['LacZ'],['Prom'],['Enh'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Day 1 reliability Value'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_reliability_all,[grace_lacz_red_d1_reliability_all tj_lacz_red_d1_reliability_all]);
[h_enh p_enh] = ttest2(arc_enh_red_d1_reliability_all,[tj_lacz_red_d1_reliability_all grace_lacz_red_d1_reliability_all]);
title([{'Red Responsive'},{['Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}]);
hold off

exportgraphics(gcf,fullfile(compfnout, ['allmice_day1_combLacZ.pdf']))

fig = figure;
subplot(2,2,1)
j=cdfplot([tj_lacz_red_d1_k_tuned_all grace_lacz_red_d1_k_tuned_all]);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
hold on
h=cdfplot(arc_prom_red_d1_k_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_k_tuned_all);
set(l, 'LineStyle', '-', 'Color', 'm', 'LineWidth',2);
legend(['LacZ- n = ' num2str(length([grace_lacz_red_d1_k_tuned_all tj_lacz_red_d1_k_tuned_all]))],...
    ['Prom- n = ' num2str(length(arc_prom_red_d1_k_tuned_all))],...
    ['Enh- n = ' num2str(length(arc_enh_red_d1_k_tuned_all))], 'Location', 'southeast')
xlim([0 30])
xlabel(['Day 1 k Value'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_k_tuned_all,[grace_lacz_red_d1_k_tuned_all tj_lacz_red_d1_k_tuned_all]);
[h_enh p_enh] = ttest2(arc_enh_red_d1_k_tuned_all,[tj_lacz_red_d1_k_tuned_all grace_lacz_red_d1_k_tuned_all]);
title([{'Red Tuned'},{['Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}]);
hold off

subplot(2,2,2)
j=cdfplot([tj_lacz_red_d1_max_tuned_all grace_lacz_red_d1_max_tuned_all]);
set(j, 'LineStyle', '-', 'Color', 'k', 'LineWidth',2);
hold on
h=cdfplot(arc_prom_red_d1_max_tuned_all);
set(h, 'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
l=cdfplot(arc_enh_red_d1_max_tuned_all);
set(l, 'LineStyle', '-', 'Color', 'm', 'LineWidth',2);
%legend(['LacZ'],['Prom'],['Enh'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Day 1 max Value'])
ylabel(['% of cells'])
[h_prom p_prom] = ttest2(arc_prom_red_d1_max_tuned_all,[grace_lacz_red_d1_max_tuned_all tj_lacz_red_d1_max_tuned_all]);
[h_enh p_enh] = ttest2(arc_enh_red_d1_max_tuned_all,[tj_lacz_red_d1_max_tuned_all grace_lacz_red_d1_max_tuned_all]);
title([{'Red Tuned'},{['Prom- p = ' num2str(chop(p_prom,2)) '; Enh- p = ' num2str(chop(p_enh,2))]}]);
hold off

figure;
subplot(2,2,1)
errorbar(0:22.5:157.5, mean(arc_prom_red_d1_respEaOri_all_sort,1),std(arc_prom_red_d1_respEaOri_all_sort,[],1)./sqrt(size(arc_prom_red_d1_respEaOri_all_sort,1)),'-or')
hold on
errorbar(0:22.5:157.5, mean(arc_enh_red_d1_respEaOri_all_sort,1),std(arc_enh_red_d1_respEaOri_all_sort,[],1)./sqrt(size(arc_enh_red_d1_respEaOri_all_sort,1)),'-om')
errorbar(0:22.5:157.5, mean([tj_lacz_red_d1_respEaOri_all_sort; grace_lacz_red_d1_respEaOri_all_sort] ,1),std([tj_lacz_red_d1_respEaOri_all_sort; grace_lacz_red_d1_respEaOri_all_sort],[],1)./sqrt(size([tj_lacz_red_d1_respEaOri_all_sort; grace_lacz_red_d1_respEaOri_all_sort],1)),'-ok')
ylim([-0.1 1.25])
ylabel('Normalized dF/F')
xlabel('Orientation (deg)')
title('Responsive')
subplot(2,2,2)
errorbar(0:22.5:157.5, mean(arc_prom_red_d1_respEaOri_tuned_all_sort,1),std(arc_prom_red_d1_respEaOri_tuned_all_sort,[],1)./sqrt(size(arc_prom_red_d1_respEaOri_tuned_all_sort,1)),'-or')
hold on
errorbar(0:22.5:157.5, mean(arc_enh_red_d1_respEaOri_tuned_all_sort,1),std(arc_enh_red_d1_respEaOri_tuned_all_sort,[],1)./sqrt(size(arc_enh_red_d1_respEaOri_tuned_all_sort,1)),'-om')
errorbar(0:22.5:157.5, mean([tj_lacz_red_d1_respEaOri_tuned_all_sort; grace_lacz_red_d1_respEaOri_tuned_all_sort] ,1),std([tj_lacz_red_d1_respEaOri_tuned_all_sort; grace_lacz_red_d1_respEaOri_tuned_all_sort],[],1)./sqrt(size([tj_lacz_red_d1_respEaOri_tuned_all_sort; grace_lacz_red_d1_respEaOri_tuned_all_sort],1)),'-ok')
ylim([-0.1 1.25])
ylabel('Normalized dF/F')
xlabel('Orientation (deg)')
title('Tuned')
subplot(2,2,3)
shadedErrorBar(0:180, mean(arc_prom_red_d1_vonMises_all_sort,1),std(arc_prom_red_d1_vonMises_all_sort,[],1)./sqrt(size(arc_prom_red_d1_vonMises_all_sort,1)),'-r');
hold on
shadedErrorBar(0:180, mean(arc_enh_red_d1_vonMises_all_sort,1),std(arc_enh_red_d1_vonMises_all_sort,[],1)./sqrt(size(arc_enh_red_d1_vonMises_all_sort,1)),'-m');
shadedErrorBar(0:180, mean([tj_lacz_red_d1_vonMises_all_sort; grace_lacz_red_d1_vonMises_all_sort] ,1),std([tj_lacz_red_d1_vonMises_all_sort; grace_lacz_red_d1_vonMises_all_sort],[],1)./sqrt(size([tj_lacz_red_d1_vonMises_all_sort; grace_lacz_red_d1_vonMises_all_sort],1)),'-k');
ylim([-0.1 1.25])
ylabel('Normalized dF/F')
xlabel('Orientation (deg)')
title('Responsive')
subplot(2,2,4)
shadedErrorBar(0:180, mean(arc_prom_red_d1_vonMises_tuned_all_sort,1),std(arc_prom_red_d1_vonMises_tuned_all_sort,[],1)./sqrt(size(arc_prom_red_d1_vonMises_tuned_all_sort,1)),'-r');
hold on
shadedErrorBar(0:180, mean(arc_enh_red_d1_vonMises_tuned_all_sort,1),std(arc_enh_red_d1_vonMises_tuned_all_sort,[],1)./sqrt(size(arc_enh_red_d1_vonMises_tuned_all_sort,1)),'-m');
shadedErrorBar(0:180, mean([tj_lacz_red_d1_vonMises_tuned_all_sort; grace_lacz_red_d1_vonMises_tuned_all_sort] ,1),std([tj_lacz_red_d1_vonMises_tuned_all_sort; grace_lacz_red_d1_vonMises_tuned_all_sort],[],1)./sqrt(size([tj_lacz_red_d1_vonMises_tuned_all_sort; grace_lacz_red_d1_vonMises_tuned_all_sort],1)),'-k');
ylim([-0.1 1.25])
ylabel('Normalized dF/F')
xlabel('Orientation (deg)')
title('Tuned')

exportgraphics(gcf,fullfile(compfnout, ['allmice_day1_tuningcurve_combLacZ.pdf']))