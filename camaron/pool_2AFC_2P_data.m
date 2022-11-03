clearvars
close all
clc

%% Load expt list
dataset = 'oriAdapt_V1';
eval(dataset);

mouse = 'i475'; % indicate mouse data to pool

filename = [mouse '_image_matching.mat'];
load(filename, "expt_list_final")

for i = 1:length(expt_list_final)
    %% Get Dataset info
    
    iexp = expt_list_final(i)
    date = expt(iexp).date;
    
    user = expt(iexp).folder;
    if user == "lindsey"
        data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey';
    elseif user == "camaron"
        data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';
    end
    
    %% For each dataset load key variables
    
    %Behavior
    run_str = ['runs-' expt(iexp).runs(1,:)];
    b_input = load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']));
    b_adapt_resp = load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_adaptResp.mat']));
    b_stim_resp = load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimResp.mat']));
    b_stim_data = load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']));
    
    %Passive
    pass_str = ['runs-' expt(iexp).pass_run];
    p_input = load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_input.mat']));
    p_adapt_resp = load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_adaptResp.mat']));
    p_stim_resp = load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_stimResp.mat']));
    p_stim_data = load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_stimData.mat']));
    
    %dir_tuning
    dir_tuning_str = ['runs-' expt(iexp).dirtuning];
    ori_tuning_and_fits = load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_tuning_str], [date '_' mouse '_' dir_tuning_str '_oriTuningAndFits.mat']));
    ori_tuning_info = load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_tuning_str], [date '_' mouse '_' dir_tuning_str '_oriTuningInfo.mat']));
    
    %% Pooling variables 
    
    %Grab stim timings across behaving and passive to check for alignment
    adaptPeriodMs(i,:) = [b_input.adapt_input.adaptPeriodMs, p_input.pass_input.adaptPeriodMs];
    dynAdaptPeriodMs(i,:) = [b_input.adapt_input.dynAdaptPeriodMs, p_input.pass_input.dynAdaptPeriodMs];
    dynAdaptFlashOffMs(i,:) = [b_input.adapt_input.dynAdaptFlashOffMs, p_input.pass_input.dynAdaptFlashOffMs];
    dates{i} = date;
    
    % Pool indices of tuned Cells
    tuned_cells(i,:) = ori_tuning_info.tunedCells;  
    ori_bins = [0 45 90 135]; % for tuned cells 
    
    
    %% Pool TCs and trial indices for conditions
    
    %Behavior 
    % TCs
    b_data_adapt_dfof_all{i} = b_adapt_resp.data_adapt_dfof; % Expt(i) = [trial_time X nCells X Trials]
    b_data_stim_dfof_all{i} = b_stim_resp.data_stim_dfof; % Expt(i) = [trial_time X nCells X Trials]
    
    % find index of trial conditions
    aGratingOri = celleqel2mat_padded(b_input.adapt_input.aGratingDirectionDeg);
    aGratingContrast = celleqel2mat_padded(b_input.adapt_input.aGratingContrast);
    aOris = unique(aGratingOri);
    naOri = length(aOris);
    ind_con = find(aGratingContrast == 1);
    ind_ori = find(aGratingOri==aOris);
    adapt_trial_ind = intersect(ind_ori,ind_con); 
    control_trial_ind = setdiff(ind_ori,ind_con); 
    
    b_adapt_trial_ind_all{i} = adapt_trial_ind; % trials when adapter is on
    b_control_trial_ind_all{i} = control_trial_ind; % trials when adapter is off
    
    % find index of responsive cells
    b_adapt_resp_ind_all(i) = b_adapt_resp.adapt_resp_ind;
    b_stim_resp_ind_all{i} = b_stim_resp.stim_resp_ind;
    
    % anyalysis window used to determine responsive cells (save these)
    base_win = b_stim_resp.base_win;
    resp_win = b_stim_resp.resp_win;
    base_win_all = b_adapt_resp.base_win_all;
    resp_win_all = b_adapt_resp.resp_win_all;
    
    % -----
    
    %Passive 
    % TCs
    p_data_adapt_dfof_all{i} = p_adapt_resp.data_adapt_dfof; % Expt(i) = [trial_time X nCells X Trials]
    p_data_stim_dfof_all{i} = p_stim_resp.data_stim_dfof; % Expt(i) = [trial_time X nCells X Trials]
    
    % find trial conditions
    aGratingOri = celleqel2mat_padded(p_input.pass_input.aGratingDirectionDeg);
    aGratingContrast = celleqel2mat_padded(p_input.pass_input.aGratingContrast);
    aOris = unique(aGratingOri);
    naOri = length(aOris);
    ind_con = find(aGratingContrast == 1);
    ind_ori = find(aGratingOri==aOris);
    adapt_trial_ind = intersect(ind_ori,ind_con); 
    control_trial_ind = setdiff(ind_ori,ind_con); 
    
    p_adapt_trial_ind_all{i} = adapt_trial_ind; % trials when adapter is on
    p_control_trial_ind_all{i} = control_trial_ind; % trials when adapter is off
    
    % find index of responsive cells
    p_adapt_resp_ind_all(i) = p_adapt_resp.adapt_resp_ind;
    p_stim_resp_ind_all{i} = p_stim_resp.stim_resp_ind;
    
    % -----
    
    % Masks (Interneurons)
    b_mask_label_all{i} =  b_adapt_resp.mask_label;
    p_mask_label_all{i} =  p_adapt_resp.mask_label;
    
    % -----
    
    adaptor_vline = [20 31 42 53];
    target_vline = 64;

end

%% 
clearvars -except adaptor_vline adaptPeriodMsb_adapt_resp_ind_all ...
    b_adapt_trial_ind_all b_control_trial_ind_all b_data_adapt_dfof_all ...
    b_data_stim_dfof_all b_mask_label_all b_stim_resp_ind_all base_win ...
    base_win_all dates dynAdaptFlashOff dynAdaptPeriodMs expt_list_final ...
    mouse ori_bins p_adapt_resp_ind_all p_adapt_trial_ind_all ...
    p_control_trial_ind_all p_data_adapt_dfof_all p_data_stim_dfof_all ...
    p_mask_label_all p_stim_resp_ind_all resp_win resp_win_all ...
    target_vline tuned_cells

%%
cd(['Z:\All_Staff\home\camaron\Analysis\2P\' mouse])
filename = [mouse '_pooled_2AFC_2P_data.mat'];
save(filename)

