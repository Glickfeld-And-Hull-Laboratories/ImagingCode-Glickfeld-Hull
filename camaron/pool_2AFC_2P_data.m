clearvars
close all
clc

%% Load expt list
dataset = 'oriAdapt_V1';
eval(dataset);

mouse = 'i475'; % indicate mouse data to pool
cd(['Z:\All_Staff\home\camaron\Analysis\2P\' mouse])
filename = [mouse '_image_matching.mat'];
load(filename, "expt_list_final", "z_pos_good_list")

% find zpos for this final list, and save for later
[C, ia, ib] = intersect(z_pos_good_list(1,:), expt_list_final); 
z_pos_expt_list_final = z_pos_good_list(:,ia);



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
    b_adapt_win_dfof_all{i} = b_adapt_resp.data_adapt_dfof; % Expt(i) = [trial_time X nCells X Trials]
    b_stim_win_dfof_all{i} = b_stim_resp.data_stim_dfof; % Expt(i) = [trial_time X nCells X Trials]
    
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

    % Break data into trial conditions (adapt when adapt on, target for both conditions)
    b_adapt_win_dfof_all_adapt{i} = b_adapt_resp.data_adapt_dfof(:, :, adapt_trial_ind); % Expt(i) = [trial_time X nCells X Adapt_Trials]
    b_stim_win_dfof_all_control{i} = b_stim_resp.data_stim_dfof(:, :, control_trial_ind); % Expt(i) = [trial_time X nCells X Control_Trials]
    b_stim_win_dfof_all_adapt{i} = b_stim_resp.data_stim_dfof(:, :, adapt_trial_ind); % Expt(i) = [trial_time X nCells X Adapt_Trials]

    % anyalysis window used to determine responsive cells (save these)
    base_win = b_stim_resp.base_win;
    resp_win = b_stim_resp.resp_win;
    base_win_all = b_adapt_resp.base_win_all;
    resp_win_all = b_adapt_resp.resp_win_all;

    % find index of responsive cells
    b_ind_sig_resp_adapt_fromTC(i) = b_adapt_resp.adapt_resp_ind; % Uses data from TC_extraction, compare with new ttest
    b_ind_sig_resp_stim_fromTC{i} = b_stim_resp.stim_resp_ind;

    b_ind_sig_resp_adapt_all{i} = find_respCells(b_adapt_win_dfof_all_adapt{i}, base_win, resp_win);
    b_ind_sig_resp_stim_control_all{i} = find_respCells(b_stim_win_dfof_all_control{i}, base_win, resp_win);
    b_ind_sig_resp_stim_adapt_all{i} = find_respCells(b_stim_win_dfof_all_adapt{i}, base_win, resp_win); % 
    % Index of sig responsive cells when adaptor is on, for use when
    % looking at cell percentages.

    
    %%
    
    %Passive 
    % TCs
    p_adapt_win_dfof_all{i} = p_adapt_resp.data_adapt_dfof; % Expt(i) = [trial_time X nCells X Trials]
    p_stim_win_dfof_all{i} = p_stim_resp.data_stim_dfof; % Expt(i) = [trial_time X nCells X Trials]
    
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

    % Break data into trial conditions (adapt when adapt on, target for both conditions)
    p_adapt_win_dfof_all_adapt{i} = p_adapt_resp.data_adapt_dfof(:, :, adapt_trial_ind); % Expt(i) = [trial_time X nCells X Adapt_Trials]
    p_stim_win_dfof_all_control{i} = p_stim_resp.data_stim_dfof(:, :, control_trial_ind); % Expt(i) = [trial_time X nCells X Control_Trials]
    p_stim_win_dfof_all_adapt{i} = p_stim_resp.data_stim_dfof(:, :, adapt_trial_ind); % Expt(i) = [trial_time X nCells X Adapt_Trials]
    
    % find index of responsive cells
    p_ind_sig_resp_adapt_fromTC(i) = p_adapt_resp.adapt_resp_ind; % 
    p_ind_sig_resp_stim_fromTC{i} = p_stim_resp.stim_resp_ind;

    p_ind_sig_resp_adapt_all{i} = find_respCells(p_adapt_win_dfof_all_adapt{i}, base_win, resp_win);
    p_ind_sig_resp_stim_control_all{i} = find_respCells(p_stim_win_dfof_all_control{i}, base_win, resp_win);
    p_ind_sig_resp_stim_adapt_all{i} = find_respCells(p_stim_win_dfof_all_adapt{i}, base_win, resp_win); % 

    
    % -----
    
    % Masks (Interneurons)  - Dont forget to use these!
    b_mask_label_all{i} =  b_adapt_resp.mask_label;
    p_mask_label_all{i} =  p_adapt_resp.mask_label;
    
    % -----
    
    adaptor_vline = [20 31 42 53];
    target_vline = 64;
   
    %% Equality test - for significance

    % 
    % The following equality tests describes if determination of
    % significantly responsive cells is equivalent across the TC extraction
    % script and find_respCells; It is for adaptor responses, it is NOT for
    % stim responses; stim responses in TC analysis include control and 
    % adapt trials, here we only use control trials.
    
    isequal(b_ind_sig_resp_adapt_fromTC{i}, b_ind_sig_resp_adapt_all{i});
    isequal(b_ind_sig_resp_stim_fromTC{i}, b_ind_sig_resp_stim_control_all{i});
    
    isequal(p_ind_sig_resp_adapt_fromTC{i}, p_ind_sig_resp_adapt_all{i});
    isequal(p_ind_sig_resp_stim_fromTC{i}, p_ind_sig_resp_stim_control_all{i});

    %% Intersection of responsive and tuned cells

    for bin = 1:length(tuned_cells(i,:))
        b_ind_tuned_sig_resp_adapt{bin} = intersect(b_ind_sig_resp_adapt_all{i}, tuned_cells{i, bin});
        b_ind_tuned_sig_resp_stim_control{bin} = intersect(b_ind_sig_resp_stim_control_all{i}, tuned_cells{i, bin});
        b_ind_tuned_sig_resp_stim_adapt{bin} = intersect(b_ind_sig_resp_stim_adapt_all{i}, tuned_cells{i, bin});


        p_ind_tuned_sig_resp_adapt{bin} = intersect(p_ind_sig_resp_adapt_all{i}, tuned_cells{i, bin});
        p_ind_tuned_sig_resp_stim_control{bin} = intersect(p_ind_sig_resp_stim_control_all{i}, tuned_cells{i, bin});
        p_ind_tuned_sig_resp_stim_adapt{bin} = intersect(p_ind_sig_resp_stim_adapt_all{i}, tuned_cells{i, bin});

    end

    b_ind_tuned_sig_resp_adapt_all{i} = b_ind_tuned_sig_resp_adapt;
    b_ind_tuned_sig_resp_stim_control_all{i} = b_ind_tuned_sig_resp_stim_control;
    b_ind_tuned_sig_resp_stim_adapt_all{i} = b_ind_tuned_sig_resp_stim_adapt;


    p_ind_tuned_sig_resp_adapt_all{i} = p_ind_tuned_sig_resp_adapt;
    p_ind_tuned_sig_resp_stim_control_all{i} = p_ind_tuned_sig_resp_stim_control;
    p_ind_tuned_sig_resp_stim_adapt_all{i} = p_ind_tuned_sig_resp_stim_adapt;


    %% Segment cells from conditional dF/F by responsivitiy and tuning

    %Responsivity 
    b_adapt_sigDuringAdapt_adapt{i} = b_adapt_win_dfof_all_adapt{i}(:,b_ind_sig_resp_adapt_all{i},:); % adapt response, on adapt trials, with cells sig resposive during adaptorON condition
    b_stim_sigDuringControl_control{i} = b_stim_win_dfof_all_control{i}(:,b_ind_sig_resp_stim_control_all{i},:); % stim response, on control trials, with cells sig during control condition
    b_stim_sigDuringControl_adapt{i} = b_stim_win_dfof_all_adapt{i}(:,b_ind_sig_resp_stim_control_all{i},:); % stim response, on adapt trials, with cells sig resposive during control condition
    b_stim_sigDuringAdapt_adapt{i} = b_stim_win_dfof_all_adapt{i}(:,b_ind_sig_resp_stim_adapt_all{i},:); % stim response, on adapt trials, with cells sig resposive during adaptorON condition


    p_adapt_sigDuringAdapt_adapt{i} = p_adapt_win_dfof_all_adapt{i}(:,p_ind_sig_resp_adapt_all{i},:);
    p_stim_sigDuringControl_control{i} = p_stim_win_dfof_all_control{i}(:,p_ind_sig_resp_stim_control_all{i},:);
    p_stim_sigDuringControl_adapt{i} = p_stim_win_dfof_all_adapt{i}(:,p_ind_sig_resp_stim_control_all{i},:);
    p_stim_sigDuringAdapt_adapt{i} = p_stim_win_dfof_all_adapt{i}(:,p_ind_sig_resp_stim_adapt_all{i},:); % stim response, on adapt trials, with cells sig resposive during adaptorON condition


    %Responsivity + tuning 

    for bin = 1:length(tuned_cells(1,:))
        b_tuned_adapt_sigDuringAdapt_adapt{i, bin} = b_adapt_win_dfof_all_adapt{i}(:,b_ind_tuned_sig_resp_adapt{bin},:); % adapt response, on adapt trials, with TUNED cells sig responsive during adaptorON condition
        b_tuned_stim_sigDuringControl_control{i, bin} = b_stim_win_dfof_all_control{i}(:,b_ind_tuned_sig_resp_stim_control{bin},:); % stim response, on control trials, with TUNED cells sig responsive during control condition
        b_tuned_stim_sigDuringControl_adapt{i, bin} = b_stim_win_dfof_all_adapt{i}(:,b_ind_tuned_sig_resp_stim_control{bin},:); % stim response, on adapt trials, with TUNED cells sig responsive during control condition
        b_tuned_stim_sigDuringAdapt_adapt{i, bin} = b_stim_win_dfof_all_adapt{i}(:,b_ind_tuned_sig_resp_stim_adapt{bin},:); % stim response, on adapt trials, with TUNED cells sig responsive during adaptorON condition


        p_tuned_adapt_sigDuringAdapt_adapt{i, bin} = p_adapt_win_dfof_all_adapt{i}(:,p_ind_tuned_sig_resp_adapt{bin},:);
        p_tuned_stim_sigDuringControl_control{i, bin} = p_stim_win_dfof_all_control{i}(:,p_ind_tuned_sig_resp_stim_control{bin},:);
        p_tuned_stim_sigDuringControl_adapt{i, bin} = p_stim_win_dfof_all_adapt{i}(:,p_ind_tuned_sig_resp_stim_control{bin},:);
        p_tuned_stim_sigDuringAdapt_adapt{i, bin} = p_stim_win_dfof_all_adapt{i}(:,p_ind_tuned_sig_resp_stim_adapt{bin},:); % stim response, on adapt trials, with TUNED cells sig responsive during adaptorON condition

    end

   
end

%% Pool mean cell activity across trials!
% Uses local helper function; 'Run Section' to use.

% All by condition

b_data_adapt_trial_mean = mean_cell_resp_by_trial(b_adapt_win_dfof_all_adapt);
b_data_stim_control_trial_mean = mean_cell_resp_by_trial(b_stim_win_dfof_all_control);
b_data_stim_adapt_trial_mean = mean_cell_resp_by_trial(b_stim_win_dfof_all_adapt);

p_data_adapt_trial_mean = mean_cell_resp_by_trial(p_adapt_win_dfof_all_adapt);
p_data_stim_control_trial_mean = mean_cell_resp_by_trial(p_stim_win_dfof_all_control);
p_data_stim_adapt_trial_mean = mean_cell_resp_by_trial(p_stim_win_dfof_all_adapt);

% Responsive by condition

b_adapt_sigDuringAdapt_adapt_trial_mean = mean_cell_resp_by_trial(b_adapt_sigDuringAdapt_adapt);
b_stim_sigDuringControl_control_trial_mean = mean_cell_resp_by_trial(b_stim_sigDuringControl_control);
b_stim_sigDuringControl_adapt_trial_mean = mean_cell_resp_by_trial(b_stim_sigDuringControl_adapt);
b_stim_sigDuringAdapt_adapt_trial_mean = mean_cell_resp_by_trial(b_stim_sigDuringAdapt_adapt);


p_adapt_sigDuringAdapt_adapt_trial_mean = mean_cell_resp_by_trial(p_adapt_sigDuringAdapt_adapt);
p_stim_sigDuringControl_control_trial_mean = mean_cell_resp_by_trial(p_stim_sigDuringControl_control);
p_stim_sigDuringControl_adapt_trial_mean = mean_cell_resp_by_trial(p_stim_sigDuringControl_adapt);
p_stim_sigDuringAdapt_adapt_trial_mean = mean_cell_resp_by_trial(p_stim_sigDuringAdapt_adapt);

% Responsive + tuned by condition 

b_tuned_adapt_sigDuringAdapt_adapt_trial_mean = mean_tuned_cell_resp_by_trial(b_tuned_adapt_sigDuringAdapt_adapt);
b_tuned_stim_sigDuringControl_control_trial_mean = mean_tuned_cell_resp_by_trial(b_tuned_stim_sigDuringControl_control);
b_tuned_stim_sigDuringControl_adapt_trial_mean = mean_tuned_cell_resp_by_trial(b_tuned_stim_sigDuringControl_adapt);
b_tuned_stim_sigDuringAdapt_adapt_trial_mean = mean_tuned_cell_resp_by_trial(b_tuned_stim_sigDuringAdapt_adapt);


p_tuned_adapt_sigDuringAdapt_adapt_trial_mean = mean_tuned_cell_resp_by_trial(p_tuned_adapt_sigDuringAdapt_adapt);
p_tuned_stim_sigDuringControl_control_trial_mean = mean_tuned_cell_resp_by_trial(p_tuned_stim_sigDuringControl_control);
p_tuned_stim_sigDuringControl_adapt_trial_mean = mean_tuned_cell_resp_by_trial(p_tuned_stim_sigDuringControl_adapt);
p_tuned_stim_sigDuringAdapt_adapt_trial_mean = mean_tuned_cell_resp_by_trial(p_tuned_stim_sigDuringAdapt_adapt);

% Reminder that experiments (cells) are at different depths...



%% Get pooled cell counts for each condition before averaging across cell (uneccessary with TC_stats function)

nAllCells = size(b_data_adapt_trial_mean,2);

b_adapt_sigDuringAdapt_nCells = size(b_adapt_sigDuringAdapt_adapt_trial_mean,2);
b_stim_sigDuringControl_nCells = size(b_stim_sigDuringControl_control_trial_mean,2);
b_stim_sigDuringAdapt_nCells = size(b_stim_sigDuringAdapt_adapt_trial_mean,2);


p_adapt_sigDuringAdapt_nCells = size(p_adapt_sigDuringAdapt_adapt_trial_mean,2);
p_stim_sigDuringControl_nCells = size(p_stim_sigDuringControl_control_trial_mean,2);
p_stim_sigDuringAdapt_nCells = size(p_stim_sigDuringAdapt_adapt_trial_mean,2);


for i = 1:length(ori_bins)
    b_tuned_adapt_sigDuringAdapt_nCells(i) = size(b_tuned_adapt_sigDuringAdapt_adapt_trial_mean{i}, 2);
    b_tuned_stim_sigDuringControl_nCells(i) = size(b_tuned_stim_sigDuringControl_control_trial_mean{i}, 2);
    b_tuned_stim_sigDuringAdapt_nCells(i) = size(b_tuned_stim_sigDuringAdapt_adapt_trial_mean{i}, 2);


    p_tuned_adapt_sigDuringAdapt_nCells(i) = size(p_tuned_adapt_sigDuringAdapt_adapt_trial_mean{i}, 2);
    p_tuned_stim_sigDuringControl_nCells(i) = size(p_tuned_stim_sigDuringControl_control_trial_mean{i}, 2);
    p_tuned_stim_sigDuringAdapt_nCells(i) = size(p_tuned_stim_sigDuringAdapt_adapt_trial_mean{i}, 2);

end


%% Cell average TC and SEM for each condition

% All by condition
b_data_adapt_stats = TC_stats(b_data_adapt_trial_mean, resp_win);
b_data_stim_control_stats  = TC_stats(b_data_stim_control_trial_mean, resp_win);
b_data_stim_adapt_stats  = TC_stats(b_data_stim_adapt_trial_mean, resp_win);

p_data_adapt_stats  = TC_stats(p_data_adapt_trial_mean, resp_win);
p_data_stim_control_stats  = TC_stats(p_data_stim_control_trial_mean, resp_win);
p_data_stim_adapt_stats  = TC_stats(p_data_stim_adapt_trial_mean, resp_win);

% Responsive by condition

b_adapt_sigDuringAdapt_adapt_stats  = TC_stats(b_adapt_sigDuringAdapt_adapt_trial_mean, resp_win);
b_stim_sigDuringControl_control_stats  = TC_stats(b_stim_sigDuringControl_control_trial_mean, resp_win);
b_stim_sigDuringControl_adapt_stats  = TC_stats(b_stim_sigDuringControl_adapt_trial_mean, resp_win);
b_stim_sigDuringAdapt_adapt_stats  = TC_stats(b_stim_sigDuringAdapt_adapt_trial_mean, resp_win);
 
p_adapt_sigDuringAdapt_adapt_stats  = TC_stats(p_adapt_sigDuringAdapt_adapt_trial_mean, resp_win);
p_stim_sigDuringControl_control_stats  = TC_stats(p_stim_sigDuringControl_control_trial_mean, resp_win);
p_stim_sigDuringControl_adapt_stats  = TC_stats(p_stim_sigDuringControl_adapt_trial_mean, resp_win);
p_stim_sigDuringAdapt_adapt_stats  = TC_stats(p_stim_sigDuringAdapt_adapt_trial_mean, resp_win);

% Responsive + tuned by condition 
for i = 1:4

    b_tuned_adapt_sigDuringAdapt_adapt_stats(i)  = TC_stats(b_tuned_adapt_sigDuringAdapt_adapt_trial_mean{i}, resp_win);
    b_tuned_stim_sigDuringControl_control_stats(i)  = TC_stats(b_tuned_stim_sigDuringControl_control_trial_mean{i}, resp_win);
    b_tuned_stim_sigDuringControl_adapt_stats(i)  = TC_stats(b_tuned_stim_sigDuringControl_adapt_trial_mean{i}, resp_win);
    b_tuned_stim_sigDuringAdapt_adapt_stats(i)  = TC_stats(b_tuned_stim_sigDuringAdapt_adapt_trial_mean{i}, resp_win);
    
    p_tuned_adapt_sigDuringAdapt_adapt_stats(i)  = TC_stats(p_tuned_adapt_sigDuringAdapt_adapt_trial_mean{i}, resp_win);
    p_tuned_stim_sigDuringControl_control_stats(i)  = TC_stats(p_tuned_stim_sigDuringControl_control_trial_mean{i}, resp_win);
    p_tuned_stim_sigDuringControl_adapt_stats(i)  = TC_stats(p_tuned_stim_sigDuringControl_adapt_trial_mean{i}, resp_win);
    p_tuned_stim_sigDuringAdapt_adapt_stats(i)  = TC_stats(p_tuned_stim_sigDuringAdapt_adapt_trial_mean{i}, resp_win);


end

%% Cell average dF/F (during response window) and SEM for each condition (uneccessary with TC_stats function)
%
% mean(b_data_adapt_stats.mean_TC(resp_win))
% OR?
% 
% win_Resp_by_cell = b_data_adapt_trial_mean(resp_win, :) % Responses in window by cell
% win_Resp_by_cell_mean = mean(win_Resp_by_cell, 1) % average response in window by cell
% avg_win_resp = mean(win_Resp_by_cell_mean) % grand average 

%% Cell adapt index for each condition





%% Most efficient way to clear uneeded variables and save the workspace?
% 
% clearvars -except adaptor_vline adaptPeriodMsb_ind_sig_resp_adapt_all ...
%     b_adapt_trial_ind_all b_control_trial_ind_all b_adapt_win_dfof_all ...
%     b_stim_win_dfof_all b_mask_label_all b_ind_sig_resp_stim_control_all base_win ...
%     base_win_all dates dynAdaptFlashOff dynAdaptPeriodMs expt_list_final ...
%     mouse ori_bins p_ind_sig_resp_adapt_all p_adapt_trial_ind_all ...
%     p_control_trial_ind_all p_adapt_win_dfof_all p_stim_win_dfof_all ...
%     p_mask_label_all p_ind_sig_resp_stim_control_all resp_win resp_win_all ...
%     target_vline tuned_cells b_adapt_win_dfof_all_adapt ...
%     b_stim_win_dfof_all_control b_stim_win_dfof_all_adapt ...
%     p_adapt_win_dfof_all_adapt p_stim_win_dfof_all_control ...
%     p_stim_win_dfof_all_adapt z_pos_expt_list_final



%%
cd(['Z:\All_Staff\home\camaron\Analysis\2P\' mouse])
filename = [mouse '_pooled_2AFC_2P_data.mat'];
save(filename)


%% Local functions

function pooled_data = mean_cell_resp_by_trial(dFoF)
% Finds mean cell activity across trials for all experiments
    pooled_data = [];
    for i = 1:length(dFoF)
        data = dFoF{i};
        mean_data = mean(data, 3, 'omitnan');
        pooled_data = cat(2, pooled_data, mean_data);
    end
end

function pooled_tuned_data = mean_tuned_cell_resp_by_trial(tuned_data)
    pooled_tuned_data = {};
    for bin = 1:4
        data_set = tuned_data(:,bin);
        pooled_tuned_data{bin} = mean_cell_resp_by_trial(data_set);
    end
end


function data_struct = TC_stats(trial_averaged_data, time_win)
    %TC mean, std, sem, cell count and more stats for a time X cell response matrix

    data_mean = mean(trial_averaged_data,2);
    data_std = std(trial_averaged_data, 0, 2);
    data_count = size(trial_averaged_data,2);
    data_sem = data_std / sqrt(data_count);

    data_struct.count_cells = data_count;
    data_struct.mean_TC = data_mean;
    data_struct.std_TC = data_std;
    data_struct.sem_TC = data_sem;

    % If response window is given return cell average DF/F and sem...

    if exist('time_win', 'var')

        win_df_by_cell = trial_averaged_data(time_win, :); % Responses in window by cell
        win_df_by_cell_mean = mean(win_df_by_cell, 1); % average response in window by cell (save this and below)
        win_df_mean = mean(win_df_by_cell_mean); % grand average 
        win_df_std = std(win_df_by_cell_mean);
        win_df_sem = win_df_std / sqrt(data_count);

        data_struct.time_win = time_win;
        data_struct.mean_df_by_cell = win_df_by_cell_mean;
        data_struct.mean_df = win_df_mean;
        data_struct.std_df = win_df_std;
        data_struct.sem_df = win_df_sem;        


    end
    
end

