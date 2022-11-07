clearvars
close all
clc

%% Load expt list
dataset = 'oriAdapt_V1';
eval(dataset);

mouse = 'i475'; % indicate mouse data to pool
cd(['Z:\All_Staff\home\camaron\Analysis\2P\' mouse])
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

    % Break data into trial conditions (adapt when adapt on, target for both conditions)
    b_data_adapt_dfof_all_adapt{i} = b_adapt_resp.data_adapt_dfof(:, :, adapt_trial_ind); % Expt(i) = [trial_time X nCells X Adapt_Trials]
    b_data_stim_dfof_all_control{i} = b_stim_resp.data_stim_dfof(:, :, control_trial_ind); % Expt(i) = [trial_time X nCells X Control_Trials]
    b_data_stim_dfof_all_adapt{i} = b_stim_resp.data_stim_dfof(:, :, adapt_trial_ind); % Expt(i) = [trial_time X nCells X Adapt_Trials]

    % anyalysis window used to determine responsive cells (save these)
    base_win = b_stim_resp.base_win;
    resp_win = b_stim_resp.resp_win;
    base_win_all = b_adapt_resp.base_win_all;
    resp_win_all = b_adapt_resp.resp_win_all;

    % find index of responsive cells
    b_adapt_resp_ind_all_TC(i) = b_adapt_resp.adapt_resp_ind; % Uses data from TC_extraction, compare with new ttest
    b_stim_resp_ind_all_TC{i} = b_stim_resp.stim_resp_ind;

    b_adapt_resp_ind_all{i} = find_respCells(b_data_adapt_dfof_all_adapt{i}, base_win, resp_win);
    b_stim_resp_ind_all{i} = find_respCells(b_data_stim_dfof_all_control{i}, base_win, resp_win);
    
    %%
    
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

    % Break data into trial conditions (adapt when adapt on, target for both conditions)
    p_data_adapt_dfof_all_adapt{i} = p_adapt_resp.data_adapt_dfof(:, :, adapt_trial_ind); % Expt(i) = [trial_time X nCells X Adapt_Trials]
    p_data_stim_dfof_all_control{i} = p_stim_resp.data_stim_dfof(:, :, control_trial_ind); % Expt(i) = [trial_time X nCells X Control_Trials]
    p_data_stim_dfof_all_adapt{i} = p_stim_resp.data_stim_dfof(:, :, adapt_trial_ind); % Expt(i) = [trial_time X nCells X Adapt_Trials]
    
    % find index of responsive cells
    p_adapt_resp_ind_all_TC(i) = p_adapt_resp.adapt_resp_ind;
    p_stim_resp_ind_all_TC{i} = p_stim_resp.stim_resp_ind;

    p_adapt_resp_ind_all{i} = find_respCells(p_data_adapt_dfof_all_adapt{i}, base_win, resp_win);
    p_stim_resp_ind_all{i} = find_respCells(p_data_stim_dfof_all_control{i}, base_win, resp_win);
    
    % -----
    
    % Masks (Interneurons)
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
    
    isequal(b_adapt_resp_ind_all_TC{i}, b_adapt_resp_ind_all{i})
    isequal(b_stim_resp_ind_all_TC{i}, b_stim_resp_ind_all{i})
    
    isequal(p_adapt_resp_ind_all_TC{i}, p_adapt_resp_ind_all{i})
    isequal(p_stim_resp_ind_all_TC{i}, p_stim_resp_ind_all{i})

    %% Intersection of responsive and tuned cells

    for bin = 1:length(tuned_cells(i,:))
        b_tuned_adapt_resp_ind{bin} = intersect(b_adapt_resp_ind_all{i}, tuned_cells{i, bin});
        b_tuned_stim_resp_ind{bin} = intersect(b_stim_resp_ind_all{i}, tuned_cells{i, bin});
        p_tuned_adapt_resp_ind{bin} = intersect(p_adapt_resp_ind_all{i}, tuned_cells{i, bin});
        p_tuned_stim_resp_ind{bin} = intersect(p_stim_resp_ind_all{i}, tuned_cells{i, bin});
    end

    b_tuned_adapt_resp_ind_all{i} = b_tuned_adapt_resp_ind;
    b_tuned_stim_resp_ind_all{i} = b_tuned_stim_resp_ind;
    p_tuned_adapt_resp_ind_all{i} = p_tuned_adapt_resp_ind;
    p_tuned_stim_resp_ind_all{i} = p_tuned_stim_resp_ind;

    %% Segment cells from conditional dF/F by responsivitiy and tuning

    %Responsivity 
    b_adapt_resp_adapt{i} = b_data_adapt_dfof_all_adapt{i}(:,b_adapt_resp_ind_all{i},:); % adapt response, on adapt trials, with cells sig resposive to adaptor
    b_stim_resp_control{i} = b_data_stim_dfof_all_control{i}(:,b_stim_resp_ind_all{i},:); % stim response, on control trials, with cells sig resposive to control stimulus
    b_stim_resp_adapt{i} = b_data_stim_dfof_all_adapt{i}(:,b_stim_resp_ind_all{i},:); % stim response, on adapt trials, with cells sig resposive to control stimulus

    p_adapt_resp_adapt{i} = p_data_adapt_dfof_all_adapt{i}(:,p_adapt_resp_ind_all{i},:);
    p_stim_resp_control{i} = p_data_stim_dfof_all_control{i}(:,p_stim_resp_ind_all{i},:);
    p_stim_resp_adapt{i} = p_data_stim_dfof_all_adapt{i}(:,p_stim_resp_ind_all{i},:);

    %Responsivity + tuning 

    for bin = 1:length(tuned_cells(1,:))
        b_tuned_adapt_resp_adapt{i, bin} = b_data_adapt_dfof_all_adapt{i}(:,b_tuned_adapt_resp_ind{bin},:); % adapt response, on adapt trials, with TUNED cells sig resposive to adaptor
        b_tuned_stim_resp_control{i, bin} = b_data_stim_dfof_all_control{i}(:,b_tuned_stim_resp_ind{bin},:); % stim response, on control trials, with TUNED cells sig resposive to control stimulus
        b_tuned_stim_resp_adapt{i, bin} = b_data_stim_dfof_all_adapt{i}(:,b_tuned_stim_resp_ind{bin},:); % stim response, on adapt trials, with TUNED cells sig resposive to control stimulus
    
        p_tuned_adapt_resp_adapt{i, bin} = p_data_adapt_dfof_all_adapt{i}(:,p_tuned_adapt_resp_ind{bin},:);
        p_tuned_stim_resp_control{i, bin} = p_data_stim_dfof_all_control{i}(:,p_tuned_stim_resp_ind{bin},:);
        p_tuned_stim_resp_adapt{i, bin} = p_data_stim_dfof_all_adapt{i}(:,p_tuned_stim_resp_ind{bin},:);
    end

   
end

%% Pool mean cell activity across trials!
% Remember that experiments are at different depths...
% Uses local helper function; 'Run Section' to use.

% All by condition

b_data_adapt_trial_mean = mean_cell_resp_by_trial(b_data_adapt_dfof_all_adapt);
b_data_stim_control_trial_mean = mean_cell_resp_by_trial(b_data_stim_dfof_all_control);
b_data_stim_adapt_trial_mean = mean_cell_resp_by_trial(b_data_stim_dfof_all_adapt);

p_data_adapt_trial_mean = mean_cell_resp_by_trial(p_data_adapt_dfof_all_adapt);
p_data_stim_control_trial_mean = mean_cell_resp_by_trial(p_data_stim_dfof_all_control);
p_data_stim_adapt_trial_mean = mean_cell_resp_by_trial(p_data_stim_dfof_all_adapt);

%% Responsive by condition

b_data_sig_adapt_trial_mean = mean_cell_resp_by_trial(b_adapt_resp_adapt);
b_data_sig_stim_control_trial_mean = mean_cell_resp_by_trial(b_stim_resp_control);
b_data_sig_stim_adapt_trial_mean = mean_cell_resp_by_trial(b_stim_resp_adapt);

p_data_sig_adapt_trial_mean = mean_cell_resp_by_trial(p_adapt_resp_adapt);
p_data_sig_stim_control_trial_mean = mean_cell_resp_by_trial(p_stim_resp_control);
p_data_sig_stim_adapt_trial_mean = mean_cell_resp_by_trial(p_stim_resp_adapt);

% Responsive + tuned by condition 

%% Responsive + tuned by condition 

b_data_sig_tuned_adapt_trial_mean = mean_tuned_cell_resp_by_trial(b_tuned_adapt_resp_adapt);
b_data_sig_tuned_stim_control_trial_mean = mean_tuned_cell_resp_by_trial(b_tuned_stim_resp_control);
b_data_sig_tuned_stim_adapt_trial_mean = mean_tuned_cell_resp_by_trial(b_tuned_stim_resp_adapt);

p_data_sig_tuned_adapt_trial_mean = mean_tuned_cell_resp_by_trial(p_tuned_adapt_resp_adapt);
p_data_sig_tuned_stim_control_trial_mean = mean_tuned_cell_resp_by_trial(p_tuned_stim_resp_control);
p_data_sig_tuned_stim_adapt_trial_mean = mean_tuned_cell_resp_by_trial(p_tuned_stim_resp_adapt);

%% Get pooled cell counts for each condition before averaging across cell

nAllCells = size(b_data_adapt_trial_mean,2);

nSigAdapt_b = size(b_data_sig_adapt_trial_mean,2);
nSigStim_control_b = size(b_data_sig_stim_control_trial_mean,2);
nSigStim_adapt_b = size(b_data_sig_stim_adapt_trial_mean,2);

nSigAdapt_p = size(p_data_sig_adapt_trial_mean,2);
nSigStim_control_p = size(p_data_sig_stim_control_trial_mean,2);
nSigStim_adapt_p = size(p_data_sig_stim_adapt_trial_mean,2);

for i = 1:length(ori_bins)
    nSigTunedAdapt_b(i) = size(b_data_sig_tuned_adapt_trial_mean{i}, 2);
    nSigTunedStim_control_b(i) = size(b_data_sig_tuned_stim_control_trial_mean{i}, 2);
    nSigTunedStim_adapt_b(i) = size(b_data_sig_tuned_stim_adapt_trial_mean{i}, 2);

    nSigTunedAdapt_p(i) = size(p_data_sig_tuned_adapt_trial_mean{i}, 2);
    nSigTunedStim_control_p(i) = size(p_data_sig_tuned_stim_control_trial_mean{i}, 2);
    nSigTunedStim_adapt_p(i) = size(p_data_sig_tuned_stim_adapt_trial_mean{i}, 2);
end




%% Cell average dF/F and SEM for each condition





%% 
clearvars -except adaptor_vline adaptPeriodMsb_adapt_resp_ind_all ...
    b_adapt_trial_ind_all b_control_trial_ind_all b_data_adapt_dfof_all ...
    b_data_stim_dfof_all b_mask_label_all b_stim_resp_ind_all base_win ...
    base_win_all dates dynAdaptFlashOff dynAdaptPeriodMs expt_list_final ...
    mouse ori_bins p_adapt_resp_ind_all p_adapt_trial_ind_all ...
    p_control_trial_ind_all p_data_adapt_dfof_all p_data_stim_dfof_all ...
    p_mask_label_all p_stim_resp_ind_all resp_win resp_win_all ...
    target_vline tuned_cells b_data_adapt_dfof_all_adapt ...
    b_data_stim_dfof_all_control b_data_stim_dfof_all_adapt ...
    p_data_adapt_dfof_all_adapt p_data_stim_dfof_all_control ...
    p_data_stim_dfof_all_adapt

%%
cd(['Z:\All_Staff\home\camaron\Analysis\2P\' mouse])
filename = [mouse '_pooled_2AFC_2P_data.mat'];
save(filename)


%% local functions

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

