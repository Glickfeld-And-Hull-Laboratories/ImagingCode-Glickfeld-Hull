clearvars
close all
clc

%% Load expt list
dataset = 'oriAdapt_V1';
eval(dataset);

mouse = 'i472'; % indicate mouse data to pool
cd(['Z:\All_Staff\home\camaron\Analysis\2P\' mouse])
filename = [mouse '_image_matching.mat'];
load(filename, "expt_list_final", "z_pos_good_list")

% find zpos for this final list, and save for later
[C, ia, ib] = intersect(z_pos_good_list(1,:), expt_list_final); 
z_pos_expt_list_final = z_pos_good_list(:,ia);


% Break data into trial conditions and get indices of cells of interest 
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

   %% find trials for each target - for target analysis 

    tGratingOri_b = b_stim_data.tGratingOri;
    tGratingOri_p = p_stim_data.tGratingOri;

    tOris_b = unique(tGratingOri_b);
    tOris_p = unique(tGratingOri_p);

        
    if ~isequal(tOris_b, tOris_p) 

        tOris = intersect(tOris_b,tOris_p);
    end

    
    
    %     Ori_trial_ind_b = cell(length(expt_list_final), length(tOris));
    %     Ori_trial_ind_p = cell(length(expt_list_final), length(tOris));
    
    
    for k = 1:length(tOris)
        Ori_trial_ind_b{i,k} = find(tGratingOri_b == (tOris(k)));
        Ori_trial_ind_p{i,k} = find(tGratingOri_p == (tOris(k)));
    
    end

    % if one of the two conditions (active, behaving) is missing target
    % orientations, use only orientations present in both




    
    
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
    b_stim_win_dfof_all_both_con{i} = b_stim_resp.data_stim_dfof; % Expt(i) = [trial_time X nCells X Adapt_Trials]


    % anyalysis window used to determine responsive cells (save these)
%     base_win = b_stim_resp.base_win;
%     resp_win = b_stim_resp.resp_win;
    base_win = [17 18 19 20 21]; % extended base win
    resp_win = [23 24 25 26 27]; % extended resp win
    base_win_all = b_adapt_resp.base_win_all;
    resp_win_all = b_adapt_resp.resp_win_all;

    % find index of responsive cells
    b_ind_sig_resp_adapt_fromTC(i) = b_adapt_resp.adapt_resp_ind; % Uses data from TC_extraction, compare with new ttest
    b_ind_sig_resp_stim_fromTC{i} = b_stim_resp.stim_resp_ind;

    b_ind_sig_resp_adapt_all{i} = find_respCells(b_adapt_win_dfof_all_adapt{i}, base_win, resp_win);
    b_ind_sig_resp_stim_control_all{i} = find_respCells(b_stim_win_dfof_all_control{i}, base_win, resp_win);
    b_ind_sig_resp_stim_adapt_all{i} = find_respCells(b_stim_win_dfof_all_adapt{i}, base_win, resp_win); % 
    b_ind_sig_resp_stim_both_con{i} = find_respCells(b_stim_win_dfof_all_both_con{i}, base_win, resp_win); % 


    [target_data_b(i,:), target_ori_data_adapt_b(i,:), target_ori_data_control_b(i,:)] = split_data_by_target_ori(b_stim_win_dfof_all_both_con{i}, Ori_trial_ind_b(i,:), adapt_trial_ind, control_trial_ind);

    
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
    p_stim_win_dfof_all_both_con{i} = p_stim_resp.data_stim_dfof; % Expt(i) = [trial_time X nCells X Adapt_Trials]

    
    % find index of responsive cells
    p_ind_sig_resp_adapt_fromTC(i) = p_adapt_resp.adapt_resp_ind; % 
    p_ind_sig_resp_stim_fromTC{i} = p_stim_resp.stim_resp_ind;

    p_ind_sig_resp_adapt_all{i} = find_respCells(p_adapt_win_dfof_all_adapt{i}, base_win, resp_win);
    p_ind_sig_resp_stim_control_all{i} = find_respCells(p_stim_win_dfof_all_control{i}, base_win, resp_win);
    p_ind_sig_resp_stim_adapt_all{i} = find_respCells(p_stim_win_dfof_all_adapt{i}, base_win, resp_win); % 
    p_ind_sig_resp_stim_both_con{i} = find_respCells(p_stim_win_dfof_all_both_con{i}, base_win, resp_win); % 


    [target_data_p(i,:), target_ori_data_adapt_p(i,:), target_ori_data_control_p(i,:)] = split_data_by_target_ori(p_stim_win_dfof_all_both_con{i}, Ori_trial_ind_p(i,:), adapt_trial_ind, control_trial_ind);

 
%%



    
    %% Intersect indices across conditions

   %Find cells responsive target during adapt and control conditions (do
    %for passive as well), then intersection across conditions (behaving, passive)
    

%     b_ind_sig_resp_adapt_stim_union = union(b_ind_sig_resp_adapt_all{i}, b_ind_sig_resp_stim_both_con{i});
%     
%     p_ind_sig_resp_adapt_stim_union = union(p_ind_sig_resp_adapt_all{i}, p_ind_sig_resp_stim_both_con{i});

    bp_ind_sig_resp_stim_intersect{i} = intersect(b_ind_sig_resp_stim_both_con{i}, p_ind_sig_resp_stim_both_con{i}); % Use intersection

%     bp_ind_sig_resp_adapt_stim_intersect{i} = intersect(b_ind_sig_resp_adapt_stim_union, p_ind_sig_resp_adapt_stim_union); % Use intersection
%     bp_ind_sig_resp_adapt_stim_union{i} = union(b_ind_sig_resp_adapt_stim_union, p_ind_sig_resp_adapt_stim_union);


    %Find intersection of cells responsive to adaptors across conditions
    %(behaving, passive), for Aix, across oris if possible...

    bp_ind_sig_resp_adapt_intersect{i} = intersect(b_ind_sig_resp_adapt_all{i}, p_ind_sig_resp_adapt_all{i}); % Use itersection
%     bp_ind_sig_resp_adapt_union{i} = union(b_ind_sig_resp_adapt_all{i}, p_ind_sig_resp_adapt_all{i});


%%

% count_items(bp_ind_sig_resp_adapt_stim_union)
% count_items(bp_ind_sig_resp_adapt_stim_intersect)
% 
% count_items(bp_ind_sig_resp_adapt_intersect)
% count_items(bp_ind_sig_resp_adapt_union)
% 
% count_items(b_ind_sig_resp_adapt_all)
% count_items(p_ind_sig_resp_adapt_all)
%%
    
    % -----
    
    % Masks (Interneurons)  - Dont forget to use these! (b and p mask label are
    % the same)
    b_mask_label_all{i} =  b_adapt_resp.mask_label;
    p_mask_label_all{i} =  p_adapt_resp.mask_label;

    ind_int{i} = find(b_mask_label_all{i}); 
    ind_pyr{i} = find(~b_mask_label_all{i}); 

    count_int(i) = length(ind_int{i});
    count_pyr(i) = length(ind_pyr{i});
    
    % index of responsive interneurons with masks 

    b_ind_sig_resp_adapt_all_interneurons{i} = intersect(b_ind_sig_resp_adapt_all{i}, ind_int{i});
    b_ind_sig_resp_stim_control_all_interneurons{i} = intersect(b_ind_sig_resp_stim_control_all{i}, ind_int{i});
    b_ind_sig_resp_stim_adapt_all_interneurons{i} = intersect(b_ind_sig_resp_stim_adapt_all{i}, ind_int{i});

    p_ind_sig_resp_adapt_all_interneurons{i} = intersect(p_ind_sig_resp_adapt_all{i}, ind_int{i});
    p_ind_sig_resp_stim_control_all_interneurons{i} = intersect(p_ind_sig_resp_stim_control_all{i}, ind_int{i});
    p_ind_sig_resp_stim_adapt_all_interneurons{i} = intersect(p_ind_sig_resp_stim_adapt_all{i}, ind_int{i});

    %--cells active across conditions 
    bp_ind_sig_resp_stim_intersect_interneurons{i} = intersect(bp_ind_sig_resp_stim_intersect{i}, ind_int{i});
    bp_ind_sig_resp_adapt_intersect_interneurons{i} = intersect(bp_ind_sig_resp_adapt_intersect{i}, ind_int{i});

    %---- pyramidal

    b_ind_sig_resp_adapt_all_pyramidal{i} = intersect(b_ind_sig_resp_adapt_all{i}, ind_pyr{i});
    b_ind_sig_resp_stim_control_all_pyramidal{i} = intersect(b_ind_sig_resp_stim_control_all{i}, ind_pyr{i});
    b_ind_sig_resp_stim_adapt_all_pyramidal{i} = intersect(b_ind_sig_resp_stim_adapt_all{i}, ind_pyr{i});

    p_ind_sig_resp_adapt_all_pyramidal{i} = intersect(p_ind_sig_resp_adapt_all{i}, ind_pyr{i});
    p_ind_sig_resp_stim_control_all_pyramidal{i} = intersect(p_ind_sig_resp_stim_control_all{i}, ind_pyr{i});
    p_ind_sig_resp_stim_adapt_all_pyramidal{i} = intersect(p_ind_sig_resp_stim_adapt_all{i}, ind_pyr{i});

    %--cells active across conditions 
    bp_ind_sig_resp_stim_intersect_pyramidal{i} = intersect(bp_ind_sig_resp_stim_intersect{i}, ind_pyr{i});
    bp_ind_sig_resp_adapt_intersect_pyramidal{i} = intersect(bp_ind_sig_resp_adapt_intersect{i}, ind_pyr{i});



   
   
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

        bp_ind_tuned_sig_resp_stim_intersect{bin} = intersect(bp_ind_sig_resp_stim_intersect{i}, tuned_cells{i, bin});
        bp_ind_tuned_sig_resp_adapt_intersect{bin} = intersect(bp_ind_sig_resp_adapt_intersect{i}, tuned_cells{i, bin});

    end

    b_ind_tuned_sig_resp_adapt_all{i} = b_ind_tuned_sig_resp_adapt;
    b_ind_tuned_sig_resp_stim_control_all{i} = b_ind_tuned_sig_resp_stim_control;
    b_ind_tuned_sig_resp_stim_adapt_all{i} = b_ind_tuned_sig_resp_stim_adapt;


    p_ind_tuned_sig_resp_adapt_all{i} = p_ind_tuned_sig_resp_adapt;
    p_ind_tuned_sig_resp_stim_control_all{i} = p_ind_tuned_sig_resp_stim_control;
    p_ind_tuned_sig_resp_stim_adapt_all{i} = p_ind_tuned_sig_resp_stim_adapt;

    bp_ind_tuned_sig_resp_stim_intersect_all = bp_ind_tuned_sig_resp_stim_intersect;
    bp_ind_tuned_sig_resp_adapt_intersect_all = bp_ind_tuned_sig_resp_adapt_intersect;




%% Segment cells from conditional dF/F by responsivitiy and tuning




    %Responsivity 
    b_adapt_sigDuringAdapt_adapt{i} = b_adapt_win_dfof_all_adapt{i}(:,b_ind_sig_resp_adapt_all{i},:); % adapt response, on adapt trials, with cells sig responsive during adaptorON condition

    b_stim_sigDuringControl_control{i} = b_stim_win_dfof_all_control{i}(:,b_ind_sig_resp_stim_control_all{i},:); % stim response, on control trials, with cells sig during control condition
    b_stim_sigDuringControl_adapt{i} = b_stim_win_dfof_all_adapt{i}(:,b_ind_sig_resp_stim_control_all{i},:); % stim response, on adapt trials, with cells sig resposive during control condition

    b_stim_sigDuringAdapt_control{i} = b_stim_win_dfof_all_control{i}(:,b_ind_sig_resp_stim_adapt_all{i},:); % stim response, on control trials, with cells sig responsive during adaptorON condition
    b_stim_sigDuringAdapt_adapt{i} = b_stim_win_dfof_all_adapt{i}(:,b_ind_sig_resp_stim_adapt_all{i},:); % stim response, on adapt trials, with cells sig responsive during adaptorON condition

    % ----

    p_adapt_sigDuringAdapt_adapt{i} = p_adapt_win_dfof_all_adapt{i}(:,p_ind_sig_resp_adapt_all{i},:);

    p_stim_sigDuringControl_control{i} = p_stim_win_dfof_all_control{i}(:,p_ind_sig_resp_stim_control_all{i},:);
    p_stim_sigDuringControl_adapt{i} = p_stim_win_dfof_all_adapt{i}(:,p_ind_sig_resp_stim_control_all{i},:);

    p_stim_sigDuringAdapt_control{i} = p_stim_win_dfof_all_control{i}(:,p_ind_sig_resp_stim_adapt_all{i},:); % stim response, on control trials, with cells sig responsive during adaptorON condition
    p_stim_sigDuringAdapt_adapt{i} = p_stim_win_dfof_all_adapt{i}(:,p_ind_sig_resp_stim_adapt_all{i},:); % stim response, on adapt trials, with cells sig responsive during adaptorON condition



    % Responsivity + tuning 

    for bin = 1:length(tuned_cells(1,:))
        b_tuned_adapt_sigDuringAdapt_adapt{i, bin} = b_adapt_win_dfof_all_adapt{i}(:,b_ind_tuned_sig_resp_adapt{bin},:); % adapt response, on adapt trials, with TUNED cells sig responsive during adaptorON condition

        b_tuned_stim_sigDuringControl_control{i, bin} = b_stim_win_dfof_all_control{i}(:,b_ind_tuned_sig_resp_stim_control{bin},:); % stim response, on control trials, with TUNED cells sig responsive during control condition
        b_tuned_stim_sigDuringControl_adapt{i, bin} = b_stim_win_dfof_all_adapt{i}(:,b_ind_tuned_sig_resp_stim_control{bin},:); % stim response, on adapt trials, with TUNED cells sig responsive during control condition

        b_tuned_stim_sigDuringAdapt_control{i, bin} = b_stim_win_dfof_all_control{i}(:,b_ind_tuned_sig_resp_stim_adapt{bin},:); % stim response, on control trials, with TUNED cells sig responsive during adaptorON condition
        b_tuned_stim_sigDuringAdapt_adapt{i, bin} = b_stim_win_dfof_all_adapt{i}(:,b_ind_tuned_sig_resp_stim_adapt{bin},:); % stim response, on adapt trials, with TUNED cells sig responsive during adaptorON condition


        p_tuned_adapt_sigDuringAdapt_adapt{i, bin} = p_adapt_win_dfof_all_adapt{i}(:,p_ind_tuned_sig_resp_adapt{bin},:);

        p_tuned_stim_sigDuringControl_control{i, bin} = p_stim_win_dfof_all_control{i}(:,p_ind_tuned_sig_resp_stim_control{bin},:);
        p_tuned_stim_sigDuringControl_adapt{i, bin} = p_stim_win_dfof_all_adapt{i}(:,p_ind_tuned_sig_resp_stim_control{bin},:);

        p_tuned_stim_sigDuringAdapt_control{i, bin} = p_stim_win_dfof_all_control{i}(:,p_ind_tuned_sig_resp_stim_adapt{bin},:); % stim response, on control trials, with TUNED cells sig responsive during adaptorON condition
        p_tuned_stim_sigDuringAdapt_adapt{i, bin} = p_stim_win_dfof_all_adapt{i}(:,p_ind_tuned_sig_resp_stim_adapt{bin},:); % stim response, on adapt trials, with TUNED cells sig responsive during adaptorON condition

    end

    % Responsivity + cell class (interneuron)
    b_adapt_sigDuringAdapt_adapt_interneurons{i} = b_adapt_win_dfof_all_adapt{i}(:,b_ind_sig_resp_adapt_all_interneurons{i},:); % adapt response, on adapt trials, with cells sig responsive during adaptorON condition

    b_stim_sigDuringControl_control_interneurons{i} = b_stim_win_dfof_all_control{i}(:,b_ind_sig_resp_stim_control_all_interneurons{i},:); % stim response, on control trials, with cells sig during control condition
    b_stim_sigDuringControl_adapt_interneurons{i} = b_stim_win_dfof_all_adapt{i}(:,b_ind_sig_resp_stim_control_all_interneurons{i},:); % stim response, on adapt trials, with cells sig resposive during control condition

    b_stim_sigDuringAdapt_control_interneurons{i} = b_stim_win_dfof_all_control{i}(:,b_ind_sig_resp_stim_adapt_all_interneurons{i},:); % stim response, on control trials, with cells sig responsive during adaptorON condition
    b_stim_sigDuringAdapt_adapt_interneurons{i} = b_stim_win_dfof_all_adapt{i}(:,b_ind_sig_resp_stim_adapt_all_interneurons{i},:); % stim response, on adapt trials, with cells sig responsive during adaptorON condition
    % ----
    p_adapt_sigDuringAdapt_adapt_interneurons{i} = p_adapt_win_dfof_all_adapt{i}(:,p_ind_sig_resp_adapt_all_interneurons{i},:);

    p_stim_sigDuringControl_control_interneurons{i} = p_stim_win_dfof_all_control{i}(:,p_ind_sig_resp_stim_control_all_interneurons{i},:);
    p_stim_sigDuringControl_adapt_interneurons{i} = p_stim_win_dfof_all_adapt{i}(:,p_ind_sig_resp_stim_control_all_interneurons{i},:);

    p_stim_sigDuringAdapt_control_interneurons{i} = p_stim_win_dfof_all_control{i}(:,p_ind_sig_resp_stim_adapt_all_interneurons{i},:); % stim response, on control trials, with cells sig responsive during adaptorON condition
    p_stim_sigDuringAdapt_adapt_interneurons{i} = p_stim_win_dfof_all_adapt{i}(:,p_ind_sig_resp_stim_adapt_all_interneurons{i},:); % stim response, on adapt trials, with cells sig responsive during adaptorON condition

    % Responsivity + cell class (pyramidal)
    b_adapt_sigDuringAdapt_adapt_pyramidal{i} = b_adapt_win_dfof_all_adapt{i}(:,b_ind_sig_resp_adapt_all_pyramidal{i},:); % adapt response, on adapt trials, with cells sig responsive during adaptorON condition

    b_stim_sigDuringControl_control_pyramidal{i} = b_stim_win_dfof_all_control{i}(:,b_ind_sig_resp_stim_control_all_pyramidal{i},:); % stim response, on control trials, with cells sig during control condition
    b_stim_sigDuringControl_adapt_pyramidal{i} = b_stim_win_dfof_all_adapt{i}(:,b_ind_sig_resp_stim_control_all_pyramidal{i},:); % stim response, on adapt trials, with cells sig resposive during control condition

    b_stim_sigDuringAdapt_control_pyramidal{i} = b_stim_win_dfof_all_control{i}(:,b_ind_sig_resp_stim_adapt_all_pyramidal{i},:); % stim response, on control trials, with cells sig responsive during adaptorON condition
    b_stim_sigDuringAdapt_adapt_pyramidal{i} = b_stim_win_dfof_all_adapt{i}(:,b_ind_sig_resp_stim_adapt_all_pyramidal{i},:); % stim response, on adapt trials, with cells sig responsive during adaptorON condition
    % ----
    p_adapt_sigDuringAdapt_adapt_pyramidal{i} = p_adapt_win_dfof_all_adapt{i}(:,p_ind_sig_resp_adapt_all_pyramidal{i},:);

    p_stim_sigDuringControl_control_pyramidal{i} = p_stim_win_dfof_all_control{i}(:,p_ind_sig_resp_stim_control_all_pyramidal{i},:);
    p_stim_sigDuringControl_adapt_pyramidal{i} = p_stim_win_dfof_all_adapt{i}(:,p_ind_sig_resp_stim_control_all_pyramidal{i},:);

    p_stim_sigDuringAdapt_control_pyramidal{i} = p_stim_win_dfof_all_control{i}(:,p_ind_sig_resp_stim_adapt_all_pyramidal{i},:); % stim response, on control trials, with cells sig responsive during adaptorON condition
    p_stim_sigDuringAdapt_adapt_pyramidal{i} = p_stim_win_dfof_all_adapt{i}(:,p_ind_sig_resp_stim_adapt_all_pyramidal{i},:); % stim response, on adapt trials, with cells sig responsive during adaptorON condition


%% NEW DATA SEGEMENTS (Intersect Passive and behaving responsive cells)

    %Responsivity

    %Adapt adapt (adapt intersect)
    b_adapt{i} = b_adapt_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_adapt_intersect{i},:);
    p_adapt{i} = p_adapt_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_adapt_intersect{i},:);

    
    %Target Adapt (target intersect)
    b_targ_adapt{i} = b_stim_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_stim_intersect{i},:);
    p_targ_adapt{i} = p_stim_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_stim_intersect{i},:);

    %Target Control (target intersect)
    b_targ_control{i} = b_stim_win_dfof_all_control{i}(:,bp_ind_sig_resp_stim_intersect{i},:);
    p_targ_control{i} = p_stim_win_dfof_all_control{i}(:,bp_ind_sig_resp_stim_intersect{i},:);
   
%     %AIX (adapt intersect)
%     b_adapt_AIX{i} = b_adapt_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_adapt_intersect{i},:);
%     p_adapt_AIX{i} = p_adapt_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_adapt_intersect{i},:);
    

    %Responsivity + tuning (few cells...)
    for bin = 1:length(tuned_cells(1,:))


        %Adapt adapt (all intersect for amplitude)
        b_adapt_AIX_tuned{i, bin} = b_adapt_win_dfof_all_adapt{i}(:,bp_ind_tuned_sig_resp_adapt_intersect_all{bin},:);
        p_adapt_AIX_tuned{i, bin} = p_adapt_win_dfof_all_adapt{i}(:,bp_ind_tuned_sig_resp_adapt_intersect_all{bin},:);

        
        %Target Adapt (amplitude)
        b_targ_adapt_tuned{i, bin} = b_stim_win_dfof_all_adapt{i}(:,bp_ind_tuned_sig_resp_stim_intersect_all{bin},:);
        p_targ_adapt_tuned{i, bin} = p_stim_win_dfof_all_adapt{i}(:,bp_ind_tuned_sig_resp_stim_intersect_all{bin},:);
    
        %Target Control (amplitude)
        b_targ_control_tuned{i, bin} = b_stim_win_dfof_all_control{i}(:,bp_ind_tuned_sig_resp_stim_intersect_all{bin},:);
        p_targ_control_tuned{i, bin} = p_stim_win_dfof_all_control{i}(:,bp_ind_tuned_sig_resp_stim_intersect_all{bin},:);
       
%         %AIX (adapt intersect)
%         b_adapt_AIX_tuned{i, bin} = b_adapt_win_dfof_all_adapt{i}(:,bp_ind_tuned_sig_resp_adapt_intersect_all{bin},:);
%         p_adapt_AIX_tuned{i, bin} = p_adapt_win_dfof_all_adapt{i}(:,bp_ind_tuned_sig_resp_adapt_intersect_all{bin},:);

    end

    %Responsivity + cell class (interneuron)

     %Adapt adapt (all intersect for amplitude)
    b_adapt_int{i} = b_adapt_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_adapt_intersect_interneurons{i},:);
    p_adapt_int{i} = p_adapt_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_adapt_intersect_interneurons{i},:);
    
    
    %Target Adapt (amplitude)
    b_targ_adapt_int{i} = b_stim_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_stim_intersect_interneurons{i},:);
    p_targ_adapt_int{i} = p_stim_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_stim_intersect_interneurons{i},:);

    %Target Control (amplitude)
    b_targ_control_int{i} = b_stim_win_dfof_all_control{i}(:,bp_ind_sig_resp_stim_intersect_interneurons{i},:);
    p_targ_control_int{i} = p_stim_win_dfof_all_control{i}(:,bp_ind_sig_resp_stim_intersect_interneurons{i},:);
   
%     %AIX (adapt intersect)
%     b_adapt_AIX_int{i} = b_adapt_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_adapt_intersect_interneurons{i},:);
%     p_adapt_AIX_int{i} = p_adapt_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_adapt_intersect_interneurons{i},:);
%     
    

    %Responsivity + cell class (pyramidal)

    %Adapt adapt (all intersect for amplitude)
    b_adapt_pyr{i} = b_adapt_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_adapt_intersect_pyramidal{i},:);
    p_adapt_pyr{i} = p_adapt_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_adapt_intersect_pyramidal{i},:);

    
    %Target Adapt (amplitude)
    b_targ_adapt_pyr{i} = b_stim_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_stim_intersect_pyramidal{i},:);
    p_targ_adapt_pyr{i} = p_stim_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_stim_intersect_pyramidal{i},:);

    %Target Control (amplitude)
    b_targ_control_pyr{i} = b_stim_win_dfof_all_control{i}(:,bp_ind_sig_resp_stim_intersect_pyramidal{i},:);
    p_targ_control_pyr{i} = p_stim_win_dfof_all_control{i}(:,bp_ind_sig_resp_stim_intersect_pyramidal{i},:);
   
%     %AIX (adapt intersect)
%     b_adapt_AIX_pyr{i} = b_adapt_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_adapt_intersect_pyramidal{i},:);
%     p_adapt_AIX_pyr{i} = p_adapt_win_dfof_all_adapt{i}(:,bp_ind_sig_resp_adapt_intersect_pyramidal{i},:);


   


end
adaptor_vline = [20 31 42 53];
target_vline = 64;

%% Pool mean cell activity across trials
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

b_stim_sigDuringAdapt_control_trial_mean = mean_cell_resp_by_trial(b_stim_sigDuringAdapt_control);
b_stim_sigDuringAdapt_adapt_trial_mean = mean_cell_resp_by_trial(b_stim_sigDuringAdapt_adapt);


p_adapt_sigDuringAdapt_adapt_trial_mean = mean_cell_resp_by_trial(p_adapt_sigDuringAdapt_adapt);

p_stim_sigDuringControl_control_trial_mean = mean_cell_resp_by_trial(p_stim_sigDuringControl_control);
p_stim_sigDuringControl_adapt_trial_mean = mean_cell_resp_by_trial(p_stim_sigDuringControl_adapt);

p_stim_sigDuringAdapt_control_trial_mean = mean_cell_resp_by_trial(p_stim_sigDuringAdapt_control);
p_stim_sigDuringAdapt_adapt_trial_mean = mean_cell_resp_by_trial(p_stim_sigDuringAdapt_adapt);

% Responsive + tuned by condition 

b_tuned_adapt_sigDuringAdapt_adapt_trial_mean = mean_tuned_cell_resp_by_trial(b_tuned_adapt_sigDuringAdapt_adapt);

b_tuned_stim_sigDuringControl_control_trial_mean = mean_tuned_cell_resp_by_trial(b_tuned_stim_sigDuringControl_control);
b_tuned_stim_sigDuringControl_adapt_trial_mean = mean_tuned_cell_resp_by_trial(b_tuned_stim_sigDuringControl_adapt);

b_tuned_stim_sigDuringAdapt_control_trial_mean = mean_tuned_cell_resp_by_trial(b_tuned_stim_sigDuringAdapt_control);
b_tuned_stim_sigDuringAdapt_adapt_trial_mean = mean_tuned_cell_resp_by_trial(b_tuned_stim_sigDuringAdapt_adapt);


p_tuned_adapt_sigDuringAdapt_adapt_trial_mean = mean_tuned_cell_resp_by_trial(p_tuned_adapt_sigDuringAdapt_adapt);

p_tuned_stim_sigDuringControl_control_trial_mean = mean_tuned_cell_resp_by_trial(p_tuned_stim_sigDuringControl_control);
p_tuned_stim_sigDuringControl_adapt_trial_mean = mean_tuned_cell_resp_by_trial(p_tuned_stim_sigDuringControl_adapt);

p_tuned_stim_sigDuringAdapt_control_trial_mean = mean_tuned_cell_resp_by_trial(p_tuned_stim_sigDuringAdapt_control);
p_tuned_stim_sigDuringAdapt_adapt_trial_mean = mean_tuned_cell_resp_by_trial(p_tuned_stim_sigDuringAdapt_adapt);



% Responsive + cell class (interneuron)
b_adapt_sigDuringAdapt_adapt_interneurons_trial_mean = mean_cell_resp_by_trial(b_adapt_sigDuringAdapt_adapt_interneurons);

b_stim_sigDuringControl_control_interneurons_trial_mean = mean_cell_resp_by_trial(b_stim_sigDuringControl_control_interneurons);
b_stim_sigDuringControl_adapt_interneurons_trial_mean = mean_cell_resp_by_trial(b_stim_sigDuringControl_adapt_interneurons);

b_stim_sigDuringAdapt_control_interneurons_trial_mean = mean_cell_resp_by_trial(b_stim_sigDuringAdapt_control_interneurons);
b_stim_sigDuringAdapt_adapt_interneurons_trial_mean = mean_cell_resp_by_trial(b_stim_sigDuringAdapt_adapt_interneurons);
% ----
p_adapt_sigDuringAdapt_adapt_interneurons_trial_mean = mean_cell_resp_by_trial(p_adapt_sigDuringAdapt_adapt_interneurons);

p_stim_sigDuringControl_control_interneurons_trial_mean = mean_cell_resp_by_trial(p_stim_sigDuringControl_control_interneurons);
p_stim_sigDuringControl_adapt_interneurons_trial_mean = mean_cell_resp_by_trial(p_stim_sigDuringControl_adapt_interneurons);

p_stim_sigDuringAdapt_control_interneurons_trial_mean = mean_cell_resp_by_trial(p_stim_sigDuringAdapt_control_interneurons);
p_stim_sigDuringAdapt_adapt_interneurons_trial_mean = mean_cell_resp_by_trial(p_stim_sigDuringAdapt_adapt_interneurons);

% Responsive + cell class (pyramidal)
b_adapt_sigDuringAdapt_adapt_pyramidal_trial_mean = mean_cell_resp_by_trial(b_adapt_sigDuringAdapt_adapt_pyramidal);

b_stim_sigDuringControl_control_pyramidal_trial_mean = mean_cell_resp_by_trial(b_stim_sigDuringControl_control_pyramidal);
b_stim_sigDuringControl_adapt_pyramidal_trial_mean = mean_cell_resp_by_trial(b_stim_sigDuringControl_adapt_pyramidal);

b_stim_sigDuringAdapt_control_pyramidal_trial_mean = mean_cell_resp_by_trial(b_stim_sigDuringAdapt_control_pyramidal);
b_stim_sigDuringAdapt_adapt_pyramidal_trial_mean = mean_cell_resp_by_trial(b_stim_sigDuringAdapt_adapt_pyramidal);
% ----
p_adapt_sigDuringAdapt_adapt_pyramidal_trial_mean = mean_cell_resp_by_trial(p_adapt_sigDuringAdapt_adapt_pyramidal);

p_stim_sigDuringControl_control_pyramidal_trial_mean = mean_cell_resp_by_trial(p_stim_sigDuringControl_control_pyramidal);
p_stim_sigDuringControl_adapt_pyramidal_trial_mean = mean_cell_resp_by_trial(p_stim_sigDuringControl_adapt_pyramidal);

p_stim_sigDuringAdapt_control_pyramidal_trial_mean = mean_cell_resp_by_trial(p_stim_sigDuringAdapt_control_pyramidal);
p_stim_sigDuringAdapt_adapt_pyramidal_trial_mean = mean_cell_resp_by_trial(p_stim_sigDuringAdapt_adapt_pyramidal);


% Reminder that experiments (cells) are at different depths...Responsive by depth

% use mean_cell_resp with tc_stats to find new segemnts...


%% Get pooled cell counts for each condition before averaging across cell (uneccessary with TC_stats function)

nAllCells = size(b_data_adapt_trial_mean,2);
nInt = sum(count_int);
nPyr = sum(count_pyr);

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
b_data_adapt_stats = TC_stats(b_data_adapt_trial_mean, base_win, resp_win);
b_data_stim_control_stats  = TC_stats(b_data_stim_control_trial_mean, base_win, resp_win);
b_data_stim_adapt_stats  = TC_stats(b_data_stim_adapt_trial_mean, base_win, resp_win);

p_data_adapt_stats  = TC_stats(p_data_adapt_trial_mean, base_win, resp_win);
p_data_stim_control_stats  = TC_stats(p_data_stim_control_trial_mean, base_win, resp_win);
p_data_stim_adapt_stats  = TC_stats(p_data_stim_adapt_trial_mean, base_win, resp_win);

% Responsive by condition

b_adapt_sigDuringAdapt_adapt_stats = TC_stats(b_adapt_sigDuringAdapt_adapt_trial_mean, base_win, resp_win);

b_stim_sigDuringControl_control_stats  = TC_stats(b_stim_sigDuringControl_control_trial_mean, base_win, resp_win);
b_stim_sigDuringControl_adapt_stats  = TC_stats(b_stim_sigDuringControl_adapt_trial_mean, base_win, resp_win);

b_stim_sigDuringAdapt_control_stats  = TC_stats(b_stim_sigDuringAdapt_control_trial_mean, base_win, resp_win);
b_stim_sigDuringAdapt_adapt_stats  = TC_stats(b_stim_sigDuringAdapt_adapt_trial_mean, base_win, resp_win);
 
p_adapt_sigDuringAdapt_adapt_stats  = TC_stats(p_adapt_sigDuringAdapt_adapt_trial_mean, base_win, resp_win);

p_stim_sigDuringControl_control_stats  = TC_stats(p_stim_sigDuringControl_control_trial_mean, base_win, resp_win);
p_stim_sigDuringControl_adapt_stats  = TC_stats(p_stim_sigDuringControl_adapt_trial_mean, base_win, resp_win);

p_stim_sigDuringAdapt_control_stats  = TC_stats(p_stim_sigDuringAdapt_control_trial_mean, base_win, resp_win);
p_stim_sigDuringAdapt_adapt_stats  = TC_stats(p_stim_sigDuringAdapt_adapt_trial_mean, base_win, resp_win);

% Responsive + tuned by condition 
for i = 1:4

    b_tuned_adapt_sigDuringAdapt_adapt_stats(i)  = TC_stats(b_tuned_adapt_sigDuringAdapt_adapt_trial_mean{i}, base_win, resp_win);

    b_tuned_stim_sigDuringControl_control_stats(i)  = TC_stats(b_tuned_stim_sigDuringControl_control_trial_mean{i}, base_win, resp_win);
    b_tuned_stim_sigDuringControl_adapt_stats(i)  = TC_stats(b_tuned_stim_sigDuringControl_adapt_trial_mean{i}, base_win, resp_win);
       
    b_tuned_stim_sigDuringAdapt_control_stats(i)  = TC_stats(b_tuned_stim_sigDuringAdapt_control_trial_mean{i}, base_win, resp_win);
    b_tuned_stim_sigDuringAdapt_adapt_stats(i)  = TC_stats(b_tuned_stim_sigDuringAdapt_adapt_trial_mean{i}, base_win, resp_win);
    
    p_tuned_adapt_sigDuringAdapt_adapt_stats(i)  = TC_stats(p_tuned_adapt_sigDuringAdapt_adapt_trial_mean{i}, base_win, resp_win);

    p_tuned_stim_sigDuringControl_control_stats(i)  = TC_stats(p_tuned_stim_sigDuringControl_control_trial_mean{i}, base_win, resp_win);
    p_tuned_stim_sigDuringControl_adapt_stats(i)  = TC_stats(p_tuned_stim_sigDuringControl_adapt_trial_mean{i}, base_win, resp_win);

    p_tuned_stim_sigDuringAdapt_control_stats(i)  = TC_stats(p_tuned_stim_sigDuringAdapt_control_trial_mean{i}, base_win, resp_win);
    p_tuned_stim_sigDuringAdapt_adapt_stats(i)  = TC_stats(p_tuned_stim_sigDuringAdapt_adapt_trial_mean{i}, base_win, resp_win);


end

% Responsive + cell class (interneuron)
b_adapt_sigDuringAdapt_adapt_interneurons_stats = TC_stats(b_adapt_sigDuringAdapt_adapt_interneurons_trial_mean, base_win, resp_win);

b_stim_sigDuringControl_control_interneurons_stats = TC_stats(b_stim_sigDuringControl_control_interneurons_trial_mean, base_win, resp_win);
b_stim_sigDuringControl_adapt_interneurons_stats = TC_stats(b_stim_sigDuringControl_adapt_interneurons_trial_mean, base_win, resp_win);

b_stim_sigDuringAdapt_control_interneurons_stats = TC_stats(b_stim_sigDuringAdapt_control_interneurons_trial_mean, base_win, resp_win);
b_stim_sigDuringAdapt_adapt_interneurons_stats = TC_stats(b_stim_sigDuringAdapt_adapt_interneurons_trial_mean, base_win, resp_win);
% ----
p_adapt_sigDuringAdapt_adapt_interneurons_stats = TC_stats(p_adapt_sigDuringAdapt_adapt_interneurons_trial_mean, base_win, resp_win);

p_stim_sigDuringControl_control_interneurons_stats = TC_stats(p_stim_sigDuringControl_control_interneurons_trial_mean, base_win, resp_win);
p_stim_sigDuringControl_adapt_interneurons_stats = TC_stats(p_stim_sigDuringControl_adapt_interneurons_trial_mean, base_win, resp_win);

p_stim_sigDuringAdapt_control_interneurons_stats = TC_stats(p_stim_sigDuringAdapt_control_interneurons_trial_mean, base_win, resp_win);
p_stim_sigDuringAdapt_adapt_interneurons_stats = TC_stats(p_stim_sigDuringAdapt_adapt_interneurons_trial_mean, base_win, resp_win);

% Responsive + cell class (pyramidal)
b_adapt_sigDuringAdapt_adapt_pyramidal_stats = TC_stats(b_adapt_sigDuringAdapt_adapt_pyramidal_trial_mean, base_win, resp_win);

b_stim_sigDuringControl_control_pyramidal_stats = TC_stats(b_stim_sigDuringControl_control_pyramidal_trial_mean, base_win, resp_win);
b_stim_sigDuringControl_adapt_pyramidal_stats = TC_stats(b_stim_sigDuringControl_adapt_pyramidal_trial_mean, base_win, resp_win);

b_stim_sigDuringAdapt_control_pyramidal_stats = TC_stats(b_stim_sigDuringAdapt_control_pyramidal_trial_mean, base_win, resp_win);
b_stim_sigDuringAdapt_adapt_pyramidal_stats = TC_stats(b_stim_sigDuringAdapt_adapt_pyramidal_trial_mean, base_win, resp_win);
% ----
p_adapt_sigDuringAdapt_adapt_pyramidal_stats = TC_stats(p_adapt_sigDuringAdapt_adapt_pyramidal_trial_mean, base_win, resp_win);

p_stim_sigDuringControl_control_pyramidal_stats = TC_stats(p_stim_sigDuringControl_control_pyramidal_trial_mean, base_win, resp_win);
p_stim_sigDuringControl_adapt_pyramidal_stats = TC_stats(p_stim_sigDuringControl_adapt_pyramidal_trial_mean, base_win, resp_win);

p_stim_sigDuringAdapt_control_pyramidal_stats = TC_stats(p_stim_sigDuringAdapt_control_pyramidal_trial_mean, base_win, resp_win);
p_stim_sigDuringAdapt_adapt_pyramidal_stats = TC_stats(p_stim_sigDuringAdapt_adapt_pyramidal_trial_mean, base_win, resp_win);

%%
% NEW CELL SEGMENT (intersections) STATS

%Adapt adapt (all intersect for amplitude)
b_adapt_stats = TC_stats(mean_cell_resp_by_trial(b_adapt), base_win, resp_win);
p_adapt_stats = TC_stats(mean_cell_resp_by_trial(p_adapt), base_win, resp_win);

%Target Adapt (amplitude)

b_targ_adapt_stats = TC_stats(mean_cell_resp_by_trial(b_targ_adapt), base_win, resp_win);
p_targ_adapt_stats = TC_stats(mean_cell_resp_by_trial(p_targ_adapt), base_win, resp_win);

%Target Control (amplitude)

b_targ_control_stats = TC_stats(mean_cell_resp_by_trial(b_targ_control), base_win, resp_win);
p_targ_control_stats = TC_stats(mean_cell_resp_by_trial(p_targ_control), base_win, resp_win);

% 
% %AIX (adapt intersect)
% 
% b_adapt_AIX_stats = TC_stats(mean_cell_resp_by_trial(b_adapt_AIX), base_win, resp_win);
% p_adapt_AIX_stats = TC_stats(mean_cell_resp_by_trial(p_adapt_AIX), base_win, resp_win);

%%  Apply cell segment and Pool data by target orientation (for target condition)
% [target_data_p(i,:), target_ori_data_adapt_p(i,:), target_ori_data_control_p(i,:)]
% [target_data_b(i,:), target_ori_data_adapt_b(i,:), target_ori_data_control_b(i,:)]
% bp_ind_sig_resp_stim_intersect - index to use

pooled_target_ori_data_all_b = mean_cell_resp_by_trial_target_ori(target_data_b, bp_ind_sig_resp_stim_intersect);
pooled_target_ori_data_all_p = mean_cell_resp_by_trial_target_ori(target_data_p, bp_ind_sig_resp_stim_intersect);

pooled_target_ori_data_adapt_b = mean_cell_resp_by_trial_target_ori(target_ori_data_adapt_b, bp_ind_sig_resp_stim_intersect);
pooled_target_ori_data_adapt_p = mean_cell_resp_by_trial_target_ori(target_ori_data_adapt_p, bp_ind_sig_resp_stim_intersect);

pooled_target_ori_data_control_b = mean_cell_resp_by_trial_target_ori(target_ori_data_control_b, bp_ind_sig_resp_stim_intersect);
pooled_target_ori_data_control_p = mean_cell_resp_by_trial_target_ori(target_ori_data_control_p, bp_ind_sig_resp_stim_intersect);



%% Get stats for target data

for i = 1:length(pooled_target_ori_data_control_b)
    b_target_ori_data_all_stats(i) = TC_stats(pooled_target_ori_data_all_b{i}, base_win, resp_win);
    p_target_ori_data_all_stats(i) = TC_stats(pooled_target_ori_data_all_p{i}, base_win, resp_win);

    b_target_ori_data_adapt_stats(i) = TC_stats(pooled_target_ori_data_adapt_b{i}, base_win, resp_win);
    p_target_ori_data_adapt_stats(i) = TC_stats(pooled_target_ori_data_adapt_p{i}, base_win, resp_win);

    b_target_ori_data_control_stats(i) = TC_stats(pooled_target_ori_data_control_b{i}, base_win, resp_win);
    p_target_ori_data_control_stats(i) = TC_stats(pooled_target_ori_data_control_p{i}, base_win, resp_win);

end

% tOris

%% Most efficient way to clear uneeded variables and save key variables from the workspace?
% Saving everything for now...
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


%% Local functions (most copied to individual scripts)

% function count_items(cell_array)
%     
%     for i = 1:length(cell_array)
%         count(i) = length(cell_array{i});
%     end
%     
%     total = sum(count);
%     disp(total)
% end