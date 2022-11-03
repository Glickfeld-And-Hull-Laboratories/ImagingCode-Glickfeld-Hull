clear all
%%

dataset = 'oriAdapt_V1_cam';
eval(dataset);
iexp = 138;
mouse = expt(iexp).mouse;
date = expt(iexp).date;
user = "lindsey";
condition = 'p';

if user == "lindsey"
    data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey';
elseif user == "camaron"
    data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';
end

run_str = ['runs']; 
run_str = [run_str '-' expt(iexp).pass_run];
load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
% input = pass_input;
load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_adaptResp.mat']))
load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimResp.mat']))

load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
dir_tuning_str = ['runs-' expt(iexp).dirtuning];
load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_tuning_str], [date '_' mouse '_' dir_tuning_str '_oriTuningAndFits.mat']))

adaptor_vline = [20 31 42 53];
target_vline = 64; % variable depending on adapter-target latency ***

nTrials = size(data_adapt_dfof,3)
nCells = size(data_adapt_dfof,2)
 
clear expt
%% Tune cells, place into bins, then get indices to pull tuned cells from TC

[max_resp prefOri] = max(vonMisesFitAllCellsAllBoots,[],1);
prefOri = squeeze(prefOri)-1;
prefOri_bootdiff = abs(prefOri(2:end,:)-prefOri(1,:));
prefOri_bootdiff(find(prefOri_bootdiff>90)) = 180-prefOri_bootdiff(find(prefOri_bootdiff>90));
ind_theta90 = find(prctile(prefOri_bootdiff,90,1)<22.5);

edges = [0 22.5:45:180 180];

[bin edges ind_bin] = histcounts(prefOri(1,:),edges);

ind_bin(find(ind_bin==5)) = 1;
bin(1) = bin(1)+bin(5);
bin(5) = [];

tunedCells = cell(1,length(bin));
for j = 1:length(bin)
    tunedCells{j} = intersect(find(ind_bin==j),ind_theta90);
end

for j = 1:length(tunedCells)
    tunedCells_bin(j) = length(cell2mat(tunedCells(j)));
end

n_cells_tuned = sum(tunedCells_bin)

ori_bins = [0 45 90 135];

clear vonMisesFitAllCellsAllBoots max_resp prefOri R_sqaure prefOri_bootdiff
%% Define indices

%CELLS
adapt_resp_ind = cell2mat(adapt_resp_ind); % cells responsive to adaptor 

% Tuning [0 45 90 135]
tuned_sig_adapt_resp_ind = cell(1,4);
tuned_sig_target_resp_ind = cell(1,4);
for ituned = 1:length(tunedCells)
    tuned_sig_adapt_resp_ind{ituned} = intersect(adapt_resp_ind, tunedCells{ituned}); % adapt responsive and tuned cells
    tuned_sig_target_resp_ind{ituned} = intersect(stim_resp_ind, tunedCells{ituned}); % target responsive and tuned cells
end

%TRIALS
ind_con = find(aGratingContrast == 1);
ind_ori = find(aGratingOri==aOris);
trial_ind_adapt = intersect(ind_ori,ind_con); % trials when adapter is on
trial_ind_control = setdiff(ind_ori,ind_con); % trials when adapter is off (control)

% Control, Adapt trials, by target; (8 targets)
target_ind_all = cell(1,8);
target_ind_control = cell(1,8);
target_ind_adapt = cell(1,8);

for itarget = 1:length(tOris)
    target_ind_all{itarget} = find(tGratingOri == tOris(itarget));
    target_ind_control{itarget} = intersect(target_ind_all{itarget}, trial_ind_control); % control trials by target
    target_ind_adapt{itarget} = intersect(target_ind_all{itarget}, trial_ind_adapt); % adapt trials by target
end

%% Single Cell TCs, Adaptor trials TCs, and Sig responding TCs + Tuned 

% data_TC_total_pop = squeeze(mean(data_adapt_dfof,2, "omitnan")); %
% average of all trials, control and adapt; unnecessary

% data_TC_adapt = mean(data_adapt_dfof(:, :, adapt_trial_ind),3, "omitnan"); % select trials with adapters 

% data_adapt_dfof_adapt_sig = data_adapt_dfof_adapt(:,adapt_resp_ind); % select cells with significant responses to adapters ** Check for erorrs in other script


% for i = 1:length(tunedCells)
%     data_adapt_dfof_adapt_sig_tuned{i} = data_TC_adapt(:, tuned_adapt_resp_ind{i}); % Tuned cells (for each ori_bin)
% end
% 
% for i = 1:length(tunedCells) % Tuned cells ONLY
%     data_adapt_dfof_adapt_tuned_only{i} = data_TC_adapt(:, tunedCells{i}); % Tuned cells (for each ori_bin)
% end

    
%% Population Averages from above (Population TCs, Adaptor trials TCs, and Sig responding TCs) / with SEM

% mean_data_adapt_dfof_total_pop = mean(data_adapt_dfof_total_pop,2, "omitnan");
% mean_data_adapt_dfof_adapt = mean(data_adapt_dfof_adapt,2,"omitnan");
% mean_data_adapt_dfof_adapt_sig = mean(data_adapt_dfof_adapt_sig,2,"omitnan");
% 
% sem_data_adapt_dfof_total_pop = std(data_adapt_dfof_total_pop,[],2, "omitnan") / sqrt(length(data_adapt_dfof_total_pop));
% sem_data_adapt_dfof_adapt = std(data_adapt_dfof_adapt,[],2,"omitnan") / sqrt(length(data_adapt_dfof_adapt));
% sem_data_adapt_dfof_adapt_sig = std(data_adapt_dfof_adapt_sig,[],2,"omitnan") / sqrt(length(data_adapt_dfof_adapt_sig));
% 
% %Tuned cells
% for i = 1:length(tunedCells)
%     mean_data_adapt_dfof_adapt_tuned{i} = mean(data_adapt_dfof_adapt_tuned{i},2,"omitnan"); % Tuned cells (for each ori_bin)
%     sem_data_adapt_dfof_adapt_tuned{i} = std(data_adapt_dfof_adapt_tuned{i},[],2,"omitnan") / sqrt(length(data_adapt_dfof_adapt_tuned{i}));
% end
% 
%% Response by adapter (adapt trials only) and by target response (control trials only)
 
amp_by_adapt_all_cells = adapt_cyc_resp(:,trial_ind_adapt,:); % select adapt ON trials
% amp_R1 = amp_by_adapt_all_cells(:,:,1); %**

% amp_R1_trial_avg = mean(amp_R1, "omitnan");
% amp_R1_cell_mean = mean(amp_R1_trial_avg);
% amp_R1_cell_sem = std(amp_R1_trial_avg) / sqrt(length(amp_R1_trial_avg));
% 
% amp_R2 = amp_by_adapt_all_cells(:,:,2);
% amp_R3 = amp_by_adapt_all_cells(:,:,3);
% amp_R4 = amp_by_adapt_all_cells(:,:,4);

targ_R_control = stim_resp(:,trial_ind_control);
targ_R_adapt = stim_resp(:,trial_ind_adapt);


%separate by target using target_ind_control
for itarget = 1:length(tOris)
    targ_R_control_by_target{:,itarget} = stim_resp(:, target_ind_control{itarget});
    target_resp_control_trial_avg{:,itarget} = mean(targ_R_control_by_target{itarget}, "omitnan");
    target_resp_control_cell_mean(itarget) = mean(target_resp_control_trial_avg{:,itarget},"omitnan");
    target_resp_control_cell_sem(itarget) =  std(target_resp_control_trial_avg{:,itarget}, "omitnan") / sqrt(length(target_resp_control_trial_avg{:,itarget}));
end

%%

% targ_0_targ_mean = [target_resp_control_cell_mean(1:4), amp_R1_cell_mean, target_resp_control_cell_mean(5:8)];
% targ_0_targ_sem = [target_resp_control_cell_sem(1:4), amp_R1_cell_sem, target_resp_control_cell_sem(5:8)];
% oris = [tOris(1:4) 0 tOris(5:8)];
% 
% figure;
% shadedErrorBar([], targ_0_targ_mean, targ_0_targ_sem)
% xticklabels(oris)


%%
save(fullfile(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\2P\' mouse '_' condition '_post_processing_update.mat']));

% save(fullfile(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\2P\' mouse ' ' condition '_post_processing.mat']))
