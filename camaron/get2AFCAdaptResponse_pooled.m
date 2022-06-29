% get2AFCAdaptResponse_Pooled
clearvars
close all
clc


%%
% Get good extractions with passive and direction tuning
dataset = 'oriAdapt_V1';
eval(dataset);

i1402_expts = 1:16;
i1403_expts = 17:29;
i475_expts = 30:79;
i472_expts = 80:125;
    
exp_list = [];

% Just behavior and direction tuning
for i = i1402_expts
    if expt(i).TCs_extracted == 1 & ~isempty(expt(i).dirtuning) 
        exp_list = [exp_list i];
    end
end

% Inlcude passive runs and direction tuning
% for i = i475_expts
%     if expt(i).TCs_extracted == 1 & ~isempty(expt(i).pass_run) & ~isempty(expt(i).dirtuning) 
%         exp_list = [exp_list i];
%     end
% end



%% Pre-procssing: Point to dataset and generate placeholders for key variables
user = "lindsey";
condition = 'b';

if user == "lindsey"
    data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey';
elseif user == "camaron"
    data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';
end

vonMisesFitAllCellsAllBoots_all = []; % For tuning
data_adapt_dfof_all = cell(1,length(exp_list)); % TC during all adapt window
mean_data_adapt_dfof_all = [];
data_adapt_dfof_adapt_all = []; % for pooled time courses of adapt window during adapt trials
% data_adapt_dfof_adapt_all_tuned = [];% cells that or responsive and also tuned
adapt_cyc_resp_all = []; % for average neural responses to adaptor (recalculate from pooled?)
mask_label_all = []; % for cell type
adapt_trial_ind_cell = []; % for adapt trial index
dates = [];


% adapt_cyc_resp_Ex = [];
% adapt_cyc_resp_In = [];

%% Loading Loop: adapt_cyc_resp_all, vonMisesFitAllCellsAllBoots_all, Pooling RESPONSIVE Cells

for i = 1:length(exp_list)
  %% Import Session info 
    iexp = exp_list(i);
    
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;

    save_loc = ['\Analysis\2P\' expt(iexp).date '_' expt(iexp).mouse];

    nrun = size(expt(iexp).runs,1);
    run_str = ['runs']; 
    for irun = 1:nrun
        run_str = [run_str '-' expt(iexp).runs(irun,:)];
    end
    % load initial mask info to be used later in loop
%     load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']))

 %%  Get Input Struct and adapt response info (necessary?? can load stimData.mat from pass or behave)
    if condition == 'b' % behavior
        for irun = 1:nrun
            run_str = ['runs']; 
            run_str = [run_str '-' expt(iexp).runs(irun,:)];
        end
        if exist(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
            load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
            input = adapt_input;
            load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_adaptResp.mat']))
            load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
        else
            continue
        end
    elseif condition == 'p' % passive
        run_str = ['runs']; 
        run_str = [run_str '-' expt(iexp).pass_run];
        if exist(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
            load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
            input = pass_input;
            load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_adaptResp.mat']))
            load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
        else
            continue
        end
    end
    
    input_temp(i) = input; 
    
    
    %Grab stim timings to check later
    adaptPeriodMs(i) = input.adaptPeriodMs;
    dynAdaptPeriodMs(i) = input.dynAdaptPeriodMs;
    dynAdaptFlashOffMs(i) = input.dynAdaptFlashOffMs;
    dates{i} = date;


    
%% Pool mean TCs and Responses to adaptor 

    aGratingOri = celleqel2mat_padded(input.aGratingDirectionDeg);
    aGratingContrast = celleqel2mat_padded(input.aGratingContrast);
    aOris = unique(aGratingOri);


    naOri = length(aOris);
    ind_con = find(aGratingContrast == 1);

%     for iOri = 1:naOri % loop not needed
    ind_ori = find(aGratingOri==aOris);
    adapt_trial_ind = intersect(ind_ori,ind_con); % trials when adapter is on
    adapt_trial_ind_cell{i} = adapt_trial_ind; 
    control_trial_ind = setdiff(ind_ori,ind_con); % trials when adapter is off
%     end

    %Take cell average responses across trials:
    clearvars data_adapt_dfof_adapt
    data_adapt_dfof_all{i} = data_adapt_dfof;
    mean_data_adapt_dfof = nanmean(mean(data_adapt_dfof,2),3);
    mean_data_adapt_dfof_all = cat(2, mean_data_adapt_dfof_all, mean_data_adapt_dfof); % pool mean tcs from each session from adapt window, (plot to look at response window). Arry is trial time X mean dF/F from experiment
    % Average again for overall mean
    data_adapt_dfof_adapt = nanmean(data_adapt_dfof(:, :, adapt_trial_ind),3); % select trials with adapters (looking back this is confusing,  to clarify data_adapt_dfof is the adapt window, CLM- 4/21/22)
    
    data_adapt_dfof_adapt_all = cat(2, data_adapt_dfof_adapt_all, data_adapt_dfof_adapt(:,adapt_resp_ind{:},:)); % select trials with significant responses to adapters 
    

    adapt_resp_ind_cell{i} = adapt_resp_ind;
    adapt_resp_cell{i} = adapt_cyc_resp(adapt_resp_ind{:},:,:); % adapt_cyc_response loaded from singleChannelTC
    if length(adapt_resp_ind{:})== 1
        adapt_cyc_resp_mean = squeeze(nanmean(adapt_resp_cell{i},2))';
    else
        adapt_cyc_resp_mean = squeeze(nanmean(adapt_resp_cell{i},2));
    end
    
    adapt_cyc_resp_all = cat(1, adapt_cyc_resp_all, adapt_cyc_resp_mean); % all cells responsive to adapter 
    
    
%     adapt_Ex  = find(mask_label(adapt_resp_ind{:}) == 0); 
%     adapt_In = find(mask_label(adapt_resp_ind{:}) == 1);
%     adapt_cyc_resp_Ex = cat(1, adapt_cyc_resp_Ex, adapt_cyc_resp_mean(adapt_Ex,:)); % Just Ex
%     adapt_cyc_resp_In = cat(1, adapt_cyc_resp_In, adapt_cyc_resp_mean(adapt_In,:)); % Just In

    

    mask_label_cell{i} = mask_label;
    mask_label_all = [mask_label_all mask_label(adapt_resp_ind{:})];
    
    
 %%  Load and Pool Cell fits
%     dir_tuning_str = 'runs-003';
    dir_tuning_str = ['runs-' expt(iexp).dirtuning];
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_tuning_str], [date '_' mouse '_' dir_tuning_str '_oriTuningAndFits.mat']))
%     load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_tuning_str], [date '_' mouse '_' dir_tuning_str '_oriTuningInfo.mat']))
    

    vonMisesFitAllCellsAllBoots_cell{i} = vonMisesFitAllCellsAllBoots;
    vonMisesFitAllCellsAllBoots_all = cat(3, vonMisesFitAllCellsAllBoots_all, vonMisesFitAllCellsAllBoots(:,:,adapt_resp_ind{:}));
    
end

n_expt_loaded = length(input_temp)
n_cells_total = size(vonMisesFitAllCellsAllBoots_all, 3)

% save(fullfile(['Y:\camaron\Analysis\2P\pooled_data_figs\' mouse '\' condition '_post_processing.mat']), 'data_adapt_dfof_all', 'mean_data_adapt_dfof_all' ,'-append')

%% Session Alignment for TARGET ANALYSIS (Passive v Behavior)

% CHECK for more than one unique adaptPeriodMS (Days must be separated by
% adaptPeriodMS)
adaptPeriodMs = adaptPeriodMs(adaptPeriodMs ~= 0);
uniqueAdaptPeriodMs = unique(adaptPeriodMs); % drop zero values! 

if length(uniqueAdaptPeriodMs) > 1
    disp('Multiple Adapt Periods found!')
    %Write code here to find the index of sessions for each unique
    %adaptPeriodMs
end



%ADD CHECK for if FB ON

%Count cells(exc,inh) for behave and passive sessions (by ori)


%% Tune cells, place into bins, then pull cells from TC

[max_resp prefOri] = max(vonMisesFitAllCellsAllBoots_all,[],1);
prefOri = squeeze(prefOri)-1;
prefOri_bootdiff = abs(prefOri(2:end,:)-prefOri(1,:));
prefOri_bootdiff(find(prefOri_bootdiff>90)) = 180-prefOri_bootdiff(find(prefOri_bootdiff>90));
ind_theta90 = find(prctile(prefOri_bootdiff,90,1)<22.5);

edges = [0 22.5:45:180 180];

[bin edges ind_bin] = histcounts(prefOri(1,:),edges);

% [bin ind_bin] = histc(prefOri(1,:),edges);
ind_bin(find(ind_bin==5)) = 1;
bin(1) = bin(1)+bin(5);
bin(5) = [];

tunedCells = cell(1,length(bin));
for j = 1:length(bin)
    tunedCells{j} = intersect(find(ind_bin==j),ind_theta90);
end

%     figure;
%     plot(edges,bin, ':o')
%     xticks(edges)
%     hold on

for j = 1:length(tunedCells)
    tunedCells_bin(j) = length(cell2mat(tunedCells(j)));
end

%     plot(edges, tunedCells_bin, ':ok')
%     legend([{'prefOri'}, {'tunedCells'}])
%     title('i472: 210625')
%     xlabel('Oris')
%     ylabel('nCells')


n_cells_tuned = sum(tunedCells_bin)
%%

%Apply mask label to df/f
for i = 1:length(tunedCells)
    ExcCells{i} = intersect(tunedCells{i}, find(mask_label_all == 0));
    InhCells{i} = intersect(tunedCells{i}, find(mask_label_all == 1));
    data_adapt_dfof_adapt_all_Exc{i} = data_adapt_dfof_adapt_all(:, ExcCells{i});
    data_adapt_dfof_adapt_all_Inh{i} = data_adapt_dfof_adapt_all(:, InhCells{i});

    data_adapt_dfof_adapt_all_tuned{i} = data_adapt_dfof_adapt_all(:, tunedCells{i});
end

%Apply tuned cells to adapt_cyc_resp
for i = 1:length(tunedCells)
%     ExcCells = intersect(tunedCells{i}, find(mask_label_all == 0));
%     InhCells = intersect(tunedCells{i}, find(mask_label_all == 1));
    adapt_resp_Exc{i} = adapt_cyc_resp_all(ExcCells{i}, :);
    adapt_resp_Inh{i} = adapt_cyc_resp_all(InhCells{i}, :);
    
    
    adapt_resp_all{i} = adapt_cyc_resp_all(tunedCells{i}, :);
end

% for i = 1:length(tunedCells)
%     adapt_resp_Ex{i} = adapt_cyc_resp_Ex(tunedCells{i},:);
%     adapt_resp_In{i} = adapt_cyc_resp_In(tunedCells{i},:);
% 
% end

%% Time Courses
clearvars TC_mean TC_sem

%Get mean and SEM for both TC's and Adapt Response
for i = 1:length(tunedCells)
    TC_mean_Exc(:,i) = nanmean(data_adapt_dfof_adapt_all_Exc{i}, 2); 
    TC_sem_Exc(:,i) =  std(data_adapt_dfof_adapt_all_Exc{i},1, 2) / sqrt(size(data_adapt_dfof_adapt_all_Exc{i},2));
    
    TC_mean_Inh(:,i) = nanmean(data_adapt_dfof_adapt_all_Inh{i}, 2); 
    TC_sem_Inh(:,i) =  std(data_adapt_dfof_adapt_all_Inh{i},1, 2) / sqrt(size(data_adapt_dfof_adapt_all_Inh{i},2));
    
    TC_mean_all(:,i) = nanmean(data_adapt_dfof_adapt_all_tuned{i}, 2); 
    TC_sem_all(:,i) =  std(data_adapt_dfof_adapt_all_tuned{i},1, 2) / sqrt(size(data_adapt_dfof_adapt_all_tuned{i},2));
end

ori_bins = [0 45 90 135];
adaptor_vline = [20 31 42 53];
target_vline = 64;

%% Plot TC's - 4 plots (4 Oris) *PUT EX and IN on same plot, diff colors (BLUE, RED)
% ori_bins = [0 45 90 135];
% [n n2] = subplotn(length(ori_bins));
% 
% adaptor_vline = [20 31 42 53];
% target_vline = 64;
% 
% figure;
% hold on
% for i = 1:length(ori_bins)
%     subplot(n, n2, i)
% %     errorbar(TC_mean_Exc(:,i), TC_sem_Exc(:,i), 'k');
%     shadedErrorBar([],TC_mean_Exc(:,i), TC_sem_Exc(:,i), 'lineProps','b')
%     hold on
%     vline([20 31 42 53])
%     vline(64, 'b:')
%     title(num2str(ori_bins(i)));
%     ylabel('dF/F')
%     xlabel('Time')
% end
% sgtitle('TCs by Ori - Exc')
% % print(fullfile([data_base '\Analysis\2P\pooled_data_figs\TC_by_ori.pdf']), '-dpdf','-bestfit')
% 
% figure;
% hold on
% for i = 1:length(ori_bins)
%     subplot(n, n2, i)
% %     errorbar(TC_mean_Inh(:,i), TC_sem_Inh(:,i), 'k');
%     shadedErrorBar([],TC_mean_Inh(:,i), TC_sem_Inh(:,i), 'lineProps','r')
%     hold on
%     vline([20 31 42 53])
%     vline(64, 'b:')
%     title(num2str(ori_bins(i)));
%     ylabel('dF/F')
%     xlabel('Time')
% end
% sgtitle('TCs by Ori - Inh')
% % print(fullfile([data_base '\Analysis\2P\pooled_data_figs\TC_by_ori.pdf']), '-dpdf','-bestfit')


%% Modulation by Adapter by Ori; data size = (4,4) = (Ori_bins, Adaptor #)

for i = 1:length(ori_bins)
    AdaptRespforCurrOri_Exc = adapt_resp_Exc{i};
    
    adapt_mean_Exc(i,:) = nanmean(AdaptRespforCurrOri_Exc, 1);
    adapt_sem_Exc(i,:) = std(AdaptRespforCurrOri_Exc,1,1) / sqrt(size(AdaptRespforCurrOri_Exc,1));
    Aix_Exc{i} = AdaptRespforCurrOri_Exc./AdaptRespforCurrOri_Exc(:,1);
    Aix_Exc_mean(i,:) = mean(Aix_Exc{i},1);
    Aix_Exc_sem(i,:) = std(Aix_Exc{i},1,1) / sqrt(size(Aix_Exc{i},1));

    
    AdaptRespforCurrOri_Inh = adapt_resp_Inh{i};
    adapt_mean_Inh(i,:) = nanmean(AdaptRespforCurrOri_Inh, 1);
    adapt_sem_Inh(i,:) = std(AdaptRespforCurrOri_Inh,1,1) ./ sqrt(size(AdaptRespforCurrOri_Inh,1));
    Aix_Inh{i} = AdaptRespforCurrOri_Inh./AdaptRespforCurrOri_Inh(:,1);
    Aix_Inh_mean(i,:) = mean(Aix_Inh{i},1);
    Aix_Inh_sem(i,:) = std(Aix_Inh{i},1,1) / sqrt(size(Aix_Inh{i},1));
    
    
    AdaptRespforCurrOri_all = adapt_resp_all{i};
    adapt_mean_all(i,:) = nanmean(AdaptRespforCurrOri_all, 1);
    adapt_sem_all(i,:) = std(AdaptRespforCurrOri_all,1,1) ./ sqrt(size(AdaptRespforCurrOri_all,1));
    Aix_all{i} = AdaptRespforCurrOri_all./AdaptRespforCurrOri_all(:,1);
    Aix_all_mean(i,:) = mean(Aix_all{i},1);
    Aix_all_sem(i,:) = std(Aix_all{i},1,1) / sqrt(size(Aix_all{i},1));
end

%% Plot Modulation by Adapter by Ori; data size = (4,4) = (Ori_bins, Adaptor #)
% 
% figure;
% hold on
% for i = 1:length(ori_bins)
%     subplot(n, n2, i)
% %     errorbar(adapt_mean_Exc(i,:), adapt_sem_Exc(i,:), 'k');
%     shadedErrorBar([],adapt_mean_Exc(i,:), adapt_sem_Exc(i,:), 'lineProps','b')
%     hold on
%     title(num2str(ori_bins(i)));
%     ylabel('dF/F')
%     xlabel('Adapt #')
% end
% sgtitle('Adapted Response by Ori - Exc')
% 
% figure;
% hold on
% for i = 1:length(ori_bins)
%     subplot(n, n2, i)
% %     errorbar(adapt_mean_Inh(i,:), adapt_sem_Inh(i,:), 'k');
%     shadedErrorBar([],adapt_mean_Inh(i,:), adapt_sem_Inh(i,:), 'lineProps','r')
%     hold on
%     title(num2str(ori_bins(i)));
%     ylabel('dF/F')
%     xlabel('Adapt #')
% end
% sgtitle('Adapted Response by Ori - Inh')
% 
% figure;
% hold on
% for i = 1:length(ori_bins)
%     subplot(n, n2, i)
% %     errorbar(Aix_Exc_mean(i,:), Aix_Exc_sem(i,:), 'k');
%     shadedErrorBar([],Aix_Exc_mean(i,:), Aix_Exc_sem(i,:), 'lineProps','b')
%     hold on
%     title(num2str(ori_bins(i)));
%     ylabel('Aix')
%     xlabel('Adapt #')
% end
% sgtitle('Aix by Ori - Exc')
% 
% figure;
% hold on
% for i = 1:length(ori_bins)
%     subplot(n, n2, i)
% %     errorbar(Aix_Inh_mean(i,:), Aix_Inh_sem(i,:), 'k');
%     shadedErrorBar([],Aix_Inh_mean(i,:), Aix_Inh_sem(i,:), 'lineProps','r')
%     hold on
%     title(num2str(ori_bins(i)));
%     ylabel('Aix')
%     xlabel('Adapt #')
% end
% sgtitle('Aix by Ori - Inh')

% save(fullfile(['Y:\camaron\Analysis\2P\pooled_data_figs\' mouse '\' condition '_post_processing.mat']), 'exp_list', ...
%     'condition', 'vonMisesFitAllCellsAllBoots_all', 'data_adapt_dfof_adapt_all','adapt_resp_cell', 'adapt_cyc_resp_all', 'mask_label_all', ...
%     'mouse', 'adaptPeriodMs', 'dynAdaptPeriodMs', 'dynAdaptFlashOffMs', 'input_temp', 'n_expt_loaded', 'n_cells_total', 'max_resp', 'prefOri', 'prefOri_bootdiff', ...
%     'ind_theta90', 'edges', 'tunedCells', 'tunedCells_bin', 'n_cells_tuned', 'ExcCells', 'InhCells', 'data_adapt_dfof_adapt_all_Exc', 'data_adapt_dfof_adapt_all_Inh',...
%     'data_adapt_dfof_adapt_all_tuned', 'adapt_resp_Exc', 'adapt_resp_Inh', 'adapt_resp_all', 'TC_mean_Exc', 'TC_sem_Exc', 'TC_mean_Inh', 'TC_sem_Inh', 'TC_mean_all', 'TC_sem_all',...
%     'adapt_mean_Exc', 'adapt_sem_Exc', 'Aix_Exc', 'Aix_Exc_mean', 'Aix_Exc_sem', 'adapt_mean_Inh', 'adapt_sem_Inh', 'Aix_Inh', 'Aix_Inh_mean', 'Aix_Inh_sem', ...
%     'adapt_mean_all', 'adapt_sem_all', 'Aix_all', 'Aix_all_mean', 'Aix_all_sem', 'ori_bins', 'adaptor_vline', 'target_vline') 


%NUKE


save(fullfile(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\2P\pooled_data_figs\' mouse '\' condition '_post_processing.mat']), 'exp_list', ...
    'condition', 'vonMisesFitAllCellsAllBoots_all', 'data_adapt_dfof_adapt_all','adapt_resp_cell', 'adapt_cyc_resp_all', 'mask_label_all', ...
    'mouse', 'adaptPeriodMs', 'dynAdaptPeriodMs', 'dynAdaptFlashOffMs', 'input_temp', 'n_expt_loaded', 'n_cells_total', 'max_resp', 'prefOri', 'prefOri_bootdiff', ...
    'ind_theta90', 'edges', 'tunedCells', 'tunedCells_bin', 'n_cells_tuned', 'ExcCells', 'InhCells', 'data_adapt_dfof_adapt_all_Exc', 'data_adapt_dfof_adapt_all_Inh',...
    'data_adapt_dfof_adapt_all_tuned', 'adapt_resp_Exc', 'adapt_resp_Inh', 'adapt_resp_all', 'TC_mean_Exc', 'TC_sem_Exc', 'TC_mean_Inh', 'TC_sem_Inh', 'TC_mean_all', 'TC_sem_all',...
    'adapt_mean_Exc', 'adapt_sem_Exc', 'Aix_Exc', 'Aix_Exc_mean', 'Aix_Exc_sem', 'adapt_mean_Inh', 'adapt_sem_Inh', 'Aix_Inh', 'Aix_Inh_mean', 'Aix_Inh_sem', ...
    'adapt_mean_all', 'adapt_sem_all', 'Aix_all', 'Aix_all_mean', 'Aix_all_sem', 'ori_bins', 'adaptor_vline', 'target_vline', 'data_adapt_dfof_all', 'mask_label_cell', 'adapt_resp_ind_cell', 'vonMisesFitAllCellsAllBoots_cell', 'adapt_trial_ind_cell', 'dates') 


