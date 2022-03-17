<<<<<<< HEAD
%%Load data
clear all; clear global; close all
clc

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
fn_analysis = fullfile(fn_base, 'home\Tierney\Analysis\2P');
fn_multi = fullfile(fn_base, 'home\Tierney\Analysis\2P\MultidayAnalysis');

ds = 'ExperimentData_TD'; % dataset info 
eval(ds)

% Baseline (1) and post-MD (2) days
day_id(2) = 11;
day_id(1) = expt(day_id(2)).matchday_baseline;

% % Baseline (1), post-MD (2), and recovery (3) days
% day_id(2) = 2;
% day_id(1) = expt(day_id(2)).matchday_baseline;
% day_id(3) = expt(day_id(2)).matchday_recovery;

nd = length(day_id);

%Specific experiment information
mouse = expt(day_id(1)).mouse;

expt(day_id(2)).multiday_time_days>0
time_str = [num2str(expt(day_id(2)).multiday_time_days) 'Days'];
fn_match = ['multiday_' time_str,'_',expt(day_id(2)).experiment];

stimData = cell(1,nd)

for id = 1:nd
    date = expt(day_id(id)).date;
    runs = expt(day_id(id)).stimruns;
    nrun = length(runs);
    run_str = catRunName(runs, nrun);
    eye_str = expt(day_id(id)).eye_str;
    contra = strcmp(expt(day_id(id)).eye_str,'Contra'); % 1 is contra eye open; 0 is ipsi eye open
    
    datemouse = [date '_' mouse];
    datemouserun = [date '_' mouse '_' run_str];
   
    fn_stimData = fullfile(fn_analysis, datemouse, datemouserun, [datemouserun '_stimData_Ori.mat']);
    temp_stimData = load(fn_stimData);
    stimData{id} = temp_stimData;
    clear temp_stimData
end
 
% Load stimData, timecourses, and input
load(fullfile(fn_multi, mouse, fn_match, ['timecourses.mat']));
load(fullfile(fn_multi, mouse, fn_match, ['input.mat']));

%% Extract variables of interest from the cell arrays loaded above

% Use this section if stimData values are different for each day; would need
% to add (id) in for loops after these variables.

% nOn = zeros(1,nd); 
% nOff = zeros(1,nd); 
% ntrials = zeros(1,nd); 
% Eyes = zeros(2,nd); 
% nEye = zeros(1,nd); 
% Dirs = cell(1,nd); 
% nDirs = zeros(1,nd); 

%  for id = 1:nd
%     nOn(id) = stimData{id}.nOn;
%     nOff(id) = stimData{id}.nOff;
%     ntrials(id) = stimData{id}.ntrials;
%     
%     tContra{id} = stimData{id}.tContra;
%     Eyes(:,id) = stimData{id}.Eyes;
%     nEye(id) = stimData{id}.nEye;
% 
%     tDir{id} = stimData{id}.tDir;
%     Dirs{id} = stimData{id}.Dirs;
%     nDirs(id) = stimData{id}.nDirs;
%  end

tContra = cell(1,nd);
tOri = cell(1,nd);

nOn = stimData{1}.nOn;
nOff = stimData{1}.nOff;

Eyes = stimData{1}.Eyes;
nEye = stimData{1}.nEye;

Oris = stimData{1}.Oris;
nOris = stimData{1}.nOris;
    
for id = 1:nd
    ntrials(id) = stimData{1}.ntrials;
    tContra{id} = stimData{1}.tContra;
    tOri{id} = stimData{1}.tOri;
end    
    
save(fullfile(fn_multi, mouse, fn_match, [mouse '_stimData_multiday.mat']), 'nOn','nOff','Eyes','nEye','Oris','nOris','ntrials','tContra','tOri')

%% looking at time courses: average across all trials

data_dfof_trial = cell(1,nd)

% In the below plots, circshift (-50) was used, so stimOn starts at 10
% frames.
figure
for id = 1:nd
    nCells = size(cellTCs_match{id},2);
    data_tc_trial{id} = reshape(cellTCs_match{id}, [nOn+nOff,ntrials(id),nCells]);
    data_f_trial{id} = mean(data_tc_trial{id}(nOff/2:nOff,:,:),1);
    data_dfof_trial{id} = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial{id}, data_f_trial{id}), data_f_trial{id});

    %looking at data with np subtracted
    tc_cell_avrg{id} = mean(data_dfof_trial{id},3);%average per cells, one row per trial
    tc_trial_avrg{id} = squeeze(mean(data_dfof_trial{id},2));%average over trials, one row per cell
    tc_cell_trial_avrg{id} = mean(tc_cell_avrg{id},2);%average over trials and cells

    subplot(2,1,id)
    plot(circshift(tc_trial_avrg{id},-50), 'LineWidth',.005);
    hold on;
    plot(circshift(tc_cell_trial_avrg{id},-50), 'LineWidth',2, 'color','k');
    hold on;
    title(['Timecourses: ' num2str(expt(id).experiment)]); % Timecourses with neuropil subtracted
    xlabel('Frames')
    ylabel('df/f')
    hold off
end

print(fullfile(fn_multi, mouse, fn_match, [mouse '_Timecourses_multiday.pdf']),'-dpdf','-bestfit')

%% get tuning data
base_win = nOff/2:nOff;
resp_win = nOff+5:nOff+nOn;

for id = 1:nd
    data_trial{id} = reshape(cellTCs_match{id}, [nOn+nOff ntrials(id) nCells]);
    data_f{id} = mean(data_trial{id}(base_win,:,:),1);
    data_dfof{id} = bsxfun(@rdivide,bsxfun(@minus,data_trial{id},data_f{id}),data_f{id});
end

resp_cell_Ori = cell(nEye,nOris,nd);
base_cell_Ori = cell(nEye,nOris,nd);
data_dfof_Ori = cell(1,nd);
h_Ori = cell(1,nd);
p_Ori = cell(1,nd);

for id = 1:nd
    data_dfof_Ori_temp = zeros(nCells,nOris,nEye,2);
    h_Ori_temp = zeros(nOris,nCells,nEye);
    p_Ori_temp = zeros(nOris,nCells,nEye);

    for iEye = 1:nEye
        ind_eye = find(tContra{id} == Eyes(iEye));
        for iOri = 1:nOris
            ind_Ori = find(tOri{id} == Oris(iOri));
            ind = intersect(ind_eye,ind_Ori);
            resp_cell_Ori{iEye,iOri,id} = squeeze(mean(data_dfof{id}(resp_win,ind,:),1));
            base_cell_Ori{iEye,iOri,id} = squeeze(mean(data_dfof{id}(base_win,ind,:),1));
            [h_Ori_temp(iOri,:,iEye), p_Ori_temp(iOri,:,iEye)] = ttest(resp_cell_Ori{iEye,iOri,id},base_cell_Ori{iEye,iOri,id},'tail','right','alpha',0.05./((nOris*nEye)-1));
            data_dfof_Ori_temp(:,iOri,iEye,1) = squeeze(mean(mean(data_dfof{id}(resp_win,ind,:),1),2));
            data_dfof_Ori_temp(:,iOri,iEye,2) = squeeze(std(mean(data_dfof{id}(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
        end
    end
    
    h_Ori{id} = h_Ori_temp;
    p_Ori{id} = p_Ori_temp;
    data_dfof_Ori{id} = data_dfof_Ori_temp;
    
    h_all_Ori{id} = squeeze(sum(h_Ori{id},1));
    resp_ind{id} = find(sum(h_all_Ori{id},2));
    ipsi_resp_ind{id} = find(h_all_Ori{id}(:,1));
    contra_resp_ind{id} = find(h_all_Ori{id}(:,2));
end 

clear h_Ori_temp p_Ori_temp data_dfof_Ori_temp

save(fullfile(fn_multi, mouse, fn_match, [mouse '_respData_multiday.mat']), 'h_Ori', 'resp_ind', 'ipsi_resp_ind', 'contra_resp_ind', 'resp_cell_Ori', 'base_cell_Ori', 'data_dfof_Ori','base_win','resp_win','data_dfof')

figure
for iEye = 1:nEye
    subplot(2,1,iEye)
%     C = linspace(1,10,length(h_all_Ori{2}(:,iEye)));
%     swarmchart((h_all_Ori{1}(:,iEye)),(h_all_Ori{2}(:,iEye)),[],C)
    swarmchart((h_all_Ori{1}(:,iEye)),(h_all_Ori{2}(:,iEye)))
    axis square
    title(['Significant orientations: ' eye_str{find(contra(1:2)==Eyes(iEye))}]);
    xlabel([num2str(expt(1).experiment)])
    ylabel([num2str(expt(2).experiment)])
    refline(1)
    hold on    
end

print(fullfile(fn_multi, mouse, fn_match, [mouse '_sigOrisScatter_multiday.pdf']),'-dpdf','-bestfit')

%% get ODI

for id = 1:nd
    [contra_resp{id} max_ind_contra{id}] = max(data_dfof_Ori{id}(:,:,find(Eyes),1),[],2);
    max_Ori_contra{id} = Oris(max_ind_contra{id});
    [ipsi_resp{id} max_ind_ipsi{id}] = max(data_dfof_Ori{id}(:,:,find(~Eyes),1),[],2); 
    max_Ori_ipsi{id} = Oris(max_ind_ipsi{id});

    real_contra_resp{id} = contra_resp{id};
    real_contra_ind{id} = find(real_contra_resp{id} < 0);
    real_contra_resp{id}(real_contra_ind{id}) = 0;

    real_ipsi_resp{id} = ipsi_resp{id};
    real_ipsi_ind{id} = find(real_ipsi_resp{id} < 0);
    real_ipsi_resp{id}(real_ipsi_ind{id}) = 0;

    ODI{id} = (real_contra_resp{id}-real_ipsi_resp{id})./(real_contra_resp{id}+real_ipsi_resp{id});
end

%scatter of max response to contra and ipsi
figure
for id = 1:nd
    subplot(2,2,id)
    contra_resp_any{id} = max(data_dfof_Ori{id}(resp_ind{id},:,find(Eyes),1),[],2);
    ipsi_resp_any{id} = max(data_dfof_Ori{id}(resp_ind{id},:,find(~Eyes),1),[],2);
%     C = linspace(1,10,length(contra_resp_any{id}(:,1)));
%     scatter(contra_resp_any{id},ipsi_resp_any{id},[],C)
    scatter(contra_resp_any{id},ipsi_resp_any{id})
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    axis square
    xlabel('Max (Contra)')
    ylabel('Max (Ipsi)')
    xlim([0.01 1])
    ylim([0.01 1])
    refline(1)
    title(['Responsive cells: ' num2str(expt(id).experiment)])
end

% ODI
for id = 1:nd
    subplot(2,2,id+2)
    ODI_any{id} = (contra_resp_any{id}-ipsi_resp_any{id})./(contra_resp_any{id}+ipsi_resp_any{id});
    histogram(ODI_any{id},[-1:0.1:1])
    [temp_counts_any{id},~] = histcounts(ODI_any{id},[-1:0.1:1])
    xlabel('ODI')
    ylabel('Cells')
    axis square
    title(['ODI: ' num2str(expt(id).experiment)])

    max_histcounts = max([temp_counts_any{id}]);
    subplot(2,2,id+2)
    axis([-1 1 0 max_histcounts+1])
end

print(fullfile(fn_multi, mouse, fn_match, [mouse '_OD_multiday.pdf']),'-dpdf','-bestfit')

save(fullfile(fn_multi, mouse, fn_match, [mouse '_ODI_multiday.mat']), 'ODI', 'real_contra_resp', 'real_ipsi_resp', 'max_Ori_contra', 'max_Ori_ipsi')

%% Population average for contra/ipsi inputs

%Get responses from individual eyes
mean_ipsi = zeros(1,nd);
mean_contra = zeros(1,nd);
 
for id = 1:nd
     mean_ipsi(id) = mean(real_ipsi_resp{id});
     mean_contra(id) = mean(real_contra_resp{id});
end

figure
mean_resp = [mean_ipsi; mean_contra]
bar(mean_resp')
axis square
xlabel('Timepoint')
set(gca,'XTick',1:3,'XTickLabel',{'Pre','Post','Rec'})
ylabel('Mean response (df/f)') % CHECK THAT THIS IS ACTUALLY DF/F
legend('Ipsi','Contra')
title(['Ipsi v. contra response'])

print(fullfile(fn_multi, mouse, fn_match, [mouse '_EyeRespHist.pdf']),'-dpdf','-bestfit')

save(fullfile(fn_multi, mouse, fn_match, [mouse '_EyeResp.mat']), 'mean_ipsi', 'mean_contra')

%% von mises

    b_ori = cell(1,nd);
    k1_ori = cell(1,nd);
    R1_ori = cell(1,nd);
    u1_ori = cell(1,nd);
    R_square_ori = cell(1,nd);
    sse_ori = cell(1,nd);
    stim_OSI = cell(1,nd);
    theta_hires = deg2rad(0:180);
    y_fit = cell(1,nd);

for id = 1:nd    
    for iEye = 1:nEye
        for iCell = 1:nCells
            data = [data_dfof_Ori{id}(iCell,:,iEye) data_dfof_Ori{id}(iCell,1,iEye)];
            theta = [deg2rad(Oris) pi];
            [b_ori{id}(iEye,iCell),k1_ori{id}(iEye,iCell),R1_ori{id}(iEye,iCell),u1_ori{id}(iEye,iCell),sse_ori{id}(iEye,iCell),R_square_ori{id}(iEye,iCell)] ...
                = miaovonmisesfit_ori(theta,data);
            [max_val{id} max_ind{id}] = max(data_dfof_Ori{id}(iCell,:,iEye),[],2);
            null_ind{id} = max_ind{id}+(nOris./2);
            null_ind{id}(find(null_ind{id}>nOris)) = null_ind{id}(find(null_ind{id}>nOris))-nOris;
            min_val{id} = data_dfof_Ori{id}(iCell,null_ind{id},iEye);
            if min_val{id}<0
                min_val{id} = 0;
            end
            stim_OSI{id}(1,iCell) = (max_val{id}-min_val{id})./(max_val{id}+min_val{id});
            y_fit{id}(:,iEye,iCell) = b_ori{id}(iEye,iCell) + R1_ori{id}(iEye,iCell) .* exp(k1_ori{id}(iEye,iCell).*(cos(2.*(theta_hires-u1_ori{id}(iEye,iCell)))-1));
        end
    end
    
    [yfit_max{id}, yfit_max_ind{id}] = max(y_fit{id}(:,:,:),[],1);
    prefOri_yfit{id} = squeeze(theta_hires(yfit_max_ind{id}));
end   
save(fullfile(fn_multi, mouse, fn_match, [mouse '_oriResp_multiday.mat']), 'data_dfof_Ori','base_win','resp_win','h_Ori', 'b_ori', 'k1_ori', 'R1_ori', 'u1_ori', 'R_square_ori', 'sse_ori','stim_OSI', 'yfit_max', 'yfit_max_ind', 'prefOri_yfit') 

%% Population average for tuning.

% Plot tuning widths in cumulative distribution plot
figure
for iEye = 1:nEye
    subplot(2,2,iEye)
    for id = 1:nd
        cdfplot(k1_ori{id}(iEye,:));
        hold on
    end
    axis square
    xlabel('Tuning width')
    ylabel('Fraction of cells')
    title(['Tuning width: ' eye_str{find(contra(1:2)==Eyes(iEye))}]);
end
legend('Pre','Post','Rec','Location','Southeast');

% Compare tuning widths of individual cells -- plot in scatterplots
% (baseline v. MD)
for iEye = 1:nEye
    subplot(2,2,iEye+2)
    scatter(k1_ori{1}(iEye,:),k1_ori{2}(iEye,:));
    axis square
    title(['Tuning width: ' eye_str{find(contra(1:2)==Eyes(iEye))}]);
    xlabel([num2str(expt(1).experiment)])
    ylabel([num2str(expt(2).experiment)])
    refline(1)
    hold on   
end

print(fullfile(fn_multi, mouse, fn_match, [mouse '_TuningPlots.pdf']),'-dpdf','-bestfit')

%% Compare dfof between pre- and post-MD

% EDIT: FINISH THIS SECTION
figure
for iEye = nEye
    
end

temp_data = data_dfof{1};
temp_max = max(reshape(temp_data,[],size(temp_data,3)),[],1);

% EDIT: PUT TEMP_MAX IN HERE INSTEAD????
[contra_resp{id} max_ind_contra{id}] = max(data_dfof_Ori{id}(:,:,find(Eyes),1),[],2);

subplot(2,1,iEye)
scatter((data_dfof_max{1}(:,iEye)),(data_dfof_max{2}(:,iEye)))
axis square
title(['Max dfof: ' eye_str{find(contra(1:2)==Eyes(iEye))}]);
xlabel([num2str(expt(1).experiment)])
ylabel([num2str(expt(2).experiment)])
refline(1)
=======
%%Load data
clear all; clear global; close all
clc

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
fn_analysis = fullfile(fn_base, 'home\Tierney\Analysis\2P');
fn_multi = fullfile(fn_base, 'home\Tierney\Analysis\2P\MultidayAnalysis');

ds = 'ExperimentData_TD'; % dataset info 
eval(ds)

% Baseline (1) and post-MD (2) days
day_id(2) = 11;
day_id(1) = expt(day_id(2)).matchday_baseline;

% % Baseline (1), post-MD (2), and recovery (3) days
% day_id(2) = 2;
% day_id(1) = expt(day_id(2)).matchday_baseline;
% day_id(3) = expt(day_id(2)).matchday_recovery;

nd = length(day_id);

%Specific experiment information
mouse = expt(day_id(1)).mouse;

expt(day_id(2)).multiday_time_days>0
time_str = [num2str(expt(day_id(2)).multiday_time_days) 'Days'];
fn_match = ['multiday_' time_str,'_',expt(day_id(2)).experiment];

stimData = cell(1,nd)

for id = 1:nd
    date = expt(day_id(id)).date;
    runs = expt(day_id(id)).stimruns;
    nrun = length(runs);
    run_str = catRunName(runs, nrun);
    eye_str = expt(day_id(id)).eye_str;
    contra = strcmp(expt(day_id(id)).eye_str,'Contra'); % 1 is contra eye open; 0 is ipsi eye open
    
    datemouse = [date '_' mouse];
    datemouserun = [date '_' mouse '_' run_str];
   
    fn_stimData = fullfile(fn_analysis, datemouse, datemouserun, [datemouserun '_stimData_Ori.mat']);
    temp_stimData = load(fn_stimData);
    stimData{id} = temp_stimData;
    clear temp_stimData
end
 
% Load stimData, timecourses, and input
load(fullfile(fn_multi, mouse, fn_match, ['timecourses.mat']));
load(fullfile(fn_multi, mouse, fn_match, ['input.mat']));

%% Extract variables of interest from the cell arrays loaded above

% Use this section if stimData values are different for each day; would need
% to add (id) in for loops after these variables.

% nOn = zeros(1,nd); 
% nOff = zeros(1,nd); 
% ntrials = zeros(1,nd); 
% Eyes = zeros(2,nd); 
% nEye = zeros(1,nd); 
% Dirs = cell(1,nd); 
% nDirs = zeros(1,nd); 

%  for id = 1:nd
%     nOn(id) = stimData{id}.nOn;
%     nOff(id) = stimData{id}.nOff;
%     ntrials(id) = stimData{id}.ntrials;
%     
%     tContra{id} = stimData{id}.tContra;
%     Eyes(:,id) = stimData{id}.Eyes;
%     nEye(id) = stimData{id}.nEye;
% 
%     tDir{id} = stimData{id}.tDir;
%     Dirs{id} = stimData{id}.Dirs;
%     nDirs(id) = stimData{id}.nDirs;
%  end

tContra = cell(1,nd);
tOri = cell(1,nd);

nOn = stimData{1}.nOn;
nOff = stimData{1}.nOff;

Eyes = stimData{1}.Eyes;
nEye = stimData{1}.nEye;

Oris = stimData{1}.Oris;
nOris = stimData{1}.nOris;
    
for id = 1:nd
    ntrials(id) = stimData{1}.ntrials;
    tContra{id} = stimData{1}.tContra;
    tOri{id} = stimData{1}.tOri;
end    
    
save(fullfile(fn_multi, mouse, fn_match, [mouse '_stimData_multiday.mat']), 'nOn','nOff','Eyes','nEye','Oris','nOris','ntrials','tContra','tOri')

%% looking at time courses: average across all trials

data_dfof_trial = cell(1,nd)

% In the below plots, circshift (-50) was used, so stimOn starts at 10
% frames.
figure
for id = 1:nd
    nCells = size(cellTCs_match{id},2);
    data_tc_trial{id} = reshape(cellTCs_match{id}, [nOn+nOff,ntrials(id),nCells]);
    data_f_trial{id} = mean(data_tc_trial{id}(nOff/2:nOff,:,:),1);
    data_dfof_trial{id} = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial{id}, data_f_trial{id}), data_f_trial{id});

    %looking at data with np subtracted
    tc_cell_avrg{id} = mean(data_dfof_trial{id},3);%average per cells, one row per trial
    tc_trial_avrg{id} = squeeze(mean(data_dfof_trial{id},2));%average over trials, one row per cell
    tc_cell_trial_avrg{id} = mean(tc_cell_avrg{id},2);%average over trials and cells

    subplot(2,1,id)
    plot(circshift(tc_trial_avrg{id},-50), 'LineWidth',.005);
    hold on;
    plot(circshift(tc_cell_trial_avrg{id},-50), 'LineWidth',2, 'color','k');
    hold on;
    title(['Timecourses: ' num2str(expt(id).experiment)]); % Timecourses with neuropil subtracted
    xlabel('Frames')
    ylabel('df/f')
    hold off
end

print(fullfile(fn_multi, mouse, fn_match, [mouse '_Timecourses_multiday.pdf']),'-dpdf','-bestfit')

%% get tuning data
base_win = nOff/2:nOff;
resp_win = nOff+5:nOff+nOn;

for id = 1:nd
    data_trial{id} = reshape(cellTCs_match{id}, [nOn+nOff ntrials(id) nCells]);
    data_f{id} = mean(data_trial{id}(base_win,:,:),1);
    data_dfof{id} = bsxfun(@rdivide,bsxfun(@minus,data_trial{id},data_f{id}),data_f{id});
end

resp_cell_Ori = cell(nEye,nOris,nd);
base_cell_Ori = cell(nEye,nOris,nd);
data_dfof_Ori = cell(1,nd);
h_Ori = cell(1,nd);
p_Ori = cell(1,nd);

for id = 1:nd
    data_dfof_Ori_temp = zeros(nCells,nOris,nEye,2);
    h_Ori_temp = zeros(nOris,nCells,nEye);
    p_Ori_temp = zeros(nOris,nCells,nEye);

    for iEye = 1:nEye
        ind_eye = find(tContra{id} == Eyes(iEye));
        for iOri = 1:nOris
            ind_Ori = find(tOri{id} == Oris(iOri));
            ind = intersect(ind_eye,ind_Ori);
            resp_cell_Ori{iEye,iOri,id} = squeeze(mean(data_dfof{id}(resp_win,ind,:),1));
            base_cell_Ori{iEye,iOri,id} = squeeze(mean(data_dfof{id}(base_win,ind,:),1));
            [h_Ori_temp(iOri,:,iEye), p_Ori_temp(iOri,:,iEye)] = ttest(resp_cell_Ori{iEye,iOri,id},base_cell_Ori{iEye,iOri,id},'tail','right','alpha',0.05./((nOris*nEye)-1));
            data_dfof_Ori_temp(:,iOri,iEye,1) = squeeze(mean(mean(data_dfof{id}(resp_win,ind,:),1),2));
            data_dfof_Ori_temp(:,iOri,iEye,2) = squeeze(std(mean(data_dfof{id}(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
        end
    end
    
    h_Ori{id} = h_Ori_temp;
    p_Ori{id} = p_Ori_temp;
    data_dfof_Ori{id} = data_dfof_Ori_temp;
    
    h_all_Ori{id} = squeeze(sum(h_Ori{id},1));
    resp_ind{id} = find(sum(h_all_Ori{id},2));
    ipsi_resp_ind{id} = find(h_all_Ori{id}(:,1));
    contra_resp_ind{id} = find(h_all_Ori{id}(:,2));
end 

clear h_Ori_temp p_Ori_temp data_dfof_Ori_temp

save(fullfile(fn_multi, mouse, fn_match, [mouse '_respData_multiday.mat']), 'h_Ori', 'resp_ind', 'ipsi_resp_ind', 'contra_resp_ind', 'resp_cell_Ori', 'base_cell_Ori', 'data_dfof_Ori','base_win','resp_win','data_dfof')

figure
for iEye = 1:nEye
    subplot(2,1,iEye)
%     C = linspace(1,10,length(h_all_Ori{2}(:,iEye)));
%     swarmchart((h_all_Ori{1}(:,iEye)),(h_all_Ori{2}(:,iEye)),[],C)
    swarmchart((h_all_Ori{1}(:,iEye)),(h_all_Ori{2}(:,iEye)))
    axis square
    title(['Significant orientations: ' eye_str{find(contra(1:2)==Eyes(iEye))}]);
    xlabel([num2str(expt(1).experiment)])
    ylabel([num2str(expt(2).experiment)])
    refline(1)
    hold on    
end

print(fullfile(fn_multi, mouse, fn_match, [mouse '_sigOrisScatter_multiday.pdf']),'-dpdf','-bestfit')

%% get ODI

for id = 1:nd
    [contra_resp{id} max_ind_contra{id}] = max(data_dfof_Ori{id}(:,:,find(Eyes),1),[],2);
    max_Ori_contra{id} = Oris(max_ind_contra{id});
    [ipsi_resp{id} max_ind_ipsi{id}] = max(data_dfof_Ori{id}(:,:,find(~Eyes),1),[],2); 
    max_Ori_ipsi{id} = Oris(max_ind_ipsi{id});

    real_contra_resp{id} = contra_resp{id};
    real_contra_ind{id} = find(real_contra_resp{id} < 0);
    real_contra_resp{id}(real_contra_ind{id}) = 0;

    real_ipsi_resp{id} = ipsi_resp{id};
    real_ipsi_ind{id} = find(real_ipsi_resp{id} < 0);
    real_ipsi_resp{id}(real_ipsi_ind{id}) = 0;

    ODI{id} = (real_contra_resp{id}-real_ipsi_resp{id})./(real_contra_resp{id}+real_ipsi_resp{id});
end

%scatter of max response to contra and ipsi
figure
for id = 1:nd
    subplot(2,2,id)
    contra_resp_any{id} = max(data_dfof_Ori{id}(resp_ind{id},:,find(Eyes),1),[],2);
    ipsi_resp_any{id} = max(data_dfof_Ori{id}(resp_ind{id},:,find(~Eyes),1),[],2);
%     C = linspace(1,10,length(contra_resp_any{id}(:,1)));
%     scatter(contra_resp_any{id},ipsi_resp_any{id},[],C)
    scatter(contra_resp_any{id},ipsi_resp_any{id})
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    axis square
    xlabel('Max (Contra)')
    ylabel('Max (Ipsi)')
    xlim([0.01 1])
    ylim([0.01 1])
    refline(1)
    title(['Responsive cells: ' num2str(expt(id).experiment)])
end

% ODI
for id = 1:nd
    subplot(2,2,id+2)
    ODI_any{id} = (contra_resp_any{id}-ipsi_resp_any{id})./(contra_resp_any{id}+ipsi_resp_any{id});
    histogram(ODI_any{id},[-1:0.1:1])
    [temp_counts_any{id},~] = histcounts(ODI_any{id},[-1:0.1:1])
    xlabel('ODI')
    ylabel('Cells')
    axis square
    title(['ODI: ' num2str(expt(id).experiment)])

    max_histcounts = max([temp_counts_any{id}]);
    subplot(2,2,id+2)
    axis([-1 1 0 max_histcounts+1])
end

print(fullfile(fn_multi, mouse, fn_match, [mouse '_OD_multiday.pdf']),'-dpdf','-bestfit')

save(fullfile(fn_multi, mouse, fn_match, [mouse '_ODI_multiday.mat']), 'ODI', 'real_contra_resp', 'real_ipsi_resp', 'max_Ori_contra', 'max_Ori_ipsi')

%% Population average for contra/ipsi inputs

%Get responses from individual eyes
mean_ipsi = zeros(1,nd);
mean_contra = zeros(1,nd);
 
for id = 1:nd
     mean_ipsi(id) = mean(real_ipsi_resp{id});
     mean_contra(id) = mean(real_contra_resp{id});
end

figure
mean_resp = [mean_ipsi; mean_contra]
bar(mean_resp')
axis square
xlabel('Timepoint')
set(gca,'XTick',1:3,'XTickLabel',{'Pre','Post','Rec'})
ylabel('Mean response (df/f)') % CHECK THAT THIS IS ACTUALLY DF/F
legend('Ipsi','Contra')
title(['Ipsi v. contra response'])

print(fullfile(fn_multi, mouse, fn_match, [mouse '_EyeRespHist.pdf']),'-dpdf','-bestfit')

save(fullfile(fn_multi, mouse, fn_match, [mouse '_EyeResp.mat']), 'mean_ipsi', 'mean_contra')

%% von mises

    b_ori = cell(1,nd);
    k1_ori = cell(1,nd);
    R1_ori = cell(1,nd);
    u1_ori = cell(1,nd);
    R_square_ori = cell(1,nd);
    sse_ori = cell(1,nd);
    stim_OSI = cell(1,nd);
    theta_hires = deg2rad(0:180);
    y_fit = cell(1,nd);

for id = 1:nd    
    for iEye = 1:nEye
        for iCell = 1:nCells
            data = [data_dfof_Ori{id}(iCell,:,iEye) data_dfof_Ori{id}(iCell,1,iEye)];
            theta = [deg2rad(Oris) pi];
            [b_ori{id}(iEye,iCell),k1_ori{id}(iEye,iCell),R1_ori{id}(iEye,iCell),u1_ori{id}(iEye,iCell),sse_ori{id}(iEye,iCell),R_square_ori{id}(iEye,iCell)] ...
                = miaovonmisesfit_ori(theta,data);
            [max_val{id} max_ind{id}] = max(data_dfof_Ori{id}(iCell,:,iEye),[],2);
            null_ind{id} = max_ind{id}+(nOris./2);
            null_ind{id}(find(null_ind{id}>nOris)) = null_ind{id}(find(null_ind{id}>nOris))-nOris;
            min_val{id} = data_dfof_Ori{id}(iCell,null_ind{id},iEye);
            if min_val{id}<0
                min_val{id} = 0;
            end
            stim_OSI{id}(1,iCell) = (max_val{id}-min_val{id})./(max_val{id}+min_val{id});
            y_fit{id}(:,iEye,iCell) = b_ori{id}(iEye,iCell) + R1_ori{id}(iEye,iCell) .* exp(k1_ori{id}(iEye,iCell).*(cos(2.*(theta_hires-u1_ori{id}(iEye,iCell)))-1));
        end
    end
    
    [yfit_max{id}, yfit_max_ind{id}] = max(y_fit{id}(:,:,:),[],1);
    prefOri_yfit{id} = squeeze(theta_hires(yfit_max_ind{id}));
end   
save(fullfile(fn_multi, mouse, fn_match, [mouse '_oriResp_multiday.mat']), 'data_dfof_Ori','base_win','resp_win','h_Ori', 'b_ori', 'k1_ori', 'R1_ori', 'u1_ori', 'R_square_ori', 'sse_ori','stim_OSI', 'yfit_max', 'yfit_max_ind', 'prefOri_yfit') 

%% Population average for tuning.

% Plot tuning widths in cumulative distribution plot
figure
for iEye = 1:nEye
    subplot(2,2,iEye)
    for id = 1:nd
        cdfplot(k1_ori{id}(iEye,:));
        hold on
    end
    axis square
    xlabel('Tuning width')
    ylabel('Fraction of cells')
    title(['Tuning width: ' eye_str{find(contra(1:2)==Eyes(iEye))}]);
end
legend('Pre','Post','Rec','Location','Southeast');

% Compare tuning widths of individual cells -- plot in scatterplots
% (baseline v. MD)
for iEye = 1:nEye
    subplot(2,2,iEye+2)
    scatter(k1_ori{1}(iEye,:),k1_ori{2}(iEye,:));
    axis square
    title(['Tuning width: ' eye_str{find(contra(1:2)==Eyes(iEye))}]);
    xlabel([num2str(expt(1).experiment)])
    ylabel([num2str(expt(2).experiment)])
    refline(1)
    hold on   
end

print(fullfile(fn_multi, mouse, fn_match, [mouse '_TuningPlots.pdf']),'-dpdf','-bestfit')

%% Compare dfof between pre- and post-MD

% EDIT: FINISH THIS SECTION
figure
for iEye = nEye
    
end

temp_data = data_dfof{1};
temp_max = max(reshape(temp_data,[],size(temp_data,3)),[],1);

% EDIT: PUT TEMP_MAX IN HERE INSTEAD????
[contra_resp{id} max_ind_contra{id}] = max(data_dfof_Ori{id}(:,:,find(Eyes),1),[],2);

subplot(2,1,iEye)
scatter((data_dfof_max{1}(:,iEye)),(data_dfof_max{2}(:,iEye)))
axis square
title(['Max dfof: ' eye_str{find(contra(1:2)==Eyes(iEye))}]);
xlabel([num2str(expt(1).experiment)])
ylabel([num2str(expt(2).experiment)])
refline(1)
>>>>>>> 56f5ff93 (Merge branch 'master' of github.com:Glickfeld-And-Hull-Laboratories/ImagingCode-Glickfeld-Hull)
hold on   