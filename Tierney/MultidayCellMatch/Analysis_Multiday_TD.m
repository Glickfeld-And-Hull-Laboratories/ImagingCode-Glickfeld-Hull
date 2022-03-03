%%Load data
clear all; clear global; close all
clc

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
fn_analysis = fullfile(fn_base, 'home\Tierney\Analysis\2P');
fn_multi = fullfile(fn_base, 'home\Tierney\Analysis\2P\MultidayAnalysis');
fn_pop = fullfile(fn_base, 'home\Tierney\Analysis\2P\PopulationAnalysis');

ds = 'ExperimentData_TD'; % dataset info 
eval(ds)

% Baseline (1), post-MD (2), and recovery (3) days
day_id(2) = 4;
day_id(1) = expt(day_id(2)).matchday_baseline;
% day_id(3) = expt(day_id(2)).matchday_recovery;

nd = length(day_id);

%Specific experiment information
mouse = expt(day_id(1)).mouse;

expt(day_id(2)).multiday_time_days>0
time_str = [num2str(expt(day_id(2)).multiday_time_days) 'Days'];
fn_match = ['multiday_' time_str,'_',expt(day_id(2)).experiment];

for id = 1:nd
    date = expt(day_id(id)).date;
    runs = expt(day_id(id)).stimruns;
    nrun = length(runs);
    contra = strcmp(expt(day_id(id)).eye_str,'Contra'); % 1 is contra eye open; 0 is ipsi eye open
    
    run_str = catRunName(runs, nrun);
    datemouse = [date '_' mouse];
    datemouserun = [date '_' mouse '_' run_str];
   
    % Load stimData, timecourses, and input
    fn_stimData = fullfile(fn_analysis, datemouse, datemouserun, [datemouserun '_stimData.mat']);
    fn_timecourses = fullfile(fn_multi, mouse, fn_match, ['timecourses.mat']);
    fn_input = fullfile(fn_multi, mouse, fn_match, ['input.mat']);
    temp_stimData = load(fn_stimData);
    temp_timecourses = load(fn_timecourses);
    temp_input = load(fn_input);
      
    stimData{id} = temp_stimData;
    timecourses{id} = temp_timecourses;
    input{id} = temp_input;
    
    clear temp_stimData
    clear temp_timecourses
    clear temp_input
end

%% Extract variables of interest from the cell arrays loaded above
 
 for id = 1:nd
     
%     nOn = input{id}.nScansOn;
%     nOff = input{id}.nScansOff;
%     ntrials = size(input{id}.tGratingDirectionDeg,2);
%      
%     tContra = celleqel2mat_padded(input{id}.contraTrialNumber); %transforms cell array into matrix (1 x ntrials)
%     Eyes = unique(tContra);
%     nEye = length(Eyes);
%     
%     tDir = celleqel2mat_padded(input{id}.tGratingDirectionDeg); %transforms cell array into matrix (1 x ntrials)
%     Dirs = unique(tDir);
%     nDirs = length(Dirs);

    nOn = stimData{id}.nOn;
    nOff = stimData{id}.nOff;
    ntrials = stimData{id}.ntrials;
    
    tContra = stimData{id}.tContra;
    Eyes = stimData{id}.Eyes;
    nEye = stimData{id}.nEye;

    tDir = stimData{id}.tDir;
    Dirs = stimData{id}.Dirs;
    nDirs = stimData{id}.nDirs;
    
    cellTCs = timecourses{id}.cellTCs_match;
    
    % I THINK SECOND ITERATION OF FOR LOOP JUST OVERWRITES FIRST. FIX!!
 end

 nCells = 48 % CHANGE THIS SO NOT HARD CODED!!!!!!

%% looking at time courses

data_tc_trial = reshape(cellTCs, [nOn+nOff,ntrials,nCells]);
data_f_trial = mean(data_tc_trial(nOff/2:nOff,:,:),1);
data_dfof_trial = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial, data_f_trial), data_f_trial);

%looking at data with np subtracted
tc_cell_avrg = mean(data_dfof_trial,3);%average pver cells, one row per trial
tc_trial_avrg = squeeze(mean(data_dfof_trial,2));%average over trials, one row per cell
tc_cell_trial_avrg = mean(tc_cell_avrg,2);%average over trials and cells

figure;
plot(tc_trial_avrg, 'LineWidth',.005,'color',[.25 .25 .25]);
hold on;
plot(tc_cell_trial_avrg, 'LineWidth',2, 'color','k');
hold on;
title('Timecourses with np subtracted');
hold off

%% make tuning curves
base_win = nOff/2:nOff;
resp_win = nOff+5:nOff+nOn;
data_trial = reshape(cellTCs{1}, [nOn+nOff ntrials nCells]);
data_f = mean(data_trial(base_win,:,:),1);
data_dfof = bsxfun(@rdivide,bsxfun(@minus,data_trial,data_f),data_f);

resp_cell_dir = cell(nEye,nDirs);
base_cell_dir = cell(nEye,nDirs);
data_dfof_dir = zeros(nCells,nDirs,nEye,2);
h_dir = zeros(nDirs,nCells,nEye);
p_dir = zeros(nDirs,nCells,nEye);
for iEye = 1:nEye
    ind_eye = find(tContra == Eyes(iEye));
    for iDir = 1:nDirs
        ind_dir = find(tDir==Dirs(iDir));
        ind = intersect(ind_eye,ind_dir);
        resp_cell_dir{iEye,iDir} = squeeze(mean(data_dfof(resp_win,ind,:),1));
        base_cell_dir{iEye,iDir} = squeeze(mean(data_dfof(base_win,ind,:),1));
        [h_dir(iDir,:,iEye), p_dir(iDir,:,iEye)] = ttest(resp_cell_dir{iEye,iDir},base_cell_dir{iEye,iDir},'tail','right','alpha',0.05./((nDirs*nEye)-1));
        data_dfof_dir(:,iDir,iEye,1) = squeeze(mean(mean(data_dfof(resp_win,ind,:),1),2));
        data_dfof_dir(:,iDir,iEye,2) = squeeze(std(mean(data_dfof(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
    end
end

h_all_dir = squeeze(sum(h_dir,1));
resp_ind = find(sum(h_all_dir,2));
ipsi_resp_ind = find(h_all_dir(:,1));
contra_resp_ind = find(h_all_dir(:,2)); 

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_respData.mat']), 'h_dir', 'resp_ind', 'ipsi_resp_ind', 'contra_resp_ind', 'resp_cell_dir', 'base_cell_dir', 'data_dfof_dir','base_win','resp_win','data_dfof')

%plot histogram of number of significant directions for contra and ipsi
figure;
[n n2] = subplotn(nEye);
for iEye = 1:nEye
subplot(n,n2,iEye)
histogram(h_all_dir(:,iEye),[0:1:nDirs])
title(eye_str{find(contra==Eyes(iEye))})
xlabel('Significant directions')
ylabel('Cells')
ylim([0 50]);
end
sgtitle([mouse ' ' date])
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_sigDirsHist.pdf']),'-dpdf','-bestfit')

%plot all cells for all direction to both contra and ipsi
start=1;
n = 1;
figure;
movegui('center')
for iCell = 1:nCells
    if start>25
        sgtitle([mouse ' ' date])
        print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_cellTuningDir' num2str(n) '.pdf']),'-dpdf','-bestfit')
        figure;movegui('center');
        start = 1;
        n = n+1;
    end
    subplot(5,5,start)
    for iEye = 1:2
        errorbar(Dirs, data_dfof_dir(iCell,:,iEye,1), data_dfof_dir(iCell,:,iEye,2), '-o')
        hold on
    end
    title(['R = ' num2str(h_all_dir(iCell,:))])
    start = start +1;
end
sgtitle([mouse ' ' date])
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_cellTuningDir' num2str(n) '.pdf']),'-dpdf','-bestfit')

[contra_resp max_ind_contra] = max(data_dfof_dir(:,:,find(Eyes),1),[],2);
max_dir_contra = Dirs(max_ind_contra);
[ipsi_resp max_ind_ipsi] = max(data_dfof_dir(:,:,find(~Eyes),1),[],2); 
max_dir_ipsi = Dirs(max_ind_ipsi);

real_contra_resp = contra_resp;
real_contra_ind = find(real_contra_resp < 0);
real_contra_resp(real_contra_ind) = 0;

real_ipsi_resp = ipsi_resp;
real_ipsi_ind = find(real_ipsi_resp < 0);
real_ipsi_resp(real_ipsi_ind) = 0;

ODI = (real_contra_resp-real_ipsi_resp)./(real_contra_resp+real_ipsi_resp);

%scatter of max response to contra and ipsi
figure;
subplot(2,3,1)
contra_resp_any = max(data_dfof_dir(resp_ind,:,find(Eyes),1),[],2);
ipsi_resp_any = max(data_dfof_dir(resp_ind,:,find(~Eyes),1),[],2);
scatter(contra_resp_any,ipsi_resp_any)
axis square
refline(1)
xlabel('Max (Contra)')
ylabel('Max (Ipsi)')
xlim([0 1])
ylim([0 1])
title('Responsive (any)')

subplot(2,3,4)
ODI_any = (contra_resp_any-ipsi_resp_any)./(contra_resp_any+ipsi_resp_any);
histogram(ODI_any,[-1:0.1:1])
[temp_counts_any,~] = histcounts(ODI_any,[-1:0.1:1])
xlabel('ODI')
ylabel('Cells')
axis square
title(['Sig cells = ' num2str(length(resp_ind))])

subplot(2,3,2)
contra_resp_contraonly = max(data_dfof_dir(contra_resp_ind,:,find(Eyes),1),[],2); 
ipsi_resp_contraonly = max(data_dfof_dir(contra_resp_ind,:,find(~Eyes),1),[],2);
scatter(contra_resp_contraonly,ipsi_resp_contraonly)
axis square
refline(1)
xlabel('Max (Contra)')
ylabel('Max (Ipsi)')
xlim([0 1])
ylim([0 1])
title('Contra responsive')

subplot(2,3,5)
ODI_contraonly = (contra_resp_contraonly-ipsi_resp_contraonly)./(contra_resp_contraonly+ipsi_resp_contraonly);
histogram(ODI_contraonly,[-1:0.1:1])
[temp_counts_contraonly,~] = histcounts(ODI_contraonly,[-1:0.1:1])
xlabel('ODI')
ylabel('Cells')
axis square
title(['Sig cells = ' num2str(length(contra_resp_ind))])

subplot(2,3,3)
contra_resp_ipsionly = max(data_dfof_dir(ipsi_resp_ind,:,find(Eyes),1),[],2);
ipsi_resp_ipsionly = max(data_dfof_dir(ipsi_resp_ind,:,find(~Eyes),1),[],2);
scatter(contra_resp_ipsionly,ipsi_resp_ipsionly)
axis square
refline(1)
xlabel('Max (Contra)')
ylabel('Max (Ipsi)')
xlim([0 1])
ylim([0 1])
title('Ipsi responsive')

subplot(2,3,6)
ODI_ipsionly = (contra_resp_ipsionly-ipsi_resp_ipsionly)./(contra_resp_ipsionly+ipsi_resp_ipsionly);
histogram(ODI_ipsionly,[-1:0.1:1])
[temp_counts_ipsionly,~] = histcounts(ODI_ipsionly,[-1:0.1:1]);
xlabel('ODI')
ylabel('Cells')
axis square
title(['Sig cells = ' num2str(length(ipsi_resp_ind))])

max_histcounts = max([temp_counts_any temp_counts_contraonly temp_counts_ipsionly]);
subplot(2,3,4)
axis([-1 1 0 max_histcounts+1])
subplot(2,3,5)
axis([-1 1 0 max_histcounts+1])
subplot(2,3,6)
axis([-1 1 0 max_histcounts+1])

print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_OD.pdf']),'-dpdf','-bestfit')

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_ODI.mat']), 'ODI', 'real_contra_resp', 'real_ipsi_resp', 'max_dir_contra', 'max_dir_ipsi')

%% von mises 

    b_ori = zeros(nEye,nCells);
    k1_ori = zeros(nEye,nCells);
    R1_ori = zeros(nEye,nCells);
    u1_ori = zeros(nEye,nCells);
    R_square_ori = zeros(nEye,nCells);
    sse_ori = zeros(nEye,nCells);
    stim_DSI = zeros(nEye,nCells);
    stim_OSI = zeros(nEye,nCells);
    theta_hires= deg2rad(0:180);
    y_fit = zeros(length(theta_hires),nEye,nCells);
    
    Oris = Dirs(1:nDirs/2);
    nOris = length(Oris);
    for iEye = 1:nEye
        data_dfof_ori(:,:,iEye)= mean(reshape(data_dfof_dir(:,:,iEye,1),[nCells nOris 2]),3);
        for iCell = 1:nCells
            data = [data_dfof_ori(iCell,:,iEye) data_dfof_ori(iCell,1,iEye)];
            theta = [deg2rad(Oris) pi];
            [b_ori(iEye,iCell),k1_ori(iEye,iCell),R1_ori(iEye,iCell),u1_ori(iEye,iCell),sse_ori(iEye,iCell),R_square_ori(iEye,iCell)] ...
                = miaovonmisesfit_ori(theta,data);
            [max_val max_ind] = max(data_dfof_ori(iCell,:,iEye),[],2);
            null_ind = max_ind+(nOris./2);
            null_ind(find(null_ind>nOris)) = null_ind(find(null_ind>nOris))-nOris;
            min_val = data_dfof_ori(iCell,null_ind,iEye);
            if min_val<0
                min_val = 0;
            end
            stim_OSI(1,iCell) = (max_val-min_val)./(max_val+min_val);
            [max_val max_ind] = max(data_dfof_dir(iCell,:,iEye,1),[],2);
            null_ind = max_ind+(nDirs./2);
            null_ind(find(null_ind>nDirs)) = null_ind(find(null_ind>nDirs))-nDirs;
            min_val = data_dfof_dir(iCell,null_ind,iEye,1);
            if min_val<0
                min_val = 0;
            end
            stim_DSI(1,iCell) = (max_val-min_val)./(max_val+min_val);
            y_fit(:,iEye,iCell) = b_ori(iEye,iCell) + R1_ori(iEye,iCell) .* exp(k1_ori(iEye,iCell).*(cos(2.*(theta_hires-u1_ori(iEye,iCell)))-1));
        end
    end
    
    [yfit_max, yfit_max_ind] = max(y_fit(:,:,:),[],1);
    prefOri_yfit = squeeze(theta_hires(yfit_max_ind));
    
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_oriResp.mat']), 'data_dfof_dir', 'data_dfof_ori','base_win','resp_win','h_dir', 'b_ori', 'k1_ori', 'R1_ori', 'u1_ori', 'R_square_ori', 'sse_ori','stim_DSI','stim_OSI', 'yfit_max', 'yfit_max_ind', 'prefOri_yfit') 