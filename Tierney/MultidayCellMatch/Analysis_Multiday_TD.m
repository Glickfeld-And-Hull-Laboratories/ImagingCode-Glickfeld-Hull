%%Load data
clear all; clear global; close all
clc

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
fn_analysis = fullfile(fn_base, 'home\Tierney\Analysis\2P');
fn_multi = fullfile(fn_base, 'home\Tierney\Analysis\2P\MultidayAnalysis');

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

stimData = cell(1,nd)

for id = 1:nd
    date = expt(day_id(id)).date;
    runs = expt(day_id(id)).stimruns;
    nrun = length(runs);
    contra = strcmp(expt(day_id(id)).eye_str,'Contra'); % 1 is contra eye open; 0 is ipsi eye open
    
    run_str = catRunName(runs, nrun);
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

% Use this section if stimData values are different
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
    
%% looking at time courses

data_dfof_trial = cell(1,nd)

for id = 1:nd
    nCells = size(cellTCs_match{id},2);
    data_tc_trial{id} = reshape(cellTCs_match{id}, [nOn+nOff,ntrials(id),nCells]);
    data_f_trial{id} = mean(data_tc_trial{id}(nOff/2:nOff,:,:),1);
    data_dfof_trial{id} = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial{id}, data_f_trial{id}), data_f_trial{id});

    %looking at data with np subtracted
    tc_cell_avrg{id} = mean(data_dfof_trial{id},3);%average pver cells, one row per trial
    tc_trial_avrg{id} = squeeze(mean(data_dfof_trial{id},2));%average over trials, one row per cell
    tc_cell_trial_avrg{id} = mean(tc_cell_avrg{id},2);%average over trials and cells

    subplot(2,1,id)
    plot(tc_trial_avrg{id}, 'LineWidth',.005,'color',[.25 .25 .25]);
    hold on;
    plot(tc_cell_trial_avrg{id}, 'LineWidth',2, 'color','k');
    hold on;
    title(['Timecourses with np subtracted day ' num2str(id)]);
    hold off
end
    
%% make tuning curves
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

% Plot histogram of number of significant directions for contra and ipsi
% figure;
% [n n2] = subplotn(nEye);
% for iEye = 1:nEye
%     subplot(n,n2,iEye)
%     histogram(h_all_dir(:,iEye),[0:1:nDirs])
%     title(eye_str{find(contra==Eyes(iEye))})
%     xlabel(['Significant directions ' num2str(id)]);
%     ylabel('Cells')
%     ylim([0 50]);
% end
% sgtitle([mouse ' ' date])
% 
% print(fullfile(fn_multi, mouse, fn_match, [datemouserun '_sigDirsHist_multiday.pdf']),'-dpdf','-bestfit')

% EDIT: MAKE THIS DATA INTO SCATTERPLOT INSTEAD. ONE FOR IPSI AND ONE FOR CONTRA (SIDE
% BY SIDE). PRE ON X AXIS AND POST ON Y AXIS, EACH DOT IS A CELL. SEE SLACK
% FROM CELINE

% EDIT: ADD TITLES, AXES TITLES, ETC.
for iEye = 1:nEye
    figure;
    swarmchart((h_all_Ori{1}(:,iEye)),(h_all_Ori{2}(:,iEye)),'k')
end

print(fullfile(fn_multi, mouse, fn_match, [mouse '_sigOrisScatter_multiday.pdf']),'-dpdf','-bestfit')

% %plot all cells for all direction to both contra and ipsi
% start=1;
% n = 1;
% figure;
% movegui('center')
% for iCell = 1:nCells
%     if start>25
%         sgtitle([mouse ' ' date])
%         print(fullfile(fn_multi, mouse, fn_match, [datemouserun '_cellTuningDir_multiday' num2str(n) '.pdf']),'-dpdf','-bestfit')
%         figure;movegui('center');
%         start = 1;
%         n = n+1;
%     end
%     subplot(5,5,start)
%     for iEye = 1:2
%         errorbar(Dirs, data_dfof_dir(iCell,:,iEye,1), data_dfof_dir(iCell,:,iEye,2), '-o')
%         hold on
%     end
%     title(['R = ' num2str(h_all_dir(iCell,:))])
%     start = start +1;
% end
% sgtitle([mouse ' ' date])
% print(fullfile(fn_multi, mouse, fn_match, [datemouserun '_cellTuningDir_multiday' num2str(n) '.pdf']),'-dpdf','-bestfit')

% EDIT: CHECK THAT THIS RUNS, ADDED FOR LOOP
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
    scatter(contra_resp_any{id},ipsi_resp_any{id})
    axis square
    xlabel('Max (Contra)')
    ylabel('Max (Ipsi)')
    xlim([0 1])
    ylim([0 1])
    refline(1)
    title(['Responsive ' num2str(id)])
    % EDIT: SAVE OUTPUT EACH DAY????? OR DOES SUBPLOT SOLVE THIS??
end

print(fullfile(fn_multi, mouse, fn_match, [datemouserun '_ResponsiveCells_multiday.pdf']),'-dpdf','-bestfit')

for id = 1:nd
    subplot(2,2,id+2) % EDIT: CAN I DO THIS??????
    ODI_any{id} = (contra_resp_any{id}-ipsi_resp_any{id})./(contra_resp_any{id}+ipsi_resp_any{id});
    histogram(ODI_any{id},[-1:0.1:1])
    [temp_counts_any{id},~] = histcounts(ODI_any{id},[-1:0.1:1])
    xlabel('ODI')
    ylabel('Cells')
    axis square
    title(['Sig cells = ' num2str(id)])

    max_histcounts = max([temp_counts_any{id}]); %EDIT: ERROR HERE FIX IT
    subplot(2,3,id+2)
    axis([-1 1 0 max_histcounts{id}+1])

    % EDIT: SAVE OUTPUT EACH DAY????? OR DOES SUBPLOT SOLVE THIS??
end

print(fullfile(fn_multi, mouse, fn_match, [datemouserun '_OD_multiday.pdf']),'-dpdf','-bestfit')

% subplot(2,3,2)
% contra_resp_contraonly = max(data_dfof_dir(contra_resp_ind,:,find(Eyes),1),[],2); 
% ipsi_resp_contraonly = max(data_dfof_dir(contra_resp_ind,:,find(~Eyes),1),[],2);
% scatter(contra_resp_contraonly,ipsi_resp_contraonly)
% axis square
% refline(1)
% xlabel('Max (Contra)')
% ylabel('Max (Ipsi)')
% xlim([0 1])
% ylim([0 1])
% title('Contra responsive')
% 
% subplot(2,3,5)
% ODI_contraonly = (contra_resp_contraonly-ipsi_resp_contraonly)./(contra_resp_contraonly+ipsi_resp_contraonly);
% histogram(ODI_contraonly,[-1:0.1:1])
% [temp_counts_contraonly,~] = histcounts(ODI_contraonly,[-1:0.1:1])
% xlabel('ODI')
% ylabel('Cells')
% axis square
% title(['Sig cells = ' num2str(length(contra_resp_ind))])
% 
% subplot(2,3,3)
% contra_resp_ipsionly = max(data_dfof_dir(ipsi_resp_ind,:,find(Eyes),1),[],2);
% ipsi_resp_ipsionly = max(data_dfof_dir(ipsi_resp_ind,:,find(~Eyes),1),[],2);
% scatter(contra_resp_ipsionly,ipsi_resp_ipsionly)
% axis square
% refline(1)
% xlabel('Max (Contra)')
% ylabel('Max (Ipsi)')
% xlim([0 1])
% ylim([0 1])
% title('Ipsi responsive')
% 
% subplot(2,3,6)
% ODI_ipsionly = (contra_resp_ipsionly-ipsi_resp_ipsionly)./(contra_resp_ipsionly+ipsi_resp_ipsionly);
% histogram(ODI_ipsionly,[-1:0.1:1])
% [temp_counts_ipsionly,~] = histcounts(ODI_ipsionly,[-1:0.1:1]);
% xlabel('ODI')
% ylabel('Cells')
% axis square
% title(['Sig cells = ' num2str(length(ipsi_resp_ind))])

% max_histcounts = max([temp_counts_any temp_counts_contraonly temp_counts_ipsionly]);
% subplot(2,3,4)
% axis([-1 1 0 max_histcounts+1])
% subplot(2,3,5)
% axis([-1 1 0 max_histcounts+1])
% subplot(2,3,6)
% axis([-1 1 0 max_histcounts+1])

save(fullfile(fn_multi, mouse, fn_match, [datemouserun '_ODI_multiday.mat']), 'ODI', 'real_contra_resp', 'real_ipsi_resp', 'max_Ori_contra', 'max_Ori_ipsi')

%% von mises 

% EDIT: HOW SHOULD I PREALLOCATE THESE? CELL ARRAYS????

    b_ori = zeros(nEye,nCells);
    k1_ori = zeros(nEye,nCells);
    R1_ori = zeros(nEye,nCells);
    u1_ori = zeros(nEye,nCells);
    R_square_ori = zeros(nEye,nCells);
    sse_ori = zeros(nEye,nCells);
    stim_OSI = zeros(nEye,nCells);
    theta_hires= deg2rad(0:180);
    y_fit = zeros(length(theta_hires),nEye,nCells);

for id = 1:nd    
    for iEye = 1:nEye
        for iCell = 1:nCells
            data = [data_dfof_ori{id}(iCell,:,iEye) data_dfof_ori{id}(iCell,1,iEye)];
            theta = [deg2rad(Oris) pi];
            [b_ori{id}(iEye,iCell),k1_ori{id}(iEye,iCell),R1_ori{id}(iEye,iCell),u1_ori{id}(iEye,iCell),sse_ori{id}(iEye,iCell),R_square_ori{id}(iEye,iCell)] ...
                = miaovonmisesfit_ori(theta,data); % EDIT: DO I NEED {ID} HERE??????
            [max_val{id} max_ind{id}] = max(data_dfof_ori{id}(iCell,:,iEye),[],2);
            null_ind{id} = max_ind{id}+(nOris./2);
            null_ind{id}(find(null_ind{id}>nOris)) = null_ind{id}(find(null_ind{id}>nOris))-nOris;
            min_val{id} = data_dfof_ori{id}(iCell,null_ind{id},iEye);
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
save(fullfile(fn_multi, mouse, fn_match, [datemouserun '_oriResp_multiday.mat']), 'data_dfof_Ori', 'data_dfof_ori','base_win','resp_win','h_Ori', 'b_ori', 'k1_ori', 'R1_ori', 'u1_ori', 'R_square_ori', 'sse_ori','stim_DSI','stim_OSI', 'yfit_max', 'yfit_max_ind', 'prefOri_yfit') 