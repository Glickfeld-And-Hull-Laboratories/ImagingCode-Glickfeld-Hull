clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsAV; %directories
eval(ds)

day_id(2) = 15;
day_id(1) = expt(day_id(2)).multiday_matchdays;
nd = size(day_id,2);

mouse = expt(day_id(1)).mouse;

if strcmp(expt(day_id(1)).data_loc,'lindsey')
    root = fullfile(rc.data,mouse);
elseif strcmp(expt(day_id(1)).data_loc,'ashley')
    root = fullfile(rc.ashleyData,mouse,'two-photon imaging');
end
fnout = fullfile(rc.lindseyAnalysis,mouse);
if expt(day_id(2)).multiday_timesincedrug_hours>0
    dart_str = [expt(day_id(2)).drug '_' num2str(expt(day_id(2)).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end
fn_multi = fullfile(rc.lindseyAnalysis,mouse,['multiday_' dart_str]);

load(fullfile(fn_multi,'timecourses.mat'))
load(fullfile(fn_multi,'multiday_alignment.mat'))
load(fullfile(fn_multi,'input.mat'))



%% stimulus props

tCon = cell(1,nd);
tDir = cell(1,nd);
tOri = cell(1,nd);

for id = 1:nd
    tCon{id} = celleqel2mat_padded(input(id).tGratingContrast);
    tDir{id} = celleqel2mat_padded(input(id).tGratingDirectionDeg);
    tOri{id} = tDir{id};
    tOri{id}(find(tDir{id}>=180)) = tDir{id}(find(tDir{id}>=180))-180;
end
oris = unique(tOri{1});
cons = unique(tCon{1});
nOri = length(oris);
nCon = length(cons);

nOn = input(1).nScansOn;
nOff = input(1).nScansOff;

%% sort TCs
nTrials = length(tCon{1});
data_dfof_match = cell(1,nd);
data_dfof_all = cell(1,nd);
fractTimeActive_all = cell(1,nd);
fractTimeActive_match = cell(1,nd);
for id = 1:nd
    nCells = size(cellTCs_match{id},2);
    nFrames = size(cellTCs_match{id},1);
    fractTimeActive_match{id} = zeros(1,nCells);
    data_trial_match = reshape(cellTCs_match{id},[nOn+nOff nTrials nCells]);
    data_f_match = mean(data_trial_match(1:nOff,:,:),1);
    data_dfof_match{id} = bsxfun(@rdivide,bsxfun(@minus,data_trial_match,data_f_match),data_f_match);
    fractTimeActive_match{id} = zeros(1,nCells);
    meansub_match = cellTCs_match{id}-mean(cellTCs_match{id},1);
    cellstd_match = std(meansub_match,[],1);
    for iCell = 1:nCells
        fractTimeActive_match{id}(:,iCell) = length(find(meansub_match(:,iCell)>3.*cellstd_match(1,iCell)))./nFrames;
    end
    nCells = size(cellTCs_all{id},2);
    data_trial_all = reshape(cellTCs_all{id},[nOn+nOff nTrials nCells]);
    data_f_all = mean(data_trial_all(1:nOff,:,:),1);
    data_dfof_all{id} = bsxfun(@rdivide,bsxfun(@minus,data_trial_all,data_f_all),data_f_all);
    fractTimeActive_all{id} = zeros(1,nCells);
    meansub_all = cellTCs_all{id}-mean(cellTCs_all{id},1);
    cellstd_all = std(meansub_all,[],1);
    for iCell = 1:nCells
        fractTimeActive_all{id}(:,iCell) = length(find(meansub_all(:,iCell)>3.*cellstd_all(1,iCell)))./nFrames;
    end
end
clear data_trial_all data_f_all data_trial_match data_f_match

%%
resp_win = [nOff+5:nOff+nOn];
base_win = [nOff/2:nOff];
resp_cell_all = cell(nd,nOri,nCon);
resp_cell_match = cell(nd,nOri,nCon);
base_cell_all = cell(nd,nOri,nCon);
base_cell_match = cell(nd,nOri,nCon);
h_all = cell(1,nd);
h_match = cell(1,nd);
resp_avg_all = cell(1,nd);
resp_avg_match = cell(1,nd);
resp_ind_all = cell(1,nd);
resp_ind_match = cell(1,nd);

for id = 1:nd
    nCells = size(data_dfof_all{id},3);
    h_all{id} = zeros(nOri,nCon,nCells);
    resp_avg_all{id} = zeros(nOri,nCon,nCells,2);
    for iOri = 1:nOri
        ind_ori = find(tOri{id} == oris(iOri));
        for iCon = 1:nCon
            ind_con = find(tCon{id} == cons(iCon));
            ind = intersect(ind_ori,ind_con);
            resp_cell_all{id,iOri,iCon} = squeeze(mean(data_dfof_all{id}(resp_win,ind,:),1));
            base_cell_all{id,iOri,iCon} = squeeze(mean(data_dfof_all{id}(base_win,ind,:),1));
            for iC = 1:size(resp_cell_all{id,iOri,iCon},2)
                h_all{id}(iOri,iCon,iC) = ttest(resp_cell_all{id,iOri,iCon}(:,iC),base_cell_all{id,iOri,iCon}(:,iC),'tail','right','alpha',0.05./((nOri.*nCon)-1));
            end
            resp_avg_all{id}(iOri,iCon,:,1) = mean(resp_cell_all{id,iOri,iCon},1);
            resp_avg_all{id}(iOri,iCon,:,2) = std(resp_cell_all{id,iOri,iCon},[],1)./sqrt(size(resp_cell_all{id,iOri,iCon},1));
        end
    end
    resp_ind_all{id} = find(sum(sum(h_all{id},1),2));

    nCells = size(data_dfof_match{id},3);
    h_match{id} = zeros(nOri,nCon,nCells);
    p_match{id} = zeros(nOri,nCon,nCells);
    resp_avg_match{id} = zeros(nOri,nCon,nCells,2);
    for iOri = 1:nOri
        ind_ori = find(tOri{id} == oris(iOri));
        for iCon = 1:nCon
            ind_con = find(tCon{id} == cons(iCon));
            ind = intersect(ind_ori,ind_con);
            resp_cell_match{id,iOri,iCon} = squeeze(mean(data_dfof_match{id}(resp_win,ind,:),1));
            base_cell_match{id,iOri,iCon} = squeeze(mean(data_dfof_match{id}(base_win,ind,:),1));
            for iC = 1:size(resp_cell_match{id,iOri,iCon},2)
                [h_match{id}(iOri,iCon,iC), p_match{id}(iOri,iCon,iC)] = ttest(resp_cell_match{id,iOri,iCon}(:,iC),base_cell_match{id,iOri,iCon}(:,iC),'tail','right','alpha',0.05./(nOri-1));
            end
            resp_avg_match{id}(iOri,iCon,:,1) = mean(resp_cell_match{id,iOri,iCon},1);
            resp_avg_match{id}(iOri,iCon,:,2) = std(resp_cell_match{id,iOri,iCon},[],1)./sqrt(size(resp_cell_match{id,iOri,iCon},1));
        end
    end
    resp_ind_match{id} = find(sum(sum(h_match{id},1),2));
end

resp_con = zeros(1,nCon);
for i = 1:nCon
    resp_con(i) = length(find(sum(h_match{2}(:,i,:),1)));
end

%%
resp_max_match = cell(1,nd);    
for id = 1:nd
    resp_max_match{id} = squeeze(max(max(resp_avg_match{id}(:,:,:,1),[],1),[],2));
end
figure; 
for id = 1:nd
    subplot(2,2,id)
    scatter(resp_max_match{1}(resp_ind_match{id}),resp_max_match{2}(resp_ind_match{id}))
    xlabel('D1- max dF/F')
    ylabel('D2- max dF/F')
    xlim([0 1])
    ylim([0 1])
    refline(1)
    title(['Day ' num2str(id) ' Responsive'])
end

print(fullfile(fn_multi,'maxResp_crossDay.pdf'),'-dpdf','-bestfit')
%%
figure;
movegui('center')
od = [1 4];
for id = 1:nd
subplot(2,2,1)
plot(id, length(intersect(find(red_ind_all{id}==0),resp_ind_all{id})),'ok')
hold on
if find(red_ind_all{id})
    plot(id, length(intersect(find(red_ind_all{id}),resp_ind_all{id})),'or')
end
ylim([0 80])
xlim([0 3])
ylabel('Visually responsive (n)')
xlabel('Session')
title('all')
subplot(2,2,2)
plot(id, length(intersect(find(red_ind_match==0),resp_ind_match{id})),'ok')
hold on
if find(red_ind_match==1)
    plot(id, length(intersect(find(red_ind_match==1),resp_ind_match{id})),'or')
end
ylim([0 80])
xlim([0 3])
ylabel('Visually responsive (n)')
xlabel('Session')
title('matched')
subplot(2,2,3)
plot(id, fractTimeActive_all{id}(:,find(red_ind_all{id}==0)),'ok')
hold on
errorbar(id+0.25, mean(fractTimeActive_all{id}(:,find(red_ind_all{id}==0)),2),std(fractTimeActive_all{id}(:,find(red_ind_all{id}==0)),[],2)./sqrt(length(find(red_ind_all{id}==0))),'ok')
if length(find(red_ind_all{id}))>0
    plot(id, fractTimeActive_all{id}(:,find(red_ind_all{id})),'or')
    errorbar(id+0.25, mean(fractTimeActive_all{id}(:,find(red_ind_all{id})),2),std(fractTimeActive_all{id}(:,find(red_ind_all{id})),[],2)./sqrt(length(find(red_ind_all{id}))),'or')
end
ylim([0 0.05])
xlim([0 3])
ylabel('Fraction time active')
xlabel('Session')
end
subplot(2,2,4)
plot(repmat(1:nd, [length(find(red_ind_match==0)) 1])', [fractTimeActive_match{1}(find(red_ind_match==0)); fractTimeActive_match{2}(find(red_ind_match==0))],'-ok')
hold on
if find(red_ind_match==1)
    plot(repmat(1:nd, [length(find(red_ind_match==1)) 1])', [fractTimeActive_match{1}(find(red_ind_match==1)); fractTimeActive_match{2}(find(red_ind_match==1))],'-or')
end
ylim([0 0.05])
xlim([0 3])
ylabel('Fraction time active')
xlabel('Session')
print(fullfile(fn_multi,'responsivity.pdf'),'-dpdf','-bestfit')

R_str = strvcat('R-','R+');
HT_str = strvcat('HT-','HT+');

nCells_match = length(red_ind_match);
figure;
movegui('center')
start = 1;
n = 1;
for iCell = 1:nCells_match
    if start>36
        print(fullfile(fn_multi,['conResp_crossDay' num2str(n) '.pdf']),'-dpdf','-bestfit')
        n = n+1;
        start = 1;
        figure;
    end
    ymax = max(max(resp_avg_match{1}(:,:,iCell,1),[],1),[],2).*1.25;
    if ymax<0.2
        ymax=0.2;
    end
    for iOri = 1:nOri
        subplot(6,6,start)
        errorbar(cons, resp_avg_match{1}(iOri,:,iCell,1), resp_avg_match{1}(iOri,:,iCell,2),'-o')
        if find(resp_ind_match{1} == iCell)
            title(R_str(2))
        end
        hold on
        ylim([-0.1 ymax])
        subplot(6,6,start+1)
        errorbar(cons, resp_avg_match{2}(iOri,:,iCell,1), resp_avg_match{2}(iOri,:,iCell,2),'-o')
        title(HT_str(red_ind_match(iCell)+1,:))
        hold on
        ylim([-0.1 ymax])
    end
    start = start+2;
end
print(fullfile(fn_multi,['conResp_crossDay' num2str(n) '.pdf']),'-dpdf','-bestfit')

%close all
figure;
movegui('center')
start = 1;
n = 1;
match_ind = find([cellImageAlign.pass]);
for iCell = 1:nCells_match
    iC = match_ind(iCell);
    if start>16
        figure;
        movegui('center')
        print(fullfile(fn_multi,['conResp_crossDay_withCell_' num2str(n) '.pdf']),'-dpdf','-bestfit')
        n = n+1;
        start = 1;
    end
    if ~isempty(find(resp_ind_match{1} == iCell)) || ~isempty(find(resp_ind_match{2} == iCell))
        ymax = max(max(resp_avg_match{1}(:,:,iCell,1),[],1),[],2).*1.25;
        if ymax<0.2
            ymax=0.2;
        end
        for iCon = 1:nCon
            subplot(4,4,start)
            errorbar(oris', resp_avg_match{1}(:,iCon,iCell,1), resp_avg_match{1}(:,iCon,iCell,2),'-o')
            hold on
            ylim([-0.1 ymax])
            title(HT_str(red_ind_match(iCell)+1,:))
            subplot(4,4,start+1)
            errorbar(oris', resp_avg_match{2}(:,iCon,iCell,1), resp_avg_match{2}(:,iCon,iCell,2),'-o')
            hold on
            ylim([-0.1 ymax])
        end
        if ~isnan(cellImageAlign(iC).r_red)
            subplot(4,4,start+2)
            imagesc(cellImageAlign(iC).d(1).red_img)
            subplot(4,4,start+3)
            imagesc(cellImageAlign(iC).d(2).red_img)
        else
            subplot(4,4,start+2)
            imagesc(cellImageAlign(iC).d(1).corr_img)
            subplot(4,4,start+3)
            imagesc(cellImageAlign(iC).d(2).corr_img)
        end
        start = start+4;
    end
end
print(fullfile(fn_multi,['conResp_crossDay_withCell_' num2str(n) '.pdf']),'-dpdf','-bestfit')


figure;
movegui('center')
ori_mat = repmat(1:nOri,[nCon 1])';
prefOri_all = cell(1,nd);
resp_avg_all_pref = cell(1,nd);
resp_avg_match_pref = cell(1,nd);
for id = 1:nd
    nC_all(id) = size(resp_avg_all{id},3);
    prefOri_all{id} = zeros(1,nC_all(id));
    resp_avg_all_pref{id} = nan(nCon,nC_all(id));
    for iC = 1:nC_all(id)
        if find(resp_ind_all{id} == iC)
            t = resp_avg_all{id}(:,:,iC,1);
            [max_val max_ind] = max(t(:));
            prefOri_all{id}(1,iC) = ori_mat(max_ind);
            resp_avg_all_pref{id}(:,iC) = squeeze(resp_avg_all{id}(ori_mat(max_ind),:,iC,1));
        end
    end
    subplot(2,2,1)
    errorbar(cons, nanmean(resp_avg_all_pref{id}(:,find(red_ind_all{id}==0)),2), nanstd(resp_avg_all_pref{id}(:,find(red_ind_all{id}==0)),[],2)./sqrt(length(find(red_ind_all{id}==0))),'-o')
    hold on
    title('all HT-')
    ylim([-0.05 0.3])
    subplot(2,2,2)
    if sum(red_ind_all{id})>0
        errorbar(cons, nanmean(resp_avg_all_pref{id}(:,find(red_ind_all{id})),2), nanstd(resp_avg_all_pref{id}(:,find(red_ind_all{id})),[],2)./sqrt(length(find(red_ind_all{id}==0))),'-o')
        hold on
    end
    title('all HT+')
    ylim([-0.05 0.3])
    nC_match = size(resp_avg_match{id},3);
    prefOri_match{id} = zeros(1,nC_match);
    resp_avg_match_pref{id} = nan(nCon,nC_match);
    for iC = 1:nC_match
        if find(resp_ind_match{1} == iC)
            t = resp_avg_match{id}(:,:,iC,1);
            [max_val max_ind] = max(t(:));
            prefOri_match{id}(1,iC) = ori_mat(max_ind);
            resp_avg_match_pref{id}(:,iC) = squeeze(resp_avg_match{id}(ori_mat(max_ind),:,iC,1));
        end
    end
    subplot(2,2,3)
    errorbar(cons, nanmean(resp_avg_match_pref{id}(:,find(red_ind_match==0)),2), nanstd(resp_avg_match_pref{id}(:,find(red_ind_match==0)),[],2)./sqrt(length(find(red_ind_match==0))),'-o')
    hold on
    title('match HT-')
    ylim([-0.05 0.3])
    subplot(2,2,4)
    if sum(red_ind_match)>0
        errorbar(cons, nanmean(resp_avg_match_pref{id}(:,find(red_ind_match)),2), nanstd(resp_avg_match_pref{id}(:,find(red_ind_match)),[],2)./sqrt(length(find(red_ind_match==0))),'-o')
    end
    hold on
    title('match HT+')
    ylim([-0.05 0.3])
end
legend({'D1','D2'})
print(fullfile(fn_multi,'conResp_crossDay_Avg.pdf'),'-dpdf','-bestfit')



figure;
movegui('center')
start = 1;
for iC = find(red_ind_match)
    subplot(6,7,start)
    scatter(resp_avg_match{1}(:,end,iC,1), resp_avg_match{1}(:,end-1,iC,1))
    hold on
    scatter(resp_avg_match{2}(:,end,iC,1), resp_avg_match{2}(:,end-1,iC,1))
    xlim([-.05 1])
    ylim([-.05 1])
    start = start+1;
    title(HT_str(2,:))
end

for iC = find(red_ind_match==0)
    subplot(6,7,start)
    scatter(resp_avg_match{1}(:,end,iC,1), resp_avg_match{1}(:,end-1,iC,1))
    hold on
    scatter(resp_avg_match{2}(:,end,iC,1), resp_avg_match{2}(:,end-1,iC,1))
    xlim([-.05 1])
    ylim([-.05 1])
    start = start+1;
    title(HT_str(1,:))
end
    
