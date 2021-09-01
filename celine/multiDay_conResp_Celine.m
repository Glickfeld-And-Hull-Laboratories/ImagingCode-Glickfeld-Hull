clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsAV; %directories
eval(ds)

day_id(2) = 60;
day_id(1) = expt(day_id(2)).multiday_matchdays;
nd = size(day_id,2);

mouse = expt(day_id(1)).mouse;

fnout = fullfile(rc.celineAnalysis,mouse);
if expt(day_id(2)).multiday_timesincedrug_hours>0
    dart_str = [expt(day_id(2)).drug '_' num2str(expt(day_id(2)).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end


fn_multi = fullfile(rc.celineAnalysis,mouse,['multiday_' dart_str]);


load(fullfile(fn_multi,'timecourses.mat'))
load(fullfile(fn_multi,'multiday_alignment.mat'))
load(fullfile(fn_multi,'input.mat'))



%% stimulus props

%tells the contrast, direction and orientation for each trial each day
tCon = cell(1,nd);
tDir = cell(1,nd);
tOri = cell(1,nd);

%find the contrasts, directions and orientations for each day
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
%make mats for all cells and for matched cells showing df/f, etc in terms
%of trials

stimStart = (nOff/2)+1; %this indicates both the perdiod to trim off the start and the stim on period after trimming
stimEnd=stimStart+nOn-1;

cellTCs_match_OG = cellTCs_match;
cellTCs_all_OG = cellTCs_all;

data_dfof_match = cell(1,nd); %make an empty array that is 1 by however many days there are (1X2 usually)
data_dfof_all = cell(1,nd);
fractTimeActive_all = cell(1,nd);
fractTimeActive_match = cell(1,nd);
for id = 1:nd %cycle through days
    nTrials = length(tDir{id})-1; %use the list of direction by trial to figure out how many trials there are
   %currently the way I center the stim on period requires me to cut out
   %one trial, hence the -1
   
    %for matched cell
    nCells = size(cellTCs_match{id},2);
    fractTimeActive_match{id} = zeros(1,nCells);
    cellTCs_match{id} = cellTCs_match{id}(stimStart:size(cellTCs_match{id},1)-stimEnd,:);
    nFrames = size(cellTCs_match{id},1);
    data_trial_match = reshape(cellTCs_match{id},[nOn+nOff nTrials nCells]);
    data_f_match = mean(data_trial_match(1: (nOff/2),:,:),1);
    data_dfof_match{id} = bsxfun(@rdivide,bsxfun(@minus,data_trial_match,data_f_match),data_f_match);
    fractTimeActive_match{id} = zeros(1,nCells);
    meansub_match = cellTCs_match{id}-mean(cellTCs_match{id},1);
    cellstd_match = std(meansub_match,[],1);
    for iCell = 1:nCells
        fractTimeActive_match{id}(:,iCell) = length(find(meansub_match(:,iCell)>3.*cellstd_match(1,iCell)))./nFrames;
    end
    
    %for all cells
    nCells = size(cellTCs_all{id},2);
    cellTCs_all{id} = cellTCs_all{id}(stimStart:size(cellTCs_all{id},1)-stimEnd,:);
    data_trial_all = reshape(cellTCs_all{id},[nOn+nOff nTrials nCells]);
    data_f_all = mean(data_trial_all(1: (nOff/2),:,:),1);
    data_dfof_all{id} = bsxfun(@rdivide,bsxfun(@minus,data_trial_all,data_f_all),data_f_all);
    fractTimeActive_all{id} = zeros(1,nCells);
    meansub_all = cellTCs_all{id}-mean(cellTCs_all{id},1);
    cellstd_all = std(meansub_all,[],1);
    for iCell = 1:nCells
        fractTimeActive_all{id}(:,iCell) = length(find(meansub_all(:,iCell)>3.*cellstd_all(1,iCell)))./nFrames;
    end
end
clear data_trial_all data_f_all data_trial_match data_f_match

%% may need to change this for centered timecourse plot




resp_win = stimStart:stimEnd;
base_win = 1: (nOff/2);

resp_cell_all = cell(nd,nOri,nCon);
resp_cell_match = cell(nd,nOri,nCon);
base_cell_all = cell(nd,nOri,nCon);
base_cell_match = cell(nd,nOri,nCon);
h_all = cell(1,nd);
h_match = cell(1,nd);
resp_avg_all = cell(1,nd); %will be the mean of the response window for all cells
resp_avg_match = cell(1,nd); %will be the mean of the response window for matched cells
resp_ind_all = cell(1,nd);
resp_ind_match = cell(1,nd);

for id = 1:nd
    tCon{id}=tCon{id}(1:nTrials); %here I need to adjust the trial info to drop the last trial
    tOri{id}=tOri{id}(1:nTrials);
    tDir{id}=tDir{id}(1:nTrials);
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
                [h_match{id}(iOri,iCon,iC), p_match{id}(iOri,iCon,iC)] = ttest(resp_cell_match{id,iOri,iCon}(:,iC),base_cell_match{id,iOri,iCon}(:,iC),'tail','right','alpha',0.05./((nOri.*nCon)-1));
            end
            resp_avg_match{id}(iOri,iCon,:,1) = mean(resp_cell_match{id,iOri,iCon},1);
            resp_avg_match{id}(iOri,iCon,:,2) = std(resp_cell_match{id,iOri,iCon},[],1)./sqrt(size(resp_cell_match{id,iOri,iCon},1));
        end
    end
    resp_ind_match{id} = find(sum(sum(h_match{id},1),2));
end

resp_ind_either = unique([resp_ind_match{1};resp_ind_match{2}]); %find the set of cells that was responsive on either day 1 or day 2

resp_con = zeros(1,nCon);
for i = 1:nCon
    resp_con(i) = length(find(sum(h_match{1}(:,i,:),1)));
end
%% get basic counts

%creating the arrays for red cells
red_d1 = find(red_ind_all{1}==1);
red_resp_d1 = intersect(find(red_ind_all{1}==1),resp_ind_all{1});
red_d2 = find(red_ind_all{2}==1);
red_resp_d2 = intersect(find(red_ind_all{2}==1),resp_ind_all{2});
red_match_respd1 = intersect(find(red_ind_match==1),resp_ind_match{1});
red_match_respd2 = intersect(find(red_ind_match==1),resp_ind_match{2});

%matched
nCells_match = length(red_ind_match); %how many matched cells total
nCells_match_respd1 = size(resp_ind_match{1},1); %how many of those responded on d1
nCells_match_respd2 = size(resp_ind_match{2},1);%how many of those responded on d2
nRed_match = sum(red_ind_match);%how many red cells matched
nRed_match_respd1 = length(red_match_respd1); %how many of the matched red cells responded on d1
nRed_match_respd2 = length(red_match_respd2); %how many of the matched red cells responded on d1

%day 1
nCells_d1 = size(data_dfof_all{1},3); %how many cells
nCells_resp_d1 = size(resp_ind_all{1},1); %how many cells responded
nRed_d1 = length(red_d1); %how many red cells
nRed_resp_d1 = length(red_resp_d1); %how many red cells responded

%day 2
nCells_d2 = size(data_dfof_all{2},3); %how many cells
nCells_resp_d2 = size(resp_ind_all{2},1); %how many cells responded
nRed_d2 = length(red_d2); %how many red cells
nRed_resp_d2 = length(red_resp_d2); %how many red cells responded

% make table of values
countsTable = table([nCells_d1; nCells_resp_d1; nRed_d1;nRed_resp_d1],[nCells_d2; nCells_resp_d2; nRed_d2;nRed_resp_d2],[nCells_match;nCells_match_respd1;nRed_match;nRed_match_respd1],[nCells_match;nCells_match_respd2;nRed_match;nRed_match_respd2],'VariableNames',{'Day 1' 'Day 2' 'Matched, day 1' 'Matched, day 2'}, 'RowNames',{'Total cells' 'Responsive cells' 'Red cells' 'Responsive red cells'})
writetable(countsTable,fullfile(fn_multi,'counts.csv'),'WriteRowNames',true)
%% looking at time courses

%day 1
tc_cell_avrg1 = mean(data_tc_trial{1},3);%average pver cells, one row per trial
tc_trial_avrg1 = squeeze(mean(data_tc_trial{1},2));%average over trials, one row per cell
%tc_trial_std = squeeze(std(data_tc_trial{1},2));%std over trials for each cell
%std over cells, after averaging over trials
tc_cell_trial_avrg1 = mean(tc_cell_avrg1,2);%average over trials and cells

figure;
plot(tc_trial_avrg1, 'LineWidth',.005,'color',[.25 .25 .25]);
hold on;
plot(tc_cell_trial_avrg1, 'LineWidth',2, 'color','k');
hold on;
vline(60,'g')
title('Day 1 timecourses');
hold off
print(fullfile(fn_multi,['timecourses_Day1']),'-dpdf');

% looking at day 2
tc_cell_avrg2 = mean(data_tc_trial{2},3);%average pver cells, one row per trial
tc_trial_avrg2 = squeeze(mean(data_tc_trial{2},2));%average over trials, one row per cell
tc_cell_trial_avrg2 = mean(tc_cell_avrg2,2);%average over trials and cells

figure;
plot(tc_trial_avrg2, 'LineWidth',.005,'color',[.25 .25 .25]);
hold on;
plot(tc_cell_trial_avrg2, 'LineWidth',2, 'color','k');
hold on;
vline(60,'g')
title('Day 2 timecourses');
hold off
print(fullfile(fn_multi,['timecourses_Day2']),'-dpdf');

%% looking at time courses matched cells only
%edit this to look at timecourses for responsive green cells, all red
%cells, and look at the 
nTrials = length(tCon{1});
data_tc_trial_match{1} = reshape(data_dfof_match{1}, [nOn+nOff,nTrials,nCells_match]);
nTrials = length(tCon{2});
data_tc_trial_match{2} = reshape(data_dfof_match{2}, [nOn+nOff,nTrials,nCells_match]);
%looking at day 1
tc_cell_match_avrg1 = mean(data_tc_trial_match{1},3);%average pver cells, one row per trial
tc_trial_match_avrg1 = squeeze(mean(data_tc_trial_match{1},2));%average over trials, one row per cell
%tc_trial_std = squeeze(std(data_tc_trial{1},2));%std over trials for each cell
%std over cells, after averaging over trials
tc_cell_trial_match_avrg1 = mean(tc_cell_match_avrg1,2);%average over trials and cells

figure;
subplot(1,2,1)
plot(tc_trial_match_avrg1(30:68,:), 'LineWidth',.005,'color',[.25 .25 .25]);
hold on;
plot(tc_cell_trial_match_avrg1(30:68,:), 'LineWidth',2, 'color','k');
hold on;
vline(30,'g')
title('Day 1 timecourses');
hold off


% looking at day 2
tc_cell_match_avrg2 = mean(data_tc_trial_match{2},3);%average pver cells, one row per trial
tc_trial_match_avrg2 = squeeze(mean(data_tc_trial_match{2},2));%average over trials, one row per cell
%tc_trial_std = squeeze(std(data_tc_trial{1},2));%std over trials for each cell
%std over cells, after averaging over trials
tc_cell_trial_match_avrg2 = mean(tc_cell_match_avrg2,2);%average over trials and cells

subplot(1,2,2)
plot(tc_trial_match_avrg2(30:68,:), 'LineWidth',.005,'color',[.25 .25 .25]);
hold on;
plot(tc_cell_trial_match_avrg2(30:68,:), 'LineWidth',2, 'color','k');
hold on;
vline(30,'g')
title('Day 2 timecourses');
hold off
print(fullfile(fn_multi,['matched timecourses']),'-dpdf');

%% plotting TC for green and red cells seperately
green_inds = 1:nCells;
green_inds = setdiff(green_inds, find(red_ind_match));


figure;
subplot(2,2,1)
plot(tc_trial_match_avrg1(30:68,green_inds), 'LineWidth',.005,'color',[.25 .25 .25]);
ylim([-.03 .05]);
vline(30,'g')
xticks(0:5:38);
xticklabels(-30:5:8);
title('Day 1 timecourses, PYR');
subplot(2,2,2)
plot(tc_trial_match_avrg1(30:68,red_ind_match), 'LineWidth',.005,'color',[.25 .25 .25]);
ylim([-.03 .05]);
vline(30,'g')
xticks(0:5:38);
xticklabels(-30:5:8);
title('Day 1 timecourses, PV');
subplot(2,2,3)
plot(tc_trial_match_avrg2(30:68,green_inds), 'LineWidth',.005,'color',[.25 .25 .25]);
ylim([-.03 .05]);
vline(30,'g')
xticks(0:5:38);
xticklabels(-30:5:8);
title('Day 2 timecourses, PYR');
subplot(2,2,4)
plot(tc_trial_match_avrg2(30:68,red_ind_match), 'LineWidth',.005,'color',[.25 .25 .25]);
ylim([-.03 .05]);
vline(30,'g')
xticks(0:5:38);
xticklabels(-30:5:8);
title('Day 2 timecourses, PV');

%% makes a scatterplot of max df/f for day 1 vs day 2, and each subplot is 
%the set of cells that were significantly responsive to at least once con/ori combo in that data
resp_max_match = cell(1,nd);    
for id = 1:nd
    resp_max_match{id} = squeeze(max(max(resp_avg_match{id}(:,:,:,1),[],1),[],2));
end
figure; movegui('center') 
for id = 1:nd
    subplot(2,2,id)
    scatter(resp_max_match{1}(resp_ind_match{id}),resp_max_match{2}(resp_ind_match{id}))
    xlabel('D1- max dF/F')
    ylabel('D2- max dF/F')
    xlim([0 .2])
    ylim([0 .2])
    refline(1)
    title(['Day ' num2str(id) ' Responsive'])
end

print(fullfile(fn_multi,'maxResp_crossDay.pdf'),'-dpdf','-bestfit')
%% comparative plots of responsivness
%this includes frac active - spontaneous activity?
figure;
movegui('center')
for id = 1:nd
subplot(2,2,1)
plot(id, length(intersect(find(red_ind_all{id}==0),resp_ind_all{id})),'ok')
hold on
if find(red_ind_all{id})
    plot(id, length(intersect(find(red_ind_all{id}),resp_ind_all{id})),'or')
end
ylim([0 100])
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
ylim([0 100])
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
%print(fullfile(fn_multi,'responsivity.pdf'),'-dpdf','-bestfit')

R_str = strvcat('R-','R+');
HT_str = strvcat('HT-','HT+');


%% tuning curves and cell images for matched cells

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
        %print(fullfile(fn_multi,['conResp_crossDay_withCell_' num2str(n) '.pdf']),'-dpdf','-bestfit')
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
            imagesc(cellImageAlign(iC).d(1).avg_img)
            subplot(4,4,start+3)
            imagesc(cellImageAlign(iC).d(2).avg_img)
        end
        start = start+4;
    end
    %print(fullfile(fn_multi,['conResp_crossDay_withCell_' num2str(n) '.pdf']),'-dpdf','-bestfit')
end

%% contrast response curves

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
    xticks(cons)
    ylim([-0.05 0.3])
    subplot(2,2,2)
    xticks(cons)
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
        if find(resp_ind_either == iC)
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
    xticks(cons)
    subplot(2,2,4)
    if sum(red_ind_match)>0
        errorbar(cons, nanmean(resp_avg_match_pref{id}(:,find(red_ind_match)),2), nanstd(resp_avg_match_pref{id}(:,find(red_ind_match)),[],2)./sqrt(length(find(red_ind_match==0))),'-o')
    end
    hold on
    title('match HT+')
    xticks(cons)
    ylim([-0.05 0.3])
end
legend({'D1','D2'})
print(fullfile(fn_multi,'conResp_crossDay_Avg.pdf'),'-dpdf','-bestfit')

%% 


%start by plotting the responses of responzive red cells on each day
%plot ori tuning at each contrast

figure;
movegui('center')
start = 1;
for iCell = 1:length(red_d1)
    iC = red_d1(iCell);
    cell_lab = iCell - 9;
    if start == length(red_d1)
        print(fullfile(fn_multi,['tuning_red_day_1_cell_' num2str(cell_lab) 'to' num2str(iCell)]),'-dpdf','-fillpage');
    elseif start>10
        print(fullfile(fn_multi,['tuning_red_day_1_cell_' num2str(cell_lab) 'to' num2str(iCell)]),'-dpdf','-fillpage');
        figure;
        movegui('center')
        n = n+1;
        start = 1;
    end
   
        ymax = max(max(resp_avg_all{1}(:,:,iCell,1),[],1),[],2).*1.25;
        if ymax<0.2
            ymax=0.2;
        end
        for iCon = 1:nCon
            subplot(5,2,start)
            errorbar(oris', resp_avg_all{1}(:,iCon,iCell,1), resp_avg_all{1}(:,iCon,iCell,2),'-o')
            hold on
            ylim([-0.1 ymax])
            if find(resp_ind_all{1} == iC)
                title('Responsive')
            end
        end

    start = start+1;
end
%%
figure;
movegui('center')
start = 1;
for iCell = 1:length(red_d2)
    iC = red_d2(iCell);
    cell_lab = iCell - 9;
    if start>10
        print(fullfile(fn_multi,['tuning_red_day_2_cell_' num2str(cell_lab) 'to' num2str(iCell)]),'-dpdf','-fillpage');
        figure;
        movegui('center')
        n = n+1;
        start = 1;
    end
   
        ymax = max(max(resp_avg_all{1}(:,:,iCell,1),[],1),[],2).*1.25;
        if ymax<0.2
            ymax=0.2;
        end
        for iCon = 1:nCon
            subplot(5,2,start)
            errorbar(oris', resp_avg_all{2}(:,iCon,iCell,1), resp_avg_all{2}(:,iCon,iCell,2),'-o')
            hold on
            ylim([-0.1 ymax])
            if find(resp_ind_all{2} == iC)
                title('Responsive')
            end
        end
    if start == length(red_d1)
        print(fullfile(fn_multi,['tuning_red_day_2_cell_' num2str(cell_lab) 'to' num2str(iCell)]),'-dpdf','-fillpage');
    end
    start = start+1;
end

%% sort TCs into day 1 and day 2, only good cells
%sort into day 1 and day 2, only want red cells
red_TCs = {cellTCs_all{1}(:,red_d1), cellTCs_all{2}(:,red_d2)};
% spontaneous activity for red cells
% I might adjust this so that it depends on how many red cells there are

pct_events_red = cell(size(red_TCs)); %I don't want pct_events to be overwritten, 
%may need to change this

for iday = 1:size(red_TCs,2) %will loop through each day
    start=1;
    figure
    sgtitle(sprintf(['Red only day ' num2str(iday)]));
    this_day_id = day_id(iday);
    d_all = red_TCs{iday}; %find the data for that day
 
    stimstart = (nOn+1):(nOn+nOff):size(d_all,1)';
    stimon = cell2mat(arrayfun(@(x) x:(x+nOn),stimstart,'unif',0));
    stimoff = setdiff(1:size(d_all,1),stimon);
    d_off = d_all(stimoff,:);
    dff = (d_off-mean(d_off,1))./mean(d_off,1);
    tt = (1:size(d_off,1))./expt(this_day_id).frame_rate./60;
    for icell = 1:size(d_all,2)
        if start>5
            cell_lab = icell - 5;
            print(fullfile(fn_multi,['red_spontaneous_dff_day'  num2str(iday) '_cell_' num2str(cell_lab) 'to' num2str(icell)]),'-dpdf','-fillpage');
            figure;
            sgtitle(sprintf(['Red only day ' num2str(iday)]));
            movegui('center')
            start = 1;
        end
        tc = dff(:,icell);
        subplot(5,1,start)
        plot(tt,tc,'k-','LineWidth',1);
        figXAxis([],'time in expt (min)',[tt(1) tt(end)],0:5:tt(end),0:5:tt(end))
        figYAxis([],'dF/F',[-1 2])
        figAxForm([],0)
        title(sprintf('Cell #%s',num2str(icell)))
        start = start+1;
    end
    
    dff_3sd = (std(dff) + mean(dff))*3;
    dff_test = dff > dff_3sd;
    pct_events_red{iday} = sum(dff_test,1)./size(dff,1);
end

%%
nc = cellfun(@(x) size(x,2),pct_events_red);
figure; 
for iday = 1:size(red_TCs,2);
    subplot(1,2,iday); hold on
    d = pct_events_red{iday};
    histogram(d,10);
    vline(mean(d));
    xlabel('pct events');
    ylabel('# cells');
    title(['Spontaneous events day ',num2str(iday)])
end
print(fullfile(fn_multi,'RedSpontaneousEvents'),'-dpdf')

%% spontaneous activity for all cells

pct_events_all = cell(size(cellTCs_all)); 
nexamplecells = 5;

for iday = 1:size(cellTCs_all,2) %will loop through each day
    this_day_id = day_id(iday);
    d_all = cellTCs_all{iday}; %find the data for that day
    
    if size(d_off,2) < nexamplecells
        n = size(d_off,2);
    else
        n = nexamplecells;
    end
    ind = randsample(size(d_all,2),n);

    figure
    sgtitle(sprintf(['Rand subset all cells day ' num2str(iday)]));
    stimstart = (nOn+1):(nOn+nOff):size(d_all,1)';
    stimon = cell2mat(arrayfun(@(x) x:(x+nOn),stimstart,'unif',0));
    stimoff = setdiff(1:size(d_all,1),stimon);
    d_off = d_all(stimoff,:);
    dff = (d_off-mean(d_off,1))./mean(d_off,1);
    tt = (1:size(d_off,1))./expt(this_day_id).frame_rate./60;
    for icell = 1:n
        tc = dff(:,ind(icell));
        subplot(nexamplecells,1,icell)
        plot(tt,tc,'k-','LineWidth',1);
        figXAxis([],'time in expt (min)',[tt(1) tt(end)],0:5:tt(end),0:5:tt(end))
        figYAxis([],'dF/F',[-1 2])
        figAxForm([],0)
        title(sprintf('Cell #%s',num2str(ind(icell))))
    end
   
    dff_3sd = (std(dff) + mean(dff))*3;
    dff_test = dff > dff_3sd;
    pct_events_all{iday} = sum(dff_test,1)./size(dff,1);
    print(fullfile(fn_multi,['all_spontaneous_dff_day'  num2str(iday) '_random_subset']),'-dpdf','-fillpage');
end
%%
nc = cellfun(@(x) size(x,2),pct_events_all);
figure; 
for iday = 1:size(cellTCs_all,2);
    subplot(1,2,iday); hold on
    d = pct_events_all{iday};
    histogram(d,10);
    vline(mean(d));
    xlabel('pct events');
    ylabel('# cells');
    title(['Spontaneous events day ',num2str(iday)])
end
print(fullfile(fn_multi,'AllSpontaneousEvents'),'-dpdf')

%% plotting max df/f and correlation matrix with masks
figure;
imagesc(dfmax{1});
colormap gray
caxis([0.05 .3])
title('max df/f day 1');
hold on
bound = cell2mat(bwboundaries(masks{2}(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
bound = cell2mat(bwboundaries(masks{1}(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','b','MarkerSize',2);
hold off

figure;
imagesc(dfmax{3});
colormap gray
caxis([0.05 .3])
title('max df/f day 2 after rotation');
hold on
bound = cell2mat(bwboundaries(masks{2}(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
bound = cell2mat(bwboundaries(masks{1}(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','b','MarkerSize',2);
hold off

%% find the subset of masks I want
mask_match_d1 = zeros(size(corrmap{1}));
for i = 1:size(match_ind,2)
   ind = match_ind(i);
   temp_mask_inds = find(masks{1}==ind);
   mask_match_d1(temp_mask_inds)=1;
end


%find red masks for day 1
mask_red_d1 = zeros(size(corrmap{1}));
for i = 1:size(red_d1,2)
   ind = red_d1(i);
   temp_mask_inds = find(masks{1}==ind);
   mask_red_d1(temp_mask_inds)=1;
end


%find red masks matched on day 2
red_match=find(red_ind_match);

mask_red_match = zeros(size(corrmap{3}));

for i = 1:size(red_match,2)
   ind = red_match(i);
   temp_mask_inds = find(masks{2}==ind);
   mask_red_match(temp_mask_inds)=1;
end


%%
figure;
imagesc(corrmap{1});
colormap gray
%caxis([0.05 .3])
title('average FOV day 1');
hold on
bound = cell2mat(bwboundaries(mask_match_d1(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','b','MarkerSize',2);
bound = cell2mat(bwboundaries(mask_red_d1(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
hold off


%%
figure;
imagesc(corrmap{3});
colormap gray
%caxis([0.05 .3])
title('average FOV day 2');
hold on
bound = cell2mat(bwboundaries(masks{2}(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','b','MarkerSize',2);
bound = cell2mat(bwboundaries(mask_red_match(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
hold off