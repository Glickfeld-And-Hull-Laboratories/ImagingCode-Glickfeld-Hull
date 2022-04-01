clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};

rc = behavConstsDART; %directories
eval(ds);

day_id = 150; %enter post-DART day
pre_day = expt(day_id).multiday_matchdays;

nd=2; %hardcoding the number of days for now

mouse = expt(day_id).mouse;

fnout = fullfile(rc.achAnalysis,mouse);
if expt(day_id).multiday_timesincedrug_hours>0
    dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end

d=string(datetime('today'));
fn_multi_analysis = fullfile(rc.achAnalysis,mouse,['multiday_' dart_str],d);
mkdir(fn_multi_analysis);
fn_multi = fullfile(rc.achAnalysis,mouse,['multiday_' dart_str]);


load(fullfile(fn_multi,'tc_keep.mat'))
load(fullfile(fn_multi,'multiday_alignment.mat'))
load(fullfile(fn_multi,'resp_keep.mat'))
load(fullfile(fn_multi,'input.mat'))
load(fullfile(fn_multi,'locomotion.mat'))

cd(fn_multi)
nKeep = size(tc_trial_avrg_keep{1},2);
clear d
% find stimulus conditions
frame_rate = input.frameImagingRateMs;
%%

%tells the contrast, direction and orientation for each trial each day
tCon_match = cell(1,nd);
tDir_match = cell(1,nd);
tOri_match = cell(1,nd);

%find the contrasts, directions and orientations for each day
for id = 1:nd
    tCon_match{id} = celleqel2mat_padded(input(id).tGratingContrast(1:nTrials(id)));
    tDir_match{id} = celleqel2mat_padded(input(id).tGratingDirectionDeg(1:nTrials(id)));
    tOri_match{id} = tDir_match{id};
    tOri_match{id}(find(tDir_match{id}>=180)) = tDir_match{id}(find(tDir_match{id}>=180))-180;
end
oris = unique(tOri_match{1});
cons = unique(tCon_match{1});
nOri = length(oris);
nCon = length(cons);

nOn = input(1).nScansOn;
nOff = input(1).nScansOff;
%% make figure with se shaded, one figure per contrast

tc_green_avrg_keep = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_keep = cell(1,nd); %same for red
tc_green_se_keep = cell(1,nd); %this will be the se across all green cells
tc_red_se_match = cell(1,nd); %same for red

for id = 1:nd
    for iCon=1:nCon
    tc_green_avrg_keep{id}(:,iCon)=nanmean(tc_trial_avrg_keep{id}(:,green_ind_keep,iCon),2);
    green_std=std(tc_trial_avrg_keep{id}(:,green_ind_keep,iCon),[],2);
    tc_green_se_keep{id}(:,iCon)=green_std/sqrt(length(green_ind_keep));
    
    tc_red_avrg_keep{id}(:,iCon)=nanmean(tc_trial_avrg_keep{id}(:,red_ind_keep,iCon),2);
    red_std=std(tc_trial_avrg_keep{id}(:,red_ind_keep,iCon),[],2);
    tc_red_se_match{id}(:,iCon)=red_std/sqrt(length(red_ind_keep));
    
    clear green_std red_std
    end
end


%creat a time axis in seconds
t=1:(size(tc_green_avrg_keep{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day

shadedErrorBar(t,tc_red_avrg_keep{2}(:,iCon),tc_red_se_match{2}(:,iCon),'r');
ylim([-.02 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_keep{2}(:,iCon),tc_green_se_keep{2}(:,iCon));
title(['Pre-DART contrast = ' num2str(cons(iCon))])
txt1 = ['HT- ' num2str(length(green_ind_keep))];
text(-1.5,0.25,txt1);
txt2 = ['HT+ ' num2str(length(red_ind_keep))];
text(-1.5,0.23,txt2,'Color','r');
ylabel('dF/F') 
xlabel('s') 

axis square


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_keep{1}(:,iCon),tc_red_se_match{1}(:,iCon),'r');
ylim([-.02 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_keep{1}(:,iCon),tc_green_se_keep{1}(:,iCon));
ylabel('dF/F') 
xlabel('s') 
title(['Post-DART contrast = ' num2str(cons(iCon))])
axis square

print(fullfile(fn_multi_analysis,[num2str(cons(iCon)) '_timecourses.pdf']),'-dpdf');
end 
clear txt1 txt2
%% make a plot of individual timecourses 

setYmin = -.2; %indicate y axes you want
setYmax = 0.8;


for iCon = 1:nCon

figure
subplot(2,2,1)
plot(t, tc_trial_avrg_keep{2}(:,green_ind_keep,iCon),'color',[0 0 0 .2])
hold on
plot(t, tc_green_avrg_keep{2}(:,iCon),'color',[0 0 0],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 1')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,2)
plot(t, tc_trial_avrg_keep{2}(:,red_ind_keep,iCon),'color',[.7 .05 .05 .2])
hold on
plot(t, tc_red_avrg_keep{2}(:,iCon),'color',[.7 .05 .05],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 1')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,3)
plot(t, tc_trial_avrg_keep{1}(:,green_ind_keep,iCon),'color',[0 0 0 .2])
hold on
plot(t, tc_green_avrg_keep{1}(:,iCon),'color',[0 0 0],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 2')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,4)
plot(t, tc_trial_avrg_keep{1}(:,red_ind_keep,iCon),'color',[.7 .05 .05 .2])
hold on
plot(t, tc_red_avrg_keep{1}(:,iCon),'color',[.7 .05 .05],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 2')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
sgtitle(['contrast = '  num2str(cons(iCon))])
hold off
print(fullfile(fn_multi_analysis,[num2str(cons(iCon)) '_indiv_timecourses.pdf']),'-dpdf');
end

%% ploting pre-and post dart average timecourse seperately for each cell
if nKeep<36
    [n n2] = subplotn(nKeep);
    tot = n.*n2;
else
    n = 6;
    n2 = 6;
    tot = 36;
end

setYmin = -.2; %indicate y axes you want
setYmax = 0.5;

figure;
movegui('center')
start = 1;
for iCell = 1:nKeep
    
    if start>tot
        figure; movegui('center')
        start = 1;
    end
    subplot(n,n2,start)
    if ismember(iCell,red_ind_keep)
        plot(t,tc_trial_avrg_keep_allCon{2}(:,iCell),'r')
        hold on
        plot(t,tc_trial_avrg_keep_allCon{1}(:,iCell),'--r')
        hold off
        ylim([setYmin setYmax]);
    else
        plot(t,tc_trial_avrg_keep_allCon{2}(:,iCell),'k')
        hold on
        plot(t,tc_trial_avrg_keep_allCon{1}(:,iCell),'--k')
        hold off
        ylim([setYmin setYmax]);
    end
    start = start+1;
end


%% time to peak
tHalfMax = cell(1,nd);

for id = 1:nd
    tHalfMaxTemp=nan(nKeep,1);

        for iCell=1:nKeep
            %pull the data for a given cell at a given contrast (pref ori)
            tempData=tc_trial_avrg_keep_allCon{id}(stimStart:stimStart+nOn,iCell);
            if resp_keep{id}(iCell)
                smoothData=smoothdata(tempData,'movmean',5) ;
                halfMax = max(smoothData(3:length(smoothData)))/2;
                tHalfMaxCell =double(min(find(smoothData>halfMax)))/double(frame_rate);
                if rem(iCell, 10) == 0 
                figure;plot(tempData)
                hold on
                plot(smoothData)
                hold on
                vline(min(find(smoothData>halfMax)))
                end
                if length(tHalfMaxCell)>0
                    tHalfMaxTemp(iCell)=tHalfMaxCell;
                end
            end
 
       end
    tHalfMax{id}=tHalfMaxTemp;
end
clear tHalfMaxCell tHalfMaxTemp tempData smoothData halfMax
% scatter for tMax



figure; movegui('center') 
subplot(1,2,1)
scatter((tHalfMax{2}(green_ind_keep)),(tHalfMax{1}(green_ind_keep)),'k')
ylabel('post-DART half-max(s)')
xlabel('pre-DART half-max(s)')
ylim([0 2.5])
xlim([0 2.5])
refline(1)
title('HT- ')
axis square


subplot(1,2,2)
scatter((tHalfMax{2}(red_ind_keep)),(tHalfMax{1}(red_ind_keep)),'MarkerEdgeColor',[.7 .05 .05])

ylabel('post-DART half-max(s)')
xlabel('pre-DART half-max(s)')
ylim([0 2.5])
xlim([0 2.5])


refline(1)
title('HT+')
axis square

sgtitle('time to half max (s), averaged over contrast')
print(fullfile(fn_multi_analysis, ['tHalfMax_crossDay.pdf']),'-dpdf','-bestfit')

%% scatterplot of max df/f for day 1 vs day 2, and each subplot is one cell type


for iCon = 1:nCon
figure; movegui('center') 
subplot(1,2,1)
scatter((resp_max_keep{2}(green_ind_keep,iCon)),(resp_max_keep{1}(green_ind_keep,iCon)),'k')
hold on
scatter((resp_max_keep{2}(intersect(green_ind_keep,find(sig_diff{iCon})),iCon)),(resp_max_keep{1}(intersect(green_ind_keep,find(sig_diff{iCon})),iCon)),'k','MarkerFaceColor','k')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])
refline(1)
title('HT- ')
axis square


subplot(1,2,2)
scatter((resp_max_keep{2}(red_ind_keep,iCon)),(resp_max_keep{1}(red_ind_keep,iCon)),'MarkerEdgeColor',[.7 .05 .05])
hold on
scatter((resp_max_keep{2}(intersect(red_ind_keep,find(sig_diff{iCon})),iCon)),(resp_max_keep{1}(intersect(red_ind_keep,find(sig_diff{iCon})),iCon)),'MarkerEdgeColor',[.7 .05 .05],'MarkerFaceColor',[.7 .05 .05])
hold off
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])

refline(1)
title('HT+')
axis square

sgtitle(num2str(cons(iCon)))
print(fullfile(fn_multi_analysis,[num2str(cons(iCon)) 'maxResp_crossDay.pdf']),'-dpdf','-bestfit')
end
%%
figure;
for iCon = 1:nCon
subplot(1,nCon,iCon)
for iRed = 1:length(red_ind_keep)
    thisCell = red_ind_keep(iRed);
    scatter((resp_max_keep{2}(thisCell,iCon)),(resp_max_keep{1}(thisCell,iCon)))
hold on
end
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])

refline(1)
axis square
title(num2str(cons(iCon)))

end
print(fullfile(fn_multi_analysis,[num2str(cons(iCon)) 'maxResp_HTCellsColored.pdf']),'-dpdf','-bestfit')
%% 
for iCon = 1:nCon
figure;
subplot(1,2,1)
scatter((resp_max_keep{2}(red_ind_keep,iCon)),dfof_max_diff_raw(red_ind_keep,iCon),'MarkerEdgeColor',[.7 .05 .05],'MarkerFaceColor',[.7 .05 .05],'MarkerFaceAlpha',.25)
lsline;
hold on
scatter((resp_max_keep{2}(intersect(red_ind_keep,find(sig_diff{iCon})),iCon)),dfof_max_diff_raw(intersect(red_ind_keep,find(sig_diff{iCon})),iCon),'MarkerFaceColor',[.7 .05 .05])
[R,p]=corrcoef((resp_max_keep{2}(red_ind_keep,iCon)),dfof_max_diff_raw(red_ind_keep,iCon),'Rows','complete');
title(['R = '  num2str(R(2)) ', p = ' num2str(p(2))])
ylabel('pre-DART - post-DART dF/F')
xlabel('pre-DART  dF/F')
axis square
hold off

subplot(1,2,2)
scatter((resp_max_keep{2}(red_ind_keep,iCon)),dfof_max_diff(red_ind_keep,iCon),'MarkerEdgeColor',[.7 .05 .05],'MarkerFaceColor',[.7 .05 .05],'MarkerFaceAlpha',.25)
lsline;
hold on
scatter((resp_max_keep{2}(intersect(red_ind_keep,find(sig_diff{iCon})),iCon)),dfof_max_diff(intersect(red_ind_keep,find(sig_diff{iCon})),iCon),'MarkerFaceColor',[.7 .05 .05])
[R,p]=corrcoef((resp_max_keep{2}(red_ind_keep,iCon)),dfof_max_diff(red_ind_keep,iCon),'Rows','complete');
title(['R = '  num2str(R(2)) ', p = ' num2str(p(2))])
ylabel('(post-pre)/(post+pre)')
xlabel('pre-DART  dF/F')
axis square
sgtitle(num2str(cons(iCon)))
print(fullfile(fn_multi_analysis,[num2str(cons(iCon)) 'd1Resp_vs_change.pdf']),'-dpdf','-bestfit')
clear R p
end


%same as plot above but comparing contrasts across days rather than days
%across contrasts. NOTE this is for two contrasts, 50% and 100%. Would need
%to change this for additional contrasts
figure;

subplot(1,2,1)
scatter(dfof_max_diff(red_ind_keep,1),dfof_max_diff(red_ind_keep,2),'MarkerEdgeColor',[.7 .05 .05],'MarkerFaceColor',[.7 .05 .05],'MarkerFaceAlpha',.25)
ylim([-1 1])
xlim([-1 1])
lsline;
[R,p]=corrcoef(dfof_max_diff(red_ind_keep,1),dfof_max_diff(red_ind_keep,2),'Rows','complete');
title(['R = '  num2str(R(2)) ', p = ' num2str(p(2))])
xlabel('Fractional change contrast 0.5')
ylabel('Fractional change contrast 1')
refline(1)
axis square

subplot(1,2,2)
scatter(dfof_max_diff_raw(red_ind_keep,1),dfof_max_diff_raw(red_ind_keep,2),'MarkerEdgeColor',[.7 .05 .05],'MarkerFaceColor',[.7 .05 .05],'MarkerFaceAlpha',.25)
ylim([-.35 .35])
xlim([-.35 .35])
lsline;
[R,p]=corrcoef(dfof_max_diff_raw(red_ind_keep,1),dfof_max_diff_raw(red_ind_keep,2),'Rows','complete');
title(['R = '  num2str(R(2)) ', p = ' num2str(p(2))])
xlabel('Raw change contrast 0.5')
ylabel('Raw change contrast 1')
refline(1)
axis square

print(fullfile(fn_multi_analysis,[num2str(cons(iCon)) 'changeVsContrast.pdf']),'-dpdf','-bestfit')
clear R p



%% fractional change in dfof

for iCon = 1:nCon
figure
x = [nanmean(dfof_max_diff(green_ind_keep,iCon)), nanmean(dfof_max_diff(red_ind_keep,iCon))];
y = [(std(dfof_max_diff(green_ind_keep,iCon)))/sqrt(length(green_ind_keep)), (std(dfof_max_diff(red_ind_keep,iCon)))/sqrt(length(red_ind_keep))];

labs =categorical({'HT-','HT+'});
bar(labs,x)                
hold on
er = errorbar(labs,x,-y,y);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
%ylim([0 .2])
hold off
title(['fractional change dfof, contrast = ', num2str(cons(iCon))])

print(fullfile(fn_multi_analysis,[num2str(cons(iCon)), '_frac_change_resp.pdf']),'-dpdf','-bestfit')

end

%%

%now I know which trials the mouse was running on, need to look at LMI
%LMI = (R_loc - R_stat)/(R_loc + R_stat)
%find the intersection of trials at the preferred orientation during
%locomotion vs stationary periods - I may not have enough trials of each
%kind

%creat a time axis in seconds
t=1:(size(tc_green_avrg_keep{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
    %make this the average for running state for each contrast
    figure;
 
        
        locr_mean = nanmean(locTCs{2,iCon}(:,red_ind_keep),2); %averaged across red cells
        locr_std = std(locTCs{2,iCon}(:,red_ind_keep),[],2);
        locr_se=locr_std/sqrt(length(red_ind_keep));
        statr_mean = nanmean(statTCs{2,iCon}(:,red_ind_keep),2);
        statr_std = std(statTCs{2,iCon}(:,red_ind_keep),[],2);
        statr_se=statr_std/sqrt(length(red_ind_keep));

        locg_mean = nanmean(locTCs{2,iCon}(:,green_ind_keep),2);
        locg_std = std(locTCs{2,iCon}(:,green_ind_keep),[],2);
        locg_se=locg_std/sqrt(length(green_ind_keep));

        statg_mean = nanmean(statTCs{2,iCon}(:,green_ind_keep),2);
        statg_std = std(statTCs{2,iCon}(:,green_ind_keep),[],2);
        statg_se=statg_std/sqrt(length(green_ind_keep));

        subplot(1,2,1)
        shadedErrorBar(t,locr_mean,locr_se,'r')
        hold on
        shadedErrorBar(t,statr_mean,statr_se,'--r')
        shadedErrorBar(t,locg_mean,locg_se,'k')
        shadedErrorBar(t,statg_mean,statg_se,'--k')
        ylim([-.05 .2]);
        title(['contrast ', num2str(cons(iCon)), ' pre-DART'])
        txt1=[num2str(locCounts{2,iCon}(1)), ' running trials']
        text(-1.5,0.18,txt1);
        txt2 = [num2str(locCounts{2,iCon}(2)), ' stationary trials']
        text(-1.5,0.16,txt2);
        axis square
        hold off   
        
 
        
        locr_mean = nanmean(locTCs{1,iCon}(:,red_ind_keep),2); %averaged across red cells
        locr_std = std(locTCs{1,iCon}(:,red_ind_keep),[],2);
        locr_se=locr_std/sqrt(length(red_ind_keep));
        statr_mean = nanmean(statTCs{1,iCon}(:,red_ind_keep),2);
        statr_std = std(statTCs{1,iCon}(:,red_ind_keep),[],2);
        statr_se=statr_std/sqrt(length(red_ind_keep));

        locg_mean = nanmean(locTCs{1,iCon}(:,green_ind_keep),2);
        locg_std = std(locTCs{1,iCon}(:,green_ind_keep),[],2);
        locg_se=locg_std/sqrt(length(green_ind_keep));

        statg_mean = nanmean(statTCs{1,iCon}(:,green_ind_keep),2);
        statg_std = std(statTCs{1,iCon}(:,green_ind_keep),[],2);
        statg_se=statg_std/sqrt(length(green_ind_keep));
        
        subplot(1,2,2)
        shadedErrorBar(t,locr_mean,locr_se,'r')
        hold on
        shadedErrorBar(t,statr_mean,statr_se,'--r')
        shadedErrorBar(t,locg_mean,locg_se,'k')
        shadedErrorBar(t,statg_mean,statg_se,'--k')
        ylim([-.05 .2]);
        title(['contrast ', num2str(cons(iCon)), ' post-DART'])
        txt1=[num2str(locCounts{1,iCon}(1)), ' running trials']
        text(-1.5,0.18,txt1);
        txt2 = [num2str(locCounts{1,iCon}(2)), ' stationary trials']
        text(-1.5,0.16,txt2);
        axis square
        hold off
        
       clear locr_mean locr_std locr_se statr_mean statr_std statr_se locg_mean locg_std locg_se statg_mean statg_std statg_se txt1 txt2
       print(fullfile(fn_multi_analysis,[num2str(cons(iCon)) '_locomotionTCs.pdf']),'-dpdf');
  
end

%% locomotion modulation index


for iCon = 1:nCon

figure; movegui('center') 
subplot(1,2,1)
scatter((LMI{2,iCon}(green_ind_keep)),(LMI{1,iCon}(green_ind_keep)),'k')
ylabel('post-DART LMI')
xlabel('pre-DART  LMI')
ylim([-1 1])
xlim([-1 1])
hline(0)
vline(0)
refline(1)
title('HT- ')
axis square


subplot(1,2,2)
scatter((LMI{2,iCon}(red_ind_keep)),(LMI{1,iCon}(red_ind_keep)),'MarkerEdgeColor',[.7 .05 .05])
ylabel('post-DART LMI')
xlabel('pre-DART  LMI')
ylim([-1 1])
xlim([-1 1])
hline(0)
vline(0)
refline(1)
title('HT+')
axis square

sgtitle(num2str(cons(iCon)))
print(fullfile(fn_multi_analysis,[num2str(cons(iCon)) '_LMI.pdf']),'-dpdf');
end

%% finding cells that are still saturated by contrast vs. not saturated pre-DART

%identify cells that are still in the rise of the contrast response
%function
risingCells = find(data_con_resp_keep{2}(:,1)<data_con_resp_keep{2}(:,2));

% make a plot of individual timecourses for rising cells
setYmin = -.1; %indicate y axes you want
setYmax = 0.8;

for iCon = 1:nCon

figure
subplot(2,2,1)
plot(t, tc_trial_avrg_keep{2}(:,intersect(risingCells,green_ind_keep),iCon),'k')
ylim([setYmin setYmax]);
title('day 1')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 


subplot(2,2,2)
plot(t, tc_trial_avrg_keep{2}(:,intersect(risingCells,red_ind_keep),iCon),'color',[.7 .05 .05])
ylim([setYmin setYmax]);
title('day 1')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 


subplot(2,2,3)
plot(t, tc_trial_avrg_keep{1}(:,intersect(risingCells,green_ind_keep),iCon),'k')
ylim([setYmin setYmax]);
title('day 2')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 


subplot(2,2,4)
plot(t, tc_trial_avrg_keep{1}(:,intersect(risingCells,red_ind_keep),iCon),'color',[.7 .05 .05])
ylim([setYmin setYmax]);
title('day 2')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
sgtitle(['contrast = '  num2str(cons(iCon))])
print(fullfile(fn_multi_analysis,[num2str(cons(iCon)) 'risingCells_indiv_timecourses.pdf']),'-dpdf');
end

%% making mask maps for various measurements

load(fullfile(fn_multi,'mask_measuremens.mat'))

% spatial distribution vs. fractional change DART effect
figure;
imagesc(keep_masks_fract_change_red)
colorbar
title('Spatial distribution of cells by fractional change from DART, HT+')
caxis([-1 1])
hold on
bound = cell2mat(bwboundaries(keep_red_masks(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','k','MarkerSize',2);
hold off
colormap bluered
axis square
print(fullfile(fn_multi_analysis,'HT+_frac_change_map.pdf'),'-dpdf','-bestfit')

figure;
imagesc(keep_masks_fract_change_green)
colorbar
title('Spatial distribution of cells by fractional change from DART, HT-')
caxis([-1 1])
hold on
bound = cell2mat(bwboundaries(keep_green_masks(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','k','MarkerSize',2);
hold off
colormap bluered
axis square
print(fullfile(fn_multi_analysis,'HT-_frac_change_map.pdf'),'-dpdf','-bestfit')


% spatial distribution vs. raw DART effect
imagesc(keep_masks_raw_change_red)
colorbar
title('Spatial distribution of cells by raw change from DART, HT+')
caxis([-.1 .1])
hold on
bound = cell2mat(bwboundaries(keep_red_masks(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','k','MarkerSize',2);
hold off
colormap bluered
axis square
print(fullfile(fn_multi_analysis,'HT+_raw_change_map.pdf'),'-dpdf','-bestfit')

figure;
imagesc(keep_masks_raw_change_green)
colorbar
title('Spatial distribution of cells by raw change from DART, HT-')
caxis([-.5 .5])
hold on
bound = cell2mat(bwboundaries(keep_green_masks(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','k','MarkerSize',2);
hold off
colormap bluered
axis square
print(fullfile(fn_multi_analysis,'HT-_raw_change_map.pdf'),'-dpdf','-bestfit')

% map of kept HT+ cells by their pre-DART df/F
imagesc(keep_masks_d1_red)
colorbar
title('Spatial distribution of cells by raw change from DART, HT+')
caxis([-.14 .14])
hold on
bound = cell2mat(bwboundaries(keep_red_masks(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','k','MarkerSize',2);
hold off
colormap bluered
axis square
print(fullfile(fn_multi_analysis,'HT+_d1_dfof_map.pdf'),'-dpdf','-bestfit')
%% plotting  ori response 
if nKeep<36
    [n n2] = subplotn(nKeep);
    tot = n.*n2;
else
    n = 4;
    n2 = 4;
    tot = 16;
end

% plot ori tuning at each contrast
figure;
movegui('center')
start = 1;

for iCell = 1:nKeep
    for id = 1:nd
        if start>tot
            figure; movegui('center')
            start = 1;
        end
        subplot(n,n2,start)

        for iCon = 1:nCon
            errorbar(oris, data_resp_keep{id}(iCell,:,iCon,1), data_resp_keep{id}(iCell,:,iCon,2),'-o')
            hold on
        end
        if ismember(iCell,red_ind_keep)
            extra_title='HT+';
        else
            extra_title='HT-';
        end
            start= start+1;
        ylim([-0.1 inf])
        title(['Day ' num2str(id) extra_title])
    end

end

%% plot contrast response averaged over orientations 
if nKeep<36
    [n n2] = subplotn(nKeep);
    tot = n.*n2;
else
    n = 5;
    n2 = 5;
    tot = 25;
end

figure;
movegui('center')
start = 1;
titleNum = 1;

for iCell = 1:nKeep
    if ismember(iCell,red_ind_keep)
        if start>tot

            figure; movegui('center')
           
            start = 1;
            
        end
        subplot(n,n2,start)
        errorbar(cons, squeeze(nanmean(data_resp_keep{2}(iCell,:,:,1),2)), squeeze(nanmean(data_resp_keep{2}(iCell,:,:,2),2)),'-')
        hold on
        errorbar(cons, squeeze(nanmean(data_resp_keep{1}(iCell,:,:,1),2)), squeeze(nanmean(data_resp_keep{1}(iCell,:,:,2),2)),'-')
        start= start+1;
        %ylim([-.05 .2])
    end     
    sgtitle('HT+')
    
end

figure;
movegui('center')
start = 1;
titleNum = 1;

for iCell = 1:nKeep
    if ~ismember(iCell,red_ind_keep)
        if start>tot
 
            figure; movegui('center')
           
            start = 1;
            
        end
        subplot(n,n2,start)
        errorbar(cons, squeeze(nanmean(data_resp_keep{2}(iCell,:,:,1),2)), squeeze(nanmean(data_resp_keep{2}(iCell,:,:,2),2)),'-')
        hold on
        errorbar(cons, squeeze(nanmean(data_resp_keep{1}(iCell,:,:,1),2)), squeeze(nanmean(data_resp_keep{1}(iCell,:,:,2),2)),'--')

        start= start+1;
        %ylim([-.05 .2])
    end     
    sgtitle('HT-')
    
end

%% plot all responses 
if nKeep<36
    [n n2] = subplotn(nKeep);
    tot = n.*n2;
else
    n = 5;
    n2 = 5;
    tot = 25;
end

figure;
movegui('center')
start = 1;

for iCell = 1:nKeep
     
        if start>tot

            figure; movegui('center')
           
            start = 1;
            
        end
        subplot(n,n2,start)
            for iCon = 1:nCon
                scatter(data_resp_keep{2}(iCell,:,iCon,1),data_resp_keep{1}(iCell,:,iCon,1))
                %the contrasts are plotted indiviudally so that the points
                %corresponding to different contrasts are different colors,
                hold on
                lsline
            end
                if ismember(iCell,red_ind_keep)
                    title('\color{red}HT+')
                else
                    title('HT-')
                end
            %max_list and upperLim are to find a single axis upper limit
            %that will work for botht eh x and y axes
            max_list = [max(max(data_resp_keep{2}(iCell,:,:,1))),max(max(data_resp_keep{1}(iCell,:,:,1)))];
            upperLim =max(max_list)+0.1;
            ylim([-.05 upperLim]);
            xlim([-.05 upperLim]);
            axis square
            hold off
            start=start+1;
                
end
%%
linCell=zeros(1,nKeep);%will provide a boolean index of which cells have a linear relationship
linCellProps = nan(4,nKeep); %this will contain the slopes and intercepts of cells that are deemed linear
data_pre_all=nan(nOri*nCon, nKeep);
count=1;
figure;
for iCell = 1:nKeep
    if ismember(iCell,red_ind_keep)
 data_pre=[];
 
 data_post=[];
  for iCon = 1:nCon
                data_pre=[data_pre,data_resp_keep{2}(iCell,:,iCon,1)]; %collects the values for day 1 over all contrasts
                data_post=[data_post,data_resp_keep{1}(iCell,:,iCon,1)]; %collects the values for day 2 over all contrasts
  end
  data_pre_all(:,iCell)=data_pre;
   [R,p]=corrcoef(data_pre,data_post); %gets the correlation data for
   %the correlation between days, pooling across contrasts
    if R > 0.5
        linCell(iCell)=1;
        linFit = polyfit(data_pre, data_post, 1);
        linCellProps(1,iCell)=linFit(1); %slope
        linCellProps(2,iCell)=linFit(2); %intercept
        linCellProps(3,iCell)=R(2);
        linCellProps(4,iCell)=p(2);
        f=polyval(linFit,data_pre);
        plot(data_pre,data_post,'ko',data_pre,f,'k-','LineWidth',.5)
        
        count=count+1;

    end
    hold on    
    end

end
% collecting info for the linear cells
        
meanFit= [mean(linCellProps(1,find(linCell))),  mean(linCellProps(2,find(linCell)))];
xRange=[min(min(data_pre_all)), max(max(data_pre_all))]
f2=polyval(meanFit,xRange);
plot(xRange,f2,'b-','LineWidth',2)
txt1 = ['y=', num2str(round(meanFit(1),2)),'x',num2str(round(meanFit(2),4))];
text(.15,0.01,txt1);
hold off
axis square
myRef = refline(1)
myRef.LineStyle = ':'

title(['\color{red}HT+ cells, all stimulus conditions, n= ' num2str(length(intersect(red_ind_keep,find(linCell))))])
xlabel('pre-DART')
ylabel('post-DART')
myRef = refline(1)
myRef.LineStyle = ':'
print(fullfile(fn_multi_analysis,'HT-_allResp.pdf'),'-dpdf','-bestfit')


% same as above for HT- cells


linCell=zeros(1,nKeep);%will provide a boolean index of which cells have a linear relationship
linCellProps = nan(4,nKeep); %this will contain the slopes and intercepts of cells that are deemed linear
data_pre_all=nan(nOri*nCon, nKeep);
count=1;
figure;
for iCell = 1:nKeep
    if ~ismember(iCell,red_ind_keep)
 data_pre=[];
 
 data_post=[];
  for iCon = 1:nCon
                data_pre=[data_pre,data_resp_keep{2}(iCell,:,iCon,1)]; %collects the values for day 1 over all contrasts
                data_post=[data_post,data_resp_keep{1}(iCell,:,iCon,1)]; %collects the values for day 2 over all contrasts
  end
  data_pre_all(:,iCell)=data_pre;
   [R,p]=corrcoef(data_pre,data_post); %gets the correlation data for
   %the correlation between days, pooling across contrasts
    if R > 0.5
        linCell(iCell)=1;
        linFit = polyfit(data_pre, data_post, 1);
        linCellProps(1,iCell)=linFit(1); %slope
        linCellProps(2,iCell)=linFit(2); %intercept
        linCellProps(3,iCell)=R(2);
        linCellProps(4,iCell)=p(2);
        f=polyval(linFit,data_pre);
        plot(data_pre,data_post,'ko',data_pre,f,'k-','LineWidth',.5)
        
        count=count+1;

    end
    hold on    
    end

end
% collecting info for the linear cells
        
meanFit= [mean(linCellProps(1,find(linCell))),  mean(linCellProps(2,find(linCell)))];
xRange=[min(min(data_pre_all)), max(max(data_pre_all))]
f2=polyval(meanFit,xRange);
plot(xRange,f2,'b-','LineWidth',2)
txt1 = ['y=', num2str(round(meanFit(1),2)),'x',num2str(round(meanFit(2),4))];
text(.2,-0.01,txt1);
hold off
axis square
myRef = refline(1)
myRef.LineStyle = ':'

title(['HT- cells, all stimulus conditions, n= ' num2str(length(intersect(green_ind_keep,find(linCell))))])
xlabel('pre-DART')
ylabel('post-DART')
myRef = refline(1)
myRef.LineStyle = ':'

print(fullfile(fn_multi_analysis,'HT+_allResp.pdf'),'-dpdf','-bestfit')


%% plot observed pref ori 
%the orientation that gave the strongest response in the observed data

fig=figure; movegui('center') 
scatter(pref_ori_keep{2}(green_ind_keep),pref_ori_keep{1}(green_ind_keep),'k','jitter', 'on', 'jitterAmount', 5)
hold on
scatter(pref_ori_keep{2}(red_ind_keep),pref_ori_keep{1}(red_ind_keep),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 5)
hold off
xlabel('Pre-DART pref ori')
ylabel('Post-DART pref ori')
% xlim([0 .4])
% ylim([0 .4])
refline(1)
title('Observed pref ori')
saveas(fig, 'obsvOri.png')

ori_diff_keep = abs(pref_ori_keep{2}-pref_ori_keep{1});
figure
boxplot(ori_diff_keep, red_keep_logical) %this isn't the best way to show this

%% extract fit orientation preference

%1 initial von M fit on real data/real orientations collected
vonMises_output = cell(1,nd);

for id = 1:nd
    vonMises_output_temp=nan(nKeep, 6);
    for i = 1:nKeep
        thisCell = data_ori_resp_keep{id}(i,:);
        [b_hat, k1_hat, R1_hat,u1_hat,sse,R_square] = miaovonmisesfit_ori(deg2rad(oris),thisCell);
        vonMises_output_temp(i,:)=[b_hat, k1_hat, R1_hat,u1_hat,sse,R_square];
    end
    vonMises_output{id}=vonMises_output_temp;

end
clear thisCell b_hat k1_hat R1_hat u1_hat sse R_square vonMises_output_temp %get rid of the last iteration of these looping variables
save(fullfile(fn_multi_analysis,'vonM_output.mat'),'vonMises_output')

%I can extract sharpness values from this structure, it is the second
%column


%2 fit with higher resolution
fit_oris = (0:1:179);

y_fits_keep = cell(1,nd); %an n cell X 180 degree matrix for each day representing the fit tuning curve
fit_pref_oris_keep = cell(1,nd); %a matrix for each day with the fit preferred ori

for id = 1:nd
    y_fits = zeros(nKeep, size(fit_oris,2));
    pref_oris = nan(nKeep, 1);
    for i = 1:nKeep
        in = vonMises_output{id}(i,:);
        b_tmp = in(1);
        k1_tmp = in(2);
        R1_tmp = in(3);
        u1_tmp = in(4);
        y_fits(i,:) = b_tmp+R1_tmp.*exp(k1_tmp.*(cos(2.*(deg2rad(fit_oris)-u1_tmp))-1));
        if size(fit_oris(find(y_fits(i,:)==max(y_fits(i,:)))),2)>1
            fit_pref_oris(i)=nan;
        else
            fit_pref_oris(i)= fit_oris(find(y_fits(i,:)==max(y_fits(i,:))));
        end
    end
    y_fits_keep{id}=y_fits;
    fit_pref_oris_keep{id}=fit_pref_oris;
clear b_tmp k1_tmp R1_tmp u1_tmp y_fits fit_pref_oris   
end

%% plots observed tuning curve with von M fit curve


figure;
movegui('center')
start = 1;

for iCell = 1:nKeep
        if start>tot
            figure; movegui('center')
            start = 1;
        end
        subplot(n,n2,start)
        for  id = 1:nd
            plot(oris, data_ori_resp_keep{id}(iCell,:))
            hold on
            plot(fit_oris,y_fits_keep{id}(iCell,:),':');
            hold on
        end
        
        if ismember(iCell,red_ind_keep)
            extra_title='HT+';
        else
            extra_title='HT-';
        end
            start= start+1;
        ylim([-0.1 inf])
        title([extra_title])
        
end


%% compare Von M fit pref ori across days
fig=figure; movegui('center') 
scatter(fit_pref_oris_keep{2}(green_ind_keep),fit_pref_oris_keep{1}(green_ind_keep),'k','jitter', 'on', 'jitterAmount', 5)
hold on
scatter(fit_pref_oris_keep{2}(red_ind_keep),fit_pref_oris_keep{1}(red_ind_keep),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 5)
hold off
xlabel('D1- pref ori')
ylabel('D2- pref ori')
% xlim([0 .4])
% ylim([0 .4])
refline(1)
title('fit pref ori')
%saveas(fig, 'fitOri.png')

%% Von Mises bootstrap
reps = 1000;

shuff_pref_ori= cell(1, nd);
for id = 1:nd
    shuff_pref_ori{id} = nan(reps,nKeep);
    for r = 1:reps
        %make random sample - this will be redone every loop
        tCon = tCon_match{id}(:,1:nTrials);
        tDir = tDir_match{id}(:,1:nTrials);
        for i = 1:nKeep
            shuff_resp=nan(1,nOri);
            for iOri = 1:nOri
                con_inds=find(tCon==(pref_con_keep{id}(i)));
                ori_inds=find(tDir == oris(iOri));
                intersect_inds=intersect(con_inds,ori_inds);
                x=round(.8*length(intersect_inds));
                trial_inds = randsample(intersect_inds, x,1); %randomly select 80% of the trials at this contrast and ori
                temp_TCs=data_trial_keep{id}(resp_win,trial_inds,i); %extract the response window of these trials for this cell
               shuff_resp(1,iOri) = mean(squeeze(mean(temp_TCs)));
            end
            %now I can use the shuffled response by ori data to calculate
            %the von M fit for this cell on this rep
            
            [b_hat, k1_hat, R1_hat,u1_hat,sse,R_square] = miaovonmisesfit_ori(deg2rad(oris),shuff_resp);
            this_fits = b_hat+R1_hat.*exp(k1_hat.*(cos(2.*(deg2rad(fit_oris)-u1_hat))-1));
            this_pref = fit_oris(min(find(this_fits==max(this_fits))));
            shuff_pref_ori{id}(r,i)=this_pref;
        end
    end
end
save(fullfile(fnout,'bootstrap_fit_ori.mat'),'shuff_pref_ori');
clear shuff_resp this_pref this_fits b_hat k1_hat R1_hat u1_hat sse R_square x temp_TCs intersect_inds ori_inds con_inds iOri i r id 
%% compare von M bootstraps to original fit 
diff_pref = cell(1,nd);
well_fit = cell(1,nd);
for id = 1:nd
    diff_pref{id}=sort(abs(shuff_pref_ori{id} - fit_pref_oris_keep{id}),1);
    well_fit{id} = find(diff_pref{id}((reps*.9),:)<45); %set to current jump size in orientation
end
%well_fit gives a list of the well-fit cell
% this compares the fit pref ori across days for cells that were well-fit on the
% pre-DART day 
[R_g p_g] = corrcoef(fit_pref_oris_keep{2},fit_pref_oris_keep{1})
figure; movegui('center') 
scatter(fit_pref_oris_keep{2}((intersect(well_fit{2},green_ind_keep))),fit_pref_oris_keep{1}((intersect(well_fit{2},green_ind_keep))),'k','jitter', 'on', 'jitterAmount', 5)
hold on
scatter(fit_pref_oris_keep{2}(intersect(well_fit{2},red_ind_keep)),fit_pref_oris_keep{1}(intersect(well_fit{2},red_ind_keep)),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 5)
hold off
xlabel('D1- pref ori')
ylabel('D2- pref ori')
refline(1)
title('fit pref ori for cells that were well-fit pre-DART')

%%

figure; movegui('center') 
subplot(1,2,1)
scatter(fit_pref_oris_keep{1}(green_ind_keep),pref_ori_keep{1}(green_ind_keep),'k','jitter', 'on', 'jitterAmount', 5)
hold on
scatter(fit_pref_oris_keep{1}(red_ind_keep),pref_ori_keep{1}(red_ind_keep),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 5)
hold off
xlabel('fit pref ori')
ylabel('raw pref ori')
title('day 1')
refline(1)

subplot(1,2,2)
scatter(fit_pref_oris_keep{2}(green_ind_keep),pref_ori_keep{2}(green_ind_keep),'k','jitter', 'on', 'jitterAmount', 5)
hold on
scatter(fit_pref_oris_keep{2}(red_ind_keep),pref_ori_keep{2}(red_ind_keep),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 5)
hold off
xlabel('fit pref ori')
ylabel('raw pref ori')
title('day 2')
refline(1)
%% example cell tcs - this is to pull out some individual example cell traces
%
cellList=[red_ind_keep]; %enter the cells you're interested in by their index wihtin the keep dataframe

place=1;
figure
for i=1:length(cellList)
    for id = 1:nd
        thisCell = cellList(i)
        %only pulling from dfof data of keep cells
        tCon=tCon_match{id}(1:nTrials);
        tDir=tDir_match{id}(1:nTrials);
        %identify the trials where ori = pref ori
        temp_ori= pref_ori_keep{id}(thisCell); %find the preferred ori of this cell and convert to degrees
        ori_inds = find(tDir==temp_ori); %these are the trials at that ori
        temp_con = pref_con_keep{id}(thisCell);%find the preferred contrast of this cell and convert to contrast value
        con_inds=find(tCon==temp_con);
        temp_trials = intersect(ori_inds, con_inds);
        temp_trials(temp_trials==160)=[]
        temp_TCs=data_trial_keep{id}(:,temp_trials,thisCell); %pulls the selected trials for the selected cell. Shape is frames X trial
        thisCellMean = mean(temp_TCs,2);
        thisCellSE=std(temp_TCs')/sqrt(length(temp_trials));
        subplot(length(cellList),nd,place)
        if red_keep_logical(thisCell)==1
            shadedErrorBar(t,thisCellMean,thisCellSE,'r')
            ylim([-.2 .5]);
        else
            shadedErrorBar(x,thisCellMean,thisCellSE,'k')
            ylim([-.2 .5]);
        end
        title(string(thisCell));
        line([0,2],[.45,.45])
    place=place+1;
    end
    
end


%%
figure;
imagesc(corrmap{3})
colormap gray
title('corr')
hold on
bound = cell2mat(bwboundaries(keep_red_masks(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
caxis([200 8000])
hold off



figure;
imagesc(dfmax{3})
colormap gray
title('dfmax')
hold on
bound = cell2mat(bwboundaries(keep_red_masks(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
caxis([-.2 1])
hold off

figure;
imagesc(fov_avg{3})
colormap gray
title('FOV avrg')
hold on
bound = cell2mat(bwboundaries(keep_red_masks(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
caxis([200 8000])
hold off


figure;
imagesc(fov_red{3})
colormap gray
title('FOV red')
hold on
bound = cell2mat(bwboundaries(keep_red_masks(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);

caxis([10 100])
hold off
