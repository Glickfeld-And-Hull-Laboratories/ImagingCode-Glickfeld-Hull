clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};

rc = behavConstsDART; %directories
eval(ds);

%      
day_id = 255; %enter post-DART day
pre_day = expt(day_id).multiday_matchdays;

nd=2; %hardcoding the number of days for now

% INDICATE THE PRE VS. POST DAYS DEPENDING ON THE ORDER OF MATCHING
pre=2;
post=1;

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
nKeep = size(tc_trial_avrg_stat{post},2);
clear d
% find stimulus conditions
frame_rate = input.frameImagingRateMs;
%

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
oris = unique(tOri_match{post});
dirs = unique(tDir_match{post});
cons = unique(tCon_match{post});
nOri = length(oris);
nCon = length(cons);

nOn = input(1).nScansOn;
nOff = input(1).nScansOff;

%% plot trial-by-trial activity in green vs red cell

trialResp=cell(1,2);
green_trialResp=cell(1,2);
red_trialResp=cell(1,2);
linCellProps = nan(6,4);

for id = 1:nd
trialResp{id} = mean(data_trial_keep{id}(stimStart:(stimStart+nOn-1),:,:),1);
green_trialResp{id}=mean(trialResp{id}(:,:,green_ind_keep),3);

red_trialResp{id}=mean(trialResp{id}(:,:,red_ind_keep),3);

end

% figure;
% subplot(2,2,1);
% plot(green_trialResp{pre});
% title('green pre');
% ylabel('trial');
% xlabel('mean dF/F over cells');
% subplot(2,2,2);
% plot(green_trialResp{post});
% title('green post');
% subplot(2,2,3);
% plot(red_trialResp{pre});
% title('red pre');
% subplot(2,2,4);
% plot(red_trialResp{post});
% title('red post');


for iCell = 1:length(red_ind_keep)
    cellID=red_ind_keep(iCell)
    thisCell_pre=mean(trialResp{pre}(:,:,cellID),3); 
    thisCell_post=mean(trialResp{post}(:,:,cellID),3); 
    figure
    scatter(green_trialResp{pre},thisCell_pre, 'MarkerFaceColor','black','MarkerEdgeColor','none','MarkerFaceAlpha', 0.5)
    hold on
    scatter(green_trialResp{pre},thisCell_post,'MarkerFaceColor','blue','MarkerEdgeColor','none','MarkerFaceAlpha', 0.5)
    h = lsline;
    set(h(1),'color','k')
    set(h(1),'color','b')
    
    [R,p]=corrcoef(green_trialResp{pre},thisCell_pre);
    txt2 = ['HT+ ' num2str(length(red_ind_keep))];
    title([num2str(cellID) ' R= ' num2str(R(2))]);
    R_p_values(1,iCell)=R(2);
    R_p_values(2,iCell)=p(2);

end

sig_corr_red = R_p_values(2,:)<0.05;
Rsq_red =  R_p_values(1,:)>.5;
%%
scatter(pref_responses_stat{pre}(red_ind_keep),R_p_values(1,:),'k')
ylabel('Max dF/F') 
xlabel('R^2') 

%
figure;
% subplot(1,2,1)
scatter(green_trialResp{1,pre},red_trialResp{1,pre},10,'MarkerEdgeColor','k')
ylabel('SOM activity')
xlabel('Pyr activity')
%title('stationary')
% ylim([-.05 .15])
% xlim([-.05 .35])
hold on
scatter(green_trialResp{1,post},red_trialResp{1,post},10,'MarkerEdgeColor','b')
hold on


idx = isnan(red_trialResp{1,pre});
linfit = polyfit(green_trialResp{1,pre}(~idx),red_trialResp{1,pre}(~idx),1);
y1 = polyval(linfit,green_trialResp{1,pre});
plot(green_trialResp{1,pre},y1,'k');
[R,p]=corrcoef(green_trialResp{1,pre}(~idx),red_trialResp{1,pre}(~idx)); 
linCellProps(1,1)=linfit(1); %slope
linCellProps(2,1)=linfit(2); %intercept
linCellProps(3,1)=R(2);
linCellProps(4,1)=p(2);
linCellProps(5,1)=min(green_trialResp{1,pre});
linCellProps(6,1)=max(green_trialResp{1,pre});

hold on
idx2 = isnan(red_trialResp{1,post});
linfit = polyfit(green_trialResp{1,post}(~idx2),red_trialResp{1,post}(~idx2),1);
y2 = polyval(linfit,green_trialResp{1,post});
plot(green_trialResp{1,post},y2,'b');
[R,p]=corrcoef(green_trialResp{1,post}(~idx2),red_trialResp{1,post}(~idx2)); 
linCellProps(1,2)=linfit(1); %slope
linCellProps(2,2)=linfit(2); %intercept
linCellProps(3,2)=R(2);
linCellProps(4,2)=p(2);
linCellProps(5,2)=min(green_trialResp{1,post});
linCellProps(6,2)=max(green_trialResp{1,post});
set(gca, 'TickDir', 'out')

% 
% subplot(1,2,2)
% scatter(green_trialResp{2,pre},red_trialResp{2,pre},10,'MarkerEdgeColor','k')
% ylabel('SOM activity')
% xlabel('Pyr activity')
% title('running')
% % ylim([-.05 .15])
% % xlim([-.05 .35])
% hold on
% scatter(green_trialResp{2,post},red_trialResp{2,post},10,'MarkerEdgeColor','b')
% hold on
% 
% idx = isnan(red_trialResp{2,pre});
% linfit = polyfit(green_trialResp{2,pre}(~idx),red_trialResp{2,pre}(~idx),1);
% y1 = polyval(linfit,green_trialResp{2,pre});
% plot(green_trialResp{2,pre},y1,'k');
% [R,p]=corrcoef(green_trialResp{2,pre}(~idx),red_trialResp{2,pre}(~idx)); 
% linCellProps(1,3)=linfit(1); %slope
% linCellProps(2,3)=linfit(2); %intercept
% linCellProps(3,3)=R(2);
% linCellProps(4,3)=p(2);
% linCellProps(5,3)=min(green_trialResp{2,pre});
% linCellProps(6,3)=max(green_trialResp{2,pre});
% 
% 
% hold on
% idx2 = isnan(red_trialResp{2,post});
% linfit1 = polyfit(green_trialResp{2,post}(~idx2),red_trialResp{2,post}(~idx2),1);
% y2 = polyval(linfit1,green_trialResp{2,post});
% plot(green_trialResp{2,post},y2,'b');
% [R,p]=corrcoef(green_trialResp{2,post}(~idx2),red_trialResp{2,post}(~idx2)); 
% linCellProps(1,4)=linfit1(1); %slope
% linCellProps(2,4)=linfit1(2); %intercept
% linCellProps(3,4)=R(2);
% linCellProps(4,4)=p(2);
% linCellProps(5,4)=min(green_trialResp{2,post});
% linCellProps(6,4)=max(green_trialResp{2,post});

sgtitle('Average response for each trial')
x0=5;
y0=5;
width=5;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca, 'TickDir', 'out')

print(fullfile(fn_multi_analysis,[ 'HT_Pyr_relationship.pdf']),'-dpdf');


%% response by condition for cells matched across all conditions
% find cells that I ahve running data for on both days
% haveRunning_pre = ~isnan(pref_responses_loc{pre});
% haveRunning_post = ~isnan(pref_responses_loc{post});
% haveRunning_both = find(haveRunning_pre.* haveRunning_post);
% green_ind_keep = intersect(haveRunning_both, green_ind_keep);
% red_ind_keep = intersect(haveRunning_both, red_ind_keep);


green_ind_keep = green_ind_keep;
red_ind_keep = red_ind_keep;

responseByCond = nan((nCon*2),4);

for iCon = 1:nCon
    if iCon == 1
        counter=1
    else
        counter=counter+2
    end
    
    responseByCond(counter,:)=[mean(pref_responses_stat{pre}(green_ind_keep,iCon), "omitnan") mean(pref_responses_stat{pre}(red_ind_keep,iCon), "omitnan") mean(pref_responses_stat{post}(green_ind_keep,iCon), "omitnan") mean(pref_responses_stat{post}(red_ind_keep,iCon), "omitnan")]
    responseByCond((counter+1),:)=[mean(pref_responses_loc{pre}(green_ind_keep,iCon), "omitnan") mean(pref_responses_loc{pre}(red_ind_keep,iCon), "omitnan") mean(pref_responses_loc{post}(green_ind_keep,iCon), "omitnan") mean(pref_responses_loc{post}(red_ind_keep,iCon), "omitnan")]

end

responseByCondProps = nan(6,2);
figure;
scatter(responseByCond(:,1),responseByCond(:,2),'k')
hold on
linfit = polyfit(responseByCond(:,1),responseByCond(:,2),1);
y1 = polyval(linfit,responseByCond(:,1));
plot(responseByCond(:,1),y1,'k');
[R,p]=corrcoef(responseByCond(:,1),responseByCond(:,2)); 
responseByCondProps(1,1)=linfit(1); %slope
responseByCondProps(2,1)=linfit(2); %intercept
responseByCondProps(3,1)=R(2);
responseByCondProps(4,1)=p(2);
responseByCondProps(5,1)=min(responseByCond(:,1));
responseByCondProps(6,1)=max(responseByCond(:,1));
hold on

scatter(responseByCond(:,3),responseByCond(:,4),'b')
hold on
linfit = polyfit(responseByCond(:,3),responseByCond(:,4),1);
y2 = polyval(linfit,responseByCond(:,3));
plot(responseByCond(:,3),y2,'b');
[R,p]=corrcoef(responseByCond(:,3),responseByCond(:,4)); 
responseByCondProps(1,2)=linfit(1); %slope
responseByCondProps(2,2)=linfit(2); %intercept
responseByCondProps(3,2)=R(2);
responseByCondProps(4,2)=p(2);
responseByCondProps(5,2)=min(responseByCond(:,3));
responseByCondProps(6,2)=max(responseByCond(:,3));

save(fullfile(fn_multi,'HT_pyr_relationship.mat'),'linCellProps','responseByCond','responseByCondProps','sig_corr_red','Rsq_red','R_p_values')


clear R p x0 y0 y1 y2 linfit
%% make figure with se shaded, averaging over contrasts and stationary vs. running

tc_green_avrg = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg = cell(1,nd); %same for red
tc_green_se = cell(1,nd); %this will be the se across all green cells
tc_red_se = cell(1,nd); %same for red

for id = 1:nd

    tc_green_avrg{id}(:)=nanmean(tc_trial_avrg_keep_allCon{id}(:,green_ind_keep),2);
    green_std=std(tc_trial_avrg_keep_allCon{id}(:,green_ind_keep),[],2);
    tc_green_se{id}(:)=green_std/sqrt(length(green_ind_keep));
    
    tc_red_avrg{id}(:)=nanmean(tc_trial_avrg_keep_allCon{id}(:,red_ind_keep),2);
    red_std=std(tc_trial_avrg_keep_allCon{id}(:,red_ind_keep),[],2);
    tc_red_se{id}(:)=red_std/sqrt(length(red_ind_keep));
    
    clear green_std red_std
    
end


%create a time axis in seconds
t=1:(size(tc_green_avrg{1},2));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day

shadedErrorBar(t,tc_red_avrg{pre},tc_red_se{pre},'r');
ylim([-.02 .3]);
hold on
shadedErrorBar(t,tc_green_avrg{pre},tc_green_se{pre});
title(['pre-DART, average'])
txt1 = ['HT- ' num2str(length(green_ind_keep))];
text(-1.5,0.25,txt1);
txt2 = ['HT+ ' num2str(length(red_ind_keep))];
text(-1.5,0.23,txt2,'Color','r');
ylabel('dF/F') 
xlabel('s') 

axis square


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg{post},tc_red_se{post},'r');
ylim([-.02 .3]);
hold on
shadedErrorBar(t,tc_green_avrg{post},tc_green_se{post});
ylabel('dF/F') 
xlabel('s') 
title(['post-DART, average'])
axis square
x0=5;
y0=5;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(fullfile(fn_multi_analysis,[ 'average_stat_timecourses.pdf']),'-dpdf');
end 
clear txt1 txt2

%% make figure with se shaded, one figure per contrast - stationary



tc_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat = cell(1,nd); %same for red
tc_green_se_stat = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat = cell(1,nd); %same for red



for id = 1:nd
    for iCon=1:nCon
        
    tc_green_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat{id}(:,green_ind_keep,iCon),2);
    green_std=nanstd(tc_trial_avrg_stat{id}(:,green_ind_keep,iCon),[],2);
    tc_green_se_stat{id}(:,iCon)=green_std/sqrt(length(green_ind_keep));
    
    tc_red_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat{id}(:,red_ind_keep,iCon),2);
    red_std=nanstd(tc_trial_avrg_stat{id}(:,red_ind_keep,iCon),[],2);
    tc_red_se_stat{id}(:,iCon)=red_std/sqrt(length(red_ind_keep));
    
    clear green_std red_std
    end
end
z=double(nOn)/double(frame_rate)

%creat a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day



%ylim([-.02 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'b');
hold on
% line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
% hold on
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['HT-',' n = ', num2str(length(green_ind_keep))])
ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon),tc_red_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_stat{post}(:,iCon),tc_red_se_stat{post}(:,iCon),'b');
%ylim([-.02 .3]);
hold on
% line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
% hold on
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
%ylabel('dF/F') 
xlabel('s') 
title(['HT+',' n = ', num2str(length(red_ind_keep))])
x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['stationary, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_stat_cellType_timecourses.pdf']),'-dpdf');
end 

clear txt1 



%% make figure with se shaded, one figure per contrast - running

tc_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_loc = cell(1,nd); %same for red
tc_green_se_loc = cell(1,nd); %this will be the se across all green cells
tc_red_se_loc = cell(1,nd); %same for red

for id = 1:nd
    for iCon=1:nCon
    tc_green_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc{id}(:,green_ind_keep,iCon),2);
    green_std=std(tc_trial_avrg_loc{id}(:,green_ind_keep,iCon),[],2);
    tc_green_se_loc{id}(:,iCon)=green_std/sqrt(length(green_ind_keep));
    
    tc_red_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc{id}(:,red_ind_keep,iCon),2);
    red_std=std(tc_trial_avrg_loc{id}(:,red_ind_keep,iCon),[],2);
    tc_red_se_loc{id}(:,iCon)=red_std/sqrt(length(red_ind_keep));
    
    clear green_std red_std
    end
end


%creat a time axis in seconds
t=1:(size(tc_green_avrg_loc{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day

shadedErrorBar(t,tc_red_avrg_loc{pre}(:,iCon),tc_red_se_loc{pre}(:,iCon),'r');
ylim([-.02 .6]);
hold on
shadedErrorBar(t,tc_green_avrg_loc{pre}(:,iCon),tc_green_se_loc{pre}(:,iCon));
title(['Running, pre-DART contrast = ' num2str(cons(iCon))])
txt1 = ['HT- ' num2str(length(green_ind_keep))];
text(-1.5,0.25,txt1);
txt2 = ['HT+ ' num2str(length(red_ind_keep))];
text(-1.5,0.23,txt2,'Color','r');
ylabel('dF/F') 
xlabel('s') 

axis square


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_loc{post}(:,iCon),tc_red_se_loc{post}(:,iCon),'r');
ylim([-.02 .6]);
hold on
shadedErrorBar(t,tc_green_avrg_loc{post}(:,iCon),tc_green_se_loc{post}(:,iCon));
ylabel('dF/F') 
xlabel('s') 
title(['Running, post-DART contrast = ' num2str(cons(iCon))])
axis square
x0=5;
y0=5;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(fullfile(fn_multi_analysis,[num2str(cons(iCon)) '_loc_timecourses.pdf']),'-dpdf');
end 
clear txt1 txt2

%% make a plot of individual timecourses 

setYmin = -.2; %indicate y axes you want
setYmax = 0.85;


for iCon = 1:nCon

figure
subplot(2,2,1)
plot(t, tc_trial_avrg_stat{pre}(:,green_ind_keep,iCon),'color',[0 0 0 .2])
hold on
plot(t, tc_green_avrg_stat{pre}(:,iCon),'color',[0 0 0],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 1')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,3)
plot(t, tc_trial_avrg_stat{pre}(:,red_ind_keep,iCon),'color',[.7 .05 .05 .2])
hold on
plot(t, tc_red_avrg_stat{pre}(:,iCon),'color',[.7 .05 .05],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 1')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,2)
plot(t, tc_trial_avrg_stat{post}(:,green_ind_keep,iCon),'color',[0 0 0 .2])
hold on
plot(t, tc_green_avrg_stat{post}(:,iCon),'color',[0 0 0],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 2')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,4)
plot(t, tc_trial_avrg_stat{post}(:,red_ind_keep,iCon),'color',[.7 .05 .05 .2])
hold on
plot(t, tc_red_avrg_stat{post}(:,iCon),'color',[.7 .05 .05],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 2')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
sgtitle(['Stationary, contrast = '  num2str(cons(iCon))])
hold off

x0=5;
y0=5;
width=6;
height=6;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fn_multi_analysis,[num2str(cons(iCon)) '_stat_indiv_timecourses.pdf']),'-dpdf');
end



for iCon = 1:nCon

figure
subplot(2,2,1)
plot(t, tc_trial_avrg_loc{pre}(:,green_ind_keep,iCon),'color',[0 0 0 .2])
hold on
plot(t, tc_green_avrg_loc{pre}(:,iCon),'color',[0 0 0],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 1')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,3)
plot(t, tc_trial_avrg_loc{pre}(:,red_ind_keep,iCon),'color',[.7 .05 .05 .2])
hold on
plot(t, tc_red_avrg_loc{pre}(:,iCon),'color',[.7 .05 .05],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 1')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,2)
plot(t, tc_trial_avrg_loc{post}(:,green_ind_keep,iCon),'color',[0 0 0 .2])
hold on
plot(t, tc_green_avrg_loc{post}(:,iCon),'color',[0 0 0],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 2')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,4)
plot(t, tc_trial_avrg_loc{post}(:,red_ind_keep,iCon),'color',[.7 .05 .05 .2])
hold on
plot(t, tc_red_avrg_loc{post}(:,iCon),'color',[.7 .05 .05],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 2')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
sgtitle(['Running, contrast = '  num2str(cons(iCon))])
hold off

x0=5;
y0=5;
width=6;
height=6;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fn_multi_analysis,[num2str(cons(iCon)) '_loc_indiv_timecourses.pdf']),'-dpdf');
end

%% ploting pre-and post dart average timecourse seperately for each cell
%this is at all contrasts, stationary

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
        plot(t,tc_trial_avrg_keep_allCon{1,pre}(:,iCell),'r')
        hold on
        plot(t,tc_trial_avrg_keep_allCon{1,post}(:,iCell),'--r')
        hold off
        ylim([setYmin setYmax]);
    else
        plot(t,tc_trial_avrg_keep_allCon{1,pre}(:,iCell),'k')
        hold on
        plot(t,tc_trial_avrg_keep_allCon{1,post}(:,iCell),'--k')
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
                %halfMax = max(smoothData(3:length(smoothData)))/2;
                halfMax = max(smoothData)/2;
                tHalfMaxCell =double(min(find(smoothData>halfMax)))/double(frame_rate);
%                 if rem(iCell, 10) == 0 
%                 figure;plot(tempData)
%                 hold on
%                 plot(smoothData)
%                 hold on
%                 vline(min(find(smoothData>halfMax)))
%                 end
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
scatter((tHalfMax{pre}(green_ind_keep)),(tHalfMax{post}(green_ind_keep)),'k','jitter', 'on', 'jitterAmount', 0.3)
ylabel('post-DART half-max(s)')
xlabel('pre-DART half-max(s)')
% ylim([0 .1])
% xlim([0 .1])
refline(1)
title('HT- ')
axis square
hold off


subplot(1,2,2)
scatter((tHalfMax{pre}(red_ind_keep)),(tHalfMax{post}(red_ind_keep)),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 0.3)
ylabel('post-DART half-max(s)')
xlabel('pre-DART half-max(s)')
% ylim([0 .1])
% xlim([0 .1])
refline(1)
title('HT+ ')
axis square
hold off

sgtitle('time to half max (s), averaged over contrast')
print(fullfile(fn_multi_analysis, ['tHalfMax_crossDay.pdf']),'-dpdf','-bestfit')

%% scatterplot of max df/f for day 1 vs day 2, and each subplot is one cell type


for iCon = 1:nCon
figure; movegui('center') 
subplot(2,2,1)
scatter((pref_responses_stat{pre}(green_ind_keep,iCon)),(pref_responses_stat{post}(green_ind_keep,iCon)),'k')
% hold on
% scatter((resp_max_keep{pre}(intersect(green_ind_keep,find(sig_diff{iCon})),iCon)),(resp_max_keep{post}(intersect(green_ind_keep,find(sig_diff{iCon})),iCon)),'k','MarkerFaceColor','k')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])
refline(1)
title('HT- stationary')
axis square


subplot(2,2,2)
scatter((pref_responses_stat{pre}(red_ind_keep,iCon)),(pref_responses_stat{post}(red_ind_keep,iCon)),'MarkerEdgeColor',[.7 .05 .05])
% hold on
% scatter((resp_max_keep{pre}(intersect(red_ind_keep,find(sig_diff{iCon})),iCon)),(resp_max_keep{post}(intersect(red_ind_keep,find(sig_diff{iCon})),iCon)),'MarkerEdgeColor',[.7 .05 .05],'MarkerFaceColor',[.7 .05 .05])
% hold off
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])

refline(1)
title('HT+ stationary')
axis square


subplot(2,2,3)
scatter((pref_responses_loc{pre}(green_ind_keep,iCon)),(pref_responses_loc{post}(green_ind_keep,iCon)),'k')
% hold on
% scatter((resp_max_keep{pre}(intersect(green_ind_keep,find(sig_diff{iCon})),iCon)),(resp_max_keep{post}(intersect(green_ind_keep,find(sig_diff{iCon})),iCon)),'k','MarkerFaceColor','k')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])
refline(1)
title('HT- running')
axis square

subplot(2,2,4)
scatter((pref_responses_loc{pre}(red_ind_keep,iCon)),(pref_responses_loc{post}(red_ind_keep,iCon)),'MarkerEdgeColor',[.7 .05 .05])
% hold on
% scatter((resp_max_keep{pre}(intersect(red_ind_keep,find(sig_diff{iCon})),iCon)),(resp_max_keep{post}(intersect(red_ind_keep,find(sig_diff{iCon})),iCon)),'MarkerEdgeColor',[.7 .05 .05],'MarkerFaceColor',[.7 .05 .05])
% hold off
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])

refline(1)
title('HT+ running')
axis square

sgtitle(num2str(cons(iCon)))
print(fullfile(fn_multi_analysis,[num2str(cons(iCon)) 'maxResp_crossDay.pdf']),'-dpdf','-bestfit')
end
%% dfof for all stimuli, running and stationary
figure; movegui('center') 
subplot(1,2,1)
scatter((pref_responses_allCon{pre}(green_ind_keep)),(pref_responses_allCon{post}(green_ind_keep)),'k')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])
refline(1)
title('HT- average')
axis square


subplot(1,2,2)
scatter((pref_responses_allCon{pre}(red_ind_keep)),(pref_responses_allCon{post}(red_ind_keep)),'MarkerEdgeColor',[.7 .05 .05])

ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])

refline(1)
title('HT+ average')
axis square

%%
figure;
for iCon = 1:nCon
subplot(1,nCon,iCon)
for iRed = 1:length(red_ind_keep)
    thisCell = red_ind_keep(iRed);
    scatter((resp_max_keep{pre}(thisCell,iCon)),(resp_max_keep{post}(thisCell,iCon)))
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
scatter((resp_max_keep{pre}(red_ind_keep,iCon)),dfof_max_diff_raw(red_ind_keep,iCon),'MarkerEdgeColor',[.7 .05 .05],'MarkerFaceColor',[.7 .05 .05],'MarkerFaceAlpha',.25)
lsline;
hold on
scatter((resp_max_keep{pre}(intersect(red_ind_keep,find(sig_diff{iCon})),iCon)),dfof_max_diff_raw(intersect(red_ind_keep,find(sig_diff{iCon})),iCon),'MarkerFaceColor',[.7 .05 .05])
[R,p]=corrcoef((resp_max_keep{pre}(red_ind_keep,iCon)),dfof_max_diff_raw(red_ind_keep,iCon),'Rows','complete');
title(['R = '  num2str(R(2)) ', p = ' num2str(p(2))])
ylabel('pre-DART - post-DART dF/F')
xlabel('pre-DART  dF/F')
axis square
hold off

subplot(1,2,2)
scatter((resp_max_keep{pre}(red_ind_keep,iCon)),dfof_max_diff(red_ind_keep,iCon),'MarkerEdgeColor',[.7 .05 .05],'MarkerFaceColor',[.7 .05 .05],'MarkerFaceAlpha',.25)
lsline;
hold on
scatter((resp_max_keep{pre}(intersect(red_ind_keep,find(sig_diff{iCon})),iCon)),dfof_max_diff(intersect(red_ind_keep,find(sig_diff{iCon})),iCon),'MarkerFaceColor',[.7 .05 .05])
[R,p]=corrcoef((resp_max_keep{pre}(red_ind_keep,iCon)),dfof_max_diff(red_ind_keep,iCon),'Rows','complete');
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

%% plots stationary and running TCs over one another - not sure I need this any more
% 
% %now I know which trials the mouse was running on, need to look at LMI
% %LMI = (R_loc - R_stat)/(R_loc + R_stat)
% %change this to use the stat and loc timecourses from above
% 
% %creat a time axis in seconds
% t=1:(size(tc_green_avrg_stat{1,1,1},1));
% t=(t-(double(stimStart)-1))/double(frame_rate);
% 
% for iCon = 1:nCon
%     %make this the average for running state for each contrast
%     figure;
%  
%         
%         locr_mean = nanmean(locTCs{pre}(:,red_ind_keep,iCon),2); %averaged across red cells
%         locr_std = std(locTCs{pre}(:,red_ind_keep,iCon),[],2);
%         locr_se=locr_std/sqrt(length(red_ind_keep));
%         statr_mean = nanmean(statTCs{pre}(:,red_ind_keep,iCon),2);
%         statr_std = std(statTCs{pre}(:,red_ind_keep,iCon),[],2);
%         statr_se=statr_std/sqrt(length(red_ind_keep));
% 
%         locg_mean = nanmean(locTCs{pre}(:,green_ind_keep,iCon),2);
%         locg_std = std(locTCs{pre}(:,green_ind_keep,iCon),[],2);
%         locg_se=locg_std/sqrt(length(green_ind_keep));
%         statg_mean = nanmean(statTCs{pre}(:,green_ind_keep,iCon),2);
%         statg_std = std(statTCs{pre}(:,green_ind_keep,iCon),[],2);
%         statg_se=statg_std/sqrt(length(green_ind_keep));
% 
%         subplot(1,2,1)
%         shadedErrorBar(t,locr_mean,locr_se,'r')
%         hold on
%         shadedErrorBar(t,statr_mean,statr_se,'--r')
%         shadedErrorBar(t,locg_mean,locg_se,'k')
%         shadedErrorBar(t,statg_mean,statg_se,'--k')
%         ylim([-.05 .2]);
%         title(['contrast ', num2str(cons(iCon)), ' pre-DART'])
%         txt1=[num2str(locCounts{pre}(1,iCon)), ' running trials']
%         text(-1.5,0.18,txt1);
%         txt2 = [num2str(locCounts{pre}(2,iCon)), ' stationary trials']
%         text(-1.5,0.16,txt2);
%         axis square
%         hold off   
%         
%  
%         
%         locr_mean = nanmean(locTCs{post}(:,red_ind_keep,iCon),2); %averaged across red cells
%         locr_std = std(locTCs{post}(:,red_ind_keep,iCon),[],2);
%         locr_se=locr_std/sqrt(length(red_ind_keep));
%         statr_mean = nanmean(statTCs{post}(:,red_ind_keep,iCon),2);
%         statr_std = std(statTCs{post}(:,red_ind_keep,iCon),[],2);
%         statr_se=statr_std/sqrt(length(red_ind_keep));
% 
%         locg_mean = nanmean(locTCs{post}(:,green_ind_keep,iCon),2);
%         locg_std = std(locTCs{post}(:,green_ind_keep,iCon),[],2);
%         locg_se=locg_std/sqrt(length(green_ind_keep));
% 
%         statg_mean = nanmean(statTCs{post}(:,green_ind_keep,iCon),2);
%         statg_std = std(statTCs{post}(:,green_ind_keep,iCon),[],2);
%         statg_se=statg_std/sqrt(length(green_ind_keep));
%         
%         subplot(1,2,2)
%         shadedErrorBar(t,locr_mean,locr_se,'r')
%         hold on
%         shadedErrorBar(t,statr_mean,statr_se,'--r')
%         shadedErrorBar(t,locg_mean,locg_se,'k')
%         shadedErrorBar(t,statg_mean,statg_se,'--k')
%         ylim([-.05 .2]);
%         title(['contrast ', num2str(cons(iCon)), ' post-DART'])
%         txt1=[num2str(locCounts{post}(1,iCon)), ' running trials']
%         text(-1.5,0.18,txt1);
%         txt2 = [num2str(locCounts{post}(2,iCon)), ' stationary trials']
%         text(-1.5,0.16,txt2);
%         axis square
%         hold off
%         
%        clear locr_mean locr_std locr_se statr_mean statr_std statr_se locg_mean locg_std locg_se statg_mean statg_std statg_se txt1 txt2
%        print(fullfile(fn_multi_analysis,[num2str(cons(iCon)) '_locomotionTCs.pdf']),'-dpdf');
%   
% end

%% locomotion modulation index


for iCon = 1:nCon

figure; movegui('center') 
subplot(1,2,1)
swarmchart((LMI{pre}(green_ind_keep,iCon)),(LMI{post}(green_ind_keep,iCon)),'k')
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
swarmchart((LMI{pre}(red_ind_keep,iCon)),(LMI{post}(red_ind_keep,iCon)),'MarkerEdgeColor',[.7 .05 .05])
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
%print(fullfile(fn_multi_analysis,[num2str(cons(iCon)) '_LMI.pdf']),'-dpdf');
end

%% finding cells that are still saturated by contrast vs. not saturated pre-DART

%identify cells that are still in the rise of the contrast response
%function
risingCells = find(data_con_resp_keep{pre}(:,1)<data_con_resp_keep{pre}(:,2));

% make a plot of individual timecourses for rising cells
setYmin = -.1; %indicate y axes you want
setYmax = 0.8;

for iCon = 1:nCon

figure
subplot(2,2,1)
plot(t, tc_trial_avrg_keep{pre}(:,intersect(risingCells,green_ind_keep),iCon),'k')
ylim([setYmin setYmax]);
title('day 1')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 


subplot(2,2,2)
plot(t, tc_trial_avrg_keep{pre}(:,intersect(risingCells,red_ind_keep),iCon),'color',[.7 .05 .05])
ylim([setYmin setYmax]);
title('day 1')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 


subplot(2,2,3)
plot(t, tc_trial_avrg_keep{post}(:,intersect(risingCells,green_ind_keep),iCon),'k')
ylim([setYmin setYmax]);
title('day 2')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 


subplot(2,2,4)
plot(t, tc_trial_avrg_keep{post}(:,intersect(risingCells,red_ind_keep),iCon),'color',[.7 .05 .05])
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
figure;
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
            errorbar(dirs, data_resp_keep{id}(iCell,:,iCon,1), data_resp_keep{id}(iCell,:,iCon,2),'-o')
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
        errorbar(cons, squeeze(nanmean(data_resp_keep{pre}(iCell,:,:,1),2)), squeeze(nanmean(data_resp_keep{pre}(iCell,:,:,2),2)),'-')
        hold on
        errorbar(cons, squeeze(nanmean(data_resp_keep{post}(iCell,:,:,1),2)), squeeze(nanmean(data_resp_keep{post}(iCell,:,:,2),2)),'-')
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
        errorbar(cons, squeeze(nanmean(data_resp_keep{pre}(iCell,:,:,1),2)), squeeze(nanmean(data_resp_keep{pre}(iCell,:,:,2),2)),'-')
        hold on
        errorbar(cons, squeeze(nanmean(data_resp_keep{post}(iCell,:,:,1),2)), squeeze(nanmean(data_resp_keep{post}(iCell,:,:,2),2)),'--')

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
                scatter(data_resp_keep{pre}(iCell,:,iCon,1),data_resp_keep{post}(iCell,:,iCon,1))
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
            max_list = [max(max(data_resp_keep{pre}(iCell,:,:,1))),max(max(data_resp_keep{post}(iCell,:,:,1)))];
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
                data_pre=[data_pre,data_resp_keep{pre}(iCell,:,iCon,1)]; %collects the values for day 1 over all contrasts
                data_post=[data_post,data_resp_keep{post}(iCell,:,iCon,1)]; %collects the values for day 2 over all contrasts
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
axis square
myRef = refline(1)
myRef.LineStyle = ':'
axis square

title(['\color{red}HT+ cells, all stimulus conditions, n= ' num2str(length(intersect(red_ind_keep,find(linCell))))])
xlabel('pre-DART')
ylabel('post-DART')

myRef = refline(1)
myRef.LineStyle = ':'
axis square
print(fullfile(fn_multi_analysis,'HT+_allResp.pdf'),'-dpdf','-bestfit')


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
                data_pre=[data_pre,data_resp_keep{pre}(iCell,:,iCon,1)]; %collects the values for day 1 over all contrasts
                data_post=[data_post,data_resp_keep{post}(iCell,:,iCon,1)]; %collects the values for day 2 over all contrasts
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
xlim(xRange);
ylim(xRange);
myRef = refline(1)
myRef.LineStyle = ':'
axis square

title(['HT- cells, all stimulus conditions, n= ' num2str(length(intersect(green_ind_keep,find(linCell))))])
xlabel('pre-DART')
ylabel('post-DART')
myRef = refline(1)
myRef.LineStyle = ':'
xlim(xRange);
ylim(xRange);
axis square

print(fullfile(fn_multi_analysis,'HT-_allResp.pdf'),'-dpdf','-bestfit')


%% plot observed pref ori 
%the orientation that gave the strongest response in the observed data

fig=figure; movegui('center') 
scatter(pref_ori_keep{pre}(green_ind_keep),pref_ori_keep{post}(green_ind_keep),'k','jitter', 'on', 'jitterAmount', 5)
hold on
scatter(pref_ori_keep{pre}(red_ind_keep),pref_ori_keep{post}(red_ind_keep),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 5)
hold off
xlabel('Pre-DART pref ori')
ylabel('Post-DART pref ori')
% xlim([0 .4])
% ylim([0 .4])
refline(1)
title('Observed pref ori')
saveas(fig, 'obsvOri.png')

ori_diff_keep = abs(pref_ori_keep{pre}-pref_ori_keep{post});
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
scatter(fit_pref_oris_keep{pre}(green_ind_keep),fit_pref_oris_keep{post}(green_ind_keep),'k','jitter', 'on', 'jitterAmount', 5)
hold on
scatter(fit_pref_oris_keep{pre}(red_ind_keep),fit_pref_oris_keep{post}(red_ind_keep),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 5)
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
[R_g p_g] = corrcoef(fit_pref_oris_keep{pre},fit_pref_oris_keep{post})
figure; movegui('center') 
scatter(fit_pref_oris_keep{pre}((intersect(well_fit{pre},green_ind_keep))),fit_pref_oris_keep{post}((intersect(well_fit{pre},green_ind_keep))),'k','jitter', 'on', 'jitterAmount', 5)
hold on
scatter(fit_pref_oris_keep{pre}(intersect(well_fit{pre},red_ind_keep)),fit_pref_oris_keep{post}(intersect(well_fit{pre},red_ind_keep)),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 5)
hold off
xlabel('D1- pref ori')
ylabel('D2- pref ori')
refline(1)
title('fit pref ori for cells that were well-fit pre-DART')

%%

figure; movegui('center') 
subplot(1,2,1)
scatter(fit_pref_oris_keep{post}(green_ind_keep),pref_ori_keep{post}(green_ind_keep),'k','jitter', 'on', 'jitterAmount', 5)
hold on
scatter(fit_pref_oris_keep{post}(red_ind_keep),pref_ori_keep{post}(red_ind_keep),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 5)
hold off
xlabel('fit pref ori')
ylabel('raw pref ori')
title('day 1')
refline(1)

subplot(1,2,2)
scatter(fit_pref_oris_keep{pre}(green_ind_keep),pref_ori_keep{pre}(green_ind_keep),'k','jitter', 'on', 'jitterAmount', 5)
hold on
scatter(fit_pref_oris_keep{pre}(red_ind_keep),pref_ori_keep{pre}(red_ind_keep),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 5)
hold off
xlabel('fit pref ori')
ylabel('raw pref ori')
title('day 2')
refline(1)
%% example cell tcs - this is to pull out some individual example cell traces
%
cellList=[35,58]; %enter the cells you're interested in by their index wihtin the keep dataframe
iCon=2
place=1;

for i=1:length(cellList)
    figure
    for id = 1:nd
        thisCell = cellList(i)

        %only pulling from dfof data of keep cells
        tCon=tCon_match{id}(1:nTrials(id));
        tDir=tDir_match{id}(1:nTrials(id));
        %identify the trials where ori = pref ori
        temp_ori= pref_dir_keep{id}(thisCell); %find the preferred ori of this cell and convert to degrees
        ori_inds = find(tDir==temp_ori); %these are the trials at that ori
        temp_con = cons(2);%find the preferred contrast of this cell and convert to contrast value
        con_inds=find(tCon==temp_con);
        temp_trials = intersect(ori_inds, con_inds);
        temp_trials(temp_trials==160)=[]
        temp_trials = intersect(temp_trials, find(~RIx{id}))
        temp_TCs=data_trial_keep{id}(:,temp_trials,thisCell); %pulls the selected trials for the selected cell. Shape is frames X trial
        thisCellMean = mean(temp_TCs,2);
        thisCellSE=std(temp_TCs')/sqrt(length(temp_trials));
        
        if id==pre
            shadedErrorBar(t,thisCellMean,thisCellSE,'k')
            ylim([-.2 .5]);
        elseif id==post
            shadedErrorBar(t,thisCellMean,thisCellSE,'b')
            ylim([-.2 .5]);
        end
        hold on
        title(string(thisCell));
        %line([0,2],[.45,.45])
    place=place+1;
    end
    hold off
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
%%
figure;
subplot(2,2,1)
imagesc(corrmap{3})
colormap gray
title('FOV avrg')
hold on
bound = cell2mat(bwboundaries(keep_green_masks(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','b','MarkerSize',2);
caxis([100 7000])
title('HT- pre-DART');
hold off

subplot(2,2,2)
imagesc(corrmap{1})
colormap gray
title('FOV avrg')
hold on
bound = cell2mat(bwboundaries(keep_green_masks(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','b','MarkerSize',2);
caxis([100 7000])
title('HT- post-DART');
hold off



subplot(2,2,3)
imagesc(fov_red{3})
colormap gray
title('FOV avrg')
hold on
bound = cell2mat(bwboundaries(keep_red_masks(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
caxis([5 40])
title('HT+ pre-DART');
hold off

subplot(2,2,4)
imagesc(fov_red{1})
colormap gray
title('FOV avrg')
hold on
bound = cell2mat(bwboundaries(keep_red_masks(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
caxis([10 50])
title('HT+ post-DART');
hold off

%% green and red FOV with masks
sz=size(fov_avg{3});
rgb = zeros(sz(1),sz(2),3);
rgb(:,:,1) = fov_red{3}./max(fov_red{3}(:));
rgb(:,:,2) = fov_avg{3}./max(fov_avg{3}(:));
figure; image(rgb);  movegui('center')
hold on
bound = cell2mat(bwboundaries(keep_masks(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','b','MarkerSize',2);
%% correlating fullTC with full wheel time
%maybe normalize the tc's

wheel_corr = cell(1,nd);


for id = 1:nd
    clean_wheel_speed{id}=wheel_speed{id};
    clean_wheel_speed{id}(find(abs(clean_wheel_speed{id})<4.884))=0;
    clean_wheel_speed{id}=downsample(clean_wheel_speed{id},10);
    clean_fullTC{id}=downsample(fullTC_keep{id},10);
    for iCell = 1:nKeep
        fullTC_keep_norm{id}(:,iCell) = clean_fullTC{id}(:,iCell) - mean(clean_fullTC{id}(:,iCell));
        wheel_corr{id}(iCell)=corr(clean_fullTC{id}(:,iCell),clean_wheel_speed{id}');
    end
end
%%


figure;
subplot(2,1,1)
plot(fullTC_keep_norm{pre}(3000:4000,red_ind_keep))
subplot(2,1,2)
plot(clean_wheel_speed{pre}(3000:4000))
sgtitle('pre-DART HT+')


figure;
subplot(2,1,1)
plot(fullTC_keep_norm{post}(3000:4000,red_ind_keep))
subplot(2,1,2)
plot(clean_wheel_speed{post}(3000:4000))
sgtitle('post-DART HT+')

% figure;
% subplot(2,1,1)
% plot(fullTC_keep_norm{pre}(:,green_ind_keep))
% subplot(2,1,2)
% plot(clean_wheel_speed{pre})
% sgtitle('pre-DART HT-')
% 
% figure;
% subplot(2,1,1)
% plot(fullTC_keep_norm{post}(:,green_ind_keep))
% subplot(2,1,2)
% plot(clean_wheel_speed{post})
% sgtitle('post-DART HT-')
% %,'Color',[0, 0, 0, 0.1]
%%


for iRed = 1:length(red_ind_keep)
    iCell = red_ind_keep(iRed)
    figure;
    scatter(clean_wheel_speed{pre},clean_fullTC{pre}(:,iCell),'k')
    hold on;
    scatter(clean_wheel_speed{post},clean_fullTC{post}(:,iCell),'b')
    ylabel("F")
    xlabel("Running speed")
    title("example HT+")
end
%%
figure;subplot(1,2,1);
scatter(wheel_corr{pre}(green_ind_keep),wheel_corr{post}(green_ind_keep));title('HT-');axis square;
%ylim([-.2 .4]);xlim([-.2 .4]);
ylabel("post-DART");xlabel("pre-DART");
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
subplot(1,2,2);
scatter(wheel_corr{pre}(red_ind_keep),wheel_corr{post}(red_ind_keep));title('HT+');axis square;
%ylim([0 .4]);xlim([0 .4]);
ylabel("post-DART");xlabel("pre-DART");
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
sgtitle("Correlation with wheel speed")