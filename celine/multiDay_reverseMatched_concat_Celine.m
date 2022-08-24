
clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc =  behavConstsDART; %directories
eval(ds);
%136 141 153 161 169 177 183 189

sess_list = [196];%enter all the sessions you want to concatenate
nSess=length(sess_list);

nd=2;%hard coding for two days per experimental session

% INDICATE THE PRE VS. POST DAYS DEPENDING ON THE ORDER OF MATCHING
pre=2;
post=1;

frame_rate = 15;
sess_title = string(sess_list(1));
for iSess = 2:nSess
    sess_title = strcat(sess_title,'_',string(sess_list(iSess)));
end
d=string(datetime('today'));
if nSess == 1
         if expt(sess_list(1)).multiday_timesincedrug_hours>0
            dart_str = [expt(sess_list(1)).drug '_' num2str(expt(sess_list(1)).multiday_timesincedrug_hours) 'Hr'];
        else
            dart_str = 'control';
        end
        
        fnout = fullfile(rc.achAnalysis,expt(sess_list(1)).mouse,['multiday_' dart_str],d);
else
    fnout= fullfile(rc.achAnalysis,strcat('concat', sess_title),d);
end
mkdir(fnout);
cd(fnout)
clear d sess_title



%% concatenating data
mice=[];
tc_trial_avrg_stat_concat=cell(1,nd);
tc_trial_avrg_loc_concat=cell(1,nd);
tc_trial_avrg_keep_allCon_stat_concat=cell(1,nd);
tc_trial_avrg_keep_allCon_loc_concat=cell(1,nd);
resp_keep_concat=cell(1,nd);
resp_max_keep_concat=cell(1,nd);
pref_responses_loc_concat=cell(1,nd);
pref_responses_stat_concat=cell(1,nd);
pref_responses_allCon__stat_concat=cell(2,nd);
pref_responses_allCon__loc_concat=cell(2,nd);
RIx_concat=cell(1,nd);
dirs_concat=[];
cons_concat=[];
green_concat=[];
red_concat=[];
redRsq_concat=[];
dfof_max_diff_concat=[];
nKeep_concat=[];
LMI_concat = cell(1,nd);
data_resp_concat = cell(1,nd);
red_fluor_concat=[];
green_fluor_concat=[];


for iSess = 1:nSess
    day_id = sess_list(iSess)
    mouse = expt(day_id).mouse;
    mice=[mice;mouse];

    if expt(day_id).multiday_timesincedrug_hours>0
        dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fn_multi = fullfile(rc.achAnalysis,mouse,['multiday_' dart_str]);

    load(fullfile(fn_multi,'tc_keep.mat'))
    load(fullfile(fn_multi,'resp_keep.mat'))
    load(fullfile(fn_multi,'input.mat'))
    load(fullfile(fn_multi,'locomotion.mat'))
    load(fullfile(fn_multi,'fluor_intensity.mat'))
    %load(fullfile(fn_multi,'HT_pyr_relationship.mat'))

    nKeep = size(tc_trial_avrg_stat{post},2);


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
    dirs = unique(tDir_match{post});
    cons = unique(tCon_match{post});

    nOn = input(1).nScansOn;
    nOff = input(1).nScansOff;
    %start conatenating
    dirs_concat = [dirs_concat;dirs]; %I will organize thisas rows so I can subsequently make sure everything matches
    cons_concat = [cons_concat;cons];
    red_concat = [red_concat, red_keep_logical];
    green_concat = [green_concat, green_keep_logical];
    nKeep_concat = [nKeep_concat,nKeep];
%   redRsq_concat = [redRsq_concat,Rsq_red];
    clear cons
for id = 1:nd
        tc_trial_avrg_keep_allCon_stat_concat{id} =cat(2,tc_trial_avrg_keep_allCon_stat_concat{id},tc_trial_avrg_keep_allCon_stat{id}(:,:));
        tc_trial_avrg_keep_allCon_loc_concat{id} =cat(2,tc_trial_avrg_keep_allCon_loc_concat{id},tc_trial_avrg_keep_allCon_loc{id}(:,:));
        tc_trial_avrg_stat_concat{id} =cat(2,tc_trial_avrg_stat_concat{id},tc_trial_avrg_stat{id}(:,:,:));
        tc_trial_avrg_loc_concat{id} =cat(2,tc_trial_avrg_loc_concat{id},tc_trial_avrg_loc{id}(:,:,:));
        resp_keep_concat{id}=cat(1,resp_keep_concat{id},resp_keep{id});
        resp_max_keep_concat{id}=cat(1,resp_max_keep_concat{id},resp_max_keep{id}(:,:));
        LMI_concat{id}=cat(1,LMI_concat{id},LMI{id}(:,:));
        pref_responses_loc_concat{id}=cat(1,pref_responses_loc_concat{id},pref_responses_loc{id}(:,:));
        pref_responses_stat_concat{id}=cat(1,pref_responses_stat_concat{id},pref_responses_stat{id}(:,:));
        pref_responses_allCon__stat_concat{1,id}=cat(2,pref_responses_allCon__stat_concat{1,id},pref_responses_allCon_stat{1,id}(:,:));
        pref_responses_allCon__stat_concat{2,id}=cat(2,pref_responses_allCon__stat_concat{2,id},pref_responses_allCon_stat{2,id}(:,:));
        pref_responses_allCon__loc_concat{1,id}=cat(2,pref_responses_allCon__loc_concat{1,id},pref_responses_allCon_loc{1,id}(:,:));
        pref_responses_allCon__loc_concat{2,id}=cat(2,pref_responses_allCon__loc_concat{2,id},pref_responses_allCon_loc{2,id}(:,:));
        RIx_concat{id}=cat(1,RIx_concat{id},sum(RIx{id}));
    end
    dfof_max_diff_concat=cat(1,dfof_max_diff_concat,dfof_max_diff);
    green_fluor_concat=cat(2,green_fluor_concat,green_fluor_keep);
    red_fluor_concat=cat(2,red_fluor_concat,red_fluor_keep);
end

clear mouse day_id nKeep iSess fn_multi cons dirs
clear explanation1 resp_keep tc_trial_avrg_keep_allCon pref_responses_allCon sig_diff pref_con_keep pref_dir_keep tDir_match tOri_match tCon_match data_trial_keep nTrials tc_trial_avrg_keep green_keep_logical red_keep_logical green_ind_keep red_ind_keep
clear LMI RIx locCounts locResp locTCs statResp statTCs wheel_tc
clear data_con_resp_keep data_dir_resp_keep data_rep_keep dfof_max_diff dfof_max_diff_raw explanation2 resp_max_keep data_resp_keep
clear red_fluor_all red_fluor_match green_fluor_match green_fluor_match red_fluor_keep green_fluor_keep
clear tc_trial_avrg_loc tc_trial_avrg_stat tc_trial_avrg_keep_allCond pref_responses_allCond
clear linCellProps responseByCond responseByCondProps sig_corr_red Rsq_red

red_ind_concat = find(red_concat);
green_ind_concat = find(green_concat);
%% check that contrasts and dirs are matched
cons=unique(cons_concat);
nCon = length(cons)
dirs=unique(dirs_concat);
nDir=length(dirs);
nKeep_total = sum(nKeep_concat);
%% find cells that I have running data for on both days
haveRunning_pre = ~isnan(pref_responses_loc_concat{pre});
haveRunning_post = ~isnan(pref_responses_loc_concat{post});
haveRunning_both = find(haveRunning_pre.* haveRunning_post);
haveRunning_green = intersect(haveRunning_both, green_ind_concat);
haveRunning_red = intersect(haveRunning_both, red_ind_concat);


%I use only these cell when plotting locomotion trials, or when directly
%comparing locomotion with stationary

%% make figure with se shaded, one figure per contrast - stationary

tc_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat = cell(1,nd); %same for red
tc_green_se_stat = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat = cell(1,nd); %same for red



for id = 1:nd
    for iCon=1:nCon
        
    tc_green_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,green_ind_concat,iCon),2);
    green_std=nanstd(tc_trial_avrg_stat_concat{id}(:,green_ind_concat,iCon),[],2);
    tc_green_se_stat{id}(:,iCon)=green_std/sqrt(length(green_ind_concat));
    
    tc_red_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,red_ind_concat,iCon),2);
    red_std=nanstd(tc_trial_avrg_stat_concat{id}(:,red_ind_concat,iCon),[],2);
    tc_red_se_stat{id}(:,iCon)=red_std/sqrt(length(red_ind_concat));
    
    clear green_std red_std
    end
end
z=double(nOn)/double(frame_rate);

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
title(['HT-',' n = ', num2str(length(green_ind_concat))])
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
title(['HT+',' n = ', num2str(length(red_ind_concat))])
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
%% shaded error tc for loc

tc_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_loc = cell(1,nd); %same for red
tc_green_se_loc = cell(1,nd); %this will be the se across all green cells
tc_red_se_loc = cell(1,nd); %same for red

for id = 1:nd
    for iCon=1:nCon
    tc_green_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,haveRunning_green,iCon),2);
    green_std=nanstd(tc_trial_avrg_loc_concat{id}(:,haveRunning_green,iCon),[],2);
    tc_green_se_loc{id}(:,iCon)=green_std/sqrt(length(haveRunning_green));
    
    tc_red_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,haveRunning_red,iCon),2);
    red_std=nanstd(tc_trial_avrg_loc_concat{id}(:,haveRunning_red,iCon),[],2);
    tc_red_se_loc{id}(:,iCon)=red_std/sqrt(length(haveRunning_red));
    
    clear green_std red_std
    end
end


%create a time axis in seconds
t=1:(size(tc_green_avrg_loc{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

z=double(nOn)/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day


%ylim([-.02 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_loc{pre}(:,iCon),tc_green_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_loc{post}(:,iCon),tc_green_se_loc{post}(:,iCon),'b');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['HT-',' n = ', num2str(length(haveRunning_green))]);
ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_loc{pre}(:,iCon),tc_red_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_loc{post}(:,iCon),tc_red_se_loc{post}(:,iCon),'b');
%ylim([-.02 .3]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['HT+',' n = ', num2str(length(haveRunning_red))])
x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['running, contrast = ' num2str(cons(iCon))])
print(fullfile(fnout,[num2str(cons(iCon)) '_loc_cellType_timecourses.pdf']),'-dpdf');
end 
clear txt1
%% scatterplot of max df/f for day 1 vs day 2, and each subplot is one cell type
green_ex_list=[]; %to highlight particular cells
red_ex_list=[];

for iCon = 1:nCon
figure; movegui('center') 
subplot(2,2,1)
scatter((pref_responses_stat_concat{pre}(green_ind_concat,iCon)),(pref_responses_stat_concat{post}(green_ind_concat,iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter((pref_responses_stat_concat{pre}(green_ex_list,iCon)),(pref_responses_stat_concat{post}(green_ex_list,iCon)),10,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter(nanmean(pref_responses_stat_concat{pre}(green_ind_concat,iCon)),nanmean(pref_responses_stat_concat{post}(green_ind_concat,iCon)),15,'r')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 1])
xlim([-.1 1])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('HT- stationary')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off


subplot(2,2,2)
scatter((pref_responses_stat_concat{pre}(red_ind_concat,iCon)),(pref_responses_stat_concat{post}(red_ind_concat,iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter((pref_responses_stat_concat{pre}(red_ex_list,iCon)),(pref_responses_stat_concat{post}(red_ex_list,iCon)),10,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter(nanmean(pref_responses_stat_concat{pre}(red_ind_concat,iCon)),nanmean(pref_responses_stat_concat{post}(red_ind_concat)),15,'r')
% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .4])
xlim([-.1 .4])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
title('HT+ stationary')
axis square
hold off

subplot(2,2,3)
scatter((pref_responses_loc_concat{pre}(haveRunning_green,iCon)),(pref_responses_loc_concat{post}(haveRunning_green,iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter(nanmean(pref_responses_loc_concat{pre}(haveRunning_green,iCon)),nanmean(pref_responses_loc_concat{post}(haveRunning_green)),15,'r')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 1.5])
xlim([-.1 1.5])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('HT- running')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off

subplot(2,2,4)
scatter((pref_responses_loc_concat{pre}(haveRunning_red,iCon)),(pref_responses_loc_concat{post}(haveRunning_red,iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter(nanmean(pref_responses_loc_concat{pre}(haveRunning_red,iCon)),nanmean(pref_responses_loc_concat{post}(haveRunning_red,iCon)),15,'r')
% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('HT+ running')
uistack(hline,'bottom');
axis square
hold off
set(gca, 'TickDir', 'out')


sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'maxResp_crossDay.pdf']),'-dpdf','-bestfit')
end
%% make a plot of individual timecourses 

setYmin = -.2; %indicate y axes you want
setYmax = 0.85;
%creat a time axis in seconds



figure
subplot(2,2,1)
plot(t, tc_trial_avrg_keep_allCon_stat_concat{pre}(:,green_ind_concat),'color',[0 0 0 .2])
hold on
plot(t, mean(tc_trial_avrg_keep_allCon_stat_concat{pre}(:,green_ind_concat),2),'color',[0 0 0],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 1')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,3)
plot(t, tc_trial_avrg_keep_allCon_stat_concat{pre}(:,red_ind_concat),'color',[.7 .05 .05 .2])
hold on
plot(t, mean(tc_trial_avrg_keep_allCon_stat_concat{pre}(:,red_ind_concat),2),'color',[.7 .05 .05],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 1')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,2)
plot(t, tc_trial_avrg_keep_allCon_stat_concat{post}(:,green_ind_concat),'color',[0 0 0 .2])
hold on
plot(t, mean(tc_trial_avrg_keep_allCon_stat_concat{post}(:,green_ind_concat),2),'color',[0 0 0],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 2')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,4)
plot(t, tc_trial_avrg_keep_allCon_stat_concat{post}(:,red_ind_concat),'color',[.7 .05 .05 .2])
hold on
plot(t, mean(tc_trial_avrg_keep_allCon_stat_concat{post}(:,red_ind_concat),2),'color',[.7 .05 .05],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 2')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
sgtitle('averaged over contrast')
hold off

x0=5;
y0=5;
width=6;
height=6;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'allCon_indiv_timecourses.pdf'),'-dpdf');

%% time to peak

tHalfMax = cell(1,nd);

for id = 1:nd
    tHalfMaxTemp=nan(nKeep_total,1);

        for iCell=1:nKeep_total
            %pull the data for a given cell at a given contrast (pref dir)
            tempData=tc_trial_avrg_keep_allCon_stat_concat{id}(stimStart:stimStart+nOn,iCell);
            if resp_keep_concat{id}(iCell)
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
%% scatter for tMax




figure; movegui('center') 
subplot(1,2,1)
scatter((tHalfMax{pre}(green_ind_concat)),(tHalfMax{post}(green_ind_concat)),'k','jitter', 'on', 'jitterAmount', 0.3)
hold on
ylabel('post-DART half-max(s)')
xlabel('pre-DART half-max(s)')
% ylim([0 2.5])
% xlim([0 2.5])
refline(1)
title('HT- ')
axis square
hold on
plot(nanmean(tHalfMax{pre}(green_ind_concat)),nanmean(tHalfMax{post}(green_ind_concat)),'b*')
hs = findobj(gca, 'Type','scatter')
NrPoints = numel(hs.XData);
txt1 = ['n = ', num2str(NrPoints)];
text(1,2.2,txt1);
hold off



subplot(1,2,2)
scatter((tHalfMax{pre}(red_ind_concat)),(tHalfMax{post}(red_ind_concat)),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 0.3)
hold on
ylabel('post-DART half-max(s)')
xlabel('pre-DART half-max(s)')
% ylim([0 2.5])
% xlim([0 2.5])
refline(1)
title('HT+')
axis square
hold on
plot(nanmean(tHalfMax{pre}(red_ind_concat)),nanmean(tHalfMax{post}(red_ind_concat)),'b*')
hs = findobj(gca, 'Type','scatter')
NrPoints = numel(hs.XData);
txt1 = ['n = ' num2str(NrPoints)];
text(1,2.2,txt1,'Color','r');
hold off

clear txt1 NrPoints

sgtitle('time to half max (s), averaged over contrast')
print(fullfile(fnout, ['tHalfMax_crossDay.pdf']),'-dpdf','-bestfit')

%% response by condition
a=mean(pref_responses_stat_concat{pre}(green_ind_concat,:), "omitnan");
b=mean(pref_responses_loc_concat{pre}(green_ind_concat,:), "omitnan");

c=mean(pref_responses_stat_concat{pre}(red_ind_concat,:), "omitnan");
d=mean(pref_responses_loc_concat{pre}(red_ind_concat,:), "omitnan");

e=mean(pref_responses_stat_concat{post}(green_ind_concat,:), "omitnan");
f=mean(pref_responses_loc_concat{post}(green_ind_concat,:), "omitnan");

g=mean(pref_responses_stat_concat{post}(red_ind_concat,:), "omitnan");
h=mean(pref_responses_loc_concat{post}(red_ind_concat,:), "omitnan");

responseByCond = horzcat([a';b'],[c';d'],[e';f'],[g';h']);
clear a b c d e f g h

figure;
scatter(responseByCond(:,1),responseByCond(:,2),'k')
hold on
scatter(responseByCond(:,3),responseByCond(:,4),'b')
hold on

h2 = lsline;
hold on
scatter(responseByCond(4:6,1),responseByCond(4:6,2),'MarkerEdgeColor', 'k','MarkerFaceColor', 'k')
hold on
scatter(responseByCond(4:6,3),responseByCond(4:6,4),'MarkerEdgeColor','b','MarkerFaceColor','b')

%% relationship of red fluorescence to dFoF diff
for id = 1:nd
pref_responses_allCon__rect{id}=pref_responses_allCon__stat_concat{1,id};
pref_responses_allCon__rect{id}(pref_responses_allCon__rect{id}<0)=0;
end

resp_diff = (pref_responses_allCon__rect{post} - pref_responses_allCon__rect{pre})./(pref_responses_allCon__rect{post} + pref_responses_allCon__rect{pre});
%resp_diff = (pref_responses_allCon__rect{post} - pref_responses_allCon__rect{pre})./(pref_responses_allCon__rect{pre});
figure;
scatter(red_fluor_concat(red_ind_concat),resp_diff(red_ind_concat),'MarkerEdgeColor', 'r','MarkerFaceColor','r')
hold on
scatter(red_fluor_concat(green_ind_concat),resp_diff(green_ind_concat),'MarkerEdgeColor', 'k','MarkerFaceColor','k')
xlabel('Red fluorescence, norm')
ylabel('[post-pre]/[post+pre]')
sgtitle('stationary, all contrasts')
hold off
print(fullfile(fnout,'norm_red_fluor.pdf'),'-dpdf','-bestfit')


%
figure
scatter(red_fluor_concat(red_ind_concat),pref_responses_allCon__stat_concat{1,pre}(red_ind_concat),'MarkerEdgeColor', 'r','MarkerFaceColor','r')
hold on
scatter(red_fluor_concat(green_ind_concat),pref_responses_allCon__stat_concat{1,pre}(green_ind_concat),'MarkerEdgeColor', 'k','MarkerFaceColor','k')
xlabel('Red fluorescence, norm')
ylabel('max response pre')
sgtitle('stationary, all contrasts')
hold off

figure
scatter(pref_responses_allCon__stat_concat{1,pre}(red_ind_concat),resp_diff(red_ind_concat),'MarkerEdgeColor', 'r','MarkerFaceColor','r')
hold on
scatter(pref_responses_allCon__stat_concat{1,pre}(green_ind_concat),resp_diff(green_ind_concat),'MarkerEdgeColor', 'k','MarkerFaceColor','k')
ylabel('[post-pre]/[post+pre]')
xlabel('max response pre')
sgtitle('stationary, all contrasts')
hold off

%% locomotion modulation index

for iCon = 1:nCon

figure; movegui('center') 
subplot(1,2,1)
scatter((LMI_concat{pre}(green_ind_concat,iCon)),(LMI_concat{post}(green_ind_concat,iCon)),'k','jitter', 'on', 'jitterAmount', 0.3)
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
scatter((LMI_concat{pre}(red_ind_concat,iCon)),(LMI_concat{post}(red_ind_concat,iCon)),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 0.3)
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
print(fullfile(fnout,[num2str(cons(iCon)) '_LMI.pdf']),'-dpdf');
end
%% scatter of stat vs. loc 

for iCon = 1:nCon

figure
subplot(2,2,1)
swarmchart(statResp_concat{pre}(green_ind_concat,iCon),locResp_concat{pre}(green_ind_concat,iCon),'k')
title('Pre-DART')
ylim([-.1 .5]);
xlim([-.1 .5]);
ylabel('Running') 
xlabel('Stationary')
axis square
refline(1)


subplot(2,2,3)
swarmchart(statResp_concat{pre}(red_ind_concat,iCon),locResp_concat{pre}(red_ind_concat,iCon),'MarkerEdgeColor',[.7 .05 .05])
title('Pre-DART')
ylim([-.1 .5]);
xlim([-.1 .5]);
ylabel('Running') 
xlabel('Stationary')
axis square
refline(1)


subplot(2,2,2)
swarmchart(statResp_concat{post}(green_ind_concat,iCon),locResp_concat{post}(green_ind_concat,iCon),'k')
title('Post-DART')
ylim([-.1 .5]);
xlim([-.1 .5]);
ylabel('Running') 
xlabel('Stationary')  
axis square
refline(1)



subplot(2,2,4)
swarmchart(statResp_concat{post}(red_ind_concat,iCon),locResp_concat{post}(red_ind_concat,iCon),'MarkerEdgeColor',[.7 .05 .05])
title('Post-DART')
ylim([-.1 .5]);
xlim([-.1 .5]);
ylabel('Running') 
xlabel('Stationary') 
axis square
refline(1)


sgtitle(['contrast = '  num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_locVSstat.pdf']),'-dpdf');
end


%% plot all responses 
if nKeep_total<36
    [n n2] = subplotn(nKeep_total);
    tot = n.*n2;
else
    n = 5;
    n2 = 5;
    tot = 25;
end

figure;
movegui('center')
start = 1;

for iCell = 1:nKeep_total
     
        if start>tot

            figure; movegui('center')
           
            start = 1;
            
        end
        subplot(n,n2,start)
            for iCon = 1:nCon
                scatter(data_resp_concat{pre}(iCell,:,iCon,1),data_resp_concat{post}(iCell,:,iCon,1))
                %the contrasts are plotted indiviudally so that the points
                %corresponding to different contrasts are different colors,
                hold on
                lsline
            end
                if ismember(iCell,red_ind_concat)
                    title('\color{red}HT+')
                else
                    title('HT-')
                end
            %max_list and upperLim are to find a single axis upper limit
            %that will work for botht eh x and y axes
            max_list = [max(max(data_resp_concat{pre}(iCell,:,:,1))),max(max(data_resp_concat{post}(iCell,:,:,1)))];
            upperLim =max(max_list)+0.1;
            ylim([-.05 upperLim]);
            xlim([-.05 upperLim]);
            axis square
            hold off
            start=start+1;
                
end
%%
linCell=zeros(1,nKeep_total);%will provide a boolean index of which cells have a linear relationship
linCellProps = nan(4,nKeep_total); %this will contain the slopes and intercepts of cells that are deemed linear
data_pre_all=nan(nDir*nCon, nKeep_total);
count=1;
figure;
for iCell = 1:nKeep_total
    if ismember(iCell,red_ind_concat)
 data_pre=[];
 
 data_post=[];
  for iCon = 1:nCon
                data_pre=[data_pre,data_resp_concat{pre}(iCell,:,iCon,1)]; %collects the values for day 1 over all contrasts
                data_post=[data_post,data_resp_concat{post}(iCell,:,iCon,1)]; %collects the values for day 2 over all contrasts
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
        plot(data_pre,f,'color',[0 0 0 .2],'LineWidth',.5)
        hold on
        scatter(data_pre,data_post,'ko')
        alpha(.3)
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
xlim([-.1 .3])
txt1 = ['y=', num2str(round(meanFit(1),2)),'x',num2str(round(meanFit(2),4))];
text(.15,0.01,txt1);
hold off

myRef = refline(1)
myRef.LineStyle = ':'

title(['\color{red}HT+ cells, all stimulus conditions, n= ' num2str(length(intersect(red_ind_concat,find(linCell))))])
xlabel('pre-DART')
ylabel('post-DART')
myRef = refline(1)
myRef.LineStyle = ':'

axis square
print(fullfile(fnout,'HT+_allResp.pdf'),'-dpdf','-bestfit')


%% same as above for HT- cells


linCell=zeros(1,nKeep_total);%will provide a boolean index of which cells have a linear relationship
linCellProps = nan(4,nKeep_total); %this will contain the slopes and intercepts of cells that are deemed linear
data_pre_all=nan(nDir*nCon, nKeep_total);
count=1;
figure;
for iCell = 1:nKeep_total
    if ~ismember(iCell,red_ind_concat)
 data_pre=[];
 
 data_post=[];
  for iCon = 1:nCon
                data_pre=[data_pre,data_resp_concat{pre}(iCell,:,iCon,1)]; %collects the values for day 1 over all contrasts
                data_post=[data_post,data_resp_concat{post}(iCell,:,iCon,1)]; %collects the values for day 2 over all contrasts
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
        plot(data_pre,f,'color',[0 0 0 .2],'LineWidth',.5)
        hold on
        scatter(data_pre,data_post,'ko')
        alpha(.3)
        %FIND A WAY TO CHANGE ALPHA FOR SCATTERS
        
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

title(['HT- cells, all stimulus conditions, n= ' num2str(length(intersect(green_ind_concat,find(linCell))))])
xlabel('pre-DART')
ylabel('post-DART')
myRef = refline(1)
myRef.LineStyle = ':'

print(fullfile(fnout,'HT-_allResp.pdf'),'-dpdf','-bestfit')
%% comparing contrast preference
%need more contrasts 

conTable = table([mean(pref_con_concat{pre}(green_ind_concat));mean(pref_con_concat{pre}(red_ind_concat))],[mean(pref_con_concat{post}(green_ind_concat));mean(pref_con_concat{post}(red_ind_concat))],'VariableNames',{'mean pref con pre' 'mean pref con post'}, 'RowNames',{'HT- cells cells'  'HT+ cells'})
writetable(conTable,fullfile(fnout,'conPref.csv'),'WriteRowNames',true)

boxplot(pref_con_concat{pre},red_concat)
boxplot(pref_con_concat{post},red_concat)

%% making tc plots for low and high Rsq cells

hi_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
low_avrg_stat = cell(1,nd); %same for red
hi_se_stat = cell(1,nd); %this will be the se across all green cells
low_se_stat = cell(1,nd); %same for red



for id = 1:nd
    for iCon=1:nCon
        
    hi_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:, (red_ind_concat(find(redRsq_concat))),iCon),2);
    green_std=nanstd(tc_trial_avrg_stat_concat{id}(:, (red_ind_concat(find(redRsq_concat))),iCon),[],2);
    hi_se_stat{id}(:,iCon)=green_std/sqrt(length( (red_ind_concat(find(redRsq_concat)))));
    
    low_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,(red_ind_concat(find(~redRsq_concat))),iCon),2);
    red_std=nanstd(tc_trial_avrg_stat_concat{id}(:,(red_ind_concat(find(~redRsq_concat))),iCon),[],2);
    low_se_stat{id}(:,iCon)=red_std/sqrt(length((red_ind_concat(find(~redRsq_concat)))));
    
    clear green_std red_std
    end
end
z=double(nOn)/double(frame_rate);

%creat a time axis in seconds
t=1:(size(hi_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for HT-



ylim([-.02 .2]);
hold on
shadedErrorBar(t,hi_avrg_stat{pre}(:,iCon),hi_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,hi_avrg_stat{post}(:,iCon),hi_se_stat{post}(:,iCon),'b');
hold on
% line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
% hold on
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['HT+ high R^2',' n = ', num2str(length((red_ind_concat(find(redRsq_concat)))))])
ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for HT+
shadedErrorBar(t,low_avrg_stat{pre}(:,iCon),low_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,low_avrg_stat{post}(:,iCon),low_se_stat{post}(:,iCon),'b');
ylim([-.02 .2]);
hold on
% line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
% hold on
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
%ylabel('dF/F') 
xlabel('s') 
title(['HT+ low R^2',' n = ', num2str(length((red_ind_concat(find(~redRsq_concat)))))])
x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['stationary, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_Rsq_timecourses.pdf']),'-dpdf');
end 

clear txt1 
%% scatterplot of max df/f for day 1 vs day 2 for high and low Rsq red cells

for iCon = 1:nCon
figure; movegui('center') 
subplot(1,2,1)
scatter((pref_responses_stat_concat{pre}((red_ind_concat(find(redRsq_concat))),iCon)),(pref_responses_stat_concat{post}((red_ind_concat(find(redRsq_concat))),iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter(nanmean(pref_responses_stat_concat{pre}((red_ind_concat(find(redRsq_concat))),iCon)),nanmean(pref_responses_stat_concat{post}((red_ind_concat(find(redRsq_concat))),iCon)),15,'r*')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])
refline(1)
title('HT+ high R^2')
axis square
set(gca, 'TickDir', 'out')



subplot(1,2,2)
scatter((pref_responses_stat_concat{pre}((red_ind_concat(find(~redRsq_concat))),iCon)),(pref_responses_stat_concat{post}((red_ind_concat(find(~redRsq_concat))),iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter(nanmean(pref_responses_stat_concat{pre}((red_ind_concat(find(~redRsq_concat))),iCon)),nanmean(pref_responses_stat_concat{post}((red_ind_concat(find(~redRsq_concat))))),15,'r*')
% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])
set(gca, 'TickDir', 'out')
refline(1)
title('HT+ low R^2')
axis square




sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'maxResp_Rsq.pdf']),'-dpdf','-bestfit')
end

%% to plot timecourses averaged over contrast

tc_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat = cell(1,nd); %same for red
tc_green_se_stat = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat = cell(1,nd); %same for red
iCon=1;


for id = 1:nd
   
        
    tc_green_avrg_stat{id}(:)=nanmean(tc_trial_avrg_keep_allCon_stat_concat{id}(:,green_ind_concat),2);
    green_std=nanstd(tc_trial_avrg_keep_allCon_stat_concat{id}(:,green_ind_concat),[],2);
    tc_green_se_stat{id}(:)=green_std/sqrt(length(green_ind_concat));
    
    tc_red_avrg_stat{id}(:)=nanmean(tc_trial_avrg_keep_allCon_stat_concat{id}(:,red_ind_concat),2);
    red_std=nanstd(tc_trial_avrg_keep_allCon_stat_concat{id}(:,red_ind_concat),[],2);
    tc_red_se_stat{id}(:)=red_std/sqrt(length(red_ind_concat));
    
    clear green_std red_std
    
end
z=double(nOn)/double(frame_rate)

%creat a time axis in seconds
t=1:(size(tc_green_avrg_stat{1},2));
t=(t-(double(stimStart)-1))/double(frame_rate);


figure
subplot(1,2,1) %

%ylim([-.02 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:),tc_green_se_stat{pre}(:),'k');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:),tc_green_se_stat{post}(:),'b');
hold on
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['HT-',' n = ', num2str(length(green_ind_concat))])
ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_stat{pre}(:),tc_red_se_stat{pre}(:),'k');
hold on
shadedErrorBar(t,tc_red_avrg_stat{post}(:),tc_red_se_stat{post}(:),'b');
%ylim([-.02 .3]);
hold on
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
%ylabel('dF/F') 
xlabel('s') 
title(['HT+',' n = ', num2str(length(red_ind_concat))])
x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['stationary, all contrasts'])

print(fullfile(fnout,[num2str(cons(iCon)) '_stat_timecourses_allCon.pdf']),'-dpdf');


clear txt1 
%% all contrast timecourse running


tc_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_loc = cell(1,nd); %same for red
tc_green_se_loc = cell(1,nd); %this will be the se across all green cells
tc_red_se_loc = cell(1,nd); %same for red
iCon=1;


for id = 1:nd
   
        
    tc_green_avrg_loc{id}(:)=nanmean(tc_trial_avrg_keep_allCon_loc_concat{id}(:,haveRunning_green),2);
    green_std=nanstd(tc_trial_avrg_keep_allCon_loc_concat{id}(:,haveRunning_green),[],2);
    tc_green_se_loc{id}(:)=green_std/sqrt(length(haveRunning_green));
    
    tc_red_avrg_loc{id}(:)=nanmean(tc_trial_avrg_keep_allCon_loc_concat{id}(:,haveRunning_red),2);
    red_std=nanstd(tc_trial_avrg_keep_allCon_loc_concat{id}(:,haveRunning_red),[],2);
    tc_red_se_loc{id}(:)=red_std/sqrt(length(haveRunning_red));
    
    clear green_std red_std
    
end
z=double(nOn)/double(frame_rate)

%creat a time axis in seconds
t=1:(size(tc_green_avrg_loc{1},2));
t=(t-(double(stimStart)-1))/double(frame_rate);


figure
subplot(1,2,1) %

%ylim([-.02 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_loc{pre}(:),tc_green_se_loc{pre}(:),'k');
hold on
shadedErrorBar(t,tc_green_avrg_loc{post}(:),tc_green_se_loc{post}(:),'b');
hold on
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['HT-',' n = ', num2str(length(haveRunning_green))])
ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_loc{pre}(:),tc_red_se_loc{pre}(:),'k');
hold on
shadedErrorBar(t,tc_red_avrg_loc{post}(:),tc_red_se_loc{post}(:),'b');
%ylim([-.02 .3]);
hold on
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
%ylabel('dF/F') 
xlabel('s') 
title(['HT+',' n = ', num2str(length(haveRunning_red))])
x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['Running, all contrasts'])

print(fullfile(fnout,[num2str(cons(iCon)) '_loc_timecourses_allCon.pdf']),'-dpdf');


clear txt1 
%% all contrast scatter
green_ex_list=[]; %to highlight particular cells
red_ex_list=[];


figure; movegui('center') 
subplot(2,2,1)
scatter((tc_trial_avrg_keep_allCon_stat_concat{pre}(green_ind_concat)),(tc_trial_avrg_keep_allCon_stat_concat{post}(green_ind_concat)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter((tc_trial_avrg_keep_allCon_stat_concat{pre}(green_ex_list)),(tc_trial_avrg_keep_allCon_stat_concat{post}(green_ex_list)),10,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter(nanmean(tc_trial_avrg_keep_allCon_stat_concat{pre}(green_ind_concat)),nanmean(tc_trial_avrg_keep_allCon_stat_concat{post}(green_ind_concat)),15,'k','MarkerFaceColor',[0 0 0])
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
% ylim([-.2 .6])
% xlim([-.2 .6])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('HT- stationary')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off


subplot(2,2,2)
scatter((tc_trial_avrg_keep_allCon_stat_concat{pre}(red_ind_concat)),(tc_trial_avrg_keep_allCon_stat_concat{post}(red_ind_concat)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter((tc_trial_avrg_keep_allCon_stat_concat{pre}(red_ex_list)),(tc_trial_avrg_keep_allCon_stat_concat{post}(red_ex_list)),10,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter(nanmean(tc_trial_avrg_keep_allCon_stat_concat{pre}(red_ind_concat)),nanmean(tc_trial_avrg_keep_allCon_stat_concat{post}(red_ind_concat)),15,'k','MarkerFaceColor',[0 0 0])
% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.2 .2])
xlim([-.2 .2])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
title('HT+ stationary')
axis square
hold off

subplot(2,2,3)
scatter((tc_trial_avrg_keep_allCon_loc_concat{pre}(haveRunning_green)),(tc_trial_avrg_keep_allCon_loc_concat{post}(haveRunning_green)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter(nanmean(tc_trial_avrg_keep_allCon_loc_concat{pre}(haveRunning_green)),nanmean(tc_trial_avrg_keep_allCon_loc_concat{post}(haveRunning_green)),15,'k','MarkerFaceColor',[0 0 0])
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.2 .6])
xlim([-.2 .6])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('HT- running')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off

subplot(2,2,4)
scatter((tc_trial_avrg_keep_allCon_loc_concat{pre}(haveRunning_red)),(tc_trial_avrg_keep_allCon_loc_concat{post}(haveRunning_red)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter(nanmean(tc_trial_avrg_keep_allCon_loc_concat{pre}(haveRunning_red)),nanmean(tc_trial_avrg_keep_allCon_loc_concat{post}(haveRunning_red)),15,'k','MarkerFaceColor',[0 0 0])
% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.2 .2])
xlim([-.2 .2])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('HT+ running')
uistack(hline,'bottom');
axis square
hold off
set(gca, 'TickDir', 'out')


print(fullfile(fnout,['allCon_maxResp_crossDay.pdf']),'-dpdf','-bestfit')
%% to choose example cells
cellList=[152]; %enter the cells you're interested in by their index wihtin the keep dataframe
iCon = 2;

place=1;
figure
for i=1:length(cellList)
    for id = 1:nd
        thisCell = cellList(i)
        
        temp_TCs=tc_trial_avrg_stat_concat{id}(:,thisCell,iCon); %pulls the selected trials for the selected cell. Shape is frames X trial
       
        subplot(length(cellList),nd,place)
        if red_concat(thisCell)==1
            plot(temp_TCs,'r')
            ylim([-.2 .5]);
        else
            plot(temp_TCs,'k')
            ylim([-.2 .5]);
        end
        title(string(thisCell));
        line([0,2],[.45,.45])
    place=place+1;
    end
    
end


c = linspace(1,10,length(cellList));
figure
scatter((pref_responses_stat_concat{pre}(cellList,iCon)),(pref_responses_stat_concat{post}(cellList,iCon)),[],c,'filled')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])
colorbar

refline(1)
axis square

