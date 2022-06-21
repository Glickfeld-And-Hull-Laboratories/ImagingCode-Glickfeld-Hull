
clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc =  behavConstsDART; %directories
eval(ds);
%131 133 138 142 163 171
sess_list = [131 133 138 142 163 171];%enter all the sessions you want to concatenate
nSess=length(sess_list);

nd=2;%hard coding for two days per experimental session

% INDICATE THE PRE VS. POST DAYS DEPENDING ON THE ORDER OF MATCHING
pre=2;
post=1;

targetCon = .5%what contrast to extract for all data - must be one that all datasets had

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
tc_trial_avrg_keep_allCond_concat=cell(1,nd);
resp_keep_concat=cell(1,nd);
resp_max_keep_concat=cell(1,nd);
pref_responses_loc_concat=cell(1,nd);
pref_responses_stat_concat=cell(1,nd);
RIx_concat=cell(1,nd);
oris_concat=[];
cons_concat=[];
green_concat=[];
red_concat=[];
dfof_max_diff_concat=[];
nKeep_concat=[];
LMI_concat = cell(1,nd);
data_resp_concat = cell(1,nd);

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
    oris = unique(tOri_match{post});
    cons = unique(tCon_match{post});
    sharedCon=find(cons==targetCon);

    nOn = input(1).nScansOn;
    nOff = input(1).nScansOff;
    %start conatenating
    oris_concat = [oris_concat,oris]; %I will organize thisas rows so I can subsequently make sure everything matches
    cons_concat = [cons_concat,cons(sharedCon)];
    red_concat = [red_concat, red_keep_logical];
    green_concat = [green_concat, green_keep_logical];
    nKeep_concat = [nKeep_concat,nKeep];
%     statCounts = [statCounts,trialCounts.preStat(1);trialCounts.postStat(1)];
%     locCounts = [locCounts,trialCounts.preLoc(1);trialCounts.postLoc(1)]
    clear cons
    
    for id = 1:nd
        tc_trial_avrg_keep_allCond_concat{id} =cat(2,tc_trial_avrg_keep_allCond_concat{id},tc_trial_avrg_keep_allCond{id}(:,:));
        tc_trial_avrg_stat_concat{id} =cat(2,tc_trial_avrg_stat_concat{id},tc_trial_avrg_stat{id}(:,:,sharedCon));
        tc_trial_avrg_loc_concat{id} =cat(2,tc_trial_avrg_loc_concat{id},tc_trial_avrg_loc{id}(:,:,sharedCon));
        resp_keep_concat{id}=cat(1,resp_keep_concat{id},resp_keep{id});
        resp_max_keep_concat{id}=cat(1,resp_max_keep_concat{id},resp_max_keep{id}(:,sharedCon));
        LMI_concat{id}=cat(1,LMI_concat{id},LMI{id}(:,sharedCon));
        pref_responses_loc_concat{id}=cat(1,pref_responses_loc_concat{id},pref_responses_loc{id}(:,sharedCon));
        pref_responses_stat_concat{id}=cat(1,pref_responses_stat_concat{id},pref_responses_stat{id}(:,sharedCon));
        RIx_concat{id}=cat(1,RIx_concat{id},sum(RIx{id}));
    end
    dfof_max_diff_concat=cat(1,dfof_max_diff_concat,dfof_max_diff(:,sharedCon));
end
%
clear mouse day_id nKeep iSess fn_multi cons oris
clear explanation1 resp_keep tc_trial_avrg_keep_allCond pref_responses_allCond sig_diff pref_con_keep pref_ori_keep tOri_match tCon_match data_trial_keep nTrials tc_trial_avrg_keep green_keep_logical red_keep_logical green_ind_keep red_ind_keep
clear LMI RIx locCounts locResp locTCs statResp statTCs wheel_tc
clear data_con_resp_keep data_ori_resp_keep data_rep_keep dfof_max_diff dfof_max_diff_raw explanation2 resp_max_keep data_resp_keep pref_responses_stat pref_responses_loc
clear tc_trial_avrg_stat tc_trial_avrg_loc
red_ind_concat = find(red_concat);
green_ind_concat = find(green_concat);
%
cons=unique(cons_concat);
nCon = length(cons)
oris=unique(oris_concat);
nOri=length(oris);
nKeep_total = sum(nKeep_concat);
mean(RIx_concat{pre})
mean(RIx_concat{post})
%% find cells that I ahve running data for on both days
haveRunning_pre = ~isnan(pref_responses_loc_concat{pre});
haveRunning_post = ~isnan(pref_responses_loc_concat{post});
haveRunning_both = find(haveRunning_pre.* haveRunning_post);
haveRunning_green = intersect(haveRunning_both, green_ind_concat);
haveRunning_red = intersect(haveRunning_both, red_ind_concat);


%% make figure with se shaded, averaging over stationary vs. running

tc_green_avrg = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg = cell(1,nd); %same for red
tc_green_se = cell(1,nd); %this will be the se across all green cells
tc_red_se = cell(1,nd); %same for red

for id = 1:nd

    tc_green_avrg{id}(:)=nanmean(tc_trial_avrg_keep_allCond_concat{id}(:,haveRunning_green),2);
    green_std=std(tc_trial_avrg_keep_allCond_concat{id}(:,haveRunning_green),[],2);
    tc_green_se{id}(:)=green_std/sqrt(length(haveRunning_green));
    
    tc_red_avrg{id}(:)=nanmean(tc_trial_avrg_keep_allCond_concat{id}(:,haveRunning_red),2);
    red_std=std(tc_trial_avrg_keep_allCond_concat{id}(:,haveRunning_red),[],2);
    tc_red_se{id}(:)=red_std/sqrt(length(haveRunning_red));
    
    clear green_std red_std
    
end


%creat a time axis in seconds
t=1:(size(tc_green_avrg{1},2));
t=(t-(double(stimStart)-1))/double(frame_rate);


figure
subplot(1,2,1) %for the first day

shadedErrorBar(t,tc_red_avrg{pre},tc_red_se{pre},'r');
ylim([-.02 .3]);
hold on
shadedErrorBar(t,tc_green_avrg{pre},tc_green_se{pre});
hold on
line([0,2],[-.01,-.01],'Color','black','LineWidth',2);
title(['pre-DART, average'])
txt1 = ['HT- ' num2str(length(green_ind_concat))];
text(-1.5,0.25,txt1);
txt2 = ['HT+ ' num2str(length(red_ind_concat))];
text(-1.5,0.23,txt2,'Color','r');
ylabel('dF/F') 
xlabel('s') 




subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg{post},tc_red_se{post},'r');
ylim([-.02 .3]);
hold on
shadedErrorBar(t,tc_green_avrg{post},tc_green_se{post});
hold on
line([0,2],[-.01,-.01],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['post-DART, average'])

x0=5;
y0=5;
width=4;
height=3;
sgtitle('average all trials')
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca, 'TickDir', 'out')
print(fullfile(fnout,[num2str(targetCon) 'average_timecourses.pdf']),'-dpdf');

clear txt1 txt2


%% make figure with se shaded, one figure per contrast - stationary

tc_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat = cell(1,nd); %same for red
tc_green_se_stat = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat = cell(1,nd); %same for red



for id = 1:nd
    for iCon=1:nCon
        
    tc_green_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,haveRunning_green,iCon),2);
    green_std=nanstd(tc_trial_avrg_stat_concat{id}(:,haveRunning_green,iCon),[],2);
    tc_green_se_stat{id}(:,iCon)=green_std/sqrt(length(haveRunning_green));
    
    tc_red_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,haveRunning_red,iCon),2);
    red_std=nanstd(tc_trial_avrg_stat_concat{id}(:,haveRunning_red,iCon),[],2);
    tc_red_se_stat{id}(:,iCon)=red_std/sqrt(length(haveRunning_red));
    
    clear green_std red_std
    end
end


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
line([0,2],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title('HT-')
txt1 = ['n = ', num2str(length(haveRunning_green))];
%txt2 = ['post ',num2str(length(find(~isnan(pref_responses_stat_concat{post}(green_ind_concat)))))];
text(-1.5,-0.03,txt1);
%text(0.75,-0.03,txt2,'Color','b');
ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon),tc_red_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_stat{post}(:,iCon),tc_red_se_stat{post}(:,iCon),'b');
%ylim([-.02 .3]);
hold on
line([0,2],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title('HT+')
txt1 = ['n = ', num2str(length(haveRunning_red))];
%txt2 = ['post ',num2str(length(find(~isnan(pref_responses_stat_concat{post}(red_ind_concat)))))];
text(-1.5,-0.03,txt1);
%text(0.75,-0.03,txt2,'Color','b');
x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['stationary, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_stat_cellType_timecourses.pdf']),'-dpdf');
end 
% 
% 
% for iCon = 1:nCon
% figure
% subplot(1,2,1) %for the first day
% 
% 
% ylim([-.02 .3]);
% hold on
% shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'k');
% hold on
% shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon),tc_red_se_stat{pre}(:,iCon),'r');
% hold on
% line([0,2],[-.01,-.01],'Color','black','LineWidth',2);
% txt1 = ['HT- ' num2str(length(green_ind_concat))];
% text(-1.5,0.25,txt1);
% txt2 = ['HT+ ' num2str(length(red_ind_concat))];
% text(-1.5,0.23,txt2,'Color','r');
% title('Post-DART')
% set(gca,'visible','off')
% 
% ylabel('dF/F') 
% xlabel('s') 
% box off
% 
% 
% subplot(1,2,2) %for the second day
% shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'k');
% hold on
% shadedErrorBar(t,tc_red_avrg_stat{post}(:,iCon),tc_red_se_stat{post}(:,iCon),'r');
% ylim([-.02 .3]);
% hold on
% line([0,2],[-.01,-.01],'Color','black','LineWidth',2);
% ylabel('dF/F') 
% xlabel('s') 
% box off
% title('Post-DART')
% sgtitle(['stationary, contrast = ' num2str(cons(iCon))])
% x0=5;
% y0=5;
% width=4;
% height=3;
% set(gcf,'units','inches','position',[x0,y0,width,height])
% 
% set(gca,'visible','off')
% print(fullfile(fnout,[num2str(cons(iCon)) '_stat_timecourses.pdf']),'-dpdf');
% end 


clear txt1 txt2

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


%creat a time axis in seconds
t=1:(size(tc_green_avrg_loc{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);
% 
% 
% for iCon = 1:nCon
% figure
% subplot(1,2,1) %for the first day
% 
% 
% ylim([-.02 .3]);
% hold on
% shadedErrorBar(t,tc_green_avrg_loc{pre}(:,iCon),tc_green_se_loc{pre}(:,iCon),'k');
% hold on
% shadedErrorBar(t,tc_red_avrg_loc{pre}(:,iCon),tc_red_se_loc{pre}(:,iCon),'r');
% hold on
% line([0,2],[-.01,-.01],'Color','black','LineWidth',2);
% txt1 = ['HT- ' num2str(length(green_ind_concat))];
% text(-1.5,0.25,txt1);
% txt2 = ['HT+ ' num2str(length(red_ind_concat))];
% text(-1.5,0.23,txt2,'Color','r');
% title('Post-DART')
% 
% 
% ylabel('dF/F') 
% xlabel('s') 
% box off
% 
% 
% subplot(1,2,2) %for the second day
% shadedErrorBar(t,tc_green_avrg_loc{post}(:,iCon),tc_green_se_loc{post}(:,iCon),'k');
% hold on
% shadedErrorBar(t,tc_red_avrg_loc{post}(:,iCon),tc_red_se_loc{post}(:,iCon),'r');
% ylim([-.02 .3]);
% hold on
% line([0,2],[-.01,-.01],'Color','black','LineWidth',2);
% ylabel('dF/F') 
% xlabel('s') 
% box off
% title('Post-DART')
% sgtitle(['running, contrast = ' num2str(cons(iCon))])
% x0=5;
% y0=5;
% width=4;
% height=3;
% set(gcf,'units','inches','position',[x0,y0,width,height])
% set(gca, 'TickDir', 'out')
% 
% print(fullfile(fnout,[num2str(cons(iCon)) '_loc_timecourses.pdf']),'-dpdf');
% end 

for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day


%ylim([-.02 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_loc{pre}(:,iCon),tc_green_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_loc{post}(:,iCon),tc_green_se_loc{post}(:,iCon),'b');
hold on
line([0,2],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title('HT-')
txt1 = ['n = ', num2str(length(haveRunning_green))];
%txt2 = ['post ',num2str(length(find(~isnan(pref_responses_loc_concat{post}(green_ind_concat)))))];
text(-1.5,-0.03,txt1);
%text(0.75,-0.03,txt2,'Color','b');
ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_loc{pre}(:,iCon),tc_red_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_loc{post}(:,iCon),tc_red_se_loc{post}(:,iCon),'b');
%ylim([-.02 .3]);
hold on
line([0,2],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title('HT+')
txt1 = ['n = ',num2str(length(haveRunning_red))];
%txt2 = ['post ', num2str()];
text(-1.5,-0.03,txt1);
%text(0.75,-0.03,txt2,'Color','b');
x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['running, contrast = ' num2str(cons(iCon))])
print(fullfile(fnout,[num2str(cons(iCon)) '_loc_cellType_timecourses.pdf']),'-dpdf');
end 
clear txt1 txt2

%% scatterplot of max df/f for day 1 vs day 2, and each subplot is one cell type
green_ex_list=[]; %to highlight particular cells
red_ex_list=[];

for iCon = 1:nCon
figure; movegui('center') 
subplot(2,2,1)
scatter((pref_responses_stat_concat{pre}(haveRunning_green,iCon)),(pref_responses_stat_concat{post}(haveRunning_green,iCon)),10,'MarkerEdgeColor',[.7 .7 .7],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter((pref_responses_stat_concat{pre}(green_ex_list,iCon)),(pref_responses_stat_concat{post}(green_ex_list,iCon)),10,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter(nanmean(pref_responses_stat_concat{pre}(haveRunning_green,iCon)),nanmean(pref_responses_stat_concat{post}(haveRunning_green,iCon)),15,'r*')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 1])
xlim([-.1 1])
refline(1)
title('HT- stationary')
axis square
set(gca, 'TickDir', 'out')



subplot(2,2,2)
scatter((pref_responses_stat_concat{pre}(haveRunning_red,iCon)),(pref_responses_stat_concat{post}(haveRunning_red,iCon)),10,'MarkerEdgeColor',[.7 .7 .7],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter((pref_responses_stat_concat{pre}(red_ex_list,iCon)),(pref_responses_stat_concat{post}(red_ex_list,iCon)),10,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter(nanmean(pref_responses_stat_concat{pre}(haveRunning_red,iCon)),nanmean(pref_responses_stat_concat{post}(haveRunning_red,iCon)),15,'r*')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 1])
xlim([-.1 1])
set(gca, 'TickDir', 'out')
refline(1)
title('HT+ stationary')
axis square


subplot(2,2,3)
scatter((pref_responses_loc_concat{pre}(haveRunning_green,iCon)),(pref_responses_loc_concat{post}(haveRunning_green,iCon)),10,'MarkerEdgeColor',[.7 .7 .7],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter(nanmean(pref_responses_loc_concat{pre}(haveRunning_green,iCon)),nanmean(pref_responses_loc_concat{post}(haveRunning_green,iCon)),15,'r*')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 1])
xlim([-.1 1])
refline(1)
title('HT- running')
axis square
set(gca, 'TickDir', 'out')

subplot(2,2,4)
scatter((pref_responses_loc_concat{pre}(haveRunning_red,iCon)),(pref_responses_loc_concat{post}(haveRunning_red,iCon)),10,'MarkerEdgeColor',[.7 .7 .7],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter(nanmean(pref_responses_loc_concat{pre}(haveRunning_red,iCon)),nanmean(pref_responses_loc_concat{post}(haveRunning_red,iCon)),15,'r*')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 1])
xlim([-.1 1])

refline(1)
title('HT+ running')
axis square

set(gca, 'TickDir', 'out')

sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'maxResp_crossDay.pdf']),'-dpdf','-bestfit')
end

%%
responseTable = table([nanmean(pref_responses_stat_concat{pre}(haveRunning_green));nanmean(pref_responses_stat_concat{post}(haveRunning_green))],[nanmean(pref_responses_stat_concat{pre}(haveRunning_red));nanmean(pref_responses_stat_concat{post}(haveRunning_red))],[nanmean(pref_responses_loc_concat{pre}(haveRunning_green));nanmean(pref_responses_loc_concat{post}(haveRunning_green))],[nanmean(pref_responses_loc_concat{pre}(haveRunning_red));nanmean(pref_responses_loc_concat{post}(haveRunning_red))],'VariableNames',{'Pyramidal cells stat'  'HT+ cells stat' 'Pyramidal cells loc'  'HT+ cells loc'}, 'RowNames',{'Pre'  'Post'})
writetable(responseTable,fullfile(fnout,[num2str(targetCon) 'responseTable.csv']),'WriteRowNames',true)

    %% 
pref_responses_stat_transform = cell(1,nd);
for id = 1:nd
    
   pref_responses_stat_transform{id} = pref_responses_stat_concat{id};
   pref_responses_stat_transform{id}(find(pref_responses_stat_transform{id}<0))=0;
   pref_responses_stat_transform{id}=log(pref_responses_stat_transform{id});
end
figure
scatter((pref_responses_stat_transform{pre}(red_ind_concat,iCon)),(pref_responses_stat_transform{post}(red_ind_concat,iCon)),10,'MarkerEdgeColor',[.4 .4 .4],'jitter', 'on', 'jitterAmount',.01)
refline(1);
hold on
scatter(nanmean(pref_responses_stat_transform{pre}(red_ind_concat,iCon)),nanmean(pref_responses_stat_transform{post}(red_ind_concat,iCon)),10,'r*')

%% calculate fract chage
raw_diff=((pref_responses_stat_concat{post}(:,iCon))-(pref_responses_stat_concat{pre}(:,iCon)));



figure
x = [nanmean(raw_diff(green_ind_concat,iCon)), nanmean(raw_diff(red_ind_concat,iCon))];
y = [(std(raw_diff(green_ind_concat,iCon)))/sqrt(length(green_ind_concat)), (std(raw_diff(red_ind_concat,iCon)))/sqrt(length(red_ind_concat))];

labs =categorical({'HT-','HT+'});
bar(labs,x)                
hold on
er = errorbar(labs,x,-y,y);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
%ylim([0 .2])
hold off
title(['Post - Pre dFoF, contrast = ', num2str(cons(iCon))])

%print(fullfile(fnout,[num2str(cons(iCon)), '_raw_change_resp.pdf']),'-dpdf','-bestfit')

%% time to peak - this is combining running, stationary, and all stimulus conditions

tHalfMax = cell(1,nd);

for id = 1:nd
    tHalfMaxTemp=nan(nKeep_total,1);

        for iCell=1:nKeep_total
            %pull the data for a given cell at a given contrast (pref ori)
            tempData=tc_trial_avrg_stat_concat{id}(stimStart:stimStart+nOn,iCell);
            if resp_keep_concat{id}(iCell)
                if mean(tempData >0)
                smoothData=smoothdata(tempData,'movmean',5) ;
                halfMax = max(smoothData(3:length(smoothData)))/2;
                tHalfMaxCell =double(min(find(smoothData>halfMax)))/double(frame_rate);
                if rem(iCell, 50) == 0 
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
 
       end
    tHalfMax{id}=tHalfMaxTemp;
end
clear tHalfMaxCell tHalfMaxTemp tempData smoothData halfMax
%% scatter for tMax



figure; movegui('center') 
subplot(1,2,1)
scatter((tHalfMax{pre}(green_ind_concat)),(tHalfMax{post}(green_ind_concat)),10,'MarkerEdgeColor',[.4 .4 .4],'jitter', 'on', 'jitterAmount',.1)
hold on
ylabel('post-DART half-max(s)')
xlabel('pre-DART half-max(s)')
ylim([0 2.5])
xlim([0 2.5])
refline(1)
title('HT- ')
axis square
hold on
scatter(nanmean(tHalfMax{pre}(green_ind_concat)),nanmean(tHalfMax{post}(green_ind_concat)),15,'r*')

hold off



subplot(1,2,2)
scatter((tHalfMax{pre}(red_ind_concat)),(tHalfMax{post}(red_ind_concat)),10,'MarkerEdgeColor',[.4 .4 .4],'jitter', 'on', 'jitterAmount',.1)
hold on
ylabel('post-DART half-max(s)')
xlabel('pre-DART half-max(s)')
ylim([0 2.5])
xlim([0 2.5])
refline(1)
title('HT+')
axis square
hold on
scatter(nanmean(tHalfMax{pre}(red_ind_concat)),nanmean(tHalfMax{post}(red_ind_concat)),15,'r*')

hold off

clear txt1 NrPoints
set(gca,'TickDir','out')

sgtitle('time to half max (s), stat')
print(fullfile(fnout,[num2str(cons(iCon)) 'stat_tHalfMax_crossDay.pdf']),'-dpdf','-bestfit')


%% locomotion modulation index


for iCon = 1:nCon

figure; movegui('center') 
subplot(1,2,1)
scatter((LMI_concat{pre}(green_ind_concat,iCon)),(LMI_concat{post}(green_ind_concat,iCon)),10,'MarkerEdgeColor',[.4 .4 .4],'jitter', 'on', 'jitterAmount',.01)
ylabel('post-DART LMI')
xlabel('pre-DART  LMI')
ylim([-1 1])
xlim([-1 1])
hline(0)
vline(0)
refline(1)
hold on
scatter(nanmean(LMI_concat{pre}(green_ind_concat,iCon)),nanmean(LMI_concat{post}(green_ind_concat,iCon)),15,'r*')
axis square
title('HT- ')



subplot(1,2,2)
scatter((LMI_concat{pre}(red_ind_concat,iCon)),(LMI_concat{post}(red_ind_concat,iCon)),10,'MarkerEdgeColor',[.4 .4 .4],'jitter', 'on', 'jitterAmount',.01)
ylabel('post-DART LMI')
xlabel('pre-DART  LMI')
ylim([-1 1])
xlim([-1 1])
hline(0)
vline(0)
refline(1)
hold on
axis square
scatter(nanmean(LMI_concat{pre}(red_ind_concat,iCon)),nanmean(LMI_concat{post}(red_ind_concat,iCon)),15,'r*')
title('HT+')


clear txt1 NrPoints
sgtitle([num2str(cons(iCon))])
set(gca,'TickDir','out')
print(fullfile(fnout,[num2str(cons(iCon)) '_LMI.pdf']),'-dpdf');

end

%% example cell tcs - this is to pull out some individual example cell traces
%
cellList=[49 93]; %enter the cells you're interested in by their index wihtin the keep dataframe


c = linspace(1,10,length(cellList));
figure
scatter((pref_responses_stat_concat{pre}(cellList)),(pref_responses_stat_concat{post}(cellList)),[],c,'filled')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])
colorbar

refline(1)
axis square

%% outputting dataframes for statistical analysis
mouseID=[];
for imouse =1:nSess
   ID_string=string(mice(imouse,1:5));
   thisID = repelem(ID_string,nKeep_concat(imouse));
   mouseID=[mouseID,thisID];
end
mouseID=mouseID';
%% 
mouseIDcol=repmat(mouseID,2,1);
days=[2 1];
day=repelem(days,nKeep_total)';
cell_type=repmat(red_concat,1,2)';
stat_resp=reshape(cell2mat(pref_responses_stat_concat),[(2*nKeep_total),1]);
loc_resp=reshape(cell2mat(pref_responses_loc_concat),[(2*nKeep_total),1]);
half_max=reshape(cell2mat(tHalfMax),[(2*nKeep_total),1]);
LMI=reshape(cell2mat(LMI_concat),[(2*nKeep_total),1]);
output = table(mouseIDcol,day,cell_type,stat_resp,loc_resp,half_max,LMI);



writetable(output,fullfile(fnout,[num2str(targetCon) '_output.csv']))

%%
stat_resp=cell2mat(pref_responses_stat_concat);
HT_ind = red_concat';
save(fullfile(fnout,'DART_dFoF_data.mat'),'stat_resp','HT_ind');
