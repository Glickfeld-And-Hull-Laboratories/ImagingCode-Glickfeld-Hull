
clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc =  behavConstsDART; %directories
eval(ds);
% 136 141 161 153 169 183 177 189
%  138 142 163 171 178 190 24 hour
% 206 210 214 atropine
sess_list = [236];%enter all the sessions you want to concatenate
nSess=length(sess_list);

nd=2;%hard coding for two days per experimental session

% INDICATE THE PRE VS. POST DAYS DEPENDING ON THE ORDER OF MATCHING
prompt = 'Which sesson was used as reference for matching: 0- baseline, 1- post-DART';
            x = input(prompt);
            switch x
                case 0
                    pre=1; %baeline session, used as reference, is in the 1st position
                    post=2;
                    "baseline used as reference"
                case 1
                  pre=2;
                  post=1; %post-DART session, used as reference, is in the 1st position  
                  "post-DART used as reference"
            end
clear x prompt

targetCon = [.25 .5 1]%what contrast to extract for all data - must be one that all datasets had

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
dirs_concat=[];
cons_concat=[];
green_concat=[];
red_concat=[];
dfof_max_diff_concat=[];
nKeep_concat=[];
LMI_concat = cell(1,nd);
data_resp_concat = cell(1,nd);
red_fluor_concat=[];
green_fluor_concat=[];
redR_concat=[];
R_values_concat=[];
wheel_corr_concat=cell(1,nd);
meanF_concat=cell(1,nd);
mean_green_concat=cell(1,nd);
norm_dir_resp_stat_concat = cell(1,nd);
norm_dir_resp_loc_concat = cell(1,nd);
pref_dir_concat=cell(1,nd);


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
    load(fullfile(fn_multi,'HT_pyr_relationship.mat'))
   
    

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
    sharedCon=find(ismember(cons, targetCon));

    nOn = input(1).nScansOn;
    nOff = input(1).nScansOff;
    %start conatenating
    dirs_concat = [dirs_concat,dirs]; %I will organize thisas rows so I can subsequently make sure everything matches
    cons_concat = [cons_concat,cons(sharedCon)];
    red_concat = [red_concat, red_keep_logical];
    green_concat = [green_concat, green_keep_logical];
    nKeep_concat = [nKeep_concat,nKeep];
    redR_concat = [redR_concat,R_red];
    R_values_concat = [R_values_concat, R_p_values];
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
        wheel_corr_concat{id}=cat(2,wheel_corr_concat{id},wheel_corr{id});
        meanF=mean(fullTC_keep{id},1);
       meanF_concat{id}=cat(2,meanF_concat{id}, meanF);
       norm_dir_resp_stat_concat{id}=cat(1,norm_dir_resp_stat_concat{id},norm_dir_resp_stat{id});
       norm_dir_resp_loc_concat{id}=cat(1,norm_dir_resp_loc_concat{id},norm_dir_resp_loc{id});
       pref_dir_concat{id}=cat(2,pref_dir_concat{id},pref_dir_keep{id});
      %mean_green_concat{id}=cat(1,mean_green_concat{id},mean_green_corr(:,id));
        clear meanF
    end
    dfof_max_diff_concat=cat(1,dfof_max_diff_concat,dfof_max_diff(:,sharedCon));
    green_fluor_concat=cat(2,green_fluor_concat,green_fluor_keep);
    red_fluor_concat=cat(2,red_fluor_concat,red_fluor_keep);
end
%
clear mouse day_id nKeep iSess fn_multi cons oris
clear explanation1 resp_keep tc_trial_avrg_keep_allCond pref_responses_allCond sig_diff pref_con_keep pref_ori_keep tOri_match tCon_match data_trial_keep nTrials tc_trial_avrg_keep green_keep_logical red_keep_logical green_ind_keep red_ind_keep
clear LMI RIx locCounts locResp locTCs statResp statTCs wheel_tc
clear data_con_resp_keep data_ori_resp_keep data_rep_keep dfof_max_diff dfof_max_diff_raw explanation2 resp_max_keep data_resp_keep pref_responses_stat pref_responses_loc
clear tc_trial_avrg_stat tc_trial_avrg_loc fullTC_keep norm_dir_resp
clear red_fluor_all red_fluor_match green_fluor_match green_fluor_match red_fluor_keep green_fluor_keep R_p_values
red_ind_concat = find(red_concat);
green_ind_concat = find(green_concat);
%
cons=unique(cons_concat);
nCon = length(cons)
dirs=unique(dirs_concat);
nDir=length(dirs);
nKeep_total = sum(nKeep_concat);
mean(RIx_concat{pre})
mean(RIx_concat{post})
%% find cells that I have running data for on both days
haveRunning_pre = ~isnan(pref_responses_loc_concat{pre});

haveRunning_post = ~isnan(pref_responses_loc_concat{post});
for iCon =1:nCon
haveRunning_both{iCon}= find(haveRunning_pre(:,iCon).* haveRunning_post(:,iCon));
haveRunning_green{iCon} = intersect(haveRunning_both{iCon}, green_ind_concat);
haveRunning_red{iCon} = intersect(haveRunning_both{iCon}, red_ind_concat);
end

clear haveRunning_pre haveRunning_post haveRunning_both


% %%alternate for times when I know I don't have enough running data
% haveRunning_green = green_ind_concat;
% haveRunning_red = red_ind_concat;


%get an array with indices for each mouse
mouseInds=cell(1,nSess);
start=1;
for iMouse = 1:nSess
    mouseInds{iMouse}=start:(start-1)+nKeep_concat(iMouse);
    start = start+nKeep_concat(iMouse)
end
clear start iMouse
% %%alternate for times when I know I don't have enough running data
% haveRunning_green = green_ind_concat;
% haveRunning_red = red_ind_concat;

%find how many haveRunning red cells I got for each mouse
cellCounts = nan(nSess,nCon);
mouseNames=[];
for iMouse = 1:nSess
    for iCon = 1:nCon
        cellCounts(iMouse, iCon,1)=length(intersect(haveRunning_red{iCon},(mouseInds{iMouse})));
        
    end
    mouseNames=[mouseNames, string(mice(iMouse,:))]
end


cellCountsGreen = nan(nSess,nCon);
mouseNames=[];
for iMouse = 1:nSess
    for iCon = 1:nCon
        cellCountsGreen(iMouse, iCon,1)=length(intersect(haveRunning_green{iCon},(mouseInds{iMouse})));
        
    end
    mouseNames=[mouseNames, string(mice(iMouse,:))]
end


cellCountTable = table(cellCounts, RowNames=mouseNames)
cellCountTableGreen = table(cellCountsGreen, RowNames=mouseNames)
%writetable(cellCountTable,fullfile(fnout,['cellCounts.csv']),'WriteRowNames',true)
clear cellCounts cellCountsGreen
%% ALL KEEP CELLS


% make figure with se shaded, one figure per contrast - stationary

tc_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat = cell(1,nd); %same for red
tc_green_se_stat = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat = cell(1,nd); %same for red



for id = 1:nd
  for iCon = 1:nCon
        
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

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

figure
subplot(1,2,1) %for the first day
ylim([-.02 .5]);
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,1),tc_green_se_stat{pre}(:,1),'g');
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,2),tc_green_se_stat{pre}(:,2),'b');
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,3),tc_green_se_stat{pre}(:,3),'c');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['Pre-DART']);

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_green_avrg_stat{post}(:,1),tc_green_se_stat{post}(:,1),'g');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:,2),tc_green_se_stat{post}(:,2),'b');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:,3),tc_green_se_stat{post}(:,3),'c');
ylim([-.02 .5]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['Post-DART'])

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['-HTP stationary contrast responses'])

print(fullfile(fnout,[num2str(cons(iCon)) 'Pyr_contrast_tc.pdf']),'-dpdf'); 




figure
subplot(1,2,1) %for the first day
ylim([-.02 .5]);
hold on
shadedErrorBar(t,tc_red_avrg_stat{pre}(:,1),tc_red_se_stat{pre}(:,1),'g');
hold on
shadedErrorBar(t,tc_red_avrg_stat{pre}(:,2),tc_red_se_stat{pre}(:,2),'b');
hold on
shadedErrorBar(t,tc_red_avrg_stat{pre}(:,3),tc_red_se_stat{pre}(:,3),'c');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['Pre-DART']);

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_stat{post}(:,1),tc_red_se_stat{post}(:,1),'g');
hold on
shadedErrorBar(t,tc_red_avrg_stat{post}(:,2),tc_red_se_stat{post}(:,2),'b');
hold on
shadedErrorBar(t,tc_red_avrg_stat{post}(:,3),tc_red_se_stat{post}(:,3),'c');
ylim([-.02 .5]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['Post-DART'])

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['+HTP stationary contrast responses'])

print(fullfile(fnout,[num2str(cons(iCon)) 'SOM_contrast_tc.pdf']),'-dpdf'); 
%% for the cells that have stationary and running

% make figure with se shaded, one figure per contrast - stationary

tc_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat = cell(1,nd); %same for red
tc_green_se_stat = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat = cell(1,nd); %same for red



for id = 1:nd
    for iCon=1:nCon
        
    tc_green_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,haveRunning_green{iCon},iCon),2);
    green_std=nanstd(tc_trial_avrg_stat_concat{id}(:,haveRunning_green{iCon},iCon),[],2);
    tc_green_se_stat{id}(:,iCon)=green_std/sqrt(length(haveRunning_green{iCon}));
    
    tc_red_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,haveRunning_red{iCon},iCon),2);
    red_std=nanstd(tc_trial_avrg_stat_concat{id}(:,haveRunning_red{iCon},iCon),[],2);
    tc_red_se_stat{id}(:,iCon)=red_std/sqrt(length(haveRunning_red{iCon}));
    
    clear green_std red_std
    end
end
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day



ylim([-.05 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['-HTP',' n = ', num2str(length(haveRunning_green{iCon}))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon),tc_red_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_stat{post}(:,iCon),tc_red_se_stat{post}(:,iCon),'b');
ylim([-.05 .3]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['+HTP',' n = ', num2str(length(haveRunning_red{iCon}))])

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['stationary, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_stat_cellType_timecourses.pdf']),'-dpdf');
end 


%% compare area between curves for all keep red
clear txt1 txt2

%now I have curves, now I want to get the integral of each curve and
%suptract the post from the pre for each contrast

AUC_red = cell(nd,iCon);

for id = 1:nd
    for iCon=1:nCon
        AUC_red{id,iCon}=trapz(tc_trial_avrg_stat_concat{id}(30:60,red_ind_concat,iCon));
    end
end

AUC_diff=cell(nCon,1);
AUC_values=nan(nCon,2);
for iCon = 1:nCon
    AUC_diff{iCon}=AUC_red{pre,iCon}-AUC_red{post,iCon};
    AUC_values(iCon,1)=mean(AUC_diff{iCon});
    AUC_values(iCon,2)=(std(AUC_diff{iCon})/sqrt(size(AUC_diff{iCon},2)));
end




%make the barchart with errobars
figure;
bar(cons,AUC_values(:,1))                
xticks(cons)
%ylim([0 2])
hold on

er = errorbar(cons,AUC_values(:,1),AUC_values(:,2),AUC_values(:,2));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

%% shaded error tc for loc

tc_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_loc = cell(1,nd); %same for red
tc_green_se_loc = cell(1,nd); %this will be the se across all green cells
tc_red_se_loc = cell(1,nd); %same for red

for id = 1:nd
    for iCon=1:nCon
    tc_green_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,haveRunning_green{iCon},iCon),2);
    green_std=nanstd(tc_trial_avrg_loc_concat{id}(:,haveRunning_green{iCon},iCon),[],2);
    tc_green_se_loc{id}(:,iCon)=green_std/sqrt(length(haveRunning_green{iCon}));
    
    tc_red_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,haveRunning_red{iCon},iCon),2);
    red_std=nanstd(tc_trial_avrg_loc_concat{id}(:,haveRunning_red{iCon},iCon),[],2);
    tc_red_se_loc{id}(:,iCon)=red_std/sqrt(length(haveRunning_red{iCon}));
    
    clear green_std red_std
    end
end


%creat a time axis in seconds
t=1:(size(tc_green_avrg_loc{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);


for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day


ylim([-.02 .5]);
hold on
shadedErrorBar(t,tc_green_avrg_loc{pre}(:,iCon),tc_green_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_loc{post}(:,iCon),tc_green_se_loc{post}(:,iCon),'b');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['-HTP',' n = ', num2str(length(haveRunning_green{iCon}))])
ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_loc{pre}(:,iCon),tc_red_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_loc{post}(:,iCon),tc_red_se_loc{post}(:,iCon),'b');
ylim([-.02 .5]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['+HTP',' n = ', num2str(length(haveRunning_red{iCon}))])
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
scatter((pref_responses_stat_concat{pre}(haveRunning_green{iCon},iCon)),(pref_responses_stat_concat{post}(haveRunning_green{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
mean_pre = nanmean(pref_responses_stat_concat{pre}(haveRunning_green{iCon},iCon));
mean_post = nanmean(pref_responses_stat_concat{post}(haveRunning_green{iCon},iCon));
scatter(mean_pre,mean_post,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
stderror_post= std(pref_responses_stat_concat{post}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
line([(mean_pre-stderror_pre),(mean_pre+stderror_pre)],[mean_post, mean_post],'Color','black','LineWidth',2);
line([mean_pre, mean_pre],[(mean_post-stderror_post),(mean_post+stderror_post)],'Color','blue','LineWidth',2);
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
limMin=min(min(pref_responses_stat_concat{pre}(haveRunning_green{iCon},iCon)),min(pref_responses_stat_concat{post}(haveRunning_green{iCon},iCon)));
limMax=max(max(pref_responses_stat_concat{pre}(haveRunning_green{iCon},iCon)),max(pref_responses_stat_concat{post}(haveRunning_green{iCon},iCon)));
ylim([limMin limMax])
xlim([limMin limMax])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('-HTP stationary')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off


subplot(2,2,2)
scatter((pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)),(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
mean_pre = nanmean(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon));
mean_post = nanmean(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon));
scatter(mean_pre,mean_post,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
stderror_post= std(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
line([(mean_pre-stderror_pre),(mean_pre+stderror_pre)],[mean_post, mean_post],'Color','black','LineWidth',2);
line([mean_pre, mean_pre],[(mean_post-stderror_post),(mean_post+stderror_post)],'Color','blue','LineWidth',2);% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
limMin=min(min(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)),min(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)));
limMax=max(max(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)),max(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)));
ylim([limMin limMax])
xlim([limMin limMax])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
title('+HTP stationary')
axis square
hold off

subplot(2,2,3)
scatter((pref_responses_loc_concat{pre}(haveRunning_green{iCon},iCon)),(pref_responses_loc_concat{post}(haveRunning_green{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
mean_pre = nanmean(pref_responses_loc_concat{pre}(haveRunning_green{iCon},iCon));
mean_post = nanmean(pref_responses_loc_concat{post}(haveRunning_green{iCon},iCon));
scatter(mean_pre,mean_post,20,'r','filled');
stderror_pre= std(pref_responses_loc_concat{pre}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
stderror_post= std(pref_responses_loc_concat{post}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
line([(mean_pre-stderror_pre),(mean_pre+stderror_pre)],[mean_post, mean_post],'Color','black','LineWidth',2);
line([mean_pre, mean_pre],[(mean_post-stderror_post),(mean_post+stderror_post)],'Color','blue','LineWidth',2);
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
limMin=min(min(pref_responses_loc_concat{pre}(haveRunning_green{iCon},iCon)),min(pref_responses_loc_concat{post}(haveRunning_green{iCon},iCon)));
limMax=max(max(pref_responses_loc_concat{pre}(haveRunning_green{iCon},iCon)),max(pref_responses_loc_concat{post}(haveRunning_green{iCon},iCon)));
ylim([limMin limMax])
xlim([limMin limMax])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('-HTP running')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off

subplot(2,2,4)
scatter((pref_responses_loc_concat{pre}(haveRunning_red{iCon},iCon)),(pref_responses_loc_concat{post}(haveRunning_red{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
mean_pre = nanmean(pref_responses_loc_concat{pre}(haveRunning_red{iCon},iCon));
mean_post = nanmean(pref_responses_loc_concat{post}(haveRunning_red{iCon},iCon));
scatter(mean_pre,mean_post,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
stderror_post= std(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
line([(mean_pre-stderror_pre),(mean_pre+stderror_pre)],[mean_post, mean_post],'Color','black','LineWidth',2);
line([mean_pre, mean_pre],[(mean_post-stderror_post),(mean_post+stderror_post)],'Color','blue','LineWidth',2);
xlabel('pre-DART  dF/F')
limMin=min(min(pref_responses_loc_concat{pre}(haveRunning_red{iCon},iCon)),min(pref_responses_loc_concat{post}(haveRunning_red{iCon},iCon)));
limMax=max(max(pref_responses_loc_concat{pre}(haveRunning_red{iCon},iCon)),max(pref_responses_loc_concat{post}(haveRunning_red{iCon},iCon)));
ylim([limMin limMax])
xlim([limMin limMax])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('+HTP running')
uistack(hline,'bottom');
axis square
hold off
set(gca, 'TickDir', 'out')



sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'maxResp_crossDay.pdf']),'-dpdf','-bestfit')
clear mean_pre mean_post stderror_post stderror_pre
end
%% t-test for the lowest contrast

[h1, p1]= ttest(pref_responses_stat_concat{pre}(haveRunning_green{1},1),pref_responses_stat_concat{post}(haveRunning_green{1},1));
[h2,p2]= ttest(pref_responses_stat_concat{pre}(haveRunning_red{1},1),pref_responses_stat_concat{post}(haveRunning_red{1},1));
[h3,p3]= ttest(pref_responses_loc_concat{pre}(haveRunning_green{1},1),pref_responses_loc_concat{post}(haveRunning_green{1},1));
[h4,p4]= ttest(pref_responses_loc_concat{pre}(haveRunning_red{1},1),pref_responses_loc_concat{post}(haveRunning_red{1},1));

%correct for four tests
p1*2
p2*2


%clear h1 p1 h2 p2 h3 p3 h4 p4
%%
responseTable = table([nanmean(pref_responses_stat_concat{pre}(haveRunning_green));nanmean(pref_responses_stat_concat{post}(haveRunning_green))],[nanmean(pref_responses_stat_concat{pre}(haveRunning_red));nanmean(pref_responses_stat_concat{post}(haveRunning_red))],[nanmean(pref_responses_loc_concat{pre}(haveRunning_green));nanmean(pref_responses_loc_concat{post}(haveRunning_green))],[nanmean(pref_responses_loc_concat{pre}(haveRunning_red));nanmean(pref_responses_loc_concat{post}(haveRunning_red))],'VariableNames',{'Pyramidal cells stat'  '+HTP cells stat' 'Pyramidal cells loc'  '+HTP cells loc'}, 'RowNames',{'Pre'  'Post'})
writetable(responseTable,fullfile(fnout,[num2str(targetCon) 'responseTable.csv']),'WriteRowNames',true)
%% polar plots of direction tuning
dirResp_for_plotting=cell(1,nd);
%adjusting the normalized response matrix for plotting
for id = 1:nd
    dirResp_for_plotting{id}=norm_dir_resp_stat_concat{id}; %copy the normalized response matrix for that day
    dirResp_for_plotting{id}(:,nDir+1)=dirResp_for_plotting{id}(:,1); %take the first element and re-paste it onto the end

end
dirs_for_plotting = deg2rad([dirs, dirs(1)]);


figure;
subplot(1,2,1)
polarplot(dirs_for_plotting,mean(dirResp_for_plotting{pre}(green_ind_concat,:),1),'k','LineWidth',2)
hold on
polarplot(dirs_for_plotting,mean(dirResp_for_plotting{post}(green_ind_concat,:),1),'b','LineWidth',2)
set(gca, 'TickDir', 'out')
title('-HTP')

subplot(1,2,2)
polarplot(dirs_for_plotting,mean(dirResp_for_plotting{pre}(red_ind_concat,:),1),'k','LineWidth',2)
hold on
polarplot(dirs_for_plotting,mean(dirResp_for_plotting{post}(red_ind_concat,:),1),'b','LineWidth',2)
set(gca, 'TickDir', 'out')
title('+HTP')

sgtitle("Normalized direction tuning, averaged over contrast")


%% linear direction plot 

dirs_for_plotting=dirs-180;

green_dir_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
red_dir_avrg_stat = cell(1,nd); %same for red
green_dir_se_stat = cell(1,nd); %this will be the se across all green cells
red_dir_se_stat = cell(1,nd); %same for red

for id = 1:nd
   
    green_dir_avrg_stat{id}=nanmean(norm_dir_resp_stat_concat{id}(green_ind_concat,:),1);
    green_std=nanstd(norm_dir_resp_stat_concat{id}(green_ind_concat,:),[],1);
    green_dir_se_stat{id}=green_std/sqrt(length(green_ind_concat));
    green_dir_avrg_stat{id}=circshift(green_dir_avrg_stat{id},4);
    green_dir_se_stat{id}=circshift(green_dir_se_stat{id},4);
    
    red_dir_avrg_stat{id}=nanmean(norm_dir_resp_stat_concat{id}(red_ind_concat,:),1);
    red_std=nanstd(norm_dir_resp_stat_concat{id}(red_ind_concat,:),[],1);
    red_dir_se_stat{id}=red_std/sqrt(length(red_ind_concat));
    red_dir_avrg_stat{id}=circshift(red_dir_avrg_stat{id},4);
    red_dir_se_stat{id}=circshift(red_dir_se_stat{id},4);
    clear green_std red_std
    
end



green_dir_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
red_dir_avrg_loc = cell(1,nd); %same for red
green_dir_se_loc = cell(1,nd); %this will be the se across all green cells
red_dir_se_loc = cell(1,nd); %same for red

for id = 1:nd
   
    green_dir_avrg_loc{id}=nanmean(norm_dir_resp_loc_concat{id}(green_ind_concat,:),1);
    green_std=nanstd(norm_dir_resp_loc_concat{id}(green_ind_concat,:),[],1);
    green_dir_se_loc{id}=green_std/sqrt(length(green_ind_concat));
    green_dir_avrg_loc{id}=circshift(green_dir_avrg_loc{id},4);
    green_dir_se_loc{id}=circshift(green_dir_se_loc{id},4);
    
    red_dir_avrg_loc{id}=nanmean(norm_dir_resp_loc_concat{id}(red_ind_concat,:),1);
    red_std=nanstd(norm_dir_resp_loc_concat{id}(red_ind_concat,:),[],1);
    red_dir_se_loc{id}=red_std/sqrt(length(red_ind_concat));
    red_dir_avrg_loc{id}=circshift(red_dir_avrg_loc{id},4);
    red_dir_se_loc{id}=circshift(red_dir_se_loc{id},4);
    clear green_std red_std
    
end



figure
subplot(2,2,1)
errorbar(dirs_for_plotting,green_dir_avrg_stat{pre},green_dir_se_stat{pre},'k')
hold on
errorbar(dirs_for_plotting,green_dir_avrg_stat{post},green_dir_se_stat{post},'b')
title('-HTP, stationary')
set(gca, 'TickDir', 'out')
axis square
box off
ylabel('dF/F')
xlabel('normalized direction')
%ylim([0 .25])

subplot(2,2,2)
errorbar(dirs_for_plotting,red_dir_avrg_stat{pre},red_dir_se_stat{pre},'k')
hold on
errorbar(dirs_for_plotting,red_dir_avrg_stat{post},red_dir_se_stat{post},'b')
title('+HTP, stationary')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('normalized direction')
%ylim([0 .25])


subplot(2,2,3)
errorbar(dirs_for_plotting,green_dir_avrg_loc{pre},green_dir_se_loc{pre},'k')
hold on
errorbar(dirs_for_plotting,green_dir_avrg_loc{post},green_dir_se_loc{post},'b')
title('-HTP, running')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('normalized direction')
%ylim([0 .25])

subplot(2,2,4)
errorbar(dirs_for_plotting,red_dir_avrg_loc{pre},red_dir_se_loc{pre},'k')
hold on
errorbar(dirs_for_plotting,red_dir_avrg_loc{post},red_dir_se_loc{post},'b')
title('+HTP, running')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('normalized direction')
%ylim([0 .25])

sgtitle("Normalized direction tuning for -HTP Putative Pyr cells")


print(fullfile(fnout,['dirTuning.pdf']),'-dpdf','-bestfit')
%% subtracted tuning curve

subtracted_dir = cell(1,nd);

for id = 1:nd
    subtracted_dir{id}=norm_dir_resp_loc_concat{id}-norm_dir_resp_stat_concat{id};

end



green_dir_avrg_subtracted = cell(1,nd); %this will be the average across all green cells - a single line
red_dir_avrg_subtracted = cell(1,nd); %same for red
green_dir_se_subtracted = cell(1,nd); %this will be the se across all green cells
red_dir_se_subtracted = cell(1,nd); %same for red

for id = 1:nd
   
    green_dir_avrg_subtracted{id}=nanmean(subtracted_dir{id}(green_ind_concat,:),1);
    green_std=nanstd(subtracted_dir{id}(green_ind_concat,:),[],1);
    green_dir_se_subtracted{id}=green_std/sqrt(length(green_ind_concat));
    green_dir_avrg_subtracted{id}=circshift(green_dir_avrg_subtracted{id},4);
    green_dir_se_subtracted{id}=circshift(green_dir_se_subtracted{id},4);
    
    red_dir_avrg_subtracted{id}=nanmean(subtracted_dir{id}(red_ind_concat,:),1);
    red_std=nanstd(subtracted_dir{id}(red_ind_concat,:),[],1);
    red_dir_se_subtracted{id}=red_std/sqrt(length(red_ind_concat));
    red_dir_avrg_subtracted{id}=circshift(red_dir_avrg_subtracted{id},4);
    red_dir_se_subtracted{id}=circshift(red_dir_se_subtracted{id},4);
    clear green_std red_std
    
end

figure;
errorbar(dirs,green_dir_avrg_subtracted{pre},green_dir_se_subtracted{pre},'k')
hold on
errorbar(dirs,green_dir_avrg_subtracted{post},green_dir_avrg_subtracted{post},'b')
title('-HTP, post')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('normalized direction')
ylim([0 .2])


%% comparing the preferred direction across days

change_pref = pref_dir_concat{pre}-pref_dir_concat{post};
change_pref=deg2rad(change_pref);
figure
subplot(1,2,1)
polarhistogram(change_pref(green_ind_concat),'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5])
title('-HTP')
set(gca, 'TickDir', 'out')

subplot(1,2,2)
polarhistogram(change_pref(red_ind_concat),'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5])
title('+HTP')
set(gca, 'TickDir', 'out')
print(fullfile(fnout,['prefDir_change.pdf']),'-dpdf','-bestfit')

%% response by condition
%finds the cells that have running and stationary for all three contrasts
green_all = intersect(haveRunning_green{1},haveRunning_green{2});
green_all = intersect(green_all, haveRunning_green{3});

red_all = intersect(haveRunning_red{1},haveRunning_red{2});
red_all = intersect(red_all, haveRunning_red{3});


a=mean(pref_responses_stat_concat{pre}(green_all,:), "omitnan");
b=mean(pref_responses_loc_concat{pre}(green_all,:), "omitnan");

c=mean(pref_responses_stat_concat{pre}(red_all,:), "omitnan");
d=mean(pref_responses_loc_concat{pre}(red_all,:), "omitnan");

e=mean(pref_responses_stat_concat{post}(green_all,:), "omitnan");
f=mean(pref_responses_loc_concat{post}(green_all,:), "omitnan");

g=mean(pref_responses_stat_concat{post}(red_all,:), "omitnan");
h=mean(pref_responses_loc_concat{post}(red_all,:), "omitnan");

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
set(gca,'TickDir','out')
xlabel('Mean Pyr dF/F')
ylabel('Mean SOM dF/F')
title('Mean response to different conditions')
print(fullfile(fnout, ['responseByCondition.pdf']),'-dpdf','-bestfit')
%% plotting vs contrast
figure;
scatter(cons,mean(pref_responses_stat_concat{pre}(red_all,:), "omitnan"),"k")
hold on
scatter(cons,mean(pref_responses_stat_concat{post}(red_all,:), "omitnan"),"b")
hold on 
scatter(cons,mean(pref_responses_loc_concat{pre}(red_all,:), "omitnan"),"k",'MarkerFaceColor', 'k')
hold on
scatter(cons,mean(pref_responses_loc_concat{post}(red_all,:), "omitnan"),"b",'MarkerFaceColor', 'b')
set(gca,'TickDir','out')

xlabel('Contrast')
ylabel('Mean SOM dF/F')
title('Mean response vs contrast')
hold off
%% calculate fract chage


for iCon = 1:nCon
raw_diff(:,iCon)=((pref_responses_stat_concat{post}(:,iCon))-(pref_responses_stat_concat{pre}(:,iCon)));
fract_diff(:,iCon)=((pref_responses_stat_concat{post}(:,iCon))-(pref_responses_stat_concat{pre}(:,iCon))./(pref_responses_stat_concat{pre}(:,iCon)));

figure
x = [nanmean(fract_diff(haveRunning_green{iCon})), nanmean(fract_diff(haveRunning_red{iCon}))];
y = [(std(fract_diff(haveRunning_green{iCon})))/sqrt(length(haveRunning_green)), (std(fract_diff(haveRunning_red{iCon})))/sqrt(length(haveRunning_red))];

labs =categorical({'+HTP','-HTP'});
bar(labs,x)                
hold on
er = errorbar(labs,x,-y,y);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
%ylim([0 .2])
hold off
title(['(Post - Pre)/Pre dFoF, contrast = ', num2str(cons(iCon))])

end
%print(fullfile(fnout,[num2str(cons(iCon)), '_raw_change_resp.pdf']),'-dpdf','-bestfit')


diff_values=nan(nCon,2);
for iCon = 1:nCon
    diff_values(iCon,1)=mean(fract_diff(red_ind_concat,iCon),1);
    diff_values(iCon,2)=(std(fract_diff(red_ind_concat,iCon))/sqrt(length(red_ind_concat)));
end




%make the barchart with errobars
figure;
bar(cons,diff_values(:,1))                
xticks(cons)
%ylim([0 2])
hold on

er = errorbar(cons,diff_values(:,1),diff_values(:,2),diff_values(:,2));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off
print(fullfile(fnout,[num2str(cons(iCon)), '_raw_change_resp.pdf']),'-dpdf','-bestfit')
%% comparing R value to change in dF/F

%averaging over contrasts -  baseline R value vs. raw change
figure;
scatter(R_values_concat(1,:),(mean(raw_diff(red_ind_concat,:),2)));
xlabel("Baseline R")
ylabel("Raw post-pre")
set(gca,'TickDir','out')
box off
title("Averaged over contrast")

%averaging over contrasts -  baseline R value vs. fractional change
figure;
scatter(R_values_concat(1,:),(mean(fract_diff(red_ind_concat,:),2)));
xlabel("Baseline R")
ylabel("Fractional post-pre")
set(gca,'TickDir','out')
box off
title("Averaged over contrast")

%only 25% contrast -  baseline R value vs. fractional change
figure;
scatter(R_values_concat(1,:),fract_diff(red_ind_concat,1));
xlabel("Baseline R")
ylabel("Fractional post-pre")
set(gca,'TickDir','out')
box off
title("25% contrast")

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
 
       end
    tHalfMax{id}=tHalfMaxTemp;
end
clear tHalfMaxCell tHalfMaxTemp tempData smoothData halfMax
%% scatter for tMax



figure; movegui('center') 
subplot(1,2,1)
scatter((tHalfMax{pre}(green_all)),(tHalfMax{post}(green_all)),10,'MarkerEdgeColor',[.4 .4 .4],'jitter', 'on', 'jitterAmount',.1)
hold on
mean_pre = nanmean(tHalfMax{pre}(green_all));
mean_post = nanmean(tHalfMax{post}(green_all));
scatter(mean_pre,mean_post,20,'r','filled');
stderror_pre= nanstd(tHalfMax{pre}(green_all)) / sqrt( length(green_all));
stderror_post= nanstd(tHalfMax{post}(green_all)) / sqrt( length(green_all));
line([(mean_pre-stderror_pre),(mean_pre+stderror_pre)],[mean_post, mean_post],'Color','black','LineWidth',2);
line([mean_pre, mean_pre],[(mean_post-stderror_post),(mean_post+stderror_post)],'Color','blue','LineWidth',2);
ylabel('post-DART half-max(s)')
xlabel('pre-DART half-max(s)')
% ylim([0 2.15])
% xlim([0 2.15])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('-HTP ')
axis square


hold off



subplot(1,2,2)
scatter((tHalfMax{pre}(red_all)),(tHalfMax{post}(red_all)),10,'MarkerEdgeColor',[.4 .4 .4],'jitter', 'on', 'jitterAmount',.1)
hold on
mean_pre = nanmean(tHalfMax{pre}(red_all));
mean_post = nanmean(tHalfMax{post}(red_all));
scatter(mean_pre,mean_post,20,'r','filled');
stderror_pre= nanstd(tHalfMax{pre}(red_all)) / sqrt( length(red_all));
stderror_post= nanstd(tHalfMax{post}(red_all)) / sqrt( length(red_all));
line([(mean_pre-stderror_pre),(mean_pre+stderror_pre)],[mean_post, mean_post],'Color','black','LineWidth',2);
line([mean_pre, mean_pre],[(mean_post-stderror_post),(mean_post+stderror_post)],'Color','blue','LineWidth',2);

%ylabel('post-DART half-max(s)')
xlabel('pre-DART half-max(s)')
ylim([0 2.15])
xlim([0 2.15])
hline=refline(1)
hline.Color = 'k';
hline.LineStyle = ':';
title('+HTP')
axis square
hold off

clear txt1 NrPoints
set(gca,'TickDir','out')

sgtitle('time to half max (s), stat')
print(fullfile(fnout,[ 'stat_tHalfMax_crossDay.pdf']),'-dpdf','-bestfit')

[h1, p1]= ttest(tHalfMax{pre}(green_all),tHalfMax{post}(green_all));
[h2,p2]= ttest(tHalfMax{pre}(red_all),tHalfMax{post}(red_all));

%correct for four tests
p1*2
p2*2
clear h1 p1 h2 p2

%% locomotion modulation index


for iCon = 1:nCon

figure; movegui('center') 
subplot(1,2,1)
scatter((LMI_concat{pre}(green_ind_concat,iCon)),(LMI_concat{post}(green_ind_concat,iCon)),10,'MarkerEdgeColor',[.4 .4 .4],'jitter', 'on', 'jitterAmount',.01)
ylabel('post-DART LMI')
xlabel('pre-DART  LMI')
ylim([-1 1])
xlim([-1 1])
% hline(0)
% vline(0)
refline(1)
hold on
scatter(nanmean(LMI_concat{pre}(green_ind_concat,iCon)),nanmean(LMI_concat{post}(green_ind_concat,iCon)),15,'r*')
axis square
title('-HTP ')



subplot(1,2,2)
scatter((LMI_concat{pre}(red_ind_concat,iCon)),(LMI_concat{post}(red_ind_concat,iCon)),10,'MarkerEdgeColor',[.4 .4 .4],'jitter', 'on', 'jitterAmount',.01)
ylabel('post-DART LMI')
xlabel('pre-DART  LMI')
ylim([-1 1])
xlim([-1 1])
% hline(0)
% vline(0)
refline(1)
hold on
axis square
scatter(nanmean(LMI_concat{pre}(red_ind_concat,iCon)),nanmean(LMI_concat{post}(red_ind_concat,iCon)),15,'r*')
title('+HTP')


clear txt1 NrPoints
sgtitle([num2str(cons(iCon))])
set(gca,'TickDir','out')
print(fullfile(fnout,[num2str(cons(iCon)) '_LMI.pdf']),'-dpdf');

end
%% compare r value to LMI
meanLMI = mean(LMI_concat{pre}(red_ind_concat,:),2);

figure
scatter(R_values_concat(1,:),meanLMI)
xlabel('r pre-DART')
ylabel('LMI pre-DART')
title('+HTP cells')
lsline()

print(fullfile(fnout,[ 'r_vs_LMI.pdf']),'-dpdf');

%% looking at R value and change in R value
figure;
subplot(2,1,1)
histogram(R_values_concat(1,:))
xlim([-1 1])
subplot(2,1,2)
histogram(R_values_concat(3,:))
xlim([-1 1])


figure;
subplot(2,1,1)
scatter(R_values_concat(1,:),R_values_concat(3,:));
xlabel("Pre R value");
ylabel("Post R value");
xlim([-.25 1])
ylim([-.25 1])
lsline
refline(1)
axis square
set(gca, 'TickDir', 'out')

subplot(2,1,2)
scatter(R_values_concat(1,:),(R_values_concat(3,:)-R_values_concat(1,:)));
xlabel("Pre R value");
ylabel("raw delta R value");
ylim([-1 1])
lsline
axis square
set(gca, 'TickDir', 'out')
print(fullfile(fnout,['changeInR.pdf']),'-dpdf');

%% making tc plots for low and high R cells
highRInds=red_ind_concat(find(redR_concat)); %find the indices of red cells with high R
lowRInds=red_ind_concat(find(~redR_concat)); %find the indices of red cells with low R

%alternative version with full TC correlation rather than stimulus response
%correlation
% highRsqInds=red_ind_concat(find(mean_green_concat{pre}>median(mean_green_concat{pre}))); %find the indices of red cells with high R
% lowRsqInds=red_ind_concat(find(mean_green_concat{pre}<median(mean_green_concat{pre}))); %find the indices of red cells with low R


%% stat high and low R
hi_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
low_avrg_stat = cell(1,nd); %same for red
hi_se_stat = cell(1,nd); %this will be the se across all green cells
low_se_stat = cell(1,nd); %same for red



for id = 1:nd
    for iCon=1:nCon
   redHigh =  intersect(haveRunning_red{iCon},highRInds); %find intersection of red cells with high Rsq 
   % and red cells that I have running data for on both days, for this contrast
    redLow =  intersect(haveRunning_red{iCon},lowRInds); %find intersection of red cells with low Rsq 
   % and red cells that I have running data for on both days, for this contrast

    hi_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:, redHigh,iCon),2);
    high_std=nanstd(tc_trial_avrg_stat_concat{id}(:, redHigh,iCon),[],2);
    hi_se_stat{id}(:,iCon)=high_std/sqrt(length(redHigh));
    
    low_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,redLow,iCon),2);
    low_std=nanstd(tc_trial_avrg_stat_concat{id}(:, redHigh,iCon),[],2);
    low_se_stat{id}(:,iCon)=low_std/sqrt(length(redLow));
    
    clear low_std high_std
    end
end
z=double(nOn)/double(frame_rate);

%creat a time axis in seconds
t=1:(size(hi_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for -HTP



ylim([-.02 .5]);
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
title(['+HTP high r',' n = ', num2str(length(redHigh))])
ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for +HTP
shadedErrorBar(t,low_avrg_stat{pre}(:,iCon),low_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,low_avrg_stat{post}(:,iCon),low_se_stat{post}(:,iCon),'b');
ylim([-.02 .5]);
hold on
% line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
% hold on
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
%ylabel('dF/F') 
xlabel('s') 
title(['+HTP low r',' n = ', num2str(length(redLow))])
x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['stationary, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) 'stat_Rsq_timecourses.pdf']),'-dpdf');
clear txt1 highRed lowRed
end 


%% scatterplot of max df/f for day 1 vs day 2 for high and low R red cells

for iCon = 1:nCon
   redHigh =  intersect(haveRunning_red{iCon},highRInds); %find intersection of red cells with high Rsq 
   % and red cells that I have running data for on both days, for this contrast
    redLow =  intersect(haveRunning_red{iCon},lowRInds); %find intersection of red cells with low Rsq 
   % and red cells that I have running data for on both days, for this contrast


figure; movegui('center') 
subplot(1,2,2)
scatter((pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)),(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
scatter((pref_responses_stat_concat{pre}(redHigh,iCon)),(pref_responses_stat_concat{post}(redHigh,iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on

mean_pre = nanmean(pref_responses_stat_concat{pre}(redHigh,iCon));
mean_post = nanmean(pref_responses_stat_concat{post}(redHigh,iCon));
scatter(mean_pre,mean_post,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(redHigh,iCon)) / sqrt( length(redHigh));
stderror_post= std(pref_responses_stat_concat{post}(redHigh,iCon)) / sqrt( length(redHigh));
line([(mean_pre-stderror_pre),(mean_pre+stderror_pre)],[mean_post, mean_post],'Color','black','LineWidth',2);
line([mean_pre, mean_pre],[(mean_post-stderror_post),(mean_post+stderror_post)],'Color','blue','LineWidth',2);


xlabel('pre-DART  dF/F')
ylim([-.1 0.4])
xlim([-.1 0.4])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('+HTP high r, stationary')
axis square
set(gca, 'TickDir', 'out')



subplot(1,2,1)
scatter((pref_responses_stat_concat{pre}(redLow,iCon)),(pref_responses_stat_concat{post}(redLow,iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
mean_pre = nanmean(pref_responses_stat_concat{pre}(redLow,iCon));
mean_post = nanmean(pref_responses_stat_concat{post}(redLow,iCon));
scatter(mean_pre,mean_post,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(redLow,iCon)) / sqrt( length(redLow));
stderror_post= std(pref_responses_stat_concat{post}(redLow,iCon)) / sqrt( length(redLow));
line([(mean_pre-stderror_pre),(mean_pre+stderror_pre)],[mean_post, mean_post],'Color','black','LineWidth',2);
line([mean_pre, mean_pre],[(mean_post-stderror_post),(mean_post+stderror_post)],'Color','blue','LineWidth',2);% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylabel('post-DART dF/F')
ylim([-.1 0.4])
xlim([-.1 0.4])
set(gca, 'TickDir', 'out')
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('+HTP low r, stationary')
axis square


sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'stat_maxResp_Rsq.pdf']),'-dpdf','-bestfit')
clear txt1 highRed lowRed
end

[h1, p1]= ttest(pref_responses_stat_concat{pre}(redLow,1),pref_responses_stat_concat{post}(redLow,1));
[h2,p2]= ttest(pref_responses_stat_concat{pre}(redHigh,1),pref_responses_stat_concat{post}(redHigh,1));

%correct for four tests
p1*2
p2*2
clear h1 p1 h2 p2

%% making tc plots for low and high R cells, running
highRInds=red_ind_concat(find(redR_concat)); %find the indices of red cells with high Rsq
lowRInds=red_ind_concat(find(~redR_concat)); %find the indices of red cells with low Rsq

hi_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
low_avrg_loc = cell(1,nd); %same for red
hi_se_loc = cell(1,nd); %this will be the se across all green cells
low_se_loc = cell(1,nd); %same for red



for id = 1:nd
    for iCon=1:nCon
   redHigh =  intersect(haveRunning_red{iCon},highRInds); %find intersection of red cells with high Rsq 
   % and red cells that I have running data for on both days, for this contrast
    redLow =  intersect(haveRunning_red{iCon},lowRInds); %find intersection of red cells with low Rsq 
   % and red cells that I have running data for on both days, for this contrast

    hi_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:, redHigh,iCon),2);
    high_std=nanstd(tc_trial_avrg_loc_concat{id}(:, redHigh,iCon),[],2);
    hi_se_loc{id}(:,iCon)=high_std/sqrt(length(redHigh));
    
    low_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,redLow,iCon),2);
    low_std=nanstd(tc_trial_avrg_loc_concat{id}(:, redHigh,iCon),[],2);
    low_se_loc{id}(:,iCon)=low_std/sqrt(length(redLow));
    
    clear low_std high_std
    end
end
z=double(nOn)/double(frame_rate);

%creat a time axis in seconds
t=1:(size(hi_avrg_loc{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for -HTP



ylim([-.02 .5]);
hold on
shadedErrorBar(t,hi_avrg_loc{pre}(:,iCon),hi_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,hi_avrg_loc{post}(:,iCon),hi_se_loc{post}(:,iCon),'b');
hold on
% line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
% hold on
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['+HTP high r',' n = ', num2str(length(redHigh))])
ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for +HTP
shadedErrorBar(t,low_avrg_loc{pre}(:,iCon),low_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,low_avrg_loc{post}(:,iCon),low_se_loc{post}(:,iCon),'b');
ylim([-.02 .5]);
hold on
% line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
% hold on
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
%ylabel('dF/F') 
xlabel('s') 
title(['+HTP low r',' n = ', num2str(length(redLow))])
x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['running, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) 'loc_Rsq_timecourses.pdf']),'-dpdf');
clear txt1 highRed lowRed
end 


%% scatterplot of max df/f for day 1 vs day 2 for high and low Rsq red cells, running

for iCon = 1:nCon
   redHigh =  intersect(haveRunning_red{iCon},highRInds); %find intersection of red cells with high Rsq 
   % and red cells that I have running data for on both days, for this contrast
    redLow =  intersect(haveRunning_red{iCon},lowRInds); %find intersection of red cells with low Rsq 
   % and red cells that I have running data for on both days, for this contrast


figure; movegui('center') 
subplot(1,2,2)
scatter((pref_responses_loc_concat{pre}(haveRunning_red{iCon},iCon)),(pref_responses_loc_concat{post}(haveRunning_red{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
scatter((pref_responses_loc_concat{pre}(redHigh,iCon)),(pref_responses_loc_concat{post}(redHigh,iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on

mean_pre = nanmean(pref_responses_loc_concat{pre}(redHigh,iCon));
mean_post = nanmean(pref_responses_loc_concat{post}(redHigh,iCon));
scatter(mean_pre,mean_post,20,'r','filled');
stderror_pre= std(pref_responses_loc_concat{pre}(redHigh,iCon)) / sqrt( length(redHigh));
stderror_post= std(pref_responses_loc_concat{post}(redHigh,iCon)) / sqrt( length(redHigh));
line([(mean_pre-stderror_pre),(mean_pre+stderror_pre)],[mean_post, mean_post],'Color','black','LineWidth',2);
line([mean_pre, mean_pre],[(mean_post-stderror_post),(mean_post+stderror_post)],'Color','blue','LineWidth',2);


xlabel('pre-DART  dF/F')
ylim([-.1 0.4])
xlim([-.1 0.4])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('+HTP high r, running')
axis square
set(gca, 'TickDir', 'out')



subplot(1,2,1)
scatter((pref_responses_loc_concat{pre}(redLow,iCon)),(pref_responses_loc_concat{post}(redLow,iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
mean_pre = nanmean(pref_responses_loc_concat{pre}(redLow,iCon));
mean_post = nanmean(pref_responses_loc_concat{post}(redLow,iCon));
scatter(mean_pre,mean_post,20,'r','filled');
stderror_pre= std(pref_responses_loc_concat{pre}(redLow,iCon)) / sqrt( length(redLow));
stderror_post= std(pref_responses_loc_concat{post}(redLow,iCon)) / sqrt( length(redLow));
line([(mean_pre-stderror_pre),(mean_pre+stderror_pre)],[mean_post, mean_post],'Color','black','LineWidth',2);
line([mean_pre, mean_pre],[(mean_post-stderror_post),(mean_post+stderror_post)],'Color','blue','LineWidth',2);% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylabel('post-DART dF/F')
ylim([-.1 0.4])
xlim([-.1 0.4])
set(gca, 'TickDir', 'out')
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('+HTP low r, running')
axis square


sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'loc_maxResp_Rsq.pdf']),'-dpdf','-bestfit')
clear txt1 highRed lowRed
end

[h1, p1]= ttest(pref_responses_loc_concat{pre}(redLow,3),pref_responses_loc_concat{post}(redLow,3));
[h2,p2]= ttest(pref_responses_loc_concat{pre}(redHigh,3),pref_responses_loc_concat{post}(redHigh,3));

%correct for four tests
p1*2
p2*2
clear h1 p1 h2 p2

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
green_col = repmat(green_fluor_concat,1,2)';
red_col = repmat(red_fluor_concat,1,2)';


output = table(mouseIDcol,day,cell_type,stat_resp,loc_resp,half_max,LMI,green_col,red_col);

writetable(output,fullfile(fnout,[num2str(targetCon) '_output.csv']))

%%
stat_resp=cell2mat(pref_responses_stat_concat);
HT_ind = red_concat';
green_intensity = green_fluor_concat';
red_intensity = red_fluor_concat';
save(fullfile(fnout,'DART_dFoF_data.mat'),'stat_resp','HT_ind','mouseID','green_intensity','red_intensity');

%% mouse by mouse means scatterplot of max df/f for day 1 vs day 2, and each subplot is one cell type
green_means_stat=cell(1,nd);
red_means_stat=cell(1,nd);
green_means_loc=cell(1,nd);
red_means_loc=cell(1,nd);

for id=1:nd
    temp_means1 = nan(nSess,nCon);
    temp_means2 = nan(nSess,nCon);
    temp_means3 = nan(nSess,nCon);
    temp_means4 = nan(nSess,nCon);
    for iMouse = 1:nSess
        for iCon = 1:nCon
            inds1 = intersect(haveRunning_green{iCon}, mouseInds{iMouse});
            temp_means1(iMouse, iCon)=mean(pref_responses_stat_concat{id}(inds1,iCon),'omitnan');

            inds2 = intersect(haveRunning_red{iCon}, mouseInds{iMouse});
            temp_means2(iMouse, iCon)=mean(pref_responses_stat_concat{id}(inds2,iCon),'omitnan');

            temp_means3(iMouse, iCon)=mean(pref_responses_loc_concat{id}(inds1,iCon),'omitnan');

            temp_means4(iMouse, iCon)=mean(pref_responses_loc_concat{id}(inds2,iCon),'omitnan');
        end
    end
    green_means_stat{id}=temp_means1;
    red_means_stat{id}=temp_means2;
    green_means_loc{id}=temp_means3;
    red_means_loc{id}=temp_means4;
end

cmap = winter(nSess); % Make colors


for iCon = 1:nCon
figure; movegui('center') 
subplot(2,2,1)
scatter((green_means_stat{pre}(:,iCon)),(green_means_stat{post}(:,iCon)),10,cmap, 'filled')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
% ylim([0 .2])
% xlim([0 .2])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('-HTP stationary')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off


subplot(2,2,2)
scatter(red_means_stat{pre}(:,iCon),red_means_stat{post}(:,iCon),10,cmap, 'filled')
hold on
% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
% ylim([0 .2])
% xlim([0 .2])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
title('+HTP stationary')
axis square
hold off

subplot(2,2,3)
scatter((green_means_loc{pre}(:,iCon)),(green_means_loc{post}(:,iCon)),10,cmap, 'filled')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
% ylim([0 .4])
% xlim([0 .4])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('-HTP running')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off

subplot(2,2,4)
scatter((red_means_loc{pre}(:,iCon)),(red_means_loc{post}(:,iCon)),10,cmap, 'filled')
hold on
% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
% ylim([0 .4])
% xlim([0 .4])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('+HTP running')
uistack(hline,'bottom');
axis square
hold off
set(gca, 'TickDir', 'out')


sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'maxResp_byMouse.pdf']),'-dpdf','-bestfit')
end
%% relationship of full timecourse to locomotion


figure; movegui('center') 
subplot(1,2,1)
scatter(wheel_corr_concat{pre}(green_ind_concat),wheel_corr_concat{post}(green_ind_concat),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
ylabel('post-DART running corr')
xlabel('pre-DART  running corr')
ylim([-.5 .5])
xlim([-.5 .5])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('-HTP')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off


subplot(1,2,2)
scatter(wheel_corr_concat{pre}(red_ind_concat),wheel_corr_concat{post}(red_ind_concat),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
% ylabel('post-DART dF/F')
xlabel('pre-DART  running corr')
ylim([-.5 .5])
xlim([-.5 .5])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
title('+HTP ')
axis square
hold off

sgtitle('Correlation of F timecourse with full wheel timecourse')
print(fullfile(fnout,['locCorr_crossDay.pdf']),'-dpdf','-bestfit')
%% 
figure
h1=cdfplot(wheel_corr_concat{pre}(red_ind_concat))
hold on
h2=cdfplot(wheel_corr_concat{pre}(green_ind_concat))
hold on
h3=cdfplot(wheel_corr_concat{post}(red_ind_concat))
hold on
h4=cdfplot(wheel_corr_concat{post}(green_ind_concat))
title("Locomotion correlation")
xlabel('Correlation of F timecourse with full wheel timecourse')
set(h1, 'LineStyle', '-', 'Color', 'k');
set(h2, 'LineStyle', '--', 'Color', 'k');
set(h3, 'LineStyle', '-', 'Color', 'b');
set(h4, 'LineStyle', '--', 'Color', 'b');
legend('+HTP','-HTP','Location','best')
hold off
print(fullfile(fnout,['locCorr_cdf.pdf']),'-dpdf','-bestfit')

%% scatterplot of average F values by day

max_lim_green=max([max(max(meanF_concat{post}(green_ind_concat))), max(max(meanF_concat{pre}(green_ind_concat)))]);
max_lim_red=max([max(max(meanF_concat{post}(red_ind_concat))), max(max(meanF_concat{pre}(red_ind_concat)))]);
min_lim=min([min(min(meanF_concat{post})), min(min(meanF_concat{pre}))]);

figure;
subplot(1,2,1);
scatter(meanF_concat{pre}(green_ind_concat),meanF_concat{post}(green_ind_concat),10,'MarkerEdgeColor',[.5 .5 .5])
title('-HTP')
xlim([min_lim max_lim_green]);
ylim([min_lim max_lim_green]);
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
xlabel('pre-DART')
ylabel('post-DART')
axis square

subplot(1,2,2);
scatter(meanF_concat{pre}(red_ind_concat),meanF_concat{post}(red_ind_concat),10,'MarkerEdgeColor',[.5 .5 .5])
title('+HTP')
xlim([min_lim max_lim_red]);
ylim([min_lim max_lim_red]);
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
xlabel('pre-DART')
axis square
sgtitle('mean raw F over all frames')

print(fullfile(fnout,['rawF_compare.pdf']),'-dpdf','-bestfit')

%% make figure with se shaded, compared stationary and loc for each day

z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_red_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
ylim([-.05 .4])
subplot(1,2,1) %for the first day
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_loc{pre}(:,iCon),tc_green_se_loc{pre}(:,iCon),'m');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['-HTP post-DART',' n = ', num2str(length(haveRunning_green{iCon}))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
ylim([-.05 .4])
shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_loc{post}(:,iCon),tc_green_se_loc{post}(:,iCon),'m');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['-HTP post-DART',' n = ', num2str(length(haveRunning_green{iCon}))])
x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['-HTP, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_locVsStat_pyr_timecourses.pdf']),'-dpdf');
end 




for iCon = 1:nCon
figure
ylim([-.05 .4])
subplot(1,2,1) %for the first day
shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon),tc_red_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_loc{pre}(:,iCon),tc_red_se_loc{pre}(:,iCon),'m');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['+HTP post-DART',' n = ', num2str(length(haveRunning_red{iCon}))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
ylim([-.05 .4])
shadedErrorBar(t,tc_red_avrg_stat{post}(:,iCon),tc_red_se_stat{post}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_loc{post}(:,iCon),tc_red_se_loc{post}(:,iCon),'m');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['+HTP post-DART',' n = ', num2str(length(haveRunning_red{iCon}))])
x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['+HTP, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_locVsStat_HT_timecourses.pdf']),'-dpdf');
end 
%% subtracted timecourses


tc_subtracted_stat = tc_trial_avrg_stat_concat{pre}-tc_trial_avrg_stat_concat{post};

    for iCon=1:nCon
        
    sub_tc_green_avrg_stat(:,iCon)=nanmean(tc_subtracted_stat(:,haveRunning_green{iCon},iCon),2);
    sub_green_std=nanstd(tc_subtracted_stat(:,haveRunning_green{iCon},iCon),[],2);
    sub_tc_green_se_stat(:,iCon)=sub_green_std/sqrt(length(haveRunning_green{iCon}));
    
    sub_tc_red_avrg_stat(:,iCon)=nanmean(tc_subtracted_stat(:,haveRunning_red{iCon},iCon),2);
    sub_red_std=nanstd(tc_subtracted_stat(:,haveRunning_red{iCon},iCon),[],2);
    sub_tc_red_se_stat(:,iCon)=sub_red_std/sqrt(length(haveRunning_red{iCon}));
    
    clear green_std red_std
    end

%%
z=double(nOn)/double(frame_rate);

%creat a time axis in seconds
t=1:(size(sub_tc_green_avrg_stat,1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure;
subplot(1,2,1) %for the first day

%ylim([-.02 .5]);
hold on
shadedErrorBar(t,sub_tc_green_avrg_stat(:,iCon),sub_tc_green_se_stat(:,iCon),'k');
hold on

line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['-HTP',' n = ', num2str(length(haveRunning_green{iCon}))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,sub_tc_red_avrg_stat(:,iCon),sub_tc_red_se_stat(:,iCon));
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['+HTP',' n = ', num2str(length(haveRunning_red{iCon}))])
x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

%sgtitle(['stationary, contrast = ' num2str(cons(iCon))])
hold on
%print(fullfile(fnout,[num2str(cons(iCon)) '_stat_cellType_timecourses.pdf']),'-dpdf');
end 
%% comparing Rsq to full tc corre

figure;scatter(mean_green_concat{pre},mean_green_concat{post})
xlabel("full TC correlation to average Pyr TC, pre")
ylabel("full TC correlation to average Pyr TC, post")

figure;scatter(mean_green_concat{pre},R_values_concat(1,:))
xlabel("full TC correlation to average Pyr TC, pre")
ylabel("Corr to average Pyr stim response")

%% timecourses for locomotion - stationary on a cell-by-cell average basis
tc_green_avrg_subtracted = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_subtracted = cell(1,nd); %same for red
tc_green_se_subtracted = cell(1,nd); %this will be the se across all green cells
tc_red_se_subtracted = cell(1,nd); %same for red

tc_subtracted = cell(1,nd);

for id = 1:nd
    for iCon=1:nCon
    tc_subtracted{id}(:,:,iCon)  = tc_trial_avrg_loc_concat{id}(:,:,iCon) - tc_trial_avrg_stat_concat{id}(:,:,iCon);
    tc_green_avrg_subtracted{id}(:,iCon)=nanmean(tc_subtracted{id}(:,haveRunning_green{iCon},iCon),2);
    green_std=nanstd(tc_subtracted{id}(:,haveRunning_green{iCon},iCon),[],2);
    tc_green_se_subtracted{id}(:,iCon)=green_std/sqrt(length(haveRunning_green{iCon}));
    
    tc_red_avrg_subtracted{id}(:,iCon)=nanmean(tc_subtracted{id}(:,haveRunning_red{iCon},iCon),2);
    red_std=nanstd(tc_subtracted{id}(:,haveRunning_red{iCon},iCon),[],2);
    tc_red_se_subtracted{id}(:,iCon)=red_std/sqrt(length(haveRunning_red{iCon}));
    
    clear green_std red_std
    end
end
 
%plotting
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_subtracted{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day



ylim([-.05 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_subtracted{pre}(:,iCon),tc_green_se_subtracted{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_subtracted{post}(:,iCon),tc_green_se_subtracted{post}(:,iCon),'b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['-HTP',' n = ', num2str(length(haveRunning_green{iCon}))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_subtracted{pre}(:,iCon),tc_red_se_subtracted{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_subtracted{post}(:,iCon),tc_red_se_subtracted{post}(:,iCon),'b');
ylim([-.05 .3]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['+HTP',' n = ', num2str(length(haveRunning_red{iCon}))])

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['loc - stat, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_subtracted_cellType_timecourses.pdf']),'-dpdf');
end 

%% subtracted mean responses

pref_resp_subtracted = cell(1,nd);

for id = 1:nd
    for iCon=1:nCon
        pref_resp_subtracted{id}(:,iCon)  = pref_responses_loc_concat{id}(:,iCon) - pref_responses_stat_concat{id}(:,iCon);
    end
end

 P = polyfit((pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon)),(pref_resp_subtracted{post}(haveRunning_green{iCon},iCon)),1);
 slope = P(1)
 intercept = P(2)
 [rho, pvalue] = corr((pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon)),(pref_resp_subtracted{post}(haveRunning_green{iCon},iCon)))

%% subtracted scatters


for iCon = 1:nCon
figure; movegui('center') 
subplot(1,2,1)
scatter((pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon)),(pref_resp_subtracted{post}(haveRunning_green{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
%hline=lsline;
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
mean_pre = nanmean(pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon));
mean_post = nanmean(pref_resp_subtracted{post}(haveRunning_green{iCon},iCon));
scatter(mean_pre,mean_post,20,'r','filled');
stderror_pre= std(pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
stderror_post= std(pref_resp_subtracted{post}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
line([(mean_pre-stderror_pre),(mean_pre+stderror_pre)],[mean_post, mean_post],'Color','black','LineWidth',2);
line([mean_pre, mean_pre],[(mean_post-stderror_post),(mean_post+stderror_post)],'Color','blue','LineWidth',2);
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
limMin=min(min(pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon)),min(pref_resp_subtracted{post}(haveRunning_green{iCon},iCon)));
limMax=max(max(pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon)),max(pref_resp_subtracted{post}(haveRunning_green{iCon},iCon)));
ylim([limMin limMax])
xlim([limMin limMax])

title('-HTP loc - stat')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off


subplot(1,2,2)
scatter((pref_resp_subtracted{pre}(haveRunning_red{iCon},iCon)),(pref_resp_subtracted{post}(haveRunning_red{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
%hline=lsline;
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
mean_pre = nanmean(pref_resp_subtracted{pre}(haveRunning_red{iCon},iCon));
mean_post = nanmean(pref_resp_subtracted{post}(haveRunning_red{iCon},iCon));
scatter(mean_pre,mean_post,20,'r','filled');
stderror_pre= std(pref_resp_subtracted{pre}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
stderror_post= std(pref_resp_subtracted{post}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
line([(mean_pre-stderror_pre),(mean_pre+stderror_pre)],[mean_post, mean_post],'Color','black','LineWidth',2);
line([mean_pre, mean_pre],[(mean_post-stderror_post),(mean_post+stderror_post)],'Color','blue','LineWidth',2);% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
limMin=min(min(pref_resp_subtracted{pre}(haveRunning_red{iCon},iCon)),min(pref_resp_subtracted{post}(haveRunning_red{iCon},iCon)));
limMax=max(max(pref_resp_subtracted{pre}(haveRunning_red{iCon},iCon)),max(pref_resp_subtracted{post}(haveRunning_red{iCon},iCon)));
ylim([limMin limMax])
xlim([limMin limMax])
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
title('+HTP loc-stat')
axis square
hold off




sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'subtracted_resp_crossDay.pdf']),'-dpdf','-bestfit')

% edges = [limMin:((limMax-limMin)/20):limMax];
% figure;
% subplot(1,2,1)
% hist(pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon),edges)
% xlim([limMin limMax])
% subplot(1,2,2)
% hist(pref_resp_subtracted{post}(haveRunning_green{iCon},iCon),edges)
% xlim([limMin limMax])
% sgtitle(num2str(cons(iCon)))
% print(fullfile(fnout,[num2str(cons(iCon)) 'subtracted_distribution.pdf']),'-dpdf','-bestfit')


clear mean_pre mean_post stderror_post stderror_pre
end

%%
[h1, p1]= ttest(pref_resp_subtracted{pre}(haveRunning_green{1},1),pref_resp_subtracted{post}(haveRunning_green{1},1));
[h2,p2]= ttest(pref_resp_subtracted{pre}(haveRunning_red{1},1),pref_resp_subtracted{post}(haveRunning_red{1},1));

%correct for four tests
p1*2
p2*2


clear h1 p1 h2 p2 

%% barcahrt
y = [mean(pref_resp_subtracted{pre}(haveRunning_green{1},1)), mean(pref_resp_subtracted{post}(haveRunning_green{1},1))]
x=[1:2];
std_pre = std(pref_resp_subtracted{pre}(haveRunning_green{1},1),[],1);
std_post = std(pref_resp_subtracted{post}(haveRunning_green{1},1),[],1);
se_pre = std_pre/sqrt(length(haveRunning_green{1}));
se_post = std_post/sqrt(length(haveRunning_green{1}));

errhigh = [y(1)+se_pre y(2)+se_post];
errlow = [y(1)-se_pre y(2)-se_post]; 

bar(x,y)

hold on

%er = errorbar(x,y,errlow,errhigh);

%% plot contrast response

conResp_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_stat = cell(1,nd); %same for red
conResp_green_se_stat = cell(1,nd); %this will be the se across all green cells
conResp_red_se_stat = cell(1,nd); %same for red



for id = 1:nd
   
        
    conResp_green_avrg_stat{id}=nanmean(pref_responses_stat_concat{id}(green_ind_concat ,:),1);
    green_std=nanstd(pref_responses_stat_concat{id}(green_ind_concat,:),1);
    conResp_green_se_stat{id}=green_std/sqrt(length(green_ind_concat));
    
    conResp_red_avrg_stat{id}=nanmean(pref_responses_stat_concat{id}(red_ind_concat,:),1);
    red_std=nanstd(pref_responses_stat_concat{id}(red_ind_concat,:),1);
    conResp_red_se_stat{id}=green_std/sqrt(length(red_ind_concat));
    
    clear green_std red_std
 
end


figure
subplot(1,2,1) %for the first day
errorbar(cons,conResp_green_avrg_stat{pre},conResp_green_se_stat{pre},'k');
hold on
errorbar(cons,conResp_green_avrg_stat{post},conResp_green_se_stat{post},'b');
title(['-HTP',' n = ', num2str(length(green_ind_concat))])
ylabel('dF/F, pref ori') 
xlabel('contrast') 
set(gca, 'TickDir', 'out')
box off

subplot(1,2,2) %for the second day
errorbar(cons,conResp_red_avrg_stat{pre},conResp_red_se_stat{pre},'k');
hold on
errorbar(cons,conResp_red_avrg_stat{post},conResp_red_se_stat{post},'b');
title(['+HTP',' n = ', num2str(length(red_ind_concat))])
ylabel('dF/F, pref ori') 
xlabel('contrast') 
set(gca, 'TickDir', 'out')
box off

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])



sgtitle(['stationary, all keep cells' ])

print(fullfile(fnout,['contrast_resposnse.pdf']),'-dpdf');



