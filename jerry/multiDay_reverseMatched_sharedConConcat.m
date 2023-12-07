
clear all; clear global; close all
clc
ds = 'DART_expt_info'; %dataset info
dataStructLabels = {'contrastxori'};

eval(ds);


sess_list = [10];%enter all the sessions you want to concatenate; ONLY ENTER THE REFERENCE DATE IF MATCHING
nSess=length(sess_list);

nd=2;%hard coding for two days per experimental session

% INDICATE THE PRE VS. POST DAYS DEPENDING ON THE ORDER OF MATCHING
prompt = 'Which sesson was used as reference for matching: 0- baseline, 1- post-DART';
            x = input(prompt);
            switch x
                case 0
                    pre=1; %baseline session, used as reference, is in the 1st position
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

if computer == 'GLNXA64'
    isilonName =  '/home/cc735@dhe.duke.edu/GlickfeldLabShare';
    base = fullfile(isilonName, '/All_Staff/home/ACh/Data/2p_data');
    
else
    isilonName = 'G:';
    base = fullfile(isilonName, '\home\ACh\Analysis\2p_analysis');
      
end


if nSess == 1
         if expt(sess_list(1)).multiday_timesincedrug_hours>0
            dart_str = [expt(sess_list(1)).drug '_' num2str(expt(sess_list(1)).multiday_timesincedrug_hours) 'Hr'];
        else
            dart_str = 'control';
        end
        
        fnout = fullfile(base,expt(sess_list(1)).mouse,['multiday_' dart_str],d);
else
    fnout= fullfile(base,strcat('concat', sess_title),d);
end
mkdir(fnout);
cd(fnout)
clear d sess_title

zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
%% concatenating data
nCon = length(targetCon)

mice={};
tc_trial_avrg_stat_concat=cell(1,nd);
tc_trial_avrg_loc_concat=cell(1,nd);
tc_trial_avrg_keep_allCond_concat=cell(1,nd);
resp_keep_concat=cell(1,nd);
resp_max_keep_concat=cell(1,nd);
pref_responses_loc_concat=cell(1,nd);
pref_responses_stat_concat=cell(1,nd);
pref_responses_allCond_concat=cell(1,nd);
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
wheel_corr_concat=cell(1,nd);
meanF_concat=cell(1,nd);
norm_dir_resp_stat_concat = cell(1,nd);
norm_dir_resp_loc_concat = cell(1,nd);
pref_dir_concat=cell(1,nd);
noiseCorr_concat = cell(1,nd);
sigCorr_concat = cell(1,nd);
pref_allTrials_stat_concat =cell(nCon,nd);
pref_allTrials_loc_concat =cell(nCon,nd);
dataTableConat=[];

for iSess = 1:nSess
    day_id = sess_list(iSess)
    mouse = expt(day_id).mouse;
    mice=[mice;mouse];

    if expt(day_id).multiday_timesincedrug_hours>0
        dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fn_multi = fullfile(base,mouse,['multiday_' dart_str]);

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

    clear cons
    
    
    for id = 1:nd
        
        tc_trial_avrg_stat_concat{id} =cat(2,tc_trial_avrg_stat_concat{id},tc_trial_avrg_stat{id}(:,:,sharedCon));
        tc_trial_avrg_loc_concat{id} =cat(2,tc_trial_avrg_loc_concat{id},tc_trial_avrg_loc{id}(:,:,sharedCon));
%        tc_trial_avrg_keep_allCond_concat{id} =cat(2,tc_trial_avrg_keep_allCond_concat{id},tc_trial_avrg_keep_allCond{id}(:,:,sharedCon));
        resp_keep_concat{id}=cat(1,resp_keep_concat{id},resp_keep{id});
        resp_max_keep_concat{id}=cat(1,resp_max_keep_concat{id},resp_max_keep{id}(:,sharedCon));
        LMI_concat{id}=cat(1,LMI_concat{id},LMI{id}(:,sharedCon));
        pref_responses_loc_concat{id}=cat(1,pref_responses_loc_concat{id},pref_responses_loc{id}(:,sharedCon));
        pref_responses_stat_concat{id}=cat(1,pref_responses_stat_concat{id},pref_responses_stat{id}(:,sharedCon));
        pref_responses_allCond_concat{id}=cat(1,pref_responses_allCond_concat{id},pref_responses_allCond{id}(:,sharedCon));
        RIx_concat{id}=cat(1,RIx_concat{id},sum(RIx{id}));
        wheel_corr_concat{id}=cat(2,wheel_corr_concat{id},wheel_corr{id});
        meanF=mean(fullTC_keep{id},1);
        meanF_concat{id}=cat(2,meanF_concat{id}, meanF);
        norm_dir_resp_stat_concat{id}=cat(1,norm_dir_resp_stat_concat{id},norm_dir_resp_stat{id});
        norm_dir_resp_loc_concat{id}=cat(1,norm_dir_resp_loc_concat{id},norm_dir_resp_loc{id});
        pref_dir_concat{id}=cat(2,pref_dir_concat{id},pref_dir_keep{id});
        noiseCorr_concat{id}=cat(2,noiseCorr_concat{id},noiseCorr{id});
        sigCorr_concat{id}=cat(2,sigCorr_concat{id},sigCorr{id});
        for i = 1:length(sharedCon)
            iCon=sharedCon(i);
            pref_allTrials_stat_concat{i,id}=[pref_allTrials_stat_concat{i,id},pref_allTrials_stat{iCon,id}];
            pref_allTrials_loc_concat{i,id}=[pref_allTrials_loc_concat{i,id},pref_allTrials_loc{iCon,id}];
        end
        clear meanF
    end
    dfof_max_diff_concat=cat(1,dfof_max_diff_concat,dfof_max_diff(:,sharedCon));
    green_fluor_concat=cat(2,green_fluor_concat,green_fluor_keep);
    red_fluor_concat=cat(2,red_fluor_concat,red_fluor_keep);
    

end
%
clear mouse day_id nKeep iSess fn_multi cons oris
clear explanation1 resp_keep tc_trial_avrg_keep_allCond pref_responses_allCond sig_diff pref_con_keep pref_ori_keep tOri_match tCon_match data_trial_keep nTrials tc_trial_avrg_keep green_keep_logical red_keep_logical green_ind_keep red_ind_keep
clear LMI RIx locCounts locResp locTCs statResp statTCs wheel_tc ttest_results_stat ttest_results_loc ttest_results_allCon_stat ttest_results_allCon_loc
clear data_con_resp_keep data_ori_resp_keep data_rep_keep dfof_max_diff dfof_max_diff_raw explanation2 resp_max_keep data_resp_keep pref_responses_stat pref_responses_loc
clear tc_trial_avrg_stat tc_trial_avrg_loc fullTC_keep norm_dir_resp sigCorr noiseCorr
clear red_fluor_all red_fluor_match green_fluor_match green_fluor_match red_fluor_keep green_fluor_keep R_p_values
red_ind_concat = find(red_concat);
green_ind_concat = find(green_concat);
%
cons=unique(cons_concat);

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

%finds the cells that have running and stationary for all three contrasts
green_all = intersect(haveRunning_green{1},haveRunning_green{2});
green_all = intersect(green_all, haveRunning_green{3});

red_all = intersect(haveRunning_red{1},haveRunning_red{2});
red_all = intersect(red_all, haveRunning_red{3});

% %for instances when I know I don't have enough running trials, just use
% all keep cells.
% 
% for iCon =1:nCon
% 
% haveRunning_green{iCon} = green_ind_concat;
% haveRunning_red{iCon} = red_ind_concat;
% end

%get an array with indices for each mouse
mouseInds=cell(1,nSess);
start=1;
for iMouse = 1:nSess
    mouseInds{iMouse}=start:(start-1)+nKeep_concat(iMouse);
    start = start+nKeep_concat(iMouse)
end
clear start iMouse


%find how many haveRunning red cells I got for each mouse
cellCountsRed = nan(nSess,nCon);
mouseNames=[];
for iMouse = 1:nSess
    for iCon = 1:nCon
        cellCountsRed(iMouse, iCon,1)=length(intersect(haveRunning_red{iCon},(mouseInds{iMouse})));
        
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


cellCountTableRed = table(cellCountsRed, RowNames=mouseNames)
cellCountTableGreen = table(cellCountsGreen, RowNames=mouseNames)
writetable(cellCountTableRed,fullfile(fnout,['cellCounts.csv']),'WriteRowNames',true)
writetable(cellCountTableGreen,fullfile(fnout,['cellCounts_Green.csv']),'WriteRowNames',true)
clear cellCountsRed cellCountsGreen
%% ALL KEEP CELLS - this is NOT using only the cells that have running trials, instead it's using all cells

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

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day



ylim([-.02 .25]);;
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
%fill([0.2 0.2 (z+.2) (z+.2)],[-.02 .3 .3 -.02],'k',FaceAlpha = 0.25,LineStyle='none')
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['-HTP',' n = ', num2str(length(green_ind_concat))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon),tc_red_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_stat{post}(:,iCon),tc_red_se_stat{post}(:,iCon),'b');
ylim([-.02 .25]);;
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
%fill([0.2 0.2 (z+.2) (z+.2)],[-.02 .3 .3 -.02],'k',FaceAlpha = 0.25,LineStyle='none')
title(['+HTP',' n = ', num2str(length(red_ind_concat))])

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['stationary, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_stat_allKeep_timecourses.pdf']),'-dpdf');
end ; 
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



ylim([-.02 .3]);;
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
%fill([0.2 0.2 (z+.2) (z+.2)],[-.02 .25 .25 -.02],'k',FaceAlpha = 0.25,LineStyle='none')
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
ylim([-.02 .3]);;
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
%fill([0.2 0.2 (z+.2) (z+.2)],[-.02 .25 .25 -.02],'k',FaceAlpha = 0.25,LineStyle='none')
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
shadedErrorBar(t,tc_green_avrg_loc{pre}(:,iCon),tc_green_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_loc{post}(:,iCon),tc_green_se_loc{post}(:,iCon),'b');
ylim([-.02 .5]);;
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
hold on
%fill([0.2 0.2 (z+.2) (z+.2)],[-.02 .45 .45 -.02],'k',FaceAlpha = 0.25,LineStyle='none')
title(['-HTP',' n = ', num2str(length(haveRunning_green{iCon}))])
ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_loc{pre}(:,iCon),tc_red_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_loc{post}(:,iCon),tc_red_se_loc{post}(:,iCon),'b');
ylim([-.02 .5]);;
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
hold on
%fill([0.2 0.2 (z+.2) (z+.2)],[-.02 .45 .45 -.02],'k',FaceAlpha = 0.25,LineStyle='none')
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
mean_pre_stat = nanmean(pref_responses_stat_concat{pre}(haveRunning_green{iCon},iCon));
mean_post_stat = nanmean(pref_responses_stat_concat{post}(haveRunning_green{iCon},iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
stderror_post= std(pref_responses_stat_concat{post}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);
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
mean_pre_stat = nanmean(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon));
mean_post_stat = nanmean(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
stderror_post= std(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);% ylabel('post-DART dF/F')
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
mean_pre_stat = nanmean(pref_responses_loc_concat{pre}(haveRunning_green{iCon},iCon));
mean_post_stat = nanmean(pref_responses_loc_concat{post}(haveRunning_green{iCon},iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_loc_concat{pre}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
stderror_post= std(pref_responses_loc_concat{post}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);
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
mean_pre_stat = nanmean(pref_responses_loc_concat{pre}(haveRunning_red{iCon},iCon));
mean_post_stat = nanmean(pref_responses_loc_concat{post}(haveRunning_red{iCon},iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
stderror_post= std(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);
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
clear mean_pre_stat mean_post_stat stderror_post stderror_pre
end



%%
responseTable = table([nanmean(pref_responses_stat_concat{pre}(green_all));nanmean(pref_responses_stat_concat{post}(green_all))],[nanmean(pref_responses_stat_concat{pre}(red_all));nanmean(pref_responses_stat_concat{post}(red_all))],[nanmean(pref_responses_loc_concat{pre}(green_all));nanmean(pref_responses_loc_concat{post}(green_all))],[nanmean(pref_responses_loc_concat{pre}(red_all));nanmean(pref_responses_loc_concat{post}(red_all))],'VariableNames',{'Pyramidal cells stat'  '+HTP cells stat' 'Pyramidal cells loc'  '+HTP cells loc'}, 'RowNames',{'Pre'  'Post'})
writetable(responseTable,fullfile(fnout,[num2str(targetCon) 'responseTable.csv']),'WriteRowNames',true)


%% linear direction plot 

dirs_for_plotting=dirs-180;

green_dir_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
red_dir_avrg_stat = cell(1,nd); %same for red
green_dir_se_stat = cell(1,nd); %this will be the se across all green cells
red_dir_se_stat = cell(1,nd); %same for red

for id = 1:nd
   % NOTE: add a cathc to only include cells that were responsive on each
   % day
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

sgtitle("Normalized direction tuning ")


print(fullfile(fnout,['dirTuning.pdf']),'-dpdf','-bestfit')

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


%% plot contrast response

conResp_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_stat = cell(1,nd); %same for red
conResp_green_se_stat = cell(1,nd); %this will be the se across all green cells
conResp_red_se_stat = cell(1,nd); %same for red



for id = 1:nd
   
        
    conResp_green_avrg_stat{id}=nanmean(pref_responses_stat_concat{id}(green_all,:),1);
    green_std=nanstd(pref_responses_stat_concat{id}(green_all,:),1);
    conResp_green_se_stat{id}=green_std/sqrt(length(green_all));
    
    conResp_red_avrg_stat{id}=nanmean(pref_responses_stat_concat{id}(red_all,:),1);
    red_std=nanstd(pref_responses_stat_concat{id}(red_all,:),1);
    conResp_red_se_stat{id}=red_std/sqrt(length(red_all));
    
    clear green_std red_std
 
end


figure
subplot(1,2,1) %for the first day
errorbar(cons,conResp_green_avrg_stat{pre},conResp_green_se_stat{pre},'k');
hold on
errorbar(cons,conResp_green_avrg_stat{post},conResp_green_se_stat{post},'b');
title(['-HTP',' n = ', num2str(length(green_all))])
ylabel('dF/F, pref ori') 
xlabel('contrast') 
set(gca, 'TickDir', 'out')
xlim([0 1])
ylim([0 .3])
xticks([.25 .5 1])
box off

subplot(1,2,2) %for the second day
errorbar(cons,conResp_red_avrg_stat{pre},conResp_red_se_stat{pre},'k');
hold on
errorbar(cons,conResp_red_avrg_stat{post},conResp_red_se_stat{post},'b');
title(['+HTP',' n = ', num2str(length(red_all))])
xlabel('contrast') 
xlim([0 1])
ylim([0 .3])
xticks([.25 .5 1])
set(gca, 'TickDir', 'out')
box off

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])



sgtitle(['stationary' ])

print(fullfile(fnout,['contrast_response.pdf']),'-dpdf');



conResp_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_loc = cell(1,nd); %same for red
conResp_green_se_loc = cell(1,nd); %this will be the se across all green cells
conResp_red_se_loc = cell(1,nd); %same for red



for id = 1:nd
   
        
    conResp_green_avrg_loc{id}=nanmean(pref_responses_loc_concat{id}(green_all ,:),1);
    green_std=nanstd(pref_responses_loc_concat{id}(green_all,:),1);
    conResp_green_se_loc{id}=green_std/sqrt(length(green_all));
    
    conResp_red_avrg_loc{id}=nanmean(pref_responses_loc_concat{id}(red_all,:),1);
    red_std=nanstd(pref_responses_loc_concat{id}(red_all,:),1);
    conResp_red_se_loc{id}=red_std/sqrt(length(red_all));
    
    clear green_std red_std
 
end


figure
subplot(1,2,1) %for the first day
errorbar(cons,conResp_green_avrg_loc{pre},conResp_green_se_loc{pre},'k');
hold on
errorbar(cons,conResp_green_avrg_loc{post},conResp_green_se_loc{post},'b');
title(['-HTP',' n = ', num2str(length(green_all))])
ylabel('dF/F, pref ori') 
xlabel('contrast') 
xlim([0 1])
ylim([0 .5])
xticks([.25 .5 1])
set(gca, 'TickDir', 'out')
box off

subplot(1,2,2) 
errorbar(cons,conResp_red_avrg_loc{pre},conResp_red_se_loc{pre},'k');
hold on
errorbar(cons,conResp_red_avrg_loc{post},conResp_red_se_loc{post},'b');
title(['+HTP',' n = ', num2str(length(red_all))])
ylim([0 .6])
xlim([0 1])
xticks([.25 .5 1])
xlabel('contrast') 
set(gca, 'TickDir', 'out')
box off

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])



sgtitle(['Running' ])

print(fullfile(fnout,['loc_contrast_resposnse.pdf']),'-dpdf');

%% response by condition



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
hold on
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
