
clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};

eval(ds);
% 259 269 emx atropine - forward matched
% 206 210 214 atropine SOM
% 201 197 emx YM90K - forward matched
% 178 190 294 %good quality SOM YM90K
%138 142 163 171 178 190 294 307 for retreat talk
%294 307 323 NES with DART
%peg 303 311 319 329
%Oct 2023 138 142 163 171 178 190 294 307 323 333

sess_list = [345];%enter all the sessions you want to concatenate
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

if computer == 'GLNXA64'
    isilonName =  '/home/cc735@dhe.duke.edu/GlickfeldLabShare';
    base = fullfile(isilonName, '/All_Staff/home/ACh/Data/2p_data');
    
else
    isilonName = 'Z:';
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
%    temp_table =readtable(fullfile(fn_multi,'dataTable.csv'));
    
    % temp_table.z_speed=zscor_xnan(temp_table.speed);
    % temp_table.z_pupil=zscor_xnan(temp_table.pupil);
    % 
    % dataTableConat=[dataTableConat; temp_table];

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



% %%
% % loop to add "b" to the end of mice IDs where I have more than one session
% % with that mouse
% %set z to be the order position of first session that is at the earlier timepoint
% z = 12;
% for iMouse = z:nSess
%     mice{iMouse}=[mice{iMouse},'b'];
% end
%% mixed model

dataTableConat.day=categorical(dataTableConat.day);
dataTableConat.mouseID=categorical(dataTableConat.mouseID);
dataTableConat.cellID=categorical(dataTableConat.cellID);
dataTableConat.cellType=categorical(dataTableConat.cellType);
dataTableConat.drug = categorical(dataTableConat.drug);
dataTableConat = sortrows(dataTableConat,"day","ascend");

%lme = fitlme(dataTableConat,'dfof_perCent~cellType+contrast+z_speed+z_pupil+day+(cellType*day)+(1|mouseID)+(1|cellID:mouseID)');
lme = fitlme(dataTableConat,'dfof~cellType+contrast+z_speed+z_pupil+day+drug+(cellType*day*drug)+(1|mouseID)+(1|cellID:mouseID)+(1|day:cellID)+(1|drug:cellID)');

%% post-hoc comaparisons for mixed model
%what is the mean df/f for each cell type on each day
vars = ["dfof"];
factors = ["drug","cellType","day"];
meanScoresByFactor = varfun(@mean, ...
                            dataTableConat, ...
                            "InputVariables",vars, ...
                            "GroupingVariables",factors)

%% find cells that I have running data for on both days
haveRunning_pre = ~isnan(pref_responses_loc_concat{pre});

haveRunning_post = ~isnan(pref_responses_loc_concat{post});
for iCon =1:nCon
haveRunning_both{iCon}= find(haveRunning_pre(:,iCon).* haveRunning_post(:,iCon));
haveRunning_green{iCon} = intersect(haveRunning_both{iCon}, green_ind_concat);
haveRunning_red{iCon} = intersect(haveRunning_both{iCon}, red_ind_concat);
end

clear haveRunning_pre haveRunning_post haveRunning_both


% %for instances when I know I don't have enough running trials, juse use
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
%% ALL KEEP CELLS

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


%% t-test for the lowest contrast

[h1, p1]= ttest(pref_responses_stat_concat{pre}(haveRunning_green{1},1),pref_responses_stat_concat{post}(haveRunning_green{1},1));
[h2,p2]= ttest(pref_responses_stat_concat{pre}(haveRunning_red{1},1),pref_responses_stat_concat{post}(haveRunning_red{1},1));
[h3,p3]= ttest(pref_responses_loc_concat{pre}(haveRunning_green{1},1),pref_responses_loc_concat{post}(haveRunning_green{1},1));
[h4,p4]= ttest(pref_responses_loc_concat{pre}(haveRunning_red{1},1),pref_responses_loc_concat{post}(haveRunning_red{1},1));

%correct for four tests
p1*4
p2*4
p3*4
p4*4


%clear h1 p1 h2 p2 h3 p3 h4 p4
%%
responseTable = table([nanmean(pref_responses_stat_concat{pre}(haveRunning_green));nanmean(pref_responses_stat_concat{post}(haveRunning_green))],[nanmean(pref_responses_stat_concat{pre}(haveRunning_red));nanmean(pref_responses_stat_concat{post}(haveRunning_red))],[nanmean(pref_responses_loc_concat{pre}(haveRunning_green));nanmean(pref_responses_loc_concat{post}(haveRunning_green))],[nanmean(pref_responses_loc_concat{pre}(haveRunning_red));nanmean(pref_responses_loc_concat{post}(haveRunning_red))],'VariableNames',{'Pyramidal cells stat'  '+HTP cells stat' 'Pyramidal cells loc'  '+HTP cells loc'}, 'RowNames',{'Pre'  'Post'})
writetable(responseTable,fullfile(fnout,[num2str(targetCon) 'responseTable.csv']),'WriteRowNames',true)


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
raw_diff_stat(:,iCon)=((pref_responses_stat_concat{post}(:,iCon))-(pref_responses_stat_concat{pre}(:,iCon)));
raw_diff_loc(:,iCon)=((pref_responses_loc_concat{post}(:,iCon))-(pref_responses_loc_concat{pre}(:,iCon)));



fract_diff_stat(:,iCon)=((pref_responses_stat_concat{post}(:,iCon))-(pref_responses_stat_concat{pre}(:,iCon)))./(pref_responses_stat_concat{pre}(:,iCon));
fract_diff_loc(:,iCon)=((pref_responses_loc_concat{post}(:,iCon))-(pref_responses_loc_concat{pre}(:,iCon)))./(pref_responses_loc_concat{pre}(:,iCon));


figure
x1 = [mean(fract_diff_stat(haveRunning_red{iCon},iCon)),mean(fract_diff_loc(haveRunning_red{iCon},iCon))];
y1 = [(std(fract_diff_stat(haveRunning_red{iCon},iCon)))/sqrt(length(haveRunning_red{iCon})),(std(fract_diff_loc(haveRunning_red{iCon},iCon)))/sqrt(length(haveRunning_red{iCon}))];

x2 = [mean(fract_diff_stat(haveRunning_green{iCon},iCon)), mean(fract_diff_loc(haveRunning_green{iCon},iCon))];
y2 = [(std(fract_diff_stat(haveRunning_green{iCon},iCon)))/sqrt(length(haveRunning_green{iCon})),(std(fract_diff_loc(haveRunning_green{iCon},iCon)))/sqrt(length(haveRunning_green{iCon}))];

labs =categorical({'stat','loc'});

subplot(1,2,1)
bar(labs,x1)                
hold on
er = errorbar(labs,x1,-y1,y1);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title("+HTP")
%ylim([-.05 .05])
ylabel('Mean')
subplot(1,2,2)
bar(labs,x2)                
hold on
er = errorbar(labs,x2,-y2,y2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
%ylim([-.05 .05])
title("-HTP")
hold off
%sgtitle(['(Post - Pre) dFoF, contrast = ', num2str(cons(iCon))])
sgtitle(['(Post - Pre)/Pre dFoF, contrast = ', num2str(cons(iCon))])
print(fullfile(fnout,[num2str(cons(iCon)), '_fract_change_resp.pdf']),'-dpdf','-bestfit')

end


diff_values=nan(nCon,2);
for iCon = 1:nCon
    diff_values(iCon,1)=mean(fract_diff_stat(red_ind_concat,iCon),1);
    diff_values(iCon,2)=(std(fract_diff_stat(red_ind_concat,iCon))/sqrt(length(red_ind_concat)));
end


 %% bar chart of fraction red cells significantly suppressed or faciltiated
norm_diff = nan(2,nCon,nKeep_total);
for i = 1:nKeep_total
    for iCon = 1:nCon
        %for stationary trials
        mean_pre_stat = mean(pref_allTrials_stat_concat{iCon,pre}{i});
        mean_post_stat=mean(pref_allTrials_stat_concat{iCon,post}{i});
        std_pre_stat = std(pref_allTrials_stat_concat{iCon,pre}{i});
        norm_diff_stat = (mean_post_stat-mean_pre_stat) / std_pre_stat;

        %for running trials
        mean_pre_loc = mean(pref_allTrials_loc_concat{iCon,pre}{i});
        mean_post_loc=mean(pref_allTrials_loc_concat{iCon,post}{i});
        std_pre_loc = std(pref_allTrials_loc_concat{iCon,pre}{i});
        norm_diff_loc = (mean_post_loc-mean_pre_loc)/ std_pre_loc;

        %putting data into matrix
        norm_diff(1,iCon,i)=norm_diff_stat; %first is stationary
        norm_diff(2,iCon,i)=norm_diff_loc; %second is running
clear mean_pre_stat mean_post_stat std_pre_stat mean_pre_loc mean_post_loc std_pre_loc 
    end 
end
%remove any infiinty values resulting from divisions by zero, and turn
%those into NANs instead
norm_diff(find(norm_diff == -Inf))=nan;
norm_diff(find(norm_diff == Inf))=nan;

%find how many cells are suppressed ( normalized diff < -1) or facilitated
%(normalized diff >1) at each contrast and behavioral state

suppressed = logical(norm_diff < -1);
facilitated = logical(norm_diff > 1);

%pull out red cells, get fractions that are suppressed vs. favilitated,
%stationary only, using all keep cells (i.e., not only the ones that have
% %running trials
fractSupp_stat = sum((squeeze(suppressed(1,:,:))),2)/nKeep_total;
fractFacil_stat = sum((squeeze(facilitated(1,:,:))),2)/nKeep_total;

% fractSupp_stat = sum((squeeze(suppressed(1,:,red_ind_concat))),2)/length(red_ind_concat);
% fractFacil_stat = sum((squeeze(facilitated(1,:,red_ind_concat))),2)/length(red_ind_concat);


% fractSupp_stat_LowR = sum((squeeze(suppressed(1,:,redLow_acrossCons))),2)/length(redLow_acrossCons);
% fractFacil_stat_LowR = sum((squeeze(facilitated(1,:,redLow_acrossCons))),2)/length(redLow_acrossCons);
% 
% fractSupp_allCond = sum(squeeze(mean(suppressed(:,:,red_ind_concat),1)),2)/length(red_ind_concat);
% fractFacil_allCond = sum(squeeze(mean(facilitated(:,:,red_ind_concat),1)),2)/length(red_ind_concat);

figure;
bar(cons,[fractSupp_stat fractFacil_stat],'stacked')
xticks([.25 .5 1])
ylabel(["Fraction HTP+ cells"]) 
xlabel('Contrast')
set(gca,'TickDir','out')
box off
title('Stationary, ')

%pull out red cells, get fractions that are suppressed vs. favilitated, for
%the red cells that have running trials at all contrasts
fractSupp_red = nan(2,3);
fractFacil_red = nan(2,3);
for iCon = 1:nCon
    fractSupp_red(:,iCon) = sum(((suppressed(:,iCon,haveRunning_red{iCon}))),3)/length(haveRunning_red{iCon});
    fractFacil_red(:,iCon) = sum(((facilitated(:,iCon,haveRunning_red{iCon}))),3)/length(haveRunning_red{iCon});
end


%averaging over contrast
supp_red = mean(fractSupp_red,2)
facil_red = mean(fractFacil_red,2)


figure;

bar([1,2],[supp_red(1) facil_red(1);supp_red(2) facil_red(2)],'stacked')
xticklabels({'Stationary','Running'})
ylabel(["Fraction HTP+ cells"]) 
set(gca,'TickDir','out')
box off
%ylim([0 .7])

norm_diff_red_allKeep = nanmedian(norm_diff(1,:,red_ind_concat),3)



%%
figure;
subplot(1,2,1)
boxchart((squeeze(mean(norm_diff(:,:,red_all),2,"omitnan")))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75]);hold on;
scatter([1 2],(squeeze(mean(norm_diff(:,:,red_all),2,"omitnan")))',"red",'jitter', 'on', 'jitterAmount',.1)
%scatter([1 2],(squeeze(mean(norm_diff(:,:,green_all),2,"omitnan")))',"black",'jitter', 'on', 'jitterAmount',.1)
%plot(squeeze(mean(norm_diff(:,:,red_all),2,"omitnan")),'--k')
ylim([-6 10])
xticklabels({'Stationary','Running'})
ylabel('(Post-Pre) / std dev Pre')
title('+HTP')
hold off
set(gca,'TickDir','out')
box off

subplot(1,2,2)
boxchart((squeeze(mean(norm_diff(:,:,green_all),2,"omitnan")))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75]);hold on;
scatter([1 2],(squeeze(mean(norm_diff(:,:,green_all),2,"omitnan")))',"black",'jitter', 'on', 'jitterAmount',.1)
%plot(squeeze(mean(norm_diff(:,:,red_all),2,"omitnan")),'--k')
ylim([-6 10])
xticklabels({'Stationary','Running'})
ylabel('(Post-Pre) / std dev Pre')
title('-HTP')
hold off
set(gca,'TickDir','out')
box off

%% 
redHigh_acrossCons = intersect(highRInds,red_ind_concat);
redHigh_acrossCons = intersect(lowRInds,red_ind_concat);

figure;
boxchart((squeeze(mean(norm_diff(:,2,haveRunning_red{2}),2,"omitnan")))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75]);hold on;
scatter([1 2],(squeeze(mean(norm_diff(:,2,redHigh_acrossCons),2,"omitnan")))',"black",'jitter', 'on', 'jitterAmount',.1)
scatter([1 2],(squeeze(mean(norm_diff(:,2,redHigh_acrossCons),2,"omitnan")))',"red",'jitter', 'on', 'jitterAmount',.1)
ylim([-8 8])
xticklabels({'Stationary','Running'})
ylabel('(Post-Pre) / std dev Pre')
title('Normalized change at 50% contrast')
hold off
set(gca,'TickDir','out')
box off

%%
figure;
subplot(1,2,1)
boxchart(squeeze(norm_diff(1,:,red_ind_concat))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75]);
hold on
scatter([1, 2, 3],squeeze(norm_diff(1,:,red_ind_concat))',"red",'jitter', 'on', 'jitterAmount',.1)
%scatter([1, 2, 3],squeeze(norm_diff(1,:,green_ind_concat))',"black",'jitter', 'on', 'jitterAmount',.1)
%plot(squeeze(norm_diff(1,[1,3],red_ind_concat)))
xticklabels({'25%','50%','100%'})
ylim([-6 10])
xlabel('Contrast')
ylabel('(Post-Pre) / std dev Pre')
title('+HTP')
hold off
set(gca,'TickDir','out')
box off

subplot(1,2,2)
boxchart(squeeze(norm_diff(1,:,green_ind_concat))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75]);
hold on
scatter([1, 2, 3],squeeze(norm_diff(1,:,green_ind_concat))',"black",'jitter', 'on', 'jitterAmount',.1)
%plot(squeeze(norm_diff(1,[1,3],red_ind_concat)))
xticklabels({'25%','50%','100%'})
ylim([-6 10])
xlabel('Contrast')
ylabel('(Post-Pre) / std dev Pre')
title('-HTP')
hold off
set(gca,'TickDir','out')
box off

%%
z=squeeze(norm_diff(1,:,red_ind_concat))';
suppFacilColors = nan(size(z));

suppFacilColors(suppFacilColors<1)=[0 0.4470 0.7410];



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
 
       end
    tHalfMax{id}=tHalfMaxTemp;
end
clear tHalfMaxCell tHalfMaxTemp tempData smoothData halfMax
% scatter for tMax



figure; movegui('center') 
subplot(1,2,1)
scatter((tHalfMax{pre}(green_all)),(tHalfMax{post}(green_all)),10,'MarkerEdgeColor',[.4 .4 .4],'jitter', 'on', 'jitterAmount',.1)
hold on
mean_pre_stat = mean(tHalfMax{pre}(green_all),"omitnan");
mean_post_stat = mean(tHalfMax{post}(green_all),"omitnan");
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= nanstd(tHalfMax{pre}(green_all)) / sqrt( length(green_all));
stderror_post= nanstd(tHalfMax{post}(green_all)) / sqrt( length(green_all));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);
ylabel('post-DART half-max(s)')
xlabel('pre-DART half-max(s)')
ylim([0 2.15])
xlim([0 2.15])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('-HTP ')
axis square


hold off



subplot(1,2,2)
scatter((tHalfMax{pre}(red_all)),(tHalfMax{post}(red_all)),10,'MarkerEdgeColor',[.4 .4 .4],'jitter', 'on', 'jitterAmount',.1)
hold on
mean_pre_stat = nanmean(tHalfMax{pre}(red_all));
mean_post_stat = nanmean(tHalfMax{post}(red_all));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= nanstd(tHalfMax{pre}(red_all)) / sqrt( length(red_all));
stderror_post= nanstd(tHalfMax{post}(red_all)) / sqrt( length(red_all));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);

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

% [h1, p1]= ttest(tHalfMax{pre}(green_all),tHalfMax{post}(green_all));
% [h2,p2]= ttest(tHalfMax{pre}(red_all),tHalfMax{post}(red_all));

%correct for four tests
% p1*2
% p2*2
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
%% compare r value to LMI
meanLMI = mean(LMI_concat{pre},2);

figure

subplot(1,2,1)
scatter(noiseCorr_concat{pre}(1,green_ind_concat),meanLMI(green_ind_concat),'MarkerEdgeColor','#47E5BC', 'LineWidth',1.25)
hold on
scatter(noiseCorr_concat{pre}(1,red_ind_concat),meanLMI(red_ind_concat),'MarkerEdgeColor','#93748A', 'LineWidth',1.25)
xlabel('Noise corr pre-DART')
ylabel('LMI pre-DART')
axis square
set(gca, 'TickDir', 'out')

subplot(1,2,2)
scatter(noiseCorr_concat{pre}(1,green_ind_concat),noiseCorr_concat{post}(1,green_ind_concat),'MarkerEdgeColor','#47E5BC', 'LineWidth',1.25)
hold on
scatter(noiseCorr_concat{pre}(1,red_ind_concat),noiseCorr_concat{post}(1,red_ind_concat),'MarkerEdgeColor','#93748A', 'LineWidth',1.25)
xlim([-.5 1])
ylim([-.5 1])
refline(1)
xlabel('Noise corr pre-DART')
ylabel('Noise corr post-DART')
axis square
set(gca, 'TickDir', 'out')

print(fullfile(fnout,[ 'r_vs_LMI.pdf']),'-dpdf');



%% stat high and low R
hi_avrg_stat = cell(1,nd);
low_avrg_stat = cell(1,nd); 
hi_se_stat = cell(1,nd); 
low_se_stat = cell(1,nd);

%fidn cells with high correlation in the baseline day
highRInds = find(noiseCorr_concat{pre}(1,:)>0.5);
lowRInds = find(noiseCorr_concat{pre}(1,:)<=0.5);

redHigh=cell(nCon,1);
redLow=cell(nCon,1);

for id = 1:nd

   for iCon=1:nCon
       redHigh{iCon} =  intersect(haveRunning_red{iCon},highRInds); %find intersection of red cells with high R 
       % and red cells that I have running data for on both days, for this contrast
        redLow{iCon} =  intersect(haveRunning_red{iCon},lowRInds); %find intersection of red cells with low R 
       % and red cells that I have running data for on both days, for this contrast


        hi_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:, redHigh{iCon},iCon),2);
        high_std=nanstd(tc_trial_avrg_stat_concat{id}(:, redHigh{iCon},iCon),[],2);
        hi_se_stat{id}(:,iCon)=high_std/sqrt(length(redHigh{iCon}));
        
        low_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,redLow{iCon},iCon),2);
        low_std=nanstd(tc_trial_avrg_stat_concat{id}(:, redLow{iCon},iCon),[],2);
        low_se_stat{id}(:,iCon)=low_std/sqrt(length(redLow{iCon}));
        
        clear low_std high_std
    end
end
z=double(nOn)/double(frame_rate);

%creat a time axis in seconds
t=1:(size(hi_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure



subplot(1,2,1) 
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
ylabel('dF/F') 
xlabel('s') 
title(['Weakly correlated',' n = ', num2str(length(redLow{iCon}))])
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) 
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
title([' Strongly correlated',' n = ', num2str(length(redHigh{iCon}))])

xlabel('s') 
set(gca,'XColor', 'none','YColor','none')

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['stationary, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) 'stat_R_timecourses.pdf']),'-dpdf');
clear txt1 highRed lowRed
end 


%% scatterplot of max df/f for day 1 vs day 2 for high and low R red cells

for iCon = 1:nCon
   redHigh =  intersect(haveRunning_red{iCon},highRInds); %find intersection of red cells with high R 
   % and red cells that I have running data for on both days, for this contrast
    redLow =  intersect(haveRunning_red{iCon},lowRInds); %find intersection of red cells with low R 
   % and red cells that I have running data for on both days, for this contrast


figure; movegui('center') 
subplot(1,2,2)
scatter((pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)),(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
scatter((pref_responses_stat_concat{pre}(redHigh,iCon)),(pref_responses_stat_concat{post}(redHigh,iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on

mean_pre_stat = nanmean(pref_responses_stat_concat{pre}(redHigh,iCon));
mean_post_stat = nanmean(pref_responses_stat_concat{post}(redHigh,iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(redHigh,iCon)) / sqrt( length(redHigh));
stderror_post= std(pref_responses_stat_concat{post}(redHigh,iCon)) / sqrt( length(redHigh));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);


xlabel('pre-DART  dF/F')
ylim([-.2 0.2])
xlim([-.2 0.2])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('+HTP high r, stationary')
axis square
set(gca, 'TickDir', 'out')



subplot(1,2,1)
scatter((pref_responses_stat_concat{pre}(redLow,iCon)),(pref_responses_stat_concat{post}(redLow,iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
mean_pre_stat = nanmean(pref_responses_stat_concat{pre}(redLow,iCon));
mean_post_stat = nanmean(pref_responses_stat_concat{post}(redLow,iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(redLow,iCon)) / sqrt( length(redLow));
stderror_post= std(pref_responses_stat_concat{post}(redLow,iCon)) / sqrt( length(redLow));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylabel('post-DART dF/F')
ylim([-.2 0.2])
xlim([-.2 0.2])
set(gca, 'TickDir', 'out')
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('+HTP low r, stationary')
axis square


sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'stat_maxResp_R.pdf']),'-dpdf','-bestfit')
clear txt1 highRed lowRed
end

[h1, p1]= ttest(pref_responses_stat_concat{pre}(redLow,1),pref_responses_stat_concat{post}(redLow,1));
[h2,p2]= ttest(pref_responses_stat_concat{pre}(redHigh,1),pref_responses_stat_concat{post}(redHigh,1));

%correct for four tests
p1*2
p2*2
clear h1 p1 h2 p2

%% making tc plots for low and high R cells, running

hi_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
low_avrg_loc = cell(1,nd); %same for red
hi_se_loc = cell(1,nd); %this will be the se across all green cells
low_se_loc = cell(1,nd); %same for red

redHigh=cell(nCon);
redLow=cell(nCon);
for id = 1:nd
    for iCon=1:nCon
   redHigh{iCon} =  intersect(haveRunning_red{iCon},highRInds); %find intersection of red cells with high R 
   % and red cells that I have running data for on both days, for this contrast
    redLow{iCon} =  intersect(haveRunning_red{iCon},lowRInds); %find intersection of red cells with low R 
   % and red cells that I have running data for on both days, for this contrast


    hi_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:, redHigh{iCon},iCon),2);
    high_std=nanstd(tc_trial_avrg_loc_concat{id}(:, redHigh{iCon},iCon),[],2);
    hi_se_loc{id}(:,iCon)=high_std/sqrt(length(redHigh{iCon}));
    
    low_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,redLow{iCon},iCon),2);
    low_std=nanstd(tc_trial_avrg_loc_concat{id}(:, redLow{iCon},iCon),[],2);
    low_se_loc{id}(:,iCon)=low_std/sqrt(length(redLow{iCon}));
    
    clear low_std high_std
    end
end
z=double(nOn)/double(frame_rate);

%creat a time axis in seconds
t=1:(size(hi_avrg_loc{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure

subplot(1,2,1) %for +HTP
shadedErrorBar(t,low_avrg_loc{pre}(:,iCon),low_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,low_avrg_loc{post}(:,iCon),low_se_loc{post}(:,iCon),'b');
ylim([-.02 .5]);
hold on
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
%ylabel('dF/F') 
xlabel('s') 
title(['Weakly correlated',' n = ', num2str(length(redLow{iCon}))])
set(gca,'XColor', 'none','YColor','none')

subplot(1,2,2) 
ylim([-.02 .5]);
hold on
shadedErrorBar(t,hi_avrg_loc{pre}(:,iCon),hi_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,hi_avrg_loc{post}(:,iCon),hi_se_loc{post}(:,iCon),'b');
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['Strongly correlated',' n = ', num2str(length(redHigh{iCon}))])
ylabel('dF/F') 
xlabel('s') 




x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['running, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) 'loc_R_timecourses.pdf']),'-dpdf');
clear txt1 highRed lowRed
end 


%% scatterplot of max df/f for day 1 vs day 2 for high and low R red cells, running

for iCon = 1:nCon
   redHigh =  intersect(haveRunning_red{iCon},highRInds); %find intersection of red cells with high R 
   % and red cells that I have running data for on both days, for this contrast
    redLow =  intersect(haveRunning_red{iCon},lowRInds); %find intersection of red cells with low R 
   % and red cells that I have running data for on both days, for this contrast


figure; movegui('center') 
subplot(1,2,2)
scatter((pref_responses_loc_concat{pre}(haveRunning_red{iCon},iCon)),(pref_responses_loc_concat{post}(haveRunning_red{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
scatter((pref_responses_loc_concat{pre}(redHigh,iCon)),(pref_responses_loc_concat{post}(redHigh,iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on

mean_pre_stat = nanmean(pref_responses_loc_concat{pre}(redHigh,iCon));
mean_post_stat = nanmean(pref_responses_loc_concat{post}(redHigh,iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_loc_concat{pre}(redHigh,iCon)) / sqrt( length(redHigh));
stderror_post= std(pref_responses_loc_concat{post}(redHigh,iCon)) / sqrt( length(redHigh));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);


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
mean_pre_stat = nanmean(pref_responses_loc_concat{pre}(redLow,iCon));
mean_post_stat = nanmean(pref_responses_loc_concat{post}(redLow,iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_loc_concat{pre}(redLow,iCon)) / sqrt( length(redLow));
stderror_post= std(pref_responses_loc_concat{post}(redLow,iCon)) / sqrt( length(redLow));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);% ylabel('post-DART dF/F')
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
print(fullfile(fnout,[num2str(cons(iCon)) 'loc_maxResp_R.pdf']),'-dpdf','-bestfit')
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
% this is for all keep cells, not only those matched across running
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
            inds1 = intersect(green_ind_concat, mouseInds{iMouse});
            temp_means1(iMouse, iCon)=mean(pref_responses_stat_concat{id}(inds1,iCon),'omitnan');

            inds2 = intersect(red_ind_concat, mouseInds{iMouse});
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
z=12;
y=z-1;


for iCon = 1:nCon
figure; movegui('center') 
subplot(2,2,1)
% scatter((green_means_stat{pre}(:,iCon)),(green_means_stat{post}(:,iCon)),10,'filled')
% hold on
scatter((green_means_stat{pre}(z:nSess,iCon)),(green_means_stat{post}(z:nSess,iCon)),10,'filled')
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
% scatter(red_means_stat{pre}(:,iCon),red_means_stat{post}(:,iCon),10, 'filled')
% hold on
scatter(red_means_stat{pre}(z:nSess,iCon),red_means_stat{post}(z:nSess,iCon),10,'filled')
% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([0 .2])
xlim([0 .2])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
title('+HTP stationary')
axis square
hold off

subplot(2,2,3)
% scatter((green_means_loc{pre}(:,iCon)),(green_means_loc{post}(:,iCon)),10, 'filled')
% hold on
scatter((green_means_loc{pre}(z:nSess,iCon)),(green_means_loc{post}(z:nSess,iCon)),10,'filled')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([0 .4])
xlim([0 .4])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('-HTP running')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off

subplot(2,2,4)
% scatter((red_means_loc{pre}(:,iCon)),(red_means_loc{post}(:,iCon)),10, 'filled')
% hold on
scatter((red_means_loc{pre}(z:nSess,iCon)),(red_means_loc{post}(z:nSess,iCon)),10,'filled')
% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([0 .4])
xlim([0 .4])
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

%% compile capture values per imaging session
%
capture = getCaptureValues_annulus(mice);

%what do the distributions look like?
figure;
subplot(2,2,1)
histogram(capture(1,:),10)
title('Fluor in FOV')
xlim([0 200])
set(gca, 'TickDir', 'out')
box off
subplot(2,2,2)
histogram(capture(2,:),10)
xlim([0 200])
title('Fluor in CNTRL')
set(gca, 'TickDir', 'out')
box off
subplot(2,2,3)
histogram(capture(3,:),10)
title('FOV / CNTRL')
set(gca, 'TickDir', 'out')
box off


% identify the sessions with highest and lowest capture
mice{find(capture(3,:)==max(capture(3,:)))}
mice{find(capture(3,:)==min(capture(3,:)))}


% do high capture vlaues have more to do with FOV being bright or control
% being dim?
figure;
subplot(1,2,1)
scatter(capture(1,:),capture(3,:))
ylabel('Capture')
xlabel('FOV intensity')
set(gca, 'TickDir', 'out')
box off
axis square

subplot(1,2,2)
scatter(capture(2,:),capture(3,:))
ylabel('Capture')
xlabel('CNTRL intensity')
set(gca, 'TickDir', 'out')
box off
axis square

%% get difference measurements per mouse

for iCon = 1:nCon
    green_stat_mouse_diff(:,iCon) = (green_means_stat{post}(:,iCon)-green_means_stat{pre}(:,iCon));
    red_stat_mouse_diff(:,iCon) = (red_means_stat{post}(:,iCon)-red_means_stat{pre}(:,iCon));
    green_loc_mouse_diff(:,iCon) = (green_means_loc{post}(:,iCon)-green_means_loc{pre}(:,iCon));
    red_loc_mouse_diff(:,iCon) = (red_means_loc{post}(:,iCon)-red_means_loc{pre}(:,iCon));

% 
%     green_stat_mouse_diff(:,iCon) = (green_means_stat{pre}(:,iCon)-green_means_stat{post}(:,iCon))./green_means_stat{pre}(:,iCon);
%     red_stat_mouse_diff(:,iCon) = (red_means_stat{pre}(:,iCon)-red_means_stat{post}(:,iCon))./red_means_stat{pre}(:,iCon);
%     green_loc_mouse_diff(:,iCon) = (green_means_loc{pre}(:,iCon)-green_means_loc{post}(:,iCon))./green_means_loc{pre}(:,iCon);
%     red_loc_mouse_diff(:,iCon) = (red_means_loc{pre}(:,iCon)-red_means_loc{post}(:,iCon))./red_means_loc{pre}(:,iCon);
end
figure
subplot(2,2,1)
histogram(green_stat_mouse_diff(:,1))
subplot(2,2,2)
histogram(green_loc_mouse_diff(:,1))
subplot(2,2,3)
histogram(red_stat_mouse_diff(:,1))
subplot(2,2,4)
histogram(red_loc_mouse_diff(:,1))
%%


for iCon=1:nCon
figure
subplot (2,2,1)
scatter(capture(3,:),green_stat_mouse_diff(:,iCon),'filled')
hold on
%scatter(capture(3,z:nSess),green_stat_mouse_diff(z:nSess,iCon),'filled')
xlabel('Capture')
ylabel('Mean dFoF response (post-pre)')
title('Pyr cells stationary')
yline(0,':')
hold off
set(gca, 'TickDir', 'out')

subplot(2,2,2)
scatter(capture(3,:),red_stat_mouse_diff(:,iCon),'filled')
hold on
%scatter(capture(3,z:nSess),red_stat_mouse_diff(z:nSess,iCon),'filled')
xlabel('Capture')
ylabel('Mean dFoF response (post-pre)')
title('SOM cells stationary')
yline(0,':')
hold off
set(gca, 'TickDir', 'out')

subplot (2,2,3)
scatter(capture(3,:),green_loc_mouse_diff(:,iCon),'filled')
hold on
%scatter(capture(3,z:nSess),green_loc_mouse_diff(z:nSess,iCon),'filled')
xlabel('Capture')
ylabel('Mean dFoF response (post-pre)')
title('Pyr cells running')
yline(0,':')
hold off
set(gca, 'TickDir', 'out')

subplot(2,2,4)
scatter(capture(3,:),red_loc_mouse_diff(:,iCon),'filled')
hold on
%scatter(capture(3,z:nSess),red_loc_mouse_diff(z:nSess,iCon),'filled')
xlabel('Capture')
ylabel('Mean dFoF response (post-pre)')
title('SOM cells running')
yline(0,':')
hold off
set(gca, 'TickDir', 'out')


sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'diffVsFractCapt_byMouse.pdf']),'-dpdf','-bestfit')
end
%% 

captureByCell = [];
for iSess = 1:nSess
    temp = repmat(capture(3,iSess),[nKeep_concat(iSess),1]);
    captureByCell=[captureByCell;temp];
end

%% scatterplot of max df/f for by capture 

cutOff=1.3;
%cutOff=median(captureByCell);
%green is above capture threshold, purple is below capture threshold
for iCon = 1:nCon

figure; movegui('center') 
subplot(2,2,1)
scatter((pref_responses_stat_concat{pre}(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff),iCon)),(pref_responses_stat_concat{post}(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff),iCon)),10,'MarkerEdgeColor',[.75 .5 .75],'MarkerFaceColor',[.75 .5 .75],'jitter', 'on', 'jitterAmount',.01)
hold on
% scatter((pref_responses_stat_concat{pre}(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff),iCon)),(pref_responses_stat_concat{post}(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff),iCon)),10,'MarkerEdgeColor',[.5 .75 .5],'MarkerFaceColor',[.5 .75 .5],'jitter', 'on', 'jitterAmount',.01)

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
scatter((pref_responses_stat_concat{pre}(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff),iCon)),(pref_responses_stat_concat{post}(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff),iCon)),10,'MarkerEdgeColor',[.75 .5 .75],'MarkerFaceColor',[.75 .5 .75],'jitter', 'on', 'jitterAmount',.01)
hold on
% scatter((pref_responses_stat_concat{pre}(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff),iCon)),(pref_responses_stat_concat{post}(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff),iCon)),10,'MarkerEdgeColor',[.5 .75 .5],'MarkerFaceColor',[.5 .75 .5],'jitter', 'on', 'jitterAmount',.01)
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
scatter((pref_responses_loc_concat{pre}(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff),iCon)),(pref_responses_loc_concat{post}(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff),iCon)),10,'MarkerEdgeColor',[.75 .5 .75],'MarkerFaceColor',[.75 .5 .75],'jitter', 'on', 'jitterAmount',.01)
hold on
% scatter((pref_responses_loc_concat{pre}(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff),iCon)),(pref_responses_loc_concat{post}(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff),iCon)),10,'MarkerEdgeColor',[.5 .75 .5],'MarkerFaceColor',[.5 .75 .5],'jitter', 'on', 'jitterAmount',.01)

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
scatter((pref_responses_loc_concat{pre}(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff),iCon)),(pref_responses_loc_concat{post}(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff),iCon)),10,'MarkerEdgeColor',[.75 .5 .75],'MarkerFaceColor',[.75 .5 .75],'jitter', 'on', 'jitterAmount',.01)
hold on
% scatter((pref_responses_loc_concat{pre}(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff),iCon)),(pref_responses_loc_concat{post}(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff),iCon)),10,'MarkerEdgeColor',[.5 .75 .5],'MarkerFaceColor',[.5 .75 .5],'jitter', 'on', 'jitterAmount',.01)

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
print(fullfile(fnout,[num2str(cons(iCon)) 'maxResp_lowCapture.pdf']),'-dpdf','-bestfit')
end
%% SOM/Pyr relationship session-by-session, with colorscale for capture

lowCapture=find(capture(3,:)<cutOff);
highCapture = find(capture(3,:)>=cutOff);

for iCon=1:nCon
figure
subplot (1,2,1)
scatter(green_stat_mouse_diff(lowCapture,iCon),red_stat_mouse_diff(lowCapture,iCon),'MarkerEdgeColor',[.75 .5 .75],'MarkerFaceColor',[.75 .5 .75])%for the session with low capture
hold on
scatter(green_stat_mouse_diff(highCapture,iCon),red_stat_mouse_diff(highCapture,iCon),'MarkerEdgeColor',[.5 .75 .5],'MarkerFaceColor',[.5 .75 .5])%for the session with high capture
xlabel('Pyr')
ylabel('SOM')
title('Stationary')
hold off
set(gca, 'TickDir', 'out')
axis square

subplot(1,2,2)
scatter(green_loc_mouse_diff(lowCapture,iCon),red_loc_mouse_diff(lowCapture,iCon),'MarkerEdgeColor',[.75 .5 .75],'MarkerFaceColor',[.75 .5 .75])%for the session with low capture
hold on
scatter(green_stat_mouse_diff(highCapture,iCon),red_stat_mouse_diff(highCapture,iCon),'MarkerEdgeColor',[.5 .75 .5],'MarkerFaceColor',[.5 .75 .5])%for the session with high capture
xlabel('Pyr')
ylabel('SOM')
title('Running')
hold off
set(gca, 'TickDir', 'out')
axis square




sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'SOM_Pyr_Capt_byMouse.pdf']),'-dpdf','-bestfit')
end

%% timecourses seperated by capture


%make figure with se shaded, one figure per contrast - stationary

tc_green_avrg_stat = cell(1,nd); 
tc_red_avrg_stat = cell(1,nd); 
tc_green_se_stat = cell(1,nd); 
tc_red_se_stat = cell(1,nd); 

tc_green_avrg_stat_b = cell(1,nd); 
tc_red_avrg_stat_b = cell(1,nd); 
tc_green_se_stat_b = cell(1,nd); 
tc_red_se_stat_b = cell(1,nd);


for id = 1:nd
    for iCon=1:nCon
        
    tc_green_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff),iCon),2);
    green_std=nanstd(tc_trial_avrg_stat_concat{id}(:,haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff),iCon),[],2);
    tc_green_se_stat{id}(:,iCon)=green_std/sqrt(length(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff)));
    
    tc_red_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff),iCon),2);
    red_std=nanstd(tc_trial_avrg_stat_concat{id}(:,haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff),iCon),[],2);
    tc_red_se_stat{id}(:,iCon)=red_std/sqrt(length(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff)));


    tc_green_avrg_stat_b{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff),iCon),2);
    green_std=nanstd(tc_trial_avrg_stat_concat{id}(:,haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff),iCon),[],2);
    tc_green_se_stat_b{id}(:,iCon)=green_std/sqrt(length(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff)));
    
    tc_red_avrg_stat_b{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff),iCon),2);
    red_std=nanstd(tc_trial_avrg_stat_concat{id}(:,haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff),iCon),[],2);
    tc_red_se_stat_b{id}(:,iCon)=red_std/sqrt(length(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff)));

    clear green_std red_std
    end
end
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(2,2,1) %for the first day

ylim([-.05 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['low capture -HTP',' n = ', num2str(length(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff)))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(2,2,2) %for the second day
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
title(['low capture +HTP',' n = ', num2str(length(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff)))])
set(gca,'XColor', 'none','YColor','none')

subplot(2,2,3) %for the first day

ylim([-.05 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_stat_b{pre}(:,iCon),tc_green_se_stat_b{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_stat_b{post}(:,iCon),tc_green_se_stat_b{post}(:,iCon),'b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['high capture -HTP',' n = ', num2str(length(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff)))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(2,2,4) %for the second day
shadedErrorBar(t,tc_red_avrg_stat_b{pre}(:,iCon),tc_red_se_stat_b{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_stat_b{post}(:,iCon),tc_red_se_stat_b{post}(:,iCon),'b');
ylim([-.05 .3]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['high capture +HTP',' n = ', num2str(length(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff)))])

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['stationary, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_stat_capt_timecourses.pdf']),'-dpdf');
end 

clear tc_green_avrg_stat tc_red_avrg_stat tc_green_se_stat tc_red_se_stat tc_green_avrg_stat_b tc_red_avrg_stat_b tc_green_se_stat_b tc_red_se_stat_b 

% make figure with se shaded, one figure per contrast - loc

tc_green_avrg_loc = cell(1,nd); 
tc_red_avrg_loc = cell(1,nd); 
tc_green_se_loc = cell(1,nd); 
tc_red_se_loc = cell(1,nd); 

tc_green_avrg_loc_b = cell(1,nd); 
tc_red_avrg_loc_b = cell(1,nd); 
tc_green_se_loc_b = cell(1,nd); 
tc_red_se_loc_b = cell(1,nd); 


for id = 1:nd
    for iCon=1:nCon
        
    tc_green_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff),iCon),2);
    green_std=nanstd(tc_trial_avrg_loc_concat{id}(:,haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff),iCon),[],2);
    tc_green_se_loc{id}(:,iCon)=green_std/sqrt(length(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff)));
    
    tc_red_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff),iCon),2);
    red_std=nanstd(tc_trial_avrg_loc_concat{id}(:,haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff),iCon),[],2);
    tc_red_se_loc{id}(:,iCon)=red_std/sqrt(length(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff)));


    tc_green_avrg_loc_b{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff),iCon),2);
    green_std=nanstd(tc_trial_avrg_loc_concat{id}(:,haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff),iCon),[],2);
    tc_green_se_loc_b{id}(:,iCon)=green_std/sqrt(length(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff)));
    
    tc_red_avrg_loc_b{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff),iCon),2);
    red_std=nanstd(tc_trial_avrg_loc_concat{id}(:,haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff),iCon),[],2);
    tc_red_se_loc_b{id}(:,iCon)=red_std/sqrt(length(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff)));
    



    clear green_std red_std
    end
end
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_loc{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);



for iCon = 1:nCon
figure
subplot(2,2,1) %for the first day

ylim([-.05 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_loc{pre}(:,iCon),tc_green_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_loc{post}(:,iCon),tc_green_se_loc{post}(:,iCon),'b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['low capture -HTP',' n = ', num2str(length(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff)))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(2,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_loc{pre}(:,iCon),tc_red_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_loc{post}(:,iCon),tc_red_se_loc{post}(:,iCon),'b');
ylim([-.05 .3]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['low capture +HTP',' n = ', num2str(length(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff)))])
set(gca,'XColor', 'none','YColor','none')

subplot(2,2,3) %for the first day

ylim([-.05 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_loc_b{pre}(:,iCon),tc_green_se_loc_b{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_loc_b{post}(:,iCon),tc_green_se_loc_b{post}(:,iCon),'b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['high capture -HTP',' n = ', num2str(length(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff)))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(2,2,4) %for the second day
shadedErrorBar(t,tc_red_avrg_loc_b{pre}(:,iCon),tc_red_se_loc_b{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_loc_b{post}(:,iCon),tc_red_se_loc_b{post}(:,iCon),'b');
ylim([-.05 .3]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['high capture +HTP',' n = ', num2str(length(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff)))])

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['running, contrast = ' num2str(cons(iCon))])


print(fullfile(fnout,[num2str(cons(iCon)) '_loc_capt_timecourses.pdf']),'-dpdf');
end 


clear tc_green_avrg_loc tc_red_avrg_loc tc_green_se_loc tc_red_se_loc tc_green_avrg_loc_b tc_red_avrg_loc_b tc_green_se_loc_b tc_red_se_loc_b 

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
% subplot(1,2,1);
% scatter(meanF_concat{pre}(green_ind_concat),meanF_concat{post}(green_ind_concat),10,'MarkerEdgeColor',[.5 .5 .5])
% title('-HTP')
% xlim([min_lim max_lim_green]);
% ylim([min_lim max_lim_green]);
% hline=refline(1);
% hline.Color = 'k';
% hline.LineStyle = ':';
% set(gca, 'TickDir', 'out')
% uistack(hline,'bottom');
% xlabel('pre-DART')
% ylabel('post-DART')
% axis square

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

subplot(1,2,1) %for the first day
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_loc{pre}(:,iCon),tc_green_se_loc{pre}(:,iCon),'m');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['-HTP pre-DART',' n = ', num2str(length(haveRunning_green{iCon}))])
ylim([-.05 .4])
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
ylim([-.05 .4])
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

subplot(1,2,1) %for the first day
shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon),tc_red_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_loc{pre}(:,iCon),tc_red_se_loc{pre}(:,iCon),'m');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['+HTP pre-DART',' n = ', num2str(length(haveRunning_red{iCon}))])
ylim([-.05 .4])
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
ylim([-.05 .4])
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

%
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
%% comparing R to full tc corre

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
mean_pre_stat = nanmean(pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon));
mean_post_stat = nanmean(pref_resp_subtracted{post}(haveRunning_green{iCon},iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
stderror_post= std(pref_resp_subtracted{post}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);
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
mean_pre_stat = nanmean(pref_resp_subtracted{pre}(haveRunning_red{iCon},iCon));
mean_post_stat = nanmean(pref_resp_subtracted{post}(haveRunning_red{iCon},iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_resp_subtracted{pre}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
stderror_post= std(pref_resp_subtracted{post}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);% ylabel('post-DART dF/F')
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


clear mean_pre_stat mean_post_stat stderror_post stderror_pre
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
std_pre_stat = std(pref_resp_subtracted{pre}(haveRunning_green{1},1),[],1);
std_post = std(pref_resp_subtracted{post}(haveRunning_green{1},1),[],1);
se_pre = std_pre_stat/sqrt(length(haveRunning_green{1}));
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
ylim([0 .4])
xticks([.25 .5 1])
set(gca, 'TickDir', 'out')
box off

subplot(1,2,2) 
errorbar(cons,conResp_red_avrg_loc{pre},conResp_red_se_loc{pre},'k');
hold on
errorbar(cons,conResp_red_avrg_loc{post},conResp_red_se_loc{post},'b');
title(['+HTP',' n = ', num2str(length(red_all))])
 ylim([0 .4])
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

%% anova on contrast response
DARTeffect = pref_responses_stat_concat{pre}(red_ind_concat,:)-pref_responses_stat_concat{post}(red_ind_concat,:);


%% time courses averaging stationary and running

tc_green_avrg_allCondition = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_allCondition = cell(1,nd); %same for red
tc_green_se_allCondition = cell(1,nd); %this will be the se across all green cells
tc_red_se_allCondition = cell(1,nd); %same for red



for id = 1:nd
    for iCon=1:nCon
        
    tc_green_avrg_allCondition{id}(:,iCon)=nanmean( tc_trial_avrg_keep_allCond_concat{id}(:,green_ind_concat,iCon),2);
    green_std=nanstd( tc_trial_avrg_keep_allCond_concat{id}(:,green_ind_concat,iCon),[],2);
    tc_green_se_allCondition{id}(:,iCon)=green_std/sqrt(length(green_ind_concat));
    
    tc_red_avrg_allCondition{id}(:,iCon)=nanmean( tc_trial_avrg_keep_allCond_concat{id}(:,red_ind_concat,iCon),2);
    red_std=nanstd( tc_trial_avrg_keep_allCond_concat{id}(:,red_ind_concat,iCon),[],2);
    tc_red_se_allCondition{id}(:,iCon)=red_std/sqrt(length(red_ind_concat));
    
    clear green_std red_std
    end
end
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_allCondition{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day



ylim([-.02 .25]);;
hold on
shadedErrorBar(t,tc_green_avrg_allCondition{pre}(:,iCon),tc_green_se_allCondition{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_allCondition{post}(:,iCon),tc_green_se_allCondition{post}(:,iCon),'b','transparent');
hold on
%line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
fill([.0 0 z z],[-.02 .3 .3 -.02],'k',FaceAlpha = 0.25,LineStyle='none')
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['-HTP',' n = ', num2str(length(green_ind_concat))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_allCondition{pre}(:,iCon),tc_red_se_allCondition{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_allCondition{post}(:,iCon),tc_red_se_allCondition{post}(:,iCon),'b');
ylim([-.02 .25]);;
hold on
%line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
fill([.0 0 z z],[-.02 .3 .3 -.02],'k',FaceAlpha = 0.25,LineStyle='none')
title(['+HTP',' n = ', num2str(length(red_ind_concat))])

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['allConditionionary, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_allCondition_allKeep_timecourses.pdf']),'-dpdf');
end ; 

%% finding how the high and low R cells are distributed over epxeriments

RbyExp = zeros(2,nSess);
for iSess = 1:nSess
    mouseIndsTemp = mouseInds{iSess};
    RbyExp(1,iSess) = length(intersect(mouseIndsTemp,highRInds));
    RbyExp(2,iSess) = length(intersect(mouseIndsTemp,lowRInds));
end

