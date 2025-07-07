
clear all; clear global; close all
clc
ds = 'DART_V1_YM90K_Celine'; %dataset info
% ds = 'DART_expt_info'
dataStructLabels = {'contrastxori'};
experimentFolder = 'VIP_YM90K';
% experimentFolder = 'PV_atropine';

rc =  behavConstsDART; %directories
eval(ds);
%285 295 300 308 324 334 DART YM90K 
% 299 289 304 312 320 330
sess_list = [2];%enter all the sessions you want to concatenate4
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

targetCon = [.125 .25 .5 1]%what contrast to extract for all data - must be one that all datasets had

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
        
        fnout = fullfile(rc.achAnalysis,experimentFolder,expt(sess_list(1)).mouse,['multiday_' dart_str],d);
else
    fnout= fullfile(rc.achAnalysis,experimentFolder,strcat('concat', sess_title),d);
end
mkdir(fnout);
cd(fnout)
clear d sess_title

nCon = length(targetCon)
nSize =5;

mice=[];
red_concat=[];
green_concat=[];
nKeep_concat=[];

tc_trial_avrg_stat_concat=cell(1,nd);
tc_trial_avrg_loc_concat=cell(1,nd);
conBySize_resp_stat_concat=cell(1,nd);
conBySize_resp_loc_concat=cell(1,nd);
h_concat=cell(1,nd);
data_resp_concat=cell(1,nd);

tc_trial_avrg_stat_largePupil_concat=cell(1,nd);
tc_trial_avrg_stat_smallPupil_concat=cell(1,nd);
resp_keep_concat=cell(1,nd);
resp_max_keep_concat=cell(1,nd);
pref_responses_loc_concat=cell(1,nd);
pref_responses_stat_concat=cell(1,nd);
pref_peak_stat_concat=cell(1,nd);
pref_peak_loc_concat=cell(1,nd);
pref_responses_stat_largePupil_concat=cell(1,nd);
pref_responses_stat_smallPupil_concat=cell(1,nd);
RIx_concat=cell(1,nd);
dirs_concat=[];
cons_concat=[];
dfof_max_diff_concat=[];
red_fluor_concat=[];
green_fluor_concat=[];
wheel_corr_concat=cell(1,nd);
meanF_concat=cell(1,nd);
norm_dir_resp_stat_concat = cell(1,nd);
norm_dir_resp_loc_concat = cell(1,nd);
pref_nonPref_stat_concat=cell(1,nd);
pref_nonPref_loc_concat=cell(1,nd);
pref_dir_concat=cell(1,nd);
noiseCorr_concat = cell(1,nd);
sigCorr_concat = cell(1,nd);
pref_allTrials_stat_concat =cell(nCon,nSize,nd);
pref_allTrials_loc_concat =cell(nCon,nSize,nd);
pref_allTrials_largePupil_concat =cell(nCon,nSize,nd);
pref_allTrials_smallPupil_concat =cell(nCon,nSize,nd);
dataTableConat=[];
drug=cell(1,nSess);
pupilMeans_concat=nan(nd,3,nSess);
motorByPupil_concat=nan(nd,2,nSess);
pupilCounts_concat=nan(nd,2,nSess);
nonPref_trial_avrg_stat_concat=cell(1,nd);
nonPref_trial_avrg_loc_concat=cell(1,nd);
data_dfof_runOnset_concat=cell(1,nd);

responseByCondProps_concat=nan(6,2,nSess);

for iSess = 1:nSess
    thisSess = sess_list(iSess);
    mouse = expt(thisSess).mouse;
    mice=[mice;mouse];
    thisDrug = expt(thisSess).drug;
    drug{iSess}=thisDrug;
    
    if expt(thisSess).multiday_timesincedrug_hours>0
        dart_str = [expt(thisSess).drug '_' num2str(expt(thisSess).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fn_multi = fullfile(rc.achAnalysis,experimentFolder,mouse,['multiday_' dart_str]);

    load(fullfile(fn_multi,'tc_keep.mat'));
    load(fullfile(fn_multi,'resp_keep.mat'));
    load(fullfile(fn_multi,'input.mat'));
    load(fullfile(fn_multi,'locomotion.mat'));
%    load(fullfile(fn_multi,'fluor_intensity.mat'));
    load(fullfile(fn_multi,'HT_pyr_relationship.mat'));
    load(fullfile(fn_multi,'pupilMeans.mat'));

    nKeep = size(tc_trial_avrg_stat{post},2);

   pupilMeans_concat(:,:,iSess)=pupilMeans;
   motorByPupil_concat(:,:,iSess)=motorByPupil;
%   pupilCounts_concat(:,:,iSess)=pupilCounts;


    %tells the contrast, direction and orientation for each trial each day
    tCon_match = cell(1,nd);
    tDir_match = cell(1,nd);
    tOri_match = cell(1,nd);
    tSize_match = cell(1,nd);

    %find the contrasts, directions and orientations for each day
    for id = 1:nd
        tCon_match{id} = celleqel2mat_padded(input(id).tGratingContrast(1:nTrials(id)));
        tDir_match{id} = celleqel2mat_padded(input(id).tGratingDirectionDeg(1:nTrials(id)));
        tOri_match{id} = tDir_match{id};
        tOri_match{id}(find(tDir_match{id}>=180)) = tDir_match{id}(find(tDir_match{id}>=180))-180;
        tSize_match{id} = celleqel2mat_padded(input(id).tGratingDiameterDeg(1:nTrials(id)));
    end
    dirs = unique(tDir_match{post});
    cons = unique(tCon_match{post});
    sharedCon=find(ismember(cons, targetCon));
    sizes = unique(tSize_match{1});

    nOn = input(1).nScansOn;
    nOff = input(1).nScansOff;
    %start conatenating
    dirs_concat = [dirs_concat,dirs]; 
    cons_concat = [cons_concat,cons(sharedCon)];
    red_concat = [red_concat, red_keep_logical];
    green_concat = [green_concat, green_keep_logical];
    nKeep_concat = [nKeep_concat,nKeep];
    % responseByCondProps_concat(:,:,iSess)=responseByCondProps;

    clear cons
    
    
    for id = 1:nd
        
        tc_trial_avrg_stat_concat{id} =cat(2,tc_trial_avrg_stat_concat{id},tc_trial_avrg_stat{id}(:,:,:,:));
        tc_trial_avrg_stat_largePupil_concat{id} = cat(2,tc_trial_avrg_stat_largePupil_concat{id},tc_trial_avrg_stat_largePupil{id}(:,:,sharedCon,:));
        tc_trial_avrg_stat_smallPupil_concat{id} = cat(2,tc_trial_avrg_stat_smallPupil_concat{id},tc_trial_avrg_stat_smallPupil{id}(:,:,sharedCon,:));
        tc_trial_avrg_loc_concat{id} =cat(2,tc_trial_avrg_loc_concat{id},tc_trial_avrg_loc{id}(:,:,sharedCon,:));
        nonPref_trial_avrg_stat_concat{id} =cat(2,nonPref_trial_avrg_stat_concat{id},nonPref_trial_avrg_stat{id}(:,:,sharedCon,:));
        nonPref_trial_avrg_loc_concat{id} =cat(2,nonPref_trial_avrg_loc_concat{id},nonPref_trial_avrg_loc{id}(:,:,sharedCon,:));
        resp_keep_concat{id}=cat(1,resp_keep_concat{id},resp_keep{id});
        resp_max_keep_concat{id}=cat(1,resp_max_keep_concat{id},resp_max_keep{id}(:,sharedCon,:));
        pref_responses_loc_concat{id}=cat(1,pref_responses_loc_concat{id},pref_responses_loc{id}(:,sharedCon,:));
        pref_responses_stat_concat{id}=cat(1,pref_responses_stat_concat{id},pref_responses_stat{id}(:,sharedCon,:));
        pref_peak_stat_concat{id}=cat(1,pref_peak_stat_concat{id},pref_peak_stat{id}(:,sharedCon,:));
        pref_peak_loc_concat{id}=cat(1,pref_peak_loc_concat{id},pref_peak_loc{id}(:,sharedCon,:));
        pref_responses_stat_largePupil_concat{id}=cat(1,pref_responses_stat_largePupil_concat{id},pref_responses_stat_largePupil{id}(:,sharedCon,:));
        pref_responses_stat_smallPupil_concat{id}=cat(1,pref_responses_stat_smallPupil_concat{id},pref_responses_stat_smallPupil{id}(:,sharedCon,:));
        RIx_concat{id}=cat(1,RIx_concat{id},sum(RIx{id}));
        wheel_corr_concat{id}=cat(2,wheel_corr_concat{id},wheel_corr{id});
        meanF=mean(fullTC_keep{id},1);
        meanF_concat{id}=cat(2,meanF_concat{id}, meanF);
        norm_dir_resp_stat_concat{id}=cat(1,norm_dir_resp_stat_concat{id},norm_dir_resp_stat{id});
        norm_dir_resp_loc_concat{id}=cat(1,norm_dir_resp_loc_concat{id},norm_dir_resp_loc{id});
        pref_dir_concat{id}=cat(2,pref_dir_concat{id},pref_dir_keep{id});
        noiseCorr_concat{id}=cat(2,noiseCorr_concat{id},noiseCorr{id});
        sigCorr_concat{id}=cat(2,sigCorr_concat{id},sigCorr{id});
        h_concat{id}=cat(1,h_concat{id},h_keep{id});
        conBySize_resp_stat_concat{id}=cat(1,conBySize_resp_stat_concat{id},conBySize_resp_stat_keep{id});
        conBySize_resp_loc_concat{id}=cat(1,conBySize_resp_loc_concat{id},conBySize_resp_loc_keep{id});
        data_resp_concat{id} = cat(1,data_resp_concat{id},data_resp_keep{id});
        data_dfof_runOnset_concat{id}=cat(2,data_dfof_runOnset_concat{id},data_dfof_runOnset_keep{id});

        for i = 1:length(sharedCon)
            iCon=sharedCon(i);
            for iSize = 1:length(sizes)
                pref_allTrials_stat_concat{i,iSize,id}=[pref_allTrials_stat_concat{i,iSize,id},pref_allTrials_stat{iCon,iSize,id}];
                pref_allTrials_loc_concat{i,iSize,id}=[pref_allTrials_loc_concat{i,iSize,id},pref_allTrials_loc{iCon,iSize,id}];
                pref_allTrials_largePupil_concat{i,iSize,id}=[pref_allTrials_largePupil_concat{i,iSize,id},pref_allTrials_largePupil{iCon,iSize,id}];
                pref_allTrials_smallPupil_concat{i,iSize,id}=[pref_allTrials_smallPupil_concat{i,iSize,id},pref_allTrials_smallPupil{iCon,iSize,id}];

            end
        end
        clear meanF i
    end
    dfof_max_diff_concat=cat(1,dfof_max_diff_concat,dfof_max_diff(:,sharedCon,:));
   % green_fluor_concat=cat(2,green_fluor_concat,green_fluor_keep);
   % red_fluor_concat=cat(2,red_fluor_concat,red_fluor_keep);
    
iSess
end


red_ind_concat = find(red_concat);
green_ind_concat = find(green_concat);
cons = targetCon;
nSize = length(sizes)
nKeep_total = sum(nKeep_concat);

clear data_resp_keep data_trial_keep green_ind_keep red_ind_keep conBySize_resp_loc_keep 
clear conBySize_resp_stat_keep tc_trial_avrg_stat tc_trial_avrg_loc conBySize_resp_loc_match 
clear conBySize_resp_stat_match red_keep_logical green_keep_logical pref_responses_stat 
clear pref_responses_loc resp_keep mouse pref_con_keep pref_dir_keep pref_size_keep 
clear dfof_max_diff dfof_max_diff_raw green_fluor_match green_fluor_keep motorByPupil nKeep noiseCorr
clear sigCorr nonPref_trial_avrg_loc nonPref_trial_avrg_stat norm_dir_resp_loc norm_dir_resp_stat 
clear pref_allTrials_largePupil pref_allTrials_smallPupil pref_allTrials_loc pref_allTrials_stat
clear pref_peak_loc pref_peak_stat pref_responses_stat_smallPupil pref_allTrials_largePupil
clear pupilMeans red_fluor_match red_fluor_keep RIx tc_trial_avrg_keep_allCond tc_trial_avrg_stat_largePupil
clear tc_trial_avrg_stat_smallPupil wheel_corr wheel_tc dart_Str expt fullTC_keep wheel_speed

% cell selection
% find cells that I have running and stationary data for on both days

haveRunning_pre=(sum(squeeze((sum(~isnan(conBySize_resp_loc_concat{pre}),2))),2)==nSize*nCon); %find cells that have a no NAN values for any size (the total of non-NAN should == the number of size)
haveRunning_post=(sum(squeeze((sum(~isnan(conBySize_resp_loc_concat{post}),2))),2)==nSize*nCon);
haveRunning_both= find(haveRunning_pre.* haveRunning_post); %find cells that meet this criteria for both days - now in indices, not logical

haveStat_pre=(sum(squeeze((sum(~isnan(conBySize_resp_stat_concat{pre}),2))),2)==nSize*nCon); %find cells that have a no NAN values for any size (the total of non-NAN should == the number of size)
haveStat_post=(sum(squeeze((sum(~isnan(conBySize_resp_stat_concat{post}),2))),2)==nSize*nCon); 
haveStat_both= find(haveStat_pre.* haveStat_post); %find cells that meet this criteria for both days - now in indices, not logical

runningCells = intersect(haveStat_both, haveRunning_both);


% %have running for large only
% 
% haveRunning_pre=squeeze((sum(~isnan(conBySize_resp_loc_concat{pre}(:,:,2)),2)))==nCon; %find cells that have a no NAN values for any size (the total of non-NAN should == the number of size)
% haveRunning_post=squeeze((sum(~isnan(conBySize_resp_loc_concat{post}(:,:,2)),2)))==nCon;
% haveRunning_both= find(haveRunning_pre.* haveRunning_post); %find cells that meet this criteria for both days - now in indices, not logical
% 
% haveStat_pre=squeeze((sum(~isnan(conBySize_resp_stat_concat{pre}(:,:,2)),2)))==nCon; %find cells that have a no NAN values for any size (the total of non-NAN should == the number of size)
% haveStat_post=squeeze((sum(~isnan(conBySize_resp_stat_concat{post}(:,:,2)),2)))==nCon;
% haveStat_both= find(haveStat_pre.* haveStat_post); %find cells that meet this criteria for both days - now in indices, not logical
% 
% runningCells = intersect(haveStat_both, haveRunning_both);


%find cells responsive to a given size
mySize = 1; %1 is the smallest size
respToSizeBothDays = cell(1,nd);
for id = 1:nd
    respToSizeBothDays{id}=sum(squeeze(sum(h_concat{id}(:,:,:,mySize),2)),2); %finding cells that responded this size, at any contrast or direction
end
respToSmall = logical(respToSizeBothDays{pre}+respToSizeBothDays{post}); %to find cells that were responsive to this size on either day

mySize = nSize; %nSize is the largest size
respToSizeBothDays = cell(1,nd);
for id = 1:nd
    respToSizeBothDays{id}=sum(squeeze(sum(h_concat{id}(:,:,:,mySize),2)),2); %finding cells that responded this size, at any contrast or direction
end
respToLarge = logical(respToSizeBothDays{pre}+respToSizeBothDays{post}); %to find cells that were responsive to this size on either day


responCriteria = cell(1,nd); %cell array that will have indices of cells that meet our response criteria on each day
for id = 1:nd
    responseCheck =sum(squeeze(sum(h_concat{id}(:,:,2:4,nSize),2)),2); %finding cells that responded large size size, within the top three contrasts
    responCriteria{id}=find(logical(responseCheck));
end

includeCells = intersect(responCriteria{pre},find(haveRunning_pre));

clear haveRunning_pre haveRunning_post haveRunning_both haveStat_both haveStat_pre haveStat_post


% to find the OSI of each cell

OSI_baseline=nan(nKeep_total,1);

dirMean = mean(squeeze(mean(data_resp_concat{pre}(:,:,:,:,1),3)),3);
dirMean(find(dirMean<0))=0;
%left with average response at each dir, averaged over size and contrast,
%for each cell, for the baseline day only
nDir = length(dirs_concat);

% if nDir > 1
% 
%     orthOrder=[3     4     1     2];
%     for iCell = 1:nKeep_total
%         [prefResp, prefDir] = max(dirMean(iCell,:));
%         orthInd = orthOrder(prefDir);
%         orthResp = dirMean(iCell,orthInd);
%         OSI_baseline(iCell)=(prefResp-orthResp)/(prefResp+orthResp);
%     end
% end

% %histogram(OSI_baseline);
% green_ind_concat = intersect(green_ind_concat, find(respToLarge));
% red_ind_concat = intersect(red_ind_concat, find(respToLarge));

runningGreen = intersect(runningCells, green_ind_concat);
runningRed= intersect(runningCells, red_ind_concat);


statGreen = green_ind_concat;
statRed = red_ind_concat;

%get an array with indices for each mouse
mouseInds=cell(1,nSess);
start=1;
for iMouse = 1:nSess
    mouseInds{iMouse}=start:(start-1)+nKeep_concat(iMouse);
    start = start+nKeep_concat(iMouse)
end
clear start iMouse


%find how many haveRunning red cells I got for each mouse
cellCounts = nan(nSess,2);
mouseNames=[];
for iMouse = 1:nSess
   cellCounts(iMouse,1)=length(intersect(runningRed',(mouseInds{iMouse})));
   cellCounts(iMouse,2)=length(intersect(statRed',(mouseInds{iMouse})));
   mouseNames=[mouseNames, string(mice(iMouse,:))];
end


cellCountsGreen = nan(nSess,2);
mouseNames=[];
for iMouse = 1:nSess
   cellCountsGreen(iMouse,1)=length(intersect(runningGreen',(mouseInds{iMouse})));
   cellCountsGreen(iMouse,2)=length(intersect(statGreen',(mouseInds{iMouse})));

   mouseNames=[mouseNames, string(mice(iMouse,:))];
end


cellCountTable = table(cellCounts, RowNames=mouseNames)
cellCountTableGreen = table(cellCountsGreen, RowNames=mouseNames)
writetable(cellCountTable,fullfile(fnout,['cellCounts.csv']),'WriteRowNames',true)
clear cellCounts cellCountsGreen

%to find the cells that have 
runningByCondition = nan(nKeep_total,nCon,nSize);

for iCon = 1:nCon
    for iSize = 1:nSize
        runningPre = ~isnan(conBySize_resp_loc_concat{pre}(:,iCon, iSize));
        runningPost = ~isnan(conBySize_resp_loc_concat{post}(:,iCon, iSize));
        runningByCondition(:,iCon,iSize)=runningPre.*runningPost;
    end
end

%%  plot shaded error bar for running onsets

tc_green_avrg_runOnset = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_runOnset = cell(1,nd); %same for red
tc_green_se_runOnset = cell(1,nd); %this will be the se across all green cells
tc_red_se_runOnset = cell(1,nd); %same for red


for id = 1:nd

        tc_green_avrg_runOnset{id}=nanmean(data_dfof_runOnset_concat{id}(:,green_ind_concat),2);
        green_std=nanstd(data_dfof_runOnset_concat{id}(:,green_ind_concat),[],2);
        tc_green_se_runOnset{id}=green_std/sqrt(length(green_ind_concat));
        
        tc_red_avrg_runOnset{id}=nanmean(data_dfof_runOnset_concat{id}(:,red_ind_concat),2);
        red_std=nanstd(data_dfof_runOnset_concat{id}(:,red_ind_concat),[],2);
        tc_red_se_runOnset{id}=red_std/sqrt(length(red_ind_concat));
        
        
        clear green_std red_st
end
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_red_avrg_runOnset{id},1));
t=(t-(double(frame_rate)-1))/double(frame_rate);

figure
    subplot(1,2,1) %for the first day

    ylim([-.01 .05]);
    hold on
    shadedErrorBar(t,tc_green_avrg_runOnset{pre},tc_green_se_runOnset{pre},'--k');
    hold on
    shadedErrorBar(t,tc_green_avrg_runOnset{post},tc_green_se_runOnset{post},'--b','transparent');
    %line([-.8,-.8],[0,.01],'Color','black','LineWidth',2);
    title(['-HTP',' n = ', num2str(length(green_ind_concat))])
    ylabel('dF/F') 
    xlabel('s') 
    box off
    
    
    subplot(1,2,2) %+HTP
    shadedErrorBar(t,tc_red_avrg_runOnset{pre},tc_red_se_runOnset{pre},'k');
    hold on
    shadedErrorBar(t,tc_red_avrg_runOnset{post},tc_red_se_runOnset{post},'b');
    ylim([-.01 .05]);
    %line([-.8,-.8],[0,.01],'Color','black','LineWidth',2);
    ylabel('dF/F') 
    xlabel('s') 
    title(['+HTP',' n = ', num2str(length(red_ind_concat))])
    
    x0=5;
    y0=5;
    width=4;
    height=3;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    box off
    sgtitle('Running onset')
print(fullfile(fnout,'runOnset_timecourse.pdf'),'-dpdf');

