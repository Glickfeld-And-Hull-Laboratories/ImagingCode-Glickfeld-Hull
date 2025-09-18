
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
sess_list = [2 4 6 16];%enter all the sessions you want to concatenate
doEye = 1; %analyze pupil info?
nSess=length(sess_list);
targetCon = [.125 .25 .5 1]%what contrast to extract for all data - must be one that all datasets had
nCon = length(targetCon)
nSize =2;
frame_rate = 15;

onlySS = 0;
onlyNOTSS = 0;
onlySmall = 0;
onlyLarge = 0;

if onlySS * onlyNOTSS == 1 | onlySmall * onlyLarge == 1
    error('ERROR: Indexing by mutually exclusive conditions.');
end


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
% fnout = "G:\home\ACh\Analysis\2p_analysis\VIP_YM90K\concat2_4_6_16\SSandRespSmallOnly";
mkdir(fnout);
cd(fnout)
clear d sess_title

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
trial_dfof_concat=cell(nSess,nd);
tSize_concat = cell(nSess,nd);
tCon_concat = cell(nSess,nd);
tRIx_concat = cell(nSess,nd);

if doEye == 1
    tc_trial_avrg_stat_largePupil_concat=cell(1,nd);
    tc_trial_avrg_stat_smallPupil_concat=cell(1,nd);
    pref_responses_stat_largePupil_concat=cell(1,nd);
    pref_responses_stat_smallPupil_concat=cell(1,nd);
    pref_allTrials_largePupil_concat =cell(nCon,nSize,nd);
    pref_allTrials_smallPupil_concat =cell(nCon,nSize,nd);
    pupilMeans_concat=nan(nd,3,nSess);
    motorByPupil_concat=nan(nd,2,nSess);
    pupilCounts_concat=nan(nd,2,nSess);
end
resp_keep_concat=cell(1,nd);
resp_max_keep_concat=cell(1,nd);
pref_responses_loc_concat=cell(1,nd);
pref_responses_stat_concat=cell(1,nd);
pref_peak_stat_concat=cell(1,nd);
pref_peak_loc_concat=cell(1,nd);
RIx_concat=cell(1,nd);
dirs_concat=[];
cons_concat=[];
dfof_max_diff_concat=[];
red_fluor_concat=[];
green_fluor_concat=[];
sSupp_concat = {};
supp_status_concat = [];
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
dataTableConat=[];
drug=cell(1,nSess);
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
    load(fullfile(fn_multi,'surr_supp.mat'));
    if doEye == 1
        load(fullfile(fn_multi,'pupilMeans.mat'));
        pupilMeans_concat(:,:,iSess)=pupilMeans;
        motorByPupil_concat(:,:,iSess)=motorByPupil;
    end
    nKeep = size(tc_trial_avrg_stat{post},2);

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
    sSupp_concat = [sSupp_concat;supp_result];
    supp_status_concat = [supp_status_concat;supp_statusArr];
    trial_dfof_concat(iSess,:) = dfofTrial_keep_2days;
    tSize_concat(iSess,:) = tSize_match;
    tCon_concat(iSess,:) = tCon_match;
    tRIx_concat(iSess,:) = RIx;
    % responseByCondProps_concat(:,:,iSess)=responseByCondProps;
    
    clear cons
    
    for id = 1:nd
        tc_trial_avrg_stat_concat{id} =cat(2,tc_trial_avrg_stat_concat{id},tc_trial_avrg_stat{id}(:,:,:,:));
        if doEye == 1
            tc_trial_avrg_stat_largePupil_concat{id} = cat(2,tc_trial_avrg_stat_largePupil_concat{id},tc_trial_avrg_stat_largePupil{id}(:,:,sharedCon,:));
            tc_trial_avrg_stat_smallPupil_concat{id} = cat(2,tc_trial_avrg_stat_smallPupil_concat{id},tc_trial_avrg_stat_smallPupil{id}(:,:,sharedCon,:));
            pref_responses_stat_largePupil_concat{id}=cat(1,pref_responses_stat_largePupil_concat{id},pref_responses_stat_largePupil{id}(:,sharedCon,:));
            pref_responses_stat_smallPupil_concat{id}=cat(1,pref_responses_stat_smallPupil_concat{id},pref_responses_stat_smallPupil{id}(:,sharedCon,:));
        end
        tc_trial_avrg_loc_concat{id} =cat(2,tc_trial_avrg_loc_concat{id},tc_trial_avrg_loc{id}(:,:,sharedCon,:));
        nonPref_trial_avrg_stat_concat{id} =cat(2,nonPref_trial_avrg_stat_concat{id},nonPref_trial_avrg_stat{id}(:,:,sharedCon,:));
        nonPref_trial_avrg_loc_concat{id} =cat(2,nonPref_trial_avrg_loc_concat{id},nonPref_trial_avrg_loc{id}(:,:,sharedCon,:));
        resp_keep_concat{id}=cat(1,resp_keep_concat{id},resp_keep{id});
        resp_max_keep_concat{id}=cat(1,resp_max_keep_concat{id},resp_max_keep{id}(:,sharedCon,:));
        pref_responses_loc_concat{id}=cat(1,pref_responses_loc_concat{id},pref_responses_loc{id}(:,sharedCon,:));
        pref_responses_stat_concat{id}=cat(1,pref_responses_stat_concat{id},pref_responses_stat{id}(:,sharedCon,:));
        pref_peak_stat_concat{id}=cat(1,pref_peak_stat_concat{id},pref_peak_stat{id}(:,sharedCon,:));
        pref_peak_loc_concat{id}=cat(1,pref_peak_loc_concat{id},pref_peak_loc{id}(:,sharedCon,:));
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
            if doEye == 1 
                pref_allTrials_largePupil_concat{i,iSize,id}=[pref_allTrials_largePupil_concat{i,iSize,id},pref_allTrials_largePupil{iCon,iSize,id}];
                pref_allTrials_smallPupil_concat{i,iSize,id}=[pref_allTrials_smallPupil_concat{i,iSize,id},pref_allTrials_smallPupil{iCon,iSize,id}];
            end
            end
        end
        clear meanF i
    end
    dfof_max_diff_concat=cat(1,dfof_max_diff_concat,dfof_max_diff(:,sharedCon,:));
   % green_fluor_concat=cat(2,green_fluor_concat,green_fluor_keep);
   % red_fluor_concat=cat(2,red_fluor_concat,red_fluor_keep);
   
iSess
end

sSupp_concat_mat = cell2mat(sSupp_concat);
sSupp_pre_ind = find(sSupp_concat_mat(:,pre));
NotSS_pre_ind = find(~sSupp_concat_mat(:,pre));

if onlySS == 0 && onlyNOTSS == 0
    red_ind_concat = find(red_concat);
    green_ind_concat = find(green_concat);
end

if onlySS == 1
    red_ind_concat = intersect(find(red_concat),sSupp_pre_ind);
    green_ind_concat = intersect(find(green_concat),sSupp_pre_ind);
elseif onlyNOTSS == 1
    red_ind_concat = intersect(find(red_concat),~sSupp_pre_ind);
    green_ind_concat = intersect(find(green_concat),~sSupp_pre_ind);
end

cons = targetCon;
nSize = length(sizes);
nKeep_total = sum(nKeep_concat);

% clear data_resp_keep data_trial_keep green_ind_keep red_ind_keep conBySize_resp_loc_keep 
% clear conBySize_resp_stat_keep tc_trial_avrg_stat tc_trial_avrg_loc conBySize_resp_loc_match 
% clear conBySize_resp_stat_match red_keep_logical green_keep_logical pref_responses_stat 
% clear pref_responses_loc resp_keep mouse pref_con_keep pref_dir_keep pref_size_keep 
% clear dfof_max_diff dfof_max_diff_raw green_fluor_match green_fluor_keep motorByPupil nKeep noiseCorr
% clear sigCorr nonPref_trial_avrg_loc nonPref_trial_avrg_stat norm_dir_resp_loc norm_dir_resp_stat 
% clear pref_allTrials_largePupil pref_allTrials_smallPupil pref_allTrials_loc pref_allTrials_stat
% clear pref_peak_loc pref_peak_stat pref_responses_stat_smallPupil pref_allTrials_largePupil
% clear pupilMeans red_fluor_match red_fluor_keep RIx tc_trial_avrg_keep_allCond tc_trial_avrg_stat_largePupil
% clear tc_trial_avrg_stat_smallPupil wheel_corr wheel_tc dart_Str expt fullTC_keep wheel_speed
% clear dfofTrial_keep_2days

% cell selection
% find cells that I have running and stationary data for on both days

haveRunning_pre=(sum(squeeze((sum(~isnan(conBySize_resp_loc_concat{pre}),2))),2)==nSize*nCon); %find cells that have a no NAN values for any size (the total of non-NAN should == the number of size)
haveRunning_post=(sum(squeeze((sum(~isnan(conBySize_resp_loc_concat{post}),2))),2)==nSize*nCon);
haveRunning_both= find(haveRunning_pre.* haveRunning_post); %find cells that meet this criteria for both days - now in indices, not logical

haveStat_pre=(sum(squeeze((sum(~isnan(conBySize_resp_stat_concat{pre}),2))),2)==nSize*nCon); %find cells that have a no NAN values for any size (the total of non-NAN should == the number of size)
haveStat_post=(sum(squeeze((sum(~isnan(conBySize_resp_stat_concat{post}),2))),2)==nSize*nCon); 
haveStat_both= find(haveStat_pre.* haveStat_post); %find cells that meet this criteria for both days - now in indices, not logical

runningCells = intersect(haveStat_both, haveRunning_both);


%% further preprocessing
%find cells responsive to a given size
mySize = 1; %1 is the smallest size
respToSizeBothDays = cell(1,nd);
for id = 1:nd
    respToSizeBothDays{id}=sum(squeeze(sum(h_concat{id}(:,:,:,mySize),2)),2); %finding cells that responded this size, at any contrast or direction
end
respToSmall = logical(respToSizeBothDays{pre}+respToSizeBothDays{post}); %to find cells that were responsive to this size on either day

if onlySmall == 1
    red_ind_concat = intersect(red_ind_concat,find(respToSmall));
    green_ind_concat = intersect(green_ind_concat,find(respToSmall));
end

mySize = nSize; %nSize is the largest size
respToSizeBothDays = cell(1,nd);
for id = 1:nd
    respToSizeBothDays{id}=sum(squeeze(sum(h_concat{id}(:,:,:,mySize),2)),2); %finding cells that responded this size, at any contrast or direction
end
respToLarge = logical(respToSizeBothDays{pre}+respToSizeBothDays{post}); %to find cells that were responsive to this size on either day

if onlyLarge == 1
    red_ind_concat = intersect(red_ind_concat,find(respToLarge));
    green_ind_concat = intersect(green_ind_concat,find(respToLarge));
end


responCriteria = cell(1,nd); %cell array that will have indices of cells that meet our response criteria on each day
for id = 1:nd
    responseCheck =sum(squeeze(sum(h_concat{id}(:,:,2:4,nSize),2)),2); %finding cells that responded large size size, within the top three contrasts
    responCriteria{id}=find(logical(responseCheck));
end

includeCells = intersect(responCriteria{pre},find(haveRunning_pre));

% clear haveRunning_pre haveRunning_post haveRunning_both haveStat_both haveStat_pre haveStat_post

% to find the OSI of each cell

OSI_baseline=nan(nKeep_total,1);

dirMean = mean(squeeze(mean(data_resp_concat{pre}(:,:,:,:,1),3)),3);
dirMean(find(dirMean<0))=0;
%left with average response at each dir, averaged over size and contrast,
%for each cell, for the baseline day only
nDir = length(unique((dirs_concat)));

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
% clear cellCounts cellCountsGreen

%to find the cells that have 
runningByCondition = nan(nKeep_total,nCon,nSize);

for iCon = 1:nCon
    for iSize = 1:nSize
        runningPre = ~isnan(conBySize_resp_loc_concat{pre}(:,iCon, iSize));
        runningPost = ~isnan(conBySize_resp_loc_concat{post}(:,iCon, iSize));
        runningByCondition(:,iCon,iSize)=runningPre.*runningPost;
    end
end
%% surround suppression scatter
%tSize_concat tCon_concat tRIx_concat trial_dfof_concat
% these are all nMice * nDay cell arrays
SSIx_concat = [];
runInclude_concat = nan(size(supp_status_concat,1),nd); % count how many cells can be included
statSSIx_concat = [];
locSSIx_concat = [];
dartIx_concat = [];

for sesh = 1:nSess
    if sesh == 1
        cellConcatCounter = 0;
    end
    nCells = size(trial_dfof_concat{sesh,1},1); % take nCells from ref day since it should be the same both days
    Mouse_SSIx = nan(nCells,nd);
    statSSIx_sesh = nan(nCells,nd);
    locSSIx_sesh = nan(nCells,nd);
    include_sm = cell(1,2);
    include_lg = cell(1,2);
    % STD_mat = nan(nCells,2);
    for id = 1:nd
        thisTdfof = trial_dfof_concat{sesh,id};
        thisNTrials = size(thisTdfof,2);
        % prep indexing matrices
        thisConKeep = tCon_concat{sesh,id}(1:thisNTrials) == max(cons); % find only con == 1
        thisSizes_sm = logical((tSize_concat{sesh,id}(1:thisNTrials) == 20) .* thisConKeep);
        thisSizes_lg = logical((tSize_concat{sesh,id}(1:thisNTrials) == 1000) .* thisConKeep);
        thisRunning = tRIx_concat{sesh,id}(1:thisNTrials);
        % check if this session had at least 3 running trial for both sizes
        % if yes, add cell idx (logical) to a matrix and calculate SSIx
        run_sm = logical(thisSizes_sm .* thisRunning);
        run_lg = logical(thisSizes_lg .* thisRunning);
        stat_sm = logical(thisSizes_sm .* ~thisRunning);
        stat_lg = logical(thisSizes_lg .* ~thisRunning);
        insertStart = cellConcatCounter + 1;
        insertEnd = cellConcatCounter + nCells;
        if sum(run_sm) >= 3 && sum(run_lg) >= 3
            thisRunInclude = true(1,nCells);
            runInclude_concat(insertStart:insertEnd,id) = thisRunInclude;
            %SSIx for stat and loc trials separately
            stat_sm_cell_avg = mean(thisTdfof(:,stat_sm),2);
            stat_lg_cell_avg = mean(thisTdfof(:,stat_lg),2);
            stat_sm_cell_avg(stat_sm_cell_avg<0) = 0;
            stat_lg_cell_avg(stat_lg_cell_avg<0) = 0;
            stat_SSIx = (stat_lg_cell_avg - stat_sm_cell_avg)./(stat_lg_cell_avg + stat_sm_cell_avg);

            loc_sm_cell_avg = mean(thisTdfof(:,run_sm),2);
            loc_lg_cell_avg = mean(thisTdfof(:,run_lg),2);
            loc_sm_cell_avg(loc_sm_cell_avg<0) = 0;
            loc_lg_cell_avg(loc_lg_cell_avg<0) = 0;
            loc_SSIx = (loc_lg_cell_avg - loc_sm_cell_avg)./(loc_lg_cell_avg + loc_sm_cell_avg);

            statSSIx_sesh(:,id) = stat_SSIx;
            locSSIx_sesh(:,id) = loc_SSIx;
        else
            thisRunInclude = false(1,nCells);
            runInclude_concat(insertStart:insertEnd,id) = thisRunInclude;
            % put nan in loc-state-depedent SSIx matrices for that session
            stat_SSIx = nan(nCells,1);
            loc_SSIx = nan(nCells,1);
            statSSIx_sesh(:,id) = stat_SSIx;
            locSSIx_sesh(:,id) = loc_SSIx;
        end
        
        % all cells
        cell_avg_dfof_sm = mean(thisTdfof(:,thisSizes_sm),2);
        cell_avg_dfof_lg = mean(thisTdfof(:,thisSizes_lg),2);

        % surround suppression index on both days
        cell_avg_dfof_sm(cell_avg_dfof_sm<0) = 0;
        cell_avg_dfof_lg(cell_avg_dfof_lg<0) = 0;
        SSIx = (cell_avg_dfof_lg - cell_avg_dfof_sm)./(cell_avg_dfof_lg + cell_avg_dfof_sm);
        % SSIx(isnan(SSIx)) = 0;
        Mouse_SSIx(:,id) = SSIx;

        include_sm{id} = stat_sm;
        include_lg{id} = stat_lg;
    end

    %calculate DART effects of small stim trials and large stim trials on
    %the pre day
    %this requires info from both days
    pre_sm = nanmean(trial_dfof_concat{sesh,pre}(:,include_sm{pre}),2);
    post_sm = nanmean(trial_dfof_concat{sesh,post}(:,include_sm{post}),2);
    pre_lg = nanmean(trial_dfof_concat{sesh,pre}(:,include_lg{pre}),2);
    post_lg = nanmean(trial_dfof_concat{sesh,post}(:,include_lg{post}),2);
    pre_sm(pre_sm<0) = 0;
    post_sm(post_sm<0) = 0;
    dartIx_sm = (post_sm - pre_sm) ./ (post_sm + pre_sm);
    pre_lg(pre_lg<0) = 0;
    post_lg(post_lg<0) = 0;
    dartIx_lg = (post_lg - pre_lg) ./ (post_lg + pre_lg);
    dartIx = [dartIx_sm dartIx_lg]; % 1st column small, 2nd column large

    dartIx_concat = [dartIx_concat;dartIx];
    SSIx_concat = [SSIx_concat;Mouse_SSIx];
    cellConcatCounter = cellConcatCounter + nCells;
    statSSIx_concat = [statSSIx_concat;statSSIx_sesh];
    locSSIx_concat = [locSSIx_concat;locSSIx_sesh];

end

%% plot pre vs post DART SSIx
SSIx_green = SSIx_concat(green_ind_concat,:);
SSIx_green_mean_pre = nanmean(SSIx_green(:,pre),1);
SSIx_green_sem_pre = nanstd(SSIx_green(:,pre))./sqrt(sum(~isnan(SSIx_green(:,pre))));
SSIx_green_mean_post = nanmean(SSIx_green(:,post),1);
SSIx_green_sem_post = nanstd(SSIx_green(:,post))./sqrt(sum(~isnan(SSIx_green(:,post))));

SSIx_red = SSIx_concat(red_ind_concat,:);
SSIx_red_mean_pre = nanmean(SSIx_red(:,pre),1);
SSIx_red_sem_pre = nanstd(SSIx_red(:,pre))./sqrt(sum(~isnan(SSIx_red(:,pre))));
SSIx_red_mean_post = nanmean(SSIx_red(:,post),1);
SSIx_red_sem_post = nanstd(SSIx_red(:,post))./sqrt(sum(~isnan(SSIx_red(:,post))));

% pyr
figure;
% subplot(2,1,2)
hold on
title(['HTP- n = ' num2str(min(sum(~isnan(SSIx_green))))])
scatter(SSIx_green(:,pre),SSIx_green(:,post))
line(-1:0.1:1,-1:0.1:1);
eb(1) = errorbar(SSIx_green_mean_pre,SSIx_green_mean_post,SSIx_green_sem_pre, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(SSIx_green_mean_pre,SSIx_green_mean_post,SSIx_green_sem_post, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'k', 'LineWidth', 1.5)
xlabel('Pre-DART')
ylabel('Post-DART')
hold off

% htp
figure;
% subplot(2,1,2)
hold on
title(['HTP+ n = ' num2str(min(sum(~isnan(SSIx_red))))])
scatter(SSIx_red(:,pre),SSIx_red(:,post))
line(-1:0.1:1,-1:0.1:1);
eb(1) = errorbar(SSIx_red_mean_pre,SSIx_red_mean_post,SSIx_red_sem_pre, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(SSIx_red_mean_pre,SSIx_red_mean_post,SSIx_red_sem_post, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'k', 'LineWidth', 1.5)
xlabel('Pre-DART')
ylabel('Post-DART')
hold off

%% plot stat vs. loc SSIx
statSSIx_green = statSSIx_concat(green_ind_concat,:);
statSSIx_green_mean_pre = nanmean(statSSIx_green(:,pre),1);
statSSIx_green_sem_pre = nanstd(statSSIx_green(:,pre))./sqrt(sum(~isnan(statSSIx_green(:,pre))));
statSSIx_green_mean_post = nanmean(statSSIx_green(:,post),1);
statSSIx_green_sem_post = nanstd(statSSIx_green(:,post))./sqrt(sum(~isnan(statSSIx_green(:,post))));

statSSIx_red = statSSIx_concat(red_ind_concat,:);
statSSIx_red_mean_pre = nanmean(statSSIx_red(:,pre),1);
statSSIx_red_sem_pre = nanstd(statSSIx_red(:,pre))./sqrt(sum(~isnan(statSSIx_red(:,pre))));
statSSIx_red_mean_post = nanmean(statSSIx_red(:,post),1);
statSSIx_red_sem_post = nanstd(statSSIx_red(:,post))./sqrt(sum(~isnan(statSSIx_red(:,post))));

locSSIx_green = locSSIx_concat(green_ind_concat,:);
locSSIx_green_mean_pre = nanmean(locSSIx_green(:,pre),1);
locSSIx_green_sem_pre = nanstd(locSSIx_green(:,pre))./sqrt(sum(~isnan(locSSIx_green(:,pre))));
locSSIx_green_mean_post = nanmean(locSSIx_green(:,post),1);
locSSIx_green_sem_post = nanstd(locSSIx_green(:,post))./sqrt(sum(~isnan(locSSIx_green(:,post))));

locSSIx_red = locSSIx_concat(red_ind_concat,:);
locSSIx_red_mean_pre = nanmean(locSSIx_red(:,pre),1);
locSSIx_red_sem_pre = nanstd(locSSIx_red(:,pre))./sqrt(sum(~isnan(locSSIx_red(:,pre))));
locSSIx_red_mean_post = nanmean(locSSIx_red(:,post),1);
locSSIx_red_sem_post = nanstd(locSSIx_red(:,post))./sqrt(sum(~isnan(locSSIx_red(:,post))));

% check how many cells are actually included (exclude nan SSIx)
boo_stat_green = ~isnan(statSSIx_green);
boo_loc_green = ~isnan(locSSIx_green);
boo_stat_red = ~isnan(statSSIx_red);
boo_loc_red = ~isnan(locSSIx_red);

% pre dart plots
N_pre_green = sum(boo_stat_green(:,pre) .* boo_loc_green(:,pre));
N_pre_red = sum(boo_stat_red(:,pre) .* boo_loc_red(:,pre));

figure;
hold on
title(['HTP- n = ' num2str(N_pre_green)])
scatter(statSSIx_green(:,pre),locSSIx_green(:,pre))
line(-1:0.1:1,-1:0.1:1);
eb(1) = errorbar(statSSIx_green_mean_pre,locSSIx_green_mean_pre,statSSIx_green_sem_pre, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(statSSIx_green_mean_pre,locSSIx_green_mean_pre,locSSIx_green_sem_pre, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'k', 'LineWidth', 1.5)
xlabel('stat')
ylabel('loc')
hold off

figure;
hold on
title(['HTP+ n = ' num2str(N_pre_red)])
scatter(statSSIx_red(:,pre),locSSIx_red(:,pre))
line(-1:0.1:1,-1:0.1:1);
eb(1) = errorbar(statSSIx_red_mean_pre,locSSIx_red_mean_pre,statSSIx_red_sem_pre, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(statSSIx_red_mean_pre,locSSIx_red_mean_pre,locSSIx_red_sem_pre, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'k', 'LineWidth', 1.5)
xlabel('stat')
ylabel('loc')
hold off

% post dart plots
N_post_green = sum(boo_stat_green(:,post) .* boo_loc_green(:,post));
N_post_red = sum(boo_stat_red(:,post) .* boo_loc_red(:,post));

figure;
hold on
title(['HTP- n = ' num2str(N_post_green)])
scatter(statSSIx_green(:,post),locSSIx_green(:,post))
line(-1:0.1:1,-1:0.1:1);
eb(1) = errorbar(statSSIx_green_mean_post,locSSIx_green_mean_post,statSSIx_green_sem_post, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(statSSIx_green_mean_post,locSSIx_green_mean_post,locSSIx_green_sem_post, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'k', 'LineWidth', 1.5)
xlabel('stat')
ylabel('loc')
hold off

figure;
hold on
title(['HTP+ n = ' num2str(N_post_red)])
scatter(statSSIx_red(:,post),locSSIx_red(:,post))
line(-1:0.1:1,-1:0.1:1);
eb(1) = errorbar(statSSIx_red_mean_post,locSSIx_red_mean_post,statSSIx_red_sem_post, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(statSSIx_red_mean_post,locSSIx_red_mean_post,locSSIx_red_sem_post, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'k', 'LineWidth', 1.5)
xlabel('stat')
ylabel('loc')
hold off

%% DART effect idx vs suppression idx
% DART effect of either large or small trials; highest contrast, only
% stationary, again pre-dart SSIx

% dartIx_concat 1st column small, 2nd column large
dIx_sm_green = dartIx_concat(green_ind_concat,1);
sm_green_avg = nanmean(dIx_sm_green);
sm_green_sem = nanstd(dIx_sm_green)./sqrt(sum(~isnan(dIx_sm_green)));

dIx_lg_green = dartIx_concat(green_ind_concat,2);
lg_green_avg = nanmean(dIx_lg_green);
lg_green_sem = nanstd(dIx_lg_green)./sqrt(sum(~isnan(dIx_lg_green)));

dIx_sm_red = dartIx_concat(red_ind_concat,1);
sm_red_avg = nanmean(dIx_sm_red);
sm_red_sem = nanstd(dIx_sm_red)./sqrt(sum(~isnan(dIx_sm_red)));

dIx_lg_red = dartIx_concat(red_ind_concat,2);
lg_red_avg = nanmean(dIx_lg_red);
lg_red_sem = nanstd(dIx_lg_red)./sqrt(sum(~isnan(dIx_lg_red)));

% calculate SSIx of the pre day as done above
SSIx_green = SSIx_concat(green_ind_concat,:);
SSIx_green_mean_pre = nanmean(SSIx_green(:,pre),1);
SSIx_green_sem_pre = nanstd(SSIx_green(:,pre))./sqrt(sum(~isnan(SSIx_green(:,pre))));

SSIx_red = SSIx_concat(red_ind_concat,:);
SSIx_red_mean_pre = nanmean(SSIx_red(:,pre),1);
SSIx_red_sem_pre = nanstd(SSIx_red(:,pre))./sqrt(sum(~isnan(SSIx_red(:,pre))));

% calculate actual nCells

n_green_sm = sum(~isnan(SSIx_green(:,pre)) .* ~isnan(dIx_sm_green));
n_green_lg = sum(~isnan(SSIx_green(:,pre)) .* ~isnan(dIx_lg_green));
n_red_sm = sum(~isnan(SSIx_red(:,pre)) .* ~isnan(dIx_sm_red));
n_red_lg = sum(~isnan(SSIx_red(:,pre)) .* ~isnan(dIx_lg_red));

% regression model
% mask1 = ~isnan(SSIx_green(:,pre)) & ~isnan(dIx_sm_green);
% model1 = polyfit(SSIx_green(mask1,pre),dIx_sm_green(mask1),1);

model1 = fitlm(SSIx_green(:,pre),dIx_sm_green);
figure
hand1 = plot(model1)

model2 = fitlm(SSIx_red(:,pre),dIx_sm_red);
figure
hand2 = plot(model2)

model3 = fitlm(SSIx_green(:,pre),dIx_lg_green);
figure
hand3 = plot(model3)

model4 = fitlm(SSIx_red(:,pre),dIx_lg_red);
figure
hand4 = plot(model4)

% plots
figure;
hold on
title(['Small Stim, HTP- n = ' num2str(n_green_sm)])
scatter(SSIx_green(:,pre),dIx_sm_green)
% line(-1:0.1:1,-1:0.1:1);
% eb(1) = errorbar(SSIx_green_mean_pre,sm_green_avg,SSIx_green_sem_pre, 'horizontal', 'LineStyle', 'none');
% eb(2) = errorbar(SSIx_green_mean_pre,sm_green_avg,sm_green_sem, 'vertical', 'LineStyle', 'none');
% set(eb, 'color', 'k', 'LineWidth', 1.5)
lsline
xlabel('SSIx')
ylabel('DART Index')
hold off

figure;
hold on
title(['Small Stim, HTP+ n = ' num2str(n_red_sm)])
scatter(SSIx_red(:,pre),dIx_sm_red)
% line(-1:0.1:1,-1:0.1:1);
% eb(1) = errorbar(SSIx_red_mean_pre,sm_red_avg,SSIx_red_sem_pre, 'horizontal', 'LineStyle', 'none');
% eb(2) = errorbar(SSIx_red_mean_pre,sm_red_avg,sm_red_sem, 'vertical', 'LineStyle', 'none');
% set(eb, 'color', 'k', 'LineWidth', 1.5)
lsline
xlabel('SSIx')
ylabel('DART Index')
hold off

figure;
hold on
title(['Large Stim, HTP- n = ' num2str(n_green_lg)])
scatter(SSIx_green(:,pre),dIx_lg_green)
% line(-1:0.1:1,-1:0.1:1);
% eb(1) = errorbar(SSIx_green_mean_pre,lg_green_avg,SSIx_green_sem_pre, 'horizontal', 'LineStyle', 'none');
% eb(2) = errorbar(SSIx_green_mean_pre,lg_green_avg,lg_green_sem, 'vertical', 'LineStyle', 'none');
% set(eb, 'color', 'k', 'LineWidth', 1.5)
lsline
xlabel('SSIx')
ylabel('DART Index')
hold off

figure;
hold on
title(['Large Stim, HTP+ n = ' num2str(n_red_lg)])
scatter(SSIx_red(:,pre),dIx_lg_red)
% line(-1:0.1:1,-1:0.1:1);
% eb(1) = errorbar(SSIx_red_mean_pre,lg_red_avg,SSIx_red_sem_pre, 'horizontal', 'LineStyle', 'none');
% eb(2) = errorbar(SSIx_red_mean_pre,lg_red_avg,lg_red_sem, 'vertical', 'LineStyle', 'none');
% set(eb, 'color', 'k', 'LineWidth', 1.5)
lsline
xlabel('SSIx')
ylabel('DART Index')
hold off


%% plot stationary timecourses for all cells

% make figure with se shaded, one figure per contrast - stationary

tc_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat = cell(1,nd); %same for red
tc_green_se_stat = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat = cell(1,nd); %same for red


for id = 1:nd
  for iCon = 1:nCon
      for iSize = 1:nSize
        
        tc_green_avrg_stat{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_stat_concat{id}(:,green_ind_concat,iCon,iSize),2);
        green_std=nanstd(tc_trial_avrg_stat_concat{id}(:,green_ind_concat,iCon,iSize),[],2);
        tc_green_se_stat{id}(:,iCon,iSize)=green_std/sqrt(length(green_ind_concat));
        
        tc_red_avrg_stat{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_stat_concat{id}(:,red_ind_concat,iCon,iSize),2);
        red_std=nanstd(tc_trial_avrg_stat_concat{id}(:,red_ind_concat,iCon,iSize),[],2);
        tc_red_se_stat{id}(:,iCon,iSize)=red_std/sqrt(length(red_ind_concat));
        
        
        clear green_std red_std
      end 
    end
end

z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
    for iSize = 1:nSize
    figure
    subplot(1,2,1) %for the first day
    
    ylim([-.02 .1]);
    hold on
    shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon,iSize),tc_green_se_stat{pre}(:,iCon,iSize),'--k');
    hold on
    shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon,iSize),tc_green_se_stat{post}(:,iCon,iSize),'--b','transparent');
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    title(['-HTP',' n = ', num2str(length(green_ind_concat))])
    
    ylabel('dF/F') 
    xlabel('s') 
    set(gca,'XColor', 'none','YColor','none')
    
    
    subplot(1,2,2) %+HTP
    shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon,iSize),tc_red_se_stat{pre}(:,iCon,iSize),'k');
    hold on
    shadedErrorBar(t,tc_red_avrg_stat{post}(:,iCon,iSize),tc_red_se_stat{post}(:,iCon,iSize),'b');
    ylim([-.02 .1]);
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    ylabel('dF/F') 
    xlabel('s') 
    title(['+HTP',' n = ', num2str(length(red_ind_concat))])
    
    x0=5;
    y0=5;
    width=4;
    height=3;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    set(gca,'XColor', 'none','YColor','none')
    
    sgtitle(['stationary, con ' num2str(cons(iCon)) ' size ' num2str(sizes(iSize))])
    
    print(fullfile(fnout,[num2str(cons(iCon)) '_' num2str(sizes(iSize)) '_stat_cellType_timecourses.pdf']),'-dpdf');
    end
end 

%% plot stationary timecourses for cells matched across behavioral state within each stim condition

% make figure with se shaded, one figure per contrast - stationary

tc_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat = cell(1,nd); %same for red
tc_green_se_stat = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat = cell(1,nd); %same for red


nGreen = nan(nCon,nSize);
nRed = nan(nCon,nSize);

for id = 1:nd
  for iCon = 1:nCon
      for iSize = 1:nSize
        theseGreenCells = intersect(green_ind_concat, find(runningByCondition(:,iCon,iSize)));
        theseRedCells = intersect(red_ind_concat, find(runningByCondition(:,iCon,iSize)));
        nGreen(iCon, iSize)=length(theseGreenCells);
        nRed(iCon, iSize)=length(theseRedCells);

        tc_green_avrg_stat{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_stat_concat{id}(:,theseGreenCells,iCon,iSize),2);
        green_std=nanstd(tc_trial_avrg_stat_concat{id}(:,theseGreenCells,iCon,iSize),[],2);
        tc_green_se_stat{id}(:,iCon,iSize)=green_std/sqrt(length(theseGreenCells));
        
        tc_red_avrg_stat{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_stat_concat{id}(:,theseRedCells,iCon,iSize),2);
        red_std=nanstd(tc_trial_avrg_stat_concat{id}(:,theseRedCells,iCon,iSize),[],2);
        tc_red_se_stat{id}(:,iCon,iSize)=red_std/sqrt(length(theseRedCells));
        
        
        clear green_std red_std
      end 
   end
end

z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
    for iSize = 1:nSize
    figure
    subplot(1,2,1) %for the first day
    
    
    
    ylim([-.02 .1]);
    hold on
    shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon,iSize),tc_green_se_stat{pre}(:,iCon,iSize),'--k');
    hold on
    shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon,iSize),tc_green_se_stat{post}(:,iCon,iSize),'--b','transparent');
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    title(['-HTP',' n = ', num2str(nGreen(iCon,iSize))])
    
    ylabel('dF/F') 
    xlabel('s') 
    set(gca,'XColor', 'none','YColor','none')
    
    
    subplot(1,2,2) %+HTP
    shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon,iSize),tc_red_se_stat{pre}(:,iCon,iSize),'k');
    hold on
    shadedErrorBar(t,tc_red_avrg_stat{post}(:,iCon,iSize),tc_red_se_stat{post}(:,iCon,iSize),'b');
    ylim([-.02 .1]);
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    ylabel('dF/F') 
    xlabel('s') 
    title(['+HTP',' n = ', num2str(nRed(iCon,iSize))])
    
    x0=5;
    y0=5;
    width=4;
    height=3;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    set(gca,'XColor', 'none','YColor','none')
    
    sgtitle(['stationary, con ' num2str(cons(iCon)) ' size ' num2str(sizes(iSize))])
    
    print(fullfile(fnout,[num2str(cons(iCon)) '_' num2str(sizes(iSize)) 'matched_stat_cellType_timecourses.pdf']),'-dpdf');
    end
end 

% plots for running trials with cell matched across behavioral state

tc_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_loc = cell(1,nd); %same for red
tc_green_se_loc = cell(1,nd); %this will be the se across all green cells
tc_red_se_loc = cell(1,nd); %same for red


for id = 1:nd
  for iCon = 1:nCon
      for iSize = 1:nSize
        theseGreenCells = intersect(green_ind_concat, find(runningByCondition(:,iCon,iSize)));
        theseRedCells = intersect(red_ind_concat, find(runningByCondition(:,iCon,iSize)));

        tc_green_avrg_loc{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_loc_concat{id}(:,theseGreenCells,iCon,iSize),2);
        green_std=nanstd(tc_trial_avrg_loc_concat{id}(:,theseGreenCells,iCon,iSize),[],2);
        tc_green_se_loc{id}(:,iCon,iSize)=green_std/sqrt(length(theseGreenCells));
        
        tc_red_avrg_loc{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_loc_concat{id}(:,theseRedCells,iCon,iSize),2);
        red_std=nanstd(tc_trial_avrg_loc_concat{id}(:,theseRedCells,iCon,iSize),[],2);
        tc_red_se_loc{id}(:,iCon,iSize)=red_std/sqrt(length(theseRedCells));
        
        clear green_std red_std
      end 
    end
end

z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_loc{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
    for iSize = 1:nSize
    figure
    subplot(1,2,1) %for the first day
    
    
    
    ylim([-.02 .2]);
    hold on
    shadedErrorBar(t,tc_green_avrg_loc{pre}(:,iCon,iSize),tc_green_se_loc{pre}(:,iCon,iSize),'--k');
    hold on
    shadedErrorBar(t,tc_green_avrg_loc{post}(:,iCon,iSize),tc_green_se_loc{post}(:,iCon,iSize),'--b','transparent');
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    title(['-HTP',' n = ', num2str(nGreen(iCon,iSize))])
    
    ylabel('dF/F') 
    xlabel('s') 
    set(gca,'XColor', 'none','YColor','none')
    
    
    subplot(1,2,2) %+HTP
    shadedErrorBar(t,tc_red_avrg_loc{pre}(:,iCon,iSize),tc_red_se_loc{pre}(:,iCon,iSize),'k');
    hold on
    shadedErrorBar(t,tc_red_avrg_loc{post}(:,iCon,iSize),tc_red_se_loc{post}(:,iCon,iSize),'b');
    ylim([-.02 .2]);
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    ylabel('dF/F') 
    xlabel('s') 
    title(['+HTP',' n = ', num2str(nRed(iCon,iSize))])
    
    x0=5;
    y0=5;
    width=4;
    height=3;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    set(gca,'XColor', 'none','YColor','none')
    
    sgtitle(['running, con ' num2str(cons(iCon)) ' size ' num2str(sizes(iSize))])
    
    print(fullfile(fnout,[num2str(cons(iCon)) '_' num2str(sizes(iSize)) 'matched_loc_cellType_timecourses.pdf']),'-dpdf');
    end
end 


%% population size tuning - averaged over contrast

%errorbar for stat resp and loc resp vs size, where error is across mice

%get average for each day

sizeResp_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
sizeResp_red_avrg_stat = cell(1,nd); %same for red
sizeResp_green_se_stat = cell(1,nd); %this will be the se across all green cells
sizeResp_red_se_stat = cell(1,nd); %same for red



for id = 1:nd
    green_data=squeeze(mean(conBySize_resp_stat_concat{id}(statGreen,:,:),2));%pulling the green cells and averaging over contrast
    sizeResp_green_avrg_stat{id}=nanmean(green_data,1);
    green_std=nanstd(green_data,1);
    sizeResp_green_se_stat{id}=green_std/sqrt(length(statGreen));
    
    
    red_data=squeeze(mean(conBySize_resp_stat_concat{id}(statRed,:,:),2));%pulling the red cells and averaging over contrast
    sizeResp_red_avrg_stat{id}=nanmean(red_data,1);
    red_std=nanstd(red_data,1);
    sizeResp_red_se_stat{id}=red_std/sqrt(length(statRed));
    
    clear green_std red_std green_data red_data
 
end



sizeResp_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
sizeResp_red_avrg_loc = cell(1,nd); %same for red
sizeResp_green_se_loc = cell(1,nd); %this will be the se across all green cells
sizeResp_red_se_loc = cell(1,nd); %same for red



for id = 1:nd
    green_data=squeeze(mean(conBySize_resp_loc_concat{id}(runningGreen,:,:),2));%pulling the green cells and averaging over contrast
    sizeResp_green_avrg_loc{id}=nanmean(green_data,1);
    green_std=nanstd(green_data,1);
    sizeResp_green_se_loc{id}=green_std/sqrt(length(runningGreen));
    
    red_data=squeeze(mean(conBySize_resp_loc_concat{id}(runningRed,:,:),2));%pulling the red cells and averaging over contrast
    sizeResp_red_avrg_loc{id}=nanmean(red_data,1);
    red_std=nanstd(red_data,1);
    sizeResp_red_se_loc{id}=red_std/sqrt(length(runningRed));
    
    clear green_std red_std green_data red_data
 
end



figure
subplot(2,2,1) %for the first day
errorbar(sizes,sizeResp_green_avrg_stat{pre},sizeResp_green_se_stat{pre},'k');
hold on
errorbar(sizes,sizeResp_green_avrg_stat{post},sizeResp_green_se_stat{post},'b');
title(['Stationary -HTP, DART',' n = ', num2str(length(statGreen))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off
ylim([-.025 .05])

subplot(2,2,2) %for the second day
errorbar(sizes,sizeResp_red_avrg_stat{pre},sizeResp_red_se_stat{pre},'k');
hold on
errorbar(sizes,sizeResp_red_avrg_stat{post},sizeResp_red_se_stat{post},'b');
title(['Stationary +HTP',' n = ', num2str(length(statRed))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off
ylim([-.025 .05])

subplot(2,2,3) %for the first day
errorbar(sizes,sizeResp_green_avrg_loc{pre},sizeResp_green_se_loc{pre},'k');
hold on
errorbar(sizes,sizeResp_green_avrg_loc{post},sizeResp_green_se_loc{post},'b');
title(['Running -HTP',' n = ', num2str(length(runningGreen))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off
ylim([0 .15])

subplot(2,2,4) %for the second day
errorbar(sizes,sizeResp_red_avrg_loc{pre},sizeResp_red_se_loc{pre},'k');
hold on
errorbar(sizes,sizeResp_red_avrg_loc{post},sizeResp_red_se_loc{post},'b');
title(['Running +HTP',' n = ', num2str(length(runningRed))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off
ylim([0 .1])

x0=5;
y0=5;
width=6;
height=6;
set(gcf,'units','inches','position',[x0,y0,width,height])



sgtitle(['population size tuning' ])

print(fullfile(fnout,['sizeTuning.pdf']),'-dpdf');

%% plotting size tuning for running vs. stationary, seperated by drug condition

figure
subplot(2,2,1) %for the first day
errorbar(sizes,sizeResp_green_avrg_stat{pre},sizeResp_green_se_stat{pre},'k');
hold on
errorbar(sizes,sizeResp_green_avrg_loc{pre},sizeResp_green_se_loc{pre},'m');
title(['Pre-DART -HTP',' n = ', num2str(length(runningGreen))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off

subplot(2,2,2) %for the second day
errorbar(sizes,sizeResp_red_avrg_stat{pre},sizeResp_red_se_stat{pre},'k');
hold on
errorbar(sizes,sizeResp_red_avrg_loc{pre},sizeResp_red_se_loc{pre},'m');
title(['Pre-DART +HTP',' n = ', num2str(length(runningRed))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off

subplot(2,2,3) %for the first day

errorbar(sizes,sizeResp_green_avrg_stat{post},sizeResp_green_se_stat{post},'k');
hold on
errorbar(sizes,sizeResp_green_avrg_loc{post},sizeResp_green_se_loc{post},'m');
title(['Post-DART -HTP',' n = ', num2str(length(runningGreen))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off

subplot(2,2,4) %for the second day
errorbar(sizes,sizeResp_red_avrg_stat{post},sizeResp_red_se_stat{post},'k');
hold on
errorbar(sizes,sizeResp_red_avrg_loc{post},sizeResp_red_se_loc{post},'m');
title(['Post-DART +HTP',' n = ', num2str(length(runningRed))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])



sgtitle(['population size tuning' ])

print(fullfile(fnout,['sizeTuningVsBehState.pdf']),'-dpdf');
%% contrast response
ymin= -0.02;
ymax=.06;
% errorbar for stat resp and loc resp vs size, where error is across mice
conResp_green_avrg_stat = cell(nSize,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_stat = cell(nSize,nd); %same for red
conResp_green_se_stat = cell(nSize,nd); %this will be the se across all green cells
conResp_red_se_stat = cell(nSize,nd); %same for red

consForPlotting = [12.5 25 50 100];
for id = 1:nd
    for iSize = 1:nSize
        green_data=conBySize_resp_stat_concat{id}(green_ind_concat,:,iSize);%pulling the green cells at this size
        conResp_green_avrg_stat{id}(iSize,:)=nanmean(green_data,1);
        green_std=nanstd(green_data,[],1);
        conResp_green_se_stat{id}(iSize,:)=green_std/sqrt(length(green_ind_concat));
        
        red_data=conBySize_resp_stat_concat{id}(red_ind_concat,:,iSize);%pulling the red cells at this size
        conResp_red_avrg_stat{id}(iSize,:)=nanmean(red_data,1);
        red_std=nanstd(red_data,[],1);
        conResp_red_se_stat{id}(iSize,:)=red_std/sqrt(length(red_ind_concat));
        
        clear green_std red_std green_data red_data
    end
end



figure
subplot(2,2,1) %for the first size, all contrasts
errorbar(consForPlotting,conResp_green_avrg_stat{pre}(1,:),conResp_green_se_stat{pre}(1,:),'--k');
hold on
errorbar(consForPlotting,conResp_green_avrg_stat{post}(1,:),conResp_green_se_stat{post}(1,:),'--b');
title(['Pyr n = ' , num2str(length(runningGreen))])
ylabel('dF/F, 20 deg') 
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax])
xlim([0 110])

subplot(2,2,3) %for the second size, all contrasts
errorbar(consForPlotting,conResp_green_avrg_stat{pre}(2,:),conResp_green_se_stat{pre}(2,:),'--k');
hold on
errorbar(consForPlotting,conResp_green_avrg_stat{post}(2,:),conResp_green_se_stat{post}(2,:),'--b');
ylabel('dF/F, Fullfield') 
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax])
xlim([0 110])

subplot(2,2,2) %for the first day
errorbar(consForPlotting,conResp_red_avrg_stat{pre}(1,:),conResp_red_se_stat{pre}(1,:),'k');
hold on
errorbar(consForPlotting,conResp_red_avrg_stat{post}(1,:),conResp_red_se_stat{post}(1,:),'b');
title(['VIP n = ' , num2str(length(runningRed))])
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax])
xlim([0 110])

subplot(2,2,4) %for the first day
errorbar(consForPlotting,conResp_red_avrg_stat{pre}(2,:),conResp_red_se_stat{pre}(2,:),'k');
hold on
errorbar(consForPlotting,conResp_red_avrg_stat{post}(2,:),conResp_red_se_stat{post}(2,:),'b');
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax])
xlim([0 110])


han=axes('visible','off'); 
han.XLabel.Visible='on';
xlabel(han,'Contrast (%)');

x0=5;
y0=5;
width=6;
height=4;
set(gcf,'units','inches','position',[x0,y0,width,height])
sgtitle('Stationary')
print(fullfile(fnout,'contrastTuning.pdf'),'-dpdf');

%for running 
ymin= 0;
ymax=.14;
% contrast response running
% errorbar for loc resp and loc resp vs size, where error is across mice
conResp_green_avrg_loc = cell(nSize,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_loc = cell(nSize,nd); %same for red
conResp_green_se_loc = cell(nSize,nd); %this will be the se across all green cells
conResp_red_se_loc = cell(nSize,nd); %same for red

consForPlotting = [12.5 25 50 100];
for id = 1:nd
    for iSize = 1:nSize
        green_data=conBySize_resp_loc_concat{id}(green_ind_concat,:,iSize);%pulling the green cells at this size
        conResp_green_avrg_loc{id}(iSize,:)=nanmean(green_data,1);
        green_std=nanstd(green_data,1);
        conResp_green_se_loc{id}(iSize,:)=green_std/sqrt(length(green_ind_concat));
        
        red_data=conBySize_resp_loc_concat{id}(red_ind_concat,:,iSize);%pulling the red cells at this size
        conResp_red_avrg_loc{id}(iSize,:)=nanmean(red_data,1);
        red_std=nanstd(red_data,1);
        conResp_red_se_loc{id}(iSize,:)=red_std/sqrt(length(red_ind_concat));
        
        % clear green_std red_std green_data red_data
    end
end

figure
subplot(2,2,1) %for the first day
errorbar(consForPlotting,conResp_green_avrg_loc{pre}(1,:),conResp_green_se_loc{pre}(1,:),'--k');
hold on
errorbar(consForPlotting,conResp_green_avrg_loc{post}(1,:),conResp_green_se_loc{post}(1,:),'--b');
title(['Pyr n = ' , num2str(length(runningGreen))])
ylabel('dF/F, 20 deg') 
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax])
xlim([0 110])

subplot(2,2,3) %for the first day
errorbar(consForPlotting,conResp_green_avrg_loc{pre}(2,:),conResp_green_se_loc{pre}(2,:),'--k');
hold on
errorbar(consForPlotting,conResp_green_avrg_loc{post}(2,:),conResp_green_se_loc{post}(2,:),'--b');
ylabel('dF/F, Fullfield') 
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax])
xlim([0 110])

subplot(2,2,2) %for the first day
errorbar(consForPlotting,conResp_red_avrg_loc{pre}(1,:),conResp_red_se_loc{pre}(1,:),'k');
hold on
errorbar(consForPlotting,conResp_red_avrg_loc{post}(1,:),conResp_red_se_loc{post}(1,:),'b');
title(['VIP n = ' , num2str(length(runningRed))])
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax])
xlim([0 110])

subplot(2,2,4) %for the first day
errorbar(consForPlotting,conResp_red_avrg_loc{pre}(2,:),conResp_red_se_loc{pre}(2,:),'k');
hold on
errorbar(consForPlotting,conResp_red_avrg_loc{post}(2,:),conResp_red_se_loc{post}(2,:),'b');
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax])
xlim([0 110])


han=axes('visible','off'); 
han.XLabel.Visible='on';
xlabel(han,'Contrast (%)');

x0=5;
y0=5;
width=6;
height=4;
set(gcf,'units','inches','position',[x0,y0,width,height])
sgtitle('Running')


print(fullfile(fnout,['contrastTuningRunning.pdf']),'-dpdf');


%% linear direction plot 
if nDir > 1
if length(dirs) ==8
    dirs_for_plotting=dirs-180;
elseif length(dirs) ==4
    dirs_for_plotting =dirs;
end

yMin = -.01;
yMax = .1;
yMax2 = .15;
iSize=2; %set to whichever size you want


green_dir_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
red_dir_avrg_stat = cell(1,nd); %same for red
green_dir_se_stat = cell(1,nd); %this will be the se across all green cells
red_dir_se_stat = cell(1,nd); %same for red

for iCon = 1:nCon
    for id = 1:nd
       
        green_dir_avrg_stat{id}=nanmean(nanmean(norm_dir_resp_stat_concat{id}(statGreen,:,iCon,iSize),4),1);
        green_std=nanstd(nanmean(norm_dir_resp_stat_concat{id}(statGreen,:,iCon,iSize),4),[],1);
        green_dir_se_stat{id}=green_std/sqrt(length(statGreen));
        green_dir_avrg_stat{id}=circshift(green_dir_avrg_stat{id},4);
        green_dir_se_stat{id}=circshift(green_dir_se_stat{id},4);
        
        red_dir_avrg_stat{id}=nanmean(nanmean(norm_dir_resp_stat_concat{id}(statRed,:,iCon,iSize),4),1);
        red_std=nanstd(nanmean(norm_dir_resp_stat_concat{id}(statRed,:,iCon,iSize),4),[],1);
        red_dir_se_stat{id}=red_std/sqrt(length(statRed));
        red_dir_avrg_stat{id}=circshift(red_dir_avrg_stat{id},4);
        red_dir_se_stat{id}=circshift(red_dir_se_stat{id},4);
        clear green_std red_std
        
    end
    
    
    
    green_dir_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
    red_dir_avrg_loc = cell(1,nd); %same for red
    green_dir_se_loc = cell(1,nd); %this will be the se across all green cells
    red_dir_se_loc = cell(1,nd); %same for red
    
    for id = 1:nd
       
        green_dir_avrg_loc{id}=nanmean(nanmean(norm_dir_resp_loc_concat{id}(runningGreen,:,iCon,iSize),4),1);
        green_std=nanstd(nanmean(norm_dir_resp_loc_concat{id}(runningGreen,:,iCon,iSize),4),[],1);
        green_dir_se_loc{id}=green_std/sqrt(length(runningGreen));
        green_dir_avrg_loc{id}=circshift(green_dir_avrg_loc{id},4);
        green_dir_se_loc{id}=circshift(green_dir_se_loc{id},4);
        
        red_dir_avrg_loc{id}=nanmean(nanmean(norm_dir_resp_loc_concat{id}(runningRed,:,iCon,iSize),4),1);
        red_std=nanstd(nanmean(norm_dir_resp_loc_concat{id}(runningRed,:,iCon,iSize),4),[],1);
        red_dir_se_loc{id}=red_std/sqrt(length(runningRed));
        red_dir_avrg_loc{id}=circshift(red_dir_avrg_loc{id},4);
        red_dir_se_loc{id}=circshift(red_dir_se_loc{id},4);
        clear green_std red_std
        
    end
    
    
    
    figure
    subplot(2,2,1)
    errorbar(dirs_for_plotting,green_dir_avrg_stat{pre},green_dir_se_stat{pre},'--k')
    hold on
    errorbar(dirs_for_plotting,green_dir_avrg_stat{post},green_dir_se_stat{post},'--b')
    title('Stationary, Pyr')
    set(gca, 'TickDir', 'out')
    xticks(dirs)
    axis square
    box off
    ylabel('dF/F')
    ylim([yMin yMax])
    xlim([-10 140])
    
    subplot(2,2,2)
    errorbar(dirs_for_plotting,red_dir_avrg_stat{pre},red_dir_se_stat{pre},'k')
    hold on
    errorbar(dirs_for_plotting,red_dir_avrg_stat{post},red_dir_se_stat{post},'b')
   title('Stationary, VIP')
    set(gca, 'TickDir', 'out')
    axis square
    box off
    xticks(dirs)
    ylim([yMin yMax])
    xlim([-10 140])
    
    subplot(2,2,3)
    errorbar(dirs_for_plotting,green_dir_avrg_loc{pre},green_dir_se_loc{pre},'--k')
    hold on
    errorbar(dirs_for_plotting,green_dir_avrg_loc{post},green_dir_se_loc{post},'--b')
    title('Running, Pyr')
    set(gca, 'TickDir', 'out')
    axis square
    box off
    xlabel('normalized direction')
    ylabel('dF/F')
    xticks(dirs)
    ylim([yMin yMax2])
    xlim([-10 140])
    
    subplot(2,2,4)
    errorbar(dirs_for_plotting,red_dir_avrg_loc{pre},red_dir_se_loc{pre},'k')
    hold on
    errorbar(dirs_for_plotting,red_dir_avrg_loc{post},red_dir_se_loc{post},'b')
    title('Running, VIP')
    set(gca, 'TickDir', 'out')
    axis square
    box off
    xlabel('normalized direction')
    xticks(dirs)
    ylim([yMin yMax2])
    xlim([-10 140])

    sgtitle(['Normalized direction tuning ',num2str(cons(iCon))])
    
    
    print(fullfile(fnout,[num2str(cons(iCon)),'dirTuning.pdf']),'-dpdf','-bestfit')
end
end
%% change in pref direction
pref_dir_change=abs(pref_dir_concat{pre}-pref_dir_concat{post});
figure;subplot(1,2,1)
polarhistogram(pref_dir_change(green_ind_concat))
title('HTP-')
subplot(1,2,2)
polarhistogram(pref_dir_change(red_ind_concat))
title('HTP+')
print(fullfile(fnout,'prefDirChange.pdf'),'-dpdf','-bestfit')

%% calculate norm_diff
norm_diff = nan(2,nCon,nSize,nKeep_total);
for i = 1:nKeep_total
    for iCon = 1:nCon
        for iSize = 1:nSize 
            %for stationary trials
            mean_pre_stat = mean(pref_allTrials_stat_concat{iCon,iSize,pre}{i});
            mean_post_stat=mean(pref_allTrials_stat_concat{iCon,iSize,post}{i});
            std_pre_stat = std(pref_allTrials_stat_concat{iCon,iSize,pre}{i});
            norm_diff_stat = (mean_post_stat-mean_pre_stat) / std_pre_stat;
    
            %for running trials
            mean_pre_loc = mean(pref_allTrials_loc_concat{iCon,iSize,pre}{i});
            mean_post_loc=mean(pref_allTrials_loc_concat{iCon,iSize,post}{i});
            std_pre_loc = std(pref_allTrials_loc_concat{iCon,iSize,pre}{i});
            norm_diff_loc = (mean_post_loc-mean_pre_loc)/ std_pre_loc;
    
            %putting data into matrix
            norm_diff(1,iCon,iSize,i)=norm_diff_stat; %first is stationary
            norm_diff(2,iCon,iSize,i)=norm_diff_loc; %second is running
%clear mean_pre_stat mean_post_stat std_pre_stat mean_pre_loc mean_post_loc std_pre_loc norn_diff_stat norm_diff_loc
        end 
    end
end
%remove any infiinty values resulting from divisions by zero, and turn
%those into NANs instead
norm_diff(find(norm_diff == -Inf))=NaN;
norm_diff(find(norm_diff == Inf))=NaN;

%% plot fraction suppressed and facilitated


norm_diff_red = norm_diff(:,:,:,red_ind_concat);
facil_red=norm_diff_red(:,:,:,:)>=1;
supp_red=norm_diff_red(:,:,:,:)<=-1;

N=length(red_ind_concat);



    facil_table_stat = squeeze(sum(facil_red(1,:,:,:),4)/N);
    supp_table_stat = squeeze(sum(supp_red(1,:,:,:),4)/N);
    
    figure;
    subplot(1,2,1)
    b=bar([1,2,3,4],[supp_table_stat(:,1),supp_table_stat(:,2)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
    b(1).FaceColor="#70D0F6"
    b(2).FaceColor="#0C8ABB"
    xticklabels({'12.5','25','50','100'})
    hold on
    title('Suppressed')
    ylim([0 .1])
    ylabel(["Fraction VIP cells"]) 
    xlabel(["Contrast"])
    set(gca,'TickDir','out')
    box off
    
    subplot(1,2,2)
    b=bar([1,2,3,4],[facil_table_stat(:,1),facil_table_stat(:,2)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
    b(1).FaceColor="#C983B1"
    b(2).FaceColor="#883367"
    xticklabels({'12.5','25','50','100'})
    hold on
    title('Facilitated')
    ylim([0 .1])
    %ylabel(["Fraction VIP cells"]) 
    xlabel(["Contrast"])
    set(gca,'TickDir','out')
    box off
    sgtitle('Stationary')
    
    x0=5;
    y0=5;
    width=3;
    height=1.75;
    set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Facil_supp_stat.pdf'),'-dpdf');
%% 
norm_diff_red = norm_diff(:,:,:,runningRed);
facil_red=norm_diff_red(:,:,:,:)>=1;
supp_red=norm_diff_red(:,:,:,:)<=-1;

N=length(red_ind_concat);

 facil_table_loc = squeeze(sum(facil_red(2,:,:,:),4)/N);
 supp_table_loc = squeeze(sum(supp_red(2,:,:,:),4)/N);
    
    figure;
    subplot(1,2,1)
    b=bar([1,2,3,4],[supp_table_loc(:,1),supp_table_loc(:,2)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
    b(1).FaceColor="#70D0F6"
    b(2).FaceColor="#0C8ABB"
    xticklabels({'12.5','25','50','100'})
    hold on
    title('Suppressed')
    ylim([0 .1])
    ylabel(["Fraction VIP cells"]) 
    xlabel(["Contrast"])
    set(gca,'TickDir','out')
    box off
    
    subplot(1,2,2)
    b=bar([1,2,3,4],[facil_table_loc(:,1),facil_table_loc(:,2)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
    b(1).FaceColor="#C983B1"
    b(2).FaceColor="#883367"
    xticklabels({'12.5','25','50','100'})
    hold on
    title('Facilitated')
    ylim([0 .1])
    %ylabel(["Fraction VIP cells"]) 
    xlabel(["Contrast"])
    set(gca,'TickDir','out')
    box off
    sgtitle('Running')
    x0=5;
    y0=5;
    width=3;
    height=1.75;
    set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Facil_supp_loc.pdf'),'-dpdf');

%%
%make a table of suppresses and facilitated cells for the cells that have
%both stationary and running within each condition
facil=norm_diff(:,:,:,:)>=1;
supp=norm_diff(:,:,:,:)<=-1;

supp_table_stat=nan(nCon,nSize);
facil_table_stat=nan(nCon,nSize);
supp_table_loc=nan(nCon,nSize);
facil_table_loc=nan(nCon,nSize);

for iCon = 1:nCon
    for iSize=1:nSize
        theseRedCells = intersect(red_ind_concat, find(runningByCondition(:,iCon,iSize)));
        supp_table_stat(iCon,iSize)=sum(supp(1,iCon,iSize,theseRedCells),4)/nRed(iCon,iSize);
        supp_table_loc(iCon,iSize)=sum(supp(2,iCon,iSize,theseRedCells),4)/nRed(iCon,iSize);

        facil_table_stat(iCon,iSize)=sum(facil(1,iCon,iSize,theseRedCells),4)/nRed(iCon,iSize);
        facil_table_loc(iCon,iSize)=sum(facil(2,iCon,iSize,theseRedCells),4)/nRed(iCon,iSize);


    end
end


figure;
subplot(2,2,1)
b=bar([1,2,3,4],[supp_table_stat(:,1),supp_table_loc(:,1)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
b(1).FaceColor="#70D0F6"
b(2).FaceColor="#0C8ABB"
xticklabels({'12.5','25','50','100'})
hold on
title('Suppressed 20-deg')
ylim([0 .6])
ylabel(["Fraction VIP cells"]) 
xlabel(["Contrast"])
set(gca,'TickDir','out')
box off

subplot(2,2,2)
b=bar([1,2,3,4],[facil_table_stat(:,1),facil_table_loc(:,1)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
b(1).FaceColor="#C983B1"
b(2).FaceColor="#883367"
xticklabels({'12.5','25','50','100'})
hold on
title('Facilitated 20-deg')
ylim([0 .6])
%ylabel(["Fraction VIP cells"]) 
xlabel(["Contrast"])
set(gca,'TickDir','out')
box off



subplot(2,2,3)
b=bar([1,2,3,4],[supp_table_stat(:,2),supp_table_loc(:,2)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
b(1).FaceColor="#70D0F6"
b(2).FaceColor="#0C8ABB"
xticklabels({'12.5','25','50','100'})
hold on
title('Suppressed fullfield')
ylim([0 .6])
ylabel(["Fraction VIP cells"]) 
xlabel(["Contrast"])
set(gca,'TickDir','out')
box off

subplot(2,2,4)
b=bar([1,2,3,4],[facil_table_stat(:,2),facil_table_loc(:,2)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
b(1).FaceColor="#C983B1"
b(2).FaceColor="#883367"
xticklabels({'12.5','25','50','100'})
hold on
title('Facilitated fullfield')
ylim([0 .6])
%ylabel(["Fraction VIP cells"]) 
xlabel(["Contrast"])
set(gca,'TickDir','out')
box off
x0=5;
y0=5;
width=4;
height=4;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'matched_facil_supp_byState.pdf'),'-dpdf');

%% norm diff for stationary vs running
for iSize = 1:nSize
    figure;
    subplot(1,2,1)
    boxchart(squeeze(norm_diff(1,:,iSize,red_ind_concat))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
    hold on
    scatter([1, 2, 3, 4],squeeze(norm_diff(1,:,iSize,red_ind_concat))',20,[.5 .15 .20], 'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)
    xticklabels({'12.5','25','50','100'})
    xlabel('Contrast(%)')
    ylabel('Normalized difference')
    ylim([-4 4])
    title('Stationary')
    hold off
    set(gca,'TickDir','out')
    box off
    x0=5;
    y0=5;
    width=1.25;
    height=2;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    
    
    subplot(1,2,2)
    boxchart(squeeze(norm_diff(2,:,iSize,red_ind_concat))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
    hold on
    scatter([1, 2, 3, 4],squeeze(norm_diff(2,:,iSize,red_ind_concat))',20,[.5 .15 .20], 'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)
    xticklabels({'12.5','25','50','100'})
    xlabel('Contrast(%)')
    %ylabel('Normalized difference')
    ylim([-10 10])
    title('Running')
    hold off
    set(gca,'TickDir','out')
    box off
    x0=5;
    y0=5;
    width=6;
    height=4;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    
    sgtitle(num2str(sizes(iSize)))

    print(fullfile(fnout,['normDiff_scatter_size', num2str(sizes(iSize)), '.pdf']),'-dpdf');
end

%% scatterplot and ttest for interneurons

%using cells matched across running and stationary for each stimulus
%condition

for iSize = 1:nSize
    for iCon = 1:nCon
        theseRedCells = intersect(red_ind_concat, find(runningByCondition(:,iCon,iSize)));
        figure; 
        subplot(1,2,1)
        scatter((conBySize_resp_stat_concat{pre}(theseRedCells,iCon,iSize)),(conBySize_resp_stat_concat{post}(theseRedCells,iCon,iSize)))
        set(gca, 'TickDir', 'out')
        axis square
        box off
        title(['Stationary con ', num2str(cons(iCon)), ' size ', num2str(sizes(iSize))])
        xlim([0 .2])
        ylim([0 .2])
        refline(1)
        [h,p]=ttest((conBySize_resp_stat_concat{pre}(theseRedCells,iCon,iSize)),(conBySize_resp_stat_concat{post}(theseRedCells,iCon,iSize)));
        if (p) < 0.05
            txt = {'p = ' num2str(round(p,3)),'n = ' num2str(length(theseRedCells))};
            text(.1,0.1,txt)
        end

        subplot(1,2,2)
        scatter((conBySize_resp_loc_concat{pre}(theseRedCells,iCon,iSize)),(conBySize_resp_loc_concat{post}(theseRedCells,iCon,iSize)))
        set(gca, 'TickDir', 'out')
        axis square
        box off
        title(['Running con ', num2str(cons(iCon)), ' size ', num2str(sizes(iSize))])
        xlim([0 .2])
        ylim([0 .2])
        refline(1)
        [h,p]=ttest((conBySize_resp_loc_concat{pre}(theseRedCells,iCon,iSize)),(conBySize_resp_loc_concat{post}(theseRedCells,iCon,iSize)));
        if (p) < 0.05
            txt = {'p = ' num2str(round(p,3)),'n = ' num2str(length(theseRedCells))};
            text(.15,0.15,txt)
        end

    end
end

%% plot large and small pupil stationary timecourses
% make figure with se shaded, one figure per contrast
%
% tc_green_avrg_small = cell(1,nd); %this will be the average across all green cells - a single line
% tc_red_avrg_small = cell(1,nd); %same for red
% tc_green_se_small = cell(1,nd); %this will be the se across all green cells
% tc_red_se_small = cell(1,nd); %same for red
% 
% 
% nGreen = nan(nCon,nSize);
% nRed = nan(nCon,nSize);
% 
% for id = 1:nd
%   for iCon = 1:nCon
%       for iSize = 1:nSize
%         tc_green_avrg_small{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_stat_smallPupil_concat{id}(:,green_ind_concat,iCon,iSize),2);
%         green_std=nanstd(tc_trial_avrg_stat_smallPupil_concat{id}(:,green_ind_concat,iCon,iSize),[],2);
%         tc_green_se_small{id}(:,iCon,iSize)=green_std/sqrt(length(green_ind_concat));
% 
%         tc_red_avrg_small{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_stat_smallPupil_concat{id}(:,red_ind_concat,iCon,iSize),2);
%         red_std=nanstd(tc_trial_avrg_stat_smallPupil_concat{id}(:,red_ind_concat,iCon,iSize),[],2);
%         tc_red_se_small{id}(:,iCon,iSize)=red_std/sqrt(length(red_ind_concat));
% 
% 
%         clear green_std red_std
%       end 
%     end
% end
% 
% z=double(nOn)/double(frame_rate);
% 
% %create a time axis in seconds
% t=1:(size(tc_green_avrg_small{1,1,1},1));
% t=(t-(double(stimStart)-1))/double(frame_rate);
% 
% for iCon = 1:nCon
%     for iSize = 1:nSize
%     figure
%     subplot(1,2,1) %for the first day
% 
% 
% 
%     ylim([-.05 .25]);
%     hold on
%     shadedErrorBar(t,tc_green_avrg_small{pre}(:,iCon,iSize),tc_green_se_small{pre}(:,iCon,iSize),'--k');
%     hold on
%     shadedErrorBar(t,tc_green_avrg_small{post}(:,iCon,iSize),tc_green_se_small{post}(:,iCon,iSize),'--b','transparent');
%     hold on
%     line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
%     hold on
%     line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
%     title(['-HTP',' n = ', num2str(nGreen(iCon,iSize))])
% 
%     ylabel('dF/F') 
%     xlabel('s') 
%     set(gca,'XColor', 'none','YColor','none')
% 
% 
%     subplot(1,2,2) %+HTP
%     shadedErrorBar(t,tc_red_avrg_small{pre}(:,iCon,iSize),tc_red_se_small{pre}(:,iCon,iSize),'k');
%     hold on
%     shadedErrorBar(t,tc_red_avrg_small{post}(:,iCon,iSize),tc_red_se_small{post}(:,iCon,iSize),'b');
%     ylim([-.05 .25]);
%     hold on
%     line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
%     hold on
%     line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
%     ylabel('dF/F') 
%     xlabel('s') 
%     title(['+HTP',' n = ', num2str(nRed(iCon,iSize))])
% 
%     x0=5;
%     y0=5;
%     width=4;
%     height=3;
%     set(gcf,'units','inches','position',[x0,y0,width,height])
%     set(gca,'XColor', 'none','YColor','none')
% 
%     sgtitle(['small pupil, con ' num2str(cons(iCon)) ' size ' num2str(sizes(iSize))])
% 
%     print(fullfile(fnout,[num2str(cons(iCon)) '_' num2str(sizes(iSize)) 'matched_small_cellType_timecourses.pdf']),'-dpdf');
%     end
% end 
% 
% 
% tc_green_avrg_large = cell(1,nd); %this will be the average across all green cells - a single line
% tc_red_avrg_large = cell(1,nd); %same for red
% tc_green_se_large = cell(1,nd); %this will be the se across all green cells
% tc_red_se_large = cell(1,nd); %same for red
% 
% 
% nGreen = nan(nCon,nSize);
% nRed = nan(nCon,nSize);
% 
% for id = 1:nd
%   for iCon = 1:nCon
%       for iSize = 1:nSize
%         tc_green_avrg_large{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_stat_largePupil_concat{id}(:,green_ind_concat,iCon,iSize),2);
%         green_std=nanstd(tc_trial_avrg_stat_largePupil_concat{id}(:,green_ind_concat,iCon,iSize),[],2);
%         tc_green_se_large{id}(:,iCon,iSize)=green_std/sqrt(length(green_ind_concat));
% 
%         tc_red_avrg_large{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_stat_largePupil_concat{id}(:,red_ind_concat,iCon,iSize),2);
%         red_std=nanstd(tc_trial_avrg_stat_largePupil_concat{id}(:,red_ind_concat,iCon,iSize),[],2);
%         tc_red_se_large{id}(:,iCon,iSize)=red_std/sqrt(length(red_ind_concat));
% 
% 
%         clear green_std red_std
%       end 
%     end
% end
% 
% z=double(nOn)/double(frame_rate);
% 
% %create a time axis in seconds
% t=1:(size(tc_green_avrg_large{1,1,1},1));
% t=(t-(double(stimStart)-1))/double(frame_rate);
% 
% for iCon = 1:nCon
%     for iSize = 1:nSize
%     figure
%     subplot(1,2,1) %for the first day
% 
% 
% 
%     ylim([-.05 .25]);
%     hold on
%     shadedErrorBar(t,tc_green_avrg_large{pre}(:,iCon,iSize),tc_green_se_large{pre}(:,iCon,iSize),'--k');
%     hold on
%     shadedErrorBar(t,tc_green_avrg_large{post}(:,iCon,iSize),tc_green_se_large{post}(:,iCon,iSize),'--b','transparent');
%     hold on
%     line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
%     hold on
%     line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
%     title(['-HTP',' n = ', num2str(nGreen(iCon,iSize))])
% 
%     ylabel('dF/F') 
%     xlabel('s') 
%     set(gca,'XColor', 'none','YColor','none')
% 
% 
%     subplot(1,2,2) %+HTP
%     shadedErrorBar(t,tc_red_avrg_large{pre}(:,iCon,iSize),tc_red_se_large{pre}(:,iCon,iSize),'k');
%     hold on
%     shadedErrorBar(t,tc_red_avrg_large{post}(:,iCon,iSize),tc_red_se_large{post}(:,iCon,iSize),'b');
%     ylim([-.05 .25]);
%     hold on
%     line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
%     hold on
%     line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
%     ylabel('dF/F') 
%     xlabel('s') 
%     title(['+HTP',' n = ', num2str(nRed(iCon,iSize))])
% 
%     x0=5;
%     y0=5;
%     width=4;
%     height=3;
%     set(gcf,'units','inches','position',[x0,y0,width,height])
%     set(gca,'XColor', 'none','YColor','none')
% 
%     sgtitle(['large pupil, con ' num2str(cons(iCon)) ' size ' num2str(sizes(iSize))])
% 
%     print(fullfile(fnout,[num2str(cons(iCon)) '_' num2str(sizes(iSize)) 'matched_large_cellType_timecourses.pdf']),'-dpdf');
%     end
% end 
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

% contrast 4 is the higher contrast
% size 2 is the large size

%% testing PCA
% 
% loc_pre = tc_trial_avrg_loc_concat{pre};
% loc_post = tc_trial_avrg_loc_concat{post};
% stat_post = tc_trial_avrg_stat_concat{post};
% stat_pre = tc_trial_avrg_stat_concat{pre};
% 
% loc_pre_red = loc_pre(:,red_ind_concat,:,:);
% loc_pre_green = loc_pre(:,green_ind_concat,:,:);
% stat_pre_green = stat_pre(:,green_ind_concat,:,:);
% stat_pre_red = stat_pre(:,red_ind_concat,:,:);
% loc_post_red = loc_post(:,red_ind_concat,:,:);
% stat_post_red = stat_post(:,red_ind_concat,:,:);
% stat_post_green = stat_post(:,green_ind_concat,:,:);
% loc_post_green = loc_post(:,green_ind_concat,:,:);
% 
% size(stat_pre_green)
% s1 = squeeze(mean(squeeze(stat_pre_green(:,:,:,1)),3));
% s2 = squeeze(mean(squeeze(stat_pre_green(:,:,:,2)),3));
% 
% s3 = squeeze(mean(squeeze(stat_post_green(:,:,:,1)),3));
% s4 = squeeze(mean(squeeze(stat_post_green(:,:,:,2)),3));
% 
% s_combined = [s1; s2];  % Concatenate observations
% [coeff, score_combined,~,~,explained_1] = pca(s_combined);
% 
% % Split scores
% score1 = score_combined(1:size(s1,1), :);
% score2 = score_combined(size(s1,1)+1:end, :);
% figure
% plot3(score1(:,1),score1(:,2),score1(:,3))
% xlabel('PC1')
% ylabel('PC2')
% zlabel('PC3')
% title('PCA Scores Plot stat Pre Green')
% hold on
% plot3(score2(:,1),score2(:,2),score2(:,3))
% legend('20deg','fullfield')
% hold off
% 
% s_combined = [s3; s4];  % Concatenate observations
% [coeff, score_combined,~,~,explained_2] = pca(s_combined);
% score1 = score_combined(1:size(s1,1), :);
% score2 = score_combined(size(s1,1)+1:end, :);
% figure
% plot3(score1(:,1),score1(:,2),score1(:,3))
% xlabel('PC1')
% ylabel('PC2')
% zlabel('PC3')
% title('PCA Scores Plot stat Post Green')
% hold on
% plot3(score2(:,1),score2(:,2),score2(:,3))
% legend('20dPost','FFPost')
% hold off
% 
% s_combined = [s1;s2;s3;s4];
% [coeff, score_combined,~,~,explained_3] = pca(s_combined);
% score1 = score_combined(1:90,:);
% score2 = score_combined(91:180,:);
% score3 = score_combined(181:270,:);
% score4 = score_combined(271:360,:);
% figure;
% plot3(score1(:,1),score1(:,2),score1(:,3))
% xlabel('PC1')
% ylabel('PC2')
% zlabel('PC3')
% title('PCA Scores Plot stat Both Green')
% hold on
% plot3(score2(:,1),score1(:,2),score1(:,3))
% plot3(score3(:,1),score1(:,2),score1(:,3))
% plot3(score4(:,1),score1(:,2),score1(:,3))
% legend('20dPre','FFPre','20dPost','FFPost')
% hold off
% % same thing for red
% size(stat_pre_red)
% s1 = squeeze(mean(squeeze(stat_pre_red(:,:,:,1)),3));
% s2 = squeeze(mean(squeeze(stat_pre_red(:,:,:,2)),3));
% 
% s3 = squeeze(mean(squeeze(stat_post_red(:,:,:,1)),3));
% s4 = squeeze(mean(squeeze(stat_post_red(:,:,:,2)),3));
% 
% s_combined = [s1; s2];  % Concatenate observations
% [coeff, score_combined,~,~,explained_1] = pca(s_combined);
% 
% % Split scores
% score1 = score_combined(1:size(s1,1), :);
% score2 = score_combined(size(s1,1)+1:end, :);
% figure
% plot3(score1(:,1),score1(:,2),score1(:,3))
% xlabel('PC1')
% ylabel('PC2')
% zlabel('PC3')
% title('PCA Scores Plot stat Pre red')
% hold on
% plot3(score2(:,1),score2(:,2),score2(:,3))
% legend('20deg','fullfield')
% hold off
% 
% s_combined = [s3; s4];  % Concatenate observations
% [coeff, score_combined,~,~,explained_2] = pca(s_combined);
% score1 = score_combined(1:size(s1,1), :);
% score2 = score_combined(size(s1,1)+1:end, :);
% figure
% plot3(score1(:,1),score1(:,2),score1(:,3))
% xlabel('PC1')
% ylabel('PC2')
% zlabel('PC3')
% title('PCA Scores Plot stat Post red')
% hold on
% plot3(score2(:,1),score2(:,2),score2(:,3))
% legend('20dPost','FFPost')
% hold off
% 
% s_combined = [s1;s2;s3;s4];
% [coeff, score_combined,~,~,explained_3] = pca(s_combined);
% score1 = score_combined(1:90,:);
% score2 = score_combined(91:180,:);
% score3 = score_combined(181:270,:);
% score4 = score_combined(271:360,:);
% figure;
% plot3(score1(:,1),score1(:,2),score1(:,3))
% xlabel('PC1')
% ylabel('PC2')
% zlabel('PC3')
% title('PCA Scores Plot stat Both red')
% hold on
% plot3(score2(:,1),score1(:,2),score1(:,3))
% plot3(score3(:,1),score1(:,2),score1(:,3))
% plot3(score4(:,1),score1(:,2),score1(:,3))
% legend('20dPre','FFPre','20dPost','FFPost')
% hold off
% 
% 
% 
% %% for loc
% size(loc_pre_green)
% s1 = squeeze(mean(squeeze(loc_pre_green(:,:,:,1)),3));
% s2 = squeeze(mean(squeeze(loc_pre_green(:,:,:,2)),3));
% 
% s3 = squeeze(mean(squeeze(loc_post_green(:,:,:,1)),3));
% s4 = squeeze(mean(squeeze(loc_post_green(:,:,:,2)),3));
% 
% s_combined = [s1; s2];  % Concatenate observations
% [coeff, score_combined,~,~,explained_1] = pca(s_combined);
% 
% % Split scores
% score1 = score_combined(1:size(s1,1), :);
% score2 = score_combined(size(s1,1)+1:end, :);
% figure
% plot3(score1(:,1),score1(:,2),score1(:,3))
% xlabel('PC1')
% ylabel('PC2')
% zlabel('PC3')
% title('PCA Scores Plot loc Pre Green')
% hold on
% plot3(score2(:,1),score2(:,2),score2(:,3))
% legend('20deg','fullfield')
% hold off
% 
% s_combined = [s3; s4];  % Concatenate observations
% [coeff, score_combined,~,~,explained_2] = pca(s_combined);
% score1 = score_combined(1:size(s1,1), :);
% score2 = score_combined(size(s1,1)+1:end, :);
% figure
% plot3(score1(:,1),score1(:,2),score1(:,3))
% xlabel('PC1')
% ylabel('PC2')
% zlabel('PC3')
% title('PCA Scores Plot loc Post Green')
% hold on
% plot3(score2(:,1),score2(:,2),score2(:,3))
% legend('20dPost','FFPost')
% hold off
% 
% s_combined = [s1;s2;s3;s4];
% [coeff, score_combined,~,~,explained_3] = pca(s_combined);
% score1 = score_combined(1:90,:);
% score2 = score_combined(91:180,:);
% score3 = score_combined(181:270,:);
% score4 = score_combined(271:360,:);
% figure;
% plot3(score1(:,1),score1(:,2),score1(:,3))
% xlabel('PC1')
% ylabel('PC2')
% zlabel('PC3')
% title('PCA Scores Plot loc Both Green')
% hold on
% plot3(score2(:,1),score1(:,2),score1(:,3))
% plot3(score3(:,1),score1(:,2),score1(:,3))
% plot3(score4(:,1),score1(:,2),score1(:,3))
% legend('20dPre','FFPre','20dPost','FFPost')
% hold off
% % same thing for red
% size(loc_pre_red)
% s1 = squeeze(mean(squeeze(loc_pre_red(:,:,:,1)),3));
% s2 = squeeze(mean(squeeze(loc_pre_red(:,:,:,2)),3));
% 
% s3 = squeeze(mean(squeeze(loc_post_red(:,:,:,1)),3));
% s4 = squeeze(mean(squeeze(loc_post_red(:,:,:,2)),3));
% 
% s_combined = [s1; s2];  % Concatenate observations
% [coeff, score_combined,~,~,explained_1] = pca(s_combined);
% 
% % Split scores
% score1 = score_combined(1:size(s1,1), :);
% score2 = score_combined(size(s1,1)+1:end, :);
% figure
% plot3(score1(:,1),score1(:,2),score1(:,3))
% xlabel('PC1')
% ylabel('PC2')
% zlabel('PC3')
% title('PCA Scores Plot loc Pre red')
% hold on
% plot3(score2(:,1),score2(:,2),score2(:,3))
% legend('20deg','fullfield')
% hold off
% 
% s_combined = [s3; s4];  % Concatenate observations
% [coeff, score_combined,~,~,explained_2] = pca(s_combined);
% score1 = score_combined(1:size(s1,1), :);
% score2 = score_combined(size(s1,1)+1:end, :);
% figure
% plot3(score1(:,1),score1(:,2),score1(:,3))
% xlabel('PC1')
% ylabel('PC2')
% zlabel('PC3')
% title('PCA Scores Plot loc Post red')
% hold on
% plot3(score2(:,1),score2(:,2),score2(:,3))
% legend('20dPost','FFPost')
% hold off
% 
% s_combined = [s1;s2;s3;s4];
% [coeff, score_combined,~,~,explained_3] = pca(s_combined);
% score1 = score_combined(1:90,:);
% score2 = score_combined(91:180,:);
% score3 = score_combined(181:270,:);
% score4 = score_combined(271:360,:);
% figure;
% plot3(score1(:,1),score1(:,2),score1(:,3))
% xlabel('PC1')
% ylabel('PC2')
% zlabel('PC3')
% title('PCA Scores Plot loc Both red')
% hold on
% plot3(score2(:,1),score1(:,2),score1(:,3))
% plot3(score3(:,1),score1(:,2),score1(:,3))
% plot3(score4(:,1),score1(:,2),score1(:,3))
% legend('20dPre','FFPre','20dPost','FFPost')
% hold off
% 
% 
