function [] = extraction_Eye(day_id,ref)

ds = 'DART_V1_YM90K_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc =  behavConstsDART; %directories
eval(ds);

%day_id = 169; %enter post-DART day
pre_day = expt(day_id).multiday_matchdays;

nd=2; %hardcoding the number of days for now
experimentFolder = 'VIP_YM90K';

%trialsToSkip
skipAction = 1; % 1- deletes trials % 2- replaces trials with NaN
skip{1} = []; % day 1 should be the reference session (whatever was identified as pre_day)
skip{2} = [];

%pupil?
if day_id == 18
    doEye = 0;
else
    doEye = 1;
end

mouse = expt(day_id).mouse;

fnout = fullfile(rc.achAnalysis,mouse);
if expt(day_id).multiday_timesincedrug_hours>0
    dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end
prompt = 'Which sesson was used as reference for matching: 0- baseline, 1- post-DART';
            x = ref;
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


fn_multi = fullfile(rc.achAnalysis,experimentFolder,mouse,['multiday_' dart_str]);

cd(fn_multi)
load(fullfile(fn_multi,'timecourses.mat'))
%load(fullfile(fn_multi,'multiday_alignment.mat'))
load(fullfile(fn_multi,'input.mat'))

[input.modifiedAction] = deal([]);
[input.modifiedFields] = deal([]);
[input.modifiedStruct] = deal([]);
[input.modifiedTime] = deal([]);
[input.modifiedTrials] = deal([]);

if ~isempty(skip{1})
    input(1) = trialDropper(input(1),skip{1},skipAction);
end

if ~isempty(skip{2})
    input(2) = trialDropper(input(2),skip{2},skipAction);
end

frame_rate = input.frameImagingRateMs;
% %% finding red fluorescence level
allDays = [day_id,pre_day];


for id = 1 %currently only doing this for the reference day, regardless of 
    % which direction the data was matched
    mouse = expt(allDays(id)).mouse;
    date = expt(allDays(id)).date;
    imgFolder = expt(allDays(id)).contrastxori_runs{1};
    fn = fullfile(rc.achAnalysis,experimentFolder, mouse,date,imgFolder);
    cd(fn);
    load(fullfile(fn,'redImage.mat'));
    load(fullfile(fn,'mask_cell.mat'));
    
    %use stackGetTimeCourses to extract the red fluorescence within each mask
    red_fluor_mask = stackGetTimeCourses(redChImg, mask_cell);
    nCells=max(max(mask_cell));
        for i = 1:nCells
            red_fluor_np(i) = stackGetTimeCourses(redChImg, mask_np(:,:,i));
        end
    %using the reference day
    %red_fluor_match=red_fluor_all(:,match_ind); %if we want to use the
    %np-subtracted red
    red_fluor_match=red_fluor_mask(:,match_ind); %to find the red within a cell, NOT no-subtracted
    z_red_fluor=zscore(red_fluor_match);
    load(fullfile(fn_multi,'multiday_alignment.mat'))
    clear red_fluor_all red_fluor_mask red_fluor_np
    % get green fluor level
    %using the reference day
    green_fluor_match=mean(cellTCs_match{1},1);   

end

if doEye == 1
    pupil=cell(1,nd);
    for id=1:nd
        pupil{id}=load(fullfile(fn,'pupil.mat'));
    end
end
%% stimulus props

nOn = input(1).nScansOn;
nOff = input(1).nScansOff;

%tells the contrast, direction and orientation for each trial each day
tCon_match = cell(1,nd);
tDir_match = cell(1,nd);
tOri_match = cell(1,nd);
tSize_match = cell(1,nd);

%find the contrasts, directions and orientations for each day
%in case of instances where the number of trails actually collected was not
%consistent with the number mWorks thinks occured, only take 1:nTrials as
%dictated above based on the number of frames recorded
for id = 1:nd
    nTrials(id) = length(input(id).tGratingContrast); %to account for times 
%when there is a disruption before the full set of trials is collcted, I'm 
%determining the number of trials each day by how many frames of data I 
%have divided by the number of frames per trial
    tCon_match{id} = celleqel2mat_padded(input(id).tGratingContrast(1:nTrials(id)));
    tDir_match{id} = celleqel2mat_padded(input(id).tGratingDirectionDeg(1:nTrials(id)));
    tOri_match{id} = tDir_match{id};
    tOri_match{id}(find(tDir_match{id}>=180)) = tDir_match{id}(find(tDir_match{id}>=180))-180;
    tSize_match{id} = celleqel2mat_padded(input(id).tGratingDiameterDeg(1:nTrials(id)));
end
oris = unique(tOri_match{1});
dirs = unique(tDir_match{1});
cons = unique(tCon_match{1});
sizes = unique(tSize_match{1});
nOri = length(oris);
nCon = length(cons);
nDir = length(dirs);
nSize = length(sizes);

%% convert to trials

data_dfof_trial_match = cell(1,nd); %make an empty array that is 1 by however many days there are (1X2 usually)

fractTimeActive_match = cell(1,nd);
cellstd_match = cell(1,nd);
for id = 1:nd %cycle through days
    mouse = expt(allDays(id)).mouse;
    date = expt(allDays(id)).date;
    imgFolder = expt(allDays(id)).contrastxori_runs{1};
    imgMatFile = [imgFolder '_000_000.mat']
    dataPath = fullfile(rc.achData, mouse, date, imgFolder);
    load(fullfile(dataPath,imgMatFile))
    [cStimOn{id} stimOffs{id}] = photoFrameFinder_Sanworks(info.frame);
    cStimOn{id}(skip{id}) = [];
    [nFrames nCells] = size(cellTCs_match{id});
    
    data_trial_match = nan(nOn+nOff,nTrials(id),nCells);
    
    for itrial = 1:nTrials(id)
      if ~isnan(cStimOn{id}(itrial)) & (cStimOn{id}(itrial)+nOn+nOff/2)<nFrames
        data_trial_match(:,itrial,:) = cellTCs_match{id}(cStimOn{id}(itrial)-nOff/2:cStimOn{id}(itrial)-1+nOn+nOff/2,:);
      end
    end

    fractTimeActive_match{id} = zeros(1,nCells);

    data_f_match = mean(data_trial_match(1: (nOff/2),:,:),1);
    data_dfof_trial_match{id} = bsxfun(@rdivide,bsxfun(@minus,data_trial_match,data_f_match),data_f_match);
    meansub_match = cellTCs_match{id}-nanmean(cellTCs_match{id},1);
    cellstd = nanstd(meansub_match,[],1);
    cellstd_match{id}=cellstd;
    for iCell = 1:nCells
        fractTimeActive_match{id}(:,iCell) = length(find(meansub_match(:,iCell)>3.*cellstd(1,iCell)))./nFrames;
    end
    
end
clear  data_f_match cellstd 

% MUST USE NANMEAN INSTEAD OF MEAN MOVING FORWARD SINCE I SET THE PADDING
% VALUES TO NAN

%% looking at wheel speed
wheel_speed = cell(1,nd); %output from wheelSpeedCalc, has whlspd for each frame

for id = 1:nd
    wheel_speed{id} = wheelSpeedCalc(input(id),32,expt(allDays(1)).wheelColor); 
    nanmean(wheel_speed{id})
end
wheel_speed_clean = cell(1,nd); %wheel_speed but jitter values are set to 0
for id = 1:nd
    wheel_speed_clean{id}=wheel_speed{id};
    wheel_speed_clean{id}(abs(wheel_speed_clean{id})<5.37)=0;
end

wheel_tc = cell(1,nd); % trial tc from whlspd_clean. RIx is calculated from this.
wheel_trial_avg= cell(1,nd);
wheel_tc_raw = cell(1,nd); % trial tc from wheel_speed but calculating with absolute values
wheel_trial_avg_raw= cell(1,nd);
RIx = cell(1,nd);

for id = 1:nd
    wheel_tc{id}=nan(nOn+nOff, nTrials(id));
    wheel_tc_raw{id}=nan(nOn+nOff, nTrials(id));
    for iTrial = 1:nTrials(id)
        if ~isnan(cStimOn{id}(iTrial)) & (cStimOn{id}(iTrial)+nOn+nOff/2)<nFrames
            wheel_tc{id}(:,iTrial) = wheel_speed_clean{id}(cStimOn{id}(iTrial)-nOff/2:cStimOn{id}(iTrial)-1+nOn+nOff/2);
            wheel_tc_raw{id}(:,iTrial) = abs(wheel_speed{id}(cStimOn{id}(iTrial)-nOff/2:cStimOn{id}(iTrial)-1+nOn+nOff/2));
        end
    end
    wheel_trial_avg{id} = mean(wheel_tc{id}(nOff/2:nOn+nOff/2,:),1);
    wheel_trial_avg_raw{id} = mean(wheel_tc_raw{id}(nOff/2:nOn+nOff/2,:),1);
    RIx{id} = wheel_trial_avg{id}>2; %~5 is the noise level in the wheel movement
    sum(RIx{id})
end


%% how stationary is stationary?
stat_wheelspd = cell(1,nd);
stat_wheelspd_avg = nan(1,nd);
stat_wheelspd_avg_noZero = nan(1,nd);

for id = 1:nd
    temp_whlspd = mean(wheel_tc_raw{id}(nOff/2:nOn+nOff/2,~RIx{id}),1);
    temp_whlspd(isnan(temp_whlspd)) = 0;
    stat_wheelspd{id} = temp_whlspd;
    stat_wheelspd_avg(id) = mean(temp_whlspd);
    stat_wheelspd_avg_noZero(id) = mean(temp_whlspd(temp_whlspd~=0));
end

% figure
% title([mouse ' stat trial wheel speed'])
% hold on
% histogram(stat_wheelspd{pre},0:0.1:2)
% histogram(stat_wheelspd{post},0:0.1:2)
% hold off
% xlim([0 2])
% ylim([0 100])
% pre_legend = ['pre n=' num2str(sum(~RIx{pre})) '; ' num2str(stat_wheelspd_avg(pre)) ';' num2str(stat_wheelspd_avg_noZero(pre))];
% post_legend = ['post n=' num2str(sum(~RIx{post})) '; ' num2str(stat_wheelspd_avg(post)) ';' num2str(stat_wheelspd_avg_noZero(post))];
% legend(pre_legend,post_legend)
% print(fullfile(fn_multi,[mouse '_stat_wheelSpd.pdf']),'-dpdf')
% print(fullfile('G:\home\jerry\reports\poster\2025\plots',[mouse '_stat_wheelSpd.pdf']),'-dpdf');
% save(fullfile(fn_multi,'wheelSpds.mat'),'stat_wheelspd_avg','stat_wheelspd_avg_noZero','stat_wheelspd');
%% extract running onsets
data_dfof_runOnset_match = cell(1,nd);
nRunOnsets=[];
for id = 1:nd
    fwdWheelClean = wheel_speed_clean{id};
    fwdWheelClean(fwdWheelClean<0)=0;
    % moveMean = movmean(fwdWheelClean,3);
    % smoothRun=smooth(fwdWheelClean,3);
    % figure;plot(fwdWheelClean(900:1300));hold on;plot(smoothRun(900:1300));hold on;plot(moveMean(900:1300));
    
    % find and refine the onset indices
    stillFrames = logical(fwdWheelClean<2); %find the stationary periods based on the sliding window
    runFrames = logical(fwdWheelClean>2); %find the stationary periods based on the sliding window

    runOnsets=(find(diff(stillFrames) == -1)+1); %find the transitions from stationary to running
    runOffsets=(find(diff(stillFrames) == 1)+1); %find the transitions from running to stationary
  
    %subset to onsets that are preceeded by 1s stillness and followed by at
    %least 1s of running
    runOnsets_clean = runOnsets;
    toDrop = [];
    for iOnset = 1:length(runOnsets_clean)
        if runOnsets_clean(iOnset) < frame_rate+1
            toDrop=[toDrop,iOnset]; %drop onsets that are less than 1 second from the start of the session, because in that case I can't test of my next criterion
        elseif runOnsets_clean(iOnset) > frame_rate+1
        testWinStill=[runOnsets_clean(iOnset)-frame_rate:runOnsets_clean(iOnset)];
        testWinRun=[runOnsets_clean(iOnset):runOnsets_clean(iOnset)+frame_rate];
            if sum(stillFrames(testWinStill))<(frame_rate*.75) |sum(runFrames(testWinRun))<(frame_rate*.75)
                toDrop=[toDrop,iOnset]; %gather a list of onsets that don't meet my criteria
            end 
        end
    end
    runOnsets_clean(toDrop)=[];
    %
    %make a vector that is 1 during the ITI and 0 when the stim is on
    ITI = ones(1,nFrames);
    for iTrial = 1:nTrials(id)
        ITI(cStimOn{id}(iTrial):cStimOn{id}(iTrial)+nOn)=0;
    end
    
    %subset the running onsets to only include those during the ITI
    ITIOnsets = runOnsets_clean;
    ITIOnsets(ITI(runOnsets_clean)==0)=[];
    
    % get neural data for run onsets 
    % extract timecourses during a window around the run onset
    onsetWin = 1; %how many seconds you want the window the be, on EITHER side of the onset
    
    data_dfof_runOnset = nan(onsetWin*2*frame_rate,nCells,length(ITIOnsets));
    runConfirmation = nan(onsetWin*2*frame_rate,length(ITIOnsets));
    for iOnset = 1:length(ITIOnsets)
        %find inds for a window around the onset
        fullWindow = [(ITIOnsets(iOnset)-(onsetWin*frame_rate)):(ITIOnsets(iOnset)+(onsetWin*frame_rate)-1)];
        bslnWindow = [(ITIOnsets(iOnset)-(frame_rate/2)):ITIOnsets(iOnset)];%the last half second before the run onset. Will use this as the baseline for df/f
        tempFull=cellTCs_match{id}(fullWindow,:);
        tempBsln=mean(cellTCs_match{id}(bslnWindow,:));
        %this is F, convert to dfof. 
        tempDfof = bsxfun(@rdivide,bsxfun(@minus,tempFull,tempBsln),tempBsln);
        data_dfof_runOnset(:,:,iOnset)=tempDfof;
        runConfirmation(:,iOnset) = fwdWheelClean(fullWindow);
    end
    
    % figure; plot(nanmean(runConfirmation,2)) %to check whether running actually does increase around the time of these onsets.
    
    data_dfof_runOnset_match{id} = mean(data_dfof_runOnset,3,'omitmissing'); %frames x cells (for all matched cells) averaged over all the onsets
    nRunOnsets=[nRunOnsets length(ITIOnsets)];
end

nRunOnsets


%% get large/small pupil trials
if doEye == 1
    for id = 1:nd
        pupil{id}.rad.stim(skip{id}) = [];
    end
    statPupilBothDays =horzcat(pupil{pre}.rad.stim(~RIx{pre}),pupil{post}.rad.stim(~RIx{post})); %combine all the pupil values for the two days
    statPupilThreshold=prctile(statPupilBothDays,50);
    %plot(statPupilBothDays); hline(statPupilThreshold); xlabel('trials, both days');ylabel('pupil diam'); hold off
    %print(fullfile(fn_multi,'pupilTraceWThreshold.pdf'),'-dpdf');
    
    pupilMeans = nan(nd,5);
    % for each day, the first column is the mean pupil size for stat trials
    % above threshold, the second column is the mean pupil size for stat trials
    % below threshold, and the third column is the mean pupil size for all
    % running trials
    PIx_stat = cell(4,nd); %pupil index for each day, first cell is inds for 
    % stationary large pupil, second cell is inds for stationary small pupil
    motorByPupil = nan(nd,2);
    statPupil = cell(3,nd);
    % 1: above thresh
    % 2: below thresh
    % 3: every stationary
    for id = 1:nd
        PIx_temp=pupil{id}.rad.stim > statPupilThreshold;
        PIx_stat{1,id}= logical(PIx_temp.*~RIx{id}); % stat, above threshold
        PIx_stat{2,id}= logical(~PIx_temp.*~RIx{id}); % stat, below threshold
        PIx_stat{3,id}= logical(PIx_temp.*RIx{id}); % running, above thresh
        PIx_stat{4,id}= logical(~PIx_temp.*RIx{id}); % running, below thresh
        pupilMeans(id,1)=mean(pupil{id}.rad.stim(PIx_stat{1,id}), 'omitmissing'); % >50%, stat
        pupilMeans(id,2)=mean(pupil{id}.rad.stim(PIx_stat{2,id}), 'omitmissing'); % <50%, stat
        pupilMeans(id,3)=mean(pupil{id}.rad.stim(RIx{id}), 'omitmissing'); %mean of all running trials
        pupilMeans(id,4)=mean(pupil{id}.rad.stim(PIx_stat{3,id}), 'omitmissing'); % >50%, loc
        pupilMeans(id,5)=mean(pupil{id}.rad.stim(PIx_stat{4,id}), 'omitmissing'); % <50%, loc
        pupilMeans(isnan(pupilMeans)) = 0;
        motorByPupil(id,1)=mean(wheel_trial_avg_raw{id}(PIx_stat{1,id}),'omitmissing'); % mean whlspd of large pupil trials
        motorByPupil(id,2)=mean(wheel_trial_avg_raw{id}(PIx_stat{2,id}),'omitmissing'); % mean whlspd of small pupil trials
        statPupil{1,id} = pupil{id}.rad.stim(PIx_stat{1,id});
        statPupil{2,id} = pupil{id}.rad.stim(PIx_stat{2,id});
        statPupil{3,id} = pupil{id}.rad.stim(~RIx{id});
    end
    save(fullfile(fn_multi,'pupilMeans.mat'),'pupilMeans','motorByPupil','statPupil');
end

end