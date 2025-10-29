clear all; clear global; 
close all
clc

%% Dataset and experiment setup
prompt = 'Enter name of instructions file: ';
instr = input(prompt, 's');
clear prompt
eval(instr);

ds = instructions.ds;
eval(ds);

rc = behavConstsDART;

day_id = str2double(instructions.session);

if length(expt) < day_id
    error('Day_id %d not valid for this dataset', day_id);
else
    pre_day = expt(day_id).multiday_matchdays;
    fprintf('Analyzing sessions: %s\n', num2str([day_id, pre_day]));
end
allDays = [day_id, pre_day];
nd = 2;
mouse = expt(day_id).mouse;
experimentFolder = expt(day_id).exptType;

if expt(day_id).multiday_timesincedrug_hours > 0
    dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end

x = instructions.refDay;
switch x
    case '1'
        pre = 1;
        post = 2;
        fprintf('Baseline used as reference\n');
    case '2'
        pre = 2;
        post = 1;
        fprintf('Post-DART used as reference\n');
end
clear x instr

stillLength = 5;
runLength = 2;

%% Load essential data
fn_multi = fullfile(rc.analysis, experimentFolder, mouse, ['multiday_' dart_str]);
cd(fn_multi)

load(fullfile(fn_multi, 'timecourses.mat'), 'cellTCs_match','red_ind_match');
load(fullfile(fn_multi, 'input.mat'));

OG_inputStructure = input;
clear input
%make this loop through the days
if exist(fullfile(fn_multi, 'correctedInputStructure.mat'), 'file')
    load("correctedInputStructure.mat");
    cStimOn = cell2mat(correctedInputStructure{id}.cStimOn);
    fprintf('Using existing corrected input structure\n');
else
    %make a corrected input structure for both days
    for id = 1:nd
        %load the data for that day
        mouse_temp = expt(allDays(id)).mouse;
        date = expt(allDays(id)).date;
        imgFolder = expt(allDays(id)).contrastxori_runs{1};
        imgMatFile = [imgFolder '_000_000.mat'];
        dataPath = fullfile(rc.data, mouse_temp, date, imgFolder);
        load(fullfile(dataPath,imgMatFile));
        %load data for that day
        if isfield(info, 'frame')
            OG_input_temp = OG_inputStructure(id);
            %check whether there is a input.cStimOn field - if not, create
            %on using nOn and nOff and add this to OG_input_temp
   
            [cStimOn, stimOffs] = photoFrameFinder_Sanworks(info.frame);
            correctedInputStructure_temp = counterValCorrect(OG_input_temp,cStimOn);
            cStimOn = cell2mat(correctedInputStructure_temp.cStimOn);
            fprintf('Creating corrected input structure with photodiode\n');
        else
            correctedInputStructure_temp = counterValCorrect_noPhotodiode(OG_inputStructure(id));
            cStimOn = cell2mat(correctedInputStructure_temp.cStimOn);
            fprintf('Creating corrected input structure without photodiode\n');
        end
    end
end

frame_rate = inputStructure(1).frameImagingRateMs;

%% Load wheel speed data for all days

wheel_speed = cell(1, nd);
%change this to use the corrected inputStrucutre
for id = 1:nd
    wheel_speed{id} = wheelSpeedCalc(inputStructure(id), 32, expt(allDays(1)).wheelColor);
end

wheel_speed_clean = cell(1, nd);
for id = 1:nd
    wheel_speed_clean{id} = wheel_speed{id};
    wheel_speed_clean{id}(abs(wheel_speed_clean{id}) < 5.37) = 0;
end

nOn = inputStructure(1).nScansOn;
nOff = inputStructure(1).nScansOff;

nTrials = zeros(1, nd);
for id = 1:nd
    nTrials(id) = size(cellTCs_match{id}, 1) / (nOn + nOff);
end

fprintf('\n=== LOADED DATA SUMMARY ===\n');
fprintf('Mouse: %s\n', mouse);
fprintf('Days: %s\n', mat2str(allDays));
fprintf('Frame rate: %.1f Hz\n', frame_rate);
fprintf('Trials per day: %s\n', mat2str(nTrials));
fprintf('Cells matched: %d\n', size(cellTCs_match{1}, 2));
fprintf('Red cells matched: %d\n', sum(red_ind_match));
if pre == 1
    ref_session = 'baseline';
else
    ref_session = 'post-DART';
end
fprintf('Reference session: %s \n', ref_session);
fprintf('ITI filtering: disabled (but ITI status will be documented)\n');
frame_rate_double = double(frame_rate);

%% looking at wheel speed
wheel_tc = cell(1,nd);
wheel_trial_avg= cell(1,nd);
wheel_tc_raw = cell(1,nd);
wheel_trial_avg_raw= cell(1,nd);
RIx = cell(1,nd);

for id = 1:nd
    wheel_tc{id}=zeros(nOn+nOff, nTrials(id));
    wheel_tc_raw{id}=zeros(nOn+nOff, nTrials(id));
    for iTrial = 1:nTrials(id)
        wheel_tc{id}(:,iTrial) = wheel_speed_clean{id}(1+((iTrial-1).*(nOn+nOff)):iTrial.*(nOn+nOff));
        wheel_tc_raw{id}(:,iTrial) = abs(wheel_speed{id}(1+((iTrial-1).*(nOn+nOff)):iTrial.*(nOn+nOff)));
    end
    wheel_trial_avg{id} = mean(wheel_tc{id}(nOff:nOn+nOff,:),1);
    wheel_trial_avg_raw{id} = mean(wheel_tc_raw{id}(nOff:nOn+nOff,:),1);
    RIx{id} = wheel_trial_avg{id}>2;
    mean(RIx{id})
end

%% Extract running onsets
speed_threshold = 5.37;
onsetWin = 5;

nRunOnsets = [];
meanSpeedAfterOnset = [];
onsetITIStatus = cell(1, nd);

for id = 1:nd
    mouse_temp = expt(allDays(id)).mouse;
    date = expt(allDays(id)).date;
    imgFolder = expt(allDays(id)).contrastxori_runs{1};
    imgMatFile = [imgFolder '_000_000.mat'];
    dataPath = fullfile(rc.data, mouse_temp, date, imgFolder);
    load(fullfile(dataPath,imgMatFile));
    %change this to use the cStimOn from the corrected input structure,
    %which will be defined above
    % if isfield(info, 'frame')
    %     [cStimOn, stimOffs] = photoFrameFinder_Sanworks(info.frame);
    %     fprintf('Getting cStimOn from photo diode\n');
    % elseif exist(fullfile(fn_multi, 'correctedInputStructure.mat'), 'file')
    %     load("correctedInputStructure.mat");
    %     cStimOn = cell2mat(correctedInputStructure{id}.cStimOn);
    %     fprintf('Getting cStimOn from existing corrected input structure\n');
    % else
    %     correctedInputStructure_temp = counterValCorrect_noPhotodiode(inputStructure(id));
    %     cStimOn = cell2mat(correctedInputStructure_temp.cStimOn);
    %     fprintf('Creating corrected input structure to get cStimOn\n');
    % end
    [nFrames, nCells] = size(cellTCs_match{id});
    
    fwdWheelClean = wheel_speed_clean{id};
    fwdWheelClean(fwdWheelClean<0) = 0;
    
    stillFrames = logical(fwdWheelClean < speed_threshold);
    runFrames = logical(fwdWheelClean >= speed_threshold);
    
    runOnsets = find(diff(stillFrames) == -1) + 1;
    runOffsets = find(diff(stillFrames) == 1) + 1;
    
    runOnsets_clean = runOnsets;
    toDrop = [];
    
    for iOnset = 1:length(runOnsets_clean)
        onsetFrame = runOnsets_clean(iOnset);
        
        if onsetFrame < frame_rate + 1 || ...
           onsetFrame <= (frame_rate * stillLength) || ...
           onsetFrame + (frame_rate * runLength) > length(stillFrames)
            toDrop = [toDrop, iOnset];
            continue;
        end
        
        testWinStill = (onsetFrame - stillLength*frame_rate):(onsetFrame - 1);
        testWinRun = (onsetFrame + 1):(onsetFrame + runLength*frame_rate);
        
        stillnessReq = sum(stillFrames(testWinStill)) == length(testWinStill);
        runningReq = sum(runFrames(testWinRun)) >= (runLength * frame_rate * 0.5);
        
        if ~stillnessReq || ~runningReq
            toDrop = [toDrop, iOnset];
        end
    end
    
    runOnsets_clean(toDrop) = [];
    
    if isempty(runOnsets_clean)
        disp(['No valid running onsets found for day ' num2str(id)]);
        data_dfof_runOnset_match{id} = nan(onsetWin*2*frame_rate, nCells);
        mean_resp_runOnset_match{id} = nan(1, nCells);
        nRunOnsets = [nRunOnsets 0];
        meanSpeedAfterOnset{id} = [];
        onsetITIStatus{id} = [];
        continue;
    end
    
    ITI = ones(1, nFrames);
    for iTrial = 1:inputStructure(id).trialSinceReset
        if cStimOn(iTrial) > 0 && cStimOn(iTrial) + inputStructure(id).nScansOn <= nFrames
            ITI(cStimOn(iTrial):cStimOn(iTrial) + inputStructure(id).nScansOn) = 0;
        end
    end
    
    finalOnsets = runOnsets_clean;
    
    itiStatusArray = cell(1, length(finalOnsets));
    for iOnset = 1:length(finalOnsets)
        onsetFrame = finalOnsets(iOnset);
        
        preOnsetFrame = max(1, onsetFrame - frame_rate);
        postOnsetFrame = min(nFrames, onsetFrame + frame_rate);
        
        if onsetFrame <= length(ITI) && ...
           ITI(preOnsetFrame) == 1 && ...
           ITI(postOnsetFrame) == 1
            itiStatusArray{iOnset} = 'ITI';
        else
            itiStatusArray{iOnset} = 'stimulus';
        end
    end
    
    onsetITIStatus{id} = itiStatusArray;
    
    if isempty(finalOnsets)
        disp(['No valid running onsets found for day ' num2str(id)]);
        data_dfof_runOnset_match{id} = nan(onsetWin*2*frame_rate, nCells);
        mean_resp_runOnset_match{id} = nan(1, nCells);
        nRunOnsets = [nRunOnsets 0];
        meanSpeedAfterOnset{id} = [];
        continue;
    end
    
    speedAfterOnset = nan(1, length(finalOnsets));
    for iOnset = 1:length(finalOnsets)
        onsetFrame = finalOnsets(iOnset);
        postOnsetWindow = onsetFrame:(onsetFrame + frame_rate - 1);
        
        postOnsetWindow = postOnsetWindow(postOnsetWindow <= length(fwdWheelClean));
        
        if ~isempty(postOnsetWindow)
            speedAfterOnset(iOnset) = mean(fwdWheelClean(postOnsetWindow), 'omitmissing');
        end
    end
    
    meanSpeedAfterOnset{id} = speedAfterOnset;
    
    data_dfof_runOnset = nan(onsetWin*2*frame_rate, nCells, length(finalOnsets));
    runConfirmation = nan(onsetWin*2*frame_rate, length(finalOnsets));
    
    for iOnset = 1:length(finalOnsets)
        onsetFrame = finalOnsets(iOnset);
        
        fullWindow = (onsetFrame - onsetWin*frame_rate):(onsetFrame + onsetWin*frame_rate - 1);
        bslnWindow = (onsetFrame - (4*frame_rate)):(onsetFrame - frame_rate/2);
        
        fullWindow = fullWindow(fullWindow > 0 & fullWindow <= nFrames);
        bslnWindow = bslnWindow(bslnWindow > 0 & bslnWindow <= nFrames);
        
        if ~isempty(fullWindow) && ~isempty(bslnWindow)
            tempFull = cellTCs_match{id}(fullWindow, :);
            tempBsln = mean(cellTCs_match{id}(bslnWindow, :));
            
            tempDfof = bsxfun(@rdivide, bsxfun(@minus, tempFull, tempBsln), tempBsln);
            
            expectedSize = onsetWin*2*frame_rate;
            if size(tempDfof, 1) == expectedSize
                data_dfof_runOnset(:, :, iOnset) = tempDfof;
            else
                tmp = nan(expectedSize, nCells);
                copySize = min(size(tempDfof, 1), expectedSize);
                tmp(1:copySize, :) = tempDfof(1:copySize, :);
                data_dfof_runOnset(:, :, iOnset) = tmp;
            end
            
            wheelWindow = fullWindow(fullWindow <= length(fwdWheelClean));
            if length(wheelWindow) == expectedSize
                runConfirmation(:, iOnset) = fwdWheelClean(wheelWindow);
            else
                tmp = nan(expectedSize, 1);
                tmp(1:length(wheelWindow)) = fwdWheelClean(wheelWindow);
                runConfirmation(:, iOnset) = tmp;
            end
        end
    end
    
    data_dfof_runOnset_match{id} = mean(data_dfof_runOnset, 3, 'omitmissing');
    nRunOnsets = [nRunOnsets length(finalOnsets)];
    mean_resp_runOnset_match{id} = mean(data_dfof_runOnset_match{id}((onsetWin*frame_rate)+1:end, :), 1, 'omitmissing');
    
    fprintf('Day %d: Found %d valid running onsets\n', id, length(finalOnsets));
    fprintf('Mean speed after onset: %.2f  %.2f cm/s\n', ...
        mean(speedAfterOnset, 'omitmissing'), std(speedAfterOnset, 'omitmissing'));
    
    itiCount = sum(strcmp(itiStatusArray, 'ITI'));
    stimCount = sum(strcmp(itiStatusArray, 'stimulus'));
    fprintf('ITI onsets: %d, Stimulus onsets: %d\n', itiCount, stimCount);
end

fprintf('\n=== SUMMARY ===\n');
fprintf('Total onsets across all days: %d\n', sum(nRunOnsets));
totalITI = 0;
totalStim = 0;
for id = 1:nd
    if ~isempty(meanSpeedAfterOnset{id})
        itiCount = sum(strcmp(onsetITIStatus{id}, 'ITI'));
        stimCount = sum(strcmp(onsetITIStatus{id}, 'stimulus'));
        totalITI = totalITI + itiCount;
        totalStim = totalStim + stimCount;
        fprintf('Day %d: %d onsets (ITI: %d, Stim: %d), mean post-onset speed = %.2f cm/s\n', ...
            id, nRunOnsets(id), itiCount, stimCount, mean(meanSpeedAfterOnset{id}, 'omitmissing'));
    else
        fprintf('Day %d: No valid onsets\n', id);
    end
end
fprintf('Overall: ITI onsets: %d, Stimulus onsets: %d\n', totalITI, totalStim);

%% save the data I want to concatenate
fn_out = fullfile(fn_multi,[num2str(stillLength),'sec_noITI']);
mkdir(fn_out)
save(fullfile(fn_out,'run_onset_analysis.mat'),'mean_resp_runOnset_match','data_dfof_runOnset_match','red_ind_match','onsetITIStatus');