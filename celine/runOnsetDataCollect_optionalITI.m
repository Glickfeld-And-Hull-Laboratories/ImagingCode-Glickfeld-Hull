clear all; clear global; 
close all
clc

%% Dataset and experiment setup
prompt = 'Enter ds (e.g., DART_V1_YM90K_Celine): ';
ds = input(prompt, 's');
clear prompt

rc = behavConstsDART;
eval(ds);

day_id = input('Enter day id ');
pre_day = expt(day_id).multiday_matchdays;
nd = 2;
mouse = expt(day_id).mouse;
experimentFolder = expt(day_id).exptType;

% Setup DART string and reference session
if expt(day_id).multiday_timesincedrug_hours > 0
    dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end

prompt = 'Which session was used as reference for matching: 0- baseline, 1- post-DART? ';
x = input(prompt);
switch x
    case 0
        pre = 1; post = 2;
        disp("baseline used as reference")
    case 1
        pre = 2; post = 1;
        disp("post-DART used as reference")
end

% Define how long the stationary period and running periods must be
stillLength = 5;
runLength = 2;

%% Load essential data
fn_multi = fullfile(rc.achAnalysis, experimentFolder, mouse, ['multiday_' dart_str]);
cd(fn_multi)

% Load matched cell timecourses and input structure
load(fullfile(fn_multi, 'timecourses.mat'), 'cellTCs_match','red_ind_match');
load(fullfile(fn_multi, 'input.mat'));

% Rename input to avoid conflict with MATLAB's input function
inputStructure = input;
clear input x prompt ds

% Get frame rate
frame_rate = inputStructure(1).frameImagingRateMs;

%% Load wheel speed data for all days
allDays = [day_id, pre_day];
wheel_speed = cell(1, nd);

for id = 1:nd
    wheel_speed{id} = wheelSpeedCalc(inputStructure(id), 32, expt(allDays(1)).wheelColor);
end

% Clean wheel speed (set small movements to zero)
wheel_speed_clean = cell(1, nd);
for id = 1:nd
    wheel_speed_clean{id} = wheel_speed{id};
    wheel_speed_clean{id}(abs(wheel_speed_clean{id}) < 5.37) = 0;
end

% Extract stimulus properties (needed for ITI calculation)
nOn = inputStructure(1).nScansOn;
nOff = inputStructure(1).nScansOff;

% Calculate number of trials for each day
nTrials = zeros(1, nd);
for id = 1:nd
    nTrials(id) = size(cellTCs_match{id}, 1) / (nOn + nOff);
end

% Summary of loaded variables
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
wheel_speed = cell(1,nd);

for id = 1:nd
    wheel_speed{id} = wheelSpeedCalc(inputStructure(id),32,expt(allDays(1)).wheelColor); 
    nanmean(wheel_speed{id})
end
wheel_speed_clean = cell(1,nd);
for id = 1:nd
    wheel_speed_clean{id}=wheel_speed{id};
    wheel_speed_clean{id}(abs(wheel_speed_clean{id})<5.37)=0;
end

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
    RIx{id} = wheel_trial_avg{id}>2; %~5 is the noise level in the wheel movement
    mean(RIx{id})
end

%% Extract running onsets
frame_rate = inputStructure.frameImagingRateMs;

speed_threshold = 5.37; % cm/s threshold for running

% Initialize variables
nRunOnsets = [];
meanSpeedAfterOnset = []; % Store mean speed after each onset
onsetITIStatus = cell(1, nd); % Store ITI status for each onset

for id = 1:nd
    mouse = expt(allDays(id)).mouse;
    date = expt(allDays(id)).date;
    imgFolder = expt(allDays(id)).contrastxori_runs{1};
    imgMatFile = [imgFolder '_000_000.mat'];
    dataPath = fullfile(rc.achData, mouse, date, imgFolder);
    load(fullfile(dataPath,imgMatFile));
    if isfield(info, 'frame')
        [cStimOn stimOffs] = photoFrameFinder_Sanworks(info.frame);
        fprintf('Getting cStimOn from photo diode\n');
    elseif exist(fullfile(fn_multi, 'correctedInputStructure.mat')) %if there is already a corrected input structure
        load("correctedInputStructure.mat");
        cStimOn = cell2mat(correctedInputStructure{id}.cStimOn); %use the corrected cStimOn from that
        fprintf('Getting cStimOn from existing corrected input structure\n');
    else %if there is no info.frame and no corrected input structure
        correctedInputStructure_temp = counterValCorrect_noPhotodiode(inputStructure(id));
        cStimOn = cell2mat(correctedInputStructure_temp.cStimOn);
        fprintf('Creating corrected input structure to get cStimOn\n');
    end
    [nFrames, nCells] = size(cellTCs_match{id});
    
    fwdWheelClean = wheel_speed_clean{id};
    fwdWheelClean(fwdWheelClean<0) = 0;
    
    % Find and refine the onset indices
    stillFrames = logical(fwdWheelClean < speed_threshold);
    runFrames = logical(fwdWheelClean >= speed_threshold); % Use >= for consistency
    
    runOnsets = find(diff(stillFrames) == -1) + 1;
    runOffsets = find(diff(stillFrames) == 1) + 1;
    
    % Subset to onsets that meet our criteria
    runOnsets_clean = runOnsets;
    toDrop = [];
    
    for iOnset = 1:length(runOnsets_clean)
        onsetFrame = runOnsets_clean(iOnset);
        
        % Check if onset is too close to start or end of session
        if onsetFrame < frame_rate + 1 || ...
           onsetFrame <= (frame_rate * stillLength) || ...
           onsetFrame + (frame_rate * runLength) > length(stillFrames)
            toDrop = [toDrop, iOnset];
            continue;
        end
        
        % Define test windows
        testWinStill = (onsetFrame - stillLength*frame_rate):(onsetFrame - 1);
        testWinRun = (onsetFrame + 1):(onsetFrame + runLength*frame_rate);
        
        % Check stillness requirement (ALL frames in stillness window should be still)
        stillnessReq = sum(stillFrames(testWinStill)) == length(testWinStill);
        
        % Check running requirement (at least 50% of frames should be running)
        runningReq = sum(runFrames(testWinRun)) >= (runLength * frame_rate * 0.5);
        
        if ~stillnessReq || ~runningReq
            toDrop = [toDrop, iOnset];
        end
    end
    
    runOnsets_clean(toDrop) = [];
    
    if isempty(runOnsets_clean)
        disp(['No valid running onsets found for day ' num2str(id)]);
        data_dfof_runOnset_match{id} = [];
        mean_resp_runOnset_match{id} = [];
        nRunOnsets = [nRunOnsets 0];
        meanSpeedAfterOnset{id} = [];
        onsetITIStatus{id} = [];
        continue;
    end
    
    % Create ITI vector for documentation (not filtering)
    ITI = ones(1, nFrames);
    for iTrial = 1:inputStructure(id).trialSinceReset
        if cStimOn(iTrial) > 0 && cStimOn(iTrial) + inputStructure(id).nScansOn <= nFrames
            ITI(cStimOn(iTrial):cStimOn(iTrial) + inputStructure(id).nScansOn) = 0;
        end
    end
    
    % Use ALL cleaned onsets (no ITI filtering)
    finalOnsets = runOnsets_clean;
    
    % Document ITI status for each onset
    itiStatusArray = cell(1, length(finalOnsets));
    for iOnset = 1:length(finalOnsets)
        onsetFrame = finalOnsets(iOnset);
        
        % Check if onset and surrounding frames are within ITI
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
    
    % Store ITI status
    onsetITIStatus{id} = itiStatusArray;
    
    if isempty(finalOnsets)
        disp(['No valid running onsets found for day ' num2str(id)]);
        data_dfof_runOnset_match{id} = [];
        mean_resp_runOnset_match{id} = [];
        nRunOnsets = [nRunOnsets 0];
        meanSpeedAfterOnset{id} = [];
        continue;
    end
    
    % Calculate mean speed in 1-second period following each onset
    speedAfterOnset = nan(1, length(finalOnsets));
    for iOnset = 1:length(finalOnsets)
        onsetFrame = finalOnsets(iOnset);
        postOnsetWindow = onsetFrame:(onsetFrame + frame_rate - 1);
        
        % Make sure we don't go beyond the data
        postOnsetWindow = postOnsetWindow(postOnsetWindow <= length(fwdWheelClean));
        
        if ~isempty(postOnsetWindow)
            speedAfterOnset(iOnset) = mean(fwdWheelClean(postOnsetWindow), 'omitmissing');
        end
    end
    
    % Store the mean speeds for this day
    meanSpeedAfterOnset{id} = speedAfterOnset;
    
    % Get neural data for run onsets 
    onsetWin = 2; % seconds on either side of onset
    
    data_dfof_runOnset = nan(onsetWin*2*frame_rate, nCells, length(finalOnsets));
    runConfirmation = nan(onsetWin*2*frame_rate, length(finalOnsets));
    
    for iOnset = 1:length(finalOnsets)
        onsetFrame = finalOnsets(iOnset);
        
        % Define windows
        fullWindow = (onsetFrame - onsetWin*frame_rate):(onsetFrame + onsetWin*frame_rate - 1);
        bslnWindow = (onsetFrame - frame_rate/2):onsetFrame;
        
        % Ensure indices are within bounds
        fullWindow = fullWindow(fullWindow > 0 & fullWindow <= nFrames);
        bslnWindow = bslnWindow(bslnWindow > 0 & bslnWindow <= nFrames);
        
        if ~isempty(fullWindow) && ~isempty(bslnWindow)
            tempFull = cellTCs_match{id}(fullWindow, :);
            tempBsln = mean(cellTCs_match{id}(bslnWindow, :));
            
            % Convert to df/f
            tempDfof = bsxfun(@rdivide, bsxfun(@minus, tempFull, tempBsln), tempBsln);
            
            % Handle size mismatches
            expectedSize = onsetWin*2*frame_rate;
            if size(tempDfof, 1) == expectedSize
                data_dfof_runOnset(:, :, iOnset) = tempDfof;
            else
                tmp = nan(expectedSize, nCells);
                copySize = min(size(tempDfof, 1), expectedSize);
                tmp(1:copySize, :) = tempDfof(1:copySize, :);
                data_dfof_runOnset(:, :, iOnset) = tmp;
            end
            
            % Handle wheel speed data
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
    
    % Store results
    data_dfof_runOnset_match{id} = mean(data_dfof_runOnset, 3, 'omitmissing');
    nRunOnsets = [nRunOnsets length(finalOnsets)];
    mean_resp_runOnset_match{id} = mean(data_dfof_runOnset_match{id}((onsetWin*frame_rate)+1:end, :), 1, 'omitmissing');
    
    % Display results for this day
    fprintf('Day %d: Found %d valid running onsets\n', id, length(finalOnsets));
    fprintf('Mean speed after onset: %.2f � %.2f cm/s\n', ...
        mean(speedAfterOnset, 'omitmissing'), std(speedAfterOnset, 'omitmissing'));
    
    % Display ITI status breakdown
    itiCount = sum(strcmp(itiStatusArray, 'ITI'));
    stimCount = sum(strcmp(itiStatusArray, 'stimulus'));
    fprintf('ITI onsets: %d, Stimulus onsets: %d\n', itiCount, stimCount);
end

% Summary across all days
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
% % Additional Plots
% onsetLocations = cell(1, nd);
% onsetReasons = cell(1, nd); % Store reason for each onset
% 
% % We need to re-run the onset detection to store onset locations and reasons
% for id = 1:nd
%     mouse = expt(allDays(id)).mouse;
%     date = expt(allDays(id)).date;
%     imgFolder = expt(allDays(id)).contrastxori_runs{1};
%     imgMatFile = [imgFolder '_000_000.mat'];
%     dataPath = fullfile(rc.achData, mouse, date, imgFolder);
%     load(fullfile(dataPath,imgMatFile));
%     [cStimOn stimOffs] = photoFrameFinder_Sanworks(info.frame);
%     [nFrames, nCells] = size(cellTCs_match{id});
% 
%     fwdWheelClean = wheel_speed_clean{id};
%     fwdWheelClean(fwdWheelClean<0) = 0;
% 
%     stillFrames = logical(fwdWheelClean < speed_threshold);
%     runFrames = logical(fwdWheelClean >= speed_threshold);
%     runOnsets = find(diff(stillFrames) == -1) + 1;
% 
%     % Store ALL original onsets and track their fate
%     allOnsets = runOnsets;
%     onsetStatus = cell(1, length(allOnsets)); % Track why each onset was kept/rejected
% 
%     % Create ITI vector for documentation
%     ITI = ones(1, nFrames);
%     for iTrial = 1:inputStructure(id).trialSinceReset
%         if cStimOn(iTrial) > 0 && cStimOn(iTrial) + inputStructure(id).nScansOn <= nFrames
%             ITI(cStimOn(iTrial):cStimOn(iTrial) + inputStructure(id).nScansOn) = 0;
%         end
%     end
% 
%     % Apply filtering and track reasons
%     for iOnset = 1:length(allOnsets)
%         onsetFrame = allOnsets(iOnset);
% 
%         % Check boundary conditions
%         if onsetFrame < frame_rate + 1
%             onsetStatus{iOnset} = 'too_close_to_start';
%             continue;
%         elseif onsetFrame <= (frame_rate * stillLength)
%             onsetStatus{iOnset} = 'insufficient_pre_data';
%             continue;
%         elseif onsetFrame + (frame_rate * runLength) > length(stillFrames)
%             onsetStatus{iOnset} = 'insufficient_post_data';
%             continue;
%         end
% 
%         % Check stillness and running requirements
%         testWinStill = (onsetFrame - stillLength*frame_rate):(onsetFrame - 1);
%         testWinRun = (onsetFrame + 1):(onsetFrame + runLength*frame_rate);
% 
%         stillnessReq = sum(stillFrames(testWinStill)) == length(testWinStill);
%         runningReq = sum(runFrames(testWinRun)) >= (runLength * frame_rate * 0.5);
% 
%         if ~stillnessReq && ~runningReq
%             onsetStatus{iOnset} = 'failed_both_criteria';
%         elseif ~stillnessReq
%             onsetStatus{iOnset} = 'insufficient_stillness';
%         elseif ~runningReq
%             onsetStatus{iOnset} = 'insufficient_running';
%         else
%             % Passed basic criteria, now document ITI status
%             % Check ITI status for documentation
%             preOnsetFrame = max(1, onsetFrame - frame_rate);
%             postOnsetFrame = min(nFrames, onsetFrame + frame_rate);
% 
%             if onsetFrame <= length(ITI) && ...
%                ITI(preOnsetFrame) == 1 && ...
%                ITI(postOnsetFrame) == 1
%                 onsetStatus{iOnset} = 'kept_ITI';
%             else
%                 onsetStatus{iOnset} = 'kept_stimulus';
%             end
%         end
%     end
% 
%     onsetLocations{id} = allOnsets;
%     onsetReasons{id} = onsetStatus;
% end
% 
% % Create output folder name (always noITI since filtering is disabled)

% 
% %%
% % 1) Plot wheel speed with ALL onsets color-coded by rejection reason and ITI status
% figure('Position', [100, 100, 1000, 600]);
% 
% % Define colors for different rejection reasons and ITI status
% colorMap = containers.Map();
% colorMap('kept_ITI') = [0, 0.8, 0]; % Green - kept onsets during ITI
% colorMap('kept_stimulus') = [0, 0.6, 0.8]; % Cyan - kept onsets during stimulus
% colorMap('too_close_to_start') = [1, 0, 1]; % Magenta
% colorMap('insufficient_pre_data') = [1, 0.5, 0]; % Orange
% colorMap('insufficient_post_data') = [0.5, 0.5, 0.5]; % Gray
% colorMap('insufficient_stillness') = [1, 0, 0]; % Red
% colorMap('insufficient_running') = [0, 0, 1]; % Blue
% colorMap('failed_both_criteria') = [0.5, 0, 0]; % Dark red
% 
% for id = 1:nd
%     subplot(nd, 1, id);
% 
%     % Plot wheel speed
%     timeVec = (1:length(wheel_speed_clean{id})) / frame_rate_double;
%     plot(timeVec, wheel_speed_clean{id}, 'k-', 'LineWidth', 0.5);
%     hold on;
% 
%     % Plot all onset locations with color coding
%     if ~isempty(onsetLocations{id})
%         onsetTimes = onsetLocations{id} / frame_rate_double;
%         yLims = [0, max(wheel_speed_clean{id})];
% 
%         for iOnset = 1:length(onsetTimes)
%             onsetTime = onsetTimes(iOnset);
%             reason = onsetReasons{id}{iOnset};
%             color = colorMap(reason);
% 
%             if startsWith(reason, 'kept')
%                 % For kept onsets, show shaded region + line
%                 endTime = onsetTime + 1; % 1 second after onset
%                 fill([onsetTime, endTime, endTime, onsetTime], ...
%                      [yLims(1), yLims(1), yLims(2), yLims(2)], ...
%                      color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%                 plot([onsetTime, onsetTime], yLims, '-', 'Color', color, 'LineWidth', 3);
%             else
%                 % For rejected onsets, show just a line
%                 plot([onsetTime, onsetTime], yLims, '-', 'Color', color, 'LineWidth', 2);
%             end
%         end
%     end
% 
%     xlabel('Time (s)');
%     ylabel('Wheel Speed (cm/s)');
%     title(['Day ' num2str(id) ' - Total: ' num2str(length(onsetLocations{id})) ', Kept: ' num2str(nRunOnsets(id))]);
%     set(gca, 'TickDir', 'out');
%     grid off;
%     box off;
% 
%     % Add legend only to first subplot
%     if id == 1
%         % Count onsets by reason for legend across all days
%         allReasons = [];
%         for dayId = 1:nd
%             if ~isempty(onsetReasons{dayId})
%                 allReasons = [allReasons, onsetReasons{dayId}];
%             end
%         end
% 
%         uniqueReasons = unique(allReasons);
%         legendHandles = [];
%         legendEntries = {};
% 
%         for iReason = 1:length(uniqueReasons)
%             reason = uniqueReasons{iReason};
%             count = sum(strcmp(allReasons, reason));
%             color = colorMap(reason);
% 
%             % Create a proper line handle for the legend
%             h = plot(NaN, NaN, '-', 'Color', color, 'LineWidth', 3);
%             legendHandles(end+1) = h;
% 
%             % Format legend entry
%             if strcmp(reason, 'kept_ITI')
%                 legendEntries{end+1} = ['kept (ITI) (' num2str(count) ')'];
%             elseif strcmp(reason, 'kept_stimulus')
%                 legendEntries{end+1} = ['kept (stimulus) (' num2str(count) ')'];
%             else
%                 legendEntries{end+1} = [strrep(reason, '_', ' ') ' (' num2str(count) ')'];
%             end
%         end
% 
%         legend(legendHandles, legendEntries, 'Location', 'best', 'FontSize', 8);
%     end
% end
% 
% sgtitle('All Detected Onsets Color-Coded by Filtering Result (ITI Status Documented)');
% 
% print(fullfile(fn_out,'runOnsetDetect.pdf'),'-dpdf');
% 
% % Print summary of rejection reasons
% fprintf('\n=== ONSET FILTERING SUMMARY ===\n');
% fprintf('ITI filtering: disabled (ITI status documented for kept onsets)\n');
% for id = 1:nd
%     fprintf('\nDay %d:\n', id);
%     reasons = onsetReasons{id};
%     uniqueReasons = unique(reasons);
% 
%     for iReason = 1:length(uniqueReasons)
%         reason = uniqueReasons{iReason};
%         count = sum(strcmp(reasons, reason));
% 
%         % Format display
%         if strcmp(reason, 'kept_ITI')
%             fprintf('  kept (ITI): %d\n', count);
%         elseif strcmp(reason, 'kept_stimulus')
%             fprintf('  kept (stimulus): %d\n', count);
%         else
%             fprintf('  %s: %d\n', reason, count);
%         end
%     end
%     fprintf('  Total detected: %d\n', length(reasons));
% end
% 
% %%
% % Initialize arrays to store table data
% dayCol = [];
% reasonCol = {};
% countCol = [];
% 
% % Collect data for table
% for id = 1:nd
%     reasons = onsetReasons{id};
%     uniqueReasons = unique(reasons);
% 
%     for iReason = 1:length(uniqueReasons)
%         reason = uniqueReasons{iReason};
%         count = sum(strcmp(reasons, reason));
% 
%         % Format reason for table
%         if strcmp(reason, 'kept_ITI')
%             formattedReason = 'kept (ITI)';
%         elseif strcmp(reason, 'kept_stimulus')
%             formattedReason = 'kept (stimulus)';
%         else
%             formattedReason = reason;
%         end
% 
%         % Append to arrays
%         dayCol = [dayCol; id];
%         reasonCol{end+1} = formattedReason;
%         countCol = [countCol; count];
%     end
% end
% 
% % Create table
% onsetSummaryTable = table(dayCol, reasonCol', countCol, ...
%     'VariableNames', {'Day', 'Reason', 'Count'});
% 
% % Display table summary
% fprintf('\n=== ONSET FILTERING SUMMARY TABLE ===\n');
% disp(onsetSummaryTable);
% 
% % Save table to file
% writetable(onsetSummaryTable, fullfile(fn_out, 'onset_filtering_summary.csv'));
% save(fullfile(fn_out, 'onset_filtering_summary.mat'), 'onsetSummaryTable');

