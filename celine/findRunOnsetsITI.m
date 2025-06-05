function [data_dfof_runOnset, runConfirmation, nRunOnsets] = findRunOnsetsITI(cellTCs, wheel_speed, frame_rate, nOn, nOff, plotResults)
% FINDRUNONSETITI - Extract neural activity around running onsets during inter-trial intervals
%
% Inputs:
%   cellTCs - Cell time courses (frames x cells)
%   wheel_speed - Forward wheel speed (frames x 1)
%   frame_rate - Sampling rate in Hz
%   nOn - Number of frames with stimulus ON in each trial
%   nOff - Number of frames with stimulus OFF (ITI) in each trial
%   plotResults - Boolean, whether to plot sanity check (default: true)
%
% Outputs:
%   data_dfof_runOnset - Neural activity around run onsets (frames x cells)
%   runConfirmation - Wheel speed around run onsets
%   nRunOnsets - Number of valid run onsets found

if nargin < 6
    plotResults = true;
end

% Get dimensions
[nFrames, nCells] = size(cellTCs);

% Clean wheel speed data - keep only forward movement
fwdWheelClean = wheel_speed;
fwdWheelClean(fwdWheelClean < 0) = 0;

% Find running and stationary periods
stillFrames = logical(fwdWheelClean < 2); % Find stationary periods based on threshold
runFrames = logical(fwdWheelClean > 2);   % Find running periods based on threshold

% Find transitions
runOnsets = (find(diff(stillFrames) == -1) + 1); % Transitions from stationary to running
runOffsets = (find(diff(stillFrames) == 1) + 1); % Transitions from running to stationary

% Subset to onsets that are preceded by 1s stillness and followed by at least 1s of running
runOnsets_clean = runOnsets;
toDrop = [];
for iOnset = 1:length(runOnsets_clean)
    if runOnsets_clean(iOnset) < frame_rate+1
        toDrop = [toDrop, iOnset]; % Drop onsets less than 1 second from start
    elseif runOnsets_clean(iOnset) > frame_rate+1 && runOnsets_clean(iOnset)+frame_rate <= length(stillFrames)
        testWinStill = [runOnsets_clean(iOnset)-frame_rate:runOnsets_clean(iOnset)];
        testWinRun = [runOnsets_clean(iOnset):runOnsets_clean(iOnset)+frame_rate];
        
        % Make sure indices are within bounds
        testWinStill = testWinStill(testWinStill > 0 & testWinStill <= length(stillFrames));
        testWinRun = testWinRun(testWinRun > 0 & testWinRun <= length(runFrames));
        
        if sum(stillFrames(testWinStill)) < (frame_rate*.75) || sum(runFrames(testWinRun)) < (frame_rate*.75)
            toDrop = [toDrop, iOnset]; % Gather onsets that don't meet criteria
        end
    else
        toDrop = [toDrop, iOnset]; % Out of bounds indices
    end
end
runOnsets_clean(toDrop) = [];

if isempty(runOnsets_clean)
    warning('No valid running onsets found');
    data_dfof_runOnset = [];
    runConfirmation = [];
    nRunOnsets = 0;
    return;
end

% Create ITI mask - 1 during ITI, 0 during stimulus
% Using fixed trial length of nOn + nOff frames
ITI = ones(1, nFrames);
trialLength = nOn + nOff;
nTrials = floor(nFrames / trialLength);

for iTrial = 0:nTrials-1
    trialStart = iTrial * trialLength + 1;
    stimStart = trialStart + nOff; % Stimulus starts after ITI
    stimEnd = stimStart + nOn - 1;
    
    % Make sure indices are within bounds
    if stimStart > 0 && stimEnd <= nFrames
        ITI(stimStart:stimEnd) = 0;
    end
end

% Subset running onsets to only include those during ITI
ITIOnsets = runOnsets_clean;
validIndices = ITIOnsets <= length(ITI);
ITIOnsets = ITIOnsets(validIndices);

if ~isempty(ITIOnsets)
    ITIOnsets(ITI(ITIOnsets) == 0) = []; % Remove onsets during stimulus
end

if isempty(ITIOnsets)
    warning('No valid ITI running onsets found');
    data_dfof_runOnset = [];
    runConfirmation = [];
    nRunOnsets = 0;
    return;
end

% Extract neural data around run onsets
onsetWin = 1; % Window size in seconds on each side of onset

data_dfof_runOnset_all = nan(onsetWin*2*frame_rate, nCells, length(ITIOnsets));
runConfirmation_all = nan(onsetWin*2*frame_rate, length(ITIOnsets));

for iOnset = 1:length(ITIOnsets)
    % Find indices for window around onset
    fullWindow = [(ITIOnsets(iOnset)-(onsetWin*frame_rate)):(ITIOnsets(iOnset)+(onsetWin*frame_rate)-1)];
    bslnWindow = [(ITIOnsets(iOnset)-(frame_rate)):(ITIOnsets(iOnset)-(frame_rate/2))]; % 0.5s baseline before running
    
    % Make sure indices are within bounds
    fullWindow = fullWindow(fullWindow > 0 & fullWindow <= nFrames);
    bslnWindow = bslnWindow(bslnWindow > 0 & bslnWindow <= nFrames);
    
    if ~isempty(fullWindow) && ~isempty(bslnWindow)
        tempFull = cellTCs(fullWindow, :);
        tempBsln = mean(cellTCs(bslnWindow, :));
        
        % Convert to dF/F: (F-F0)/F0
        tempDfof = bsxfun(@rdivide, bsxfun(@minus, tempFull, tempBsln), tempBsln);
        
        % Handle case when window size doesn't match expected size
        if size(tempDfof, 1) == onsetWin*2*frame_rate
            data_dfof_runOnset_all(:, :, iOnset) = tempDfof;
        else
            % Pad or truncate to expected size
            tmp = nan(onsetWin*2*frame_rate, nCells);
            tmp(1:min(size(tempDfof, 1), onsetWin*2*frame_rate), :) = tempDfof(1:min(size(tempDfof, 1), onsetWin*2*frame_rate), :);
            data_dfof_runOnset_all(:, :, iOnset) = tmp;
        end
        
        % Handle wheelspeed data similarly
        if length(fullWindow) == onsetWin*2*frame_rate && all(fullWindow <= length(fwdWheelClean))
            runConfirmation_all(:, iOnset) = fwdWheelClean(fullWindow);
        else
            tmp = nan(onsetWin*2*frame_rate, 1);
            validIndices = fullWindow(fullWindow <= length(fwdWheelClean));
            tmp(1:length(validIndices)) = fwdWheelClean(validIndices);
            runConfirmation_all(:, iOnset) = tmp;
        end
    end
end

% Average across all onsets
data_dfof_runOnset = mean(data_dfof_runOnset_all, 3, 'omitmissing');
runConfirmation = mean(runConfirmation_all, 2, 'omitmissing');
nRunOnsets = length(ITIOnsets);

% Plot sanity check
if plotResults && nRunOnsets > 0
    figure; 
    plot((-(onsetWin*frame_rate):(onsetWin*frame_rate)-1)/frame_rate, runConfirmation, 'LineWidth', 1.5);
    hold on;
    plot([0 0], ylim, 'k--');
    xlabel('Time from onset (s)');
    ylabel('Wheel speed');
    title(['Average wheel speed around run onset (n = ', num2str(nRunOnsets), ' onsets)']);
    set(gca, 'TickDir', 'out');
    grid off;
    box off;
end

fprintf('Number of valid ITI running onsets: %d\n', nRunOnsets);
end
