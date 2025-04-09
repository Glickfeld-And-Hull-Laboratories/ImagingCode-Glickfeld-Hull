%% EPSC Detection in Whole Cell Recordings (.abf format)
% This script detects excitatory postsynaptic currents (EPSCs) in 
% whole cell voltage clamp recordings from .abf files.
% Handles multiple sweeps (time × channel × sweep data structure).
% Modified to prevent overlapping events and start detection after 200 ms.

clear all; close all; clc;

%% Parameters for EPSC detection
params.minAmplitude = 5;       % Minimum amplitude for EPSC detection (pA)
params.maxAmplitude = 100;      % Maximum amplitude for EPSC detection (pA)
params.minSlope = 5;            % Minimum slope for EPSC onset (pA/ms)
params.minWidth = 2;            % Minimum width of EPSC (ms)
params.baseline = 'pre';        % Method for baseline calculation ('pre', 'local', or 'median')
params.smoothWindow = 0.2;      % Window size for smoothing (ms)
params.channel = 1;             % Channel to analyze (default: 1)
params.inverted = true;         % Set to true for inward currents (negative deflections)
params.conversionFactor = 1;    % Some recordings had units other than pA 
params.plotStyle.tickDir = 'out'; % Direction of tick marks
params.plotStyle.grid = 'off';  % Grid off
params.plotStyle.box = 'off';   % Box off
params.startTime = 200;         % Start detection after this time (ms)

%% Step 1: Import ABF file
% Using MATLAB's built-in abfload function
% If not available, please download from the MathWorks File Exchange
% or use the function provided at the end of this script

% Define the filename directly or use dialog if not provided
% To use directly, specify the full path to the file
fileNameInput = '24110000'; % Enter your filename here, e.g., 'C:\Data\recording.abf'
fileNameInput = [fileNameInput,'.abf']
if isempty(fileNameInput)
    % Use dialog if no filename is provided
    [fileName, filePath] = uigetfile('*.abf', 'Select ABF file');
    if fileName == 0
        error('No file selected');
    end
    fullPath = fullfile(filePath, fileName);
else
    % Use the direct filename input
    fullPath = fileNameInput;
    [filePath, fileName, fileExt] = fileparts(fullPath);
    fileName = [fileName, fileExt]; % Add extension back to filename
end

try
    [data, si, h] = abfload(fullPath);
catch
    error(['Could not load ABF file. Make sure abfload.m is in your path. ' ...
           'Alternative: try using importdata for simple .abf files.']);
end

% Extract data and convert sampling interval to ms
samplingInterval = si/1000; % Convert to ms
timeVector = (0:size(data,1)-1)' * samplingInterval;

% Check data dimensions
[numTimePoints, numChannels, numSweeps] = size(data);
fprintf('File loaded: %d time points, %d channels, %d sweeps\n', numTimePoints, numChannels, numSweeps);
fprintf('File: %s\n', fileName);

% Apply conversion factor
data = data * params.conversionFactor;

% Calculate the index corresponding to the start time
startIdx = find(timeVector >= params.startTime, 1);
if isempty(startIdx)
    error('Start time exceeds recording duration');
end
fprintf('Starting analysis at %.2f ms (index %d)\n', timeVector(startIdx), startIdx);

% Calculate effective duration for frequency calculation (total time - start time)
effectiveDuration = (timeVector(end) - params.startTime) / 1000; % in seconds
fprintf('Effective duration for rate calculation: %.2f seconds\n', effectiveDuration);

%% Step 2: Process each sweep
% Initialize results structure for each sweep
allResults = struct('sweep', {}, 'epscs', {}, 'stats', {});

% Process each sweep individually
for sweepIdx = 1:numSweeps
    fprintf('\nProcessing sweep %d of %d...\n', sweepIdx, numSweeps);
    
    % Extract current data for this sweep
    currentData = data(:, params.channel, sweepIdx);
    
    % If inverted (inward currents are negative), flip the sign
    if params.inverted
        currentData = -currentData; % Invert so EPSCs appear as positive deflections
    end
    
    % Apply median filter to remove artifacts
    currentFiltered = medfilt1(currentData, 5);
    
    % Smooth the data
    smoothWindow = ceil(params.smoothWindow / samplingInterval);
    currentSmoothed = smooth(currentFiltered, smoothWindow);
    
    % Calculate the first derivative (slope)
    currentSlope = diff(currentSmoothed) / samplingInterval;
    % Ensure currentSlope is a column vector and pad to match original length
    if size(currentSlope, 2) > 1
        % If currentSlope is a row vector or matrix, transpose it
        currentSlope = currentSlope';
    end
    currentSlope = [NaN(1, size(currentSlope, 2))'; currentSlope]; % Pad to match original length

    %% Step 3: Detect potential EPSCs for this sweep
    % Find local maxima in the slope (potential EPSC onsets), starting after the specified start time
    [~, allPotentialOnsets] = findpeaks(currentSlope, 'MinPeakHeight', params.minSlope, 'MinPeakDistance', 10);
    
    % Filter out onsets that occur before the start time
    potentialOnsets = allPotentialOnsets(allPotentialOnsets >= startIdx);
    
    % Initialize arrays for EPSC detection
    epscs = struct('onset', {}, 'peak', {}, 'offset', {}, 'amplitude', {}, 'slope', {}, 'halfWidth', {}, 'area', {});
    epscCount = 0;
    
    % Keep track of detected event intervals
    detectedIntervals = []; % Will store [onset, offset] pairs
    
    % For each potential onset, validate if it's an EPSC
    for i = 1:length(potentialOnsets)
        onsetIdx = potentialOnsets(i);
        
        % Skip if too close to the end of recording
        if onsetIdx + 100 > length(currentSmoothed)
            continue;
        end
        
        % Check if this onset falls within any previously detected event
        overlaps = false;
        for j = 1:size(detectedIntervals, 1)
            if onsetIdx >= detectedIntervals(j, 1) && onsetIdx <= detectedIntervals(j, 2)
                overlaps = true;
                break;
            end
        end
        
        % Skip this potential event if it overlaps with an existing one
        if overlaps
            continue;
        end
        
        % Determine baseline
        baselineStart = 2.5/samplingInterval; % baseline will start 25ms before the event
        baselineEnd = 0.5/samplingInterval; % baseline will end 5 ms before the event
        if strcmp(params.baseline, 'pre')
            if onsetIdx < baselineStart
                continue; % Skip if too close to the beginning
            end
            baselineWindow = max(1, onsetIdx-baselineStart):onsetIdx-baselineEnd;
            baseline = mean(currentSmoothed(baselineWindow));
        elseif strcmp(params.baseline, 'local')
            baseline = currentSmoothed(onsetIdx);
        else % median
            baselineWindow = max(1, onsetIdx-50):min(length(currentSmoothed), onsetIdx+50);
            baseline = median(currentSmoothed(baselineWindow));
        end
        
        % Find peak within a reasonable window after onset, here more than
        % 1ms after onset
        [peakValue, peakOffset] = max(currentSmoothed(onsetIdx:min(onsetIdx+(1/samplingInterval), length(currentSmoothed))));
        peakIdx = onsetIdx + peakOffset - 1;
        
        % Calculate amplitude
        amplitude = peakValue - baseline;
        
        % Check if amplitude meets threshold
        if amplitude < params.minAmplitude || amplitude > params.maxAmplitude
            continue;
        end
        
        % Find offset (return to baseline or next EPSC onset)
        offsetIdx = peakIdx;
        while offsetIdx < length(currentSmoothed) && ...
              currentSmoothed(offsetIdx) > (baseline + 0.2*amplitude)
            offsetIdx = offsetIdx + 1;
        end
        
        % Calculate half-width
        halfAmp = baseline + amplitude/2;
        
        % Find first crossing of half-amplitude
        halfStartIdx = onsetIdx;
        while halfStartIdx < peakIdx && currentSmoothed(halfStartIdx) < halfAmp
            halfStartIdx = halfStartIdx + 1;
        end
        
        % Find second crossing of half-amplitude
        halfEndIdx = peakIdx;
        while halfEndIdx < offsetIdx && currentSmoothed(halfEndIdx) > halfAmp
            halfEndIdx = halfEndIdx + 1;
        end
        
        % Calculate half-width in ms
        halfWidth = (halfEndIdx - halfStartIdx) * samplingInterval;
        
        % Check if width meets threshold
        if halfWidth < params.minWidth
            continue;
        end
        
        % Calculate area under the curve (charge transfer)
        area = trapz(timeVector(onsetIdx:offsetIdx), currentSmoothed(onsetIdx:offsetIdx) - baseline) * samplingInterval;
        
        % Store the EPSC
        epscCount = epscCount + 1;
        epscs(epscCount).onset = onsetIdx;
        epscs(epscCount).peak = peakIdx;
        epscs(epscCount).offset = offsetIdx;
        epscs(epscCount).amplitude = amplitude;
        epscs(epscCount).slope = currentSlope(onsetIdx);
        epscs(epscCount).halfWidth = halfWidth;
        epscs(epscCount).area = area;
        epscs(epscCount).time = timeVector(onsetIdx); % Store actual time in ms
        if exist('tau', 'var')
            epscs(epscCount).tau = tau;
        end
        
        % Add this event to the list of detected intervals
        detectedIntervals(end+1, :) = [onsetIdx, offsetIdx];
    end
    
    fprintf('Detected %d EPSCs in sweep %d (after %.1f ms)\n', epscCount, sweepIdx, params.startTime);
    
    % Store results for this sweep
    allResults(sweepIdx).sweep = sweepIdx;
    allResults(sweepIdx).epscs = epscs;
    
    % Calculate sweep statistics
    if epscCount > 0
        allResults(sweepIdx).stats.meanAmplitude = mean([epscs.amplitude]);
        allResults(sweepIdx).stats.meanHalfWidth = mean([epscs.halfWidth]);
        % Calculate frequency based on effective duration (after start time)
        allResults(sweepIdx).stats.frequency = epscCount/effectiveDuration; % in Hz
        if isfield(epscs, 'tau') && ~all(isnan([epscs.tau]))
            validTaus = [epscs.tau];
            allResults(sweepIdx).stats.meanTau = mean(validTaus(~isnan(validTaus)));
        end
    else
        allResults(sweepIdx).stats = struct();
        allResults(sweepIdx).stats.meanAmplitude = NaN;
        allResults(sweepIdx).stats.meanHalfWidth = NaN;
        allResults(sweepIdx).stats.frequency = 0;
    end
    
    % Store start time information
    allResults(sweepIdx).analysisStartTime = params.startTime;
    allResults(sweepIdx).effectiveDuration = effectiveDuration;
end
