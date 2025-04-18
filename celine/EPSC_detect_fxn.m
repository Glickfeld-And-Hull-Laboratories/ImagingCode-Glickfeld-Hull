function [] = EPSC_detect_fxn(filename,outpath)
%% EPSC Detection in Whole Cell Recordings (.abf format)
% This script detects excitatory postsynaptic currents (EPSCs) in 
% whole cell voltage clamp recordings from .abf files.
% Handles multiple sweeps (time  channel  sweep data structure).
% Modified to prevent overlapping events and start detection after 200 ms.

close all; clc;
cd('\home\celine\Data\patch_data')
%% Parameters for EPSC detection
params.minAmplitude = 5;       % Minimum amplitude for EPSC detection (pA)
params.maxAmplitude = 100;      % Maximum amplitude for EPSC detection (pA)
params.minSlope = 2;            % Minimum slope for EPSC onset (pA/ms)
params.minWidth = 1;            % Minimum width of EPSC (ms)
params.baseline = 'pre';        % Method for baseline calculation ('pre', 'local', or 'median')
params.smoothWindow = 1;       % Window size for smoothing (ms)
params.channel = 1;             % Channel to analyze (default: 1)
params.inverted = true;         % Set to true for inward currents (negative deflections)
params.plotStyle.tickDir = 'out'; % Direction of tick marks
params.plotStyle.grid = 'off';  % Grid off
params.plotStyle.box = 'off';   % Box off
params.startTime = 350;         % Start detection after this time (ms)

%% Step 1: Import ABF file
% Using MATLAB's built-in abfload function
% If not available, please download from the MathWorks File Exchange
% or use the function provided at the end of this script

% Define the filename directly or use dialog if not provided
% To use directly, specify the full path to the file

filePath = [filename,'.abf']
[data, si, h] = abfload(filePath);

if h.recChUnits{1,1}=='pA'
    params.conversionFactor = 1;  
elseif h.recChUnits{1,1}=='nA'
    params.conversionFactor = 1000;
else
    error('Units are not nA or pA')
end

% Extract data and convert sampling interval to ms
samplingInterval = si/1000; % Convert to ms
samplesPerMs = 1/samplingInterval;
timeVector = (0:size(data,1)-1)' * samplingInterval;
Fs = samplesPerMs*1000;
% Check data dimensions
[numTimePoints, numChannels, numSweeps] = size(data);
fprintf('File loaded: %d time points, %d channels, %d sweeps\n', numTimePoints, numChannels, numSweeps);
fprintf('File: %s\n', filename);

% Apply conversion factor
data = data * params.conversionFactor;
%data = data * h.fInstrumentScaleFactor(1);

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
firstPlotMade=false
% Process each sweep individually
for sweepIdx = 1:numSweeps
    fprintf('\nProcessing sweep %d of %d...\n', sweepIdx, numSweeps);
    
    % Extract current data for this sweep
    currentDataRaw = data(:, params.channel, sweepIdx);
%% frequency filtering
    cutoffFreq = 1300; % Cutoff frequency for lowpass filter (1.3 kHz)
    notchFreq = 60;    % Frequency to remove (60 Hz)
    order = 4;         % Filter order for lowpass filter
    
    % Design and apply lowpass filter
    Wn_low = cutoffFreq/(Fs/2);
    [b_low, a_low] = butter(order, Wn_low, 'low');
    lowpassData = filtfilt(b_low, a_low, currentDataRaw);
    
    % Design notch filter to remove 60 Hz
    notchQ = 35;       % Quality factor - higher Q = narrower notch
    wo = notchFreq/(Fs/2);  % Normalized frequency
    bw = wo/notchQ;         % Bandwidth
    [b_notch, a_notch] = iirnotch(wo, bw);
    
    % Apply notch filter to the already lowpass filtered data
    currentData = filtfilt(b_notch, a_notch, lowpassData);
    % the median filter, smoothing, etc will move on with the current data
    % that has done through the frequency filtering.

    
%% Check for baseline drift
% Find the index corresponding to the start time
startTimeIdx = find(timeVector >= params.startTime, 1);
if isempty(startTimeIdx)
    error('Start time exceeds recording duration');
end

% Only analyze from the start time onwards
filteredDataForAnalysis = currentDataRaw(startTimeIdx:end);
timeVectorForAnalysis = timeVector(startTimeIdx:end);

% Get the start portion after the startTime
startPortion = filteredDataForAnalysis(1:min(round(0.1*length(filteredDataForAnalysis)), 1000));
endPortion = filteredDataForAnalysis(max(1, end-min(round(0.1*length(filteredDataForAnalysis)), 1000)):end);

% Calculate mean values
startMean = mean(startPortion);
endMean = mean(endPortion);

% Check if baseline drift exceeds threshold
baselineDriftCheck = abs(startMean - endMean) <= 5; % pA

% Also check if any extended portion differs from baseline by more than 5 pA
extendedDeviationCheck = true;
windowSize = round(200 / samplingInterval); % 200 ms window converted to samples
baselineThreshold = 5; % pA

if windowSize < length(filteredDataForAnalysis)
    % Use a sliding window approach to check for extended deviations
    for i = 1:windowSize:length(filteredDataForAnalysis)-windowSize
        windowMean = mean(filteredDataForAnalysis(i:i+windowSize-1));
        if abs(windowMean - startMean) > baselineThreshold
            extendedDeviationCheck = false;
            break;
        end
    end
end

% Skip this sweep if either check fails
if ~baselineDriftCheck || ~extendedDeviationCheck
    fprintf('Skipping sweep %d due to baseline instability after start time (%.2f ms)\n', sweepIdx, params.startTime);
    
    % Store empty results for this sweep with NaN values
    allResults(sweepIdx).sweep = sweepIdx;
    allResults(sweepIdx).epscs = struct('onset', {}, 'peak', {}, 'offset', {}, 'amplitude', {}, 'slope', {}, 'halfWidth', {}, 'area', {});
    allResults(sweepIdx).stats = struct();
    allResults(sweepIdx).stats.meanAmplitude = NaN;
    allResults(sweepIdx).stats.meanHalfWidth = NaN;
    allResults(sweepIdx).stats.frequency = NaN;
    allResults(sweepIdx).analysisStartTime = params.startTime;
    allResults(sweepIdx).effectiveDuration = effectiveDuration;
    allResults(sweepIdx).baselineStable = false;
    
    % Skip to the next sweep
    continue;
else
    allResults(sweepIdx).baselineStable = true;
end

    
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

    %figure;plot(currentDataRaw*-1);hold on;plot(currentData);hold on;plot(currentFiltered);hold on;plot(currentSmoothed)

    %% Step 3: Detect potential EPSCs for this sweep
    % Find local maxima in the slope (potential EPSC onsets), starting after the specified start time
    [~, allPotentialOnsets] = findpeaks(currentSlope, 'MinPeakHeight', params.minSlope, 'MinPeakDistance', 20);
    
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
        baselineStart = 10*samplesPerMs; % baseline will start 20ms before the event
        baselineEnd = 5*samplesPerMs; % baseline will end 5 ms before the event
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
        
        %% find peaks
        % Find peak based on derivative zero-crossing (slope reversal)
        % We already have currentSlope calculated from earlier
        
        % Find where the derivative changes from positive to negative
        % (indicating a switch from rising to falling)
        zeroIdx = []; % Will store indices where slope crosses zero
        for j = onsetIdx+2:min(onsetIdx+(15*samplesPerMs), length(currentSlope)-1)
            % Check if slope crosses from positive to negative
            if currentSlope(j-1) > 0 && currentSlope(j) <= 0
                zeroIdx = [zeroIdx, j];
            end
        end
        
        % If we found zero crossings, use the first one after onset
        if ~isempty(zeroIdx)
            % Find the first zero crossing
            firstZero = zeroIdx(1);
            
            % Look for maximum value within a few samples of the zero crossing
            % (2 samples should capture the true peak)
            searchStart = max(1, firstZero-2);
            searchEnd = min(length(currentSmoothed), firstZero+2);
            [peakValue, peakOffset] = max(currentSmoothed(searchStart:searchEnd));
            peakIdx = searchStart + peakOffset - 1;
        else
            % Fallback: if no zero crossing found, use older method
            [peakValue, peakOffset] = max(currentSmoothed(onsetIdx+(.5*samplesPerMs):min(onsetIdx+(10*samplesPerMs), length(currentSmoothed))));
            peakIdx = onsetIdx + peakOffset + floor(0.5*samplesPerMs) - 1;
        end
        %%
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
        
        % Add this event to the list of detected intervals
        detectedIntervals(end+1, :) = [onsetIdx, offsetIdx];
    end
    
    % Calculate frequency for this sweep (events per second)
    frequency = epscCount / effectiveDuration;
    
    fprintf('Detected %d EPSCs in sweep %d (after %.1f ms), Frequency: %.2f Hz\n', ...
        epscCount, sweepIdx, params.startTime, frequency);
    
    % Store results for this sweep
    allResults(sweepIdx).sweep = sweepIdx;
    allResults(sweepIdx).epscs = epscs;
    
    % Calculate sweep statistics
    if epscCount > 0
        allResults(sweepIdx).stats.meanAmplitude = mean([epscs.amplitude]);
        allResults(sweepIdx).stats.meanHalfWidth = mean([epscs.halfWidth]);
        % Calculate frequency based on effective duration (after start time)
        allResults(sweepIdx).stats.frequency = frequency; % in Hz
    else
        allResults(sweepIdx).stats = struct();
        allResults(sweepIdx).stats.meanAmplitude = NaN;
        allResults(sweepIdx).stats.meanHalfWidth = NaN;
        allResults(sweepIdx).stats.frequency = 0;
    end
    
    % Store start time information
    allResults(sweepIdx).analysisStartTime = params.startTime;
    allResults(sweepIdx).effectiveDuration = effectiveDuration;

    % Plot the first sweep that has events
if epscCount > 0 && ~firstPlotMade
    % Create a new figure for this sweep
    figure('Position', [100, 100, 1200, 800]);
    
    % Plot the smoothed signal
    ax1 = subplot(2,1,1);
    plot(timeVector, currentSmoothed, 'k');
    hold on;
    
    % Mark onsets, peaks, and offsets
    onsetTimes = timeVector([epscs.onset]);
    peakTimes = timeVector([epscs.peak]);
    offsetTimes = timeVector([epscs.offset]);
    
    plot(onsetTimes, currentSmoothed([epscs.onset]), 'go', 'MarkerSize', 5, 'LineWidth', 1.5);
    plot(peakTimes, currentSmoothed([epscs.peak]), 'ro', 'MarkerSize', 5, 'LineWidth', 1.5);
    plot(offsetTimes, currentSmoothed([epscs.offset]), 'bo', 'MarkerSize', 5, 'LineWidth', 1.5);
    
    % Add frequency information to the plot title
    title(sprintf('Detected EPSCs in Sweep %d (%d events, %.2f Hz)', sweepIdx, epscCount, frequency));
    
    % Set axis labels and style
    xlabel('Time (ms)');
    if params.inverted
        ylabel('Current (pA, inverted)');
    else
        ylabel('Current (pA)');
    end
    set(ax1, 'TickDir', params.plotStyle.tickDir, 'Box', params.plotStyle.box);
    grid(params.plotStyle.grid);
    
    % Plot the slope in the second subplot
    ax2 = subplot(2,1,2);
    plot(timeVector, currentSlope, 'b');
    hold on;
    
    % Mark EPSC onsets on slope trace
    plot(onsetTimes, currentSlope([epscs.onset]), 'go', 'MarkerSize', 5, 'LineWidth', 1.5);
    
    title('Current Slope (First Derivative)');
    xlabel('Time (ms)');
    ylabel('Slope (pA/ms)');
    set(ax2, 'TickDir', params.plotStyle.tickDir, 'Box', params.plotStyle.box);
    grid(params.plotStyle.grid);
    
    % Link the x-axes of the two plots
    linkaxes([ax1, ax2], 'x');
    sgtitle(filename);
    
    % Save the figure
    savefig(fullfile(outpath, [filename '_sweep' num2str(sweepIdx) '.fig']));
    
    % Add a flag to indicate that we've already made the first plot
    firstPlotMade = true;
end
end


%% Step 5: Display EPSC statistics and summaries
% Summary statistics across all sweeps
totalEvents = 0;
allAmplitudes = [];
allHalfWidths = [];
allAreas = [];
allTaus = [];
sweepsWithEvents = 0;
validSweeps = 0; % Count of sweeps with valid baseline

for s = 1:numSweeps
    % Check if this sweep had a stable baseline
    if isfield(allResults(s), 'baselineStable') && allResults(s).baselineStable
        validSweeps = validSweeps + 1;
        
        epscCount = length(allResults(s).epscs);
        if epscCount > 0
            totalEvents = totalEvents + epscCount;
            allAmplitudes = [allAmplitudes, [allResults(s).epscs.amplitude]];
            allHalfWidths = [allHalfWidths, [allResults(s).epscs.halfWidth]];
            allAreas = [allAreas, [allResults(s).epscs.area]];

            sweepsWithEvents = sweepsWithEvents + 1;
        end
    end
end

% Calculate overall frequency across all valid sweeps
totalDuration = validSweeps * effectiveDuration; % Total recording time in seconds
overallFrequency = totalEvents / totalDuration;

% Display overall summary
fprintf('\n===== EPSC ANALYSIS SUMMARY =====\n');
fprintf('Total sweeps: %d\n', numSweeps);
fprintf('Valid sweeps (stable baseline): %d\n', validSweeps);
fprintf('Sweeps with EPSCs: %d (%.1f%% of valid sweeps)\n', sweepsWithEvents, 100*sweepsWithEvents/max(1,validSweeps));
fprintf('Total EPSCs detected: %d\n', totalEvents);
fprintf('Average EPSCs per sweep: %.2f\n', totalEvents/max(1,validSweeps));
fprintf('Overall frequency: %.2f Hz\n', overallFrequency);

if totalEvents > 0
    % Calculate overall stats
    fprintf('\nEPSC Characteristics (mean � std):\n');
    fprintf('Amplitude: %.2f � %.2f pA\n', mean(allAmplitudes), std(allAmplitudes));
    fprintf('Half-width: %.2f � %.2f ms\n', mean(allHalfWidths), std(allHalfWidths));
    fprintf('Area (charge): %.2f � %.2f pAms\n', mean(allAreas), std(allAreas));

    % Plot histogram of EPSC amplitudes
    figure('Position', [100, 500, 800, 400]);

    subplot(1,2,1);
    histogram(allAmplitudes, min(20, ceil(sqrt(length(allAmplitudes)))));
    title('EPSC Amplitude Distribution');
    xlabel('Amplitude (pA)');
    ylabel('Count');
    set(gca, 'TickDir', params.plotStyle.tickDir, 'Box', params.plotStyle.box);
    grid(params.plotStyle.grid);

    subplot(1,2,2);
    histogram(allHalfWidths, min(20, ceil(sqrt(length(allHalfWidths)))));
    title('EPSC Half-Width Distribution');
    xlabel('Half-Width (ms)');
    ylabel('Count');
    set(gca, 'TickDir', params.plotStyle.tickDir, 'Box', params.plotStyle.box);
    grid(params.plotStyle.grid);
    print(fullfile(outpath,[filename '_histograms.pdf']),'-dpdf');

    % Create summary table for each sweep
    sweepSummary = table();
    sweepSummary.Sweep = (1:numSweeps)';

    for s = 1:numSweeps
        sweepSummary.EPSCCount(s) = length(allResults(s).epscs);
        
        % Set baselineStable column
        if isfield(allResults(s), 'baselineStable')
            sweepSummary.BaselineStable(s) = allResults(s).baselineStable;
        else
            sweepSummary.BaselineStable(s) = false;
        end

        if sweepSummary.EPSCCount(s) > 0
            sweepSummary.MeanAmplitude(s) = allResults(s).stats.meanAmplitude;
            sweepSummary.MeanHalfWidth(s) = allResults(s).stats.meanHalfWidth;
            sweepSummary.Frequency(s) = allResults(s).stats.frequency;
        else
            sweepSummary.MeanAmplitude(s) = NaN;
            sweepSummary.MeanHalfWidth(s) = NaN;
            if sweepSummary.BaselineStable(s)
                sweepSummary.Frequency(s) = 0;
            else
                sweepSummary.Frequency(s) = NaN;
            end
        end
    end

    % Display the summary table
    disp(sweepSummary);

    % Save results automatically 
    saveFile = fullfile(outpath, [filename '_epsc_results.mat']);
    save(saveFile, 'allResults', 'sweepSummary', 'params', 'overallFrequency');
    fprintf('Results saved to %s\n', saveFile);

    % Export event times to text file automatically
    exportFile = fullfile(outpath, [filename '_epsc_events.txt']);
    fid = fopen(exportFile, 'w');

    fprintf(fid, 'Sweep\tEvent\tOnset(ms)\tPeak(ms)\tAmplitude(pA)\tHalfWidth(ms)\n');

    for s = 1:numSweeps
        for e = 1:length(allResults(s).epscs)
            fprintf(fid, '%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
                s, e, ...
                timeVector(allResults(s).epscs(e).onset), ...
                timeVector(allResults(s).epscs(e).peak), ...
                allResults(s).epscs(e).amplitude, ...
                allResults(s).epscs(e).halfWidth);
        end
    end

    fclose(fid);
    fprintf('Event times exported to %s\n', exportFile);

    % Create and save a frequency summary file
    freqFile = fullfile(outpath, [filename '_frequency_summary.txt']);
    fid = fopen(freqFile, 'w');
    
    fprintf(fid, 'Sweep\tFrequency(Hz)\tEvent_Count\tBaseline_Stable\n');
    for s = 1:numSweeps
        if isfield(allResults(s), 'baselineStable')
            baselineStable = allResults(s).baselineStable;
        else
            baselineStable = false;
        end
        
        if baselineStable && isfield(allResults(s), 'stats') && isfield(allResults(s).stats, 'frequency')
            frequency = allResults(s).stats.frequency;
        else
            frequency = NaN;
        end
        
        eventCount = length(allResults(s).epscs);
        
        fprintf(fid, '%d\t%.2f\t%d\t%d\n', s, frequency, eventCount, baselineStable);
    end
    
    % Write overall summary
    fprintf(fid, '\n===== SUMMARY =====\n');
    fprintf(fid, 'Overall Frequency: %.2f Hz\n', overallFrequency);
    fprintf(fid, 'Total Events: %d\n', totalEvents);
    fprintf(fid, 'Valid Sweeps: %d\n', validSweeps);
    fprintf(fid, 'Sweeps With Events: %d\n', sweepsWithEvents);
    
    fclose(fid);
    fprintf('Frequency summary exported to %s\n', freqFile);

  
else
    disp('No EPSCs detected in any sweep.');
end

%% Create a frequency plot across sweeps
figure('Position', [100, 100, 900, 400]);

% Extract frequencies and event counts for plotting
sweepNumbers = 1:numSweeps;
frequencies = zeros(1, numSweeps);
eventCounts = zeros(1, numSweeps);
stableBaseline = false(1, numSweeps);

for s = 1:numSweeps
    if isfield(allResults(s), 'baselineStable')
        stableBaseline(s) = allResults(s).baselineStable;
    end
    
    if stableBaseline(s) && isfield(allResults(s), 'stats') && isfield(allResults(s).stats, 'frequency')
        frequencies(s) = allResults(s).stats.frequency;
        eventCounts(s) = length(allResults(s).epscs);
    else
        frequencies(s) = NaN;
        eventCounts(s) = 0;
    end
end

% Plot frequencies across sweeps
subplot(2,1,1);
bar(sweepNumbers, frequencies);
hold on;

% Add overall frequency reference line
if validSweeps > 0
    plot([0.5, numSweeps+0.5], [overallFrequency, overallFrequency], 'r--', 'LineWidth', 1.5);
    text(numSweeps/2, overallFrequency*1.1, sprintf('Overall: %.2f Hz', overallFrequency), ...
        'HorizontalAlignment', 'center', 'Color', 'r');
end

% Mark unstable baseline sweeps
for s = 1:numSweeps
    if ~stableBaseline(s)
        % Mark unstable sweeps with a light red bar
        bar(s, 0, 'FaceColor', [1 0.8 0.8], 'EdgeColor', [0.8 0 0]);
    end
end

title('EPSC Frequency by Sweep');
xlabel('Sweep Number');
ylabel('Frequency (Hz)');
set(gca, 'TickDir', params.plotStyle.tickDir, 'Box', params.plotStyle.box);
grid(params.plotStyle.grid);
xlim([0.5, numSweeps+0.5]);

% Plot event counts across sweeps
subplot(2,1,2);
bar(sweepNumbers, eventCounts);

% Mark unstable baseline sweeps
hold on;
for s = 1:numSweeps
    if ~stableBaseline(s)
        % Mark unstable sweeps with a light red bar
        bar(s, 0, 'FaceColor', [1 0.8 0.8], 'EdgeColor', [0.8 0 0]);
    end
end

title('EPSC Count by Sweep');
xlabel('Sweep Number');
ylabel('Event Count');
set(gca, 'TickDir', params.plotStyle.tickDir, 'Box', params.plotStyle.box);
grid(params.plotStyle.grid);
xlim([0.5, numSweeps+0.5]);


% Save the frequency plot
print(fullfile(outpath,[filename '_frequency_plot.pdf']),'-dpdf');

end