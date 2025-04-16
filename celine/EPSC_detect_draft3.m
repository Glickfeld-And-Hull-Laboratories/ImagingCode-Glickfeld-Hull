%% EPSC Detection in Whole Cell Recordings (.abf format)
% This script detects excitatory postsynaptic currents (EPSCs) in 
% whole cell voltage clamp recordings from .abf files.
% Handles multiple sweeps (time  channel  sweep data structure).
% Modified to prevent overlapping events and start detection after 200 ms.

clear all; close all; clc;
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
params.conversionFactor = 1000;    % Some recordings had units other than 1pA 
params.plotStyle.tickDir = 'out'; % Direction of tick marks
params.plotStyle.grid = 'off';  % Grid off
params.plotStyle.box = 'off';   % Box off
params.startTime = 350;         % Start detection after this time (ms)
params.remove60Hz = true;       %option to remove 60Hs noise
%% Step 1: Import ABF file
% Using MATLAB's built-in abfload function
% If not available, please download from the MathWorks File Exchange
% or use the function provided at the end of this script

% Define the filename directly or use dialog if not provided
% To use directly, specify the full path to the file
fileName = '23719004'; % Enter your filename here, e.g., 'C:\Data\recording.abf'
filePath = [fileName,'.abf']
[data, si, h] = abfload(filePath);

% Extract data and convert sampling interval to ms
samplingInterval = si/1000; % Convert to ms
samplesPerMs = 1/samplingInterval;
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

%% Step 1.5: Remove 60 Hz noise
% Add this parameter to the parameters section:
% params.remove60Hz = true;       % Flag to apply 60 Hz notch filter

if params.remove60Hz
    fprintf('Applying 60 Hz notch filter to remove line noise...\n');
    
    % Calculate sampling frequency in Hz
    Fs = 1000 / samplingInterval; % Convert from ms to Hz
    
    % Design a notch filter
    % Calculate the normalized frequency (60 Hz / (Fs/2))
    w0 = 60/(Fs/2);
    bw = w0/35; % Bandwidth of the notch filter (adjust as needed)
    
    % Create notch filter
    [b, a] = iirnotch(w0, bw);
    
    % Create a second-order notch filter for the 120 Hz harmonic
    w1 = 120/(Fs/2);
    [b2, a2] = iirnotch(w1, bw);
    
    % Store original data for comparison
    originalData = data;
    
    % Apply filters to each sweep and channel
    fprintf('Applying notch filter to all sweeps...\n');
    for sweepIdx = 1:numSweeps
        for chanIdx = 1:numChannels
            % Apply 60 Hz filter
            data(:, chanIdx, sweepIdx) = filtfilt(b, a, data(:, chanIdx, sweepIdx));
            
            % Apply 120 Hz harmonic filter
            data(:, chanIdx, sweepIdx) = filtfilt(b2, a2, data(:, chanIdx, sweepIdx));
        end
    end
    
    % Plot comparison for first sweep to demonstrate filter effect
    figure('Position', [100, 500, 1200, 600]);
    
    % Plot in time domain
    subplot(2,2,1);
    plot(timeVector, originalData(:, params.channel, 1), 'b', 'LineWidth', 1);
    hold on;
    plot(timeVector, data(:, params.channel, 1), 'r', 'LineWidth', 1);
    title('Filter Effect on Raw Trace (Sweep 1)');
    xlabel('Time (ms)');
    ylabel('Current (pA)');
    legend('Original', 'Filtered');
    set(gca, 'TickDir', params.plotStyle.tickDir, 'Box', params.plotStyle.box);
    grid(params.plotStyle.grid);
    
    % Plot zoom of small segment
    subplot(2,2,2);
    segment = 1000:2000; % Adjust as needed based on your data
    plot(timeVector(segment), originalData(segment, params.channel, 1), 'b', 'LineWidth', 1);
    hold on;
    plot(timeVector(segment), data(segment, params.channel, 1), 'r', 'LineWidth', 1);
    title('Zoomed Segment');
    xlabel('Time (ms)');
    ylabel('Current (pA)');
    set(gca, 'TickDir', params.plotStyle.tickDir, 'Box', params.plotStyle.box);
    grid(params.plotStyle.grid);
    
    % Plot in frequency domain
    subplot(2,2,[3,4]);
    
    % Calculate FFT of original and filtered data
    L = length(timeVector);
    NFFT = 2^nextpow2(L);
    Y_orig = fft(originalData(:, params.channel, 1), NFFT)/L;
    Y_filt = fft(data(:, params.channel, 1), NFFT)/L;
    f = Fs/2*linspace(0,1,NFFT/2+1);
    
    % Plot single-sided amplitude spectrum
    plot(f, 2*abs(Y_orig(1:NFFT/2+1)), 'b');
    hold on;
    plot(f, 2*abs(Y_filt(1:NFFT/2+1)), 'r');
    title('Single-Sided Amplitude Spectrum');
    xlabel('Frequency (Hz)');
    ylabel('|Y(f)|');
    xlim([0, 200]); % Focus on frequencies up to 200 Hz
    legend('Original', 'Filtered');
    
    % Add markers for 60 Hz and 120 Hz
    line([60, 60], ylim, 'Color', 'k', 'LineStyle', '--');
    line([120, 120], ylim, 'Color', 'k', 'LineStyle', '--');
    text(60, max(2*abs(Y_orig(1:NFFT/2+1)))*0.8, '60 Hz', 'HorizontalAlignment', 'center');
    text(120, max(2*abs(Y_orig(1:NFFT/2+1)))*0.6, '120 Hz', 'HorizontalAlignment', 'center');
    
    set(gca, 'TickDir', params.plotStyle.tickDir, 'Box', params.plotStyle.box);
    grid(params.plotStyle.grid);
    
    fprintf('60 Hz noise filtering complete.\n');
end

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
        for i = onsetIdx+2:min(onsetIdx+(15*samplesPerMs), length(currentSlope)-1)
            % Check if slope crosses from positive to negative
            if currentSlope(i-1) > 0 && currentSlope(i) <= 0
                zeroIdx = [zeroIdx, i];
            end
        end
        
        % If we found zero crossings, use the first one after onset
        if ~isempty(zeroIdx)
            % Find the first zero crossing
            firstZero = zeroIdx(1);
            
            % Look for maximum value within a few samples of the zero crossing
            % (2 samples should capture the true peak)
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

%% Step 4: Plot the results
% Figure for event detection
figure('Position', [100, 100, 1200, 800]);

% Select the first sweep with events (or first sweep if none have events)
sweepToPlot = 20;
for i = 1:length(allResults)
    if ~isempty(allResults(i).epscs)
        sweepToPlot = i;
        break;
    end
end

% Get current data for selected sweep
currentData = data(:, params.channel, sweepToPlot);
if params.inverted
    currentData = -currentData; % Invert if looking at inward currents
    yLabel = 'Current (pA, inverted)';
else
    yLabel = 'Current (pA)';
end

% Apply preprocessing to match detection algorithm
currentFiltered = medfilt1(currentData, 5);
smoothWindow = ceil(params.smoothWindow / samplingInterval);
currentSmoothed = smooth(currentFiltered, smoothWindow);
currentSlope = diff(currentSmoothed)' / samplingInterval;
currentSlope = [0; currentSlope];

% Plot the current trace
ax1 = subplot(2,1,1);
plot(timeVector, currentSmoothed, 'k');
hold on;

% Get events for this sweep
epscs = allResults(sweepToPlot).epscs;
epscCount = length(epscs);

% Mark detected EPSCs
if epscCount > 0
    onsetTimes = timeVector([epscs.onset]);
    peakTimes = timeVector([epscs.peak]);
    offsetTimes = timeVector([epscs.offset]);
    
    % Mark onsets, peaks, and offsets
    plot(onsetTimes, currentSmoothed([epscs.onset]), 'go', 'MarkerSize', 5, 'LineWidth', 1.5);
    plot(peakTimes, currentSmoothed([epscs.peak]), 'ro', 'MarkerSize', 5, 'LineWidth', 1.5);
    plot(offsetTimes, currentSmoothed([epscs.offset]), 'bo', 'MarkerSize', 5, 'LineWidth', 1.5);
end

title(sprintf('Detected EPSCs in Sweep %d (%d events)', sweepToPlot, epscCount));
xlabel('Time (ms)');
ylabel(yLabel);
set(ax1, 'TickDir', params.plotStyle.tickDir, 'Box', params.plotStyle.box);
grid(params.plotStyle.grid);

% Plot the slope
ax2 = subplot(2,1,2);
plot(timeVector, currentSlope, 'b');
hold on;

% Mark EPSC onsets on slope trace
if epscCount > 0
    plot(onsetTimes, currentSlope([epscs.onset]), 'go', 'MarkerSize', 5, 'LineWidth', 1.5);
end

title('Current Slope (First Derivative)');
xlabel('Time (ms)');
ylabel('Slope (pA/ms)');
set(ax2, 'TickDir', params.plotStyle.tickDir, 'Box', params.plotStyle.box);
grid(params.plotStyle.grid);

% Link the x-axes of the two plots
linkaxes([ax1, ax2], 'x');
%%
% Highlight a random EPSC if any are detected
if epscCount > 0
    % Create a new figure for detailed view of a single EPSC
    figure('Position', [500, 100, 800, 400]);
    
    % Select a random EPSC to highlight
    epscIdx = randi(epscCount);
    
    % Extract event times
    eventOnset = epscs(epscIdx).onset;
    eventPeak = epscs(epscIdx).peak;
    eventOffset = epscs(epscIdx).offset;
    
    % Create window around the event (20ms)
    windowStart = max(1, eventOnset - round(20/samplingInterval));
    windowEnd = min(length(timeVector), eventOffset + round(20/samplingInterval));
    
    % Plot the event
    plot(timeVector(windowStart:windowEnd), currentSmoothed(windowStart:windowEnd), 'k', 'LineWidth', 1.5);
    hold on;
    
    % Mark key points
    plot(timeVector(eventOnset), currentSmoothed(eventOnset), 'go', 'MarkerSize', 8, 'LineWidth', 2);
    plot(timeVector(eventPeak), currentSmoothed(eventPeak), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    plot(timeVector(eventOffset), currentSmoothed(eventOffset), 'bo', 'MarkerSize', 8, 'LineWidth', 2);
    
    % Plot baseline
    if strcmp(params.baseline, 'pre')
        baselineWindow = max(1, eventOnset-30):eventOnset-5;
        baseline = mean(currentSmoothed(baselineWindow));
    else
        baseline = currentSmoothed(eventOnset);
    end
    
    line([timeVector(windowStart), timeVector(windowEnd)], [baseline, baseline], 'Color', 'r', 'LineStyle', '--');
    
    % Add half-width markers
    halfAmp = baseline + epscs(epscIdx).amplitude/2;
    line([timeVector(windowStart), timeVector(windowEnd)], [halfAmp, halfAmp], 'Color', 'g', 'LineStyle', ':');
    
    % Add labels
    title(sprintf('EPSC Event #%d from Sweep %d', epscIdx, sweepToPlot));
    xlabel('Time (ms)');
    ylabel(yLabel);
    legend('Current', 'Onset', 'Peak', 'Offset', 'Baseline', 'Half-amplitude');
    set(gca, 'TickDir', params.plotStyle.tickDir, 'Box', params.plotStyle.box);
    grid(params.plotStyle.grid);
    
    % Add event properties as text
    amplitude = epscs(epscIdx).amplitude;
    halfWidth = epscs(epscIdx).halfWidth;
    area = epscs(epscIdx).area;
    
    textStr = sprintf('Amplitude: %.1f pA\nHalf-width: %.2f ms\nArea: %.1f pAms', amplitude, halfWidth, area);
    if isfield(epscs, 'tau') && ~isnan(epscs(epscIdx).tau)
        textStr = [textStr sprintf('\nDecay : %.2f ms', epscs(epscIdx).tau)];
    end
    
    text(timeVector(windowStart) + 5, max(currentSmoothed(windowStart:windowEnd)), textStr, 'VerticalAlignment', 'top');
end


%% Step 5: Display EPSC statistics and summaries
% Summary statistics across all sweeps
totalEvents = 0;
allAmplitudes = [];
allHalfWidths = [];
allAreas = [];
allTaus = [];
sweepsWithEvents = 0;

for s = 1:numSweeps
    epscCount = length(allResults(s).epscs);
    if epscCount > 0
        totalEvents = totalEvents + epscCount;
        allAmplitudes = [allAmplitudes, [allResults(s).epscs.amplitude]];
        allHalfWidths = [allHalfWidths, [allResults(s).epscs.halfWidth]];
        allAreas = [allAreas, [allResults(s).epscs.area]];
        
        if isfield(allResults(s).epscs, 'tau')
            taus = [allResults(s).epscs.tau];
            allTaus = [allTaus, taus(~isnan(taus))];
        end
        
        sweepsWithEvents = sweepsWithEvents + 1;
    end
end

% Display overall summary
fprintf('\n===== EPSC ANALYSIS SUMMARY =====\n');
fprintf('Total sweeps: %d\n', numSweeps);
fprintf('Sweeps with EPSCs: %d (%.1f%%)\n', sweepsWithEvents, 100*sweepsWithEvents/numSweeps);
fprintf('Total EPSCs detected: %d\n', totalEvents);
fprintf('Average EPSCs per sweep: %.2f\n', totalEvents/numSweeps);

if totalEvents > 0
    % Calculate overall stats
    fprintf('\nEPSC Characteristics (mean  std):\n');
    fprintf('Amplitude: %.2f  %.2f pA\n', mean(allAmplitudes), std(allAmplitudes));
    fprintf('Half-width: %.2f  %.2f ms\n', mean(allHalfWidths), std(allHalfWidths));
    fprintf('Area (charge): %.2f  %.2f pAms\n', mean(allAreas), std(allAreas));
    
    if ~isempty(allTaus)
        fprintf('Decay time constant (tau): %.2f  %.2f ms\n', mean(allTaus), std(allTaus));
    end

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

    % Create summary table for each sweep
    sweepSummary = table();
    sweepSummary.Sweep = (1:numSweeps)';

    for s = 1:numSweeps
        sweepSummary.EPSCCount(s) = length(allResults(s).epscs);

        if sweepSummary.EPSCCount(s) > 0
            sweepSummary.MeanAmplitude(s) = allResults(s).stats.meanAmplitude;
            sweepSummary.MeanHalfWidth(s) = allResults(s).stats.meanHalfWidth;
            sweepSummary.Frequency(s) = allResults(s).stats.frequency;
        else
            sweepSummary.MeanAmplitude(s) = NaN;
            sweepSummary.MeanHalfWidth(s) = NaN;
            sweepSummary.Frequency(s) = 0;
        end
    end

    % Display sweep summary
    disp('Sweep-by-Sweep Summary:');
    disp(sweepSummary);

    % Plot event frequency across sweeps
    figure('Position', [100, 100, 900, 400]);

    subplot(1,2,1);
    bar(sweepSummary.Sweep, sweepSummary.EPSCCount);
    title('EPSC Count by Sweep');
    xlabel('Sweep Number');
    ylabel('Number of EPSCs');
    set(gca, 'TickDir', params.plotStyle.tickDir, 'Box', params.plotStyle.box);
    grid(params.plotStyle.grid);

    subplot(1,2,2);
    bar(sweepSummary.Sweep, sweepSummary.Frequency);
    title('EPSC Frequency by Sweep');
    xlabel('Sweep Number');
    ylabel('Frequency (Hz)');
    set(gca, 'TickDir', params.plotStyle.tickDir, 'Box', params.plotStyle.box);
    grid(params.plotStyle.grid);

    % Save results automatically rather than asking
    
    saveFile = [fileName '_epsc_results.mat'];
    % save(saveFile, 'allResults', 'sweepSummary', 'params');
    % fprintf('Results saved to %s\n', saveFile);
    % 
    % % Export event times to text file automatically
    % exportFile = [fileName '_epsc_events.txt'];
    % fid = fopen(exportFile, 'w');

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

    % Create a stack of aligned EPSCs automatically if not too many
    if totalEvents > 0 && totalEvents <= 500
        % Setup figure for stack plot
        figure('Position', [300, 200, 1000, 800]);

        % Define time window around peak for alignment (20ms)
        windowSize = round(20/samplingInterval);

        % Collect aligned traces
        alignedTraces = [];
        colorCodes = [];

        for s = 1:numSweeps
            for e = 1:length(allResults(s).epscs)
                peakIdx = allResults(s).epscs(e).peak;

                % Define window around peak
                startIdx = max(1, peakIdx - windowSize);
                endIdx = min(length(timeVector), peakIdx + windowSize);

                % Extract and align trace
                currentData = data(:, params.channel, s);
                if params.inverted
                    currentData = -currentData;
                end

                if endIdx - startIdx == 2*windowSize
                    trace = currentData(startIdx:endIdx);
                    alignedTraces = [alignedTraces; trace'];
                    colorCodes = [colorCodes; s];
                end
            end
        end

        % Create time vector centered at 0 (peak)
        alignedTime = (-windowSize:windowSize) * samplingInterval;

        % Plot aligned traces
        subplot(2,1,1);
        plot(alignedTime, alignedTraces', 'Color', [0.7 0.7 0.7]);
        hold on;
        plot(alignedTime, mean(alignedTraces), 'k', 'LineWidth', 2);

        title('Aligned EPSCs');
        xlabel('Time from Peak (ms)');
        ylabel(yLabel);
        set(gca, 'TickDir', params.plotStyle.tickDir, 'Box', params.plotStyle.box);
        grid(params.plotStyle.grid);

        % Plot heatmap of aligned traces
        subplot(2,1,2);
        imagesc(alignedTime, 1:size(alignedTraces,1), alignedTraces);
        colormap(jet);
        colorbar;

        title('EPSC Heatmap');
        xlabel('Time from Peak (ms)');
        ylabel('Event Number');
        set(gca, 'TickDir', params.plotStyle.tickDir, 'Box', params.plotStyle.box);
    end
else
    disp('No EPSCs detected in any sweep.');
end


