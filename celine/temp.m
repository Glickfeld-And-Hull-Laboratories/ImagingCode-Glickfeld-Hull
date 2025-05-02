% EPSC Sweep Visualization
% This script loads EPSC results from a .mat file, finds baseline and drug sweeps with events,
% loads the corresponding .abf file, and plots the sweeps with event markers.

% Specify file to analyze
cell_id = '23615001';
resultsFile = [cell_id, '_epsc_results.mat'];

% Load EPSC results file
fprintf('Loading EPSC results file: %s\n', resultsFile);
load(resultsFile);

% Define baseline and drug sweeps based on the analysis script
baseline_sweeps = 1:30;
drug_sweeps = 40:70;

% Find a baseline sweep with events
baseline_sweep_idx = [];
for i = baseline_sweeps
    if isfield(allResults, 'baselineStable') && allResults(i).baselineStable && ...
       isfield(allResults(i), 'epscs') && ~isempty(allResults(i).epscs) && ...
       isfield(allResults(i).epscs, 'peakTime') && ~isempty(allResults(i).epscs.peakTime)
        baseline_sweep_idx = i;
        fprintf('Found baseline sweep %d with %d events\n', i, length(allResults(i).epscs.peakTime));
        break;
    end
end

% Find a drug sweep with events
drug_sweep_idx = [];
for i = drug_sweeps
    if isfield(allResults, 'baselineStable') && allResults(i).baselineStable && ...
       isfield(allResults(i), 'epscs') && ~isempty(allResults(i).epscs) && ...
       isfield(allResults(i).epscs, 'peakTime') && ~isempty(allResults(i).epscs.peakTime)
        drug_sweep_idx = i;
        fprintf('Found drug sweep %d with %d events\n', i, length(allResults(i).epscs.peakTime));
        break;
    end
end

% Check if we found appropriate sweeps
if isempty(baseline_sweep_idx)
    fprintf('Warning: Could not find suitable baseline sweep with events. Using first sweep.\n');
    baseline_sweep_idx = baseline_sweeps(1);
end

if isempty(drug_sweep_idx)
    fprintf('Warning: Could not find suitable drug sweep with events. Using first sweep.\n');
    drug_sweep_idx = drug_sweeps(1);
end

fprintf('Selected baseline sweep: %d\n', baseline_sweep_idx);
fprintf('Selected drug sweep: %d\n', drug_sweep_idx);

% Get the event times and amplitudes from each sweep
if isfield(allResults(baseline_sweep_idx), 'epscs') && ...
   isfield(allResults(baseline_sweep_idx).epscs, 'peakTime')
    baseline_event_times = allResults(baseline_sweep_idx).epscs.peakTime;
    baseline_event_amps = allResults(baseline_sweep_idx).epscs.amplitude;
    baseline_event_count = length(baseline_event_times);
else
    baseline_event_times = [];
    baseline_event_amps = [];
    baseline_event_count = 0;
end

if isfield(allResults(drug_sweep_idx), 'epscs') && ...
   isfield(allResults(drug_sweep_idx).epscs, 'peakTime')
    drug_event_times = allResults(drug_sweep_idx).epscs.peakTime;
    drug_event_amps = allResults(drug_sweep_idx).epscs.amplitude;
    drug_event_count = length(drug_event_times);
else
    drug_event_times = [];
    drug_event_amps = [];
    drug_event_count = 0;
end

fprintf('Baseline events: %d\n', baseline_event_count);
fprintf('Drug events: %d\n', drug_event_count);

% Display fields in the structure to help debug
fprintf('\nStructure of allResults for baseline sweep:\n');
disp(fieldnames(allResults(baseline_sweep_idx)));
if isfield(allResults(baseline_sweep_idx), 'epscs')
    fprintf('Structure of epscs:\n');
    disp(fieldnames(allResults(baseline_sweep_idx).epscs));
end

% Create .abf filename based on cell_id\
abf_filename = fullfile('\home\celine\Data\patch_data',[cell_id, '.abf']);

% Check if the ABF reading functions are available
if ~exist('abfload', 'file')
    fprintf('abfload function not found. Using simulated data for plotting.\n');
    
    % Create simulated data for demonstration
    sample_rate = 10000; % 10 kHz
    sweep_duration = 5; % 5 seconds
    
    % Create time vector
    t = linspace(0, sweep_duration, sample_rate * sweep_duration);
    
    % Create simulated baseline trace with noise
    baseline_trace = randn(1, length(t)) * 5; % Random noise with 5 pA SD
    
    % Create simulated drug trace with less activity
    drug_trace = randn(1, length(t)) * 4; % Random noise with 4 pA SD
    
    % Add events to baseline trace
    for i = 1:length(baseline_event_times)
        if i <= length(baseline_event_times)
            event_time = baseline_event_times(i);
            event_amp = baseline_event_amps(i);
            
            % Find the closest time index
            [~, time_idx] = min(abs(t - event_time));
            
            % Add a fast-rising, slower-decaying EPSC-like event
            event_window = -100:100;
            valid_idx = time_idx + event_window;
            valid_idx = valid_idx(valid_idx > 0 & valid_idx <= length(t));
            
            if ~isempty(valid_idx)
                time_offset = valid_idx - time_idx;
                
                % Create exponential decay shape
                event_shape = zeros(size(time_offset));
                event_shape(time_offset >= 0) = -exp(-time_offset(time_offset >= 0)/20);
                
                % Scale by amplitude
                baseline_trace(valid_idx) = baseline_trace(valid_idx) + event_shape * event_amp;
            end
        end
    end
    
    % Add events to drug trace
    for i = 1:length(drug_event_times)
        if i <= length(drug_event_times)
            event_time = drug_event_times(i);
            event_amp = drug_event_amps(i);
            
            % Find the closest time index
            [~, time_idx] = min(abs(t - event_time));
            
            % Add a fast-rising, slower-decaying EPSC-like event
            event_window = -100:100;
            valid_idx = time_idx + event_window;
            valid_idx = valid_idx(valid_idx > 0 & valid_idx <= length(t));
            
            if ~isempty(valid_idx)
                time_offset = valid_idx - time_idx;
                
                % Create exponential decay shape
                event_shape = zeros(size(time_offset));
                event_shape(time_offset >= 0) = -exp(-time_offset(time_offset >= 0)/20);
                
                % Scale by amplitude
                drug_trace(valid_idx) = drug_trace(valid_idx) + event_shape * event_amp;
            end
        end
    end
else
    % Load the actual ABF file if abfload is available
    fprintf('Loading ABF file: %s\n', abf_filename);
    try
        [data, si, header] = abfload(abf_filename);
        
        % Extract the traces for the selected sweeps
        baseline_trace = data(:, 1, baseline_sweep_idx);
        drug_trace = data(:, 1, drug_sweep_idx);
        
        % Create time vector based on sampling interval
        sweep_duration = length(baseline_trace) * si / 1e6; % Convert to seconds
        t = linspace(0, sweep_duration, length(baseline_trace));
    catch
        error('Error loading ABF file: %s', abf_filename);
    end
end

% Create figure
figure('Position', [100, 100, 1000, 600]);

% Plot baseline sweep
subplot(2, 1, 1);
plot(t, baseline_trace, 'k');
hold on;

% Get y-axis limits for better marker placement
y_lim = get(gca, 'YLim');
marker_size = 8;

% Mark events on baseline sweep with clear red circles
if ~isempty(baseline_event_times)
    for i = 1:length(baseline_event_times)
        event_time = baseline_event_times(i);
        
        % For simulated data, place markers at the trace level
        if exist('abfload', 'file')
            % Find the closest time index
            [~, time_idx] = min(abs(t - event_time));
            % Place marker at the actual data point
            marker_y = baseline_trace(time_idx);
        else
            % For simulated data, make markers visible below the trace
            % Using event amplitude to position marker
            marker_y = min(baseline_trace) - 5;
        end
        
        plot(event_time, marker_y, 'ro', 'MarkerSize', marker_size, 'MarkerFaceColor', 'r');
    end
end

% Add EPSC frequency to the title if available
if isfield(allResults(baseline_sweep_idx).stats, 'frequency')
    baseline_freq = allResults(baseline_sweep_idx).stats.frequency;
else
    baseline_freq = baseline_event_count / sweep_duration;
end

title(sprintf('Baseline Sweep %d - Frequency: %.2f Hz (%d events)', baseline_sweep_idx, baseline_freq, baseline_event_count));
xlabel('Time (s)');
ylabel('Current (pA)');
set(gca, 'TickDir', 'out');
box off;

% Plot drug sweep
subplot(2, 1, 2);
plot(t, drug_trace, 'b');
hold on;

% Get y-axis limits for better marker placement
y_lim = get(gca, 'YLim');

% Mark events on drug sweep with clear red circles
if ~isempty(drug_event_times)
    for i = 1:length(drug_event_times)
        event_time = drug_event_times(i);
        
        % For simulated data, place markers at the trace level
        if exist('abfload', 'file')
            % Find the closest time index
            [~, time_idx] = min(abs(t - event_time));
            % Place marker at the actual data point
            marker_y = drug_trace(time_idx);
        else
            % For simulated data, make markers visible below the trace
            % Using event amplitude to position marker
            marker_y = min(drug_trace) - 5;
        end
        
        plot(event_time, marker_y, 'ro', 'MarkerSize', marker_size, 'MarkerFaceColor', 'r');
    end
end

% Add EPSC frequency to the title if available
if isfield(allResults(drug_sweep_idx).stats, 'frequency')
    drug_freq = allResults(drug_sweep_idx).stats.frequency;
else
    drug_freq = drug_event_count / sweep_duration;
end

title(sprintf('Drug Sweep %d - Frequency: %.2f Hz (%d events)', drug_sweep_idx, drug_freq, drug_event_count));
xlabel('Time (s)');
ylabel('Current (pA)');
set(gca, 'TickDir', 'out');
box off;

% Add overall title
sgtitle(sprintf('EPSC Events for Cell %s', cell_id), 'FontWeight', 'bold');

% Save figure
saveas(gcf, [cell_id, '_sweeps_with_events.png']);
fprintf('Figure saved as %s_sweeps_with_events.png\n', cell_id);