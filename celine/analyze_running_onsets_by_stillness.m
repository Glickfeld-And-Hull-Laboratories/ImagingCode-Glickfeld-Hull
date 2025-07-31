function [results] = analyze_running_onsets_by_stillness(day_ids, dataset_choice)
% ANALYZE_RUNNING_ONSETS_BY_STILLNESS - Comprehensive analysis of neural responses to running onsets
% categorized by preceding stillness duration
%
% Inputs:
%   day_ids - vector of experiment day IDs to analyze
%   dataset_choice - 0 for DART_expt_info, 1 for DART_V1_YM90K_Celine
%
% Outputs:
%   results - structure containing concatenated results across all experiments

%% Initialize
clear global;
close all;

% Setup dataset
switch dataset_choice
    case 0
        ds = 'DART_expt_info';
    case 1
        ds = 'DART_V1_YM90K_Celine';
end

rc = behavConstsDART;
eval(ds);

% Analysis parameters
stillLength_min = 2; % minimum stillness before onset
runLength = 2; % minimum running after onset
speed_threshold = 5.37; % cm/s threshold for running
onsetWin = 2; % seconds on either side of onset for timecourse
baseline_win = 0.5; % seconds before onset for baseline

% Stillness categorization bins (in seconds)
stillness_bins = [2, 5, 20, Inf]; % short: 2-5s, medium: 5-20s, long: 20+s
bin_names = {'short', 'medium', 'long'};
n_bins = length(bin_names);

% Initialize results storage
all_timecourses = {}; % cell array for each experiment
all_mean_responses = {};
all_correlations = {};
all_htp_status = {}; % track HTP+ vs HTP- cells
experiment_info = {};

fprintf('Starting analysis of %d experiments...\n', length(day_ids));

%% Main analysis loop
for exp_idx = 1:length(day_ids)
    day_id = day_ids(exp_idx);
    fprintf('\n=== Processing Experiment %d (day_id: %d) ===\n', exp_idx, day_id);
    
    try
        % Get experiment info
        pre_day = expt(day_id).multiday_matchdays;
        nd = 2;
        mouse = expt(day_id).mouse;
        experimentFolder = expt(day_id).exptType;
        
        % Setup DART string
        if expt(day_id).multiday_timesincedrug_hours > 0
            dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
        else
            dart_str = 'control';
        end
        
        % Load data
        fn_multi = fullfile(rc.achAnalysis, experimentFolder, mouse, ['multiday_' dart_str]);
        
        if ~exist(fn_multi, 'dir')
            fprintf('Warning: Directory not found for day %d, skipping...\n', day_id);
            continue;
        end
        
        % Load essential data
        load(fullfile(fn_multi, 'timecourses.mat'), 'cellTCs_match','red_ind_match');
        load(fullfile(fn_multi, 'input.mat'));
        load(fullfile(fn_multi, 'multiday_alignment.mat'), 'match_ind');
        
        inputStructure = input;
        clear input;
        
        frame_rate = double(inputStructure(1).frameImagingRateMs);
        allDays = [day_id, pre_day];
        
        % Initialize experiment-specific storage
        exp_timecourses = cell(1, nd);
        exp_mean_responses = cell(1, nd);
        exp_correlations = cell(1, nd);
        
        %% Process each day
        for day_idx = 1:nd
            fprintf('Processing day %d/%d...\n', day_idx, nd);
            
            % Calculate wheel speed
            wheel_speed_raw = wheelSpeedCalc(inputStructure(day_idx), 32, expt(allDays(1)).wheelColor);
            wheel_speed_clean = wheel_speed_raw;
            wheel_speed_clean(abs(wheel_speed_clean) < speed_threshold) = 0;
            wheel_speed_clean(wheel_speed_clean < 0) = 0; % Only forward running
            
            % Get neural data
            neural_data = cellTCs_match{day_idx};
            [nFrames, nCells] = size(neural_data);
            
            % Load stimulus timing for ITI calculation
            mouse_day = expt(allDays(day_idx)).mouse;
            date_day = expt(allDays(day_idx)).date;
            imgFolder = expt(allDays(day_idx)).contrastxori_runs{1};
            imgMatFile = [imgFolder '_000_000.mat'];
            dataPath = fullfile(rc.achData, mouse_day, date_day, imgFolder);
            
            if ~exist(fullfile(dataPath, imgMatFile), 'file')
                fprintf('Warning: Image data file not found for day %d, skipping...\n', day_idx);
                continue;
            end
            
            load(fullfile(dataPath, imgMatFile));
            [cStimOn, stimOffs] = photoFrameFinder_Sanworks(info.frame);
            
            % %% Calculate cross-correlations between wheel speed and neural activity
            % fprintf('Calculating cross-correlations...\n');
            % correlations = zeros(nCells, 1);
            % 
            % % Ensure same length
            % min_length = min(length(wheel_speed_clean), nFrames);
            % wheel_for_corr = wheel_speed_clean(1:min_length);
            % neural_for_corr = neural_data(1:min_length, :);
            % 
            % for cell_idx = 1:nCells
            %     try
            %         cell_trace = neural_for_corr(:, cell_idx);
            % 
            %         % Remove NaN and Inf values for correlation
            %         valid_idx = ~isnan(cell_trace) & ~isnan(wheel_for_corr) & ...
            %                    ~isinf(cell_trace) & ~isinf(wheel_for_corr);
            % 
            %         if sum(valid_idx) > 100 % require at least 100 valid points
            %             wheel_valid = wheel_for_corr(valid_idx);
            %             cell_valid = cell_trace(valid_idx);
            % 
            %             % Ensure both are column vectors
            %             wheel_valid = wheel_valid(:);
            %             cell_valid = cell_valid(:);
            % 
            %             % Check for sufficient variance (avoid constant signals)
            %             if std(wheel_valid) > 1e-10 && std(cell_valid) > 1e-10
            %                 correlations(cell_idx) = corr(wheel_valid, cell_valid);
            %             else
            %                 correlations(cell_idx) = NaN;
            %             end
            %         else
            %             correlations(cell_idx) = NaN;
            %         end
            %     catch ME
            %         fprintf('Warning: Correlation failed for cell %d: %s\n', cell_idx, ME.message);
            %         correlations(cell_idx) = NaN;
            %     end
            % end
            % 
            %% Find running onsets and categorize by stillness duration
            fprintf('Finding and categorizing running onsets...\n');
            
            % Create ITI mask
            ITI = ones(1, nFrames);
            for iTrial = 1:inputStructure(day_idx).trialSinceReset
                if cStimOn(iTrial) > 0 && cStimOn(iTrial) + inputStructure(day_idx).nScansOn <= nFrames
                    ITI(cStimOn(iTrial):cStimOn(iTrial) + inputStructure(day_idx).nScansOn) = 0;
                end
            end
            
            % Find potential onsets
            stillFrames = wheel_speed_clean < speed_threshold;
            runOnsets = find(diff(stillFrames) == -1) + 1;
            
            % Categorize onsets by stillness duration
            valid_onsets = [];
            stillness_durations = [];
            stillness_categories = [];
            
            for onset_idx = 1:length(runOnsets)
                onsetFrame = runOnsets(onset_idx);
                
                % Check if we have enough data around onset
                if onsetFrame <= onsetWin*frame_rate || ...
                   onsetFrame + onsetWin*frame_rate > nFrames || ...
                   onsetFrame > length(ITI)
                    continue;
                end
                
                % Check ITI requirement
                pre_onset_check = max(1, onsetFrame - frame_rate);
                post_onset_check = min(nFrames, onsetFrame + frame_rate);
                if ITI(pre_onset_check) ~= 1 || ITI(post_onset_check) ~= 1
                    continue;
                end
                
                % Calculate duration of stillness before onset
                stillness_start = onsetFrame - 1;
                while stillness_start > 1 && stillFrames(stillness_start)
                    stillness_start = stillness_start - 1;
                end
                stillness_duration = (onsetFrame - stillness_start - 1) / frame_rate;
                
                % Check minimum stillness requirement
                if stillness_duration < stillLength_min
                    continue;
                end
                
                % Check running requirement after onset
                testWinRun = onsetFrame:(onsetFrame + runLength*frame_rate);
                testWinRun = testWinRun(testWinRun <= length(stillFrames));
                runningReq = sum(~stillFrames(testWinRun)) >= (length(testWinRun) * 0.5);
                
                if ~runningReq
                    continue;
                end
                
                % Categorize by stillness duration
                bin_category = find(stillness_duration >= stillness_bins(1:end-1) & ...
                                  stillness_duration < stillness_bins(2:end));
                if isempty(bin_category)
                    bin_category = n_bins; % longest bin
                end
                
                valid_onsets(end+1) = onsetFrame;
                stillness_durations(end+1) = stillness_duration;
                stillness_categories(end+1) = bin_category;
            end
            
            fprintf('Found %d valid onsets\n', length(valid_onsets));
            
            %% Extract neural responses for each onset
            if isempty(valid_onsets)
                exp_timecourses{day_idx} = nan(n_bins, nCells, onsetWin*2*frame_rate);
                exp_mean_responses{day_idx} = nan(n_bins, nCells);
            else
                % Initialize storage for this day
                onset_timecourses = nan(length(valid_onsets), nCells, onsetWin*2*frame_rate);
                onset_mean_responses = nan(length(valid_onsets), nCells);
                
                for onset_idx = 1:length(valid_onsets)
                    onsetFrame = valid_onsets(onset_idx);
                    
                    % Define windows
                    full_window = (onsetFrame - onsetWin*frame_rate):(onsetFrame + onsetWin*frame_rate - 1);
                    baseline_window = (onsetFrame - baseline_win*frame_rate):(onsetFrame - 1);
                    response_window = onsetFrame:(onsetFrame + 2*frame_rate - 1); % 2 seconds after onset
                    
                    % Ensure windows are within bounds
                    full_window = full_window(full_window > 0 & full_window <= nFrames);
                    baseline_window = baseline_window(baseline_window > 0 & baseline_window <= nFrames);
                    response_window = response_window(response_window > 0 & response_window <= nFrames);
                    
                    if length(full_window) == onsetWin*2*frame_rate && ~isempty(baseline_window)
                        % Get neural data
                        neural_full = neural_data(full_window, :);
                        neural_baseline = mean(neural_data(baseline_window, :), 1);
                        neural_response = neural_data(response_window, :);
                        
                        % Calculate df/f
                        neural_dfof = (neural_full - neural_baseline) ./ neural_baseline;
                        response_dfof = (neural_response - neural_baseline) ./ neural_baseline;
                        
                        % Store
                        onset_timecourses(onset_idx, :, :) = neural_dfof';
                        onset_mean_responses(onset_idx, :) = mean(response_dfof, 1);
                    end
                end
                
                %% Average across onsets in each stillness bin
                day_timecourses = nan(n_bins, nCells, onsetWin*2*frame_rate);
                day_mean_responses = nan(n_bins, nCells);
                
                for bin_idx = 1:n_bins
                    bin_onsets = stillness_categories == bin_idx;
                    if sum(bin_onsets) > 0
                        day_timecourses(bin_idx, :, :) = squeeze(mean(onset_timecourses(bin_onsets, :, :), 1, 'omitnan'));
                        day_mean_responses(bin_idx, :) = mean(onset_mean_responses(bin_onsets, :), 1, 'omitnan');
                    end
                    fprintf('Bin %s: %d onsets\n', bin_names{bin_idx}, sum(bin_onsets));
                end
                
                exp_timecourses{day_idx} = day_timecourses;
                exp_mean_responses{day_idx} = day_mean_responses;
            end
            
            exp_correlations{day_idx} = correlations;
        end
        
        %% Store experiment results
        all_timecourses{exp_idx} = exp_timecourses;
        all_mean_responses{exp_idx} = exp_mean_responses;
        all_correlations{exp_idx} = exp_correlations;
        all_htp_status{exp_idx} = red_ind_match; % HTP+ indicator
        
        experiment_info{exp_idx} = struct(...
            'day_id', day_id, ...
            'mouse', mouse, ...
            'dart_str', dart_str, ...
            'nCells', nCells, ...
            'nHTP_positive', sum(red_ind_match), ...
            'frame_rate', frame_rate);
        
        fprintf('Experiment %d completed successfully\n', exp_idx);
        
    catch ME
        fprintf('Error processing experiment %d: %s\n', exp_idx, ME.message);
        continue;
    end
end

%% Concatenate results across experiments
fprintf('\n=== Concatenating results across experiments ===\n');

% Initialize concatenated matrices
concat_timecourses = [];
concat_mean_responses = [];
concat_correlations = [];
concat_htp_status = [];
experiment_labels = [];

total_cells = 0;
for exp_idx = 1:length(all_timecourses)
    if isempty(all_timecourses{exp_idx})
        continue;
    end
    
    exp_cells = size(all_timecourses{exp_idx}{1}, 2);
    total_cells = total_cells + exp_cells;
end

% Pre-allocate concatenated arrays
concat_timecourses = cell(1, nd);
concat_mean_responses = cell(1, nd);
concat_correlations = cell(1, nd);

for day_idx = 1:nd
    concat_timecourses{day_idx} = nan(n_bins, total_cells, onsetWin*2*frame_rate);
    concat_mean_responses{day_idx} = nan(n_bins, total_cells);
    concat_correlations{day_idx} = nan(total_cells, 1);
end

concat_htp_status = false(total_cells, 1);
experiment_labels = zeros(total_cells, 1);

% Fill concatenated arrays
cell_start_idx = 1;
for exp_idx = 1:length(all_timecourses)
    if isempty(all_timecourses{exp_idx})
        continue;
    end
    
    exp_cells = length(all_htp_status{exp_idx});
    cell_end_idx = cell_start_idx + exp_cells - 1;
    
    for day_idx = 1:nd
        if ~isempty(all_timecourses{exp_idx}{day_idx})
            concat_timecourses{day_idx}(:, cell_start_idx:cell_end_idx, :) = all_timecourses{exp_idx}{day_idx};
            concat_mean_responses{day_idx}(:, cell_start_idx:cell_end_idx) = all_mean_responses{exp_idx}{day_idx};
            concat_correlations{day_idx}(cell_start_idx:cell_end_idx) = all_correlations{exp_idx}{day_idx};
        end
    end
    
    concat_htp_status(cell_start_idx:cell_end_idx) = all_htp_status{exp_idx};
    experiment_labels(cell_start_idx:cell_end_idx) = exp_idx;
    
    cell_start_idx = cell_end_idx + 1;
end

%% Package results
results = struct();
results.timecourses = concat_timecourses;
results.mean_responses = concat_mean_responses;
results.correlations = concat_correlations;
results.htp_status = concat_htp_status;
results.experiment_labels = experiment_labels;
results.experiment_info = experiment_info;
results.analysis_params = struct(...
    'stillness_bins', stillness_bins, ...
    'bin_names', {bin_names}, ...
    'speed_threshold', speed_threshold, ...
    'onsetWin', onsetWin, ...
    'baseline_win', baseline_win, ...
    'frame_rate', frame_rate);

%% Summary statistics
fprintf('\n=== FINAL SUMMARY ===\n');
fprintf('Total experiments processed: %d\n', length(experiment_info));
fprintf('Total cells: %d\n', total_cells);
fprintf('HTP+ cells: %d (%.1f%%)\n', sum(concat_htp_status), 100*mean(concat_htp_status));
fprintf('HTP- cells: %d (%.1f%%)\n', sum(~concat_htp_status), 100*mean(~concat_htp_status));

for day_idx = 1:nd
    fprintf('\nDay %d:\n', day_idx);
    for bin_idx = 1:n_bins
        valid_cells = ~isnan(concat_mean_responses{day_idx}(bin_idx, :));
        fprintf('  %s stillness: %d cells with valid responses\n', ...
            bin_names{bin_idx}, sum(valid_cells));
    end
end

end

%% Helper function to run the analysis
function run_analysis_example()
    % Example usage:
    day_ids = [1, 2, 3, 4, 5]; % Replace with your actual day IDs
    dataset_choice = 0; % 0 for DART_expt_info, 1 for DART_V1_YM90K_Celine
    
    results = analyze_running_onsets_by_stillness(day_ids, dataset_choice);
    
    % Save results
    save('running_onset_analysis_results.mat', 'results');
    
    % Example of how to access results:
    % results.timecourses{day_idx}(bin_idx, cell_idx, time_point)
    % results.mean_responses{day_idx}(bin_idx, cell_idx)
    % results.correlations{day_idx}(cell_idx)
    % results.htp_status(cell_idx) % true for HTP+, false for HTP-
end