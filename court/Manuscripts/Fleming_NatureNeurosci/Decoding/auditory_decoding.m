% This script load the granual cell data for the auditory condition. This
% Code was writted by GDF in Fall of 2022; 
% cleaned and commented by GDF 05-2023.

animal = 1; % which animal to use, 1, 2, 3, or 'all'
if ~strcmp(animal,'all')
    animal_folder = ['animal-',num2str(animal)];
else
    animal_folder = 'animal-all';
end

cell_down_sample_factor = 1; %This is the fraction of cells that will be used.

% verbose flags
verbose = 1; % plots stuff
print_verbose = 0; % export key figures

% This block defines parameters used for decoding. 
% resp_flag and decode_flag are two key flags that determine what features
% of the responses will be used for decoding. resp_flag defines the
% amplitude, the time, or both to be used for decoding. decode_flag defines
% the initial response, peak response or a combo of the two to be used for
% decoding.
resp_flag = 1;  % 1 = amplitude; 2 = time; 3 = both
decode_flag = 2; % 1 = initial; 2 = peak; 3 = integrated
shuffle_flag = 0; % shuffle trials for each cell.
non_zero_flag = 1; % only keep cells that have a non_zero response under control conditions;

% define parameters of the decoding
smoothing_scale = 3; % frames to smooth over
num_sum = 3;    % number of frames over which to sum the initial response

% range of frames overwhich to define a response.
stim_start = 28; % frame that the stimulus starts (some temporal slop here to -2)
stim_end = 65;  % frame that the stimulus ends + 5

%% Load Data
if ~exist('control_all_cells', 'var') || ~exist('dart_all_cells', 'var')
    cd /Users/gfield/Dropbox/Work/Documents/Projects/Granule-Cell-Project/data/ % needs to be edited by user to point to data
    load tonecelltrials_control
    control_all_cells = allcells;
    load tonecelltrials_dart
    dart_all_cells = allcells;
end
clear allcells

% got to the right directory for writing data (this will need to be edited
% by user for their directory structure)
cd /Users/gfield/Dropbox/Work/Documents/Projects/Granule-Cell-Project/code/
cd(animal_folder)

% get the number of conditions present 
num_conds = length(control_all_cells{1});

%% INSPECT the structure of the Data 

% get some initial stuff
control_num_cells = length(control_all_cells);
dart_num_cells = length(dart_all_cells);

% control condition
trial_matrix = zeros(control_num_cells, length(control_all_cells{1}));
for gc = 1:control_num_cells
    for stim_cond = 1:length(control_all_cells{1})
        trial_matrix(gc, stim_cond) = size(control_all_cells{gc}{stim_cond}, 2);
    end
end
figure(1); clf;
plot(trial_matrix)
legend('1 kHz, 68 dB', '1 kHz, 72 dB', '5 kHz, 68 dB', '5 kHz 72 dB',...
        '10 kHz, 68 dB', '10 kHz, 72 dB',...
        'location', 'northeast')
legend('boxoff')
title('CONTROL')
ylabel('trials')
xlabel('cell number')

% DART condition
trial_matrix = zeros(dart_num_cells, length(dart_all_cells{1}));
for gc = 1:dart_num_cells
    for stim_cond = 1:length(dart_all_cells{1})
        trial_matrix(gc, stim_cond) = size(dart_all_cells{gc}{stim_cond}, 2);
    end
end
figure(2); clf
plot(trial_matrix)
legend('1 kHz, 68 dB', '1 kHz, 72 dB', '5 kHz, 68 dB', '5 kHz 72 dB',...
        '10 kHz, 68 dB', '10 kHz, 72 dB',...
        'location', 'northeast')
legend('boxoff')
title('DART')
ylabel('trials')
xlabel('cell number')

%% Define cell ranges and trial numbers according to animals:

% Specify which group of cells to use
switch animal
    case 1
        % defines the range of cells and number of trials for animal 1
        bg_cell = 1;
        end_cell = 329;
        control_num_trials = 23; % minimum number of trials per condition, control and DART
        dart_num_trials = 22;

    case 2
        bg_cell = 331;
        end_cell = 412;
        control_num_trials = 11;
        dart_num_trials = 23;

    case 3
        bg_cell = 413;
        end_cell = 729;
        control_num_trials = 12;
        dart_num_trials = 16;

    case 'all'
        bg_cell = 1;
        end_cell = 729;
        control_num_trials = 11;
        dart_num_trials = 16;
end
decode_cell_num = end_cell - bg_cell + 1;

%% process the individual traces for control condition

% define filter
ca_filter = [zeros(1,smoothing_scale), ones(1,smoothing_scale), zeros(1,smoothing_scale)];
ca_filter = ca_filter ./ norm(ca_filter);

% preallocate cell arrays
all_init_resps = cell(control_num_cells,1);
all_peak_resps = cell(control_num_cells,1);
all_integrate_resps = cell(control_num_cells,1);  

for gc = 1:control_num_cells
    
    % preallocate
    temp_all_init_resps = cell(1,num_conds);
    temp_all_peak_resps = cell(1,num_conds);
    temp_all_integrate_resps = cell(1,num_conds);
    
    for stim_cond = 1:length(control_all_cells{gc})
    
        % extract the responses for a single cell under a stim condition
        temp_gc = cell2mat(control_all_cells{gc}{stim_cond});

        % preallocate     
        resp_vec1 = zeros(size(temp_gc, 2), 2);
        resp_vec2 = zeros(size(temp_gc, 2), 2);
        resp_vec3 = zeros(size(temp_gc, 2), 2);
        
        %put some edge buffering in 
        bb = repmat(temp_gc(1,:), [smoothing_scale, 1]);
        ee = repmat(temp_gc(end,:), [smoothing_scale, 1]);
        padded_resps = [bb; temp_gc; ee];
        
        % smooth data by taking boxcar average
        clear filtered_signals
        filtered_signals = zeros(size(temp_gc));
        for trial = 1:size(temp_gc, 2)
            temp_resp = conv(padded_resps(:,trial), ca_filter, 'same') ./ sum(ca_filter);        
            filtered_signals(:,trial) = temp_resp(smoothing_scale+1:end-smoothing_scale);        
        end          

        % after smoothing get mean and SD of baseline activity across
        % trials
        baseline_resp = filtered_signals(1:25,:);
        baseline_noise_threshold = 2*std(baseline_resp(:));

        thresholded_signals = filtered_signals;
        thresholded_signals(filtered_signals < baseline_noise_threshold) = 0;
        
        % extract some different possible response quantifications:
        for trial = 1:size(temp_gc, 2)
            
            % cut out the response in the time_window of interest
            temp_trial = thresholded_signals(stim_start:stim_end, trial);
            
         % 1.) time and amplitude of initial response
            temp_index = find(temp_trial > baseline_noise_threshold);
            if isempty(temp_index)
                resp_vec1(trial,:) = [0 0];
            else
                if (temp_index(1) + num_sum) > length(temp_trial)
                    temp_end = length(temp_trial);
                else
                    temp_end = temp_index(1)+num_sum;
                end
                resp_vec1(trial,:) = [temp_index(1), sum(temp_trial(temp_index(1):temp_end))];
            end
                
        % 2.) time and amplitude of peak response
            temp_index = find(temp_trial > baseline_noise_threshold);
            if isempty(temp_index)
                resp_vec2(trial,:) = [0 0];
            else
                [temp_max, max_index] = max(thresholded_signals(temp_index));
                resp_vec2(trial,:) = [max_index temp_max];
            end
           
        % 3.) integrate the calcium signal across all frames that exceed
           % the threshold
            temp_index = find(temp_trial > baseline_noise_threshold);
            if isempty(temp_index)
                resp_vec3(trial,:) = [0 0];
            else
                resp_vec3(trial,:) = [temp_index(1), sum(temp_trial(temp_index))];
            end            
        end  
        
        temp_all_init_resps{stim_cond} = resp_vec1;
        temp_all_peak_resps{stim_cond} = resp_vec2;
        temp_all_integrate_resps{stim_cond} = resp_vec3;
    end
    
    all_init_resps{gc} = temp_all_init_resps;
    all_peak_resps{gc} = temp_all_peak_resps;
    all_integrate_resps{gc} = temp_all_integrate_resps;  
    
end


%% Inspect the responses
all_responses = zeros(decode_cell_num, length(control_all_cells{1}), control_num_trials, 2);

if decode_flag == 1
    resps_to_decode = all_init_resps;
elseif decode_flag == 2
    resps_to_decode = all_peak_resps;
elseif decode_flag == 3
    resps_to_decode = all_integrate_resps;
else
    error('decode_flag is not a recognized value, must be 1,2,3')
end

for stim_cond = 1:length(control_all_cells{1})
    for gc = bg_cell:end_cell
        temp_resps = resps_to_decode{gc}{stim_cond};
        % Shuffle responses if required
        if shuffle_flag
            rand_ind = randperm(control_num_trials);
            all_responses(gc-bg_cell+1, stim_cond, 1:control_num_trials,:) = temp_resps(rand_ind,:);
        else
            all_responses(gc-bg_cell+1, stim_cond, 1:control_num_trials,:) = temp_resps(1:control_num_trials,:);
        end
    end
end

if resp_flag == 1
    distilled_resps = squeeze(all_responses(:,:,:,2));
    distilled_resps = permute(distilled_resps, [2 , 3, 1]);
    distilled_resps = reshape(distilled_resps, [length(control_all_cells{1})*control_num_trials, decode_cell_num]);
elseif resp_flag == 2
    distilled_resps = squeeze(all_responses(:,:,:,1));
    distilled_resps = permute(distilled_resps, [2 , 3, 1]);
    distilled_resps = reshape(distilled_resps, [length(control_all_cells{1})*control_num_trials, decode_cell_num]);
elseif resp_flag == 3
    distilled_resps = permute(all_responses, [2 , 3, 4, 1]);
    distilled_resps = reshape(distilled_resps, [length(control_all_cells{1})*control_num_trials, decode_cell_num*2]);
else
    error('resp_flag is not a recognized value, must be 1,2,3')
end

% only keep responses that are non_zero if flag is set to 1;
if non_zero_flag
    temp_check = zeros(size(distilled_resps,2), 1);
    for gc = 1:size(distilled_resps, 2)
        temp_check(gc) = any(distilled_resps(:,gc));
    end
    distilled_resps = distilled_resps(:,logical(temp_check));
end
[~, pc_weights, eigenvalz] = pca(distilled_resps); 


if verbose
    figure(3); clf;
    plot(cumsum(eigenvalz ./ sum(eigenvalz)), 'ko')
    xlabel('PC')
    ylabel('cumulative variance')

    % plot for each category
    figure(4); clf;
    pc_dim_1 = 1;
    pc_dim_2 = 2;
    plot(pc_weights(1:control_num_trials,pc_dim_1), pc_weights(1:control_num_trials,2), '^r')
    hold on
    plot(pc_weights(control_num_trials + 1:2*control_num_trials,pc_dim_1), pc_weights(control_num_trials + 1:2*control_num_trials,pc_dim_2), 'vk')
    plot(pc_weights(control_num_trials*2 + 1: 3*control_num_trials,pc_dim_1), pc_weights(control_num_trials*2 + 1: 3*control_num_trials,pc_dim_2), '+r')
    plot(pc_weights(control_num_trials*3 + 1: 4*control_num_trials,pc_dim_1), pc_weights(control_num_trials*3 + 1: 4*control_num_trials,pc_dim_2), 'xk')
    plot(pc_weights(control_num_trials*4 + 1: 5*control_num_trials,pc_dim_1), pc_weights(control_num_trials*4 + 1: 5*control_num_trials,pc_dim_2), 'or')
    plot(pc_weights(control_num_trials*5 + 1: 6*control_num_trials,pc_dim_1), pc_weights(control_num_trials*5 + 1: 6*control_num_trials,pc_dim_2), '*k')
    hold off
    title('Control PC weigths')
    xlabel('PC 1')
    ylabel('PC 2')
    legend('1 kHz, 68 dB', '1 kHz, 72 dB', '5 kHz, 68 dB', '5 kHz 72 dB',...
        '10 kHz, 68 dB', '10 kHz, 72 dB',...
        'location', 'northeast')
    legend('boxoff')

    if print_verbose
        f = gcf;
        exportgraphics(f, 'Control_scatter.pdf', 'ContentType', 'vector')
    end

    Y = tsne(distilled_resps);
    figure(5); clf;
    plot(Y(1:control_num_trials,1), Y(1:control_num_trials,2), '^r')
    hold on
    plot(Y(control_num_trials + 1:2*control_num_trials,1), Y(control_num_trials + 1:2*control_num_trials,2), 'vk')
    plot(Y(control_num_trials*2 + 1: 3*control_num_trials,1), Y(control_num_trials*2 + 1: 3*control_num_trials,2), '+r')
    plot(Y(control_num_trials*3 + 1: 4*control_num_trials,1), Y(control_num_trials*3 + 1: 4*control_num_trials,2), 'xk')
    plot(Y(control_num_trials*4 + 1: 5*control_num_trials,1), Y(control_num_trials*4 + 1: 5*control_num_trials,2), 'or')
    plot(Y(control_num_trials*5 + 1: 6*control_num_trials,1), Y(control_num_trials*5 + 1: 6*control_num_trials,2), '*k')
    title('tsne: Control conditions')
    hold off
    legend('1 kHz, 68 dB', '1 kHz, 72 dB', '5 kHz, 68 dB', '5 kHz 72 dB',...
        '10 kHz, 68 dB', '10 kHz, 72 dB',...
        'location', 'northeast')
    legend('boxoff')

     if print_verbose
        f = gcf;
        exportgraphics(f, 'Control t_sne.pdf', 'ContentType', 'vector')
    end
end

%% decode

observations_dists = pdist(pc_weights); 
dists_matrix = squareform(observations_dists);

condition_hits = 0;
frequency_hits = 0;
amplitude_hits = 0;

for trial = 1:control_num_trials * length(control_all_cells{gc})
    temp_row = dists_matrix(trial,:);
    max_val = max(temp_row);
    temp_row(temp_row == 0) = max_val;
    [min_val, min_ind] = min(temp_row);
    
    
    trial_group = ceil(trial ./ control_num_trials);
    clossest_neighbor_group = ceil(min_ind ./ control_num_trials);
    
    % accumulate hits for same condition
    if trial_group == clossest_neighbor_group
        condition_hits = condition_hits +1;
    end
    
    % accumulate hits for same frequency
    if trial_group == 1 || trial_group == 2
        if ~isempty(intersect(clossest_neighbor_group, [1 2]))
            frequency_hits = frequency_hits + 1;
        end
    elseif trial_group == 3 || trial_group == 4
        if ~isempty(intersect(clossest_neighbor_group, [3 4]))
            frequency_hits = frequency_hits + 1;
        end
    elseif trial_group == 5 || trial_group == 6
        if ~isempty(intersect(clossest_neighbor_group, [5 6]))
            frequency_hits = frequency_hits + 1;
        end
    end

    % accumulate hits for same amplitude
    if mod(trial_group, 2) == mod(clossest_neighbor_group,2)
        amplitude_hits = amplitude_hits + 1;
    end  
end
    
control_condition_hits = condition_hits ./ (control_num_trials * length(control_all_cells{gc}));
control_frequency_hits = frequency_hits ./ (control_num_trials * length(control_all_cells{gc}));
control_amplitude_hits = amplitude_hits ./ (control_num_trials * length(control_all_cells{gc}));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process DART Data

% process the individual traces
num_cells = length(dart_all_cells);

ca_filter = [zeros(1,smoothing_scale), ones(1,smoothing_scale), zeros(1,smoothing_scale)];
ca_filter = ca_filter ./ norm(ca_filter);

% preallocate
all_init_resps = cell(num_cells,1);
all_peak_resps = cell(num_cells,1);
all_integrate_resps = cell(num_cells,1);  

for gc = 1:num_cells
    
    %preallocate
    temp_all_init_resps = cell(1,num_conds);
    temp_all_peak_resps = cell(1,num_conds);
    temp_all_integrate_resps = cell(1,num_conds);
    
    for stim_cond = 1:length(dart_all_cells{gc})
    
        % extract the responses for a single cell under a stim condition
        temp_gc = cell2mat(dart_all_cells{gc}{stim_cond});

        % preallocate     
        resp_vec1 = zeros(size(temp_gc, 2), 2);
        resp_vec2 = zeros(size(temp_gc, 2), 2);
        resp_vec3 = zeros(size(temp_gc, 2), 2);
        
        %put some edge buffering in 
        bb = repmat(temp_gc(1,:), [smoothing_scale, 1]);
        ee = repmat(temp_gc(end,:), [smoothing_scale, 1]);
        padded_resps = [bb; temp_gc; ee];
        

        % smooth data by taking boxcar average
        clear filtered_signals
        filtered_signals = zeros(size(temp_gc));
        for trial = 1:size(temp_gc, 2)
            temp_resp = conv(padded_resps(:,trial), ca_filter, 'same') ./ sum(ca_filter);        
            filtered_signals(:,trial) = temp_resp(smoothing_scale+1:end-smoothing_scale);        
        end          

        % after smoothing get mean and SD of baseline activity across
        % trials
        baseline_resp = filtered_signals(1:25,:);
        baseline_noise_threshold = 2*std(baseline_resp(:));

        thresholded_signals = filtered_signals;
        thresholded_signals(filtered_signals < baseline_noise_threshold) = 0;
        
        % extract some different possible response quantifications:
        for trial = 1:size(temp_gc, 2)
            
            % cut out the response in the time_window of interest
            temp_trial = thresholded_signals(stim_start:stim_end, trial);
            
         % 1.) time and amplitude of initial response
            temp_index = find(temp_trial > baseline_noise_threshold);
            if isempty(temp_index)
                resp_vec1(trial,:) = [0 0];
            else
                if (temp_index(1) + num_sum) > length(temp_trial)
                    temp_end = length(temp_trial);
                else
                    temp_end = temp_index(1)+num_sum;
                end
                resp_vec1(trial,:) = [temp_index(1), sum(temp_trial(temp_index(1):temp_end))];
            end
                
        % 2.) time and amplitude of peak response
            temp_index = find(temp_trial > baseline_noise_threshold);
            if isempty(temp_index)
                resp_vec2(trial,:) = [0 0];
            else
                [temp_max, max_index] = max(thresholded_signals(temp_index));
                resp_vec2(trial,:) = [max_index temp_max];
            end
           
        % 3.) integrate the calcium signal across all frames that exceed
           % the threshold
            temp_index = find(temp_trial > baseline_noise_threshold);
            if isempty(temp_index)
                resp_vec3(trial,:) = [0 0];
            else
                resp_vec3(trial,:) = [temp_index(1), sum(temp_trial(temp_index))];
            end            
        end  
        
        temp_all_init_resps{stim_cond} = resp_vec1;
        temp_all_peak_resps{stim_cond} = resp_vec2;
        temp_all_integrate_resps{stim_cond} = resp_vec3;
    end
    
    all_init_resps{gc} = temp_all_init_resps;
    all_peak_resps{gc} = temp_all_peak_resps;
    all_integrate_resps{gc} = temp_all_integrate_resps;  
    
end

%% Inspect the responses
all_responses = zeros(decode_cell_num, length(dart_all_cells{1}), dart_num_trials, 2);

if decode_flag == 1
    resps_to_decode = all_init_resps;
elseif decode_flag == 2
    resps_to_decode = all_peak_resps;
elseif decode_flag == 3
    resps_to_decode = all_integrate_resps;
else
    error('that is not a recognized value for decode flag, must be 1,2,3')
end

for stim_cond = 1:length(dart_all_cells{1})
    for gc = bg_cell:end_cell
        temp_resps = resps_to_decode{gc}{stim_cond};
        % Shuffle responses if requires
        if shuffle_flag
            rand_ind = randperm(dart_num_trials);
            all_responses(gc-bg_cell+1, stim_cond, 1:dart_num_trials,:) = temp_resps(rand_ind,:);
        else
            all_responses(gc-bg_cell+1, stim_cond, 1:dart_num_trials,:) = temp_resps(1:dart_num_trials,:);
        end
    end
end

if resp_flag == 1
    distilled_resps = squeeze(all_responses(:,:,:,2));
    distilled_resps = permute(distilled_resps, [2 , 3, 1]);
    distilled_resps = reshape(distilled_resps, [length(dart_all_cells{1})*dart_num_trials, decode_cell_num]);
elseif resp_flag == 2
    distilled_resps = squeeze(all_responses(:,:,:,1));
    distilled_resps = permute(distilled_resps, [2 , 3, 1]);
    distilled_resps = reshape(distilled_resps, [length(dart_all_cells{1})*dart_num_trials, decode_cell_num]);
elseif resp_flag == 3
    distilled_resps = permute(all_responses, [2 , 3, 4, 1]);
    distilled_resps = reshape(distilled_resps, [length(dart_all_cells{1})*dart_num_trials, decode_cell_num*2]);
else
    error('resp_flag is not a recognized value, must be 1,2,3')
end
   
% removed cells that didn't response significantly during the control
% condition if flag is positive.
if non_zero_flag
    distilled_resps = distilled_resps(:,logical(temp_check));
end
[~, pc_weights, eigenvalz] = pca(distilled_resps); 


if verbose
    figure(6)
    plot(cumsum(eigenvalz ./ sum(eigenvalz)), 'ko')
    xlabel('PC')
    ylabel('cumulative variance')
    
    % plot for each category
    figure(7); clf;
    pc_dim_1 = 1;
    pc_dim_2 = 2;
    plot(pc_weights(1:dart_num_trials,1), pc_weights(1:dart_num_trials,2), '^r')
    hold on
    plot(pc_weights(dart_num_trials + 1:2*dart_num_trials,1), pc_weights(dart_num_trials + 1:2*dart_num_trials,2), 'vk')
    plot(pc_weights(dart_num_trials*2 + 1: 3*dart_num_trials,1), pc_weights(dart_num_trials*2 + 1: 3*dart_num_trials,2), '+r')
    plot(pc_weights(dart_num_trials*3 + 1: 4*dart_num_trials,1), pc_weights(dart_num_trials*3 + 1: 4*dart_num_trials,2), 'xk')
    plot(pc_weights(dart_num_trials*4 + 1: 5*dart_num_trials,1), pc_weights(dart_num_trials*4 + 1: 5*dart_num_trials,2), 'or')
    plot(pc_weights(dart_num_trials*5 + 1: 6*dart_num_trials,1), pc_weights(dart_num_trials*5 + 1: 6*dart_num_trials,2), '*k')
    hold off
    title('DART PC weights')
    xlabel('PC 1')
    ylabel('PC 2')
    legend('1 kHz, 68 dB', '1 kHz, 72 dB', '5 kHz, 68 dB', '5 kHz 72 dB',...
        '10 kHz, 68 dB', '10 kHz, 72 dB',...
        'location', 'northeast')
    legend('boxoff')

    if print_verbose
        f = gcf;
        exportgraphics(f, 'DART-scatter.pdf', 'ContentType', 'vector')
    end

    Y = tsne(distilled_resps);
    figure(8); clf;
    plot(Y(1:dart_num_trials,1), Y(1:dart_num_trials,2), '^r')
    hold on
    plot(Y(dart_num_trials + 1:2*dart_num_trials,1), Y(dart_num_trials + 1:2*dart_num_trials,2), 'vk')
    plot(Y(dart_num_trials*2 + 1: 3*dart_num_trials,1), Y(dart_num_trials*2 + 1: 3*dart_num_trials,2), '+r')
    plot(Y(dart_num_trials*3 + 1: 4*dart_num_trials,1), Y(dart_num_trials*3 + 1: 4*dart_num_trials,2), 'xk')
    plot(Y(dart_num_trials*4 + 1: 5*dart_num_trials,1), Y(dart_num_trials*4 + 1: 5*dart_num_trials,2), 'or')
    plot(Y(dart_num_trials*5 + 1: 6*dart_num_trials,1), Y(dart_num_trials*5 + 1: 6*dart_num_trials,2), '*k')
    title('tsne: DART conditions')
    hold off
    legend('1 kHz, 68 dB', '1 kHz, 72 dB', '5 kHz, 68 dB', '5 kHz 72 dB',...
        '10 kHz, 68 dB', '10 kHz, 72 dB',...
        'location', 'northeast')
    legend('boxoff')

    if print_verbose
        f = gcf;
        exportgraphics(f, 'DART-t_sne.pdf', 'ContentType', 'vector')
    end
end

%% decode

observations_dists = pdist(pc_weights); 
dists_matrix = squareform(observations_dists);

condition_hits = 0;
frequency_hits = 0;
amplitude_hits = 0;

for trial = 1:dart_num_trials * length(dart_all_cells{gc})
    temp_row = dists_matrix(trial,:);
    max_val = max(temp_row);
    temp_row(temp_row == 0) = max_val;
    [min_val, min_ind] = min(temp_row);
    
    
    trial_group = ceil(trial ./ dart_num_trials);
    clossest_neighbor_group = ceil(min_ind ./ dart_num_trials);
    
    % accumulate hits for same condition
    if trial_group == clossest_neighbor_group
        condition_hits = condition_hits +1;
    end
    
    % accumulate hits for same frequency
    if trial_group == 1 || trial_group == 2
        if ~isempty(intersect(clossest_neighbor_group, [1 2]))
            frequency_hits = frequency_hits + 1;
        end
    elseif trial_group == 3 || trial_group == 4
        if ~isempty(intersect(clossest_neighbor_group, [3 4]))
            frequency_hits = frequency_hits + 1;
        end
    elseif trial_group == 5 || trial_group == 6
        if ~isempty(intersect(clossest_neighbor_group, [5 6]))
            frequency_hits = frequency_hits + 1;
        end
    end

    % accumulate hits for same amplitude
    if mod(trial_group, 2) == mod(clossest_neighbor_group,2)
        amplitude_hits = amplitude_hits + 1;
    end  
end
    
dart_condition_hits = condition_hits ./ (dart_num_trials * length(dart_all_cells{gc}));
dart_frequency_hits = frequency_hits ./ (dart_num_trials * length(dart_all_cells{gc}));
dart_amplitude_hits = amplitude_hits ./ (dart_num_trials * length(dart_all_cells{gc}));

condition_performance = [control_condition_hits dart_condition_hits];
frequenc_performance = [control_frequency_hits dart_frequency_hits];
amplitude_performance = [control_amplitude_hits dart_amplitude_hits];

%% plot stuff and calculate p-values

control_trial_total = control_num_trials * num_conds;
dart_trial_total = dart_num_trials * num_conds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(10); clf;
subplot(1,3,1)
X = [1 2];
bar(X, condition_performance, 'FaceColor', [0.8 0.8 0.8])
axis([0 3 0 1])
xticklabels({'control', 'DART'})
hold on
plot([0 3], [1/6 1/6], 'k--')

%zeta = 1-(0.05/2); % confidence interval
ci_c_condition = sqrt(control_condition_hits*(1-control_condition_hits)/(control_num_trials * num_conds));
ci_d_condition = sqrt(dart_condition_hits*(1-dart_condition_hits)/(dart_num_trials * num_conds));

er = errorbar(X,condition_performance,[ci_c_condition, ci_d_condition]);   
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylabel('fraction correct')

p_hat = (control_trial_total * control_condition_hits + dart_trial_total *dart_condition_hits) ./ (control_trial_total+dart_trial_total);
test_val = (control_condition_hits - dart_condition_hits) ./ sqrt(p_hat*(1-p_hat)*(1/control_trial_total + 1/dart_trial_total));
condition_performance_p_val = 1-normcdf(test_val);
xlabel(num2str(condition_performance_p_val))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,2)
X = [1 2];
bar(X, frequenc_performance, 'FaceColor', [0.8 0.8 0.8])
axis([0 3 0 1])
xticklabels({'control', 'DART'})
hold on
plot([0 3], [1/3 1/3], 'k--')

%zeta = 1-(0.05/2); % confidence interval
ci_c_frequency = sqrt(control_frequency_hits*(1-control_frequency_hits)/(control_num_trials* num_conds));
ci_d_frequency = sqrt(dart_frequency_hits*(1-dart_frequency_hits)/(dart_num_trials * num_conds));

er = errorbar(X,frequenc_performance,[ci_c_frequency, ci_d_frequency]);   
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

p_hat = (control_trial_total * control_frequency_hits + dart_trial_total *dart_frequency_hits) ./ (control_trial_total+dart_trial_total);
test_val = (control_frequency_hits - dart_frequency_hits) ./ sqrt(p_hat*(1-p_hat)*(1/control_trial_total + 1/dart_trial_total));
frequency_performance_p_val = 1-normcdf(test_val);
xlabel(num2str(frequency_performance_p_val))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,3)
X = [1 2];
bar(X, amplitude_performance, 'FaceColor', [0.8 0.8 0.8])
axis([0 3 0 1])
xticklabels({'control', 'DART'})
hold on
plot([0 3], [1/2 1/2], 'k--')

%zeta = 1-(0.05/2); % confidence interval
ci_c_amplitude = sqrt(control_amplitude_hits*(1-control_amplitude_hits)/(control_num_trials * num_conds));
ci_d_amplitude = sqrt(dart_amplitude_hits*(1-dart_amplitude_hits)/(dart_num_trials * num_conds));

er = errorbar(X,amplitude_performance,[ci_c_amplitude, ci_d_amplitude]); 
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

p_hat = (control_trial_total * control_amplitude_hits + dart_trial_total *dart_amplitude_hits) ./ (control_trial_total+dart_trial_total);
test_val = (control_amplitude_hits - dart_amplitude_hits) ./ sqrt(p_hat*(1-p_hat)*(1/control_trial_total + 1/dart_trial_total));
amplitude_performance_p_val = 1-normcdf(test_val);
xlabel(num2str(amplitude_performance_p_val))

subplot(1,3,2)
title(animal_folder)
if print_verbose
    f = gcf;
    exportgraphics(f, 'classification-performance.pdf', 'ContentType', 'vector')
end

%% computing probability of chance performance
normcdf(1/6, control_condition_hits, ci_c_condition)
normcdf(1/3, control_frequency_hits, ci_c_frequency)
normcdf(1/2, control_amplitude_hits, ci_c_amplitude)


