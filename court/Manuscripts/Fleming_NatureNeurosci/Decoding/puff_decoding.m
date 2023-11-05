
if ~exist('control_all_cells', 'var')
    cd /Users/gfield/Dropbox/Work/Documents/Projects/Granule-Cell-Project/data/
    load puffcelltrials_control
    control_all_cells = allcells;
end
clear allcells
cd ../code/

% get number of cells and number of conditions
num_cells = length(control_all_cells);
num_conds = length(control_all_cells{1});

trial_tracker = zeros(num_cells, num_conds);
for gc = 1:num_cells
    for s_cond = 1:num_conds
       trial_tracker(gc, s_cond) = size(control_all_cells{gc}{s_cond}, 2);
    end
end
if 1
    figure
    plot(trial_tracker)
end

%% process the individual traces
num_cells = length(control_all_cells);
smoothing_scale = 3;
num_sum = 1;    % number of frames over which to sum the initial response
stim_start = 27; % frame that the stimulus starts
stim_end = 37;  % frame that the stimulus ends
verbose = 1;
resp_flag = 2;  % 1 = amplitude; 2 = time; 3 = both
decode_flag = 2 % 1 = initial; 2 = peak; 3 = integrated

% 275 for control puff and 275 for dart puff
control_trunc_cells = round(275/1); % only use cells for which there are more than 15 trials
control_num_trials = round(15/1);
dart_trunc_cells = round(275/1); % only use cells for which there are more than 15 trials
dart_num_trials = round(15/1);

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
all_responses = zeros(control_trunc_cells, length(control_all_cells{1}), control_num_trials, 2);

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
    for gc = 1:control_trunc_cells
        temp_resps = resps_to_decode{gc}{stim_cond};
        all_responses(gc, stim_cond, 1:control_num_trials,:) = temp_resps(1:control_num_trials,:);
    end
end

if resp_flag == 1
    distilled_resps = squeeze(all_responses(:,:,:,2));
    distilled_resps = permute(distilled_resps, [2 , 3, 1]);
    distilled_resps = reshape(distilled_resps, [length(control_all_cells{1})*control_num_trials, control_trunc_cells]);
    [eigenvecz, pc_weights, eigenvalz] = pca(distilled_resps);
elseif resp_flag == 2
    distilled_resps = squeeze(all_responses(:,:,:,1));
    distilled_resps = permute(distilled_resps, [2 , 3, 1]);
    distilled_resps = reshape(distilled_resps, [length(control_all_cells{1})*control_num_trials, control_trunc_cells]);
    [eigenvecz, pc_weights, eigenvalz] = pca(distilled_resps);
elseif resp_flag == 3
    distilled_resps = permute(all_responses, [2 , 3, 4, 1]);
    distilled_resps = reshape(distilled_resps, [length(control_all_cells{1})*control_num_trials, control_trunc_cells*2]);
    [eigenvecz, pc_weights, eigenvalz] = pca(distilled_resps); 
else
    error('resp_flag is not a recognized value, must be 1,2,3')
end
    
    
if verbose
    figure(1); clf;
    plot(cumsum(eigenvalz ./ sum(eigenvalz)), 'ko')
    xlabel('PC')
    ylabel('cumulative variance')

    % plot for each category
    figure(2); clf;
    plot(pc_weights(1:control_num_trials,1), pc_weights(1:control_num_trials,2), '*r')
    hold on
    plot(pc_weights(control_num_trials + 1:2*control_num_trials,1), pc_weights(control_num_trials + 1:2*control_num_trials,2), '*r')
    plot(pc_weights(control_num_trials*2 + 1: 3*control_num_trials,1), pc_weights(control_num_trials*2 + 1: 3*control_num_trials,2), '*b')
    hold off
    title('PC plot')

    Y = tsne(distilled_resps);
    figure(3); clf;
    plot(Y(1:control_num_trials,1), Y(1:control_num_trials,2), '*r')
    hold on
    plot(Y(control_num_trials + 1:2*control_num_trials,1), Y(control_num_trials + 1:2*control_num_trials,2), '*r')
    plot(Y(control_num_trials*2 + 1: 3*control_num_trials,1), Y(control_num_trials*2 + 1: 3*control_num_trials,2), '*b')
    title('tsne: Control conditions')
    hold off
end

%% decode

observations_dists = pdist(pc_weights); 
dists_matrix = squareform(observations_dists);
if verbose
    figure(4); clf; imagesc(dists_matrix)
end

condition_hits = 0;

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
    
end
    
control_condition_hits = condition_hits ./ (control_num_trials * length(control_all_cells{gc}));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DART condition

if ~exist('dart_all_cells', 'var')
    cd /Users/gfield/Dropbox/Work/Documents/Projects/Granule-Cell-Project/data/
    load puffcelltrials_dart
    dart_all_cells = allcells;
end
clear allcells
cd ../code/

% get number of cells and number of conditions
num_cells = length(dart_all_cells);
num_conds = length(dart_all_cells{1});

trial_tracker = zeros(num_cells, num_conds);
for gc = 1:num_cells
    for s_cond = 1:num_conds
       trial_tracker(gc, s_cond) = size(dart_all_cells{gc}{s_cond}, 2);
    end
end
 if 1
     plot(trial_tracker)
 end

%% process the individual traces
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
all_responses = zeros(dart_trunc_cells, length(dart_all_cells{1}), dart_num_trials, 2);

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
    for gc = 1:dart_trunc_cells
        temp_resps = resps_to_decode{gc}{stim_cond};
        all_responses(gc, stim_cond, 1:dart_num_trials,:) = temp_resps(1:dart_num_trials,:);
    end
end

if resp_flag == 1
    distilled_resps = squeeze(all_responses(:,:,:,2));
    distilled_resps = permute(distilled_resps, [2 , 3, 1]);
    distilled_resps = reshape(distilled_resps, [length(dart_all_cells{1})*dart_num_trials, dart_trunc_cells]);
    [eigenvecz, pc_weights, eigenvalz] = pca(distilled_resps);
elseif resp_flag == 2
    distilled_resps = squeeze(all_responses(:,:,:,1));
    distilled_resps = permute(distilled_resps, [2 , 3, 1]);
    distilled_resps = reshape(distilled_resps, [length(dart_all_cells{1})*dart_num_trials, dart_trunc_cells]);
    [eigenvecz, pc_weights, eigenvalz] = pca(distilled_resps);
elseif resp_flag == 3
    distilled_resps = permute(all_responses, [2 , 3, 4, 1]);
    distilled_resps = reshape(distilled_resps, [length(dart_all_cells{1})*dart_num_trials, dart_trunc_cells*2]);
    [eigenvecz, pc_weights, eigenvalz] = pca(distilled_resps); 
else
    error('resp_flag is not a recognized value, must be 1,2,3')
end
   

if verbose
    figure(5)
    plot(cumsum(eigenvalz ./ sum(eigenvalz)), 'ko')
    xlabel('PC')
    ylabel('cumulative variance')
    
    % plot for each category
    figure(6); clf;
    plot(pc_weights(1:dart_num_trials,1), pc_weights(1:dart_num_trials,2), '*r')
    hold on
    plot(pc_weights(dart_num_trials + 1:2*dart_num_trials,1), pc_weights(dart_num_trials + 1:2*dart_num_trials,2), '*r')
    plot(pc_weights(dart_num_trials*2 + 1: 3*dart_num_trials,1), pc_weights(dart_num_trials*2 + 1: 3*dart_num_trials,2), '*b')
    hold off
    title('PC weights (DART)')

    Y = tsne(distilled_resps);
    figure(7); clf;
    plot(Y(1:dart_num_trials,1), Y(1:dart_num_trials,2), '*r')
    hold on
    plot(Y(dart_num_trials + 1:2*dart_num_trials,1), Y(dart_num_trials + 1:2*dart_num_trials,2), '*r')
    plot(Y(dart_num_trials*2 + 1: 3*dart_num_trials,1), Y(dart_num_trials*2 + 1: 3*dart_num_trials,2), '*b')
 
    title('tsne: DART conditions')
    hold off
end

%% decode

observations_dists = pdist(pc_weights); 
dists_matrix = squareform(observations_dists);
if verbose
    figure(8); clf; imagesc(dists_matrix)
end

condition_hits = 0;

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
    
end
    
dart_condition_hits = condition_hits ./ (dart_num_trials * length(dart_all_cells{gc}));

condition_performance = [control_condition_hits dart_condition_hits]

%% Generate plots

figure(10); clf;
%X = categorical({'control','DART'});
X = [1 2];
bar(X, condition_performance, 'FaceColor', [0.8 0.8 0.8])
axis([0 3 0 1])
xticklabels({'control', 'DART'})
hold on
plot([0 3], [1/3 1/3], 'k--')

%zeta = 1-(0.05/2); % confidence interval
ci_c = sqrt(control_condition_hits*(1-control_condition_hits)/control_num_trials);
ci_d = sqrt(dart_condition_hits*(1-dart_condition_hits)/dart_num_trials);

er = errorbar(X,condition_performance,[ci_c, ci_d]); % dividing by 2, to make it ~SEM    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('puff discrimination')
ylabel('fraction correct')


p_hat = (45 * control_condition_hits + 45 *dart_condition_hits) ./ 90;
test_val = (control_condition_hits - dart_condition_hits) ./ sqrt(p_hat*(1-p_hat)*(2/45));
p_val = normcdf(test_val)
