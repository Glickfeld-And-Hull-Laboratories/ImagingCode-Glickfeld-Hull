
if ~exist('control_all_cells', 'var')
    cd /Users/gfield/Dropbox/Work/Documents/Projects/Granule-Cell-Project/data/
    load tonecelltrials_control
    control_all_cells = allcells;
end
cd ../code/
clear allcells


%% process the individual traces
num_cells = length(control_all_cells);
num_conds = length(control_all_cells{1});
smoothing_scale = 3;
num_sum = 3;    % number of frames over which to sum the initial response
stim_start = 28; % frame that the stimulus starts
stim_end = 65;  % frame that the stimulus ends
verbose = 1;
resp_flag = 3;  % 1 = amplitude; 2 = time; 3 = both
decode_flag = 1; % 1 = initial; 2 = peak; 3 = integrated

% 328 for control auditory and 319 for dart auditory
%control_trunc_cells = round(328/1); % only use cells for which there are more than 23 trials
%dart_trunc_cells = round(319/1); % only use cells for which there are more than 22 trials

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

%% some code for looking at how many cells/trials can be included
trial_matrix = zeros(num_cells, length(control_all_cells{1}));
for gc = 1:num_cells
    for stim_cond = 1:length(control_all_cells{1})
        trial_matrix(gc, stim_cond) = size(control_all_cells{gc}{stim_cond}, 2);
    end
end
figure(1); clf;
plot(trial_matrix)

%% Inspect the responses
bg_cell = 413;
end_cell = 729;
%bg_cell = 330;
%end_cell = 412;
%bg_cell = 1;
%end_cell = 320;
temp_num_cells = end_cell-bg_cell+1;
control_num_trials = 12;
dart_num_trials = 12;
all_responses = zeros(temp_num_cells, length(control_all_cells{1}), control_num_trials, 2);

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
        all_responses(gc-bg_cell+1, stim_cond, 1:control_num_trials,:) = temp_resps(1:control_num_trials,:);
    end
end

if resp_flag == 1
    distilled_resps = squeeze(all_responses(:,:,:,2));
    distilled_resps = permute(distilled_resps, [2 , 3, 1]);
    distilled_resps = reshape(distilled_resps, [length(control_all_cells{1})*control_num_trials, temp_num_cells]);
    [eigenvecz, pc_weights, eigenvalz] = pca(distilled_resps);
elseif resp_flag == 2
    distilled_resps = squeeze(all_responses(:,:,:,1));
    distilled_resps = permute(distilled_resps, [2 , 3, 1]);
    distilled_resps = reshape(distilled_resps, [length(control_all_cells{1})*control_num_trials, temp_num_cells]);
    [eigenvecz, pc_weights, eigenvalz] = pca(distilled_resps);
elseif resp_flag == 3
    distilled_resps = permute(all_responses, [2 , 3, 4, 1]);
    distilled_resps = reshape(distilled_resps, [length(control_all_cells{1})*control_num_trials, temp_num_cells*2]);
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
    plot(pc_weights(control_num_trials*3 + 1: 4*control_num_trials,1), pc_weights(control_num_trials*3 + 1: 4*control_num_trials,2), '*b')
    plot(pc_weights(control_num_trials*4 + 1: 5*control_num_trials,1), pc_weights(control_num_trials*4 + 1: 5*control_num_trials,2), '*k')
    plot(pc_weights(control_num_trials*5 + 1: 6*control_num_trials,1), pc_weights(control_num_trials*5 + 1: 6*control_num_trials,2), '*k')
    hold off
    title('PC plot')

end

%% decode

observations_dists = pdist(pc_weights); 
dists_matrix = squareform(observations_dists);
if verbose
    figure(4); clf; imagesc(dists_matrix)
end

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
        if ~isempty(intersect([ clossest_neighbor_group], [3 4]))
            frequency_hits = frequency_hits + 1;
        end
    elseif trial_group == 5 || trial_group == 5
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
%% DART condition
if ~exist('dart_all_cells', 'var')
    cd /Users/gfield/Dropbox/Work/Documents/Projects/Granule-Cell-Project/data/
    load tonecelltrials_dart
    dart_all_cells = allcells;
end
cd ../code/
clear allcells

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

%% some code for looking at how many cells/trials can be included
trial_matrix = zeros(num_cells, length(dart_all_cells{1}));
for gc = 1:num_cells
    for stim_cond = 1:length(dart_all_cells{1})
        trial_matrix(gc, stim_cond) = size(dart_all_cells{gc}{stim_cond}, 2);
    end
end
plot(trial_matrix)

%% Inspect the responses
all_responses = zeros(temp_num_cells, length(dart_all_cells{1}), dart_num_trials, 2);

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
        all_responses(gc-bg_cell+1, stim_cond, 1:dart_num_trials,:) = temp_resps(1:control_num_trials,:);
    end
end


if resp_flag == 1
    distilled_resps = squeeze(all_responses(:,:,:,2));
    distilled_resps = permute(distilled_resps, [2 , 3, 1]);
    distilled_resps = reshape(distilled_resps, [length(dart_all_cells{1})*dart_num_trials, temp_num_cells]);
    [eigenvecz, pc_weights, eigenvalz] = pca(distilled_resps);
elseif resp_flag == 2
    distilled_resps = squeeze(all_responses(:,:,:,1));
    distilled_resps = permute(distilled_resps, [2 , 3, 1]);
    distilled_resps = reshape(distilled_resps, [length(dart_all_cells{1})*dart_num_trials, temp_num_cells]);
    [eigenvecz, pc_weights, eigenvalz] = pca(distilled_resps);
elseif resp_flag == 3
    distilled_resps = permute(all_responses, [2 , 3, 4, 1]);
    distilled_resps = reshape(distilled_resps, [length(dart_all_cells{1})*dart_num_trials, temp_num_cells*2]);
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
    plot(pc_weights(dart_num_trials*3 + 1: 4*dart_num_trials,1), pc_weights(dart_num_trials*3 + 1: 4*dart_num_trials,2), '*b')
    plot(pc_weights(dart_num_trials*4 + 1: 5*dart_num_trials,1), pc_weights(dart_num_trials*4 + 1: 5*dart_num_trials,2), '*k')
    plot(pc_weights(dart_num_trials*5 + 1: 6*dart_num_trials,1), pc_weights(dart_num_trials*5 + 1: 6*dart_num_trials,2), '*k')
    hold off
    title('PC weights (DART)')

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
        if ~isempty(intersect([ clossest_neighbor_group], [3 4]))
            frequency_hits = frequency_hits + 1;
        end
    elseif trial_group == 5 || trial_group == 5
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


condition_performance = [control_condition_hits dart_condition_hits]
frequenc_performance = [control_frequency_hits dart_frequency_hits]
amplitude_performance = [control_amplitude_hits dart_amplitude_hits]

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
condition_performance_p_val = 1-normcdf(test_val)
xlabel(num2str(condition_performance_p_val))

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
frequency_performance_p_val = 1-normcdf(test_val)
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
amplitude_performance_p_val = 1-normcdf(abs(test_val))
xlabel(num2str(amplitude_performance_p_val));

%f = gcf;
%exportgraphics(f, 'classification-performance.pdf', 'ContentType', 'vector')

%% computing probability of chance performance

%normcdf(1/6, control_condition_hits, ci_c_condition)
%normcdf(1/3, control_frequency_hits, ci_c_frequency)
normcdf(1/2, control_amplitude_hits, ci_c_amplitude)


