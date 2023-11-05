cd ~/Desktop/forGreg/
load tonecelltrials_dart

% scan through and look at some responses of individual cells
num_cells = length(allcells);
condition_choice = 1;
for gc = 1:num_cells
    test_gc = cell2mat(allcells{gc}{condition_choice});
    figure(1); clf;
    imagesc(test_gc')
    hold on
    plot([30 30], [1 25], 'w');
    hold off
    xticks([30 60 90 120])
    xticklabels({'1', '2', '3', '4'})
    xlabel('seconds')
    ylabel('trials')  
    pause
end

%% Scan through and look at trials across cells for particular conditions

condition_choice = 1;

[num_frames, num_trials] = size(allcells{1}{condition_choice});

temp_cell_resp = zeros(num_cells, num_trials, num_frames);

for  gc = 1:num_cells
    [temp_trials, temp_frames] = size(cell2mat(allcells{gc}{condition_choice})');
    temp_cell_resp(gc, 1:temp_trials, 1:temp_frames) = cell2mat(allcells{gc}{condition_choice})';
end

for trials = 1:num_trials
    imagesc(squeeze(temp_cell_resp(:,trials,:)));
    hold on
    plot([30 30], [1 729], 'w');
    hold off
    xticks([30 60 90 120])
    xticklabels({'1', '2', '3', '4'})
    xlabel('seconds')
    ylabel('cells')  
    pause
end

%% process the individual traces
num_cells = length(allcells);
condition_choice = 1;
smoothing_scale = 3;
num_sum = 3;    % number of frames over which to sum the initial response
stim_start = 28; % frame that the stimulus starts
stim_end = 65;  % frame that the stimulus ends

ca_filter = [zeros(1,smoothing_scale), ones(1,smoothing_scale), zeros(1,smoothing_scale)];
ca_filter = ca_filter ./ norm(ca_filter);


% preallocate
all_init_resps = cell(num_cells,1);
all_peak_resps = cell(num_cells,1);
all_integrate_resps = cell(num_cells,1);  

for gc = 1:num_cells
    
    %preallocate
    temp_all_init_resps = cell(1,6);
    temp_all_peak_resps = cell(1,6);
    temp_all_integrate_resps = cell(1,6);
    
    for stim_cond = 1:length(allcells{gc})
    
        % extract the responses for a single cell under a stim condition
        temp_gc = cell2mat(allcells{gc}{stim_cond});

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


trial_matrix = zeros(num_cells, length(allcells{1}));
for gc = 1:num_cells
    for stim_cond = 1:length(allcells{1})
        trial_matrix(gc, stim_cond) = size(allcells{gc}{stim_cond}, 2);
    end
end

plot(trial_matrix)

%% Inspect the responses

% 328 for control auditory and 319 for dart auditory
trunc_cells = 319; % only use cells for which there are more than 20 trials
num_trials = 22;
all_responses = zeros(trunc_cells, length(allcells{1}), num_trials, 2);
resps_to_decode = all_integrate_resps;

for stim_cond = 1:length(allcells{1})
    for gc = 1:trunc_cells
        temp_resps = resps_to_decode{gc}{stim_cond};
        all_responses(gc, stim_cond, 1:num_trials,:) = temp_resps(1:num_trials,:);
    end
end

amps_only = squeeze(all_responses(:,:,:,2));
amps_only = permute(amps_only, [2 , 3, 1]);
amps_only = reshape(amps_only, [length(allcells{1})*num_trials, trunc_cells]);
[eigenvecz, pc_weights, eigenvalz] = pca(amps_only);

plot(cumsum(eigenvalz ./ sum(eigenvalz)), 'ko')

% plot for each category
figure(1); clf;
plot(pc_weights(1:num_trials,1), pc_weights(1:num_trials,2), '*r')
hold on
plot(pc_weights(num_trials + 1:2*num_trials,1), pc_weights(num_trials + 1:2*num_trials,2), '*r')
plot(pc_weights(num_trials*2 + 1: 3*num_trials,1), pc_weights(num_trials*2 + 1: 3*num_trials,2), '*b')
plot(pc_weights(num_trials*3 + 1: 4*num_trials,1), pc_weights(num_trials*3 + 1: 4*num_trials,2), '*b')
plot(pc_weights(num_trials*4 + 1: 5*num_trials,1), pc_weights(num_trials*4 + 1: 5*num_trials,2), '*k')
plot(pc_weights(num_trials*5 + 1: 6*num_trials,1), pc_weights(num_trials*5 + 1: 6*num_trials,2), '*k')
hold off

Y = tsne(amps_only);
figure(2); clf;
plot(Y(1:num_trials,1), Y(1:num_trials,2), '*r')
hold on
plot(Y(num_trials + 1:2*num_trials,1), Y(num_trials + 1:2*num_trials,2), '*r')
plot(Y(num_trials*2 + 1: 3*num_trials,1), Y(num_trials*2 + 1: 3*num_trials,2), '*b')
plot(Y(num_trials*3 + 1: 4*num_trials,1), Y(num_trials*3 + 1: 4*num_trials,2), '*b')
plot(Y(num_trials*4 + 1: 5*num_trials,1), Y(num_trials*4 + 1: 5*num_trials,2), '*k')
plot(Y(num_trials*5 + 1: 6*num_trials,1), Y(num_trials*5 + 1: 6*num_trials,2), '*k')

title('tsne: DART conditions')
hold off


%% decode

observations_dists = pdist(pc_weights);
 
dists_matrix = squareform(observations_dists);
figure(3)
imagesc(dists_matrix)

condition_hits = 0;
frequency_hits = 0;
amplitude_hits = 0;

for trial = 1:num_trials * length(allcells{gc})
    temp_row = dists_matrix(trial,:);
    max_val = max(temp_row);
    temp_row(temp_row == 0) = max_val;
    [min_val, min_ind] = min(temp_row);
    
    
    trial_group = ceil(trial ./ num_trials);
    clossest_neighbor_group = ceil(min_ind ./ num_trials);
    
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
    
condition_hits = condition_hits ./ (num_trials * length(allcells{gc}))
frequency_hits = frequency_hits ./ (num_trials * length(allcells{gc}))
amplitude_hits = amplitude_hits ./ (num_trials * length(allcells{gc}))



%%
for gc = 1:trunc_cells   
    temp_cell = filtered_gc{gc};
    
    for trial = 1:size(temp_cell,2)
        
        windowed_resp = temp_cell(28:63,trial);
        resp_indices = find(windowed_resp > baseline_noise(gc));
        new_resp(trial) = sum(windowed_resp(resp_indices));        
    end    
    all_resps{gc} = new_resp;    
end


temp_val = zeros(6,num_cells);
for gc = 1:num_cells
    for stim_cond = 1:length(allcells{1})
        temp_val(stim_cond,gc) = size(allcells{gc}{stim_cond}, 2);
    end
end

plot(temp_val')
xlabel('cells')
ylabel('trials')
axis([0 800 0 30])


%% start playing with decoding
cd ~/Desktop/forGreg/
load puffcelltrials_control

Control_Condition = allcells;
clear allcells;

num_cells = length(Control_Condition);
num_stim = length(Control_Condition{1});
[num_frames, num_trials] = size(Control_Condition{1}{1});
all_responses = zeros(num_cells, num_stim, num_trials, num_frames);

for stim = 1:num_stim
    for  gc = 1:num_cells
        [temp_trials, temp_frames] = size(cell2mat(Control_Condition{gc}{stim})');
        all_responses(gc, stim, 1:temp_trials, 1:temp_frames) = cell2mat(Control_Condition{gc}{stim})';
    end
end

collapsed_responses = sum(all_responses(:,:,:,28:34), 4);
collapsed_responses = reshape(collapsed_responses, num_cells, []);
[eigenvecs, pc_weights, eigenvals] = pca(collapsed_responses');
var_explained_vals = eigenvals./sum(eigenvals);
figure(4)
semilogy(var_explained_vals, 'ko')
axis([0 20 0.001 1])
title('control data')

figure(1); clf;
plot(pc_weights(1:25,1), pc_weights(1:25,2), 'g*')
hold on
plot(pc_weights(26:50,1), pc_weights(26:50,2), 'r*')
plot(pc_weights(51:75,1), pc_weights(51:75,2), 'b*')
legend('10 psi', '15 psi', '20 psi')
hold off

figure(2); clf;
plot3(pc_weights(1:25,1), pc_weights(1:25,2), pc_weights(1:25,3), 'g*')
hold on
plot3(pc_weights(26:50,1), pc_weights(26:50,2), pc_weights(26:50,3), 'r*')
plot3(pc_weights(51:75,1), pc_weights(51:75,2), pc_weights(51:75,3), 'b*')
hold off

Y = tsne(collapsed_responses');
figure(3); clf;
plot(Y(1:25,1), Y(1:25,2), 'g*')
hold on
plot(Y(26:50,1), Y(26:50,2), 'r*')
plot(Y(51:75,1), Y(51:75,2), 'b*')
title('tsne: control conditions')
hold off

%%
cd ~/Desktop/forGreg/
load puffcelltrials_dart

DART_Condition = allcells;
clear allcells;

num_cells = length(DART_Condition);
num_stim = length(DART_Condition{1});
[num_frames, num_trials] = size(DART_Condition{1}{1});
all_responses = zeros(num_cells, num_stim, num_trials, num_frames);

for stim = 1:num_stim
    for  gc = 1:num_cells
        [temp_trials, temp_frames] = size(cell2mat(DART_Condition{gc}{stim})');
        all_responses(gc, stim, 1:temp_trials, 1:temp_frames) = cell2mat(DART_Condition{gc}{stim})';
    end
end

collapsed_responses = sum(all_responses(:,:,:,28:34), 4);
collapsed_responses = reshape(collapsed_responses, num_cells, []);
[eigenvecs, pc_weights, eigenvals] = pca(collapsed_responses');
var_explained_vals = eigenvals./sum(eigenvals);
figure(8)
semilogy(var_explained_vals, 'ko')
axis([0 20 0.001 1])
title('DART')

figure(5); clf;
plot(pc_weights(1:25,1), pc_weights(1:25,2), 'g*')
hold on
plot(pc_weights(26:50,1), pc_weights(26:50,2), 'r*')
plot(pc_weights(51:75,1), pc_weights(51:75,2), 'b*')
hold off

figure(6); clf;
plot3(pc_weights(1:25,1), pc_weights(1:25,2), pc_weights(1:25,3), 'g*')
hold on
plot3(pc_weights(26:50,1), pc_weights(26:50,2), pc_weights(26:50,3), 'r*')
plot3(pc_weights(51:75,1), pc_weights(51:75,2), pc_weights(51:75,3), 'b*')
hold off

Y = tsne(collapsed_responses');
figure(7); clf;
plot(Y(1:25,1), Y(1:25,2), 'g*')
hold on
plot(Y(26:50,1), Y(26:50,2), 'r*')
plot(Y(51:75,1), Y(51:75,2), 'b*')
title('tsne: DART conditions')
hold off



