%% This script compares discrimination performance between measured uni-sensory vs. 
% measured multisensory granule cell responses to discrimination
% performance between measured uni-sensory vs. simuluated multi-sensory
% responses produced by linear combinations of uni-sensory, This analysis
% is presented in Figure 5.

% Load Multi sensory Data
cd /Users/gfield/Dropbox/Work/Documents/Projects/Granule-Cell-Project/data/
load combocelltrials_dart
control_all_cells = allcells;
clear allcells
cd ../code/

% get number of cells and number of conditions
num_cells = length(control_all_cells);
num_conds = length(control_all_cells{1});


%% View the number of trials per condition. 

trial_tracker = zeros(num_cells, num_conds);
for gc = 1:num_cells
    for s_cond = 1:num_conds
       trial_tracker(gc, s_cond) = size(control_all_cells{gc}{s_cond}, 2);
    end
end
figure(1); clf;
plot(trial_tracker)
sum_unimodal_trials = sum(trial_tracker(:,1:3), 2);
min_trials = min(trial_tracker(:));
sum_multimodal_trials = sum(trial_tracker(:,4:5), 2);
min_sum_unimodal_trial_num = min(sum(trial_tracker(:,1:3), 2));
min_sum_multimodal_trail_num = min(sum(trial_tracker(:,4:5), 2));
legend('1 kHz, 68 dB', '10 kHz, 68 dB', '10 psi', ...
        '1 kHz 68 dB + 10psi', '10 kHz 68 dB + 10psi',...
        'Location', 'northwest')
legend('boxoff')


%% Process calcium signals 
% this will filter the signals and threshold to reduce noise.
conditions_to_process = [1,2,3,4,5];
all_resps = process_calcium_signals(control_all_cells, conditions_to_process);

% get responses for integrating the calcium signals
temp_resps = all_resps.all_integrate_resps;

%% Sort Responses

% sort unimodal resps into one group
unimodal_resps = zeros(num_cells, min_sum_unimodal_trial_num);
for gc = 1:num_cells
    long_buffer_counter = 1;
    for u_cond = 1:3 % unimodal conditions        
        short_buffer_counter = 1;
        temp_buffer = temp_resps{gc}{u_cond}(:,1);
        buffer_length = length(temp_buffer);
        while long_buffer_counter <= min_sum_unimodal_trial_num && short_buffer_counter <= buffer_length
        
            unimodal_resps(gc, long_buffer_counter) = temp_buffer(short_buffer_counter);
            long_buffer_counter = long_buffer_counter + 1;
            short_buffer_counter = short_buffer_counter +1;
            
        end
    end
end

% sort multimodal resps into one group
multimodal_resps = zeros(num_cells, min_sum_multimodal_trail_num);
for gc = 1:num_cells
    long_buffer_counter = 1;
    for u_cond = 4:5 % multimodal conditions
        
        short_buffer_counter = 1;
        temp_buffer = temp_resps{gc}{u_cond}(:,1);
        buffer_length = length(temp_buffer);

        while long_buffer_counter <= min_sum_multimodal_trail_num && short_buffer_counter <= buffer_length
        
            multimodal_resps(gc, long_buffer_counter) = temp_buffer(short_buffer_counter);
            long_buffer_counter = long_buffer_counter + 1;
            short_buffer_counter = short_buffer_counter +1;
            
        end
    end
end


%% discriminate unimodal from multimodal responses

% compute distances between all responses
num_trials = min([size(unimodal_resps, 2), size(multimodal_resps, 2)]);
observations_dists = pdist([unimodal_resps(:,1:num_trials), multimodal_resps(:,1:num_trials)]'); 
dists_matrix = squareform(observations_dists);

% compute correct identification of unimodal and multimodal responses
correct_hits = 0;
for trial = 1:size(dists_matrix, 1)
    temp_row = dists_matrix(trial,:);
    max_val = max(temp_row);
    temp_row(temp_row == 0) = max_val;
    [min_val, min_ind] = min(temp_row);
    
    
    if trial <= 30 && min_ind <= 30
        correct_hits = correct_hits + 1;
    elseif trial > 30 && min_ind > 30
        correct_hits = correct_hits + 1;
    end
 
end    
p_correct = correct_hits ./ size(dists_matrix, 1);


%% discriminate unimodal from SIMULATED multimodal responses

rng(1) % seed randomnumber generator to porduce repeateable output.
num_sims = 200; % set the number of simulations
sim_p_correct = zeros(num_sims,1);
for sim_counter = 1:num_sims
    sim_resps = zeros(num_cells, num_trials);
    for trial = 1:num_trials

        % simulate a multimodal response by combining random 2 unimodal
        for gc = 1:num_cells
            % randomly select auditory condition 1 versus 2
            aud_cond = ceil(2 * rand(1));

            % simulate a population of made up responses    
            % figure out how many trials there were for this cel
            aud_num_trials = size(control_all_cells{gc}{aud_cond},2);
            %randomly select a trial from aud condition 1
            aud_trial_index = ceil(aud_num_trials * rand(1));

            % randomly select a trial from puff condition
            num_trials_puff = size(control_all_cells{gc}{3},2);
            puff_trial_index = ceil(num_trials_puff * rand(1));

            % add responses from rand aud 1 to rand puff 1
            sim_combo_resp1 = cell2mat(control_all_cells{gc}{aud_cond}(:,aud_trial_index)) +...
                              cell2mat(control_all_cells{gc}{3}(:,puff_trial_index));    

            sim_all_cells{gc}{1} = sim_combo_resp1; 
        end

        all_resps = process_auditory_signals(sim_all_cells, 1);
        for gc = 1:num_cells
            sim_resps(gc, trial) = all_resps.all_integrate_resps{gc}{1}(1);
        end
    end
    
    % compute distance of simulated response from other responses
    observations_dists = pdist([unimodal_resps(:,1:num_trials), sim_resps(:,1:num_trials)]'); 
    dists_matrix = squareform(observations_dists);

    correct_hits = 0;

    for trial = 1:size(dists_matrix, 1)
        temp_row = dists_matrix(trial,:);
        max_val = max(temp_row);
        temp_row(temp_row == 0) = max_val;
        [min_val, min_ind] = min(temp_row);

        if trial <= 30 && min_ind <= 30
            correct_hits = correct_hits + 1;
        elseif trial > 30 && min_ind > 30
            correct_hits = correct_hits + 1;
        end

    end

    sim_p_correct(sim_counter) = correct_hits ./ size(dists_matrix, 1);
end
 

histogram(sim_p_correct, [0:0.01:1])
hold on 
plot([p_correct p_correct], [0 20], 'k-', 'LineWidth', 4)
hold off

    
mean(sim_p_correct)
std(sim_p_correct) ./ sqrt(length(sim_p_correct))
z_score = (p_correct - mean(sim_p_correct)) ./ std(sim_p_correct)                  
    

normpdf(z_score, 0, 1)
