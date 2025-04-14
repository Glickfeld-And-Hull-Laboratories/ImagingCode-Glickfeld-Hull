%% Model for contrast by drug condition in the stationary state
% instructions to Claude: 
% I have the outcome variable pref_responses_stat_concat that is a two cell
% array. Each cell corresponds to a different drug condition, where
% cell 1 is the "DART" condition and cell 2 is the "control" condition. 
% Within each cell is a n neuron by 3 contrast matrix of visual responses
% at different visual % contrasts (25%, 50%, and 100%). I have two types of
% cells. the logical % array red_concat has 1's for SST cells and 0's for
% pyramidal cells.mouseID is a n neuron by 1 string array with the identifier of the mouse each neurons came from.
% Seperately for each cell type, I want to run a mixed
% effects model for neural response within each cell type, with drug condition and 
% contrast as the fixed effects and neuron and mouse as random effects.

% Separate data by cell type
sst_indices = find(red_concat == 1); % SST cells
pyr_indices = find(red_concat == 0); % Pyramidal cells

% Create tables for mixed-effects modeling
% Initialize empty arrays to store the data
responses = [];
drug_condition = {};
contrast_levels = [];
neuron_id = [];
mouse_id = {};
cell_type = {};

% Loop through drug conditions
for drug = 1:2
    if drug == 1
        drug_name = 'DART';
    else
        drug_name = 'Control';
    end
    
    % Get responses for this drug condition
    drug_responses = pref_responses_stat_concat{drug};
    
    % Loop through neurons
    for neuron = 1:size(drug_responses, 1)
        % Loop through contrast levels
        for contrast_idx = 1:3
            if contrast_idx == 1
                contrast = 25;
            elseif contrast_idx == 2
                contrast = 50;
            else
                contrast = 100;
            end
            
            % Add this data point
            responses = [responses; drug_responses(neuron, contrast_idx)];
            drug_condition = [drug_condition; drug_name];
            contrast_levels = [contrast_levels; contrast];
            neuron_id = [neuron_id; neuron];
            mouse_id = [mouse_id; mouseID(neuron)];
            
            if ismember(neuron, sst_indices)
                cell_type = [cell_type; 'SST'];
            else
                cell_type = [cell_type; 'PYR'];
            end
        end
    end
end

% Create one data table with all information
data_table = table(responses, drug_condition, contrast_levels, neuron_id, mouse_id, cell_type);

% Split data by cell type
sst_data = data_table(strcmp(data_table.cell_type, 'SST'), :);
pyr_data = data_table(strcmp(data_table.cell_type, 'PYR'), :);

% Convert categorical variables with proper reference levels
% Set Control as the reference level for drug_condition
sst_data.drug_condition = categorical(sst_data.drug_condition, {'Control', 'DART'});
sst_data.mouse_id = categorical(sst_data.mouse_id);

pyr_data.drug_condition = categorical(pyr_data.drug_condition, {'Control', 'DART'});
pyr_data.mouse_id = categorical(pyr_data.mouse_id);

% Fit mixed-effects models for each cell type
% For SST cells
lme_sst_stationary = fitlme(sst_data, 'responses ~ drug_condition * contrast_levels + (1|mouse_id) + (1|neuron_id)', 'FitMethod', 'REML');
lme_sst_stationary_noMouse = fitlme(sst_data, 'responses ~ drug_condition * contrast_levels +  (1|neuron_id)', 'FitMethod', 'REML');

% For Pyramidal cells
lme_pyr_stationary = fitlme(pyr_data, 'responses ~ drug_condition * contrast_levels + (1|mouse_id) + (1|neuron_id)', 'FitMethod', 'REML');
lme_pyr_stationary_noMouse = fitlme(pyr_data, 'responses ~ drug_condition * contrast_levels + (1|neuron_id)', 'FitMethod', 'REML');

exportLMEResultsTable(lme_sst_stationary, lme_pyr_stationary);
%% log likelihood ratio for random effect of mouse
% Extract log-likelihood values
ll_test(lme_sst_stationary,lme_sst_stationary_noMouse)

ll_test(lme_pyr_stationary,lme_pyr_stationary_noMouse)

%% model SST cells by noise correlation

% Separate data by cell type and noise correlation
sst_indices = find(red_concat == 1); % SST cells

% Create tables for mixed-effects modeling
% Initialize empty arrays to store the data
responses = [];
drug_condition = {};
contrast_levels = [];
neuron_id = [];
mouse_id = {};
noise_corr_group = {};

% Loop through drug conditions
for drug = 1:2
    if drug == 1
        drug_name = 'DART';
    else
        drug_name = 'Control';
    end
    
    % Get responses for this drug condition
    drug_responses = pref_responses_stat_concat{drug};
    
    % Loop through SST neurons only
    for i = 1:length(sst_indices)
        neuron = sst_indices(i);
        
        % Determine noise correlation group
        if ismember(neuron, redHigh)
            nc_group = 'High';
        elseif ismember(neuron, redLow)
            nc_group = 'Low';
        else
            continue; % Skip if not in either high or low group
        end
        
        % Loop through contrast levels
        for contrast_idx = 1:3
            if contrast_idx == 1
                contrast = 25;
            elseif contrast_idx == 2
                contrast = 50;
            else
                contrast = 100;
            end
            
            % Add this data point
            responses = [responses; drug_responses(neuron, contrast_idx)];
            drug_condition = [drug_condition; drug_name];
            contrast_levels = [contrast_levels; contrast];
            neuron_id = [neuron_id; neuron];
            mouse_id = [mouse_id; mouseID(neuron)];
            noise_corr_group = [noise_corr_group; nc_group];
        end
    end
end

% Create one data table with all information
sst_data_by_nc = table(responses, drug_condition, contrast_levels, neuron_id, mouse_id, noise_corr_group);

% Convert categorical variables with proper reference levels
sst_data_by_nc.drug_condition = categorical(sst_data_by_nc.drug_condition, {'Control', 'DART'});
sst_data_by_nc.mouse_id = categorical(sst_data_by_nc.mouse_id);
sst_data_by_nc.noise_corr_group = categorical(sst_data_by_nc.noise_corr_group, {'Low', 'High'});

% Split data by noise correlation group
sst_high_data = sst_data_by_nc(sst_data_by_nc.noise_corr_group == 'High', :);
sst_low_data = sst_data_by_nc(sst_data_by_nc.noise_corr_group == 'Low', :);

% Fit mixed-effects models for each noise correlation group
% For high noise correlation SST cells
lme_sst_high = fitlme(sst_high_data, 'responses ~ drug_condition * contrast_levels + (1|mouse_id) + (1|neuron_id)', 'FitMethod', 'REML');

% For low noise correlation SST cells
lme_sst_low = fitlme(sst_low_data, 'responses ~ drug_condition * contrast_levels + (1|mouse_id) + (1|neuron_id)', 'FitMethod', 'REML');

% Also fit a model with noise correlation as an interaction term
lme_sst_nc_interaction = fitlme(sst_data_by_nc, 'responses ~ drug_condition * contrast_levels * noise_corr_group + (1|mouse_id) + (1|neuron_id)', 'FitMethod', 'REML');

% Export results
fprintf('===== HIGH NOISE CORRELATION SST CELLS =====\n');
disp(lme_sst_high);

fprintf('\n===== LOW NOISE CORRELATION SST CELLS =====\n');
disp(lme_sst_low);

fprintf('\n===== INTERACTION MODEL WITH NOISE CORRELATION =====\n');
disp(lme_sst_nc_interaction);

exportLMEResultsTable(lme_sst_high, lme_sst_low);




%% Model for contrast X drug condition X behavioral state
% pref_responses_loc_concat is another cell array outcome variable in the
% same format at pref_responses_stat_concat. pref_responses_stat_concat has
% neural responses when mice are stationary and pref_responses_loc_concat
% has responses when mice are running, so these measure two different
% behavioral states. For each cell type, I want to run a mixed effects
% model for neural response with behavioral state, drug condition and 
% contrast as the fixed effects and neuron and mouse as random effects.
% Combine stationary and running data into a single table
% Initialize empty arrays to store the data
responses = [];
drug_condition = {};
contrast_levels = [];
behavioral_state = {};
neuron_id = [];
mouse_id = {};
cell_type = {};

% 1. STATIONARY STATE MODEL
% Create tables for mixed-effects modeling - stationary state
% Initialize empty arrays to store the data
responses_stat = [];
drug_condition_stat = {};
contrast_levels_stat = [];
neuron_id_stat = [];
mouse_id_stat = {};
cell_type_stat = {};

% Loop through drug conditions
for drug = 1:2
    if drug == 1
        drug_name = 'DART';
    else
        drug_name = 'Control';
    end
    
    % Get responses for this drug condition
    drug_responses = pref_responses_stat_concat{drug};
    
    % Loop through neurons
    for neuron = 1:size(drug_responses, 1)
        % Loop through contrast levels
        for contrast_idx = 1:3
            if contrast_idx == 1
                contrast = 25;
            elseif contrast_idx == 2
                contrast = 50;
            else
                contrast = 100;
            end
            
            % Add this data point
            responses_stat = [responses_stat; drug_responses(neuron, contrast_idx)];
            drug_condition_stat = [drug_condition_stat; drug_name];
            contrast_levels_stat = [contrast_levels_stat; contrast];
            neuron_id_stat = [neuron_id_stat; neuron];
            mouse_id_stat = [mouse_id_stat; mouseID(neuron)];
            
            % Assign cell type based on new index lists
            if ismember(neuron, red_all)
                cell_type_stat = [cell_type_stat; 'SST'];
            elseif ismember(neuron, green_all)
                cell_type_stat = [cell_type_stat; 'PYR'];
            else
                cell_type_stat = [cell_type_stat; 'Other']; % In case neuron is in neither list
            end
        end
    end
end

% Create one data table with all information
data_table_stat = table(responses_stat, drug_condition_stat, contrast_levels_stat, neuron_id_stat, mouse_id_stat, cell_type_stat);

% Split data by cell type
sst_data_stat = data_table_stat(strcmp(data_table_stat.cell_type_stat, 'SST'), :);
pyr_data_stat = data_table_stat(strcmp(data_table_stat.cell_type_stat, 'PYR'), :);

% Convert categorical variables with proper reference levels
% Set Control as the reference level for drug_condition
sst_data_stat.drug_condition_stat = categorical(sst_data_stat.drug_condition_stat, {'Control', 'DART'});
sst_data_stat.mouse_id_stat = categorical(sst_data_stat.mouse_id_stat);

pyr_data_stat.drug_condition_stat = categorical(pyr_data_stat.drug_condition_stat, {'Control', 'DART'});
pyr_data_stat.mouse_id_stat = categorical(pyr_data_stat.mouse_id_stat);

% Fit mixed-effects models for each cell type - stationary state
% For SST cells
lme_sst_stat = fitlme(sst_data_stat, 'responses_stat ~ drug_condition_stat * contrast_levels_stat + (1|mouse_id_stat) + (1|neuron_id_stat)', 'FitMethod', 'REML');
lme_sst_stat_reduced = fitlme(sst_data_stat, 'responses_stat ~ drug_condition_stat * contrast_levels_stat + (1|neuron_id_stat)', 'FitMethod', 'REML');

% For Pyramidal cells
lme_pyr_stat = fitlme(pyr_data_stat, 'responses_stat ~ drug_condition_stat * contrast_levels_stat + (1|mouse_id_stat) + (1|neuron_id_stat)', 'FitMethod', 'REML');
lme_pyr_stat_reduced = fitlme(pyr_data_stat, 'responses_stat ~ drug_condition_stat * contrast_levels_stat + (1|neuron_id_stat)', 'FitMethod', 'REML');

exportLMEResultsTable(lme_sst_stat,lme_pyr_stat)


ll_test(lme_sst_stat,lme_sst_stat_reduced)

ll_test(lme_pyr_stat,lme_pyr_stat_reduced)


%% 2. RUNNING STATE MODEL
% Initialize empty arrays to store the data
responses_run = [];
drug_condition_run = {};
contrast_levels_run = [];
neuron_id_run = [];
mouse_id_run = {};
cell_type_run = {};

% Loop through drug conditions
for drug = 1:2
    if drug == 1
        drug_name = 'DART';
    else
        drug_name = 'Control';
    end
    
    % Get responses for this drug condition (running state)
    drug_responses = pref_responses_loc_concat{drug};
    
    % Loop through neurons
    for neuron = 1:size(drug_responses, 1)
        % Loop through contrast levels
        for contrast_idx = 1:3
            if contrast_idx == 1
                contrast = 25;
            elseif contrast_idx == 2
                contrast = 50;
            else
                contrast = 100;
            end
            
            % Add this data point
            responses_run = [responses_run; drug_responses(neuron, contrast_idx)];
            drug_condition_run = [drug_condition_run; drug_name];
            contrast_levels_run = [contrast_levels_run; contrast];
            neuron_id_run = [neuron_id_run; neuron];
            mouse_id_run = [mouse_id_run; mouseID(neuron)];
            
            % Assign cell type based on new index lists
            if ismember(neuron, red_all)
                cell_type_run = [cell_type_run; 'SST'];
            elseif ismember(neuron, green_all)
                cell_type_run = [cell_type_run; 'PYR'];
            else
                cell_type_run = [cell_type_run; 'Other']; % In case neuron is in neither list
            end
        end
    end
end

% Create one data table with all information
data_table_run = table(responses_run, drug_condition_run, contrast_levels_run, neuron_id_run, mouse_id_run, cell_type_run);

% Split data by cell type
sst_data_run = data_table_run(strcmp(data_table_run.cell_type_run, 'SST'), :);
pyr_data_run = data_table_run(strcmp(data_table_run.cell_type_run, 'PYR'), :);

% Convert categorical variables with proper reference levels
% Set Control as the reference level for drug_condition
sst_data_run.drug_condition_run = categorical(sst_data_run.drug_condition_run, {'Control', 'DART'});
sst_data_run.mouse_id_run = categorical(sst_data_run.mouse_id_run);

pyr_data_run.drug_condition_run = categorical(pyr_data_run.drug_condition_run, {'Control', 'DART'});
pyr_data_run.mouse_id_run = categorical(pyr_data_run.mouse_id_run);

% Fit mixed-effects models for each cell type - running state
% For SST cells
lme_sst_run = fitlme(sst_data_run, 'responses_run ~ drug_condition_run * contrast_levels_run + (1|mouse_id_run) + (1|neuron_id_run)', 'FitMethod', 'REML');
lme_sst_run_reduced = fitlme(sst_data_run, 'responses_run ~ drug_condition_run * contrast_levels_run + (1|neuron_id_run)', 'FitMethod', 'REML');

% For Pyramidal cells
lme_pyr_run = fitlme(pyr_data_run, 'responses_run ~ drug_condition_run * contrast_levels_run + (1|mouse_id_run) + (1|neuron_id_run)', 'FitMethod', 'REML');
lme_pyr_run_reduced = fitlme(pyr_data_run, 'responses_run ~ drug_condition_run * contrast_levels_run + (1|neuron_id_run)', 'FitMethod', 'REML');
exportLMEResultsTable(lme_sst_run,lme_pyr_run)

ll_test(lme_sst_run,lme_sst_run_reduced)

ll_test(lme_pyr_run,lme_pyr_run_reduced)

%% 3. COMBINED MODEL (STATIONARY AND RUNNING)
% Initialize empty arrays to store the data
responses = [];
drug_condition = {};
contrast_levels = [];
behavioral_state = {};
neuron_id = [];
mouse_id = {};
cell_type = {};

% Define behavioral states
behav_states = {'Stationary', 'Running'};

% Process data from both behavioral states
for state_idx = 1:2
    if state_idx == 1
        state_name = 'Stationary';
        state_data = pref_responses_stat_concat;
    else
        state_name = 'Running';
        state_data = pref_responses_loc_concat;
    end
    
    % Loop through drug conditions
    for drug = 1:2
        if drug == 1
            drug_name = 'DART';
        else
            drug_name = 'Control';
        end
        
        % Get responses for this drug condition
        drug_responses = state_data{drug};
        
        % Loop through neurons
        for neuron = 1:size(drug_responses, 1)
            % Loop through contrast levels
            for contrast_idx = 1:3
                if contrast_idx == 1
                    contrast = 25;
                elseif contrast_idx == 2
                    contrast = 50;
                else
                    contrast = 100;
                end
                
                % Add this data point
                responses = [responses; drug_responses(neuron, contrast_idx)];
                drug_condition = [drug_condition; drug_name];
                contrast_levels = [contrast_levels; contrast];
                behavioral_state = [behavioral_state; state_name];
                neuron_id = [neuron_id; neuron];
                mouse_id = [mouse_id; mouseID(neuron)];
                
                % Assign cell type based on new index lists
                if ismember(neuron, red_all)
                    cell_type = [cell_type; 'SST'];
                elseif ismember(neuron, green_all)
                    cell_type = [cell_type; 'PYR'];
                else
                    cell_type = [cell_type; 'Other']; % In case neuron is in neither list
                end
            end
        end
    end
end

% Create one data table with all information
data_table = table(responses, drug_condition, contrast_levels, behavioral_state, neuron_id, mouse_id, cell_type);

% Split data by cell type
sst_data = data_table(strcmp(data_table.cell_type, 'SST'), :);
pyr_data = data_table(strcmp(data_table.cell_type, 'PYR'), :);

% Convert categorical variables with proper reference levels
% For SST cells
sst_data.drug_condition = categorical(sst_data.drug_condition, {'Control', 'DART'});  % Set Control as reference
sst_data.behavioral_state = categorical(sst_data.behavioral_state, {'Stationary', 'Running'});  % Set Stationary as reference
sst_data.mouse_id = categorical(sst_data.mouse_id);

% For Pyramidal cells
pyr_data.drug_condition = categorical(pyr_data.drug_condition, {'Control', 'DART'});  % Set Control as reference
pyr_data.behavioral_state = categorical(pyr_data.behavioral_state, {'Stationary', 'Running'});  % Set Stationary as reference
pyr_data.mouse_id = categorical(pyr_data.mouse_id);

% Fit mixed-effects models for each cell type with three-way interaction
% For SST cells
lme_sst_full = fitlme(sst_data, 'responses ~ drug_condition * contrast_levels * behavioral_state + (1|mouse_id) + (1|neuron_id)', 'FitMethod', 'REML');
lme_sst_reduced = fitlme(sst_data, 'responses ~ drug_condition * contrast_levels * behavioral_state +  (1|neuron_id)', 'FitMethod', 'REML');

% For Pyramidal cells
lme_pyr_full = fitlme(pyr_data, 'responses ~ drug_condition * contrast_levels * behavioral_state + (1|mouse_id) + (1|neuron_id)', 'FitMethod', 'REML');
lme_pyr_reduced = fitlme(pyr_data, 'responses ~ drug_condition * contrast_levels * behavioral_state  + (1|neuron_id)', 'FitMethod', 'REML');
exportLMEResultsTable(lme_sst_full,lme_pyr_full)

ll_test(lme_sst_full,lme_sst_reduced)

ll_test(lme_pyr_full,lme_pyr_reduced)