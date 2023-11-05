function all_resps = process_auditory_signals(allcells, condition_flag, varargin)


% process the variable inputs
p = inputParser;
p.addParameter('smoothing_scale', 3, @isnumeric)
p.addParameter('num_sum', 3, @isnumeric)
p.addParameter('stim_start', 20, @isnumeric)
p.addParameter('stim_end', 56, @isnumeric)
p.addParameter('verbose', false, @islogical)
p.parse(varargin{:});

% dole out the inputs to variables for convenience
smoothing_scale = p.Results.smoothing_scale;
num_sum = p.Results.num_sum;
stim_start = p.Results.stim_start;
stim_end = p.Results.stim_end;
verbose = p.Results.verbose;

% make filter for smoothing the calcium traces
ca_filter = [zeros(1,smoothing_scale), ones(1,smoothing_scale), zeros(1,smoothing_scale)];
ca_filter = ca_filter ./ norm(ca_filter);

% figure out the number of cells and conditions that need to be looped
num_conditions = length(condition_flag);
num_cells = length(allcells);

% loop over cells
for gc = 1:num_cells
    
    %preallocate
    temp_all_init_resps = cell(1,num_conditions);
    temp_all_peak_resps = cell(1,num_conditions);
    temp_all_integrate_resps = cell(1,num_conditions);
    
    for stim_cond = 1:num_conditions
    
        % extract the responses for a single cell under a stim condition
        if condition_flag == 1
            temp_gc = cell2mat(allcells{gc}(1));
        else
            temp_gc = cell2mat(allcells{gc}{condition_flag(stim_cond)});
        end
        
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

all_resps.all_init_resps = all_init_resps;
all_resps.all_peak_resps = all_peak_resps;
all_resps.all_integrate_resps = all_integrate_resps;


