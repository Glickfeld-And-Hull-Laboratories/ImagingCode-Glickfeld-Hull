function slopes = calculateContrastSlopes(norm_diff)
    % calculateContrastSlopes - Calculate slopes from 25% to 100% contrast for each neuron
    %
    % Input:
    %   norm_diff - 3D array where:
    %       1st dimension: behavioral state (1=stationary, 2=running)
    %       2nd dimension: contrast (1=25%, 2=50%, 3=100%)
    %       3rd dimension: neuron index
    %
    % Output:
    %   slopes - 2D array where:
    %       1st dimension: behavioral state (1=stationary, 2=running)
    %       2nd dimension: neuron index
    
    % Get dimensions
    [num_states, ~, num_neurons] = size(norm_diff);
    
    % Preallocate the slopes array
    slopes = zeros(num_states, num_neurons);
    
    % Calculate slope using only the 25% and 100% contrast points
    for state = 1:num_states
        for neuron = 1:num_neurons
            % Get values at 25% (index 1) and 100% (index 3)
            y_start = norm_diff(state, 1, neuron);
            y_end = norm_diff(state, 3, neuron);
            
            % Check if either point is NaN
            if isnan(y_start) || isnan(y_end)
                slopes(state, neuron) = NaN;
            else
                % Calculate rise over run from 25% to 100%
                rise = y_end - y_start;
                run = 100 - 25;  % contrast difference
                slopes(state, neuron) = rise / run;
            end
        end
    end
end