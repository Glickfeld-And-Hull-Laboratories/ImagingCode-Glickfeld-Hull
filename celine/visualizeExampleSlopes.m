function visualizeExampleSlopes(norm_diff, slopes, numExamples)
    % visualizeExampleSlopes - Visualize the calculated slopes for example neurons
    %
    % Input:
    %   norm_diff - 3D array of normalized difference values
    %   slopes - 2D array of calculated slopes
    %   numExamples - Number of example neurons to visualize (default: 4)
    
    % If numExamples is not provided, default to 4
    if nargin < 3
        numExamples = 4;
    end
    
    % Get dimensions
    [num_states, ~, num_neurons] = size(norm_diff);
    
    % Find non-NaN neurons (those with valid slopes in both states)
    valid_neurons = find(~isnan(slopes(1,:)) & ~isnan(slopes(2,:)));
    
    % If no valid neurons found, exit
    if isempty(valid_neurons)
        disp('No neurons with valid slopes in both states found.');
        return;
    end
    
    % Select numExamples random neurons from valid ones
    if length(valid_neurons) <= numExamples
        example_neurons = valid_neurons;
    else
        % Pick random indices without replacement
        rng(42); % For reproducibility
        rand_indices = randperm(length(valid_neurons), numExamples);
        example_neurons = valid_neurons(rand_indices);
    end
    
    % Contrast values (x-axis) - now only using endpoints
    contrast_values = [25, 50, 100];
    
    % Set up the figure
    figure('Position', [100, 100, 1200, 800]);
    
    % Create a subplot for each example neuron
    for i = 1:length(example_neurons)
        neuron = example_neurons(i);
        
        % Create a subplot for this neuron
        subplot(2, ceil(length(example_neurons)/2), i);
        
        % Plot data and slope lines for both states
        for state = 1:num_states
            % Extract data for this neuron/state
            y_values = squeeze(norm_diff(state, :, neuron));
            
            % Determine color based on state
            if state == 1
                color = 'b'; % Blue for stationary
                state_name = 'Stationary';
            else
                color = 'r'; % Red for running
                state_name = 'Running';
            end
            
            % Plot all data points (including the middle point)
            scatter(contrast_values, y_values, 50, color, 'filled', 'MarkerFaceAlpha', 0.7);
            hold on;
            
            % Highlight the endpoints that were used for slope calculation
            scatter([contrast_values(1), contrast_values(3)], [y_values(1), y_values(3)], 80, color, 'LineWidth', 2);
            
            % Plot the slope line connecting only the 25% and 100% points
            plot([contrast_values(1), contrast_values(3)], [y_values(1), y_values(3)], color, 'LineWidth', 2, 'DisplayName', sprintf('%s (Slope: %.4f)', state_name, slopes(state, neuron)));
        end
        
        % Add labels and title
        xlabel('Contrast (%)');
        ylabel('Normalized Difference');
        title(sprintf('Neuron %d', neuron));
        legend('Location', 'best');
        grid on;
        
        % Set reasonable axis limits
        ylim_current = ylim;
        ylim_padding = (ylim_current(2) - ylim_current(1)) * 0.2;
        ylim([ylim_current(1) - ylim_padding, ylim_current(2) + ylim_padding]);
        xlim([20, 105]);
    end
    
    % Add overall title
    sgtitle('Example Neurons: Contrast Responses and 25%-to-100% Slopes');
end