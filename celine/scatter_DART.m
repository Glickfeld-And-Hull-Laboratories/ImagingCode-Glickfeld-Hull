function scatter_DART(dataset, pre, post, interneuron_indices, pyramidal_indices)
% scatter_DART Creates scatter plots for neural responses and saves to PDF
%
% INPUTS:
%   dataset - Cell array containing neural response data
%   pre - Index for pre condition in dataset
%   post - Index for post condition in dataset
%   interneuron_indices - Indices for interneurons
%   pyramidal_indices - Indices for pyramidal cells
%
% This function creates scatter plots of pre vs. post neural responses
% for both interneurons and pyramidal cells, and saves each figure as a PDF.

% Get number of contrast conditions
[~, nCon] = size(dataset{pre});

% Loop through each contrast condition
for iCon = 1:nCon
    % Create a new figure
    fig = figure;
    
    % Extract data for this contrast condition
    pre_interneuron = dataset{pre}(interneuron_indices, iCon);
    post_interneuron = dataset{post}(interneuron_indices, iCon);
    pre_pyramidal = dataset{pre}(pyramidal_indices, iCon);
    post_pyramidal = dataset{post}(pyramidal_indices, iCon);
    
    % Create interneuron plot with its own axis limits
    subplot(1, 2, 1);
    scatter(pre_interneuron, post_interneuron, 'MarkerFaceColor','black','MarkerEdgeColor','white','MarkerFaceAlpha', 0.5);
    hold on;
    
    % Calculate axis limits for interneurons (same for x and y)
    int_min = round(min([pre_interneuron; post_interneuron]),1);
    int_max = round(max([pre_interneuron; post_interneuron]),1);
    int_axis_lim = [int_min int_max];
    
    % Add unity line and set axis properties
    plot(int_axis_lim, int_axis_lim, 'k--');
    xlim(int_axis_lim);
    ylim(int_axis_lim);
    axis square
    xlabel('Control', 'FontSize', 10);
    ylabel('YM90KDART', 'FontSize', 10);
    title('Interneurons');
    % Ensure ticks at min and max values
    xticks([int_axis_lim(1), (int_axis_lim(1)+int_axis_lim(2))/2, int_axis_lim(2)]);
    yticks([int_axis_lim(1), (int_axis_lim(1)+int_axis_lim(2))/2, int_axis_lim(2)]);
    set(gca, 'TickDir', 'out', 'FontSize', 8)
    grid off
    box off
    
    % Create pyramidal cell plot with its own axis limits
    subplot(1, 2, 2);
    scatter(pre_pyramidal, post_pyramidal, 'MarkerFaceColor','black','MarkerEdgeColor','white','MarkerFaceAlpha', 0.5);
    hold on;
    
    % Calculate axis limits for pyramidal cells (same for x and y)
    pyr_min = round(min([pre_pyramidal; post_pyramidal]),1);
    pyr_max = round(max([pre_pyramidal; post_pyramidal]),1);
    pyr_axis_lim = [pyr_min pyr_max];
    
    % Add unity line and set axis properties
    plot(pyr_axis_lim, pyr_axis_lim, 'k--');
    xlim(pyr_axis_lim);
    ylim(pyr_axis_lim);
    axis square
    xlabel('Pre Response', 'FontSize', 10);
    ylabel('Post Response', 'FontSize', 10);
    title('Pyramidal Cells');
    % Ensure ticks at min and max values
    xticks([pyr_axis_lim(1), (pyr_axis_lim(1)+pyr_axis_lim(2))/2, pyr_axis_lim(2)]);
    yticks([pyr_axis_lim(1), (pyr_axis_lim(1)+pyr_axis_lim(2))/2, pyr_axis_lim(2)]);
    set(gca, 'TickDir', 'out', 'FontSize', 8)
    grid off
    box off
    
    % Add contrast condition information
    sgtitle(['Contrast Condition: ' num2str(iCon)]);
    
    % Set the figure size to make the plots 0.75 inches square
    % First set the figure size in inches
    set(fig, 'Units', 'inches');
    
    % Calculate needed figure width (2 plots × 0.75 inches plus some spacing)
    total_width = 0.75 * 2 + 0.5; % 0.75 inches per plot plus 0.5 inch for spacing/margins
    total_height = 0.75 + 0.5; % 0.75 inches for plot plus 0.5 inch for margins and title
    
    % Set figure size
    set(fig, 'Position', [1, 1, total_width, total_height]);
    
    % Adjust the subplot positions to ensure they are exactly 0.75 inches square
    pos1 = get(subplot(1,2,1), 'Position');
    pos2 = get(subplot(1,2,2), 'Position');
    
    % Calculate new positions (make width and height 0.75/total_width and 0.75/total_height of figure)
    subplot_width = 0.75/total_width;
    subplot_height = 0.75/total_height;
    
    % Position for first subplot (left side)
    set(subplot(1,2,1), 'Position', [0.1, 0.2, subplot_width, subplot_height]);
    
    % Position for second subplot (right side)
    set(subplot(1,2,2), 'Position', [0.1 + subplot_width + 0.1, 0.2, subplot_width, subplot_height]);
    
    % Save figure as PDF
    filename = sprintf('scatter_contrast_%d.pdf', iCon);
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, filename, '-dpdf', '-bestfit');
    
    fprintf('Saved figure for contrast %d as %s\n', iCon, filename);
end
end