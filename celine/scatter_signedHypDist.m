function scatter_signedHypDist(dataset, pre, post, interneuron_indices, pyramidal_indices, filename)
% scatter_signedHypDist Creates scatter plots with signed hypotenuse probability distributions for neural responses and saves to PDF
%
% INPUTS:
%   dataset - Cell array containing neural response data
%   pre - Index for pre condition in dataset
%   post - Index for post condition in dataset
%   interneuron_indices - Indices for interneurons
%   pyramidal_indices - Indices for pyramidal cells
%   filename - Base filename for saving PDFs
%
% This function creates scatter plots of pre vs. post neural responses
% for both interneurons and pyramidal cells, along with probability density plots of post-pre differences,
% and saves each figure as a PDF.

% Get number of contrast conditions
[~, nCon] = size(dataset{pre});

% FIRST PASS: Find global min/max values across all contrasts for consistent axes
int_min_global = Inf;
int_max_global = -Inf;
pyr_min_global = Inf;
pyr_max_global = -Inf;

% Calculate global axis limits for interneurons and pyramidal cells
for iCon = 1:nCon
    % Extract data for this contrast condition
    pre_interneuron = dataset{pre}(interneuron_indices, iCon);
    post_interneuron = dataset{post}(interneuron_indices, iCon);
    pre_pyramidal = dataset{pre}(pyramidal_indices, iCon);
    post_pyramidal = dataset{post}(pyramidal_indices, iCon);
    
    % Update global min/max for interneurons
    int_min_raw = min([pre_interneuron; post_interneuron]);
    int_min_global = min(int_min_global, int_min_raw);
    
    int_max_raw = max([pre_interneuron; post_interneuron]);
    int_max_global = max(int_max_global, int_max_raw);
    
    % Update global min/max for pyramidal cells
    pyr_min_raw = min([pre_pyramidal; post_pyramidal]);
    pyr_min_global = min(pyr_min_global, pyr_min_raw);
    
    pyr_max_raw = max([pre_pyramidal; post_pyramidal]);
    pyr_max_global = max(pyr_max_global, pyr_max_raw);
end

% Round down global min to integer with 1 decimal place
int_min_global = floor(int_min_global * 10) / 10;
pyr_min_global = floor(pyr_min_global * 10) / 10;

% Round up global max to integer with 1 decimal place and double it
int_max_global = ceil(int_max_global * 10) / 10 * 2;
pyr_max_global = ceil(pyr_max_global * 10) / 10 * 2;

% Set consistent axis limits
int_axis_lim = [int_min_global int_max_global];
pyr_axis_lim = [pyr_min_global pyr_max_global];

% Calculate axis limits for signed hypotenuse distributions based on global ranges
int_leg_length = int_axis_lim(2) - int_axis_lim(1);
pyr_leg_length = pyr_axis_lim(2) - pyr_axis_lim(1);

% Calculate hypotenuse of right triangle with leg length = x-axis range of scatterplot
int_hyp_bound = int_leg_length / 2; % Half of diagonal length for consistent scaling
pyr_hyp_bound = pyr_leg_length / 2; % Half of diagonal length for consistent scaling

% Set symmetric axis limits for signed hypotenuse
hyp_int_axis_lim = [-int_hyp_bound int_hyp_bound];
hyp_pyr_axis_lim = [-pyr_hyp_bound pyr_hyp_bound];

% SECOND PASS: Create plots using consistent axis limits
for iCon = 1:nCon
    % Create a new figure with a 2x2 grid layout with explicit size in inches
    fig = figure('Units', 'inches', 'Position', [1 1 4 4]);
    
    % Extract data for this contrast condition
    pre_interneuron = dataset{pre}(interneuron_indices, iCon);
    post_interneuron = dataset{post}(interneuron_indices, iCon);
    pre_pyramidal = dataset{pre}(pyramidal_indices, iCon);
    post_pyramidal = dataset{post}(pyramidal_indices, iCon);
    
    % Calculate signed hypotenuse (direction-aware Euclidean distance) between pre and post responses
    % Calculate magnitude (standard hypotenuse)
    hyp_interneuron = sqrt(pre_interneuron.^2 + post_interneuron.^2);
    hyp_pyramidal = sqrt(pre_pyramidal.^2 + post_pyramidal.^2);
    
    % Add sign based on whether pre > post (reversed from previous version)
    signed_hyp_interneuron = hyp_interneuron .* sign(pre_interneuron - post_interneuron);
    signed_hyp_pyramidal = hyp_pyramidal .* sign(pre_pyramidal - post_pyramidal);
    
    % Calculate x-values for signed hypotenuse distributions 
    x_smooth_hyp_int = linspace(hyp_int_axis_lim(1), hyp_int_axis_lim(2), 100);
    x_smooth_hyp_pyr = linspace(hyp_pyr_axis_lim(1), hyp_pyr_axis_lim(2), 100);
    
    % Compute smoothed distributions for signed hypotenuse
    hyp_int_density = ksdensity_safe(signed_hyp_interneuron, x_smooth_hyp_int);
    hyp_pyr_density = ksdensity_safe(signed_hyp_pyramidal, x_smooth_hyp_pyr);
    
    % Normalize densities to peak at 1.0 for better comparison
    hyp_int_density = hyp_int_density / max(hyp_int_density);
    hyp_pyr_density = hyp_pyr_density / max(hyp_pyr_density);
    
    % Define colors for scatter plots (match these colors with distributions)
    interneuron_color = 'black';
    pyramidal_color = 'black';
    
    % Interneuron scatter - top left
    subplot(2, 2, 1);
    % Use smaller marker size (15) with white edge and black fill
    scatter(pre_interneuron, post_interneuron, 15, 'MarkerFaceColor', interneuron_color, 'MarkerEdgeColor', 'white', 'MarkerFaceAlpha', 0.5, 'LineWidth', 0.5);
    hold on;
    
    % Add unity line and set axis properties
    plot(int_axis_lim, int_axis_lim, 'k--');
    xlim(int_axis_lim);
    ylim(int_axis_lim);
    axis square
    % xlabel('Pre Response', 'FontSize', 10);
    % ylabel('Post Response', 'FontSize', 10);
    % title('Interneurons');
    
    % Create exactly 3 tick marks at min, middle, max with one decimal place
    % Ensure middle tick is exactly halfway between min and max
    int_middle = (int_axis_lim(1) + int_axis_lim(2)) / 2;
    x_ticks = [int_axis_lim(1), int_middle, int_axis_lim(2)];
    y_ticks = x_ticks;
    
    xticks(x_ticks);
    yticks(y_ticks);
    xticklabels(compose('%.1f', x_ticks));
    yticklabels(compose('%.1f', y_ticks));
    set(gca, 'TickDir', 'out', 'FontSize', 8, 'XTickLabelRotation', 0, 'YTickLabelRotation', 0)
    box off
    
    % Interneuron signed hypotenuse distribution - top right
    subplot(2, 2, 2);
    area(x_smooth_hyp_int, hyp_int_density, 'FaceColor', interneuron_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    xlim(hyp_int_axis_lim);
    ylim([0 1.1]);  % Fixed y-axis for normalized density
    % title('Interneurons (Signed Hypotenuse)');
    % xlabel('Signed Magnitude (+ if Pre > Post)', 'FontSize', 10);
    % ylabel('Normalized Density', 'FontSize', 10);
    % Add vertical line at zero
    hold on;
    plot([0 0], [0 1.1], 'k--');
    % Use equally spaced ticks with zero in the middle
    x_ticks = [hyp_int_axis_lim(1), 0, hyp_int_axis_lim(2)];
    xticks(x_ticks);
    xticklabels(compose('%.1f', x_ticks));
    set(gca, 'TickDir', 'out', 'FontSize', 8, 'XTickLabelRotation', 0, 'YTickLabelRotation', 0)
    box off
    
    % Pyramidal cell scatter - bottom left
    subplot(2, 2, 3);
    % Use smaller marker size (15) with white edge and black fill
    scatter(pre_pyramidal, post_pyramidal, 15, 'MarkerFaceColor', pyramidal_color, 'MarkerEdgeColor', 'white', 'MarkerFaceAlpha', 0.5, 'LineWidth', 0.5);
    hold on;
    
    % Add unity line and set axis properties
    plot(pyr_axis_lim, pyr_axis_lim, 'k--');
    xlim(pyr_axis_lim);
    ylim(pyr_axis_lim);
    axis square
    % xlabel('Pre Response', 'FontSize', 10);
    % ylabel('Post Response', 'FontSize', 10);
    % title('Pyramidal Cells');
    
    % Create exactly 3 tick marks at min, middle, max with one decimal place
    % Ensure middle tick is exactly halfway between min and max
    pyr_middle = (pyr_axis_lim(1) + pyr_axis_lim(2)) / 2;
    x_ticks = [pyr_axis_lim(1), pyr_middle, pyr_axis_lim(2)];
    y_ticks = x_ticks;
    
    xticks(x_ticks);
    yticks(y_ticks);
    xticklabels(compose('%.1f', x_ticks));
    yticklabels(compose('%.1f', y_ticks));
    set(gca, 'TickDir', 'out', 'FontSize', 8, 'XTickLabelRotation', 0, 'YTickLabelRotation', 0)
    box off
    
    % Pyramidal signed hypotenuse distribution - bottom right
    subplot(2, 2, 4);
    area(x_smooth_hyp_pyr, hyp_pyr_density, 'FaceColor', pyramidal_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    xlim(hyp_pyr_axis_lim);
    ylim([0 1.1]);  % Fixed y-axis for normalized density
    % title('Pyramidal Cells (Signed Hypotenuse)');
    % xlabel('Signed Magnitude (+ if Pre > Post)', 'FontSize', 10);
    % ylabel('Normalized Density', 'FontSize', 10);
    % Add vertical line at zero
    hold on;
    plot([0 0], [0 1.1], 'k--');
    % Use equally spaced ticks with zero in the middle
    x_ticks = [hyp_pyr_axis_lim(1), 0, hyp_pyr_axis_lim(2)];
    xticks(x_ticks);
    xticklabels(compose('%.1f', x_ticks));
    set(gca, 'TickDir', 'out', 'FontSize', 8, 'XTickLabelRotation', 0, 'YTickLabelRotation', 0)
    box off
    
    % Set scatterplot size to 0.75 inch square
    ax1 = subplot(2, 2, 1);
    ax1.Units = 'inches';
    pos1 = ax1.Position;
    pos1(3:4) = [0.75 0.75]; % width and height in inches
    ax1.Position = pos1;
    
    ax3 = subplot(2, 2, 3);
    ax3.Units = 'inches';
    pos3 = ax3.Position;
    pos3(3:4) = [0.75 0.75]; % width and height in inches
    ax3.Position = pos3;
    
    % Adjust PDF plots size accordingly
    ax2 = subplot(2, 2, 2);
    ax2.Units = 'inches';
    pos2 = ax2.Position;
    pos2(3:4) = [sqrt(2)*0.75 0.25]; % Scale width proportionally to the scatterplot
    ax2.Position = pos2;
    
    ax4 = subplot(2, 2, 4);
    ax4.Units = 'inches';
    pos4 = ax4.Position;
    pos4(3:4) = [sqrt(2)*0.75 0.25]; % Scale width proportionally to the scatterplot
    ax4.Position = pos4;
    
    % Add contrast condition information
    sgtitle(['Contrast Condition: ' num2str(iCon)]);
    
    % Save figure as PDF with explicit size control
    fullfilename = sprintf([filename '_scatter_signed_hyp_contrast_%d.pdf'], iCon);
    
    % Set paper size to match figure size
    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperSize', [fig.Position(3) fig.Position(4)]);
    set(fig, 'PaperPositionMode', 'manual');
    set(fig, 'PaperPosition', [0 0 fig.Position(3) fig.Position(4)]);
    
    % Save with high resolution
    print(fig, fullfilename, '-dpdf', '-r300');
    
    fprintf('Saved figure for contrast %d as %s\n', iCon, fullfilename);
end
end

function y = ksdensity_safe(data, xi)
% Helper function for safe kernel density estimation
% Falls back to histogram-based approach if ksdensity fails

try
    % Try to compute kernel density estimate
    y = ksdensity(data, xi);
catch
    % If ksdensity fails (e.g., insufficient unique points), fallback to histogram
    [counts, edges] = histcounts(data, 10, 'Normalization', 'probability');
    centers = (edges(1:end-1) + edges(2:end))/2;
    y = interp1(centers, counts, xi, 'linear', 0);
    % Normalize to ensure it integrates to 1
    y = y / trapz(xi, y);
end
end