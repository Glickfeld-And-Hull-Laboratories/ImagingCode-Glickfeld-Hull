function plotSizeResponseBootstrap(data, red_ind, green_ind, contrast, sizes, varargin)
% PLOTSIZERESPONSEBOOTSTRAP Plot size tuning with random 2/3 subsampling across 5 iterations

p = inputParser;
addRequired(p, 'data', @iscell);
addRequired(p, 'red_ind', @isnumeric);
addRequired(p, 'green_ind', @isnumeric);
addRequired(p, 'contrast', @isnumeric);
addRequired(p, 'sizes', @isnumeric);
addParameter(p, 'UseDashedLines', [false, true], @(x) islogical(x) && length(x)==2);
addParameter(p, 'Titles', {'Dataset 1', 'Dataset 2'}, @(x) iscell(x) && length(x)==2);
addParameter(p, 'YLabel', 'dF/F', @ischar);
addParameter(p, 'Colors1', {'k', 'b'}, @(x) iscell(x) && length(x)==2);
addParameter(p, 'Colors2', {'k', 'b'}, @(x) iscell(x) && length(x)==2);
parse(p, data, red_ind, green_ind, contrast, sizes, varargin{:});

use_dashed_lines = p.Results.UseDashedLines;
titles = p.Results.Titles;
y_label = p.Results.YLabel;
colors1 = p.Results.Colors1;
colors2 = p.Results.Colors2;

pre = 2;
post = 1;
n_iterations = 5;
n_red = round(length(red_ind) * 2/3);
n_green = round(length(green_ind) * 2/3);

figure('Units', 'inches', 'Position', [2, 2, 8, 3*n_iterations]);

% Get highest contrast index (last one)
contrast_idx = size(data{pre}, 2);

fprintf('Random subsampling: selecting %d/%d red cells and %d/%d green cells\n', ...
    n_red, length(red_ind), n_green, length(green_ind));

% Initialize storage for selected indices
selected_cells = table('Size', [n_iterations 3], ...
    'VariableTypes', {'double', 'cell', 'cell'}, ...
    'VariableNames', {'Iteration', 'Red_Cells', 'Green_Cells'});

for iter = 1:n_iterations
    % Random subsample 2/3 of each cell type
    red_subset = red_ind(randperm(length(red_ind), n_red));
    green_subset = green_ind(randperm(length(green_ind), n_green));
    
    % Store selected indices
    selected_cells.Iteration(iter) = iter;
    selected_cells.Red_Cells{iter} = red_subset;
    selected_cells.Green_Cells{iter} = green_subset;
    
    fprintf('Iteration %d: Red cells [%s...] Green cells [%s...]\n', ...
        iter, num2str(red_subset(1:min(5,n_red))), num2str(green_subset(1:min(5,n_green))));
    
    % Extract data for highest contrast
    red_pre = squeeze(data{pre}(red_subset, contrast_idx, :));
    red_post = squeeze(data{post}(red_subset, contrast_idx, :));
    green_pre = squeeze(data{pre}(green_subset, contrast_idx, :));
    green_post = squeeze(data{post}(green_subset, contrast_idx, :));
    
    % Calculate mean and SEM
    red_pre_mean = mean(red_pre, 1, 'omitnan');
    red_pre_se = std(red_pre, 1, 'omitnan') / sqrt(n_red);
    red_post_mean = mean(red_post, 1, 'omitnan');
    red_post_se = std(red_post, 1, 'omitnan') / sqrt(n_red);
    
    green_pre_mean = mean(green_pre, 1, 'omitnan');
    green_pre_se = std(green_pre, 1, 'omitnan') / sqrt(n_green);
    green_post_mean = mean(green_post, 1, 'omitnan');
    green_post_se = std(green_post, 1, 'omitnan') / sqrt(n_green);
    
    % Determine y-axis limits for this iteration
    all_vals = [red_pre_mean - red_pre_se, red_pre_mean + red_pre_se, ...
                red_post_mean - red_post_se, red_post_mean + red_post_se, ...
                green_pre_mean - green_pre_se, green_pre_mean + green_pre_se, ...
                green_post_mean - green_post_se, green_post_mean + green_post_se];
    ymin = min(all_vals);
    ymax = max(all_vals);
    padding = 0.1 * (ymax - ymin);
    ylim_range = [ymin - padding, ymax + padding];
    
    % Plot red cells
    subplot(n_iterations, 2, (iter-1)*2 + 1);
    line_style = '-';
    if use_dashed_lines(1)
        line_style = '--';
    end
    errorbar(sizes, red_pre_mean, red_pre_se, [line_style, 'o'], ...
        'Color', colors1{1}, 'LineWidth', 1.5, 'MarkerSize', 6);
    hold on;
    errorbar(sizes, red_post_mean, red_post_se, [line_style, 'o'], ...
        'Color', colors1{2}, 'LineWidth', 1.5, 'MarkerSize', 6);
    
    xlim([min(sizes)*0.8, max(sizes)*1.2]);
    ylim(ylim_range);
    xticks(sizes);
    set(gca, 'TickDir', 'out', 'box', 'off');
    grid off;
    
    if iter == 1
        title([titles{1}, ' (Contrast: ', num2str(contrast), ') n = ', num2str(n_red)], 'FontWeight', 'normal');
    end
    if iter == n_iterations
        xlabel('Size (deg)');
    end
    ylabel(y_label);
    
    % Plot green cells
    subplot(n_iterations, 2, (iter-1)*2 + 2);
    line_style = '-';
    if use_dashed_lines(2)
        line_style = '--';
    end
    errorbar(sizes, green_pre_mean, green_pre_se, [line_style, 'o'], ...
        'Color', colors2{1}, 'LineWidth', 1.5, 'MarkerSize', 6);
    hold on;
    errorbar(sizes, green_post_mean, green_post_se, [line_style, 'o'], ...
        'Color', colors2{2}, 'LineWidth', 1.5, 'MarkerSize', 6);
    
    xlim([min(sizes)*0.8, max(sizes)*1.2]);
    ylim(ylim_range);
    xticks(sizes);
    set(gca, 'TickDir', 'out', 'box', 'off');
    grid off;
    
    if iter == 1
        title([titles{2}, ' (Contrast: ', num2str(contrast), ') n = ', num2str(n_green)], 'FontWeight', 'normal');
    end
    if iter == n_iterations
        xlabel('Size (deg)');
    end
end

% Save selected cells table
save('bootstrap_selected_cells.mat', 'selected_cells');
fprintf('\nSaved selected cell indices to bootstrap_selected_cells.mat\n');

end