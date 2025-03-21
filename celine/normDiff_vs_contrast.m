% Create a new figure
figure;
weights = bsln_std(1,1,red_ind_concat); %a matrix of the baseline standard deviations that I'll use to wieght values in the regression
% Extract x and y data for regression
x_data = noiseCorr_OG_concat{pre}(1,red_ind_concat);
y_data = norm_diff(1,1,red_ind_concat);

% Check dimensions and ensure they match
x_size = size(x_data);
y_size = size(y_data);

fprintf('Size of x_data: [%s]\n', num2str(x_size));
fprintf('Size of y_data: [%s]\n', num2str(y_size));

% Make sure both are column vectors
x_data = x_data(:);
y_data = y_data(:);
weights = weights(:); % Ensure weights is also a column vector

% Check if they have different lengths (but don't trim)
if length(x_data) ~= length(y_data)
    error('ERROR: x_data and y_data have different lengths (%d vs %d). This will cause indexing errors.', length(x_data), length(y_data));
end

% Remove NaN values for regression analysis
valid_idx = ~isnan(x_data) & ~isnan(y_data) & ~isnan(weights) & ~isinf(x_data) & ~isinf(y_data) & ~isinf(weights);
x_valid = x_data(valid_idx);
y_valid = y_data(valid_idx);
weights_valid = weights(valid_idx);

% Check if we have enough valid data points
if length(x_valid) < 2
    error('Not enough valid data points for regression analysis');
end

% Create scatter plot with sizes based on weights
scatter(x_data, y_data, 50*weights/max(weights), 'filled', 'MarkerFaceAlpha', 0.7);

% Hold the plot to add the regression line
hold on;

% Compute weighted linear regression on valid data
% Prepare design matrix X with column of ones for intercept
X = [ones(length(x_valid), 1), x_valid];

% Calculate weighted least squares parameters: (X'WX)^-1 X'Wy
W = diag(weights_valid);
beta = (X' * W * X) \ (X' * W * y_valid);

% Extract slope and intercept
slope = beta(2);
intercept = beta(1);

% Create a sequence of x values for smoother line plotting
x_range = linspace(min(x_valid), max(x_valid), 100);
y_fit = intercept + slope * x_range;

% Plot regression line
plot(x_range, y_fit, 'r-', 'LineWidth', 2);

% Calculate weighted R-squared
y_mean = sum(weights_valid .* y_valid) / sum(weights_valid); % Weighted mean
y_pred = intercept + slope * x_valid;
weighted_SS_total = sum(weights_valid .* (y_valid - y_mean).^2);
weighted_SS_resid = sum(weights_valid .* (y_valid - y_pred).^2);
Rsq = 1 - weighted_SS_resid / weighted_SS_total;

% Calculate standard error of the slope
n = length(x_valid);
dof = n - 2; % Degrees of freedom
MSE = weighted_SS_resid / dof;
var_beta = MSE * inv(X' * W * X);
se_slope = sqrt(var_beta(2,2));

% Calculate t-statistic and p-value for the slope
t_stat = slope / se_slope;
slope_p_value = 2 * (1 - tcdf(abs(t_stat), dof));

% Calculate weighted correlation coefficient
x_mean = sum(weights_valid .* x_valid) / sum(weights_valid);
weighted_cov = sum(weights_valid .* (x_valid - x_mean) .* (y_valid - y_mean)) / sum(weights_valid);
weighted_var_x = sum(weights_valid .* (x_valid - x_mean).^2) / sum(weights_valid);
weighted_var_y = sum(weights_valid .* (y_valid - y_mean).^2) / sum(weights_valid);
weighted_corr = weighted_cov / sqrt(weighted_var_x * weighted_var_y);
p_value = slope_p_value; % For weighted regression, use slope p-value

% Display statistics on the plot
equation = sprintf('Slope = %.4f', slope);
r_squared = sprintf('R^2 = %.4f', Rsq);
corr_p_val = sprintf('Corr p = %.4g', p_value);
reg_p_val = sprintf('Reg p = %.4g', slope_p_value);
n_points = sprintf('n = %d', n);

% Create text box with statistics
dim = [.2 .02 .3 .3];
str = {equation, r_squared, corr_p_val, reg_p_val, n_points};
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'white');

% Print statistics to command window
fprintf('\nWeighted Regression Statistics:\n');
fprintf('Beta: %.4f\n', slope);
fprintf('Intercept: %.4f\n', intercept);
fprintf('R-squared: %.4f\n', Rsq);
fprintf('Correlation (weighted): %.4f\n', weighted_corr);
fprintf('Regression p-value: %.6f\n', slope_p_value);
fprintf('Sample size: %d\n', n);

% Calculate prediction intervals (95%)
% Note: This is an approximation for weighted regression
t_crit = tinv(0.975, dof);
pred_se = sqrt(MSE * (1 + 1/n + ((x_range - x_mean).^2) / sum((x_valid - x_mean).^2)));
plot(x_range, y_fit + t_crit * pred_se, 'r--', 'LineWidth', 1);
plot(x_range, y_fit - t_crit * pred_se, 'r--', 'LineWidth', 1);

% Customize appearance
set(gca, 'TickDir', 'out');
box off;
xlabel('R value');
ylabel('Norm Diff 25%');
title('Weighted Regression: R value vs Normalized Difference 25%');

% Add legend
legend('Data points (size by weight)', 'Regression line', '95% Prediction interval', 'Location', 'best');

hold off;

%%
% Script to save MATLAB data structures in formats that can be easily uploaded
% This includes saving as CSV files and a combined JSON file

% Ensure output directory exists
if ~exist('export_data', 'dir')
    mkdir('export_data');
end

% Part 1: Save each structure as separate CSV files
% For 3D matrices bsln_std and norm_diff:
% Save each behavioral state and contrast combination separately

% For bsln_std
for state = 1:2
    for contrast = 1:3
        % Extract the slice for this state and contrast
        data_slice = bsln_std(state, contrast, :);
        % Reshape to a column vector
        data_slice = reshape(data_slice, [], 1);
        % Create filename
        filename = sprintf('export_data/bsln_std_state%d_contrast%d.csv', state, contrast);
        % Write to CSV
        writematrix(data_slice, filename);
    end
end

% For norm_diff
for state = 1:2
    for contrast = 1:3
        % Extract the slice for this state and contrast
        data_slice = norm_diff(state, contrast, :);
        % Reshape to a column vector
        data_slice = reshape(data_slice, [], 1);
        % Create filename
        filename = sprintf('export_data/norm_diff_state%d_contrast%d.csv', state, contrast);
        % Write to CSV
        writematrix(data_slice, filename);
    end
end

% For noiseCorr_OG_concat
for cell_idx = 1:2
    % Extract the entire cell
    cell_data = noiseCorr_OG_concat{cell_idx};
    % Save each state row
    for state = 1:2
        % Extract the row for this state
        data_row = cell_data(state, :);
        % Create filename
        filename = sprintf('export_data/noiseCorr_cell%d_state%d.csv', cell_idx, state);
        % Write to CSV
        writematrix(data_row, filename);
    end
end

% Part 2: Save as MAT file with version 7 (more compatible)
save('export_data/all_data_v7.mat', 'bsln_std', 'norm_diff', 'noiseCorr_OG_concat', '-v7');

% Part 3: Save a JSON representation of the data (for even more compatibility)
% Note: This requires a custom function to convert MATLAB arrays to JSON
% First, create a structure to hold the data
data_struct = struct();

% Create a simpler structure for bsln_std
bsln_std_struct = struct();
for state = 1:2
    for contrast = 1:3
        fieldname = sprintf('state%d_contrast%d', state, contrast);
        bsln_std_struct.(fieldname) = reshape(bsln_std(state, contrast, :), [], 1);
    end
end
data_struct.bsln_std = bsln_std_struct;

% Create a simpler structure for norm_diff
norm_diff_struct = struct();
for state = 1:2
    for contrast = 1:3
        fieldname = sprintf('state%d_contrast%d', state, contrast);
        norm_diff_struct.(fieldname) = reshape(norm_diff(state, contrast, :), [], 1);
    end
end
data_struct.norm_diff = norm_diff_struct;

% Create a simpler structure for noiseCorr_OG_concat
noiseCorr_struct = struct();
for cell_idx = 1:2
    cell_struct = struct();
    for state = 1:2
        fieldname = sprintf('state%d', state);
        cell_struct.(fieldname) = noiseCorr_OG_concat{cell_idx}(state, :);
    end
    noiseCorr_struct.(sprintf('cell%d', cell_idx)) = cell_struct;
end
data_struct.noiseCorr_OG_concat = noiseCorr_struct;

% Save the structure as JSON
json_str = jsonencode(data_struct);
fid = fopen('export_data/all_data.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

disp('All data has been exported to the export_data directory.');
disp('Files saved:');
disp(' - CSV files for each slice of the data');
disp(' - all_data_v7.mat (MATLAB format v7)');
disp(' - all_data.json (JSON format)');

% Part 4: Create a README file explaining the data
readme_content = [
    'Data Export Readme\n', ...
    '==================\n\n', ...
    'This directory contains exported data from MATLAB structures:\n\n', ...
    '1. bsln_std: 3D matrix [behavioral state (2), stimulus contrast (3), neuron]\n', ...
    '   - state1 = stationary, state2 = running\n', ...
    '   - contrast1 = 25%, contrast2 = 50%, contrast3 = 100%\n\n', ...
    '2. norm_diff: 3D matrix [behavioral state (2), stimulus contrast (3), neuron]\n', ...
    '   - state1 = stationary, state2 = running\n', ...
    '   - contrast1 = 25%, contrast2 = 50%, contrast3 = 100%\n\n', ...
    '3. noiseCorr_OG_concat: Cell array (2 cells), each containing a 2D matrix [behavioral state (2), neuron]\n', ...
    '   - cell1 = experimental, cell2 = baseline\n', ...
    '   - state1 = stationary, state2 = running\n\n', ...
    'File formats:\n', ...
    '- CSV files: Each slice of the data is saved as a separate CSV file\n', ...
    '- MAT file: all_data_v7.mat (MATLAB format v7 for better compatibility)\n', ...
    '- JSON file: all_data.json (JSON format for non-MATLAB use)\n'
];

fid = fopen('export_data/README.txt', 'w');
fprintf(fid, '%s', readme_content);
fclose(fid);

disp(' - README.txt (explanation of the data)');
%% 
% I want to perform and plot weighted regressions with norm_diff as the dependent variable and
% the baseline day of noiseCorr_OG_concat as the independent variable, in the stationary state,
% one regression for each contrast. Use bsln_std for the weights.

% Weighted Regression Analysis for Neural Data
% Analyze relation between noise correlations and normalized difference
% With higher weight for smaller standard deviations and robustness to NaNs
mySample = red_ind_concat;
% Weighted Regression Analysis for Neural Data
% Analyze relation between noise correlations and normalized difference
% With higher weight for smaller standard deviations and robustness to NaNs

% Extract data (stationary state only)
stationary_norm_diff = squeeze(norm_diff(1,:,mySample)); % [contrast, neuron]
stationary_bsln_std = squeeze(bsln_std(1,:,mySample));   % [contrast, neuron]
baseline_noise_corr = noiseCorr_OG_concat{2}(1,mySample); % [1, neuron]

% Number of contrasts and neurons
[num_contrasts, num_neurons] = size(stationary_norm_diff);

% Store results
slopes = zeros(num_contrasts, 1);
intercepts = zeros(num_contrasts, 1);
r_squared = zeros(num_contrasts, 1);
p_values = zeros(num_contrasts, 1);
std_errors = zeros(num_contrasts, 1);

% Create figure
figure('Position', [100, 100, 1200, 400]);

% Analyze each contrast
for contrast = 1:num_contrasts
    % Extract data for this contrast
    y = stationary_norm_diff(contrast, :)';  % neural response difference
    x = baseline_noise_corr';                 % noise correlations
    
    % Calculate weights - higher weight for smaller std (inverse variance)
    % But with moderation to prevent extremely high weights
    raw_weights = 1 ./ (stationary_bsln_std(contrast, :)' .^ 2);
    
    % Apply a cap to prevent excessive weighting
    % Method 1: Use a percentile-based cap (e.g., 95th percentile)
    weight_cap = prctile(raw_weights(~isnan(raw_weights) & isfinite(raw_weights)), 95);
    weights = min(raw_weights, weight_cap);
    
    % Alternative method (commented out): Apply a square root transformation to compress the range
    % weights = sqrt(raw_weights);
    
    % Handle NaNs and missing data
    valid_idx = ~isnan(x) & ~isnan(y) & ~isnan(weights) & isfinite(x) & isfinite(y) & isfinite(weights);
    
    % Check if we have enough valid data points
    if sum(valid_idx) < 3
        warning('Not enough valid data points for contrast %d. Skipping.', contrast);
        title(sprintf('Contrast %d: Insufficient Data', contrast));
        continue;
    end
    
    % Use only valid data
    x_valid = x(valid_idx);
    y_valid = y(valid_idx);
    weights_valid = weights(valid_idx);
    
    % Create design matrix for regression
    X = [ones(length(x_valid), 1), x_valid];
    
    % Perform weighted regression
    [b, stats] = lscov(X, y_valid, weights_valid);
    
    % Store results
    intercepts(contrast) = b(1);
    slopes(contrast) = b(2);
    std_errors(contrast) = stats(2);
    
    % Calculate weighted R-squared
    y_hat = b(1) + b(2) * x_valid;
    weighted_mean = sum(weights_valid .* y_valid) / sum(weights_valid);
    weighted_SST = sum(weights_valid .* (y_valid - weighted_mean).^2);
    weighted_SSE = sum(weights_valid .* (y_valid - y_hat).^2);
    r_squared(contrast) = 1 - weighted_SSE / weighted_SST;
    
    % Calculate p-value for slope
    t_stat = b(2) / stats(2);
    p_values(contrast) = 2 * (1 - tcdf(abs(t_stat), length(x_valid) - 2));
    
    % Plot regression
    subplot(1, 3, contrast);
    
    % Create a normalized weight for visualization (bubble size)
    % Keep the relationship correct: smaller std = larger bubble
    % Use log scale for better visualization distribution
    log_weights = log10(weights_valid);
    max_log_weight = max(log_weights);
    min_log_weight = min(log_weights);
    range_log_weight = max_log_weight - min_log_weight;
    
    if range_log_weight > 0
        normalized_weights = 30 + 120 * (log_weights - min_log_weight) / range_log_weight;
    else
        normalized_weights = 80 * ones(size(weights_valid)); % Default size if all weights are equal
    end
    
    % Ensure minimum and maximum sizes are reasonable
    normalized_weights = max(normalized_weights, 20);  % Minimum point size
    normalized_weights = min(normalized_weights, 200); % Maximum point size
    
    % Create scatter plot with size based on weights
    % Use log scale for color mapping to better visualize the range
    scatter(x_valid, y_valid, normalized_weights, log10(weights_valid), 'filled', 'MarkerFaceAlpha', 0.7);
    colormap(gca, 'parula');
    hold on;
    
    % Add regression line
    x_range = linspace(min(x_valid), max(x_valid), 100);
    y_line = b(1) + b(2) * x_range;
    plot(x_range, y_line, 'r-', 'LineWidth', 2);
    
    % Labels and title
    contrasts = [25, 50, 100];
    title(sprintf('Contrast: %d%%', contrasts(contrast)), 'FontSize', 12);
    xlabel('Noise Correlation (Baseline)', 'FontSize', 11);
    ylabel('Normalized Difference', 'FontSize', 11);
    
    % % Set axis limits with some padding
    % x_range = max(x_valid) - min(x_valid);
    % y_range = max(y_valid) - min(y_valid);
    % x_padding = 0.15 * x_range;
    % y_padding = 0.15 * y_range;
    % xlim([min(x_valid) - x_padding, max(x_valid) + x_padding]);
    % ylim([min(y_valid) - y_padding, max(y_valid) + y_padding]);

    xlim([-.2 1])
    ylim([-8 8])

    hline(0)
    
    % Add regression equation and stats
    text_x = min(x_valid);
    text_y_start = max(y_valid) - 0.1 * y_range;
    text_y_step = 0.1 * y_range;
    
    text(text_x, text_y_start, sprintf('y = %.3fx + %.3f', b(2), b(1)), 'FontSize', 10);
    text(text_x, text_y_start - text_y_step, sprintf('R� = %.3f', r_squared(contrast)), 'FontSize', 10);
    text(text_x, text_y_start - 2*text_y_step, sprintf('p = %.4f', p_values(contrast)), 'FontSize', 10);
    
    % Add sample size
    text(text_x, text_y_start - 3*text_y_step, sprintf('n = %d', sum(valid_idx)), 'FontSize', 10);
    
    % Add colorbar for weights (blue = higher weight = lower std)
    c = colorbar;
    c.Label.String = 'log10(Weight)';
    c.Label.FontSize = 9;
    
    % Add text indicating weighting approach
    annotation('textbox', [0.01, 0.01, 0.4, 0.05], 'String', ...
        'Weights capped at 95th percentile to prevent excessive influence', ...
        'FontSize', 8, 'FitBoxToText', 'on', 'EdgeColor', 'none');
    
    
    % Set box style and tick direction
    box off;
    set(gca, 'TickDir', 'out');
end

% Adjust spacing
sgtitle('Weighted Regression: Noise Correlation vs. Normalized Difference (Stationary State)', 'FontSize', 14);
set(gcf, 'Color', 'w');

% Create a results table
contrast_labels = {'25%', '50%', '100%'};
results_table = table(contrast_labels', slopes, std_errors, intercepts, r_squared, p_values, ...
    'VariableNames', {'Contrast', 'Slope', 'StdError', 'Intercept', 'RSquared', 'PValue'});

% Display results table
disp('Weighted Regression Results:');
disp(results_table);

% Save figure
saveas(gcf, 'weighted_regression_plot.png');
fprintf('Figure saved as weighted_regression_plot.png\n');

% Optional: Save results to CSV
writetable(results_table, 'weighted_regression_results.csv');
fprintf('Results saved to weighted_regression_results.csv\n');

% Additional analysis: Test if slopes are significantly different across contrasts
fprintf('\nComparing slopes across contrast levels:\n');

% Pairwise comparison of slopes (z-test)
for i = 1:num_contrasts
    for j = i+1:num_contrasts
        % Calculate Z-statistic for difference between slopes
        z_diff = (slopes(i) - slopes(j)) / sqrt(std_errors(i)^2 + std_errors(j)^2);
        p_diff = 2 * (1 - normcdf(abs(z_diff)));
        fprintf('Contrast %s vs %s: Difference = %.3f, z = %.3f, p = %.4f\n', ...
            contrast_labels{i}, contrast_labels{j}, slopes(i) - slopes(j), z_diff, p_diff);
    end
end

% Additional figure: Compare slopes across contrasts
figure('Position', [100, 500, 500, 400]);
bar(1:num_contrasts, slopes);
hold on;

% Add error bars
errorbar(1:num_contrasts, slopes, std_errors, 'k.', 'LineWidth', 1.5);

% Add labels and title
xlabel('Contrast', 'FontSize', 12);
ylabel('Regression Slope', 'FontSize', 12);
title('Comparison of Regression Slopes Across Contrasts', 'FontSize', 14);
xticks(1:num_contrasts);
xticklabels(contrast_labels);
grid on;
box off;
set(gca, 'TickDir', 'out');

% Save comparison figure
saveas(gcf, 'slope_comparison.png');
fprintf('Slope comparison figure saved as slope_comparison.png\n');