%% Figure S3
clear all
cons=[.25 .5 1] %the visual stimulus contrasts used
nCon = length(cons);
pre=2; post=1; %indexing the two imaging days
load('FigureS3AB_data.mat')
load('FigureS3_data.mat')

%% Figure S3A - Distribution of correlation values
figure;
histogram(control_noiseCorr(1,sst_ind),"NumBins",11)
%% Figure S3B -  Normalized difference vs. noise correlation scatter
%This plots the scatters for all contrasts; 25% is in Figure 4, 50% and
%100% are in Figure S3

% Create structures to store results
slopes = zeros(nCon, 1);
intercepts = zeros(nCon, 1);
r_squared = zeros(nCon, 1);
p_values = zeros(nCon, 1);
std_errors = zeros(nCon, 1);

% Create figure with fixed size
figure('Units', 'inches', 'Position', [0, 0, 5, 2], 'PaperUnits', 'inches', 'PaperSize', [5, 2]);

% Use tiledlayout
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'Unweighted Regression: Noise Correlation vs. Normalized Difference', 'FontSize', 10);

% Analyze each iCon
for iCon = 1:nCon
    % Extract data for this iCon
    y = stationary_norm_diff(iCon, :)';  % neural response difference
    x = control_noiseCorr(1,sst_ind)';                 % noise correlations
    
    % Handle NaNs and missing data
    valid_idx = ~isnan(x) & ~isnan(y) & isfinite(x) & isfinite(y);
    
    % Check if we have enough valid data points
    if sum(valid_idx) < 3
        warning('Not enough valid data points for iCon %d. Skipping.', iCon);
        continue;
    end
    
    % Use only valid data
    x_valid = x(valid_idx);
    y_valid = y(valid_idx);
    
    % Perform unweighted linear regression
    [p, stats] = polyfit(x_valid, y_valid, 1);
    
    % Store results
    intercepts(iCon) = p(2);
    slopes(iCon) = p(1);
    
    % Calculate standard error of the slope
    yfit = polyval(p, x_valid);
    residuals = y_valid - yfit;
    SSE = sum(residuals.^2);
    n = length(x_valid);
    
    % Standard error of the regression
    SE_regression = sqrt(SSE / (n-2));
    
    % Sum of squares of x deviations
    SS_x = sum((x_valid - mean(x_valid)).^2);
    
    % Standard error of the slope
    std_errors(iCon) = SE_regression / sqrt(SS_x);
    
    % Calculate R-squared
    SS_total = sum((y_valid - mean(y_valid)).^2);
    r_squared(iCon) = 1 - SSE/SS_total;
    
    % Calculate p-value for slope
    t_stat = slopes(iCon) / std_errors(iCon);
    p_values(iCon) = 2 * (1 - tcdf(abs(t_stat), length(x_valid) - 2));
    
    % Create a nexttile for this iCon
    ax = nexttile;
    
    % Create scatter plot with smaller grey points and white outlines
    scatter(x_valid, y_valid, 15, 'MarkerFaceColor','black','MarkerEdgeColor', 'white', 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.5);
    hold on;
    
    % Add regression line
    x_range = linspace(min(x_valid), max(x_valid), 100);
    y_line = p(1) * x_range + p(2);
    plot(x_range, y_line, 'r-', 'LineWidth', 2);
    
    % Labels and title
    title(sprintf('Contrast: %d%%', cons(iCon)), 'FontSize', 8);
    xlabel('Noise Correlation (Baseline)', 'FontSize', 7);
    ylabel('Normalized Difference', 'FontSize', 7);
    
    % Set axis limits
    xlim([-.2 1]);
    ylim([-8 8]);

    % Add horizontal line at y=0
    hline(0);
    
    % Add regression equation and stats - using smaller font size
    % Calculate text position dynamically
    y_range = range(ylim);
    text_x = -0.15;  % Fixed position to align across plots
    text_y_start = 7;  % Fixed position near top
    text_y_step = 1;  % Fixed step size
    
    text(text_x, text_y_start, sprintf('y = %.3fx + %.3f', slopes(iCon), intercepts(iCon)), 'FontSize', 6);
    text(text_x, text_y_start - text_y_step, sprintf('R^2 = %.3f', r_squared(iCon)), 'FontSize', 6);
    text(text_x, text_y_start - 2*text_y_step, sprintf('p = %.4f', p_values(iCon)), 'FontSize', 6);
    text(text_x, text_y_start - 3*text_y_step, sprintf('n = %d', sum(valid_idx)), 'FontSize', 6);
    
    % Set box style and tick direction
    box off;
    set(gca, 'TickDir', 'out');
    grid off;
    
    % Make the axis ticks and numbers smaller
    set(gca, 'FontSize', 6);
    
    % Make sure the plots have equal sizes
    axis square;
end

% Create a results table
contrast_labels = {'25%', '50%', '100%'};
results_table = table(contrast_labels', slopes, std_errors, intercepts, r_squared, p_values, ...
    'VariableNames', {'Contrast', 'Slope', 'StdError', 'Intercept', 'RSquared', 'PValue'});

% Display results table
disp('Unweighted Regression Results:');
disp(results_table);


%% Figure S3C timecourses at 25% contrast
%identify the strongly and weakly correlated SST cells
highRInds = find(control_noiseCorr_PEG(1,:)>0.5);
lowRInds = find(control_noiseCorr_PEG(1,:)<=0.5);
sstLow=intersect(lowRInds, sst_ind_PEG);
sstHigh=intersect(highRInds, sst_ind_PEG);

plotNeuralTimecourse_singleContrast(tc_trial_avrg_stat_concat_PEG, tc_trial_avrg_stat_concat_PEG, ...
    sstLow, sstHigh, 1,...
    'UseDashedLines', [false, false], ...
    'Colors1', {'k', 'c'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'c'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'Weak corr. SST', 'Strong corr. SST'}, ...
    'StimStart', 31);
%% Figure S3D - contrast response plots
plotContrastResponse(pref_responses_stat_concat_PEG, pref_responses_stat_concat_PEG, ...
    sstLow, sstHigh, cons, ...
    'UseDashedLines', [false, false], ...  % Dashed lines for the right plot
    'Titles', {'Weak corr. SST', 'Strong corr. SST'}, ...
    'YLabel', 'dF/F','Colors', {'k', 'c'});
%% Figure S3E -  Normalized difference vs. noise correlation scatter
%This plots the scatters for all contrasts; 25% is in Figure 4, 50% and
%100% are in Figure S3

% Create structures to store results
slopes = zeros(nCon, 1);
intercepts = zeros(nCon, 1);
r_squared = zeros(nCon, 1);
p_values = zeros(nCon, 1);
std_errors = zeros(nCon, 1);

% Create figure with fixed size
figure('Units', 'inches', 'Position', [0, 0, 5, 2], 'PaperUnits', 'inches', 'PaperSize', [5, 2]);

% Use tiledlayout
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'Unweighted Regression: Noise Correlation vs. Normalized Difference', 'FontSize', 10);

% Analyze each iCon
for iCon = 1:nCon
    % Extract data for this iCon
    y = stationary_norm_diff_PEG(iCon, :)';  % neural response difference
    x = control_noiseCorr_PEG(1,sst_ind_PEG)';                 % noise correlations
    
    % Handle NaNs and missing data
    valid_idx = ~isnan(x) & ~isnan(y) & isfinite(x) & isfinite(y);
    
    % Check if we have enough valid data points
    if sum(valid_idx) < 3
        warning('Not enough valid data points for iCon %d. Skipping.', iCon);
        continue;
    end
    
    % Use only valid data
    x_valid = x(valid_idx);
    y_valid = y(valid_idx);
    
    % Perform unweighted linear regression
    [p, stats] = polyfit(x_valid, y_valid, 1);
    
    % Store results
    intercepts(iCon) = p(2);
    slopes(iCon) = p(1);
    
    % Calculate standard error of the slope
    yfit = polyval(p, x_valid);
    residuals = y_valid - yfit;
    SSE = sum(residuals.^2);
    n = length(x_valid);
    
    % Standard error of the regression
    SE_regression = sqrt(SSE / (n-2));
    
    % Sum of squares of x deviations
    SS_x = sum((x_valid - mean(x_valid)).^2);
    
    % Standard error of the slope
    std_errors(iCon) = SE_regression / sqrt(SS_x);
    
    % Calculate R-squared
    SS_total = sum((y_valid - mean(y_valid)).^2);
    r_squared(iCon) = 1 - SSE/SS_total;
    
    % Calculate p-value for slope
    t_stat = slopes(iCon) / std_errors(iCon);
    p_values(iCon) = 2 * (1 - tcdf(abs(t_stat), length(x_valid) - 2));
    
    % Create a nexttile for this iCon
    ax = nexttile;
    
    % Create scatter plot with smaller grey points and white outlines
    scatter(x_valid, y_valid, 15, 'MarkerFaceColor','black','MarkerEdgeColor', 'white', 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.5);
    hold on;
    
    % Add regression line
    x_range = linspace(min(x_valid), max(x_valid), 100);
    y_line = p(1) * x_range + p(2);
    plot(x_range, y_line, 'r-', 'LineWidth', 2);
    
    % Labels and title
    title(sprintf('Contrast: %d%%', cons(iCon)), 'FontSize', 8);
    xlabel('Noise Correlation (Baseline)', 'FontSize', 7);
    ylabel('Normalized Difference', 'FontSize', 7);
    
    % Set axis limits
    xlim([-.2 1]);
    ylim([-8 8]);

    % Add horizontal line at y=0
    hline(0);
    
    % Add regression equation and stats - using smaller font size
    % Calculate text position dynamically
    y_range = range(ylim);
    text_x = -0.15;  % Fixed position to align across plots
    text_y_start = 7;  % Fixed position near top
    text_y_step = 1;  % Fixed step size
    
    text(text_x, text_y_start, sprintf('y = %.3fx + %.3f', slopes(iCon), intercepts(iCon)), 'FontSize', 6);
    text(text_x, text_y_start - text_y_step, sprintf('R^2 = %.3f', r_squared(iCon)), 'FontSize', 6);
    text(text_x, text_y_start - 2*text_y_step, sprintf('p = %.4f', p_values(iCon)), 'FontSize', 6);
    text(text_x, text_y_start - 3*text_y_step, sprintf('n = %d', sum(valid_idx)), 'FontSize', 6);
    
    % Set box style and tick direction
    box off;
    set(gca, 'TickDir', 'out');
    grid off;
    
    % Make the axis ticks and numbers smaller
    set(gca, 'FontSize', 6);
    
    % Make sure the plots have equal sizes
    axis square;
end

% Create a results table
contrast_labels = {'25%', '50%', '100%'};
results_table = table(contrast_labels', slopes, std_errors, intercepts, r_squared, p_values, ...
    'VariableNames', {'Contrast', 'Slope', 'StdError', 'Intercept', 'RSquared', 'PValue'});

% Display results table
disp('Unweighted Regression Results:');
disp(results_table);

