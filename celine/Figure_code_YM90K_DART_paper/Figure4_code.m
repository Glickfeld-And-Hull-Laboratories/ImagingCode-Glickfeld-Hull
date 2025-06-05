%% Figure 4
cons=[.25 .5 1] %the visual stimulus contrasts used
nCon = length(cons);
pre=2; post=1; %indexing the two imaging days
load("Figure4A_data.mat");
load("Figure4_data.mat");
%% Figure 4A - example cell scatters
% Example cell 1 = cell 44
% Example cell 2  = cell 49
% Both from i2052
for iCell = 1:size(example_sst_meanSub,2)
    thisCell = example_sst_meanSub(:,iCell);
    [R,p]=corrcoef(example_pyr_meanSub,thisCell,'rows','complete');
        figure
        scatter(example_pyr_meanSub,thisCell, 'MarkerFaceColor','black','MarkerEdgeColor','none','MarkerFaceAlpha', 0.5)
        hold on
        h = lsline;
        title([' R= ' num2str(R(2))]);
        ylabel('Mean-subtracted activity for individual SST cell')
        xlabel('Mean-subtracted activity for all Pyr cells')
        box off;
        set(gca, 'TickDir', 'out');
        grid off;
end
%% Figure 4B - timecourses at 25% contrast
%identify the strongly and weakly correlated SST cells
highRInds = find(control_noiseCorr(1,:)>0.5);
lowRInds = find(control_noiseCorr(1,:)<=0.5);
sstLow=intersect(lowRInds, sst_ind);
sstHigh=intersect(highRInds, sst_ind);

plotNeuralTimecourse_singleContrast(tc_trial_avrg_stat_concat, tc_trial_avrg_stat_concat, ...
    sstLow, sstHigh, 1,...
    'UseDashedLines', [false, false], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'Weak corr. SST', 'Strong corr. SST'}, ...
    'StimStart', 31);
%% Figure 4C - contrast response plots
plotContrastResponse(pref_responses_stat_concat, pref_responses_stat_concat, ...
    sstLow, sstHigh, cons, ...
    'UseDashedLines', [false, false], ...  % Dashed lines for the right plot
    'Titles', {'Weak corr. SST', 'Strong corr. SST'}, ...
    'YLabel', 'dF/F');
%% Figure 4C - ANOVA 
w = table(categorical([1 1 1 2 2 2 ].'), categorical([1 2 3 1 2 3].'), 'VariableNames', {'DART', 'contrast'}); % within-design
rm_SST_low = fitrm(SST_low_stat_dfof, 'd1c1-d2c3 ~ 1', 'WithinDesign', w)
ranova(rm_SST_low, 'withinmodel', 'DART*contrast')

rm_SST_high = fitrm(SST_high_stat_dfof, 'd1c1-d2c3 ~ 1', 'WithinDesign', w)
ranova(rm_SST_high, 'withinmodel', 'DART*contrast')

%% Figure 4D - Normalized difference vs. noise correlation scatter
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

%% Figure 4E - Comparing slopes 
% Additional figure: Compare slopes across contrasts
figure('Position', [100, 500, 500, 400]);
bar(1:nCon, slopes);
ylim([-2,2])
hold on;

% Add error bars
errorbar(1:nCon, slopes, std_errors, 'k.', 'LineWidth', 1.5);

% Add labels and title
xlabel('Contrast', 'FontSize', 12);
ylabel('Regression Slope', 'FontSize', 12);
title('Comparison of Regression Slopes Across Contrasts', 'FontSize', 14);
xticks(1:nCon);
xticklabels(contrast_labels);
grid on;
box off;
set(gca, 'TickDir', 'out');
%% Figure 4E - Test linear relationship between beta values and contrast using linear model

% Extract the beta values (slopes) and contrast levels
beta_values = slopes;  % These are already calculated in your original code
contrast_levels = [25, 50, 100]';  % Make it a column vector

% Fit linear model
mdl = fitlm(contrast_levels, beta_values);
fprintf('Linear model summary:\n');
disp(mdl);

