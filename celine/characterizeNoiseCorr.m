%% instructions given to Claude
% pref_responses_loc_concat{pre} is a N neuron by 3 contrast matrix of
% neural responses during running and pref_responses_stat_concat{pre} is
% the same during stationary periods. I want to make a new matrix of
% locomotion modulation index for each neuron at each contrast, where LMI =
% (running response - stationary response) / (running response + stationary
% response). I also want a matrix of the mean LMI for each cell, averaging
% over contrasts. Next, do the same thing for arousal modulation index (AMI)
% using the matrices pref_responses_stat_largePupil_concat{pre} for neural
% responses when arousal is high and
% pref_responses_stat_smallPupil_concat{pre} for responses when arousal is
% low.

% Use pref_responses_stat_concat{pre} to find the slope of the line that
% best fits the change of response over contrast for each cell. The
% contrast values at 25%, 50%, and 100%.
%% Calculate LMI
% Initialize matrix for LMI values
LMI = zeros(size(pref_responses_loc_concat{pre}));

% Calculate LMI for each neuron at each contrast
for neuron = 1:size(pref_responses_loc_concat{pre}, 1)
    for contrast = 1:size(pref_responses_loc_concat{pre}, 2)
        % Get responses and set negative values to zero
        run_resp = max(0, pref_responses_loc_concat{pre}(neuron, contrast));
        stat_resp = max(0, pref_responses_stat_concat{pre}(neuron, contrast));
        
        % Calculate LMI = (running - stationary) / (running + stationary)
        if (run_resp + stat_resp) ~= 0
            LMI(neuron, contrast) = (run_resp - stat_resp) / (run_resp + stat_resp);
        else
            LMI(neuron, contrast) = 0; % If both responses are zero
        end
    end
end

% Calculate mean LMI for each neuron (averaging across contrasts)
meanLMI = mean(LMI, 2);

% check LMI distribution
%figure; histogram(meanLMI)

%% Calculate Arousal Modulation Index AMI
% Initialize matrix for AMI values
AMI = zeros(size(pref_responses_stat_largePupil_concat{pre}));

% Calculate AMI for each neuron at each contrast
for neuron = 1:size(pref_responses_stat_largePupil_concat{pre}, 1)
    for contrast = 1:size(pref_responses_stat_largePupil_concat{pre}, 2)
        % Get responses and set negative values to zero
        high_arousal = max(0, pref_responses_stat_largePupil_concat{pre}(neuron, contrast));
        low_arousal = max(0, pref_responses_stat_smallPupil_concat{pre}(neuron, contrast));
        
        % Calculate AMI = (high arousal - low arousal) / (high arousal + low arousal)
        if (high_arousal + low_arousal) ~= 0
            AMI(neuron, contrast) = (high_arousal - low_arousal) / (high_arousal + low_arousal);
        else
            AMI(neuron, contrast) = 0; % If both responses are zero
        end
    end
end

% Calculate mean AMI for each neuron (averaging across contrasts)
meanAMI = mean(AMI, 2);
% check LMI distribution
%figure; histogram(meanAMI)

%% Calculate contrast response slope
% best fit line is found with a linear regression
% Define contrast values (25%, 50%, 100%)
contrast_values = [25, 50, 100];

% Initialize vector to store slopes for each neuron
response_slopes = zeros(size(pref_responses_stat_concat{pre}, 1), 1);

% Calculate slope for each neuron
for neuron = 1:size(pref_responses_stat_concat{pre}, 1)
    % Extract responses for this neuron across all contrasts
    responses = pref_responses_stat_concat{pre}(neuron, :);
    
    % Use polyfit to find the line of best fit
    % First parameter is degree 1 for a linear fit
    p = polyfit(contrast_values, responses, 1);
    
    % The first coefficient is the slope
    response_slopes(neuron) = p(1);
end

% Optional: Visualize fits for a few example neurons
figure;
for i = 1:min(5, size(pref_responses_stat_concat{pre}, 1))
    subplot(min(5, size(pref_responses_stat_concat{pre}, 1)), 1, i);
    
    % Plot the actual data points
    plot(contrast_values, pref_responses_stat_concat{pre}(i, :), 'o');
    hold on;
    
    % Plot the fitted line
    p = polyfit(contrast_values, pref_responses_stat_concat{pre}(i, :), 1);
    fitted_line = polyval(p, [0 100]);
    plot([0 100], fitted_line, '-');
    
    title(['Neuron ' num2str(i) ', Slope = ' num2str(p(1))]);
    xlabel('Contrast (%)');
    ylabel('Response');
    set(gca, 'TickDir', 'out');
    box off;
    hold off;
end
% Initialize vector to store half-max contrast values for each neuron
half_max_contrasts = zeros(size(pref_responses_stat_concat{pre}, 1), 1);

% Calculate half-max contrast for each neuron
for neuron = 1:size(pref_responses_stat_concat{pre}, 1)
    % Extract responses for this neuron across all contrasts
    responses = pref_responses_stat_concat{pre}(neuron, :);
    
    % Use polyfit to find the line of best fit
    p = polyfit(contrast_values, responses, 1);
    
    % Slope and y-intercept from the fit
    slope = p(1);
    intercept = p(2);
    
    % Find the maximum predicted response at 100% contrast
    max_response = slope * 100 + intercept;
    
    % Calculate half of the maximum response
    half_max = max_response / 2;
    
    % Find the contrast value where response equals half_max
    % Using the equation: half_max = slope * contrast + intercept
    % Rearranging: contrast = (half_max - intercept) / slope
    if slope ~= 0
        half_max_contrast = (half_max - intercept) / slope;
        
        % Store the result, but check if it's a reasonable value
        if half_max_contrast > 0 && half_max_contrast <= 100
            half_max_contrasts(neuron) = half_max_contrast;
        else
            % If the value is outside [0, 100], mark as NaN or another special value
            half_max_contrasts(neuron) = NaN;
        end
    else
        % If slope is zero, the half-max contrast is undefined
        half_max_contrasts(neuron) = NaN;
    end
end

% Optional: Update the visualization to show half-max points
figure;
for i = 1:min(5, size(pref_responses_stat_concat{pre}, 1))
    subplot(min(5, size(pref_responses_stat_concat{pre}, 1)), 1, i);
    
    % Plot the actual data points
    plot(contrast_values, pref_responses_stat_concat{pre}(i, :), 'o');
    hold on;
    
    % Calculate fitted line parameters
    p = polyfit(contrast_values, pref_responses_stat_concat{pre}(i, :), 1);
    x_range = [0 100];
    fitted_line = polyval(p, x_range);
    
    % Plot the fitted line
    plot(x_range, fitted_line, '-');
    
    % Calculate and plot half-max point
    max_response = p(1) * 100 + p(2);
    half_max = max_response / 2;
    
    if ~isnan(half_max_contrasts(i))
        % Plot half-max point
        plot(half_max_contrasts(i), half_max, 'rx', 'MarkerSize', 10);
        % Add a horizontal line at half-max response
        plot([0, 100], [half_max, half_max], 'r--');
        % Add a vertical line at half-max contrast
        plot([half_max_contrasts(i), half_max_contrasts(i)], [0, max_response], 'r--');
    end
    
    title(['Neuron ' num2str(i) ', C_{50} = ' num2str(half_max_contrasts(i))]);
    xlabel('Contrast (%)');
    ylabel('Response');
    set(gca, 'TickDir', 'out');
    box off;
    hold off;
end
%% Exploring relationship with OSI added as an independent variable
% Add OSI to the individual correlations and the regression

mySample = red_ind_concat;

% Define the state (1 = stationary)
state_idx = 1;

% Create a figure with 5 subplots (one for each dependent variable, including OSI)
figure('Position', [100 100 600 1200]);

% Loop through dependent variables, now including OSI
dep_vars = {meanLMI, response_slopes, pref_responses_stat_concat{pre}(:,3), OSI_concat{pre}};
dep_var_names = {'Mean LMI', 'Response Slopes', 'Response Amplitude', 'Orientation Selectivity (OSI)'};

for dep_var_idx = 1:length(dep_vars)
    % Calculate subplot position
    subplot(5, 1, dep_var_idx);
    
    % Extract noise correlations for stationary state
    noise_corr = noiseCorr_OG_concat{pre}(state_idx, :);
    
    % Get the current dependent variable
    dep_var = dep_vars{dep_var_idx};
    
    % First filter by mySample indices
    % Make sure the indices are valid
    valid_indices = mySample(mySample <= length(noise_corr) & mySample <= length(dep_var));
    
    % Extract data for valid indices only
    subset_noise_corr = noise_corr(valid_indices);
    subset_dep_var = dep_var(valid_indices)';
    
    % Then remove any NaN values
    valid_data = ~isnan(subset_noise_corr) & ~isnan(subset_dep_var);
    x_data = subset_noise_corr(valid_data);
    y_data = subset_dep_var(valid_data);
    
    % Check if we have enough data points
    if length(x_data) < 2
        warning(['Not enough valid data points for ' dep_var_names{dep_var_idx}]);
        text(0.5, 0.5, 'Insufficient data', 'FontSize', 14, 'HorizontalAlignment', 'center');
        set(gca, 'TickDir', 'out');
        box off;
        xlabel('Noise Correlation (Stationary)');
        ylabel(dep_var_names{dep_var_idx});
        title(['Noise Corr vs. ' dep_var_names{dep_var_idx}], 'FontSize', 10);
        continue;
    end
    
    % Scatter plot
    scatter(x_data, y_data, 50, 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    
    % Compute linear regression
    [p, S] = polyfit(x_data, y_data, 1);
    x_fit = linspace(min(x_data), max(x_data), 100);
    y_fit = polyval(p, x_fit);
    
    % Add regression line
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
    
    % Calculate R^2 and p-value
    [r, p_val] = corrcoef(x_data, y_data);
    r_squared = r(1,2)^2;
    
    % Add R^2 and p-value text
    if p_val(1,2) < 0.001
        p_text = 'p < 0.001';
    else
        p_text = ['p = ' num2str(p_val(1,2), '%.3f')];
    end
    
    % Position the text to better fit a narrower plot - adjust positioning to bottom right
    text(min(x_data) + 0.6 * range(x_data), min(y_data) + 0.15 * range(y_data), ...
        {['R^2 = ' num2str(r_squared, '%.3f')], p_text, ['y = ' num2str(p(1), '%.3f') 'x + ' num2str(p(2), '%.3f')]}, ...
        'FontSize', 9, 'BackgroundColor', [0.9 0.9 0.9, 0.7]);
    
    % Labels and title
    xlabel('Noise Correlation (Stationary)');
    ylabel(dep_var_names{dep_var_idx});
    title(['Relationship Between Noise Correlations and ' dep_var_names{dep_var_idx}]);
    
    % Customize appearance
    set(gca, 'TickDir', 'out');
    box off;
end

% Adjust spacing and use a more compact title
sgtitle('Noise Correlations vs. Neural Measures (SST)', 'FontSize', 12);
set(gcf, 'Color', 'w');
print('-dpdf', 'scattersVsNoiseCorr.pdf', '-bestfit');

%% Update the analysis of stationary_norm_diff_25 to include OSI
mySample = red_ind_concat;
% Extract data (stationary state only, 25% only)
stationary_norm_diff_25 = squeeze(norm_diff(1,1,mySample)); 

% Create a figure for individual correlations with OSI included
figure('Position', [100 100 800 1000]);

% Define independent variables including OSI
raw_ind_vars = {meanLMI(mySample), response_slopes(mySample), pref_responses_stat_concat{pre}(mySample,3), noiseCorr_OG_concat{pre}(1,mySample)', OSI_concat{pre}(mySample)};
ind_var_names = {'Mean LMI', 'Response Slopes', 'Response Amplitude', 'Noise Correlation', 'Orientation Selectivity (OSI)'};

% Z-score the independent variables
ind_vars = cell(size(raw_ind_vars));
for i = 1:length(raw_ind_vars)
    valid_data = ~isnan(raw_ind_vars{i});
    temp_data = raw_ind_vars{i};
    temp_data(valid_data) = zscore(temp_data(valid_data));
    ind_vars{i} = temp_data;
end

% Loop through independent variables for correlation plots
for var_idx = 1:length(ind_vars)
    subplot(3, 2, var_idx);
    
    % Extract data
    x_data = ind_vars{var_idx};
    y_data = stationary_norm_diff_25;
    
    % Remove NaN values
    valid_data = ~isnan(x_data) & ~isnan(y_data);
    x_clean = x_data(valid_data);
    y_clean = y_data(valid_data);
    
    % Check if we have enough data points
    if length(x_clean) < 2
        warning(['Not enough valid data points for ' ind_var_names{var_idx}]);
        text(0.5, 0.5, 'Insufficient data', 'FontSize', 14, 'HorizontalAlignment', 'center');
        continue;
    end
    
    % Scatter plot
    scatter(x_clean, y_clean, 50, 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    
    % Compute linear regression
    [p, S] = polyfit(x_clean, y_clean, 1);
    x_fit = linspace(min(x_clean), max(x_clean), 100);
    y_fit = polyval(p, x_fit);
    
    % Add regression line
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
    
    % Calculate R^2 and p-value
    [r, p_val] = corrcoef(x_clean, y_clean);
    r_squared = r(1,2)^2;
    
    % Add R^2 and p-value text
    if p_val(1,2) < 0.001
        p_text = 'p < 0.001';
    else
        p_text = ['p = ' num2str(p_val(1,2), '%.3f')];
    end
    
    text(min(x_clean) + 0.6 * range(x_clean), min(y_clean) + 0.15 * range(y_clean), ...
        {['R^2 = ' num2str(r_squared, '%.3f')], p_text, ['y = ' num2str(p(1), '%.3f') 'x + ' num2str(p(2), '%.3f')]}, ...
        'FontSize', 9, 'BackgroundColor', [0.9 0.9 0.9, 0.7]);
    
    % Labels and title
    xlabel([ind_var_names{var_idx} ' (z-scored)']);
    ylabel('Normalized Difference (25% contrast)');
    title(['Relationship with ' ind_var_names{var_idx}]);
    
    % Apply your preferred plot settings
    set(gca, 'TickDir', 'out');
    box off;
end

% Adjust spacing between subplots
sgtitle('Correlations with Normalized Difference at 25% Contrast', 'FontSize', 12);
set(gcf, 'Color', 'w');
print('-dpdf', 'scattersVsNormDiff.pdf', '-bestfit');


%% Multiple regression analysis with OSI included

% Define your sample - using the mySample variable from your original code
mySample = red_ind_concat;

% Define independent variables and their names with OSI added
raw_ind_vars = {meanLMI(mySample), response_slopes(mySample), pref_responses_stat_concat{pre}(mySample,3), noiseCorr_OG_concat{pre}(1,mySample)', OSI_concat{pre}(mySample)};
ind_var_names = {'MeanLMI', 'ResponseSlopes', 'ResponseAmplitude', 'NoiseCorrelation', 'OSI'};

% Dependent variable
y = stationary_norm_diff_25; % Your normalized difference at 25% contrast

% Z-score the independent variables
ind_vars = cell(size(raw_ind_vars));
for i = 1:length(raw_ind_vars)
    valid_data = ~isnan(raw_ind_vars{i});
    temp_data = raw_ind_vars{i};
    temp_data(valid_data) = zscore(temp_data(valid_data));
    ind_vars{i} = temp_data;
end

% Prepare data matrix for regression
X = zeros(length(mySample), length(ind_vars));
for i = 1:length(ind_vars)
    X(:, i) = ind_vars{i};
end

% Identify rows with all valid data (no NaN values)
valid_rows = true(length(mySample), 1);
for i = 1:length(ind_vars)
    valid_rows = valid_rows & ~isnan(X(:, i));
end
valid_rows = valid_rows & ~isnan(y);

% Filter data to remove rows with NaN values
X_clean = X(valid_rows, :);
y_clean = y(valid_rows);

% Check if we have enough data for regression
if size(X_clean, 1) < length(ind_vars) + 1
    warning('Not enough valid data points for multiple regression');
    disp(['Only ' num2str(size(X_clean, 1)) ' complete cases available, but need at least ' num2str(length(ind_vars) + 1)]);
else
    % Create a table with appropriate variable names
    tbl = array2table(X_clean, 'VariableNames', ind_var_names);
    tbl.y = y_clean;

    % Fit linear model
    mdl = fitlm(tbl, 'y ~ 1 + MeanLMI + ResponseSlopes + ResponseAmplitude + NoiseCorrelation + OSI');

    % Display model results
    disp(mdl);

    % Create a figure with appropriate size
    figure('Units', 'inches', 'Position', [1 1 7 3]);
    
    % Get standardized coefficients (betas)
    beta_values = mdl.Coefficients.Estimate(2:end); % Skip intercept
    beta_errors = mdl.Coefficients.SE(2:end);
    beta_pvalues = mdl.Coefficients.pValue(2:end);
    coef_names = mdl.CoefficientNames(2:end);
    
    % Create the axes with 1 inch height and width
    ax = axes('Units', 'inches', 'Position', [1.5 1 2 1]);
    
    % Create the bar chart with error bars
    b = bar(beta_values);
    hold on;
    errorbar(1:length(beta_values), beta_values, beta_errors, '.k');
    
    % Add a horizontal line at y=0
    plot([0.5, length(beta_values)+0.5], [0, 0], 'k--');
    
    % Add p-values above each bar
    for i = 1:length(beta_values)
        if beta_pvalues(i) < 0.05
            % Format p-value with appropriate precision
            if beta_pvalues(i) < 0.001
                p_text = 'p<0.001';
            else
                p_text = ['p=' num2str(beta_pvalues(i), '%.3f')];
            end
        else
            p_text = 'n.s.';
        end
        
        % Position text above or below bar depending on bar direction
        if beta_values(i) >= 0
            text_pos = beta_values(i) + beta_errors(i) + 0.05*max(abs(beta_values));
        else
            text_pos = beta_values(i) - beta_errors(i) - 0.15*max(abs(beta_values));
        end
        
        text(i, text_pos, p_text, 'HorizontalAlignment', 'center', 'FontSize', 8);
    end
    
    % Format plot
    set(gca, 'XTick', 1:length(coef_names), 'XTickLabel', coef_names, 'XTickLabelRotation', 45);
    set(gca, 'FontSize', 8);
    ylabel('Standardized Coefficient (Beta)', 'FontSize', 8);
    
    % Apply user preferences
    set(gca, 'TickDir', 'out');
    grid off;
    box off;
    title('Standardized Coefficients (All Variables)', 'FontSize', 9);
    print('-dpdf', 'NoisecorrMixedModel_withOSI_1.pdf', '-bestfit');
end

%% Regression model without amplitude, but with OSI

% Check if we have enough data for regression
if size(X_clean, 1) < length(ind_vars) + 1
    warning('Not enough valid data points for multiple regression');
    disp(['Only ' num2str(size(X_clean, 1)) ' complete cases available, but need at least ' num2str(length(ind_vars) + 1)]);
else
    % Create a table with appropriate variable names
    tbl = array2table(X_clean, 'VariableNames', ind_var_names);
    tbl.y = y_clean;

    % Fit linear model (excluding ResponseAmplitude)
    mdl = fitlm(tbl, 'y ~ 1 + MeanLMI + ResponseSlopes + NoiseCorrelation + OSI');

    % Display model results
    disp(mdl);

    % Create a figure with appropriate size
    figure('Units', 'inches', 'Position', [1 1 6 3]);
    
    % Get standardized coefficients (betas)
    beta_values = mdl.Coefficients.Estimate(2:end); % Skip intercept
    beta_errors = mdl.Coefficients.SE(2:end);
    beta_pvalues = mdl.Coefficients.pValue(2:end);
    coef_names = mdl.CoefficientNames(2:end);
    
    % Create the axes with 1 inch height and width
    ax = axes('Units', 'inches', 'Position', [1.5 1 1.5 1]);
    
    % Create the bar chart with error bars
    b = bar(beta_values);
    hold on;
    errorbar(1:length(beta_values), beta_values, beta_errors, '.k');
    
    % Add a horizontal line at y=0
    plot([0.5, length(beta_values)+0.5], [0, 0], 'k--');
    
    % Add p-values above each bar
    for i = 1:length(beta_values)
        if beta_pvalues(i) < 0.05
            % Format p-value with appropriate precision
            if beta_pvalues(i) < 0.001
                p_text = 'p<0.001';
            else
                p_text = ['p=' num2str(beta_pvalues(i), '%.3f')];
            end
        else
            p_text = 'n.s.';
        end
        
        % Position text above or below bar depending on bar direction
        if beta_values(i) >= 0
            text_pos = beta_values(i) + beta_errors(i) + 0.05*max(abs(beta_values));
        else
            text_pos = beta_values(i) - beta_errors(i) - 0.15*max(abs(beta_values));
        end
        
        text(i, text_pos, p_text, 'HorizontalAlignment', 'center', 'FontSize', 8);
    end
    
    % Format plot
    set(gca, 'XTick', 1:length(coef_names), 'XTickLabel', coef_names, 'XTickLabelRotation', 45);
    set(gca, 'FontSize', 8);
    ylabel('Standardized Coefficient (Beta)', 'FontSize', 8);
    
    % Apply user preferences
    set(gca, 'TickDir', 'out');
    grid off;
    box off;
    title('Standardized Coefficients (Excluding Amplitude)', 'FontSize', 9);
    print('-dpdf', 'NoisecorrMixedModel_withOSI_2.pdf', '-bestfit');
end