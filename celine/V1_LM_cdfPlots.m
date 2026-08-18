%load('V1_noiseCorrCompare.mat')
%load('LM_noiseCorrCompare.mat')

V1_mean_pref_red = mean_pref_resp_V1(red_ind_V1);
V1_noiseCorr_red = noiseCorr_concat_V1(red_ind_V1)';
V1_normDiff_red = squeeze(norm_diff_lowConFF(red_ind_V1));
LM_mean_pref_red = mean_pref_resp_LM(red_ind_LM);
LM_noiseCorr_red = noiseCorr_concat_LM(red_ind_LM)';
LM_normDiff_red = squeeze(norm_diff_lowConFF_LM(red_ind_LM));

% Define custom colors
V1_color = '#E09F3E';
LM_color = '#5B9279';

figure;
[f1, xi1] = ksdensity(V1_noiseCorr_red);
plot(xi1, f1, 'Color', V1_color, 'LineWidth', 2);
xlabel('Noise Correlation');
ylabel('Probability Density');
hold on;
[f2, xi2] = ksdensity(LM_noiseCorr_red);
plot(xi2, f2, 'Color', LM_color, 'LineWidth', 2);
legend('V1', 'LM', 'Location', 'best');
set(gca, 'TickDir', 'out');
grid off;
box off;
sgtitle('SST cells')

% Define optimal bin boundaries (6 bins, splitting original bin 2)
bin_edges = [-0.05, 0, 0.025, 0.05, 0.1, 0.2, 0.5];

% Define custom colors
V1_color = '#E09F3E';
LM_color = '#5B9279';
V1_color_rgb = [224 159 62]/255;
LM_color_rgb = [91 146 121]/255;

figure('Position', [100, 100, 800, 600])

% Plot CDFs for each bin comparing V1 and LM
for i = 1:6
    subplot(2, 3, i);
    
    % Create masks for current bin for both datasets
    if i == 6
        mask_V1 = V1_mean_pref_red >= bin_edges(i) & V1_mean_pref_red <= bin_edges(i+1);
        mask_LM = LM_mean_pref_red >= bin_edges(i) & LM_mean_pref_red <= bin_edges(i+1);
    else
        mask_V1 = V1_mean_pref_red >= bin_edges(i) & V1_mean_pref_red < bin_edges(i+1);
        mask_LM = LM_mean_pref_red >= bin_edges(i) & LM_mean_pref_red < bin_edges(i+1);
    end
    
    % Plot V1 CDF
    data_V1 = V1_noiseCorr_red(mask_V1);
    if ~isempty(data_V1)
        [f_V1, x_V1] = ecdf(data_V1);
        plot(x_V1, f_V1, 'Color', V1_color, 'LineWidth', 1.5);
        hold on;
        text(0.95, 0.15, sprintf('V1: %d', sum(mask_V1)), 'Units', 'normalized', 'Color', V1_color, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
        text(0.95, 0.05, sprintf('LM: %d', sum(mask_LM)), 'Units', 'normalized', 'Color', LM_color, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    end
    
    % Plot LM CDF
    data_LM = LM_noiseCorr_red(mask_LM);
    if ~isempty(data_LM)
        [f_LM, x_LM] = ecdf(data_LM);
        plot(x_LM, f_LM, 'Color', LM_color, 'LineWidth', 1.5);
    end
    
    title(sprintf('Bin %d: [%.3f, %.3f]', i, bin_edges(i), bin_edges(i+1)));
    if i == 1
        legend('V1', 'LM', 'Location', 'best');
    end
    set(gca, 'TickDir', 'out');
    grid off; box off;
    hold off;
end
sgtitle('SST cells')
print('SST_CDF_plots.pdf', '-dpdf', '-bestfit')

%%
% Create grouped bar chart for SST cells
figure('Position', [100, 100, 600, 400]);
V1_means = [];
V1_sems = [];
LM_means = [];
LM_sems = [];

for i = 1:6
    if i == 6
        mask_V1 = V1_mean_pref_red >= bin_edges(i) & V1_mean_pref_red <= bin_edges(i+1);
        mask_LM = LM_mean_pref_red >= bin_edges(i) & LM_mean_pref_red <= bin_edges(i+1);
    else
        mask_V1 = V1_mean_pref_red >= bin_edges(i) & V1_mean_pref_red < bin_edges(i+1);
        mask_LM = LM_mean_pref_red >= bin_edges(i) & LM_mean_pref_red < bin_edges(i+1);
    end
    
    data_V1 = V1_noiseCorr_red(mask_V1);
    data_LM = LM_noiseCorr_red(mask_LM);
    
    V1_means(i) = mean(data_V1);
    V1_sems(i) = std(data_V1) / sqrt(length(data_V1));
    LM_means(i) = mean(data_LM);
    LM_sems(i) = std(data_LM) / sqrt(length(data_LM));
end

x = 1:6;
bar_data = [V1_means; LM_means]';
b = bar(x, bar_data, 'grouped');
b(1).FaceColor = V1_color_rgb;
b(2).FaceColor = LM_color_rgb;

hold on;
errorbar(x - 0.15, V1_means, V1_sems, 'k.', 'LineWidth', 1);
errorbar(x + 0.15, LM_means, LM_sems, 'k.', 'LineWidth', 1);
bin_labels = {sprintf('[%.3f,%.3f]', bin_edges(1), bin_edges(2)), ...
              sprintf('[%.3f,%.3f]', bin_edges(2), bin_edges(3)), ...
              sprintf('[%.3f,%.3f]', bin_edges(3), bin_edges(4)), ...
              sprintf('[%.3f,%.3f]', bin_edges(4), bin_edges(5)), ...
              sprintf('[%.3f,%.3f)]', bin_edges(5), bin_edges(6)), ...
              sprintf('[%.3f,%.3f]', bin_edges(6), bin_edges(7))};
set(gca, 'XTickLabel', bin_labels);
xtickangle(45);

xlabel('Mean Response Range');
ylabel('Mean Noise Correlation');
title('SST cells');
legend('V1', 'LM', 'Location', 'best');
set(gca, 'TickDir', 'out');
grid off; box off;
print('SST_bar_chart.pdf', '-dpdf', '-bestfit')
%%
% Define the green subsets
total_cells_V1 = length(mean_pref_resp_V1);
green_ind_V1 = setdiff(1:total_cells_V1, red_ind_V1);

total_cells_LM = length(mean_pref_resp_LM);
green_ind_LM = setdiff(1:total_cells_LM, red_ind_LM);

V1_mean_pref_green = mean_pref_resp_V1(green_ind_V1);
V1_noiseCorr_green = noiseCorr_concat_V1(green_ind_V1)';
V1_normDiff_green = mean_norm_diff_V1(green_ind_V1);
LM_mean_pref_green = mean_pref_resp_LM(green_ind_LM);
LM_noiseCorr_green = noiseCorr_concat_LM(green_ind_LM)';
LM_normDiff_green = mean_norm_diff_LM(green_ind_LM);

figure;
[f1, xi1] = ksdensity(V1_noiseCorr_green);
plot(xi1, f1, 'Color', V1_color, 'LineWidth', 2);
xlabel('Noise Correlation');
ylabel('Probability Density');
hold on;
[f2, xi2] = ksdensity(LM_noiseCorr_green);
plot(xi2, f2, 'Color', LM_color, 'LineWidth', 2);
legend('V1', 'LM', 'Location', 'best');
set(gca, 'TickDir', 'out');
grid off;
box off;
sgtitle('Pyr cells')

figure('Position', [100, 100, 800, 600])

% Plot CDFs for each bin comparing V1 and LM
for i = 1:6
    subplot(2, 3, i);
    
    % Create masks for current bin for both datasets
    if i == 6
        mask_V1 = V1_mean_pref_green >= bin_edges(i) & V1_mean_pref_green <= bin_edges(i+1);
        mask_LM = LM_mean_pref_green >= bin_edges(i) & LM_mean_pref_green <= bin_edges(i+1);
    else
        mask_V1 = V1_mean_pref_green >= bin_edges(i) & V1_mean_pref_green < bin_edges(i+1);
        mask_LM = LM_mean_pref_green >= bin_edges(i) & LM_mean_pref_green < bin_edges(i+1);
    end
    
    % Plot V1 CDF
    data_V1 = V1_noiseCorr_green(mask_V1);
    if ~isempty(data_V1)
        [f_V1, x_V1] = ecdf(data_V1);
        plot(x_V1, f_V1, 'Color', V1_color, 'LineWidth', 1.5);
        hold on;
        text(0.95, 0.15, sprintf('V1: %d', sum(mask_V1)), 'Units', 'normalized', 'Color', V1_color, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
        text(0.95, 0.05, sprintf('LM: %d', sum(mask_LM)), 'Units', 'normalized', 'Color', LM_color, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    end
    
    % Plot LM CDF
    data_LM = LM_noiseCorr_green(mask_LM);
    if ~isempty(data_LM)
        [f_LM, x_LM] = ecdf(data_LM);
        plot(x_LM, f_LM, 'Color', LM_color, 'LineWidth', 1.5);
    end
    xlim([-.2 1])
    title(sprintf('Bin %d: [%.3f, %.3f]', i, bin_edges(i), bin_edges(i+1)));
    if i == 1
        legend('V1', 'LM', 'Location', 'best');
    end
    set(gca, 'TickDir', 'out');
    grid off; box off;
    hold off;
end
sgtitle('Pyr cells')
print('Pyr_CDF_plots.pdf', '-dpdf', '-bestfit')


%%
% Create grouped bar chart for Pyr cells
figure('Position', [100, 100, 600, 400]);
V1_means = [];
V1_sems = [];
LM_means = [];
LM_sems = [];

for i = 1:6
    if i == 6
        mask_V1 = V1_mean_pref_green >= bin_edges(i) & V1_mean_pref_green <= bin_edges(i+1);
        mask_LM = LM_mean_pref_green >= bin_edges(i) & LM_mean_pref_green <= bin_edges(i+1);
    else
        mask_V1 = V1_mean_pref_green >= bin_edges(i) & V1_mean_pref_green < bin_edges(i+1);
        mask_LM = LM_mean_pref_green >= bin_edges(i) & LM_mean_pref_green < bin_edges(i+1);
    end
    
    data_V1 = V1_noiseCorr_green(mask_V1);
    data_LM = LM_noiseCorr_green(mask_LM);
    
    V1_means(i) = mean(data_V1);
    V1_sems(i) = std(data_V1) / sqrt(length(data_V1));
    LM_means(i) = mean(data_LM);
    LM_sems(i) = std(data_LM) / sqrt(length(data_LM));
end

x = 1:6;
bar_data = [V1_means; LM_means]';
b = bar(x, bar_data, 'grouped');
b(1).FaceColor = V1_color_rgb;
b(2).FaceColor = LM_color_rgb;

hold on;
errorbar(x - 0.15, V1_means, V1_sems, 'k.', 'LineWidth', 1);
errorbar(x + 0.15, LM_means, LM_sems, 'k.', 'LineWidth', 1);

bin_labels = {sprintf('[%.3f,%.3f]', bin_edges(1), bin_edges(2)), ...
              sprintf('[%.3f,%.3f]', bin_edges(2), bin_edges(3)), ...
              sprintf('[%.3f,%.3f]', bin_edges(3), bin_edges(4)), ...
              sprintf('[%.3f,%.3f]', bin_edges(4), bin_edges(5)), ...
              sprintf('[%.3f,%.3f)]', bin_edges(5), bin_edges(6)), ...
              sprintf('[%.3f,%.3f]', bin_edges(6), bin_edges(7))};
set(gca, 'XTickLabel', bin_labels);
xtickangle(45);

xlabel('Mean Response Range');
ylabel('Mean Noise Correlation');
title('Pyr cells');
legend('V1', 'LM', 'Location', 'best');
set(gca, 'TickDir', 'out');
grid off; box off;
print('Pyr_bar_chart.pdf', '-dpdf', '-bestfit')
%% 
% Extract the data
noiseCorr_data = noiseCorr_concat{pre}(1,:);
norm_diff_data = squeeze(norm_diff_concat(1,:,:,:));

% Define colors
LM_color_rgb = [91 146 121]/255;

% SST cells figure
figure('Position', [100, 100, 900, 900]);

for contrast = 1:3
    for size = 1:3
        subplot(3, 3, (contrast-1)*3 + size);
        
        % Extract normalized difference for this contrast/size combination for SST cells
        norm_diff_subset = squeeze(norm_diff_data(contrast, size, red_ind_concat));
        noiseCorr_subset = noiseCorr_data(red_ind_concat);
        
        % Create scatterplot
        scatter(noiseCorr_subset, norm_diff_subset, 30, LM_color_rgb, 'filled', 'MarkerFaceAlpha', 0.6);
        
        xlabel('Noise Correlation');
        ylabel('Normalized Difference');
        title(sprintf('Contrast %d, Size %d', contrast, size));
        
        set(gca, 'TickDir', 'out');
        grid off;
        box off;
    end
end
sgtitle('SST cells');

% Pyr cells figure
figure('Position', [100, 100, 900, 900]);

for contrast = 1:3
    for size = 1:3
        subplot(3, 3, (contrast-1)*3 + size);
        
        % Extract normalized difference for this contrast/size combination for Pyr cells
        norm_diff_subset = squeeze(norm_diff_data(contrast, size, green_ind_concat));
        noiseCorr_subset = noiseCorr_data(green_ind_concat);
        
        % Create scatterplot
        scatter(noiseCorr_subset, norm_diff_subset, 30, LM_color_rgb, 'filled', 'MarkerFaceAlpha', 0.6);
        
        xlabel('Noise Correlation');
        ylabel('Normalized Difference');
        title(sprintf('Contrast %d, Size %d', contrast, size));
        
        set(gca, 'TickDir', 'out');
        grid off;
        box off;
    end
end
sgtitle('Pyr cells');