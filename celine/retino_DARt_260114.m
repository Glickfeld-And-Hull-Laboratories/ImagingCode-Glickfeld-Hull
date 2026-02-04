%retino_DARt_260114
figure;
data = [ret_distance_keep_concat{pre}; ret_distance_keep_concat{post}];
group = [ones(length(ret_distance_keep_concat{1}), 1); 2*ones(length(ret_distance_keep_concat{2}), 1)];

boxplot(data, group, 'Labels', {'Pre', 'Post'});
hold on;
scatter(ones(size(ret_distance_keep_concat{pre})), ret_distance_keep_concat{pre}, 20, 'k', 'filled', 'jitter', 'on', 'jitterAmount', 0.1);
scatter(2*ones(size(ret_distance_keep_concat{post})), ret_distance_keep_concat{post}, 20, 'k', 'filled', 'jitter', 'on', 'jitterAmount', 0.1);

set(gca, 'TickDir', 'out');
grid off;
box off;
ylabel('Distance');

figure;
subplot(2,1,1)
hist(ret_distance_keep_concat{pre})
xlim([0 12])
title('ret distance pre')
subplot(2,1,2)
hist(ret_distance_keep_concat{post})
xlim([0 12])
title('ret distance post')
saveas(gcf, sprintf('ret_distance_distribution.pdf'));
%% 
%designating close cells as those with distance <6 on both days
closeCellsPre = find(ret_distance_keep_concat{pre}<6);
closeCellsPost = find(ret_distance_keep_concat{post}<6);
closeCellsBothDays = intersect(closeCellsPre,closeCellsPost);
close_red = intersect(closeCellsBothDays,red_ind_concat);
close_green = intersect(closeCellsBothDays,green_ind_concat);

plotNeuralTimecourse(tc_trial_avrg_stat_concat, tc_trial_avrg_stat_concat, ...
    close_red, close_green,'DayOrder',matchDrx, ...
    'UseDashedLines', [false, true], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'HTP+ ret <= 7.1', 'HTP- ret <= 5'}, ...
    'StimStart', 31);

plotSizeResponse(pref_responses_stat_concat, pref_responses_stat_concat, ...
    close_red, close_green, targetCon,targetSize,'DayOrder',matchDrx, ...
    'UseDashedLines', [false, true], ...  % Dashed lines for the right plot
    'Titles', {'HTP+', 'HTP-'}, ...
    'YLabel', 'dF/F');
sgtitle('Stationary')
saveas(gcf, sprintf('closeCell_stationary_size_response.pdf'));
%%
farCellsPre = find(ret_distance_keep_concat{pre}>6);
farCellsPost = find(ret_distance_keep_concat{post}>6);
farCellsBothDays = intersect(farCellsPre,farCellsPost);
far_red = intersect(farCellsBothDays,red_ind_concat);
far_green = intersect(farCellsBothDays,green_ind_concat);

plotNeuralTimecourse(tc_trial_avrg_stat_concat, tc_trial_avrg_stat_concat, ...
    far_red, far_green,'DayOrder',matchDrx, ...
    'UseDashedLines', [false, true], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'HTP+ ret <= 7.1', 'HTP- ret <= 7.1'}, ...
    'StimStart', 31);
%% look at delta distance

delta_ret_distance = ret_distance_keep_concat{post}-ret_distance_keep_concat{pre};
data = squeeze(norm_diff_concat(1,3,:,:));
num_sizes = size(data,1);

keep = all(data <= 5, 1);
data = data(:, keep);
delta_ret_distance_clean = delta_ret_distance(keep);

y_lim = [min(data(:)), max(data(:))];
x_lim = [-10, 10];

figure;
for i = 1:num_sizes
    subplot(1,num_sizes,i)
    scatter(delta_ret_distance_clean, data(i,:))
    hold on
    
    mdl = fitlm(delta_ret_distance_clean, data(i,:));
    r_squared = mdl.Rsquared.Ordinary;
    p_value = mdl.Coefficients.pValue(2);
    
    x_fit = linspace(min(delta_ret_distance_clean), max(delta_ret_distance_clean), 100)';
    [y_fit, y_ci] = predict(mdl, x_fit);
    
    fill([x_fit; flipud(x_fit)], [y_ci(:,1); flipud(y_ci(:,2))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5)
    
    xlabel('Delta Ret Distance')
    ylabel('Norm Diff')
    title(sprintf('Size %d, R² = %.2f, p = %.2f', i, r_squared, p_value))
    ylim(y_lim)
    set(gca, 'TickDir', 'out')
    grid off
    box off
    hold off
end

set(gcf, 'PaperOrientation', 'landscape')
set(gcf, 'PaperPosition', [0.25 2 10.5 4])
print('normDiff_vsRet', '-dpdf')
%% plot responses for cells with little change in ret distance
stableRet = find(delta_ret_distance<2.5);

stableRed = intersect(stableRet,red_ind_concat);
stableGreen = intersect(stableRet,green_ind_concat);

plotNeuralTimecourse(tc_trial_avrg_stat_concat, tc_trial_avrg_stat_concat, ...
    stableRed, stableGreen,'DayOrder',matchDrx, ...
    'UseDashedLines', [false, true], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'HTP+ ret <= 7.1', 'HTP- ret <= 5'}, ...
    'StimStart', 31);

plotSizeResponse(pref_responses_stat_concat, pref_responses_stat_concat, ...
    stableRed, stableGreen, targetCon,targetSize,'DayOrder',matchDrx, ...
    'UseDashedLines', [false, true], ...  % Dashed lines for the right plot
    'Titles', {'HTP+', 'HTP-'}, ...
    'YLabel', 'dF/F');
sgtitle('Stationary')
saveas(gcf, sprintf('stableCell_stationary_size_response.pdf'));