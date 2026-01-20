%retino_DARt_260114
figure;
data = [ret_distance_keep_concat{1}; ret_distance_keep_concat{2}];
group = [ones(length(ret_distance_keep_concat{1}), 1); 2*ones(length(ret_distance_keep_concat{2}), 1)];

boxplot(data, group, 'Labels', {'Pre', 'Post'});
hold on;
scatter(ones(size(ret_distance_keep_concat{1})), ret_distance_keep_concat{1}, 20, 'k', 'filled', 'jitter', 'on', 'jitterAmount', 0.1);
scatter(2*ones(size(ret_distance_keep_concat{2})), ret_distance_keep_concat{2}, 20, 'k', 'filled', 'jitter', 'on', 'jitterAmount', 0.1);

set(gca, 'TickDir', 'out');
grid off;
box off;
ylabel('Distance');

%% 
%designating close cells as those with distance <8 on both days
closeCellsPre = find(ret_distance_keep_concat{pre}<8);
closeCellsPost = find(ret_distance_keep_concat{post}<8);
closeCellsBothDays = intersect(closeCellsPre,closeCellsPost);
close_red = intersect(closeCellsBothDays,red_ind_concat);
close_green = intersect(closeCellsBothDays,green_ind_concat);

plotNeuralTimecourse(tc_trial_avrg_stat_concat, tc_trial_avrg_stat_concat, ...
    close_red, close_green,'DayOrder',matchDrx, ...
    'UseDashedLines', [false, true], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'HTP+ ret <= 7.1', 'HTP- ret <= 7.1'}, ...
    'StimStart', 31);
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

% Remove neurons with any value > 5
keep = all(data <= 5, 1);
data = data(:, keep);
delta_ret_distance_clean = delta_ret_distance(keep);

y_lim = [min(data(:)), max(data(:))];
%y_lim=[]
x_lim=[-20, 20]
figure;
for i = 1:num_sizes
    subplot(1,num_sizes,i)
    scatter(delta_ret_distance_clean, data(i,:))
    xlabel('Delta Ret Distance')
    ylabel('Norm Diff')
    title(['Size ' num2str(i)])
    ylim(y_lim)
    xlim(x_lim)
    set(gca, 'TickDir', 'out')
    grid off
    box off
end

set(gcf, 'PaperOrientation', 'landscape')
set(gcf, 'PaperPosition', [0.25 2 10.5 4])
print('normDiff_vsRet', '-dpdf')
