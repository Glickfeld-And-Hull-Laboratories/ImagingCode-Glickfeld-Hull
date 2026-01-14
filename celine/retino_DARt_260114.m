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