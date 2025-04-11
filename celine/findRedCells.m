function high_red_cells = findRedCells(redChImg, mask_cell, percentile)
% findHighRedCells - Identifies cells with red fluorescence above a specified percentile
%
% INPUTS:
%   id              - Index to use in allDays
%   percentile      - Percentile threshold (0-100) to select cells


% Extract red fluorescence for all cells
red_fluor_mask = stackGetTimeCourses(redChImg, mask_cell);

% If match_ind exists in workspace, use it, otherwise use all cells
if exist('match_ind', 'var')
    red_fluor = red_fluor_mask(:, match_ind);
else
    red_fluor = red_fluor_mask;
end

% Calculate mean fluorescence for each cell
mean_red_fluor = mean(red_fluor, 1);

% Calculate the threshold based on the specified percentile
threshold = prctile(mean_red_fluor, percentile);

% Find cells above threshold
cells_above_threshold = find(mean_red_fluor > threshold);

% If we're using matched indices, convert back to original indices
if exist('match_ind', 'var')
    high_red_cells = match_ind(cells_above_threshold);
else
    high_red_cells = cells_above_threshold;
end

% Optional: Visualize the results
figure;
histogram(mean_red_fluor);
hold on;
xline(threshold, 'r--', ['p' num2str(percentile) ' threshold'], 'LineWidth', 2);
xlabel('Mean Red Fluorescence');
ylabel('Number of Cells');
title(['Cells Above ' num2str(percentile) 'th Percentile Red Fluorescence']);
set(gca, 'TickDir', 'out');
grid off;
box off;

fprintf('Found %d cells (%.1f%%) above the %.1f percentile threshold.\n', ...
    length(high_red_cells), 100*length(high_red_cells)/length(mean_red_fluor), percentile);
end