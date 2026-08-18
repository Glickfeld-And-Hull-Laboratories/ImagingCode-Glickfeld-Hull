function plotPeakVsLargeBoxplot(id, h_largeVsPeak_keep, prefDir_keep, stat_resp_keep, ...
    data_dfof_trial_keep, tCon_match, tSize_match, tDir_match, ...
    RIx, nTrials, dirs, cons, sizes, resp_win, varargin)

% Parse optional inputs
p = inputParser;
addParameter(p, 'nCells', 6, @isnumeric);
addParameter(p, 'passColor', [0.2 0.6 0.8], @isnumeric);
addParameter(p, 'failColor', [0.8 0.2 0.2], @isnumeric);
parse(p, varargin{:});
n_cells_plot = p.Results.nCells;
pass_color = p.Results.passColor;
fail_color = p.Results.failColor;

% Get cells that pass and fail the test
pass_cells = find(h_largeVsPeak_keep{id} == 1);
fail_cells = find(h_largeVsPeak_keep{id} == 0);

% Select subset for plotting
n_plot_pass = min(n_cells_plot, length(pass_cells));
n_plot_fail = min(n_cells_plot, length(fail_cells));

% Extract trial info
tCon = tCon_match{id}(:, 1:nTrials(id));
tSize = tSize_match{id}(:, 1:nTrials(id));
tDir = tDir_match{id}(:, 1:nTrials(id));
stat_inds = find(~RIx{id});

% Create figure
figure;

pos = 1;
positions = [];
labels_plot = {};
all_cells = [pass_cells(1:n_plot_pass), fail_cells(1:n_plot_fail)];
cell_types = [ones(1, n_plot_pass), zeros(1, n_plot_fail)]; % 1=pass, 0=fail

% First pass: collect all data to avoid multiple boxplot calls
data_matrix = [];
group_matrix = [];
color_vector = [];

for idx = 1:length(all_cells)
    iCell = all_cells(idx);
    is_pass = cell_types(idx);
    color = pass_color * is_pass + fail_color * (1 - is_pass);
    
    % Get preferred direction and find peak size
    pref_dir = prefDir_keep{id}(iCell);
    resp_by_size = squeeze(stat_resp_keep{id}(iCell, pref_dir, end, :));
    [~, peakSize] = max(resp_by_size);
    
    % Get trial indices for preferred direction and highest contrast
    ind_dir = find(tDir == dirs(pref_dir));
    ind_con = find(tCon == cons(end));
    ind_temp = intersect(ind_dir, ind_con);
    
    % Get responses for peak size and large size
    ind_large = intersect(intersect(ind_temp, find(tSize == sizes(end))), stat_inds);
    ind_peak = intersect(intersect(ind_temp, find(tSize == sizes(peakSize))), stat_inds);
    
    if ~isempty(ind_large) && ~isempty(ind_peak)
        resp_large = squeeze(nanmean(data_dfof_trial_keep{id}(resp_win, ind_large, iCell), 1));
        resp_peak = squeeze(nanmean(data_dfof_trial_keep{id}(resp_win, ind_peak, iCell), 1));
        
        % Store data in matrix format for single boxplot call
        max_trials = max(length(resp_peak(:)), length(resp_large(:)));
        
        % Pad with NaNs to make equal length
        peak_padded = [resp_peak(:); NaN(max_trials - length(resp_peak(:)), 1)];
        large_padded = [resp_large(:); NaN(max_trials - length(resp_large(:)), 1)];
        
        data_matrix = [data_matrix, peak_padded, large_padded];
        
        % Track positions and labels
        positions = [positions, pos, pos+1];
        labels_plot{end+1} = sprintf('Peak\nCell %d', idx);
        labels_plot{end+1} = sprintf('Large\nCell %d', idx);
        color_vector = [color_vector; color; color];
        pos = pos + 3; % Space between cell groups
    end
end

% Single boxplot call with matrix input
boxplot(data_matrix, 'positions', positions, 'symbol', '', 'widths', 0.6);

% Color the boxes after creation
h_boxes = findobj(gca, 'Tag', 'Box');
for i = 1:min(length(h_boxes), size(color_vector, 1))
    set(h_boxes(end-i+1), 'Color', color_vector(i, :), 'LineWidth', 1.5);
end

set(gca, 'TickDir', 'out');
grid off;
box off;
set(gca, 'XTick', positions, 'XTickLabel', labels_plot);
xtickangle(45);
ylabel('\Delta F/F');
title(sprintf('Peak vs Large Size Responses: Day %d\n(Blue = No Suppression, Red = Has Suppression)', id));

% Add legend
hold on;
h1 = plot(NaN, NaN, 's', 'MarkerFaceColor', pass_color, 'MarkerEdgeColor', pass_color);
h2 = plot(NaN, NaN, 's', 'MarkerFaceColor', fail_color, 'MarkerEdgeColor', fail_color);
legend([h1, h2], {'No Suppression', 'Has Suppression'}, 'Location', 'best');

% Print summary
fprintf('Day %d Summary:\n', id);
fprintf(' Cells with no suppression: %d/%d plotted\n', n_plot_pass, length(pass_cells));
fprintf(' Cells with suppression: %d/%d plotted\n', n_plot_fail, length(fail_cells));

end