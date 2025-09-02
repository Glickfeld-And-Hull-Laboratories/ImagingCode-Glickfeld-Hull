function plotPeakVsLargeBoxplot(id, h_largeVsPeak_keep, prefDir_keep, stat_resp_keep, ...
                                data_dfof_trial_keep, tCon_match, tSize_match, tDir_match, ...
                                RIx, nTrials, dirs, cons, sizes, resp_win, varargin)
% plotPeakVsLargeBoxplot - Plot peak vs large size responses for cells that pass/fail t-test
%
% Inputs:
%   id - Day index to plot
%   h_largeVsPeak_keep - Cell array of hypothesis test results
%   prefDir_keep - Cell array of preferred directions
%   stat_resp_keep - Cell array of stationary responses
%   data_dfof_trial_keep - Cell array of trial data
%   tCon_match, tSize_match, tDir_match - Trial condition arrays
%   RIx - Running index for each day
%   nTrials - Number of trials per day
%   dirs, cons, sizes - Stimulus parameter vectors
%   resp_win - Response window indices
%   varargin - Optional parameters:
%       'nCells' - Number of cells to plot per group (default: 10)
%       'passColor' - Color for passing cells (default: [0.2 0.6 0.8])
%       'failColor' - Color for failing cells (default: [0.8 0.2 0.2])

% Parse optional inputs
p = inputParser;
addParameter(p, 'nCells', 10, @isnumeric);
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
pass_subset = pass_cells(1:n_plot_pass);
fail_subset = fail_cells(1:n_plot_fail);

% Extract trial info
tCon = tCon_match{id}(:, 1:nTrials(id));
tSize = tSize_match{id}(:, 1:nTrials(id));
tDir = tDir_match{id}(:, 1:nTrials(id));
stat_inds = find(~RIx{id});

% Collect responses for plotting
all_responses = [];
all_labels = {};

for cell_type = 1:2
    if cell_type == 1
        cells_to_plot = pass_subset;
        type_name = 'Pass';
    else
        cells_to_plot = fail_subset;
        type_name = 'Fail';
    end
    
    for i = 1:length(cells_to_plot)
        iCell = cells_to_plot(i);
        
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
            
            % Ensure column vectors
            resp_large = resp_large(:);
            resp_peak = resp_peak(:);
            
            % Add to data arrays
            all_responses = [all_responses; resp_peak; resp_large];
            all_labels = [all_labels; repmat({sprintf('Peak_%s_C%d', type_name, i)}, length(resp_peak), 1); ...
                         repmat({sprintf('Large_%s_C%d', type_name, i)}, length(resp_large), 1)];
        end
    end
end

% Create figure
figure;
positions = [];
labels_plot = {};

pos = 1;
for cell_type = 1:2
    if cell_type == 1
        cells_to_plot = pass_subset;
        type_name = 'Pass';
        color = pass_color;
    else
        cells_to_plot = fail_subset;
        type_name = 'Fail';
        color = fail_color;
    end
    
    for i = 1:length(cells_to_plot)
        % Peak response
        peak_idx = strcmp(all_labels, sprintf('Peak_%s_C%d', type_name, i));
        if any(peak_idx)
            boxplot(all_responses(peak_idx), 'positions', pos, 'colors', color, 'symbol', '');
            positions = [positions, pos];
            labels_plot{end+1} = sprintf('Peak\n%s C%d', type_name, i);
            pos = pos + 1;
        end
        
        % Large response
        large_idx = strcmp(all_labels, sprintf('Large_%s_C%d', type_name, i));
        if any(large_idx)
            hold on;
            boxplot(all_responses(large_idx), 'positions', pos, 'colors', color, 'symbol', '');
            positions = [positions, pos];
            labels_plot{end+1} = sprintf('Large\n%s C%d', type_name, i);
            pos = pos + 2; % Extra space between cells
        end
    end
    pos = pos + 3; % Extra space between pass/fail groups
end

set(gca, 'TickDir', 'out');
grid off;
box off;
set(gca, 'XTick', positions, 'XTickLabel', labels_plot);
xtickangle(45);
ylabel('\Delta F/F');
title(sprintf('Peak vs Large Size Responses: Day %d (Pass: %d cells, Fail: %d cells)', ...
              id, n_plot_pass, n_plot_fail));

% Add legend
h1 = plot(NaN, NaN, 's', 'MarkerFaceColor', pass_color, 'MarkerEdgeColor', pass_color);
h2 = plot(NaN, NaN, 's', 'MarkerFaceColor', fail_color, 'MarkerEdgeColor', fail_color);
legend([h1, h2], {'Pass t-test', 'Fail t-test'}, 'Location', 'best');

% Print summary
fprintf('Day %d Summary:\n', id);
fprintf('  Cells passing t-test: %d/%d plotted\n', n_plot_pass, length(pass_cells));
fprintf('  Cells failing t-test: %d/%d plotted\n', n_plot_fail, length(fail_cells));

end