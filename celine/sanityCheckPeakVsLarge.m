function sanityCheckPeakVsLarge(id, h_largeVsPeak_keep, prefDir_keep, stat_resp_keep, ...
    data_dfof_trial_keep, tCon_match, tSize_match, tDir_match, ...
    RIx, nTrials, dirs, cons, sizes, resp_win, varargin)

% Parse optional inputs
p = inputParser;
addParameter(p, 'nExamples', 3, @isnumeric);
parse(p, varargin{:});
n_examples = p.Results.nExamples;

% Get cells that pass and fail the test
pass_cells = find(h_largeVsPeak_keep{id} == 1);
fail_cells = find(h_largeVsPeak_keep{id} == 0);

% Select examples
n_pass = min(n_examples, length(pass_cells));
n_fail = min(n_examples, length(fail_cells));

if n_pass == 0 && n_fail == 0
    fprintf('No cells found for day %d\n', id);
    return;
end

% Extract trial info
tCon = tCon_match{id}(:, 1:nTrials(id));
tSize = tSize_match{id}(:, 1:nTrials(id));
tDir = tDir_match{id}(:, 1:nTrials(id));
stat_inds = find(~RIx{id});

figure;

subplot_idx = 1;

% Plot pass cells
for i = 1:n_pass
    iCell = pass_cells(i);
    
    % Get preferred direction and find peak size
    pref_dir = prefDir_keep{id}(iCell);
    resp_by_size = squeeze(stat_resp_keep{id}(iCell, pref_dir, end, :));
    [~, peakSize] = max(resp_by_size);
    
    % Get trial indices
    ind_dir = find(tDir == dirs(pref_dir));
    ind_con = find(tCon == cons(end));
    ind_temp = intersect(ind_dir, ind_con);
    
    % Get responses
    ind_large = intersect(intersect(ind_temp, find(tSize == sizes(end))), stat_inds);
    ind_peak = intersect(intersect(ind_temp, find(tSize == sizes(peakSize))), stat_inds);
    
    if ~isempty(ind_large) && ~isempty(ind_peak)
        resp_large = squeeze(nanmean(data_dfof_trial_keep{id}(resp_win, ind_large, iCell), 1));
        resp_peak = squeeze(nanmean(data_dfof_trial_keep{id}(resp_win, ind_peak, iCell), 1));
        
        subplot(2, max(n_pass, n_fail), subplot_idx);
        
        plot(1, resp_peak, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        hold on;
        plot(2, resp_large, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Peak', 'Large'});
        set(gca, 'TickDir', 'out');
        grid off;
        box off;
        ylabel('\Delta F/F');
        title(sprintf('Pass Cell %d', i), 'Color', 'b');
        xlim([0.5 2.5]);
        
        fprintf('Pass Cell %d: Peak=%.3f±%.3f, Large=%.3f±%.3f\n', ...
            i, mean(resp_peak), std(resp_peak), mean(resp_large), std(resp_large));
        
        subplot_idx = subplot_idx + 1;
    end
end

% Plot fail cells
for i = 1:n_fail
    iCell = fail_cells(i);
    
    % Get preferred direction and find peak size
    pref_dir = prefDir_keep{id}(iCell);
    resp_by_size = squeeze(stat_resp_keep{id}(iCell, pref_dir, end, :));
    [~, peakSize] = max(resp_by_size);
    
    % Get trial indices
    ind_dir = find(tDir == dirs(pref_dir));
    ind_con = find(tCon == cons(end));
    ind_temp = intersect(ind_dir, ind_con);
    
    % Get responses
    ind_large = intersect(intersect(ind_temp, find(tSize == sizes(end))), stat_inds);
    ind_peak = intersect(intersect(ind_temp, find(tSize == sizes(peakSize))), stat_inds);
    
    if ~isempty(ind_large) && ~isempty(ind_peak)
        resp_large = squeeze(nanmean(data_dfof_trial_keep{id}(resp_win, ind_large, iCell), 1));
        resp_peak = squeeze(nanmean(data_dfof_trial_keep{id}(resp_win, ind_peak, iCell), 1));
        
        subplot(2, max(n_pass, n_fail), max(n_pass, n_fail) + i);
        
        plot(1, resp_peak, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        hold on;
        plot(2, resp_large, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Peak', 'Large'});
        set(gca, 'TickDir', 'out');
        grid off;
        box off;
        ylabel('\Delta F/F');
        title(sprintf('Fail Cell %d', i), 'Color', 'r');
        xlim([0.5 2.5]);
        
        fprintf('Fail Cell %d: Peak=%.3f±%.3f, Large=%.3f±%.3f\n', ...
            i, mean(resp_peak), std(resp_peak), mean(resp_large), std(resp_large));
    end
end

sgtitle(sprintf('Peak vs Large Size Responses: Day %d Sanity Check', id));

fprintf('\nSanity Check Summary for Day %d:\n', id);
fprintf('Pass cells should show Peak > Large responses\n');
fprintf('Fail cells should show Peak ≈ Large responses\n');

end