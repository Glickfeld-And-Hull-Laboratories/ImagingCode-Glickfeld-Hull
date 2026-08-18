function plotSizeResponse_contrastOverlay_byCondition_singleDay(statData, locData, cell_indices1, cell_indices2, contrasts, sizes, runningByCondition, varargin)
% PLOTSIZERESPONSE_CONTRASTOVERLAY_BYCONDITION_SINGLEDAY  Size tuning overlay comparing stationary vs running
%   Shows all contrast conditions overlaid, comparing stationary vs running for condition-matched cells
%
%   plotSizeResponse_contrastOverlay_byCondition_singleDay(statData, locData, cell_indices1, cell_indices2, contrasts, sizes, runningByCondition)
%   plotSizeResponse_contrastOverlay_byCondition_singleDay(..., 'Name', Value)
%
% Inputs:
%   statData - Stationary response data (nCells x nContrasts x nSizes)
%   locData - Running response data (nCells x nContrasts x nSizes)
%   cell_indices1 - Indices for population 1 (e.g., red cells)
%   cell_indices2 - Indices for population 2 (e.g., green cells)
%   contrasts - Array of contrast values
%   sizes - Array of size values
%   runningByCondition - nCells x nContrasts x nSizes logical matrix
%
% Optional parameters:
%   'Titles'      - {title1, title2}  Default: {'Population 1', 'Population 2'}
%   'YLabel'      - y-axis label      Default: 'dF/F'
%   'XLabel'      - x-axis label      Default: 'Size (deg)'
%   'XLim'        - [xmin xmax]       Default: auto from sizes
%   'FigureSize'  - [w h] inches      Default: [8, 4]
%   'ShowLegend'  - logical           Default: true

p = inputParser;
addRequired(p, 'statData');
addRequired(p, 'locData');
addRequired(p, 'cell_indices1');
addRequired(p, 'cell_indices2');
addRequired(p, 'contrasts');
addRequired(p, 'sizes');
addRequired(p, 'runningByCondition');
addParameter(p, 'Titles',     {'Population 1', 'Population 2'});
addParameter(p, 'YLabel',     'dF/F');
addParameter(p, 'XLabel',     'Size (deg)');
addParameter(p, 'XLim',       []);
addParameter(p, 'FigureSize', [8, 4]);
addParameter(p, 'ShowLegend', true);
parse(p, statData, locData, cell_indices1, cell_indices2, contrasts, sizes, runningByCondition, varargin{:});
r = p.Results;

nCon  = length(contrasts);
nSize = length(sizes);

if isempty(r.XLim), r.XLim = [min(sizes)*0.8, max(sizes)*1.2]; end

% Create grayscale colormap: lightest to darkest as contrast increases
gray_values = linspace(0.7, 0.0, nCon);
colors = repmat(gray_values', 1, 3);

% Compute condition-specific mean ± SE for each population
avg_stat1 = nan(nCon, nSize);  se_stat1 = nan(nCon, nSize);
avg_loc1  = nan(nCon, nSize);  se_loc1  = nan(nCon, nSize);
avg_stat2 = nan(nCon, nSize);  se_stat2 = nan(nCon, nSize);
avg_loc2  = nan(nCon, nSize);  se_loc2  = nan(nCon, nSize);

n_cells1 = nan(nCon, nSize);
n_cells2 = nan(nCon, nSize);

for iCon = 1:nCon
    for iSize = 1:nSize
        % Find cells that have running data for this specific condition
        cells_with_running = find(runningByCondition(:, iCon, iSize));
        
        % Population 1
        cells1_this_condition = intersect(cell_indices1, cells_with_running);
        if ~isempty(cells1_this_condition)
            stat_vals1 = statData(cells1_this_condition, iCon, iSize);
            loc_vals1 = locData(cells1_this_condition, iCon, iSize);
            
            avg_stat1(iCon, iSize) = mean(stat_vals1, 'omitnan');
            avg_loc1(iCon, iSize) = mean(loc_vals1, 'omitnan');
            
            se_stat1(iCon, iSize) = std(stat_vals1, 'omitnan') / sqrt(length(cells1_this_condition));
            se_loc1(iCon, iSize) = std(loc_vals1, 'omitnan') / sqrt(length(cells1_this_condition));
            
            n_cells1(iCon, iSize) = length(cells1_this_condition);
        end
        
        % Population 2
        cells2_this_condition = intersect(cell_indices2, cells_with_running);
        if ~isempty(cells2_this_condition)
            stat_vals2 = statData(cells2_this_condition, iCon, iSize);
            loc_vals2 = locData(cells2_this_condition, iCon, iSize);
            
            avg_stat2(iCon, iSize) = mean(stat_vals2, 'omitnan');
            avg_loc2(iCon, iSize) = mean(loc_vals2, 'omitnan');
            
            se_stat2(iCon, iSize) = std(stat_vals2, 'omitnan') / sqrt(length(cells2_this_condition));
            se_loc2(iCon, iSize) = std(loc_vals2, 'omitnan') / sqrt(length(cells2_this_condition));
            
            n_cells2(iCon, iSize) = length(cells2_this_condition);
        end
    end
end

% Shared y-axis limits
allVals = [avg_stat1-se_stat1; avg_stat1+se_stat1; avg_loc1-se_loc1; avg_loc1+se_loc1; ...
           avg_stat2-se_stat2; avg_stat2+se_stat2; avg_loc2-se_loc2; avg_loc2+se_loc2];
ymin = min(allVals(:), [], 'omitnan'); ymax = max(allVals(:), [], 'omitnan');
pad  = 0.1 * (ymax - ymin);
ymin = ymin - pad; ymax = ymax + pad;

figure('Units', 'inches', 'Position', [5, 5, r.FigureSize(1), r.FigureSize(2)]);

% Population 1
subplot(1, 2, 1);
hold on;
legend_entries = cell(nCon*2, 1);
legend_count = 0;

for iCon = 1:nCon
    % Stationary (solid line)
    if any(~isnan(avg_stat1(iCon,:)))
        errorbar(sizes, avg_stat1(iCon,:), se_stat1(iCon,:), '-o', ...
            'Color', colors(iCon,:), 'LineWidth', 1.5, 'MarkerSize', 6, ...
            'MarkerFaceColor', colors(iCon,:));
        legend_count = legend_count + 1;
        legend_entries{legend_count} = sprintf('%.2f stat', contrasts(iCon));
    end
    
    % Running (dashed line)
    if any(~isnan(avg_loc1(iCon,:)))
        errorbar(sizes, avg_loc1(iCon,:), se_loc1(iCon,:), '--o', ...
            'Color', colors(iCon,:), 'LineWidth', 1.5, 'MarkerSize', 6, ...
            'MarkerFaceColor', 'none');
        legend_count = legend_count + 1;
        legend_entries{legend_count} = sprintf('%.2f run', contrasts(iCon));
    end
end

xlim(r.XLim); ylim([ymin ymax]); xticks(sizes);
set(gca, 'TickDir', 'out', 'Box', 'off'); grid off;
ylabel(r.YLabel); xlabel(r.XLabel);

% Show cell count range
n_range1 = [min(n_cells1(:), [], 'omitnan'), max(n_cells1(:), [], 'omitnan')];
if n_range1(1) == n_range1(2)
    n_str1 = sprintf('n=%d', n_range1(1));
else
    n_str1 = sprintf('n=%d-%d', n_range1(1), n_range1(2));
end
title(sprintf('%s %s', r.Titles{1}, n_str1), 'FontWeight', 'normal');

if r.ShowLegend
    legend(legend_entries(1:legend_count), 'Location', 'best', 'Box', 'off');
end

% Population 2
subplot(1, 2, 2);
hold on;
legend_count = 0;

for iCon = 1:nCon
    % Stationary (solid line)
    if any(~isnan(avg_stat2(iCon,:)))
        errorbar(sizes, avg_stat2(iCon,:), se_stat2(iCon,:), '-o', ...
            'Color', colors(iCon,:), 'LineWidth', 1.5, 'MarkerSize', 6, ...
            'MarkerFaceColor', colors(iCon,:));
        legend_count = legend_count + 1;
    end
    
    % Running (dashed line)
    if any(~isnan(avg_loc2(iCon,:)))
        errorbar(sizes, avg_loc2(iCon,:), se_loc2(iCon,:), '--o', ...
            'Color', colors(iCon,:), 'LineWidth', 1.5, 'MarkerSize', 6, ...
            'MarkerFaceColor', 'none');
        legend_count = legend_count + 1;
    end
end

xlim(r.XLim); ylim([ymin ymax]); xticks(sizes);
set(gca, 'TickDir', 'out', 'Box', 'off'); grid off;
xlabel(r.XLabel);

% Show cell count range
n_range2 = [min(n_cells2(:), [], 'omitnan'), max(n_cells2(:), [], 'omitnan')];
if n_range2(1) == n_range2(2)
    n_str2 = sprintf('n=%d', n_range2(1));
else
    n_str2 = sprintf('n=%d-%d', n_range2(1), n_range2(2));
end
title(sprintf('%s %s', r.Titles{2}, n_str2), 'FontWeight', 'normal');

if r.ShowLegend
    legend(legend_entries(1:legend_count), 'Location', 'best', 'Box', 'off');
end

end
