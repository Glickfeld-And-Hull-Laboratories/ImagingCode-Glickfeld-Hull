function plotSizeResponse_byCondition_singleDay(statData, locData, cell_indices1, cell_indices2, contrasts, sizes, runningByCondition, varargin)
% PLOTSIZERESPONSE_BYCONDITION_SINGLEDAY Plot size response curves using cells with running data for each condition
%   Single-day version - compares stationary vs running for same cells by condition
% 
% Usage:
%   plotSizeResponse_byCondition_singleDay(statData, locData, cell_indices1, cell_indices2, contrasts, sizes, runningByCondition) 
%   plotSizeResponse_byCondition_singleDay(..., 'Name', Value, ...)
%
% Inputs:
%   statData - Stationary response data (nCells x nContrasts x nSizes)
%   locData - Running response data (nCells x nContrasts x nSizes)
%   cell_indices1 - Indices for population 1 (e.g., red cells)
%   cell_indices2 - Indices for population 2 (e.g., green cells)
%   contrasts - Array of contrast values
%   sizes - Array of size values
%   runningByCondition - nCells x nContrasts x nSizes logical matrix indicating which cells have running data
%
% Optional Name-Value Pairs:
%   'Colors1' - Colors for [stationary, running] in population 1
%               Default: {'k', 'r'}
%   'Colors2' - Colors for [stationary, running] in population 2
%               Default: {'k', 'r'}
%   'XLim' - X axis limits [min, max]
%            Default: [min(sizes)*0.8, max(sizes)*1.2]
%   'Titles' - Cell array with titles for each population {'title1', 'title2'}
%              Default: {'Population 1', 'Population 2'}
%   'YLabel' - Label for y-axis
%              Default: 'dF/F'
%   'XLabel' - Label for x-axis
%              Default: '' (no x-label)
%   'FigureSize' - Size of the figure [width, height] in inches
%                 Default: [8, 3*nContrasts]

p = inputParser;
addRequired(p, 'statData');
addRequired(p, 'locData');
addRequired(p, 'cell_indices1');
addRequired(p, 'cell_indices2');
addRequired(p, 'contrasts');
addRequired(p, 'sizes');
addRequired(p, 'runningByCondition');
addParameter(p, 'Colors1', {'k', 'r'});
addParameter(p, 'Colors2', {'k', 'r'});
addParameter(p, 'XLim', []);
addParameter(p, 'Titles', {'Population 1', 'Population 2'});
addParameter(p, 'YLabel', 'dF/F');
addParameter(p, 'XLabel', '');
addParameter(p, 'FigureSize', []);

parse(p, statData, locData, cell_indices1, cell_indices2, contrasts, sizes, runningByCondition, varargin{:});

colors1 = p.Results.Colors1;
colors2 = p.Results.Colors2;
xlim_range = p.Results.XLim;
titles = p.Results.Titles;
y_label = p.Results.YLabel;
x_label = p.Results.XLabel;
figure_size = p.Results.FigureSize;

if isempty(xlim_range)
    xlim_range = [min(sizes)*0.8, max(sizes)*1.2];
end

nContrasts = length(contrasts);
nSizes = length(sizes);

if isempty(figure_size)
    figure_size = [8, 3*nContrasts];
end

% Initialize arrays for responses and standard errors
sizeResp_stat1_avrg = nan(nContrasts, nSizes);
sizeResp_loc1_avrg = nan(nContrasts, nSizes);
sizeResp_stat2_avrg = nan(nContrasts, nSizes);
sizeResp_loc2_avrg = nan(nContrasts, nSizes);

sizeResp_stat1_se = nan(nContrasts, nSizes);
sizeResp_loc1_se = nan(nContrasts, nSizes);
sizeResp_stat2_se = nan(nContrasts, nSizes);
sizeResp_loc2_se = nan(nContrasts, nSizes);

n_cells1 = nan(nContrasts, nSizes);
n_cells2 = nan(nContrasts, nSizes);

% Calculate responses for each condition, using only cells that have running data for that condition
for iContrast = 1:nContrasts
    for iSize = 1:nSizes
        % Find cells that have running data for this specific condition
        cells_with_running = find(runningByCondition(:, iContrast, iSize));
        
        % Population 1 (e.g., red cells)
        cells1_this_condition = intersect(cell_indices1, cells_with_running);
        if ~isempty(cells1_this_condition)
            stat_vals1 = statData(cells1_this_condition, iContrast, iSize);
            loc_vals1 = locData(cells1_this_condition, iContrast, iSize);
            
            sizeResp_stat1_avrg(iContrast, iSize) = mean(stat_vals1, 'omitnan');
            sizeResp_loc1_avrg(iContrast, iSize) = mean(loc_vals1, 'omitnan');
            
            sizeResp_stat1_se(iContrast, iSize) = std(stat_vals1, 'omitnan') / sqrt(length(cells1_this_condition));
            sizeResp_loc1_se(iContrast, iSize) = std(loc_vals1, 'omitnan') / sqrt(length(cells1_this_condition));
            
            n_cells1(iContrast, iSize) = length(cells1_this_condition);
        end
        
        % Population 2 (e.g., green cells)
        cells2_this_condition = intersect(cell_indices2, cells_with_running);
        if ~isempty(cells2_this_condition)
            stat_vals2 = statData(cells2_this_condition, iContrast, iSize);
            loc_vals2 = locData(cells2_this_condition, iContrast, iSize);
            
            sizeResp_stat2_avrg(iContrast, iSize) = mean(stat_vals2, 'omitnan');
            sizeResp_loc2_avrg(iContrast, iSize) = mean(loc_vals2, 'omitnan');
            
            sizeResp_stat2_se(iContrast, iSize) = std(stat_vals2, 'omitnan') / sqrt(length(cells2_this_condition));
            sizeResp_loc2_se(iContrast, iSize) = std(loc_vals2, 'omitnan') / sqrt(length(cells2_this_condition));
            
            n_cells2(iContrast, iSize) = length(cells2_this_condition);
        end
    end
end

% Find global y-axis limits
all_vals = [sizeResp_stat1_avrg(:) - sizeResp_stat1_se(:); 
            sizeResp_stat1_avrg(:) + sizeResp_stat1_se(:);
            sizeResp_loc1_avrg(:) - sizeResp_loc1_se(:); 
            sizeResp_loc1_avrg(:) + sizeResp_loc1_se(:);
            sizeResp_stat2_avrg(:) - sizeResp_stat2_se(:); 
            sizeResp_stat2_avrg(:) + sizeResp_stat2_se(:);
            sizeResp_loc2_avrg(:) - sizeResp_loc2_se(:); 
            sizeResp_loc2_avrg(:) + sizeResp_loc2_se(:)];

ymin = min(all_vals, [], 'omitnan');
ymax = max(all_vals, [], 'omitnan');
padding = 0.1 * (ymax - ymin);
ymin = ymin - padding;
ymax = ymax + padding;

% Create figure
figure('Units', 'inches', 'Position', [5, 5, figure_size(1), figure_size(2)]);

for iContrast = 1:nContrasts
    % Population 1 plot
    subplot(nContrasts, 2, (iContrast-1)*2 + 1);
    
    % Stationary data
    errorbar(sizes, sizeResp_stat1_avrg(iContrast, :), sizeResp_stat1_se(iContrast, :), ...
        'o', 'Color', colors1{1}, 'LineWidth', 1.5, 'MarkerSize', 6, 'LineStyle', 'none');
    hold on;
    
    % Running data
    errorbar(sizes, sizeResp_loc1_avrg(iContrast, :), sizeResp_loc1_se(iContrast, :), ...
        'o', 'Color', colors1{2}, 'LineWidth', 1.5, 'MarkerSize', 6, 'LineStyle', 'none');
    
    if iContrast == 1
        ylabel(y_label);
    end
    xlim(xlim_range);
    ylim([ymin, ymax]);
    xticks(sizes);
    set(gca, 'TickDir', 'out', 'box', 'off');
    grid off;
    if iContrast == nContrasts && ~isempty(x_label)
        xlabel(x_label);
    end
    
    % Show range of cell counts for this population
    n_range1 = [min(n_cells1(iContrast, :), [], 'omitnan'), max(n_cells1(iContrast, :), [], 'omitnan')];
    if n_range1(1) == n_range1(2)
        n_str1 = sprintf('n=%d', n_range1(1));
    else
        n_str1 = sprintf('n=%d-%d', n_range1(1), n_range1(2));
    end
    
    title([titles{1}, ' (Contrast: ', num2str(contrasts(iContrast)), ') ', n_str1], 'FontWeight', 'normal');
    
    % Population 2 plot
    subplot(nContrasts, 2, (iContrast-1)*2 + 2);
    
    % Stationary data
    errorbar(sizes, sizeResp_stat2_avrg(iContrast, :), sizeResp_stat2_se(iContrast, :), ...
        'o', 'Color', colors2{1}, 'LineWidth', 1.5, 'MarkerSize', 6, 'LineStyle', 'none');
    hold on;
    
    % Running data
    errorbar(sizes, sizeResp_loc2_avrg(iContrast, :), sizeResp_loc2_se(iContrast, :), ...
        'o', 'Color', colors2{2}, 'LineWidth', 1.5, 'MarkerSize', 6, 'LineStyle', 'none');
    
    xlim(xlim_range);
    ylim([ymin, ymax]);
    xticks(sizes);
    set(gca, 'TickDir', 'out', 'box', 'off');
    grid off;
    if iContrast == nContrasts && ~isempty(x_label)
        xlabel(x_label);
    end
    
    % Show range of cell counts for this population
    n_range2 = [min(n_cells2(iContrast, :), [], 'omitnan'), max(n_cells2(iContrast, :), [], 'omitnan')];
    if n_range2(1) == n_range2(2)
        n_str2 = sprintf('n=%d', n_range2(1));
    else
        n_str2 = sprintf('n=%d-%d', n_range2(1), n_range2(2));
    end
    
    title([titles{2}, ' (Contrast: ', num2str(contrasts(iContrast)), ') ', n_str2], 'FontWeight', 'normal');
end

end