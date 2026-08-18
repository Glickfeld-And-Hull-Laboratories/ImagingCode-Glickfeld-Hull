function plotSizeResponse_conditionMatched_singleDay(locData, cell_indices1, cell_indices2, contrasts, sizes, runningByCondition, varargin)
% PLOTSIZERESPONSE_CONDITIONMATCHED_SINGLEDAY  Running-only size tuning using condition-matched cells
%   Shows only running responses, but each condition uses only cells that have running data for that condition
%
%   plotSizeResponse_conditionMatched_singleDay(locData, cell_indices1, cell_indices2, contrasts, sizes, runningByCondition)
%   plotSizeResponse_conditionMatched_singleDay(..., 'Name', Value)
%
% Inputs:
%   locData - Running response data (nCells x nContrasts x nSizes)
%   cell_indices1 - Indices for population 1 (e.g., red cells)
%   cell_indices2 - Indices for population 2 (e.g., green cells)
%   contrasts - Array of contrast values
%   sizes - Array of size values
%   runningByCondition - nCells x nContrasts x nSizes logical matrix indicating which cells have running data
%
% Optional Name-Value Pairs:
%   'Colors'      - {color1, color2}  Default: {'k', 'b'}
%   'Titles'      - {title1, title2}  Default: {'Population 1', 'Population 2'}
%   'YLabel'      - y-axis label      Default: 'dF/F'
%   'XLabel'      - x-axis label      Default: ''
%   'XLim'        - [xmin xmax]       Default: auto from sizes
%   'FigureSize'  - [w h] inches      Default: [8, 3*nCon]

p = inputParser;
addRequired(p, 'locData');
addRequired(p, 'cell_indices1');
addRequired(p, 'cell_indices2');
addRequired(p, 'contrasts');
addRequired(p, 'sizes');
addRequired(p, 'runningByCondition');
addParameter(p, 'Colors',     {'k', 'b'});
addParameter(p, 'Titles',     {'Population 1', 'Population 2'});
addParameter(p, 'YLabel',     'dF/F');
addParameter(p, 'XLabel',     '');
addParameter(p, 'XLim',       []);
addParameter(p, 'FigureSize', []);
parse(p, locData, cell_indices1, cell_indices2, contrasts, sizes, runningByCondition, varargin{:});
r = p.Results;

nCon  = length(contrasts);
nSize = length(sizes);

if isempty(r.XLim),       r.XLim = [min(sizes)*0.8, max(sizes)*1.2]; end
if isempty(r.FigureSize), r.FigureSize = [8, 3*nCon]; end

% Compute condition-specific mean ± SE for each population
avg1 = nan(nCon, nSize);  se1 = nan(nCon, nSize);  n1 = nan(nCon, nSize);
avg2 = nan(nCon, nSize);  se2 = nan(nCon, nSize);  n2 = nan(nCon, nSize);

for iCon = 1:nCon
    for iSize = 1:nSize
        % Find cells that have running data for this specific condition
        cells_with_running = find(runningByCondition(:, iCon, iSize));
        
        % Population 1
        cells1_this_condition = intersect(cell_indices1, cells_with_running);
        if ~isempty(cells1_this_condition)
            vals1 = locData(cells1_this_condition, iCon, iSize);
            avg1(iCon, iSize) = mean(vals1, 'omitnan');
            se1(iCon, iSize) = std(vals1, 'omitnan') / sqrt(length(cells1_this_condition));
            n1(iCon, iSize) = length(cells1_this_condition);
        end
        
        % Population 2  
        cells2_this_condition = intersect(cell_indices2, cells_with_running);
        if ~isempty(cells2_this_condition)
            vals2 = locData(cells2_this_condition, iCon, iSize);
            avg2(iCon, iSize) = mean(vals2, 'omitnan');
            se2(iCon, iSize) = std(vals2, 'omitnan') / sqrt(length(cells2_this_condition));
            n2(iCon, iSize) = length(cells2_this_condition);
        end
    end
end

% Shared y-axis limits
allVals = [avg1-se1; avg1+se1; avg2-se2; avg2+se2];
ymin = min(allVals(:), [], 'omitnan'); ymax = max(allVals(:), [], 'omitnan');
pad  = 0.1 * (ymax - ymin);
ymin = ymin - pad; ymax = ymax + pad;

figure('Units', 'inches', 'Position', [5, 5, r.FigureSize(1), r.FigureSize(2)]);

for iCon = 1:nCon
    subplot(nCon, 2, (iCon-1)*2 + 1);
    errorbar(sizes, avg1(iCon,:), se1(iCon,:), '-o', ...
        'Color', r.Colors{1}, 'LineWidth', 1.5, 'MarkerSize', 6);
    xlim(r.XLim); ylim([ymin ymax]); xticks(sizes);
    set(gca, 'TickDir', 'out', 'Box', 'off'); grid off;
    if iCon == 1, ylabel(r.YLabel); end
    if iCon == nCon && ~isempty(r.XLabel), xlabel(r.XLabel); end
    
    % Show cell count range for this contrast
    n_range1 = [min(n1(iCon, :), [], 'omitnan'), max(n1(iCon, :), [], 'omitnan')];
    if n_range1(1) == n_range1(2)
        n_str1 = sprintf('n=%d', n_range1(1));
    else
        n_str1 = sprintf('n=%d-%d', n_range1(1), n_range1(2));
    end
    title(sprintf('%s  Con: %g  %s', r.Titles{1}, contrasts(iCon), n_str1), 'FontWeight', 'normal');

    subplot(nCon, 2, (iCon-1)*2 + 2);
    errorbar(sizes, avg2(iCon,:), se2(iCon,:), '-o', ...
        'Color', r.Colors{2}, 'LineWidth', 1.5, 'MarkerSize', 6);
    xlim(r.XLim); ylim([ymin ymax]); xticks(sizes);
    set(gca, 'TickDir', 'out', 'Box', 'off'); grid off;
    if iCon == nCon && ~isempty(r.XLabel), xlabel(r.XLabel); end
    
    % Show cell count range for this contrast
    n_range2 = [min(n2(iCon, :), [], 'omitnan'), max(n2(iCon, :), [], 'omitnan')];
    if n_range2(1) == n_range2(2)
        n_str2 = sprintf('n=%d', n_range2(1));
    else
        n_str2 = sprintf('n=%d-%d', n_range2(1), n_range2(2));
    end
    title(sprintf('%s  Con: %g  %s', r.Titles{2}, contrasts(iCon), n_str2), 'FontWeight', 'normal');
end
end
