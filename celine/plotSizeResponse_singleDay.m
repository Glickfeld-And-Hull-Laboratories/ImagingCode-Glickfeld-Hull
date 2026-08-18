function plotSizeResponse_singleDay(data, cell_indices1, cell_indices2, contrasts, sizes, varargin)
% PLOTSIZERESPONSE_SINGLEDAY  Size tuning curves for two cell populations.
%
%   plotSizeResponse_singleDay(data, cell_indices1, cell_indices2, contrasts, sizes)
%   plotSizeResponse_singleDay(..., 'Name', Value)
%
% Inputs:
%   data          - nCells x nCon x nSize response amplitude array
%   cell_indices1 - row indices into data for population 1 (e.g. red cells)
%   cell_indices2 - row indices into data for population 2 (e.g. green cells)
%   contrasts     - 1 x nCon contrast values
%   sizes         - 1 x nSize size values
%
% Optional parameters:
%   'Colors'      - {color1, color2}  Default: {'k', 'b'}
%   'Titles'      - {title1, title2}  Default: {'Population 1', 'Population 2'}
%   'YLabel'      - y-axis label      Default: 'dF/F'
%   'XLabel'      - x-axis label      Default: ''
%   'XLim'        - [xmin xmax]       Default: auto from sizes
%   'FigureSize'  - [w h] inches      Default: [8, 3*nCon]

p = inputParser;
addRequired(p, 'data');
addRequired(p, 'cell_indices1');
addRequired(p, 'cell_indices2');
addRequired(p, 'contrasts');
addRequired(p, 'sizes');
addParameter(p, 'Colors',     {'k', 'b'});
addParameter(p, 'Titles',     {'Population 1', 'Population 2'});
addParameter(p, 'YLabel',     'dF/F');
addParameter(p, 'XLabel',     '');
addParameter(p, 'XLim',       []);
addParameter(p, 'FigureSize', []);
parse(p, data, cell_indices1, cell_indices2, contrasts, sizes, varargin{:});
r = p.Results;

nCon  = length(contrasts);
nSize = length(sizes);

if isempty(r.XLim),       r.XLim = [min(sizes)*0.8, max(sizes)*1.2]; end
if isempty(r.FigureSize), r.FigureSize = [8, 3*nCon]; end

% Compute mean ± SE for each population x contrast x size
avg1 = zeros(nCon, nSize);  se1 = zeros(nCon, nSize);
avg2 = zeros(nCon, nSize);  se2 = zeros(nCon, nSize);
for iCon = 1:nCon
    d1 = squeeze(data(cell_indices1, iCon, :));
    d2 = squeeze(data(cell_indices2, iCon, :));
    avg1(iCon,:) = mean(d1, 1, 'omitnan');
    avg2(iCon,:) = mean(d2, 1, 'omitnan');
    se1(iCon,:)  = std(d1, 0, 1, 'omitnan') / sqrt(length(cell_indices1));
    se2(iCon,:)  = std(d2, 0, 1, 'omitnan') / sqrt(length(cell_indices2));
end

% Shared y-axis limits
allVals = [avg1-se1; avg1+se1; avg2-se2; avg2+se2];
ymin = min(allVals(:)); ymax = max(allVals(:));
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
    title(sprintf('%s  Con: %g  n=%d', r.Titles{1}, contrasts(iCon), length(cell_indices1)), 'FontWeight', 'normal');

    subplot(nCon, 2, (iCon-1)*2 + 2);
    errorbar(sizes, avg2(iCon,:), se2(iCon,:), '-o', ...
        'Color', r.Colors{2}, 'LineWidth', 1.5, 'MarkerSize', 6);
    xlim(r.XLim); ylim([ymin ymax]); xticks(sizes);
    set(gca, 'TickDir', 'out', 'Box', 'off'); grid off;
    if iCon == nCon && ~isempty(r.XLabel), xlabel(r.XLabel); end
    title(sprintf('%s  Con: %g  n=%d', r.Titles{2}, contrasts(iCon), length(cell_indices2)), 'FontWeight', 'normal');
end
