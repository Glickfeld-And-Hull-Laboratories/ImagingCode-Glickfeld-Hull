function plotSizeResponse_contrastOverlay_singleDay(data, cell_indices1, cell_indices2, contrasts, sizes, varargin)
% PLOTSIZERESPONSE_CONTRASTOVERLAY_SINGLEDAY  Size tuning with all contrasts overlaid
%   Shows all contrast conditions on the same plot using grayscale where darker = higher contrast.
%   Marker size scales with the number of cells contributing to each size point.
%
%   plotSizeResponse_contrastOverlay_singleDay(data, cell_indices1, cell_indices2, contrasts, sizes)
%   plotSizeResponse_contrastOverlay_singleDay(..., 'Name', Value)
%
% Inputs:
%   data          - nCells x nCon x nSize response amplitude array
%   cell_indices1 - row indices into data for population 1 (e.g. red cells)
%   cell_indices2 - row indices into data for population 2 (e.g. green cells)
%   contrasts     - 1 x nCon contrast values
%   sizes         - 1 x nSize size values
%
% Optional parameters:
%   'Titles'      - {title1, title2}  Default: {'Population 1', 'Population 2'}
%   'YLabel'      - y-axis label      Default: 'dF/F'
%   'XLabel'      - x-axis label      Default: 'Size (deg)'
%   'XLim'        - [xmin xmax]       Default: auto from sizes
%   'FigureSize'  - [w h] inches      Default: [8, 4]
%   'ShowLegend'  - logical           Default: true
%   'MarkerRange' - [minSz maxSz]     Marker size range in points. Default: [4, 14]

p = inputParser;
addRequired(p, 'data');
addRequired(p, 'cell_indices1');
addRequired(p, 'cell_indices2');
addRequired(p, 'contrasts');
addRequired(p, 'sizes');
addParameter(p, 'Titles',      {'Population 1', 'Population 2'});
addParameter(p, 'YLabel',      'dF/F');
addParameter(p, 'XLabel',      'Size (deg)');
addParameter(p, 'XLim',        []);
addParameter(p, 'FigureSize',  [8, 4]);
addParameter(p, 'ShowLegend',  true);
addParameter(p, 'MarkerRange', [4, 14]);
parse(p, data, cell_indices1, cell_indices2, contrasts, sizes, varargin{:});
r = p.Results;

nCon  = length(contrasts);
nSize = length(sizes);

if isempty(r.XLim), r.XLim = [min(sizes)*0.8, max(sizes)*1.2]; end

gray_values = linspace(0.7, 0.0, nCon);
colors = repmat(gray_values', 1, 3);

% Compute mean ± SE and N for each population x contrast x size
avg1 = zeros(nCon, nSize);  se1 = zeros(nCon, nSize);  n1 = zeros(nCon, nSize);
avg2 = zeros(nCon, nSize);  se2 = zeros(nCon, nSize);  n2 = zeros(nCon, nSize);
for iCon = 1:nCon
    d1 = squeeze(data(cell_indices1, iCon, :));
    d2 = squeeze(data(cell_indices2, iCon, :));
    avg1(iCon,:) = mean(d1, 1, 'omitnan');
    avg2(iCon,:) = mean(d2, 1, 'omitnan');
    se1(iCon,:)  = std(d1, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(d1), 1));
    se2(iCon,:)  = std(d2, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(d2), 1));
    n1(iCon,:)   = sum(~isnan(d1), 1);
    n2(iCon,:)   = sum(~isnan(d2), 1);
end

% Map N -> marker size independently per population
mlo = r.MarkerRange(1);
mhi = r.MarkerRange(2);
makeScaler = @(n_mat) deal( ...
    min(n_mat(n_mat > 0)), ...
    max(n_mat(:)));
[nMin1, nMax1] = makeScaler(n1);
[nMin2, nMax2] = makeScaler(n2);
scaleN1 = @(n) mlo + (mhi - mlo) * (n - nMin1) / max(nMax1 - nMin1, 1);
scaleN2 = @(n) mlo + (mhi - mlo) * (n - nMin2) / max(nMax2 - nMin2, 1);

allVals = [avg1-se1; avg1+se1; avg2-se2; avg2+se2];
ymin = min(allVals(:)); ymax = max(allVals(:));
pad  = 0.1 * (ymax - ymin);
ymin = ymin - pad; ymax = ymax + pad;

figure('Units', 'inches', 'Position', [5, 5, r.FigureSize(1), r.FigureSize(2)]);

legend_entries = arrayfun(@(c) sprintf('%.2f', c), contrasts, 'UniformOutput', false);

for iPop = 1:2
    subplot(1, 2, iPop);
    hold on;

    if iPop == 1
        avg = avg1; se = se1; n = n1; idx = cell_indices1;
        scaleN = scaleN1; nMin = nMin1; nMax = nMax1;
    else
        avg = avg2; se = se2; n = n2; idx = cell_indices2;
        scaleN = scaleN2; nMin = nMin2; nMax = nMax2;
    end

    for iCon = 1:nCon
        c = colors(iCon,:);
        for iSz = 1:nSize
            if n(iCon, iSz) == 0, continue; end
            ms = scaleN(n(iCon, iSz));
            errorbar(sizes(iSz), avg(iCon,iSz), se(iCon,iSz), '-o', ...
                'Color', c, 'LineWidth', 1.5, 'MarkerSize', ms, ...
                'MarkerFaceColor', c, 'CapSize', 4);
        end
        % Overlay connecting line (no markers, so per-contrast legend handle is clean)
        plot(sizes, avg(iCon,:), '-', 'Color', c, 'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
    end

    xlim(r.XLim); ylim([ymin ymax]); xticks(sizes);
    set(gca, 'TickDir', 'out', 'Box', 'off'); grid off;
    xlabel(r.XLabel);
    if iPop == 1
        ylabel(r.YLabel);
    end
    title(sprintf('%s (max n=%d)', r.Titles{iPop}, max(n(:))), 'FontWeight', 'normal');

    if r.ShowLegend
        hCon = gobjects(nCon, 1);
        for iCon = 1:nCon
            hCon(iCon) = plot(nan, nan, '-o', 'Color', colors(iCon,:), ...
                'MarkerFaceColor', colors(iCon,:), 'MarkerSize', 6, 'LineWidth', 1.5);
        end
        nLevels = unique(n(n > 0));
        hN = gobjects(length(nLevels), 1);
        nLabels = cell(length(nLevels), 1);
        for k = 1:length(nLevels)
            ms = scaleN(nLevels(k));
            hN(k) = plot(nan, nan, 'o', 'Color', [0.4 0.4 0.4], ...
                'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerSize', ms, 'LineWidth', 1);
            nLabels{k} = sprintf('n=%d', nLevels(k));
        end
        legend([hCon; hN], [legend_entries(:); nLabels(:)], ...
            'Location', 'best', 'Box', 'off');
    end
end

end