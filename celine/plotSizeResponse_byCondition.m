function plotSizeResponse_byCondition(data1, data2, cell_indices1, cell_indices2, contrasts, sizes, runningByCondition, varargin)

p = inputParser;
addRequired(p, 'data1');
addRequired(p, 'data2');
addRequired(p, 'cell_indices1');
addRequired(p, 'cell_indices2');
addRequired(p, 'contrasts');
addRequired(p, 'sizes');
addRequired(p, 'runningByCondition');
addParameter(p, 'Colors1', {'k', 'b'});
addParameter(p, 'Colors2', {'k', 'b'});
addParameter(p, 'XLim', []);
addParameter(p, 'Titles', {'Dataset 1', 'Dataset 2'});
addParameter(p, 'YLabel', 'dF/F');
addParameter(p, 'XLabel', '');
addParameter(p, 'FigureSize', []);
addParameter(p, 'DayOrder', 'r');

parse(p, data1, data2, cell_indices1, cell_indices2, contrasts, sizes, runningByCondition, varargin{:});

colors1     = p.Results.Colors1;
colors2     = p.Results.Colors2;
xlim_range  = p.Results.XLim;
titles      = p.Results.Titles;
y_label     = p.Results.YLabel;
x_label     = p.Results.XLabel;
figure_size = p.Results.FigureSize;

if isempty(xlim_range)
    xlim_range = [min(sizes)*0.8, max(sizes)*1.2];
end

nContrasts = length(contrasts);
nSizes     = length(sizes);
nd         = size(data1, 2);

if isempty(figure_size)
    figure_size = [8, 3*nContrasts];
end

if strcmp(p.Results.DayOrder, 'f')
    pre = 1; post = 2;
else
    pre = 2; post = 1;
end

sizeResp_data1_avrg = nan(nd, nContrasts, nSizes);
sizeResp_data2_avrg = nan(nd, nContrasts, nSizes);
sizeResp_data1_se   = nan(nd, nContrasts, nSizes);
sizeResp_data2_se   = nan(nd, nContrasts, nSizes);
n_cells_data1       = nan(nd, nContrasts, nSizes);
n_cells_data2       = nan(nd, nContrasts, nSizes);

for id = 1:nd
    for iContrast = 1:nContrasts
        for iSize = 1:nSizes
            cells_with_running = find(runningByCondition(:, iContrast, iSize));

            cells_data1 = intersect(cell_indices1, cells_with_running);
            if ~isempty(cells_data1)
                data1_vals = data1{id}(cells_data1, iContrast, iSize);
                sizeResp_data1_avrg(id, iContrast, iSize) = mean(data1_vals, 'omitnan');
                sizeResp_data1_se(id, iContrast, iSize)   = std(data1_vals, 'omitnan') / sqrt(length(cells_data1));
                n_cells_data1(id, iContrast, iSize)        = length(cells_data1);
            end

            cells_data2 = intersect(cell_indices2, cells_with_running);
            if ~isempty(cells_data2)
                data2_vals = data2{id}(cells_data2, iContrast, iSize);
                sizeResp_data2_avrg(id, iContrast, iSize) = mean(data2_vals, 'omitnan');
                sizeResp_data2_se(id, iContrast, iSize)   = std(data2_vals, 'omitnan') / sqrt(length(cells_data2));
                n_cells_data2(id, iContrast, iSize)        = length(cells_data2);
            end
        end
    end
end

ymin = min([sizeResp_data1_avrg(:) - sizeResp_data1_se(:); sizeResp_data2_avrg(:) - sizeResp_data2_se(:)], [], 'omitnan');
ymax = max([sizeResp_data1_avrg(:) + sizeResp_data1_se(:); sizeResp_data2_avrg(:) + sizeResp_data2_se(:)], [], 'omitnan');
padding = 0.1 * (ymax - ymin);
ymin = ymin - padding;
ymax = ymax + padding;

figure('Units', 'inches', 'Position', [5, 5, figure_size(1), figure_size(2)]);

for iContrast = 1:nContrasts
    subplot(nContrasts, 2, (iContrast-1)*2 + 1);

    errorbar(sizes, squeeze(sizeResp_data1_avrg(pre, iContrast, :)), squeeze(sizeResp_data1_se(pre, iContrast, :)), ...
        'o', 'Color', colors1{1}, 'LineWidth', 1.5, 'MarkerSize', 6, 'LineStyle', 'none');
    hold on;
    errorbar(sizes, squeeze(sizeResp_data1_avrg(post, iContrast, :)), squeeze(sizeResp_data1_se(post, iContrast, :)), ...
        'o', 'Color', colors1{2}, 'LineWidth', 1.5, 'MarkerSize', 6, 'LineStyle', 'none');

    if iContrast == 1, ylabel(y_label); end
    if iContrast == nContrasts && ~isempty(x_label), xlabel(x_label); end
    xlim(xlim_range); ylim([ymin, ymax]); xticks(sizes);
    set(gca, 'TickDir', 'out', 'box', 'off'); grid off;

    n1 = n_cells_data1([pre post], iContrast, :);
    title(sprintf('%s (Contrast: %g, N=%d–%d)', titles{1}, contrasts(iContrast), min(n1(:),[],'omitnan'), max(n1(:),[],'omitnan')), 'FontWeight', 'normal');

    subplot(nContrasts, 2, (iContrast-1)*2 + 2);

    errorbar(sizes, squeeze(sizeResp_data2_avrg(pre, iContrast, :)), squeeze(sizeResp_data2_se(pre, iContrast, :)), ...
        'o', 'Color', colors2{1}, 'LineWidth', 1.5, 'MarkerSize', 6, 'LineStyle', 'none');
    hold on;
    errorbar(sizes, squeeze(sizeResp_data2_avrg(post, iContrast, :)), squeeze(sizeResp_data2_se(post, iContrast, :)), ...
        'o', 'Color', colors2{2}, 'LineWidth', 1.5, 'MarkerSize', 6, 'LineStyle', 'none');

    if iContrast == nContrasts && ~isempty(x_label), xlabel(x_label); end
    xlim(xlim_range); ylim([ymin, ymax]); xticks(sizes);
    set(gca, 'TickDir', 'out', 'box', 'off'); grid off;

    n2 = n_cells_data2([pre post], iContrast, :);
    title(sprintf('%s (Contrast: %g, N=%d–%d)', titles{2}, contrasts(iContrast), min(n2(:),[],'omitnan'), max(n2(:),[],'omitnan')), 'FontWeight', 'normal');
end

end