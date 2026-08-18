function plotNeuralTimecourse_singleDay(tc, cell_indices1, cell_indices2, varargin)
% PLOTNEURALTIMECOURSE_SINGLEDAY  Timecourses for two cell populations.
%   One figure per stimulus size; rows = contrasts; left col = pop1, right = pop2.
%
%   plotNeuralTimecourse_singleDay(tc, cell_indices1, cell_indices2)
%   plotNeuralTimecourse_singleDay(..., 'Name', Value)
%
% Inputs:
%   tc            - nFrames x nCells x nCon x nSize trial-averaged timecourse
%   cell_indices1 - indices for population 1
%   cell_indices2 - indices for population 2
%
% Optional parameters:
%   'Colors'       - {color1, color2}  Default: {'k', 'b'}
%   'Titles'       - {title1, title2}  Default: {'Population 1', 'Population 2'}
%   'StimDuration' - stimulus duration in seconds  Default: 2
%   'FrameRate'    - Hz                Default: 15
%   'StimStart'    - frame index of stimulus onset  Default: floor(nFrames/3)
%   'FigureSize'   - [w h] inches      Default: [4, 9]

p = inputParser;
addRequired(p, 'tc');
addRequired(p, 'cell_indices1');
addRequired(p, 'cell_indices2');
addParameter(p, 'Colors',       {'k', 'b'});
addParameter(p, 'Titles',       {'Population 1', 'Population 2'});
addParameter(p, 'StimDuration', 2);
addParameter(p, 'FrameRate',    15);
addParameter(p, 'StimStart',    []);
addParameter(p, 'FigureSize',   [4, 9]);
parse(p, tc, cell_indices1, cell_indices2, varargin{:});
r = p.Results;

[nFrames, ~, nCon, nSize] = size(tc);
if isempty(r.StimStart), r.StimStart = floor(nFrames/3); end

% Pre-compute mean ± SE for all sizes to find global y limits
avg1 = cell(nCon, nSize); se1 = cell(nCon, nSize);
avg2 = cell(nCon, nSize); se2 = cell(nCon, nSize);
for iSize = 1:nSize
    for iCon = 1:nCon
        d1 = squeeze(tc(:, cell_indices1, iCon, iSize));
        d2 = squeeze(tc(:, cell_indices2, iCon, iSize));
        avg1{iCon,iSize} = mean(d1, 2, 'omitnan');
        avg2{iCon,iSize} = mean(d2, 2, 'omitnan');
        se1{iCon,iSize}  = std(d1, 0, 2, 'omitnan') / sqrt(length(cell_indices1));
        se2{iCon,iSize}  = std(d2, 0, 2, 'omitnan') / sqrt(length(cell_indices2));
    end
end

allMin = cellfun(@(a,s) min(a-s), [avg1(:);avg2(:)], [se1(:);se2(:)]);
allMax = cellfun(@(a,s) max(a+s), [avg1(:);avg2(:)], [se1(:);se2(:)]);
ymin = min(allMin) - 0.1*range([min(allMin) max(allMax)]);
ymax = max(allMax) + 0.1*range([min(allMin) max(allMax)]);
pad  = 0.1 * (ymax - ymin);

t = ((1:nFrames) - (double(r.StimStart) - 1)) / double(r.FrameRate);
z = r.StimDuration;

for iSize = 1:nSize
    figure('Units', 'inches', 'Position', [5+(iSize-1)*0.5, 0, r.FigureSize(1), r.FigureSize(2)]);

    for iCon = 1:nCon
        % Population 1
        subplot(nCon, 2, (iCon-1)*2 + 1);
        shadedErrorBar(t, avg1{iCon,iSize}, se1{iCon,iSize}, r.Colors{1});
        ylim([ymin ymax]);
        hold on;
        line([0 z], [ymin+0.1*pad, ymin+0.1*pad], 'Color', 'k', 'LineWidth', 2);
        set(gca, 'TickDir', 'out', 'XColor', 'none', 'YColor', 'none', 'Box', 'off'); grid off;
        if iCon == 1
            title(sprintf('%s  Size %d  n=%d', r.Titles{1}, iSize, length(cell_indices1)));
            line([-1.8 -1.8], [ymin+2*pad, ymin+2*pad+0.05], 'Color', 'k', 'LineWidth', 2);
        end

        % Population 2
        subplot(nCon, 2, (iCon-1)*2 + 2);
        shadedErrorBar(t, avg2{iCon,iSize}, se2{iCon,iSize}, r.Colors{2});
        ylim([ymin ymax]);
        hold on;
        line([0 z], [ymin+0.1*pad, ymin+0.1*pad], 'Color', 'k', 'LineWidth', 2);
        set(gca, 'TickDir', 'out', 'XColor', 'none', 'YColor', 'none', 'Box', 'off'); grid off;
        if iCon == 1
            title(sprintf('%s  Size %d  n=%d', r.Titles{2}, iSize, length(cell_indices2)));
            line([-1.8 -1.8], [ymin+2*pad, ymin+2*pad+0.05], 'Color', 'k', 'LineWidth', 2);
        end
    end
end
