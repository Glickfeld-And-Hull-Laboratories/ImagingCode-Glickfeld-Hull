function plotNeuralTimecourse_byCondition_singleDay(tc_stat, tc_loc, cell_indices1, cell_indices2, runningByCondition, varargin)
% PLOTNEURALTIMECOURSE_BYCONDITION_SINGLEDAY  Timecourses comparing stationary vs running for condition-matched cells
%   One figure per stimulus size; rows = contrasts; compares stationary vs running for cells that have both
%
%   plotNeuralTimecourse_byCondition_singleDay(tc_stat, tc_loc, cell_indices1, cell_indices2, runningByCondition)
%   plotNeuralTimecourse_byCondition_singleDay(..., 'Name', Value)
%
% Inputs:
%   tc_stat           - nFrames x nCells x nCon x nSize stationary timecourses
%   tc_loc            - nFrames x nCells x nCon x nSize running timecourses
%   cell_indices1     - indices for population 1 (e.g., red cells)
%   cell_indices2     - indices for population 2 (e.g., green cells)
%   runningByCondition - nCells x nCon x nSize logical array indicating which cells have running data
%
% Optional parameters:
%   'Colors1'      - {stat_color, run_color} for population 1  Default: {'k', 'r'}
%   'Colors2'      - {stat_color, run_color} for population 2  Default: {'k', 'r'}
%   'Titles'       - {title1, title2}  Default: {'Population 1', 'Population 2'}
%   'StimDuration' - stimulus duration in seconds  Default: 2
%   'FrameRate'    - Hz                Default: 15
%   'StimStart'    - frame index of stimulus onset  Default: floor(nFrames/3)
%   'FigureSize'   - [w h] inches      Default: [4, 9]

p = inputParser;
addRequired(p, 'tc_stat');
addRequired(p, 'tc_loc');
addRequired(p, 'cell_indices1');
addRequired(p, 'cell_indices2');
addRequired(p, 'runningByCondition');
addParameter(p, 'Colors1',      {'k', 'r'});
addParameter(p, 'Colors2',      {'k', 'r'});
addParameter(p, 'Titles',       {'Population 1', 'Population 2'});
addParameter(p, 'StimDuration', 2);
addParameter(p, 'FrameRate',    15);
addParameter(p, 'StimStart',    []);
addParameter(p, 'FigureSize',   [4, 9]);
parse(p, tc_stat, tc_loc, cell_indices1, cell_indices2, runningByCondition, varargin{:});
r = p.Results;

[nFrames, ~, nCon, nSize] = size(tc_stat);
if isempty(r.StimStart), r.StimStart = floor(nFrames/3); end

% Pre-compute mean ± SE for all sizes and conditions to find global y limits
avg_stat1 = cell(nCon, nSize); se_stat1 = cell(nCon, nSize);
avg_loc1  = cell(nCon, nSize); se_loc1  = cell(nCon, nSize);
avg_stat2 = cell(nCon, nSize); se_stat2 = cell(nCon, nSize);
avg_loc2  = cell(nCon, nSize); se_loc2  = cell(nCon, nSize);

n_cells1 = nan(nCon, nSize);
n_cells2 = nan(nCon, nSize);

for iSize = 1:nSize
    for iCon = 1:nCon
        % Find cells that have running data for this specific condition
        cells_with_running = find(runningByCondition(:, iCon, iSize));
        
        % Population 1 - cells that are in population 1 AND have running data for this condition
        cells1_this_condition = intersect(cell_indices1, cells_with_running);
        if ~isempty(cells1_this_condition)
            d_stat1 = squeeze(tc_stat(:, cells1_this_condition, iCon, iSize));
            d_loc1  = squeeze(tc_loc(:,  cells1_this_condition, iCon, iSize));
            
            avg_stat1{iCon,iSize} = mean(d_stat1, 2, 'omitnan');
            avg_loc1{iCon,iSize}  = mean(d_loc1,  2, 'omitnan');
            
            se_stat1{iCon,iSize} = std(d_stat1, 0, 2, 'omitnan') / sqrt(length(cells1_this_condition));
            se_loc1{iCon,iSize}  = std(d_loc1,  0, 2, 'omitnan') / sqrt(length(cells1_this_condition));
            
            n_cells1(iCon, iSize) = length(cells1_this_condition);
        else
            % No cells for this condition - fill with NaN
            avg_stat1{iCon,iSize} = nan(nFrames, 1);
            avg_loc1{iCon,iSize}  = nan(nFrames, 1);
            se_stat1{iCon,iSize}  = nan(nFrames, 1);
            se_loc1{iCon,iSize}   = nan(nFrames, 1);
            n_cells1(iCon, iSize) = 0;
        end
        
        % Population 2 - cells that are in population 2 AND have running data for this condition
        cells2_this_condition = intersect(cell_indices2, cells_with_running);
        if ~isempty(cells2_this_condition)
            d_stat2 = squeeze(tc_stat(:, cells2_this_condition, iCon, iSize));
            d_loc2  = squeeze(tc_loc(:,  cells2_this_condition, iCon, iSize));
            
            avg_stat2{iCon,iSize} = mean(d_stat2, 2, 'omitnan');
            avg_loc2{iCon,iSize}  = mean(d_loc2,  2, 'omitnan');
            
            se_stat2{iCon,iSize} = std(d_stat2, 0, 2, 'omitnan') / sqrt(length(cells2_this_condition));
            se_loc2{iCon,iSize}  = std(d_loc2,  0, 2, 'omitnan') / sqrt(length(cells2_this_condition));
            
            n_cells2(iCon, iSize) = length(cells2_this_condition);
        else
            % No cells for this condition - fill with NaN
            avg_stat2{iCon,iSize} = nan(nFrames, 1);
            avg_loc2{iCon,iSize}  = nan(nFrames, 1);
            se_stat2{iCon,iSize}  = nan(nFrames, 1);
            se_loc2{iCon,iSize}   = nan(nFrames, 1);
            n_cells2(iCon, iSize) = 0;
        end
    end
end

% Find global y-axis limits
allMin = cellfun(@(a,s) min(a-s), [avg_stat1(:);avg_loc1(:);avg_stat2(:);avg_loc2(:)], ...
                                  [se_stat1(:);se_loc1(:);se_stat2(:);se_loc2(:)], 'UniformOutput', false);
allMax = cellfun(@(a,s) max(a+s), [avg_stat1(:);avg_loc1(:);avg_stat2(:);avg_loc2(:)], ...
                                  [se_stat1(:);se_loc1(:);se_stat2(:);se_loc2(:)], 'UniformOutput', false);

allMin = cellfun(@(x) min(x), allMin);
allMax = cellfun(@(x) max(x), allMax);

ymin = min(allMin) - 0.1*range([min(allMin) max(allMax)]);
ymax = max(allMax) + 0.1*range([min(allMin) max(allMax)]);
pad  = 0.1 * (ymax - ymin);

t = ((1:nFrames) - (double(r.StimStart) - 1)) / double(r.FrameRate);
z = r.StimDuration;

for iSize = 1:nSize
    figure('Units', 'inches', 'Position', [5+(iSize-1)*0.5, 0, r.FigureSize(1), r.FigureSize(2)]);

    for iCon = 1:nCon
        % Population 1 - stationary vs running
        subplot(nCon, 2, (iCon-1)*2 + 1);
        
        % Plot stationary
        if any(~isnan(avg_stat1{iCon,iSize}))
            shadedErrorBar(t, avg_stat1{iCon,iSize}, se_stat1{iCon,iSize}, r.Colors1{1});
            hold on;
        end
        
        % Plot running
        if any(~isnan(avg_loc1{iCon,iSize}))
            shadedErrorBar(t, avg_loc1{iCon,iSize}, se_loc1{iCon,iSize}, r.Colors1{2});
            hold on;
        end
        
        ylim([ymin ymax]);
        line([0 z], [ymin+0.1*pad, ymin+0.1*pad], 'Color', 'k', 'LineWidth', 2);
        set(gca, 'TickDir', 'out', 'XColor', 'none', 'YColor', 'none', 'Box', 'off'); grid off;
        
        if iCon == 1
            n_range1 = [min(n_cells1(:, iSize)), max(n_cells1(:, iSize))];
            if n_range1(1) == n_range1(2)
                n_str1 = sprintf('n=%d', n_range1(1));
            else
                n_str1 = sprintf('n=%d-%d', n_range1(1), n_range1(2));
            end
            title(sprintf('%s  Size %d  %s', r.Titles{1}, iSize, n_str1));
            line([-1.8 -1.8], [ymin+2*pad, ymin+2*pad+0.05], 'Color', 'k', 'LineWidth', 2);
        end

        % Population 2 - stationary vs running
        subplot(nCon, 2, (iCon-1)*2 + 2);
        
        % Plot stationary
        if any(~isnan(avg_stat2{iCon,iSize}))
            shadedErrorBar(t, avg_stat2{iCon,iSize}, se_stat2{iCon,iSize}, r.Colors2{1});
            hold on;
        end
        
        % Plot running
        if any(~isnan(avg_loc2{iCon,iSize}))
            shadedErrorBar(t, avg_loc2{iCon,iSize}, se_loc2{iCon,iSize}, r.Colors2{2});
            hold on;
        end
        
        ylim([ymin ymax]);
        line([0 z], [ymin+0.1*pad, ymin+0.1*pad], 'Color', 'k', 'LineWidth', 2);
        set(gca, 'TickDir', 'out', 'XColor', 'none', 'YColor', 'none', 'Box', 'off'); grid off;
        
        if iCon == 1
            n_range2 = [min(n_cells2(:, iSize)), max(n_cells2(:, iSize))];
            if n_range2(1) == n_range2(2)
                n_str2 = sprintf('n=%d', n_range2(1));
            else
                n_str2 = sprintf('n=%d-%d', n_range2(1), n_range2(2));
            end
            title(sprintf('%s  Size %d  %s', r.Titles{2}, iSize, n_str2));
            line([-1.8 -1.8], [ymin+2*pad, ymin+2*pad+0.05], 'Color', 'k', 'LineWidth', 2);
        end
    end
end
end
