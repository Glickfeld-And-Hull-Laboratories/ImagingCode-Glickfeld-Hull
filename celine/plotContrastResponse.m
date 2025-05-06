function plotContrastResponse(data1, data2, cell_indices1, cell_indices2, contrasts, varargin)
% PLOTCONTRASTRESPONSE Plot contrast response curves for different cell populations
% 
% Usage:
%   plotContrastResponse(data1, data2, cell_indices1, cell_indices2, contrasts) 
%   plotContrastResponse(data1, data2, cell_indices1, cell_indices2, contrasts, 'Name', Value, ...)
%
% Inputs:
%   data1 - First dataset (e.g., pref_responses_stat_concat)
%   data2 - Second dataset (can be the same as data1)
%   cell_indices1 - Indices of cells to use for first dataset (e.g., red_ind_concat)
%   cell_indices2 - Indices of cells to use for second dataset (e.g., green_ind_concat)
%   contrasts - Array of contrast values (e.g., cons)
%
% Optional Name-Value Pairs:
%   'UseDashedLines' - Logical [dash1, dash2] indicating whether to use dashed lines
%                     Default: [false, false]
%   'Colors1' - Colors for pre and post in first dataset [pre_color, post_color]
%               Default: {'k', 'b'}
%   'Colors2' - Colors for pre and post in second dataset [pre_color, post_color]
%               Default: {'k', 'b'}
%   'YLim' - Y axis limits [min, max]
%            Default: [0, 0.18]
%   'XLim' - X axis limits [min, max]
%            Default: [0, 1.2]
%   'Titles' - Cell array with titles for each dataset {'title1', 'title2'}
%              Default: {'Dataset 1', 'Dataset 2'}
%   'YLabel' - Label for y-axis
%              Default: 'dF/F'
%   'XLabel' - Label for x-axis (optional)
%              Default: '' (no x-label)
%   'FigureSize' - Size of the figure [width, height] in inches
%                 Default: [6, 4]

% Parse inputs
p = inputParser;
addRequired(p, 'data1', @iscell);
addRequired(p, 'data2', @iscell);
addRequired(p, 'cell_indices1', @isnumeric);
addRequired(p, 'cell_indices2', @isnumeric);
addRequired(p, 'contrasts', @isnumeric);
addParameter(p, 'UseDashedLines', [false, false], @(x) islogical(x) && length(x)==2);
addParameter(p, 'Colors1', {'k', 'b'}, @(x) iscell(x) && length(x)==2);
addParameter(p, 'Colors2', {'k', 'b'}, @(x) iscell(x) && length(x)==2);
addParameter(p, 'YLim', [0, 0.18], @(x) isnumeric(x) && length(x)==2);
addParameter(p, 'XLim', [0, 1.2], @(x) isnumeric(x) && length(x)==2);
addParameter(p, 'Titles', {'Dataset 1', 'Dataset 2'}, @(x) iscell(x) && length(x)==2);
addParameter(p, 'YLabel', 'dF/F', @ischar);
addParameter(p, 'XLabel', '', @ischar);
addParameter(p, 'FigureSize', [6, 4], @(x) isnumeric(x) && length(x)==2);

parse(p, data1, data2, cell_indices1, cell_indices2, contrasts, varargin{:});

% Extract parameters
use_dashed_lines = p.Results.UseDashedLines;
colors1 = p.Results.Colors1;
colors2 = p.Results.Colors2;
ylim_range = p.Results.YLim;
xlim_range = p.Results.XLim;
titles = p.Results.Titles;
y_label = p.Results.YLabel;
x_label = p.Results.XLabel;
figure_size = p.Results.FigureSize;

% Get number of days
nd = size(data1, 2);
pre = 2;    % Index for pre-treatment day
post = 1;   % Index for post-treatment day

% Initialize arrays for average and standard error
conResp_data1_avrg = cell(1, nd);
conResp_data2_avrg = cell(1, nd);
conResp_data1_se = cell(1, nd);
conResp_data2_se = cell(1, nd);

% Calculate average and standard error for each day
for id = 1:nd
    % First dataset
    conResp_data1_avrg{id} = mean(data1{id}(cell_indices1, :), 1, 'omitnan');
    data1_std = std(data1{id}(cell_indices1, :), 1, 'omitnan');
    conResp_data1_se{id} = data1_std / sqrt(length(cell_indices1));
    
    % Second dataset
    conResp_data2_avrg{id} = mean(data2{id}(cell_indices2, :), 1, 'omitnan');
    data2_std = std(data2{id}(cell_indices2, :), 1, 'omitnan');
    conResp_data2_se{id} = data2_std / sqrt(length(cell_indices2));
end

% Create figure
figure('Units', 'inches', 'Position', [5, 5, figure_size(1), figure_size(2)]);

% Plot first dataset
subplot(1, 2, 1);

% Choose line style based on UseDashedLines parameter
line_style1_pre = '-';
line_style1_post = '-';
marker_style1_pre = 'o';
marker_style1_post = 'o';
if use_dashed_lines(1)
    line_style1_pre = '--';
    line_style1_post = '--';
end

% Plot with error bars
errorbar(contrasts, conResp_data1_avrg{pre}, conResp_data1_se{pre}, ...
    [line_style1_pre, colors1{1}], 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
errorbar(contrasts, conResp_data1_avrg{post}, conResp_data1_se{post}, ...
    [line_style1_post, colors1{2}], 'LineWidth', 1.5, 'MarkerSize', 6);

% Set axis properties
ylabel(y_label);
xlim(xlim_range);
ylim(ylim_range);
xticks(contrasts);
set(gca, 'TickDir', 'out', 'box', 'off');
if ~isempty(x_label)
    xlabel(x_label);
end

% Add title
title([titles{1}, ' n = ', num2str(length(cell_indices1))], 'FontWeight', 'normal');

% Plot second dataset
subplot(1, 2, 2);

% Choose line style based on UseDashedLines parameter
line_style2_pre = '-';
line_style2_post = '-';
marker_style2_pre = 'o';
marker_style2_post = 'o';
if use_dashed_lines(2)
    line_style2_pre = '--';
    line_style2_post = '--';
end

% Plot with error bars
errorbar(contrasts, conResp_data2_avrg{pre}, conResp_data2_se{pre}, ...
    [line_style2_pre, colors2{1}], 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
errorbar(contrasts, conResp_data2_avrg{post}, conResp_data2_se{post}, ...
    [line_style2_post, colors2{2}], 'LineWidth', 1.5, 'MarkerSize', 6);

% Set axis properties
xlim(xlim_range);
ylim(ylim_range);
xticks(contrasts);
set(gca, 'TickDir', 'out', 'box', 'off');
if ~isempty(x_label)
    xlabel(x_label);
end

% Add title
title([titles{2}, ' n = ', num2str(length(cell_indices2))], 'FontWeight', 'normal');

end