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
%   'Colors' - Colors for all lines. Can be specified as:
%              - Single color string/RGB: applies to all lines (e.g., 'r' or [1 0 0])
%              - Cell array {pre_color, post_color}: applies to both datasets
%              - Cell array {{data1_pre, data1_post}, {data2_pre, data2_post}}: specific colors for each
%              Default: {{'k', 'b'}, {'k', 'b'}}
%   'Colors1' - Colors for pre and post in first dataset [pre_color, post_color]
%               (Overrides 'Colors' for dataset 1 if specified)
%               Default: {'k', 'b'}
%   'Colors2' - Colors for pre and post in second dataset [pre_color, post_color]
%               (Overrides 'Colors' for dataset 2 if specified)
%               Default: {'k', 'b'}
%   'UseDashedLines' - Logical [dash1, dash2] indicating whether to use dashed lines
%                     Default: [false, false]
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
addParameter(p, 'Colors', {}, @(x) validateColors(x));
addParameter(p, 'UseDashedLines', [false, false], @(x) islogical(x) && length(x)==2);
addParameter(p, 'Colors1', {'k', 'b'}, @(x) iscell(x) && length(x)==2);
addParameter(p, 'Colors2', {'k', 'b'}, @(x) iscell(x) && length(x)==2);
addParameter(p, 'XLim', [0, 1.2], @(x) isnumeric(x) && length(x)==2);
addParameter(p, 'Titles', {'Dataset 1', 'Dataset 2'}, @(x) iscell(x) && length(x)==2);
addParameter(p, 'YLabel', 'dF/F', @ischar);
addParameter(p, 'XLabel', '', @ischar);
addParameter(p, 'FigureSize', [6, 4], @(x) isnumeric(x) && length(x)==2);

parse(p, data1, data2, cell_indices1, cell_indices2, contrasts, varargin{:});

% Extract parameters
use_dashed_lines = p.Results.UseDashedLines;
colors_input = p.Results.Colors;
colors1_input = p.Results.Colors1;
colors2_input = p.Results.Colors2;
xlim_range = p.Results.XLim;
titles = p.Results.Titles;
y_label = p.Results.YLabel;
x_label = p.Results.XLabel;
figure_size = p.Results.FigureSize;

% Process color inputs
[colors1, colors2] = processColors(colors_input, colors1_input, colors2_input);

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

% Find global min and max for y-axis scaling
ymin = Inf;
ymax = -Inf;

% Check dataset 1
for id = 1:nd
    % Calculate min considering error bars
    temp_min = min(conResp_data1_avrg{id} - conResp_data1_se{id});
    if ~isnan(temp_min) && temp_min < ymin
        ymin = temp_min;
    end
    
    % Calculate max considering error bars
    temp_max = max(conResp_data1_avrg{id} + conResp_data1_se{id});
    if ~isnan(temp_max) && temp_max > ymax
        ymax = temp_max;
    end
end

% Check dataset 2
for id = 1:nd
    % Calculate min considering error bars
    temp_min = min(conResp_data2_avrg{id} - conResp_data2_se{id});
    if ~isnan(temp_min) && temp_min < ymin
        ymin = temp_min;
    end
    
    % Calculate max considering error bars
    temp_max = max(conResp_data2_avrg{id} + conResp_data2_se{id});
    if ~isnan(temp_max) && temp_max > ymax
        ymax = temp_max;
    end
end

% Add some padding to the y-axis limits (10% of the range)
padding = 0.1 * (ymax - ymin);
ymin = max(0, ymin - padding); % Ensure we don't go below zero if data is all positive
ymax = ymax + padding;

% Create figure
figure('Units', 'inches', 'Position', [5, 5, figure_size(1), figure_size(2)]);

% Plot first dataset
subplot(1, 2, 1);

% Choose line style based on UseDashedLines parameter
line_style1_pre = '-';
line_style1_post = '-';
if use_dashed_lines(1)
    line_style1_pre = '--';
    line_style1_post = '--';
end

% Plot with error bars
errorbar(contrasts, conResp_data1_avrg{pre}, conResp_data1_se{pre}, ...
    [line_style1_pre, 'o'], 'Color', colors1{1}, 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
errorbar(contrasts, conResp_data1_avrg{post}, conResp_data1_se{post}, ...
    [line_style1_post, 'o'], 'Color', colors1{2}, 'LineWidth', 1.5, 'MarkerSize', 6);

% Set axis properties
ylabel(y_label);
xlim(xlim_range);
ylim([ymin, ymax]);
xticks(contrasts);
set(gca, 'TickDir', 'out', 'box', 'off');
grid off;
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
if use_dashed_lines(2)
    line_style2_pre = '--';
    line_style2_post = '--';
end

% Plot with error bars
errorbar(contrasts, conResp_data2_avrg{pre}, conResp_data2_se{pre}, ...
    [line_style2_pre, 'o'], 'Color', colors2{1}, 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
errorbar(contrasts, conResp_data2_avrg{post}, conResp_data2_se{post}, ...
    [line_style2_post, 'o'], 'Color', colors2{2}, 'LineWidth', 1.5, 'MarkerSize', 6);

% Set axis properties
xlim(xlim_range);
ylim([ymin, ymax]);
xticks(contrasts);
set(gca, 'TickDir', 'out', 'box', 'off');
grid off;
if ~isempty(x_label)
    xlabel(x_label);
end

% Add title
title([titles{2}, ' n = ', num2str(length(cell_indices2))], 'FontWeight', 'normal');

end

% Helper function to validate color inputs
function isValid = validateColors(colors)
    isValid = true;
    if isempty(colors)
        return;
    end
    
    % Single color (string or RGB)
    if ischar(colors) || (isnumeric(colors) && length(colors) == 3)
        return;
    end
    
    % Cell array
    if iscell(colors)
        if length(colors) == 2
            % Could be {pre, post} or {{data1_pre, data1_post}, {data2_pre, data2_post}}
            if iscell(colors{1})
                % Nested format
                if length(colors{1}) == 2 && length(colors{2}) == 2
                    return;
                end
            else
                % Simple {pre, post} format
                return;
            end
        end
    end
    
    isValid = false;
end

% Helper function to process color inputs
function [colors1, colors2] = processColors(colors_input, colors1_input, colors2_input)
    % Default colors
    colors1 = colors1_input;
    colors2 = colors2_input;
    
    if isempty(colors_input)
        return;
    end
    
    % Single color for all lines
    if ischar(colors_input) || (isnumeric(colors_input) && length(colors_input) == 3)
        colors1 = {colors_input, colors_input};
        colors2 = {colors_input, colors_input};
        return;
    end
    
    % Cell array format
    if iscell(colors_input)
        if length(colors_input) == 2
            if iscell(colors_input{1})
                % Nested format: {{data1_pre, data1_post}, {data2_pre, data2_post}}
                colors1 = colors_input{1};
                colors2 = colors_input{2};
            else
                % Simple format: {pre_color, post_color} - applies to both datasets
                colors1 = colors_input;
                colors2 = colors_input;
            end
        end
    end
end