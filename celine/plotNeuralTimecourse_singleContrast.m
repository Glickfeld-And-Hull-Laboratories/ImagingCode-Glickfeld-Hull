function plotNeuralTimecourse_singleContrast(data1, data2, cell_indices1, cell_indices2, contrastIndex, varargin)
% PLOTNEURALTIMECOURSE_SINGLECONTRAST Plot timecourses of neural activity for a single contrast
% 
% Usage:
%   plotNeuralTimecourse_singleContrast(data1, data2, cell_indices1, cell_indices2, contrastIndex)
%   plotNeuralTimecourse_singleContrast(data1, data2, cell_indices1, cell_indices2, contrastIndex, 'Name', Value, ...)
%
% Inputs:
%   data1 - First dataset (e.g., tc_trial_avrg_stat_concat)
%   data2 - Second dataset (can be the same as data1)
%   cell_indices1 - Indices of cells to use for first dataset (e.g., red_ind_concat)
%   cell_indices2 - Indices of cells to use for second dataset (e.g., green_ind_concat)
%   contrastIndex - Index of the contrast to plot (e.g., 1, 2, or 3)
%
% Optional Name-Value Pairs:
%   'UseDashedLines' - Logical vector [dash1, dash2] indicating whether to use dashed lines
%                       Default: [false, false]
%   'Colors1' - Colors for pre and post in first dataset [pre_color, post_color]
%                Default: {'k', 'b'}
%   'Colors2' - Colors for pre and post in second dataset [pre_color, post_color]
%                Default: {'k', 'b'}
%   'YLim1' - Y axis limits for first dataset [min, max]
%                Default: [-0.02, 0.17]
%   'YLim2' - Y axis limits for second dataset [min, max]
%                Default: [-0.02, 0.25]
%   'Titles' - Cell array with titles for each dataset {'title1', 'title2'}
%                Default: {'Dataset 1', 'Dataset 2'}
%   'StimDuration' - Duration of stimulus in seconds
%                Default: 2
%   'FrameRate' - Frame rate in Hz
%                Default: 15
%   'StimStart' - Frame at which stimulus starts
%                Default: []
%   'FigureSize' - Size of the figure [width, height] in inches
%                Default: [4, 3]

% Parse inputs
p = inputParser;
addRequired(p, 'data1', @iscell);
addRequired(p, 'data2', @iscell);
addRequired(p, 'cell_indices1', @isnumeric);
addRequired(p, 'cell_indices2', @isnumeric);
addRequired(p, 'contrastIndex', @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'UseDashedLines', [false, false], @(x) islogical(x) && length(x)==2);
addParameter(p, 'Colors1', {'k', 'b'}, @(x) iscell(x) && length(x)==2);
addParameter(p, 'Colors2', {'k', 'b'}, @(x) iscell(x) && length(x)==2);
addParameter(p, 'YLim1', [-0.02, 0.17], @(x) isnumeric(x) && length(x)==2);
addParameter(p, 'YLim2', [-0.02, 0.25], @(x) isnumeric(x) && length(x)==2);
addParameter(p, 'Titles', {'Dataset 1', 'Dataset 2'}, @(x) iscell(x) && length(x)==2);
addParameter(p, 'StimDuration', 2, @isnumeric);
addParameter(p, 'FrameRate', 15, @isnumeric);
addParameter(p, 'StimStart', [], @isnumeric);
addParameter(p, 'FigureSize', [4, 3], @(x) isnumeric(x) && length(x)==2);
addParameter(p, 'DayOrder', 'r', @(x) ischar(x) && (strcmp(x,'f') || strcmp(x,'r')));

parse(p, data1, data2, cell_indices1, cell_indices2, contrastIndex, varargin{:});

% Extract parameters
use_dashed_lines = p.Results.UseDashedLines;
colors1 = p.Results.Colors1;
colors2 = p.Results.Colors2;
ylim1 = p.Results.YLim1;
ylim2 = p.Results.YLim2;
titles = p.Results.Titles;
stim_duration = p.Results.StimDuration;
frame_rate = p.Results.FrameRate;
stimStart = p.Results.StimStart;
figure_size = p.Results.FigureSize;
iCon = contrastIndex; % Using the direct variable instead of p.Results

% Get parameters from the data
nd = size(data1, 2);  % Number of days

% Set pre and post indices based on DayOrder
if strcmp(p.Results.DayOrder, 'f')
    pre = 1;
    post = 2;
else
    pre = 2;
    post = 1;
end

nCon = size(data1{pre}, 3);  % Number of conditions (contrasts)

% Check if contrast index is valid
if iCon > nCon
    error('Contrast index exceeds the number of available contrasts (%d)', nCon);
end

% If stimStart is not provided, use a reasonable default
if isempty(stimStart)
    stimStart = round(size(data1{pre}, 1) / 3);
end

% Initialize arrays for average and standard error
tc_data1_avrg = cell(1, nd);
tc_data2_avrg = cell(1, nd);
tc_data1_se = cell(1, nd);
tc_data2_se = cell(1, nd);

% Calculate average and standard error for each day for the selected contrast
for id = 1:nd
    % First dataset
    tc_data1_avrg{id} = nanmean(data1{id}(:, cell_indices1, iCon), 2);
    data1_std = nanstd(data1{id}(:, cell_indices1, iCon), [], 2);
    tc_data1_se{id} = data1_std / sqrt(length(cell_indices1));
    
    % Second dataset
    tc_data2_avrg{id} = nanmean(data2{id}(:, cell_indices2, iCon), 2);
    data2_std = nanstd(data2{id}(:, cell_indices2, iCon), [], 2);
    tc_data2_se{id} = data2_std / sqrt(length(cell_indices2));
end

% Create time axis in seconds
t = 1:(size(tc_data1_avrg{1}, 1));
t = (t - (double(stimStart) - 1)) / double(frame_rate);

% Create a stimulus marker
z = stim_duration;

% Create figure
figure('Units', 'inches', 'Position', [5, 0, figure_size(1), figure_size(2)]);

% Plot in a 1x2 layout (side by side)
% Plot first dataset
subplot(1, 2, 1);

% Choose line style based on UseDashedLines parameter
line_style1_pre = '';
line_style1_post = '';
if use_dashed_lines(1)
    line_style1_pre = '--';
    line_style1_post = '--';
end

% Plot with shaded error bars
shadedErrorBar(t, tc_data1_avrg{pre}, tc_data1_se{pre}, [line_style1_pre, colors1{1}]);
hold on;
shadedErrorBar(t, tc_data1_avrg{post}, tc_data1_se{post}, [line_style1_post, colors1{2}]);

ylim(ylim1);
hold on;
line([0, z], [-0.01, -0.01], 'Color', 'black', 'LineWidth', 2);
title([titles{1}, ' n = ', num2str(length(cell_indices1))]);
line([-1.8, -1.8], [0.01, 0.06], 'Color', 'black', 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'box', 'off');
grid off;
axis off
% Plot second dataset
subplot(1, 2, 2);

% Choose line style based on UseDashedLines parameter
line_style2_pre = '-';
line_style2_post = '-';
if use_dashed_lines(2)
    line_style2_pre = '--';
    line_style2_post = '--';
end

% Plot with shaded error bars
shadedErrorBar(t, tc_data2_avrg{pre}, tc_data2_se{pre}, [line_style2_pre, colors2{1}]);
hold on;
shadedErrorBar(t, tc_data2_avrg{post}, tc_data2_se{post}, [line_style2_post, colors2{2}]);

ylim(ylim2);
hold on;
line([0, z], [-0.01, -0.01], 'Color', 'black', 'LineWidth', 2);
title([titles{2}, ' n = ', num2str(length(cell_indices2))]);
line([-1.8, -1.8], [0.01, 0.06], 'Color', 'black', 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'box', 'off');
grid off;
axis off


end