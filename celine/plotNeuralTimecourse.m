function plotNeuralTimecourse(data1, data2, cell_indices1, cell_indices2, varargin)
% PLOTNEURALTIMECOURSE Plot timecourses of neural activity for different stimulus conditions
% 
% Usage:
%   plotNeuralTimecourse(data1, data2, cell_indices1, cell_indices2) 
%   plotNeuralTimecourse(data1, data2, cell_indices1, cell_indices2, 'Name', Value, ...)
%
% Inputs:
%   data1 - First dataset (e.g., tc_trial_avrg_stat_concat)
%   data2 - Second dataset (can be the same as data1)
%   cell_indices1 - Indices of cells to use for first dataset (e.g., red_ind_concat)
%   cell_indices2 - Indices of cells to use for second dataset (e.g., green_ind_concat)
%
% Optional Name-Value Pairs:
%   'UseDashedLines' - Logical vector [dash1, dash2] indicating whether to use dashed lines
%                       Default: [false, false]
%   'Colors1' - Colors for pre and post in first dataset [pre_color, post_color]
%                Default: {'k', 'b'}
%   'Colors2' - Colors for pre and post in second dataset [pre_color, post_color]
%                Default: {'k', 'b'}
%   'Titles' - Cell array with titles for each dataset {'title1', 'title2'}
%                Default: {'Dataset 1', 'Dataset 2'}
%   'StimDuration' - Duration of stimulus in seconds
%                Default: 2
%   'FrameRate' - Frame rate in Hz
%                Default: 15
%   'StimStart' - Frame at which stimulus starts
%                Default: []
%   'FigureSize' - Size of the figure [width, height] in inches
%                Default: [4, 9]

% Parse inputs
p = inputParser;
addRequired(p, 'data1', @iscell);
addRequired(p, 'data2', @iscell);
addRequired(p, 'cell_indices1', @isnumeric);
addRequired(p, 'cell_indices2', @isnumeric);
addParameter(p, 'UseDashedLines', [false, false], @(x) islogical(x) && length(x)==2);
addParameter(p, 'Colors1', {'k', 'b'}, @(x) iscell(x) && length(x)==2);
addParameter(p, 'Colors2', {'k', 'b'}, @(x) iscell(x) && length(x)==2);
addParameter(p, 'Titles', {'Dataset 1', 'Dataset 2'}, @(x) iscell(x) && length(x)==2);
addParameter(p, 'StimDuration', 2, @isnumeric);
addParameter(p, 'FrameRate', 15, @isnumeric);
addParameter(p, 'StimStart', [], @isnumeric);
addParameter(p, 'FigureSize', [4, 9], @(x) isnumeric(x) && length(x)==2);

parse(p, data1, data2, cell_indices1, cell_indices2, varargin{:});

% Extract parameters
use_dashed_lines = p.Results.UseDashedLines;
colors1 = p.Results.Colors1;
colors2 = p.Results.Colors2;
titles = p.Results.Titles;
stim_duration = p.Results.StimDuration;
frame_rate = p.Results.FrameRate;
stimStart = p.Results.StimStart;
figure_size = p.Results.FigureSize;

% Get parameters from the data
nd = size(data1, 2);  % Number of days
pre = 2;              % Index for pre-treatment day
post = 1;             % Index for post-treatment day
nCon = size(data1{pre}, 3);  % Number of contrasts
nSizes = size(data1{pre}, 4);  % Number of sizes

% If stimStart is not provided, use a reasonable default
if isempty(stimStart)
    stimStart = round(size(data1{pre}, 1) / 3);
end

% Find global min and max for y-axis scaling across ALL sizes and conditions
ymin_global = Inf;
ymax_global = -Inf;

for iSize = 1:nSizes
    % Initialize arrays for average and standard error
    tc_data1_avrg = cell(1, nd);
    tc_data2_avrg = cell(1, nd);
    tc_data1_se = cell(1, nd);
    tc_data2_se = cell(1, nd);
    
    % Calculate average and standard error for each day and condition
    for id = 1:nd
        for iCon = 1:nCon
            % First dataset
            tc_data1_avrg{id}(:, iCon) = nanmean(data1{id}(:, cell_indices1, iCon, iSize), 2);
            data1_std = nanstd(data1{id}(:, cell_indices1, iCon, iSize), [], 2);
            tc_data1_se{id}(:, iCon) = data1_std / sqrt(length(cell_indices1));
            
            % Second dataset
            tc_data2_avrg{id}(:, iCon) = nanmean(data2{id}(:, cell_indices2, iCon, iSize), 2);
            data2_std = nanstd(data2{id}(:, cell_indices2, iCon, iSize), [], 2);
            tc_data2_se{id}(:, iCon) = data2_std / sqrt(length(cell_indices2));
        end
    end
    
    % Update global min/max for this size
    for id = 1:nd
        for iCon = 1:nCon
            % Check dataset 1
            temp_min = min(tc_data1_avrg{id}(:, iCon) - tc_data1_se{id}(:, iCon));
            if ~isnan(temp_min) && temp_min < ymin_global
                ymin_global = temp_min;
            end
            temp_max = max(tc_data1_avrg{id}(:, iCon) + tc_data1_se{id}(:, iCon));
            if ~isnan(temp_max) && temp_max > ymax_global
                ymax_global = temp_max;
            end
            
            % Check dataset 2
            temp_min = min(tc_data2_avrg{id}(:, iCon) - tc_data2_se{id}(:, iCon));
            if ~isnan(temp_min) && temp_min < ymin_global
                ymin_global = temp_min;
            end
            temp_max = max(tc_data2_avrg{id}(:, iCon) + tc_data2_se{id}(:, iCon));
            if ~isnan(temp_max) && temp_max > ymax_global
                ymax_global = temp_max;
            end
        end
    end
end

% Add padding to the global y-axis limits (10% of the range)
padding = 0.1 * (ymax_global - ymin_global);
ymin_final = ymin_global - padding;
ymax_final = ymax_global + padding;

% Now plot each size with the same y-axis limits
for iSize = 1:nSizes
    % Recalculate data for this size (since we need it for plotting)
    tc_data1_avrg = cell(1, nd);
    tc_data2_avrg = cell(1, nd);
    tc_data1_se = cell(1, nd);
    tc_data2_se = cell(1, nd);
    
    for id = 1:nd
        for iCon = 1:nCon
            % First dataset
            tc_data1_avrg{id}(:, iCon) = nanmean(data1{id}(:, cell_indices1, iCon, iSize), 2);
            data1_std = nanstd(data1{id}(:, cell_indices1, iCon, iSize), [], 2);
            tc_data1_se{id}(:, iCon) = data1_std / sqrt(length(cell_indices1));
            
            % Second dataset
            tc_data2_avrg{id}(:, iCon) = nanmean(data2{id}(:, cell_indices2, iCon, iSize), 2);
            data2_std = nanstd(data2{id}(:, cell_indices2, iCon, iSize), [], 2);
            tc_data2_se{id}(:, iCon) = data2_std / sqrt(length(cell_indices2));
        end
    end
    
    % Create time axis in seconds
    t = 1:(size(tc_data1_avrg{1}, 1));
    t = (t - (double(stimStart) - 1)) / double(frame_rate);
    
    % Create a stimulus marker
    z = stim_duration;
    
    % Create figure
    figure('Units', 'inches', 'Position', [5 + (iSize-1)*0.5, 0, figure_size(1), figure_size(2)]);
    
    % Define subplot positions
    positions = [1, 2; 3, 4; 5, 6];
    
    % Plot each condition
    for iCon = 1:nCon
        p1 = positions(iCon, 1);
        p2 = positions(iCon, 2);
        
        % Plot first dataset
        subplot(3, 2, p1);
        
        % Choose line style based on UseDashedLines parameter
        line_style1_pre = '';
        line_style1_post = '';
        if use_dashed_lines(1)
            line_style1_pre = '--';
            line_style1_post = '--';
        end
        
        % Plot with shaded error bars
        shadedErrorBar(t, tc_data1_avrg{pre}(:, iCon), tc_data1_se{pre}(:, iCon), [line_style1_pre, colors1{1}]);
        hold on;
        shadedErrorBar(t, tc_data1_avrg{post}(:, iCon), tc_data1_se{post}(:, iCon), [line_style1_post, colors1{2}]);
        
        % Apply global y-axis limits to all plots
        ylim([ymin_final, ymax_final]);
        hold on;
        line([0, z], [ymin_final + 0.1*padding, ymin_final + 0.1*padding], 'Color', 'black', 'LineWidth', 2);
        
        if iCon == 1
            title([titles{1}, ' Size ', num2str(iSize), ' n = ', num2str(length(cell_indices1))]);
            % Vertical calibration bar: 5% df/f
            y_bar_height = 0.05;
            y_bar_start = ymin_final + 2*padding;
            line([-1.8, -1.8], [y_bar_start, y_bar_start + y_bar_height], 'Color', 'black', 'LineWidth', 2);
        end
        set(gca, 'TickDir', 'out', 'XColor', 'none', 'YColor', 'none', 'box', 'off');
        grid off;
        
        % Plot second dataset
        subplot(3, 2, p2);
        
        % Choose line style based on UseDashedLines parameter
        line_style2_pre = '-';
        line_style2_post = '-';
        if use_dashed_lines(2)
            line_style2_pre = '--';
            line_style2_post = '--';
        end
        
        % Plot with shaded error bars
        shadedErrorBar(t, tc_data2_avrg{pre}(:, iCon), tc_data2_se{pre}(:, iCon), [line_style2_pre, colors2{1}] );
        hold on;
        shadedErrorBar(t, tc_data2_avrg{post}(:, iCon), tc_data2_se{post}(:, iCon), [line_style2_post,colors2{2}]);
        
        % Apply global y-axis limits to all plots
        ylim([ymin_final, ymax_final]);
        hold on;
        line([0, z], [ymin_final + 0.1*padding, ymin_final + 0.1*padding], 'Color', 'black', 'LineWidth', 2);
        
        if iCon == 1
            title([titles{2}, ' Size ', num2str(iSize), ' n = ', num2str(length(cell_indices2))]);
            % Vertical calibration bar: 5% df/f
            y_bar_height = 0.05;
            y_bar_start = ymin_final + 2*padding;
            line([-1.8, -1.8], [y_bar_start, y_bar_start + y_bar_height], 'Color', 'black', 'LineWidth', 2);
        end
        set(gca, 'TickDir', 'out', 'XColor', 'none', 'YColor', 'none', 'box', 'off');
        grid off;
    end
end
end