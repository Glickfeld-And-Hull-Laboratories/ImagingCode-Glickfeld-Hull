function h = plotWithErrorBars(data, varargin)
% PLOTWITHERRORBARS Creates a line plot with error bars calculated from raw data
%
% Usage:
%   h = plotWithErrorBars(data)
%   h = plotWithErrorBars(data, 'ErrorType', 'sem')
%   h = plotWithErrorBars(data, 'Dimension', 1, 'XValues', 1:10)
%
% Inputs:
%   data - Raw data matrix/array to plot
%
% Optional Name-Value Pairs:
%   'XValues'     - Vector of x-axis values (default: 1:size(data, otherDim))
%   'ErrorType'   - String, 'sem' (default) or 'std' to specify error bar type
%   'Dimension'   - Dimension to average over (default: 1)
%   'LineColor'   - RGB triplet for line color, default [0 0 0]
%   'ErrorColor'  - RGB triplet for error bar color, default is same as LineColor
%   'LineWidth'   - Width of the line, default 1.5
%   'MarkerSize'  - Size of markers, default 6
%   'MarkerType'  - Type of marker, default 'o'
%   'FigHandle'   - Handle to figure, default is current figure
%   'AxHandle'    - Handle to axes, default is current axes
%
% Output:
%   h - Struct containing handles to the plot elements

% Parse inputs
p = inputParser;
addRequired(p, 'data', @isnumeric);
addParameter(p, 'XValues', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'ErrorType', 'sem', @(x) ismember(lower(x), {'sem', 'std'}));
addParameter(p, 'Dimension', 1, @(x) isnumeric(x) && isscalar(x) && (x == 1 || x == 2));
addParameter(p, 'LineColor', [0 0 0], @(x) isnumeric(x) && length(x) == 3);
addParameter(p, 'ErrorColor', [], @(x) isempty(x) || (isnumeric(x) && length(x) == 3));
addParameter(p, 'LineWidth', 1.5, @isnumeric);
addParameter(p, 'MarkerSize', 6, @isnumeric);
addParameter(p, 'MarkerType', 'o', @ischar);
addParameter(p, 'FigHandle', [], @(x) isempty(x) || ishandle(x));
addParameter(p, 'AxHandle', [], @(x) isempty(x) || ishandle(x));

parse(p, data, varargin{:});
opts = p.Results;

% Setup figure and axes if not provided
if isempty(opts.FigHandle)
    opts.FigHandle = gcf;
end

if isempty(opts.AxHandle)
    opts.AxHandle = gca;
end

% Set error color same as line color if not specified
if isempty(opts.ErrorColor)
    opts.ErrorColor = opts.LineColor;
end

% Determine averaging dimension and other dimension
dim = opts.Dimension;
otherDim = 3 - dim; % If dim=1, otherDim=2; if dim=2, otherDim=1

% Create x values if not provided
if isempty(opts.XValues)
    x = 1:size(data, otherDim);
else
    x = opts.XValues;
    if length(x) ~= size(data, otherDim)
        error('Length of XValues must match size(data, %d)', otherDim);
    end
end

% Calculate mean and error
y_mean = mean(data, dim);
if strcmpi(opts.ErrorType, 'sem')
    % Standard error of the mean
    err = std(data, 0, dim) ./ sqrt(size(data, dim));
else
    % Standard deviation
    err = std(data, 0, dim);
end

% Ensure proper orientation for plotting
if dim == 2
    y_mean = y_mean(:)';
    err = err(:)';
end

% Create the plot
hold(opts.AxHandle, 'on');

% Plot the line
h.line = plot(opts.AxHandle, x, y_mean, ...
    'Color', opts.LineColor, ...
    'LineWidth', opts.LineWidth, ...
    'Marker', opts.MarkerType, ...
    'MarkerSize', opts.MarkerSize, ...
    'MarkerFaceColor', opts.LineColor, ...
    'MarkerEdgeColor', opts.LineColor);

% Add error bars
h.errorbar = errorbar(opts.AxHandle, x, y_mean, err, 'LineStyle', 'none', ...
    'Color', opts.ErrorColor, ...
    'LineWidth', opts.LineWidth/2);

% Set the properties according to preferences
set(opts.AxHandle, 'TickDir', 'out');
grid(opts.AxHandle, 'off');
box(opts.AxHandle, 'off');

% Add buffer to axes to prevent error bars from overlapping with axes
ax_xlim = get(opts.AxHandle, 'XLim');
ax_ylim = get(opts.AxHandle, 'YLim');
x_range = max(x) - min(x);
y_range = max(y_mean + err) - min(y_mean - err);
x_buffer = 0.05 * x_range;
y_buffer = 0.05 * y_range;

% Set new limits with buffer
if x_range > 0
    set(opts.AxHandle, 'XLim', [min(x) - x_buffer, max(x) + x_buffer]);
end
if y_range > 0
    set(opts.AxHandle, 'YLim', [min(y_mean - err) - y_buffer, max(y_mean + err) + y_buffer]);
end

% Set the figure current
figure(opts.FigHandle);

end