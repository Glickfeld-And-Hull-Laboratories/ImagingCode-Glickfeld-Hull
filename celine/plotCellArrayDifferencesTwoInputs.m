function plotCellArrayDifferencesTwoInputs(cellArray1, cellArray2)
    % plotCellArrayDifferencesSkipTwo Processes two cell arrays, skips the 
    % first two numeric values in each, computes differences, converts the 
    % first cell array from microseconds to milliseconds, and plots the results.
    %
    % Usage:
    %   plotCellArrayDifferencesSkipTwo(cellArray1, cellArray2)
    %
    % Inputs:
    %   cellArray1 - A cell array with numeric values (in microseconds)
    %   cellArray2 - A cell array with numeric values

    % Validate inputs
    if ~iscell(cellArray1) || ~iscell(cellArray2)
        error('Both inputs must be cell arrays.');
    end

    % Extract numeric values from cellArray1 and convert to milliseconds
    values1 = [];
    for i = 1:numel(cellArray1)
        if isnumeric(cellArray1{i})
            values1 = [values1; cellArray1{i}(:)];
        else
            warning('Non-numeric element in cellArray1 at index %d ignored.', i);
        end
    end
    values1 = values1 / 1000; % Convert to milliseconds

    % Extract numeric values from cellArray2
    values2 = [];
    for i = 1:numel(cellArray2)
        if isnumeric(cellArray2{i})
            values2 = [values2; cellArray2{i}(:)];
        else
            warning('Non-numeric element in cellArray2 at index %d ignored.', i);
        end
    end

    % Ensure there are enough values for difference calculation
    if numel(values1) < 4 || numel(values2) < 4
        error('Both cell arrays must have enough numeric values to compute differences after skipping the first two items.');
    end

    % Skip the first two numeric values in each array
    values1 = values1(1:end);
    values2 = values2(1:end);

    % Compute differences
    differences1 = diff(values1);
    differences2 = diff(values2);

    % Plot the differences
    figure;
    subplot(2, 1, 1);
    plot(differences1, 'LineWidth', .25);
    xlabel('Index');
    ylabel('Difference (ms)');
    title('Differences counterTimes');
    %hline(round((1/frameRate)*1000),'r') %to add a refernce line at the
    %expected value, based on frame rate
    %ylim([20 120])
    grid on;

    subplot(2, 1, 2);
    plot(differences2,  'LineWidth', .25);
    xlabel('Index');
    ylabel('Difference');
    title('Differences counterValues');
    grid on;
    %ylim([-2 2])
end
