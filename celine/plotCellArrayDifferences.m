function plotCellArrayDifferences(cellArray)
    % plotCellArrayDifferencesExtracts values from a cell array, skips 
    % the first pair of numeric items, computes differences between the remaining 
    % consecutive elements, and plots the differences.
    %
    % Usage:
    %   plotCellArrayDifferencesSkipFirst(cellArray)
    %
    % Input:
    %   cellArray - A cell array containing numeric values or arrays.

    % Validate input
    if ~iscell(cellArray)
        error('Input must be a cell array.');
    end

    % Extract numeric values from the cell array
    values = [];
    for i = 1:numel(cellArray)
        if isnumeric(cellArray{i})
            values = [values; cellArray{i}(:)];
        else
            warning('Non-numeric element in cell array at index %d ignored.', i);
        end
    end

    % Check if there are enough values for difference calculation
    if numel(values) < 3
        error('Not enough numeric values to compute differences after skipping the first pair.');
    end

    % Skip the first pair
    values = values(3:end);

    % Compute differences between consecutive values
    differences = diff(values);

    % Plot the differences
    figure;
    plot(differences, '-o', 'LineWidth', 1.5);
    xlabel('Index');
    ylabel('Difference');
    %title('Differences Between Consecutive Values (Skipping First Item)');
    grid on;
end
