function [R, p] = calculateCorrelation(thisCell, otherCells)
    % Check if otherCells is empty or has no valid data
    if isempty(otherCells) || size(otherCells, 2) == 0
        R = NaN;
        p = NaN;
        return;
    end
    
    % Calculate mean of other cells
    otherCellsMean = mean(otherCells, 2, "omitnan");
    
    % Remove NaN values for correlation calculation
    validIdx = ~isnan(otherCellsMean) & ~isnan(thisCell);
    
    if sum(validIdx) > 2  % Need at least 3 data points for correlation
        [corrMatrix, pMatrix] = corrcoef(otherCellsMean(validIdx), thisCell(validIdx));
        R = corrMatrix(1, 2);
        p = pMatrix(1, 2);
    else
        R = NaN;
        p = NaN;
    end
end