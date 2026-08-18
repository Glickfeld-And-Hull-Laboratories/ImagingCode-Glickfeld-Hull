function plotCorrelation(otherCellsData, thisCellData, R_value, cellNum, corrType, dayID, pre)
    % Create correlation plot
    otherCellsMean = mean(otherCellsData, 2, "omitnan");
    
    figure;
    % Color based on experimental day
    if dayID == pre
        markerColor = 'black';
    else
        markerColor = 'blue';
    end
    
    scatter(otherCellsMean, thisCellData, 'MarkerFaceColor', markerColor, ...
           'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
    hold on;
    lsline;  % Add regression line
    
    title([corrType ' Cell ' num2str(cellNum) ', R = ' num2str(R_value, '%.3f')]);
    xlabel('Mean of Other Cells');
    ylabel('This Cell');
    grid on;
end
