function tempFig = setFigure

tempFig = figure; %create figure

ax = gca; %modify axes
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.01 0.01];
ax.FontName = 'Helvetica';
ax.Box = 'off';
ax.LineWidth = 1;

set(0,'DefaultFigureWindowStyle','normal')
end