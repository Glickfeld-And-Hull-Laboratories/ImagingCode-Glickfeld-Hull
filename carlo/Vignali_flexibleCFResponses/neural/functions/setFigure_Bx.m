function tempFig = setFigure_Bx

tempFig = figure; %create figure

ax = gca; %modify axes
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.01 0.01];
ax.FontName = 'Helvetica';
ax.YLim = [0 10];
ax.XLim = [-1 4];
ax.Box = 'off';
%ax.Position = [.1 .1 .8 .4];
end