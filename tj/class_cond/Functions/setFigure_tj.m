function tempFig = setFigure

tempFig = figure; %create figure

ax = gca; %modify axes
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.01 0.01];
ax.FontName = 'Helvetica';
%ax.Box = 'off';

%set(0,'DefaultFigureWindowStyle','docked')
end
