function  RotEncodTimeGrid(TimeGridA,TimeGridB, RotEncod)

TS2 = RotEncod;

title_ = 'RotEncod';

TS2 = TimeGridUnit(TimeGridA, TimeGridB, TS2); %select ticks in time grid

TimeTotal = 0;
AllSpikes = 0;

FR = zeros(length(TimeGridB),1);
for f = 1:length(FR)
    TSWin = (TS2(TS2 > TimeGridA(f) & TS2 < TimeGridB(f))); 
    AllSpikes = AllSpikes + length(TSWin);
    FR(f) = length(TS2(TS2 > TimeGridA(f) & TS2 < TimeGridB(f))) /(TimeGridB(f)-TimeGridA(f));
    FR(f) = ((FR(f)/2000)*2*pi*7)/(TimeGridB(1)-TimeGridA(1)); % convert ticks/bin to running rate based on particular metrics of the wheen (2000 ticks/rev, radius = 7 cm)
    TimeTotal = TimeTotal + (TimeGridB(f)-TimeGridA(f));
end

%FRate = AllSpikes/TimeTotal;

%plot(clusts);
figure;
bar (FR, 'k');
%histogram (TS2, 'Binwidth', bwidth, 'Facecolor', [0 0 0], 'Linestyle', 'none', 'FaceAlpha',1, 'Normalization', 'countdensity')
box off;
ax = gca; 
ax.TickDir = 'out';
ax.FontName = 'Calibri'; 'FixedWidth';
ax.FontSize = 18;

%title(title_);
title(title_);
end