function FRate = FRstructTimeGridTimeLimit(TimeGridA,TimeGridB, TimeLim, struct, unit)

n = find([struct.unitID] == unit); %% n changes to index in struct pointing to specified unit
TS2 = [struct(n).timestamps];  %% Make vector, TimeStamps2, that has timestamps from unit.
length(TS2)
title_ = [struct(n).unitID];
title_ = num2str(title_);
titleTr_ = inputname(1);
title_ = strcat('Hz over time is', titleTr_);

TS2 = TS2(TS2 < TimeLim(2)); %time limit timestamps
TS2 = TS2(TS2 > TimeLim(1));

TS2 = TimeGridUnit(TimeGridA, TimeGridB, TS2);

TimeTotal = 0;
AllSpikes = 0;

FR = zeros(length(TimeGridB),1);
for f = 1:length(FR)
    TSWin = (TS2(TS2 > TimeGridA(f) & TS2 < TimeGridB(f)));
    AllSpikes = AllSpikes + length(TSWin);
    FR(f) = length(TS2(TS2 > TimeGridA(f) & TS2 < TimeGridB(f))) /(TimeGridB(f)-TimeGridA(f));
    TimeTotal = TimeTotal + (TimeGridB(f)-TimeGridA(f));
end

FRate = AllSpikes/TimeTotal;

%plot(clusts);
figure;
bar (FR, 1,'k');
%histogram (TS2, 'Binwidth', bwidth, 'Facecolor', [0 0 0], 'Linestyle', 'none', 'FaceAlpha',1, 'Normalization', 'countdensity')
box off;
ax = gca; 
ax.TickDir = 'out';
ax.FontName = 'Calibri'; 'FixedWidth';
ax.FontSize = 18;

%title(title_);
title([num2str(unit),   ' Hz over time']);
end