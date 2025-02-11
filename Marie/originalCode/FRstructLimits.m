function FR = FRstructLimits(struct, bwidth, TimeWindow, unit)
% TimeWindow = [tmin, tmax]

n = find([struct.unitID] == unit); %% n changes to index in struct pointing to specified unit
TS2 = [struct(n).timestamps];  %% Make vector, TimeStamps2, that has timestamps from unit.
length(TS2)
title_ = [struct(n).unitID];
title_ = num2str(title_);
titleTr_ = inputname(1);
title_ = strcat('Hz over time is', titleTr_);

TS2 = TS2( (TS2 >= TimeWindow(1) ) & (TS2 <= TimeWindow(2)));

range = TimeWindow;                       % designate x-axis range in sec
L = (length(TS2));                           % L = number of spikes
start = TS2(1);
stop = TS2(end);

FR = length(TS2)/(TimeWindow(2)-TimeWindow(1));

%plot(clusts);
figure;
histogram (TS2, 'BinLimits', range, 'Binwidth', bwidth, 'Facecolor', [0 0 0], 'Linestyle', 'none', 'FaceAlpha',1, 'Normalization', 'countdensity')
%histogram (TS2, 'Binwidth', bwidth, 'Facecolor', [0 0 0], 'Linestyle', 'none', 'FaceAlpha',1, 'Normalization', 'countdensity')
box off;
ax = gca; 
ax.TickDir = 'out';
ax.FontName = 'Calibri'; 'FixedWidth';
ax.FontSize = 18;

%title(title_);
title([num2str(unit),   ' Hz over time']);
end