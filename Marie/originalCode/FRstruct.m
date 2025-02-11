function FRstruct(unit, struct, bwidth)

n = find([struct.unitID] == unit); %% n changes to index in struct pointing to specified unit
TS2 = [struct(n).timestamps];  %% Make vector, TimeStamps2, that has timestamps from unit.
length(TS2)
title_ = [struct(n).unitID];
title_ = num2str(title_);
titleTr_ = inputname(1);
title_ = strcat('Hz over time is', titleTr_);

range = [0, TS2(end)];                       % designate x-axis range in sec
L = (length(TS2));                           % L = number of spikes
start = TS2(1);
stop = TS2(end);

%plot(clusts);
%figure;
%histogram (TS2, 'BinLimits', range, 'Binwidth', bwidth, 'Facecolor', [0 0 0], 'Linestyle', 'none', 'FaceAlpha',.5, 'Normalization', 'probability')
[N, edges] = histcounts(TS2, 'BinLimits', range, 'Binwidth', bwidth);

bwidth = edges(2)-edges(1);
N = (N)*(1/bwidth); % convert bars to firing rates
edges = edges(1:(length(edges)-1)); % remove last traiiling edge so sizes of N and edges match)
edges = edges + (.5*bwidth); % shift edges so bars are aligned to left side of each bin.

%hold off
%figure
%hold on
bar (edges, N, 1, 'k');


%histogram (TS2, 'Binwidth', bwidth, 'Facecolor', [0 0 0], 'Linestyle', 'none', 'FaceAlpha',1, 'Normalization', 'countdensity')
box off;
ax = gca; 
ax.TickDir = 'out';
ax.FontName = 'Calibri'; 'FixedWidth';
ax.FontSize = 18;

%title(title_);
title([num2str(unit),   ' Hz over time']);
end