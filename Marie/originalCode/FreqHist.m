function N = FreqHist(N, edges, L1, color)

% call like this
% [N, edges] = histcounts(useDeltas, 'BinLimits', range, 'Binwidth', bwidth);
% FreqHist(N, edges, length(TS1), 'k');

bwidth = edges(2)-edges(1);
N = (N/L1)*(1/bwidth); % convert bars to firing rates
edges = edges(1:(length(edges)-1)); % remove last traiiling edge so sizes of N and edges match)
edges = edges + (.5*bwidth); % shift edges so bars are aligned to left side of each bin.

%hold off
%figure
%hold on
bar (edges, N, 1, color);



end