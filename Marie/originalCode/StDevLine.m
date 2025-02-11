function [meanLine, stdevLine] = StDevLine(N, edges)
index0 = find(edges == 0)-1;
Prestim = N(1:index0);
stdevLine = std(Prestim);
meanLine = mean(Prestim);

end
