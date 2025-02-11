function meanLine = AddStDevLines(meanLine, stdevLine)
%pass N (histogram counts) and edges (with trailing edge already removed to
%match counts)
%with hold on on figure of interest
%plots green and yellow lines showing stdevs away from mean


UpLine = meanLine + 2*stdevLine;
yline(UpLine, 'g');
DownLine = meanLine - 2*stdevLine;
if DownLine >0
yline(DownLine, 'y');
end
end