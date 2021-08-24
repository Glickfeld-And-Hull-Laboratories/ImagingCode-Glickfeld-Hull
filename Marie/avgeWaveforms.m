function AvgWvF = avgeWaveforms(WvTraces)
%title_ = inputname(1);
sz = size(WvTraces);
numrows = sz(1);
numcolumns = sz(1,2);
AvgWvF = (sum(WvTraces, 2))/(numcolumns);

%figure
%plot(AvgWvF, 'r')
%title(title_);
%xticks([0 30 60 90 120 150])
%xticklabels({'0' '.001','.002','.003', '.004', '.005'})
%box off;
%ax.TickDir = 'out'
%ax = gca; 
%ax.TickDir = 'out';
%ax.FontName = 'Calibri'; 'FixedWidth';
%ax.FontSize = 18;
end