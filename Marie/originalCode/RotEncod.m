function  [Avg, reporter] = RotEncod(RotEncod, bwidth) % gives the average ticks/bin histogram for the rotatary encoder data

TS2 = RotEncod; %ticks for rotary encoder

title_ = 'RotEncod';



TimeTotal = 0;
AllSpikes = 0;

EndIter = floor(RotEncod(end)/bwidth); %breaks the data stream into bins, chopping off the end if it doesn't fill a full bin.
for f = 1:(EndIter- 1)
    TSWin = TS2((TS2 <(bwidth*f + bwidth)) & (TS2 > (bwidth *f))); % checks the number of ticks in each bin
    AllSpikes = AllSpikes + length(TSWin); % adds the number of ticks to a total count
    count(f) = length(TSWin); %adds the number of ticks to a vector, FR = each index is the number of ticks in the respective bin
    
    %FR(f) = ((FR(f)/2000)*2*pi*7)/bwidth;
   
    runRate(f) = (((count(f)/2000)*2*pi*7)/bwidth); % convert ticks/bin to running rate based on particular metrics of the wheen (2000 ticks/rev, radius = 7 cm)
    reporter = runRate;
end

Avg = AllSpikes/(EndIter*bwidth); % average running rate

%plot(clusts);
figure;
bar (runRate, 'k');
%histogram (TS2, 'Binwidth', bwidth, 'Facecolor', [0 0 0], 'Linestyle', 'none', 'FaceAlpha',1, 'Normalization', 'countdensity')
box off;
ax = gca; 
ax.TickDir = 'out';
ax.FontName = 'Calibri'; 'FixedWidth';
ax.FontSize = 18;

%title(title_);
title(title_);
end