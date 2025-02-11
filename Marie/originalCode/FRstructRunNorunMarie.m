function [FRate, FRrun, FRstop, reporter] = FRstructRunNorunMarie(ifrun, struct, burstLim, unit)

% takes ifrun, a 2-column matrix where the first column is the start time
% of every trial, starting at 0. The second column has a 1 if the trial is
% stationary and a 2 if the trial is locomotion. Written to analyze
% Shuyangs running data. 5/17/21 MEH


n = find([struct.unitID] == unit); %% n changes to index in struct pointing to specified unit
TS2 = [struct(n).timestamps];  %% Make vector, TimeStamps2, that has timestamps from unit.
length(TS2)
title_ = [struct(n).unitID];
title_ = num2str(title_);
titleTr_ = inputname(1);
title_ = strcat('Hz over time is', titleTr_);

TS2 = RemoveBurst(TS2, burstLim);

for i = 1:length(ifrun)-1
    ifrun(i, 3) = ifrun(i,1) + 5; % how much to stay away from trial change time?
    ifrun(i, 4) = ifrun(i+1, 1) - 5;
end


TimeTotal = 0;
AllSpikes = 0;
AllSpikesRun = 0;
TimeRun = 0;
AllSpikesStop = 0;
TimeStop = 0;

%FR = zeros(length(TimeGridB),1);
for i = 1:length(ifrun)
    TSWin = (TS2(TS2 > ifrun(i,3) & TS2 < ifrun(i,4)));
    ifrun(i,6) = length(TSWin);
    AllSpikes = AllSpikes + length(TSWin);
    if ifrun(i,2)==1
        AllSpikesStop = AllSpikesStop + length(TSWin);
        TimeStop = TimeStop + (ifrun(i,4) - ifrun(i,3));
        ifrun(i,7) = AllSpikesStop;
    end
    if ifrun(i,2)==2
        AllSpikesRun = AllSpikesRun + length(TSWin);
        TimeRun = TimeRun + (ifrun(i,4) - ifrun(i,3));
    end
    ifrun(i,5) =  length(TSWin) /(ifrun(i,4) - ifrun(i,3));
    TimeTotal = TimeTotal + (ifrun(i,4) - ifrun(i,3));
end

reporter = ifrun;
FRate = AllSpikes/TimeTotal;
FRrun = AllSpikesRun/TimeRun;
FRstop = AllSpikesStop/TimeStop;
%FR = ifrun(:,5);

allspikesrunSpikes = (AllSpikesRun)
allspikesstopSpikes = (AllSpikesStop)

ifrunRUN = ifrun((ifrun(:,2) == 2),:);
ifrunSTOP = ifrun((ifrun(:,2) == 1),:);


%plot(clusts);
figure;
bar (ifrunRUN(:,1), ifrunRUN(:,5), 1,'g');
hold on
bar (ifrunSTOP(:,1), ifrunSTOP(:,5), 1, 'r');

%histogram (TS2, 'Binwidth', bwidth, 'Facecolor', [0 0 0], 'Linestyle', 'none', 'FaceAlpha',1, 'Normalization', 'countdensity')
box off;
ax = gca; 
ax.TickDir = 'out';
ax.FontName = 'Calibri'; 'FixedWidth';
ax.FontSize = 18;

%title(title_);
title([num2str(unit),   ' Hz over time']);
end