% created to make cross-correlograms on data presented as structs. unit1
% and unit2 are particular units in a list of many units. struct contains 2
% fields= .unitID, which is a string identifying which unit, and
% .timestamps, which is a vector of timestamps where that unit fires.
% Adapted from iter_crosscorr
% MEH 3/3/21

function reporter = xCorrStructRunNorun(ifrun, unit1, unit2, struct, xmin, xmax, bwidth, limMin, limMax) %times in seconds, n unit of interest for struct (.unitID= string, .timestamps= vector 

range = [xmin, xmax];                       % designate x-axis range in sec, xmin should be negative
trange = abs(xmin);
if trange < xmax
    trange = xmax;
end

% if you are passing a particular unit, n, instead of using n/for loop to cycle
%through many:
Iunit1 = find([struct.unitID] == unit1);        % find index for units of interest
Iunit2 = find([struct.unitID] == unit2);
Unit1 = [struct(Iunit1).timestamps];         %% Unit 1 is timestamps where unit one fires
Unit2 = [struct(Iunit2).timestamps]; 
title1 = [struct(Iunit1).unitID];           % collect unitIDs as strings for titling the graph
title2 = [struct(Iunit2).unitID];

Unit1 = Unit1((limMin < Unit1) & (Unit1 < limMax)); %limit analysis to time window specified
Unit2 = Unit2((limMin < Unit2) & (Unit2 < limMax));


TS1 = Unit1; % just rename timestamps
TS2 = Unit2;


for i = 1:length(ifrun)-1
    ifrun(i, 3) = ifrun(i,1) + 5; % how much to stay away from trial change time?
    ifrun(i, 4) = ifrun(i+1, 1) - 5;
end



TS2run = [];
TS2stop = [];


for i = 1:length(ifrun)     % make two timestamp vectors, one from the central portion of running trial and one from the central portion of no-run trials for TS2
    TSWin = (TS2(TS2 > ifrun(i,3) & TS2 < ifrun(i,4))); % check if the timestamp is in a section we want to analyze (not right by the change between running and not running)
    % ifrun(i,6) = length(TSWin);
    %AllSpikes = AllSpikes + length(TSWin);
    if ifrun(i,2)==1
        TS2stop = [TS2stop TSWin.'];
    end
    if ifrun(i,2)==2
        TS2run = [TS2run TSWin.'];
    end
   
end

TS1run = [];
TS1stop = [];

for i = 1:length(ifrun)     % and repeat for TS1
    TSWin = (TS1(TS1 > ifrun(i,3) & TS1 < ifrun(i,4))); % check if the timestamp is in a section we want to analyze (not right by the change between running and not running)
    % ifrun(i,6) = length(TSWin);
    %AllSpikes = AllSpikes + length(TSWin);
    if ifrun(i,2)==1
        TS1stop = [TS1stop TSWin.'];
    end
    if ifrun(i,2)==2
        TS1run = [TS1run TSWin.'];
    end
   
end


 
for m =1:2 %run all the cross-correlograms twice, once in the stop condition and once in the run condition
    if m == 1
        Unit1 = TS1stop.'; %The first time we use TSstop
        Unit2 = TS2stop.';
    end
    if m == 2
        Unit1 = TS1run.';
        Unit2 = TS2run.';
    end



L1 = (length(Unit1));            % L = number of spikes
%disp(L1)
L2 = (length(Unit2));

useDeltas = [];                                     %be sure useDeltas is reset after m = 1;
k = 1;                                              %create counter for output vector index
for j = 1:L1                                           % for every element in the first spiketime vector
    for i=1:L2                                    % for every element in the second (or mirror) spiketime vector
       test = Unit2(i,:)- Unit1(j,:);             % get difference between spiketimes
       %if test == 0
       %    test= nan;                               % eliminate zeros
       %end
       if ((test <= trange) && (test >= -trange))    % Check if difference is in histogram range
           useDeltas (k,1) = test;                  % If yes, add difference to vector that will create histogram
           k = k+1;                                 % update index
       end
           
    end  
end

if exist ('useDeltas', 'var')
    
    %strTitle = [num2str(title1) ' vs ' num2str(title2)];
title_ = [num2str(title1) ' & ' num2str(title2) ' from ' num2str(limMin) ' to ' num2str(limMax)];

if m == 1
    color = 'r';
    figure
    
    [Nnorun, edgesNorun] = histcounts(useDeltas, 'BinLimits', range, 'Binwidth', bwidth);
    FreqHist(Nnorun, edgesNorun, length(Unit1), color);
    
    title([title_ ' NoRun']);
    box off;
    %ax.TickDir = 'out'
    ax = gca; 
    ax.TickDir = 'out';
    ax.FontName = 'Calibri';
    ax.FontSize = 18;


end
if m ==2
    color = 'g';
    figure
    
    [Nrun, edgesRun] = histcounts(useDeltas, 'BinLimits', range, 'Binwidth', bwidth);
    FreqHist(Nrun, edgesRun, length(Unit1), color);
   
    title([title_ ' Run']);
     box off;
    %ax.TickDir = 'out'
    ax = gca; 
    ax.TickDir = 'out';
    ax.FontName = 'Calibri';
    ax.FontSize = 18;
  
    
    
    figure %% the second time through, also make a figure that has both graphs
hold on

TranspFreqHist(Nrun, edgesRun, length(TS1run), 1, 'g');
TranspFreqHist(Nnorun, edgesNorun, length(TS1stop), .5, 'r');
reporter = [length(TS1stop); length(TS1run); length(TS2stop); length(TS2run)]; 
 
%figure
%histogram (useDeltas, 'BinLimits', range, 'Binwidth', bwidth, 'Facecolor', 'k', 'Linestyle', 'none', 'Facealpha', 1);
title(title_ );
box off;
%ax.TickDir = 'out'
ax = gca; 
ax.TickDir = 'out';
ax.FontName = 'Calibri';
%ax.FontName = 'FixedWidth';
ax.FontSize = 18;
%length(L1)/bwidth;

hold off

end
    


%tiledlayout(flow);
else
    fprintf('%d has no spikes in window\n',title2)
    
end

end

end
