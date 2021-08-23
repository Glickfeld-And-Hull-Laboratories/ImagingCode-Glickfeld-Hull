% Creates histogram using a vector of triggers and a struct containing unitID
% and timstamps. Modified from iter_crosscor. MEH 2/25/21
%
% go to folder of interest!!!  & call this function
%
%
% TS1 = TimeStampOne, a vector of timestamps that will be used to trigger
% the histogram.
%
% n = unit of interest for the histogram, then index pointing to specified
% unit.
%
% struct = structure of units (field: unitID, as a string) and associated vectors containing
% timestamps for each unit (field: timstamps).
% xmin, smax, and bwidth for the histogram in seconds

function  reporter = OneUnitHistStructEmptyTrials(TS1, n, struct, xmin, xmax, bwidth)

range = [xmin, xmax];                       % designate x-axis range in sec

trange = abs(xmin);                         % define a range that is used to find timestamps for histogram
if trange < xmax
    trange = xmax;
end


n = find([struct.unitID] == n); %% n changes to index in struct pointing to specified unit
TS2 = [struct(n).timestamps];  %% Make vector, TimeStamps2, that has timestamps from unit.
title_ = [struct(n).unitID];
title_ = num2str(title_);
titleTr_ = inputname(1);
title_ = strcat(titleTr_, ',', '  ', title_);

ONE = TS1;
TWO = TS2;

TS2 = ONE;
TS1 = TWO;

L1 = (length(TS1));                           % L1 = number of triggers
L2 = (length(TS2));                            %L2 = number of spikes
%deltaT = tall(zeros(L*L,1));
%column = zeros(100)

reporter = 0;
k = 1;                                              %create counter for output vector index
for  i=1:L2                                           % for every element in the first spiketime vector
     state = 0;                                     % state machine to see if there are any spikes in the window for this trigger
    for j = 1:L1                                  % for every element in the second (or mirror) spiketime vector
       test = TS2(i,:)- TS1(j,:);                  % get difference between spike and trigger
       if test == 0
           test= nan;                               % eliminate zeros
       end
       if ((test <= trange) && (test >= -trange))    % Check if difference is in histogram range
           useDeltas (k,1) = test;                  % If yes, add difference to vector that will create histogram
           k = k + 1;                               % update index
           
           state = 1;                              % yes there are spikes in this window
       end
           
    end
        if state == 1
         reporter = reporter + 1  
         b = i
    
        end
    if (state == 0 )                          % add in a blank trial if there were no hits for this TS1 element
        useDeltas (k,1) = NaN;
        k = k+1;
    end
    
end


if exist ('useDeltas')                              % As long as there are spikes in the range, create the histogram.

figure

histogram (useDeltas, 'BinLimits', range, 'Binwidth', bwidth, 'Facecolor', 'k', 'Linestyle', 'none', 'Facealpha', 1); %, 'Normalization', 'countdensity'
xline(0,'b');
title(title_);
box off;
%ax.TickDir = 'out'
ax = gca; 
ax.TickDir = 'out';
ax.FontName = 'Calibri'; 'FixedWidth';
ax.FontSize = 18;
%length(TS1)
yticklabels(yticks/(length(TS1))/bwidth);


[N, edges] = histcounts(useDeltas, 'BinLimits', range, 'Binwidth', bwidth);
FreqHist(N, edges, length(TS1));
title(title_);
box off;
ax.TickDir = 'out'
ax = gca; 
ax.TickDir = 'out';
ax.FontName = 'Calibri'; 'FixedWidth';
ax.FontSize = 18;
xline(0, 'b');

%tiledlayout(flow);
else
    fprintf('%s has no spikes in window\n',title_) % Alert user if no spikes in range.
    
end


end