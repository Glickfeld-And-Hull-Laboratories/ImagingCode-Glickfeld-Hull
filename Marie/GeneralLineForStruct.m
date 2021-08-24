% Makes a histogram using a vector of timestamps (trigger) and a structure
% of unit names (.unitID= string, .timestamps= vectors of timestamps for each unit). n is not the unit, but is the index in the structure pointing to one unit. 
% This code is primarily useful with HistLoop to create histograms for all
% the units in the structure.
%
% Similar to OneUnitHistStruct, but instead of the user specifying the unit
% it goes through a whole list.
%
 %TS1 = TimeStampOne, a vector of timestamps that will be used to trigger
% the histogram.
%
% n = index for marching through units.
%
% struct = structure of units (field: unitID, as a string) and associated vectors containing
% timestamps for each unit (field: timstamps).
% xmin, smax, and bwidth for the histogram in seconds
%
% go to folder of interest!!!  & set these variables

function  [meanline, N, edges, stdevLine] = GeneralHistForStruct(TS1, n, struct, xmin, xmax, bwidth, color)  

range = [xmin, xmax];                       % designate x-axis range in sec

trange = abs(xmin);                         % set range to check if spike is in window
if trange < xmax
    trange = xmax;
end


TS2 = [struct(n).timestamps];                   % extract vector of timestamps that corresponds to unit at index n
title_ = [struct(n).unitID];                    % extract name/number of unit at index n (string that will be used in title of histogram)
title_ = num2str(title_);                      
titleTr_ = inputname(1);                         % Extract name of trigger (as a string, to be used in title)
title_ = strcat(titleTr_, ',', '  ', title_);   % make the title (can't get a space to show up after the comma)

L1 = (length(TS1));                           % L1 = number of triggers
L2 = (length(TS2));                             % L2 = number of timestamps

k = 1;                                              %create counter for output vector index
for j = 1:L1                                           % for every trigger
    state = 0;                                      % will check and see if there are any spikes for this trigger
    for i=1:L2                                     % for every element in the spiketime vector
       test = TS2(i,:)- TS1(j,:);             % get difference between spiketimes
       if test == 0
           test= nan;                               % eliminate zeros
       end
       if ((test <= trange) && (test >= -trange))    % Check if difference is in histogram range
           useDeltas (k,1) = test;                  % If yes, add difference to vector that will create histogram
           k = k+1;                                 % update index
           state = 1;                               % confirm spikes for this trigger
       end
           
    end
    if (state == 0)
        useDeltas (k,1) = NaN;              % if there are no spikes for this trigger, add a NaN to mark the trial.
    end
end


if exist ('useDeltas', 'var')                              % Check that spikes exist in the histogram window and create histogram.
    
%figure

%histogram (useDeltas, 'BinLimits', range, 'Binwidth', bwidth, 'Facecolor', color, 'Linestyle', 'none', 'Facealpha', 1)  %'Normalization', 'probability', 
[N, edges] = histcounts(useDeltas, 'BinLimits', range, 'Binwidth', bwidth);
N = FreqLine(N, edges, length(TS1), color);
[meanline, stdevLine] = StDevLine(N, edges);
%AddStDevLines(meanline, stdevLine);
reporter = N;
reporter2 = edges;

xline(0,'b');
title(title_)
box off
%ax.TickDir = 'out'
ax = gca; 
ax.TickDir = 'out';
ax.FontName = 'Times New Roman';
ax.FontSize = 10;

%tiledlayout(flow);
else
    fprintf('%s has no spikes in window\n',title_)  % If no spikes exist in the histogram window, alert user.
    
end


end