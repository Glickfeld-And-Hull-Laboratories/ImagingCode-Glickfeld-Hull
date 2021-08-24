% created to make cross-correlograms on data presented as structs. unit1
% and unit2 are particular units in a list of many units. struct contains 2
% fields= .unitID, which is a string identifying which unit, and
% .timestamps, which is a vector of timestamps where that unit fires.
% Adapted from iter_crosscorr
% MEH 3/3/21

function reporter = xCorrTimeGrid(TimeGridA, TimeGridB, struct, xmin, xmax, bwidth, unit1, unit2) %times in seconds, n unit of interest for struct (.unitID= string, .timestamps= vector 

range = [xmin, xmax];                       % designate x-axis range in sec, xmin should be negative
trange = abs(xmin);
if trange < xmax;
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

Unit1 = TimeGridUnit(TimeGridA, TimeGridB, Unit1);
Unit2 = TimeGridUnit(TimeGridA, TimeGridB, Unit2);

%counts = xcorrBusiness(Unit1, Unit2, .001, -.02, .02);

L1 = (length(Unit1));           % L = number of spikes
%disp(L1)
L2 = (length(Unit2));

k = 1;                              %create counter for output vector index
TO = 0;
for j = 1:L1                                           % for every element in the first spiketime vector
    for i=1:L2                                    % for every element in the second (or mirror) spiketime vector
       test = Unit2(i,:)- Unit1(j,:);             % get difference between spiketimes
       if test == 0;
           test= nan;                               % eliminate zeros
           
       end
       if ((test <= trange) && (test >= -trange));    % Check if difference is in histogram range
           useDeltas (k,1) = test;                  % If yes, add difference to vector that will create histogram
           k = k+1;                                 % update index
       end
           
    end  
end

if exist ('useDeltas', 'var')
%reporter = useDeltas    
% figure
%reporter = length(useDeltas);
% histogram (useDeltas, 'BinLimits', range, 'Binwidth', bwidth, 'Facecolor', 'k', 'Linestyle', 'none', 'Facealpha', 1); %'Normalization', 'Probability'
 
% strTitle = [num2str(title1) ' vs ' num2str(title2)];
% title([num2str(title1) ' & ' num2str(title2)]);
% box off;
% %ax.TickDir = 'out'
% ax = gca; 
% ax.TickDir = 'out';
% ax.FontName = 'Calibri';
% %ax.FontName = 'FixedWidth';
% ax.FontSize = 18;
% %length(L1)/bwidth;
% %yticklabels(yticks/L1/bwidth);


[N, edges] = histcounts(useDeltas, 'BinLimits', range, 'Binwidth', bwidth);

FreqHist(N, edges, L1, 'k');
reporter = TO;
title([num2str(title1) ' & ' num2str(title2)]);
box off;
%ax.TickDir = 'out'
ax = gca; 
ax.TickDir = 'out';
ax.FontName = 'Calibri'; 'FixedWidth';
ax.FontSize = 18;
xline(0, 'b');

%tiledlayout(flow);
else
    fprintf('%d has no spikes in window\n',title2)
    
end


end
