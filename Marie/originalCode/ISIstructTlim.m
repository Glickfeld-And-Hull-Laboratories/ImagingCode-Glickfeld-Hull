% go to folder of interest!!!  & set these variables
function ISIstructTlim(struct, unit, tmin, tlim, range, bwidth)

n = find([struct.unitID] == unit); %% n changes to index in struct pointing to specified unit
unitTS = [struct(n).timestamps];  %% Make vector, TimeStamps2, that has timestamps from unit.                    % create a vector of timestamps for cluster of interest
unitTS = unitTS((tmin < unitTS) & (unitTS < tlim)); %%%% adjust here to change tlim structure
%selectCL_ts = cast(selectCL_ts,'uint32');      % recast ts as a double instead of uint64
L1 = (length(unitTS));                           % L = number of spikes




k = 1;                                              %create counter for output vector index
for i=1:(L1-1)                                           % for every element in the first spiketime vector
       test = unitTS(i+1,:)- unitTS(i,:);             % get difference between spiketimes
       useDeltas (k,1) = test;                  % If yes, add difference to vector that will create histogram
           k = k+1;                                 % update index
     
           
    
end
figure
histogram (useDeltas, 'BinLimits', range, 'Binwidth', bwidth, 'Facecolor', 'k', 'Linestyle', 'none', 'Facealpha', 1) %'Normalization', 'probability', 
title([' ISI ' num2str(unit)])
box off
ax = gca; 
ax.TickDir = 'out';
ax.FontName = 'Calibri', 'FixedWidth';
ax.FontSize = 18;
%tiledlayout(flow);
end


