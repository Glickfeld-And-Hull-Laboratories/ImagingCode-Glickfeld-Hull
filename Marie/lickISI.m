allLicks = [];
for i = 1:109
allLicks = [allLicks JuiceLicksStateMach{i,2}];
end
clust1 = allLicks.';
L1 = (length(clust1));                           % L = number of spikes
clear('useDeltas');



k = 1;                                              %create counter for output vector index
for i=1:(L1-1)                                           % for every element in the first spiketime vector
       test = clust1(i+1,:)- clust1(i,:);             % get difference between spiketimes
       useDeltas (k,1) = test;                  % If yes, add difference to vector that will create histogram
           k = k+1;                                 % update index
     
           
    
end
figure
histogram (useDeltas, 'BinLimits', [-0.025, 0.3], 'Binwidth', 0.001, 'Facecolor', 'k', 'Linestyle', 'none', 'Facealpha', 1) %'Normalization', 'probability', 

box off
ax = gca; 
ax.TickDir = 'out';
ax.FontName = 'Calibri', 'FixedWidth';
ax.FontSize = 18;