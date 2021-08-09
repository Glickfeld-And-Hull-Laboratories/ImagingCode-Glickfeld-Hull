function [UnitStruct] = CreateGoodStruct(cluster, SpikeTimes, clustersAllST) % pass only good clusters to get struct of good units & timestamps
clear k
clear j
clear log1
clear k1

cl_ts = [clustersAllST SpikeTimes];                              % create matrices of cluster, timestamps

% find vectors of spiketimes of each unit and make structure together with
% the name/number of the unit
for j = 1:length(cluster)
    unit = cluster(j);
    log1 = clustersAllST(:,1) == unit;                     % create logical vector for cluster of interest THIS IS WHERE YOU ARE!!!!!!!!!!! somehow causing problems when it switches to double digits. maybe convert srting back to double somehow???
    k1 = find(log1);  % create vector of indices for clusters of interest
    unit_timestamps = cl_ts(k1,2);                      % create a vector of timestamps for cluster of interest
    %strUnit = string(unit);
    %eval('strUnit = unit_timestamps');
    UnitStruct(j).unitID = unit;
    UnitStruct(j).timestamps = unit_timestamps;
    clear k1
    clear log1
end



end