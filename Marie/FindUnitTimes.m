function ClusterTimes = FindUnitTimes(unit, SpikeTimes, clusterID)

log1 = clusterID(:,1) == unit;                     % create logical vector for first cluster of interest
k1 = find(log1);                                % create vector of indices for clusters of interest

ClusterTimes = SpikeTimes(k1);                      % create matrices of timestamps for clustesr of interest
end