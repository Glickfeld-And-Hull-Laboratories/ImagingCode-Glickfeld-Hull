% Imports KS (and other) data to begin matlab analysis.
%
% Gets .imec file that was used for KS and phy to parse metadata.
%
% Returns vectors of unit name/numbers (as scalars) that are classified as
% good, noise, and mua in py (GoodUnits, NoiseUnits, MUnits).
%
% Returns structure of all units that contains unit (as a string( in field
% unitID and vectors of timestamps for corresponding units (as vectors) in
% field .timstamps (AllUnitStruct).
%
% Also returns structure of good units that contains unit (as a string) in
% field .unitID and vectors of timestamps for corresponding units (as
% vectors) in field .timestamps (GoodUnitStruct).
%
%

function [cluster_struct, NoiseUnits, GoodUnits, MUnits, AllUnitStruct, GoodUnitStruct, MultiUnitStruct] = ImportKSdata (); % Add LaserStim to left side if wanted

    SpikeTimes = SampToSec();         % use SampToSec funtion to get spiketimes and convert them to seconds using exact sampling rate
    
    clustersAllST = double(readNPY('spike_clusters.npy'));           % For every spiketime, the cluster associated with it
    
    fileID = fopen('cluster_group.tsv');                             % Get sorting data: which units are good, noise, and mu
    C = textscan(fileID, '%s %s', 'HeaderLines', 1);
    cluster = C{1};
    group = C{2};
    cluster_struct = struct('cluster', cluster, 'group', group);
    
    [NoiseUnits, GoodUnits, MUnits] = ReadInClusters(cluster_struct);
    
    [AllUnitStruct] = CreateUnivStruct(cluster, SpikeTimes, clustersAllST);
    
    cellGoodUnits = {GoodUnits};
    
    [GoodUnitStruct] = CreateGoodStruct(GoodUnits, SpikeTimes, clustersAllST);
    [MultiUnitStruct] = CreateGoodStruct(MUnits, SpikeTimes, clustersAllST);
    
   % fid = fopen('LaserStim.txt');
   % LaserStim = fscanf(fid, '%f');
  %  fclose(fid);
    
    
    
end
    