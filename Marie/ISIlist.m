% go to folder of interest!!!  & set these variables
% creat unitIDmatrix= list of clusterIDs of interest
range = [-.100, .900];                       % designate x-axis range in sec

bwidth = [.001];                               %designate bin size in sec
%%%%

ncell = length(unitIDmatrix);
m = 1;
for m = 1:ncell
    clearvars -except range bwidth trange unitIDmatrix m ncell
    cluster1 = unitIDmatrix(m, 1);                                    % designate cluster of interest 



ts = readNPY('spike_times.npy');           %convert all spike times from python to matlab
dts = double(ts)/30000.0000;              % create vector of all spike timestamps
cl = double(readNPY('spike_clusters.npy'));           % create vector designating cluster for each ts
if size(dts) ~= size(cl);                       % check that vector sizes are equal
       printf('warning: spike_times and spike_clusters are unequal')
end
log = cl(:,1) == cluster1;                     % create logical vector for  cluster of interest
k = find(log);                                % create vector of indices for clusters of interest
cl_ts = [cl dts];                              % create matrices of cluster, timestamps
clone_ts = cl_ts(k,:);                      % create matrices of timestamps for clustesr of interest


clust1 = dts(k,:);                      % create a vector of timestamps for cluster of interest

%selectCL_ts = cast(selectCL_ts,'uint32');      % recast ts as a double instead of uint64
L1 = (length(clust1));                           % L = number of spikes




k = 1;                                              %create counter for output vector index
for i=1:(L1-1)                                           % for every element in the first spiketime vector
       test = clust1(i+1,:)- clust1(i,:);             % get difference between spiketimes
       useDeltas (k,1) = test;                  % If yes, add difference to vector that will create histogram
           k = k+1;                                 % update index
     
           
    
end
figure(m)
histogram (useDeltas, 'BinLimits', range, 'Binwidth', bwidth, 'Normalization','probability')
title([' ISI ' num2str(cluster1)])
m= m+1
end
%tiledlayout(flow);


