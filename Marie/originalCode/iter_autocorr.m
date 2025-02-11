% go to folder of interest!!!  & set these variables
cluster = 87;                                 % designate the cluster of interest
range = [-.200, .200];                       % designate x-axis range in sec
trange = .200;
bwidth = [.001];                               %designate bin size in sec
%%%%

ts = readNPY('spike_times.npy');           %convert all spike times from python to matlab
dts = double(ts)/30000.0000;              % create vector of all spike timestamps
cl = double(readNPY('spike_clusters.npy'));           % create vector designating cluster for each ts
if size(dts) ~= size(cl);                       % check that vector sizes are equal
       printf('warning: spike_times and spike_clusters are unequal')
end
log = cl(:,1) == cluster;                     % create logical vector for cluster of interest
k = find(log);                                % create vector of indices for cluster of interest
cl_ts = [cl dts];                              % create matrix of cluster, timestamp
onecl_ts = cl_ts(k,:);                      % create matrix of timestamps for cluster of interest
clusts = dts(k,:);                      % create a vector of timestamps for cluster of interest
%selectCL_ts = cast(selectCL_ts,'uint32');      % recast ts as a double instead of uint64
L = (length(clusts));                           % L = number of spikes
%deltaT = tall(zeros(L*L,1));
%column = zeros(100)


k = 1;                                              %create counter for output vector index
for i=1:L                                           % for every element in the first spiketime vector
    for j = 1:L                                     % for every element in the second (or mirror) spiketime vector
       test = clusts(i,:)- clusts(j,:);             % get difference between spiketimes
       if test == 0;
           test= nan;                               % eliminate zeros
       end
       if ((test <= trange) & (test >= -trange));    % Check if difference is in histogram range
           useDeltas (k,1) = test;                  % If yes, add difference to vector that will create histogram
           k = k+1;                                 % update index
       end
           
    end  
end
histogram (useDeltas, 'BinLimits', range, 'Binwidth', bwidth)

