% go to folder of interest!!!  & set these variables
cluster1 = 240;                                 % designate first cluster of interest
cluster2 = 239;                                  %designate second cluster of interest
range = [-.050, .050];                       % designate x-axis range in sec
trange = .50;
bwidth = [.001];                               %designate bin size in sec
%%%%

ts = readNPY('spike_times.npy');           %convert all spike times from python to matlab
dts = double(ts)/30000.0000;              % create vector of all spike timestamps
cl = double(readNPY('spike_clusters.npy'));           % create vector designating cluster for each ts
if size(dts) ~= size(cl);                       % check that vector sizes are equal
       printf('warning: spike_times and spike_clusters are unequal')
end
log1 = cl(:,1) == cluster1;                     % create logical vector for first cluster of interest
log2 = cl(:,1) == cluster2;                  %create logical vector for second cluster of interest
k1 = find(log1);                                % create vector of indices for clusters of interest
k2 = find(log2);
cl_ts = [cl dts];                              % create matrices of cluster, timestamps
clone_ts = cl_ts(k1,:);                      % create matrices of timestamps for clustesr of interest
cltwo_ts = cl_ts(k2, :);

clust1 = dts(k1,:);                      % create a vector of timestamps for cluster of interest
clust2 = dts(k2,:);
%selectCL_ts = cast(selectCL_ts,'uint32');      % recast ts as a double instead of uint64
L1 = (length(clust1));                           % L = number of spikes
L2 = (length(clust2));
%deltaT = tall(zeros(L*L,1));
%column = zeros(100)


k = 1;                                              %create counter for output vector index
for i=1:L2                                           % for every element in the first spiketime vector
    for j = 1:L1                                     % for every element in the second (or mirror) spiketime vector
       test = clust2(i,:)- clust1(j,:);             % get difference between spiketimes
       if test == 0;
           test= nan;                               % eliminate zeros
       end
       if ((test <= 0.01) & (test > 0));    % Check if there is a violation
           violations (k,1) = clust1(j,:);                  % If yes, add clust 1 timestamp to vector that will point out violations
           k = k+1;                                 % update index
       end
           
    end  
end
%figure(1)
%histogram (useDeltas, 'BinLimits', range, 'Binwidth', bwidth)
%title('211 & 395')

