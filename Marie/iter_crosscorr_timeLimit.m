% go to folder of interest!!!  & set these variables
cluster1 = 42;                                 % designate first cluster of interest
cluster2 = 42;                                  %designate second cluster of interest
                %set histogram variables
range = [-.100, .100];                       % designate x-axis range in sec
trange = .500;
bwidth = [.001];                               %designate bin size in sec

                %set time Limit variables
Lbound = 0
Ubound = inf
                %set title etc. at end
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
%some stuff to limit time window
log3 = (Lbound < clust1(:,1) < Ubound);      % create logical vector for  when timestamp is between bounds
log4 = (Lbound < clust2(:,1) < Ubound);
k1 = find(log3);                             % use logical vector to find indexes when timestamp is between bounds
k2 = find(log4);
clust1Lim = clust1(k1, :);                   % create vector of timestamps within bounds
clust2Lim = clust2(k2, :);

%selectCL_ts = cast(selectCL_ts,'uint32');      % recast ts as a double instead of uint64
L1 = (length(clust1Lim));                           % L = number of spikes
L2 = (length(clust2Lim));
%deltaT = tall(zeros(L*L,1));
%column = zeros(100)


k = 1;                                              %create counter for output vector index
i = 1;
j =1;
for i=1:L2                                           % for every element in the first spiketime vector
    for j = 1:L1                                     % for every element in the second (or mirror) spiketime vector
       test = clust2Lim(i,:)- clust1Lim(j,:);             % get difference between spiketimes
       if test == 0;
           test= nan;                               % eliminate zeros
       end
       if ((test <= trange) & (test >= -trange));    % Check if difference is in histogram range
           useDeltas (k,1) = test;                  % If yes, add difference to vector that will create histogram
           k = k+1;                                 % update index
       end
           
    end  
end
figure(1)
histogram (useDeltas, 'BinLimits', range, 'Binwidth', bwidth)
title('42')
