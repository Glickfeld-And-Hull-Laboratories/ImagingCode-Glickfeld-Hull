% go to folder of interest!!!
cluster = 13;                                 % designate the cluster of interest
range = [-.100, .100];                       % designate x-axis range in sec
bwidth = [.001];                               %designate bin size in sec

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
L = (length(clusts));
%deltaT = tall(zeros(L*L,1));
%column = zeros(100)


k = 0
for i=1:L
    for j = 1:L
       deltaT(i,j)=clusts(i,:)- clusts(j,:);
       if deltaT(i,j) == 0;
           deltaT(i, j)= nan;
       end
       if ((deltaT(i,j) <= range) & (delta(i,j) >= -range));
           useDeltas (k,1) = deltaT(i,j);
           k = k+1;
       end
           
    end  
end
% histogram (deltaT, 'BinLimits', range, 'Binwidth', bwidth)

