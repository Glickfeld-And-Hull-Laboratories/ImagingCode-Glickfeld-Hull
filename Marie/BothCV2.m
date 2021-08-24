% go to folder of interest!!!  & set these variables
clear
cluster1 = 82;                                 % designate cluster of interest 

range = [-.200, .800];                       % designate x-axis range in sec

bwidth = [.010];                               %designate bin size in sec
window = [300, 360];                                     % designate 60 sec time window for stats
%%%%

ts = readNPY('spike_times.npy');           %convert all spike times from python to matlab
dts = double(ts)/30000.0000;              % create vector of all spike timestamps
cl = double(readNPY('spike_clusters.npy'));           % create vector designating cluster for each ts
if size(dts) ~= size(cl)                         % check that vector sizes are equal
       printf('warning: spike_times and spike_clusters are unequal')
end
logical = cl(:,1) == cluster1;                     % create logical vector for  cluster of interest
k = find(logical);                                % create vector of indices for clusters of interest
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


%CV2
L2 = length(useDeltas);
m = 1;
for j=1:(L2-1)
    CV2(m) = 2*(abs(useDeltas(j+1)-useDeltas(j)))/(useDeltas(j+1)+ useDeltas(j));
    m = m+1;
end
allCV2 = mean(CV2)

clear k i test useDeltas j L2 m CV2 

k = 1;
for i=1:(L1-1)                                           % for every element in the first spiketime vector
        if((window(1)<= clust1(i)) & (clust1(i)<= window(2)));
         
            test = clust1(i+1,:)- clust1(i,:);             % get difference between spiketimes
            useDeltas (k,1) = test;                  % If yes, add difference to vector that will create histogram
            %reporter(k,:) = k; 
             k = k+1;    
        end
                                   % update index
     
           
   
end
%CV2
L2 = length(useDeltas);
m = 1;
for j=1:(L2-1)
    CV2(m) = 2*(abs(useDeltas(j+1)-useDeltas(j)))/(useDeltas(j+1)+ useDeltas(j));
    m = m+1;
end
windowCV2 = mean(CV2)




