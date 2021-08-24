% go to folder of interest!!!  & set these variables
clear
cluster1 = 1253;                                 % designate cluster of interest 

range = [-.200, .800];                       % designate x-axis range in sec

bwidth = [.010];                               %designate bin size in sec
window = [400, 460];                                     % designate 60 sec time window for stats
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
%logical60 = clust1(:,1) == ((logical(1,1) <= clust1(:,1)) & ( clust1(:,1)<= logical(2,1);
%k2 = find(logical60);
%selectTS = clust1(k2,:);

%selectCL_ts = cast(selectCL_ts,'uint32');      % recast ts as a double instead of uint64
L1 = (length(clust1));                           % L = number of spikes




k = 1;                                              %create counter for output vector index
for i=1:(L1-1)                                           % for every element in the first spiketime vector
        if((window(1)<= clust1(i)) & (clust1(i)<= window(2)));
         
            test = clust1(i+1,:)- clust1(i,:);             % get difference between spiketimes
            useDeltas (k,1) = test;                  % If yes, add difference to vector that will create histogram
            %reporter(k,:) = k; 
             k = k+1;    
        end
                                   % update index
     
           
   
end
figure
histogram (useDeltas, 'BinLimits', range, 'Binwidth', bwidth)
title([' ISI ' num2str(cluster1)])

%average Frequency
avgHz = (L1/(clust1(L1)-clust1(1)))
windowAvgHz = (length(useDeltas)/60)

%median ISI
medISI = median(useDeltas)

%CVlog ISI
msecISI = 100*useDeltas;
logISI = log(msecISI);
CVlog = std(logISI)/mean(logISI)

%CV2
L2 = length(useDeltas);
m = 1;
for j=1:(L2-1)
    CV2(m) = 2*(abs(useDeltas(j+1)-useDeltas(j)))/(useDeltas(j+1)+ useDeltas(j));
    m = m+1;
end
meanCV2 = mean(CV2)

% Fifth percentile of the ISI distribution
fvpctISI = prctile(useDeltas,5)

% print values on chart
gtext(['avg Hz = ' num2str(avgHz)])
gtext(['median ISI = ' num2str(medISI)])
gtext(['CVlog = ' num2str(CVlog)])
gtext(['CV2 = ' num2str(meanCV2)])
gtext(['fifthISI = ' num2str(fvpctISI)])


    



%tiledlayout(flow);


