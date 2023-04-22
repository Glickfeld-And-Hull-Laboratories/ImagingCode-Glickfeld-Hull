function [FR, spikeTimes] = SpikeRate(dataTrain,sampleFreq)
%Takes a train of data (eg the output from cleanDataCeline) and the
%sampling frequency and returns the spike rats as well as the indices of
%the spikes

%% get the second derivative of the data train
data_diff = diff(diff(dataTrain));

spikeTimes=[];
for i = 2:length(data_diff)
    if data_diff(i) > 3 
        if data_diff(i-1)<3
            spikeTimes=[spikeTimes,i];
        end
    end
end

SpikeNum = length(spikeTimes);
timeSec = (length(dataTrain)-1)/sampleFreq;
FR=SpikeNum/timeSec;
spikeTimes=spikeTimes;

end