function [baseline, drug] = eventsPerMinDRUG(dataIN,nSweeps,sweepLength,drugStart,drugBuffer,drugDuration)
    inputFile = [dataIN,'.csv'];
    fullData =readtable(inputFile);
    data=fullData(:,1);
    data=table2array(data);
    counts=zeros(nSweeps,1);
    for iSweep = 1:nSweeps
        if length(find(data==iSweep))>0
            counts(iSweep)=length(find(data==iSweep));
        end
    end

    counts = (counts./sweepLength)*60;

    baseline = mean(counts(1:drugStart));

    drug = mean(counts(drugBuffer:(drugBuffer+drugDuration)));


end

