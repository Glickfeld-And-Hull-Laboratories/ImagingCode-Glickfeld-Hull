function [counts] = spikesPerSweep(inputFile,nSweeps)
    %UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here
    %read in file as "data"
    fullData =readtable(inputFile);
    data=fullData(:,1);
    data=table2array(data);
%     counts = groupcounts(data);
    counts=zeros(nSweeps,1);
    for iSweep = 1:nSweeps
        if length(find(data==iSweep))>0
            counts(iSweep)=length(find(data==iSweep));
        end
    end
    figure;
    plot(counts);
    set(gca, 'TickDir', 'out')
    box off
    xlabel("sweep")
    ylabel("# spikes")
end

