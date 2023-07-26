function [counts] = spikesPerSweep(dataIN,nSweeps)
    inputFile = [dataIN,'.csv']
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
    ylabel("# events")
    %ylim([0 50])
    print(fullfile([dataIN 'events_per_sweep.pdf']),'-dpdf','-bestfit')

end

