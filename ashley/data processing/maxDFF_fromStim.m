function maxDFF_stim = maxDFF_fromStim(data,stimOnFr,stimID,nBaseFr,nStimFr)

    baseFr = cellfun(@(x) (x-nBaseFr+1):x,num2cell(stimOnFr),'unif',0);
    stimFr = cellfun(@(x) (x):(nStimFr+x-1),num2cell(stimOnFr),'unif',0);
    
    F = cellfun(@(x) mean(data(:,:,x),3),baseFr,'unif',0);
    dFF = cellfun(@(x,y) (mean(data(:,:,x),3)-y)./y,stimFr,F,'unif',0);
    clear F data
    
    stims = unique(stimID);
    maxDFF_stim = nan(size(dFF{1},1),size(dFF{1},2),length(stims));
    for istim = 1:length(stims)
        ind = stimID == stims(istim);
        maxDFF_stim(:,:,istim) = max(reshape(cell2mat(dFF(ind)),...
            [size(dFF{1},1),size(dFF{1},2),sum(ind)]),[],3);
    end    
end