function [isResponsive_any,isResponsive_eaStim] = findRespCells_anyContrast(data,frameRateHz)
    
    nBaseFr = round(frameRateHz); % 1 second
    nStimFr = round(frameRateHz); % 1 second
    
    isResponsive_eaStim = cell2mat(cellfun(@(x) ttest(squeeze(mean(x(1:nBaseFr,:,:),1)),...
        squeeze(mean(x((nBaseFr+1):(nStimFr+nBaseFr),:,:),1)),...
        'dim',2,'tail','left'),data,'unif',0));
    
    allData = [];
    for istim = 1:length(data)
        allData = cat(3,allData,data{istim});
    end
    isResponsive_any = sum(isResponsive_eaStim,2)>0 | ...
        ttest(squeeze(mean(allData(1:nBaseFr,:,:),1)),...
        squeeze(mean(allData((nBaseFr+1):(nStimFr+nBaseFr),:,:),1)),...
        'dim',2,'tail','left');    
end