function [avgResp,avgRespErr,tc,tcErr,contrasts,tc_stimSort] = getContrastResp(data,mw,frameRateHz)
    
    nBaseFr = round(frameRateHz); % 1 second
    nStimFr = round(frameRateHz); % 1 second
        
    [nfr,nc] = size(data);
    stimID = celleqel2mat_padded(mw.tGratingContrast);
    contrasts = unique(stimID);
    nstim = length(contrasts);
    off = mw.nScansOff;
    on = mw.nScansOn;
    
    nt = floor(nfr./(off+on));
    if nt > length(stimID)
        stimID = stimID(1:nt);
    end
    
    trialStartFrame = (off+1):(off+on):((off+on)*nt);
    
    trialFrameInd = arrayfun(@(x) (x-nBaseFr):(x+nStimFr-1),...
        trialStartFrame,'unif',0);
    
    data_eaTrial = cellfun(@(x) data(x,:),trialFrameInd,'unif',0);
    dff_eaTrial = cellfun(@(x) ...
        (x-(mean(x(1:nBaseFr,:),1)))./mean(x(1:nBaseFr,:),1),...
        data_eaTrial,'unif',0);
    
    tc_stimSort = cell(1,nstim);
    for istim = 1:nstim
        ind = stimID == contrasts(istim);
        tc_stimSort{istim} = reshape(cell2mat(...
            dff_eaTrial(ind)),[nBaseFr+nStimFr,nc,sum(ind)]);
    end
    
    tc = reshape(cell2mat(...
        cellfun(@(x) mean(x,3),tc_stimSort,'unif',0)),...
        [nBaseFr+nStimFr,nc,nstim]);
    tcErr = reshape(cell2mat(...
        cellfun(@(x) ste(x,3),tc_stimSort,'unif',0)),...
        [nBaseFr+nStimFr,nc,nstim]);
    avgResp = cell2mat(cellfun(...
        @(x) mean(mean(x((nBaseFr+1):(nBaseFr+nStimFr),:,:),1),3),...
        tc_stimSort,'unif',0)');
    avgRespErr = cell2mat(cellfun(...
        @(x) ste(mean(x((nBaseFr+1):(nBaseFr+nStimFr),:,:),1),3),...
        tc_stimSort,'unif',0)');
end