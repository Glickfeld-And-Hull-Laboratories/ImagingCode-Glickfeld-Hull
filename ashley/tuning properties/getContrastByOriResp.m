function [avgResp,avgRespErr,contrasts,orientations,tc_stimSort] = getContrastByOriResp(data,mw,frameRateHz)
    
    nBaseFr = round(frameRateHz); % 1 second
    nStimFr = round(frameRateHz); % 1 second
        
    [nfr,nc] = size(data);
    conID = celleqel2mat_padded(mw.tGratingContrast);
    oriID = celleqel2mat_padded(mw.tGratingDiameterDeg);
    contrasts = unique(conID);
    orientations = unique(oriID);
    ncon = length(contrasts);
    nori = length(orientations);
    off = mw.nScansOff;
    on = mw.nScansOn;
    
    nt = floor(nfr./(off+on));
    if nt > length(conID)
        conID = conID(1:nt);
        oriID = oriID(1:nt);
    end
    
    trialStartFrame = (off+1):(off+on):((off+on)*nt);
    
    trialFrameInd = arrayfun(@(x) (x-nBaseFr):(x+nStimFr-1),...
        trialStartFrame,'unif',0);
    
    data_eaTrial = cellfun(@(x) data(x,:),trialFrameInd,'unif',0);
    dff_eaTrial = cellfun(@(x) ...
        (x-(mean(x(1:nBaseFr,:),1)))./mean(x(1:nBaseFr,:),1),...
        data_eaTrial,'unif',0);
    
    tc_stimSort = cell(nori,ncon);
    for icon = 1:ncon
        for iori = 1:nori
            ind = conID == contrasts(icon) & oriID == orientations(iori);
            tc_stimSort{iori,icon} = reshape(cell2mat(...
                dff_eaTrial(ind)),[nBaseFr+nStimFr,nc,sum(ind)]);
        end
    end
    avgResp = cell2mat(cellfun(...
        @(x) mean(mean(x((nBaseFr+1):(nBaseFr+nStimFr),:,:),1),3),...
        tc_stimSort,'unif',0)');
    avgRespErr = cell2mat(cellfun(...
        @(x) ste(mean(x((nBaseFr+1):(nBaseFr+nStimFr),:,:),1),3),...
        tc_stimSort,'unif',0)');
end