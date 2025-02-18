function [avgResp,avgRespErr,contrasts,sizes,tc_stimSort] = getEaTrialResp_visstimret(data,mw,frameRateHz)
    
    nBaseFr = round(frameRateHz); % 1 second
    nStimFr = round(frameRateHz); % 1 second
%         
%     [nfr,nc] = size(data);
%     conID = celleqel2mat_padded(mw.tGratingContrast);
%     sizeID = celleqel2mat_padded(mw.tGratingDiameterDeg);
%     contrasts = unique(conID);
%     sizes = unique(sizeID);
%     ncon = length(contrasts);
%     nsize = length(sizes);
    off = mw.nScansOff;
    on = mw.nScansOn;
    
    nt = floor(nfr./(off+on));
    if nt > length(conID)
        conID = conID(1:nt);
        sizeID = sizeID(1:nt);
    end
    
    trialStartFrame = (off+1):(off+on):((off+on)*nt);
    
    trialFrameInd = arrayfun(@(x) (x-nBaseFr):(x+nStimFr-1),...
        trialStartFrame,'unif',0);
    
    data_eaTrial = cellfun(@(x) data(x,:),trialFrameInd,'unif',0);
    dff_eaTrial = cellfun(@(x) ...
        (x-(mean(x(1:nBaseFr,:),1)))./mean(x(1:nBaseFr,:),1),...
        data_eaTrial,'unif',0);
%     
%     tc_stimSort = cell(nsize,ncon);
%     for icon = 1:ncon
%         for isize = 1:nsize
%             ind = conID == contrasts(icon) & sizeID == sizes(isize);
%             tc_stimSort{isize,icon} = reshape(cell2mat(...
%                 dff_eaTrial(ind)),[nBaseFr+nStimFr,nc,sum(ind)]);
%         end
%     end
    
%     tc = reshape(cell2mat(...
%         cellfun(@(x) mean(x,3),tc_stimSort,'unif',0)),...
%         [nBaseFr+nStimFr,nc,ncon]);
%     tcErr = reshape(cell2mat(...
%         cellfun(@(x) ste(x,3),tc_stimSort,'unif',0)),...
%         [nBaseFr+nStimFr,nc,ncon]);
    avgResp = cell2mat(cellfun(...
        @(x) mean(mean(x((nBaseFr+1):(nBaseFr+nStimFr),:,:),1),3),...
        tc_stimSort,'unif',0)');
    avgRespErr = cell2mat(cellfun(...
        @(x) ste(mean(x((nBaseFr+1):(nBaseFr+nStimFr),:,:),1),3),...
        tc_stimSort,'unif',0)');
end