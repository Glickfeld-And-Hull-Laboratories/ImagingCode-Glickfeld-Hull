function [dff_eaTrial,avgResp_eaTrial] = getEaTrialResp_visstimret(data,mw,frameRateHz)
    
    nBaseFr = round(frameRateHz); % 1 second
    nStimFr = round(frameRateHz.*0.5); % 0.5 second
%         
    nfr = cellfun(@(x) size(x,1),data);
%     conID = celleqel2mat_padded(mw.tGratingContrast);
%     sizeID = celleqel2mat_padded(mw.tGratingDiameterDeg);
%     contrasts = unique(conID);
%     sizes = unique(sizeID);
%     ncon = length(contrasts);
%     nsize = length(sizes);
    off = mw{1}.nScansOff;
    on = mw{1}.nScansOn;
    
    nt = cellfun(@(x) length(x.tGratingContrast),mw);
    nt_allframes = floor(nfr./double(off+on));
    if any(nt > nt_allframes)
        nt(nt>nt_allframes) = nt_allframes(nt>nt_allframes);
    end
    
    trialStartFrame = arrayfun(@(x) (off+1):(off+on):((off+on)*x),nt,'unif',0);
    
    data_eaTrial = [];
    for irun = 1:length(data)
        trialFrameInd = arrayfun(@(x) (x-nBaseFr):(x+nStimFr-1),...
            trialStartFrame{irun},'unif',0);
        d = cellfun(@(x) data{irun}(x,:),trialFrameInd,'unif',0);
        data_eaTrial = cat(2,data_eaTrial,d);
    end
    
    dff_eaTrial = cellfun(@(x) ...
        (x-(mean(x(1:nBaseFr,:),1)))./mean(x(1:nBaseFr,:),1),...
        data_eaTrial,'unif',0);

    avgResp_eaTrial = cell2mat(cellfun(...
        @(x) mean(mean(x((nBaseFr+1):(nBaseFr+nStimFr),:,:),1),3),...
        dff_eaTrial,'unif',0)');
end