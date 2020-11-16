clear all
close all
ds = 'FSAV_attentionV1'; % 'FSAV_V1_100ms_naive'  'FSAV_V1_naive_GCaMP6m'  'FSAV_attentionV1'   'FSAV_attentionV1_noAttn'
cellsOrDendrites = 1;
doLoadPreviousAnalysis = false;
analysisDate = '200210';
attnAnalysisDate = '191211';
doDecoding = true;
%%
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms

eval(ds)
titleStr = ds(6:end);
mice = unique({expt.SubNum});
mouse_str = ['i' strjoin(mice,'_i')];

if cellsOrDendrites == 1
    load(fullfile(rc.caOutputDir,ds,...
        [mouse_str '_trOutcomeStruct_cells' ds(5:end) '.mat']));
    fnout = fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
        [titleStr '_' datestr(now,'yymmdd') '_']); 
elseif cellsOrDendrites == 2
    load(fullfile(rc.caOutputDir,ds,...
        [mouse_str '_trOutcomeStruct_cells' ds(5:end) '.mat']));
    fnout = fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
        [titleStr '_den_' datestr(now,'yymmdd') '_']); 
end

if strcmp(ds,'FSAV_attentionV1')
    load(fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
        'V1_100ms_naive_anticipation_imgStats'))
    imgStats_naive = imgStats;
    clear imgStats
end

%%
nav = 2;
nBaselineFr = mouse(1).expt(1).info.preAlignFrames;
nFrames1s = frameRateHz;
nexp = size(expt,2);
nCycles = 8;
% lateCycles = 5:nCycles;
lateWinFr = (45:88)+nBaselineFr;
firstWinFr = (3:44)+nBaselineFr;
minTargetRT = (nVisDelayFr_target+respwin_target(1)-nBaselineFr)./frameRateHz*1000;

oriBinSize = 45;
orientations = 0:oriBinSize:(180-oriBinSize);
oriBinEdges = [0, (oriBinSize/2):oriBinSize:(180-(oriBinSize/2)), 180];
nOri = length(orientations);

nMovWin = 15;
movWinLabelFr = 30:(30+nMovWin-1);
movWinLabelFr_target = 30:(30+nMovWin-1);
% movWinLabelMs = 

% minCellN_SIFRmatch = 36;

trOutType = {'h';'m';'fa';'cr';'yes';'no'};
% trOutTypeName = {'H-All';'H-HT';'H-ET';'M-All';'M-HT';'M-ET';'FA';'CR'};
% trOutTypeName = {'H';'M';'FA';'CR'};
%% pool experiment data
if doLoadPreviousAnalysis
    load([fnout(1:end-7) analysisDate '_imgAnalysisData'])
else    
    rng(0)
    antiDataExpt = struct;
    oriTuningExpt = struct;
    targetDataExpt = struct;
    decodeDataExpt = struct;
    decodeDataExpt.av(visualTrials).name = 'Visual';
    decodeDataExpt.av(auditoryTrials).name = 'Auditory';
    decodeDataExpt.catch.name = 'Catch';
    nTargets = 2; %sum(unique(targetInd) > 1);
    nTrialsPerExpt = nan(1,nexp);
    rewExptInd = nan(1,nexp);
    for imouse = 1:size(mouse,2)
        for iexp = 1:size(mouse(imouse).expt,2)
            if imouse == 1 & iexp == 1
                exptN = 1;
            else
                exptN = exptN+1;
            end

            d = mouse(imouse).expt(iexp);
            exptID = strcmp({expt.SubNum},mouse(imouse).mouse_name) & strcmp({expt.date},d.date);
            if strcmp(ds,'FSAV_attentionV1')
                if expt(exptID).catchRew == 1
                    rewExptInd(exptN) = true;
                else
                    rewExptInd(exptN) = false;
                end
            else
                rewExptInd(exptN) = false;
            end
            nTrialsPerExpt(exptN) = length(d.av(visualTrials).align(1).outcome) + ...
                length(d.av(auditoryTrials).align(1).outcome);

            cycLengthFr = d.info.cycTimeFrames;
            nCycLongTC = ceil(longTrialLengthFr./cycLengthFr);

            maxCycles = min([max(d.av(visualTrials).align(alignStart).nCycles),...
                max(d.av(auditoryTrials).align(alignStart).nCycles)]);
            
            cycRespEaTrial = cell(1,2);
            outcomeEaTrial = cell(1,2);
            stimEaTrial = cell(1,2);
            for iav = 1:2
                de = d.av(iav).align(alignStart);
                outs = de.outcome;
                outs(strcmp(outs,'failure')) = {'fa'};
                outs(strcmp(outs,'success')) = {'h'};
                outs(strcmp(outs,'ignore')) = {'m'};
                outcomeEaTrial{iav} = outs;
                n = length(de.nCycles);
                trN = 0;
                cycRespEaTrial{iav} = cell(1,n);
                stimEaTrial{iav} = nan(1,n);
                for itr = 1:n   
                    r = nan(de.nCycles(itr),size(de.respTC,2));
                    for icyc = 1:de.nCycles(itr)
                        cycStartOffset = ((icyc-1).*cycLengthFr)+nBaselineFr;
                        tc = de.respTC(...
                            (cycStartOffset-nBaselineFr+1):(cycStartOffset+nFrames1s),...
                            :,itr);
                        r(icyc,:) = mean(tc(respwin,:),1) - mean(tc(basewin_0,:),1);
                    end                      
                    if ~strcmp(de.outcome(itr),'failure')
                        trN = trN+1;
                        cycStartOffset = ((icyc).*cycLengthFr)+nBaselineFr;
                        tc = d.av(iav).align(alignTarget).respTC(:,:,trN);
                        r(de.nCycles(itr)+1,:) = ...
                            mean(tc(respwin_target,:),1) - mean(tc(basewin_0_target,:),1);
                        if iav == 1
                            stimEaTrial{iav}(itr) = d.av(iav).align(alignTarget).ori(trN);
                        else
                            stimEaTrial{iav}(itr) = d.av(iav).align(alignTarget).amp(trN);
                        end
                    end
                    cycRespEaTrial{iav}{itr} = r; 
                end
            end
            decodeDataExpt(exptN).cycRespEaTrial = cycRespEaTrial;
            decodeDataExpt(exptN).outcomeEaTrial = outcomeEaTrial;
            decodeDataExpt(exptN).stimEaTrial = stimEaTrial;
            
            
            % catch trials
            if length(d.av(1).align) == 5
                iav = 1;
                de = d.av(iav).align(5);
                outs = de.outcome;
                outs(strcmp(outs,'failure')) = {'fa'};
                outs(strcmp(outs,'success')) = {'h'};
                outs(strcmp(outs,'ignore')) = {'m'};
                outcomeEaTrial = outs;
                catch_outs = de.catchOutcome;
                catch_outs(cellfun(@isempty,catch_outs)) = {nan};
                catch_outs(strcmp(catch_outs,'FA')) = {'h'};
                catch_outs(strcmp(catch_outs,'CR')) = {'m'};
                catchOutcomeEaTrial = catch_outs;
                catchOriEaTr = de.catchOri;
                catchOriEaTr(cellfun(@isnan,catch_outs)) = nan;
                
                n = length(de.nCycles);
                cycRespEaTrial = cell(1,n);
                stimEaTrial = nan(1,n);
                for itr = 1:n   
                    r = nan(de.nCycles(itr),size(de.respTC,2));
                    for icyc = 1:de.nCycles(itr)
                        cycStartOffset = ((icyc-1).*cycLengthFr)+nBaselineFr;
                        tc = de.respTC(...
                            (cycStartOffset-nBaselineFr+1):(cycStartOffset+nFrames1s),...
                            :,itr);
                        r(icyc,:) = mean(tc(respwin,:),1) - mean(tc(basewin_0,:),1);
                    end                      
                    if ~strcmp(outs(itr),'fa')
                        trN = trN+1;
                        cycStartOffset = ((icyc).*cycLengthFr)+nBaselineFr;
                        tc = de.respTC(:,:,itr);
                        r(de.nCycles(itr)+1,:) = ...
                            mean(tc(respwin_target,:),1) - mean(tc(basewin_0_target,:),1);
                        
                        stimEaTrial(itr) = de.amp(itr);
                        
                    end
                    cycRespEaTrial{itr} = r; 
                end
                cCyc = de.catchCycle;
                cCyc(cellfun(@isnan,catch_outs)) = nan;
                
                decodeDataExpt(exptN).catch.cycRespEaTrial = cycRespEaTrial;
                decodeDataExpt(exptN).catch.catchOutcomeEaTrial = catchOutcomeEaTrial;
                decodeDataExpt(exptN).catch.catchStimEaTrial = catchOriEaTr;
                decodeDataExpt(exptN).catch.outcomeEaTrial = outcomeEaTrial;
                decodeDataExpt(exptN).catch.stimEaTrial = stimEaTrial;
                decodeDataExpt(exptN).catch.catchCycN = cCyc;
                
            end
            

            tc_AV = [];
    %         cycTC = cell(2,nCycles);
            cycTC = cell(2,maxCycles);
            longTC = cell(2,1);
            for iav = 1:2
                dd = d.av(iav);
                trOut = [];
                trStim = [];
                trResp = [];
                trTC = [];
                trResp_movWin = [];
                nCyc_dc = [];
                for ialign = 2:4
                    ddd = dd.align(ialign);
    %                 if iav 
    %                     fprintf('%s/%s auditory trials < 200 ms\n', ...
    %                         num2str(sum(ddd.reactTime < 200)),num2str(length(ddd.reactTime)))
    %                 end
                    tc = ddd.respTC;
                    if ialign == 4
                        tc_temp = circshift(tc,-1,1);
                        tc_temp(end,:,:) = nan;
                        trTC = cat(3,trTC,tc);
                        trResp = cat(2,trResp,squeeze(mean(tc(respwin_target,:,:),1)...
                            - mean(tc(basewin_0_target,:,:),1)));
                    else
                        trResp = cat(2,trResp,squeeze(mean(tc(respwin,:,:),1)...
                            - mean(tc(basewin_0,:,:),1)));
                        trTC = cat(3,trTC,circshift(tc,-1,1)- mean(tc(basewin_0_target,:,:),1));
                    end
                    nfr = length(respwin)-1;
                    nMovWin = 15;
                    trResp_movWin_temp = nan(size(tc,2),size(tc,3),nMovWin);
                    if ialign == 4
                        tempwinstart = 28;
                    else
                        tempwinstart = 29;
                    end
                    for iwin = 1:nMovWin
                        tempwin = tempwinstart+(iwin-1):(tempwinstart+(iwin-1)+nfr);
                        if ialign == 4
                            trResp_movWin_temp(:,:,iwin) = squeeze(mean(tc(tempwin,:,:),1)) - ...
                                squeeze(mean(tc(basewin_0_target,:,:),1));
                        else
                            trResp_movWin_temp(:,:,iwin) = squeeze(mean(tc(tempwin,:,:),1)) - ...
                                squeeze(mean(tc(basewin_0,:,:),1));
                        end
                    end
                    trResp_movWin = cat(2,trResp_movWin,trResp_movWin_temp);
                    if ialign == alignCR
                        trOut = cat(2,trOut,repmat({'cr'},[1,length(ddd.outcome)]));
                    else
                        trOut = cat(2,trOut,ddd.outcome);
                    end
                    if isempty(ddd.ori) & isempty(ddd.amp)
                        trStim = cat(2,trStim,zeros(1,length(ddd.outcome)));
                    elseif iav == 1 | strcmp(ds,'FSAV_V1_audControl')
                        trStim = cat(2,trStim,ddd.ori);
                    elseif iav == 2
                        trStim = cat(2,trStim,ddd.amp);
                    end
                    nCyc_dc = cat(2,nCyc_dc,ddd.nCycles);
                end
                trOut(strcmp(trOut,'success')) = {'h'};
                trOut(strcmp(trOut,'ignore')) = {'m'};
                trOut(strcmp(trOut,'failure')) = {'fa'};

                decodeDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
                decodeDataExpt(exptN).av(iav).outcome = trOut;
                decodeDataExpt(exptN).av(iav).stim = trStim;
                decodeDataExpt(exptN).av(iav).resp = trResp;
                decodeDataExpt(exptN).av(iav).movWinResp = trResp_movWin;
                decodeDataExpt(exptN).av(iav).tc = trTC(1:(nBaselineFr*2),:,:);
                decodeDataExpt(exptN).av(iav).nCyc = nCyc_dc;    
                if iav == 1
                    binEdges = oriBins;
                elseif iav == 2
                    binEdges = ampBins;
                end
                trStimID = discretize(trStim,binEdges);
                rng(0)
                trOutTC = cell(1,length(trOutType));
                
                if strcmp(ds,'FSAV_attentionV1') | strcmp(ds,'FSAV_attentionV1_noAttn')
                    hitInd = strcmp(trOut,'h');
                    missInd = strcmp(trOut,'m');
                    matchInd = getMatchedOutcomeTrialIndex(...
                        trOut,trStimID,minTrN);
                    trOutTC = {trTC(:,:,matchInd & hitInd),trTC(:,:,matchInd & missInd),...
                        trTC(:,:,strcmp(trOut,'fa')),trTC(:,:,strcmp(trOut,'cr')),...
                        cat(3,trTC(:,:,matchInd & hitInd),trTC(:,:,strcmp(trOut,'fa'))),...
                        cat(3,trTC(:,:,matchInd & missInd),trTC(:,:,strcmp(trOut,'cr')))};
                end
                targetDataExpt(exptN).av(iav).trOutTC = trOutTC;

                de = d.av(iav).align(alignStart);
                if strcmp(ds,'FSAV_attentionV1')|strcmp(ds,'FSAV_attentionV1_noAttn')
                    hits = strcmp(de.outcome,'success');
                elseif strcmp(ds,'FSAV_V1_100ms_naive')
                    hits = true(1,length(de.outcome));
                end
                misses = strcmp(de.outcome,'ignore');
                tc_AV = cat(3,tc_AV,de.respTC(:,:,hits)); % make hits or misses
    %             for icyc = 1:nCycles
                for icyc = 1:maxCycles
                    tc = de.respTC(:,:,de.nCycles >= icyc & (hits | misses));
                    cycStartOffset = ((icyc-1).*cycLengthFr)+nBaselineFr;
                    cycTC{iav,icyc} = tc(...
                        (cycStartOffset-nBaselineFr+1):(cycStartOffset+nFrames1s),:,:);
                end
                longTC{iav} = de.respTC(1:(longTrialLengthFr+nBaselineFr),:,...
                    de.nCycles >= nCycLongTC & (hits | misses));

                de = d.av(iav).align(alignTarget);
                if iav == 1 | strcmp(ds,'FSAV_V1_audControl')
                    binEdges = oriBins;
                    targets = de.ori;
                elseif iav == 2
                    binEdges = ampBins;
                    targets = de.amp;
                end
                tc = de.respTC(1:(nBaselineFr*2),:,:);
                targetInd = discretize(targets,binEdges);
                targetDataExpt(exptN).av(iav).tc = cell(1,2);
                for itar = 1:nTargets
                    ind = targetInd == itar+1;
                    targetDataExpt(exptN).av(iav).tc{itar} = tc(:,:,ind);
                    targetDataExpt(exptN).av(iav).targets{itar} = targets(ind);
                end
            end
            antiDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
            antiDataExpt(exptN).mouse = mouse(imouse).mouse_name;
            antiDataExpt(exptN).exptCycLengthFr = cycLengthFr;
            antiDataExpt(exptN).longTC = longTC;
            antiDataExpt(exptN).cycTC = cycTC;
            targetDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
            if size(d.av(visualTrials).align,2) > 4
                r = squeeze(mean(...
                    d.av(visualTrials).align(5).respTC(respwin_target,:,:),1) - ...
                    mean(d.av(visualTrials).align(5).respTC(basewin_0_target,:,:),1));
                trOut =  d.av(visualTrials).align(5).catchOutcome;
                ind = cellfun(@(x) ~isempty(x),trOut);
                trOut = trOut(ind);
                r = r(:,ind);
                trOut(strcmp(trOut,'FA')) = {'h'};
                trOut(strcmp(trOut,'CR')) = {'m'};
                decodeDataExpt(exptN).av(visualTrials).catchResp = r;
                decodeDataExpt(exptN).av(visualTrials).catchOutcome = trOut;
                decodeDataExpt(exptN).av(visualTrials).catchStim = ...
                    d.av(visualTrials).align(5).catchOri;
            end
            eaCycSI = cellfun(@(x,y) ...
                getSelectivityIndex(squeeze(mean(x(respwin,:,:),1)-mean(x(basewin_0,:,:),1))',...
                squeeze(mean(y(respwin,:,:),1)-mean(y(basewin_0,:,:),1))'),...
                cycTC(visualTrials,:),...
                cycTC(auditoryTrials,:),'unif',0);
            antiDataExpt(exptN).eaCycSI = eaCycSI;
            
            do = d.oriTuning;
            if size(do.oriResp,2) == 8
                oriTuningExpt(exptN).oriResp = do.oriResp(:,1:2:8);
                oriTuningExpt(exptN).oriRespErr = do.oriRespSem(:,1:2:8);
            elseif size(do.oriResp,2) ~= 4
                error('error in orientations used for passivie tuning')
            else
                oriTuningExpt(exptN).oriResp = do.oriResp;
                oriTuningExpt(exptN).oriRespErr = do.oriRespSem;
            end
            oriTuningExpt(exptN).fit = do.oriFit;
            oriTuningExpt(exptN).isTuned = do.oriFitReliability < tuningReliabilityThresh;
            [~,oriTuningExpt(exptN).fitPeak] = max(do.oriFit,[],1);
            [~,oriPref] = histc(oriTuningExpt(exptN).fitPeak,oriBinEdges);
            oriPref(oriPref == length(orientations)+1 | oriPref == length(oriBinEdges) == 1) = 1;
            oriTuningExpt(exptN).oriPref = oriPref;
            oriTuningExpt(exptN).tuningReliability = do.oriFitReliability;
        end
    end

    shortCycExptInd = cell2mat({antiDataExpt.exptCycLengthFr}) == 11;
    isShortCycExpt = [];
    respCellsExpt = struct;
    for iexp = 1:nexp
        respCellsExpt(iexp).mouse = antiDataExpt(iexp).mouse;
        firstTCAV = cat(3,antiDataExpt(iexp).cycTC{visualTrials,1},...
            antiDataExpt(iexp).cycTC{auditoryTrials,1});
        fprintf('%s: %s first stim\n',antiDataExpt(iexp).exptName,num2str(size(firstTCAV,3)))

        lateTC_vis = [];
        lateTC_aud = [];
        lateCycles = lateCyclesMin:size(antiDataExpt(iexp).cycTC,2);
        for icyc = 1:length(lateCycles)
            lateTC_vis = cat(3,lateTC_vis,...
                antiDataExpt(iexp).cycTC{visualTrials,lateCycles(icyc)});
            lateTC_aud = cat(3,lateTC_aud,...
                antiDataExpt(iexp).cycTC{auditoryTrials,lateCycles(icyc)});
        end
        lateTCAV = cat(3,lateTC_vis,lateTC_aud);
        fprintf('%s: %s late stim\n',antiDataExpt(iexp).exptName,num2str(size(lateTCAV,3)))
        
        longTCAV = cat(3,antiDataExpt(iexp).longTC{visualTrials,1},...
            antiDataExpt(iexp).longTC{auditoryTrials,1});
        fprintf('%s: %s long trials\n',antiDataExpt(iexp).exptName,num2str(size(longTCAV,3)))
        
        % noise correlation
        
        rv = squeeze(mean(lateTC_vis(respwin,:,:),1)-mean(lateTC_vis(basewin_0,:,:),1))';
        ra = squeeze(mean(lateTC_aud(respwin,:,:),1)-mean(lateTC_aud(basewin_0,:,:),1))';
        rav = cat(1,rv,ra);
%         rv_rect = rv;
%         rv_rect(rv<0) = nan;
%         ra_rect = ra;
%         ra_rect(ra<0) = nan;
%         rav_rect = cat(1,rv_rect,ra_rect);
%         rv_rect = rv-min(rv)+1;
%         ra_rect = ra-min(ra)+1;
%         rav_rect = cat(1,rv,ra)-min(cat(1,rv,ra))+1;
        
        rsc_vis = corrcoef(rv);
        rsc_aud = corrcoef(ra);
        rsc_all = corrcoef(rav);
        
        nc = size(rv,2);
        gm_vis = nan(nc);
        gm_aud = nan(nc);
        gm_all = nan(nc);
        for icell1 = 1:nc
            rv1 = rv(:,icell1);
            ra1 = ra(:,icell1);
            rav1 = rav(:,icell1);
            for icell2 = 1:nc
                rv2 = rv(:,icell2);
                ra2 = ra(:,icell2);
                rav2 = rav(:,icell2);
                d = [rv1;rv2];
                gm_vis(icell1,icell2) = geomean(d-min(d)+1)+min(d)-1;
                d = [ra1;ra2];
                gm_aud(icell1,icell2) = geomean(d-min(d)+1)+min(d)-1;
                d = [rav1;rav2];
                gm_all(icell1,icell2) = geomean(d-min(d)+1)+min(d)-1;
            end
        end
        
        respCellsExpt(iexp).gm = {gm_vis,gm_aud,gm_all};
        respCellsExpt(iexp).rsc = {rsc_vis,rsc_aud,rsc_all};
        % randomly select trials for responsive cells test so that each
        % experiment test the same number of trials
        rng(0)
        if size(lateTCAV,3) > minTrN_lateresp
            ind = randsample(size(firstTCAV,3),minTrN_firstresp);
            firstRespCells = ttest(...
                squeeze(mean(firstTCAV(respwin,:,ind),1)),...
                squeeze(mean(firstTCAV(basewin_0,:,ind),1)),...
                'dim',2,'tail','right','alpha',cellGroupsAlpha);
            ind = randsample(size(lateTCAV,3),minTrN_lateresp);
            lateCycRespCells = ttest(...
                squeeze(mean(lateTCAV(respwin,:,ind),1)),...
                squeeze(mean(lateTCAV(basewin_0,:,ind),1)),...
                'dim',2,'tail','right','alpha',cellGroupsAlpha);
        else
            firstRespCells = ttest(...
                squeeze(mean(firstTCAV(respwin,:,:),1)),...
                squeeze(mean(firstTCAV(basewin_0,:,:),1)),...
                'dim',2,'tail','right','alpha',cellGroupsAlpha);
            lateCycRespCells = ttest(...
                squeeze(mean(lateTCAV(respwin,:,:),1)),...
                squeeze(mean(lateTCAV(basewin_0,:,:),1)),...
                'dim',2,'tail','right','alpha',cellGroupsAlpha);
        end        
        if size(longTCAV,3) > minTrN_latewin
            ind = randsample(size(longTCAV,3),minTrN_latewin);
            lateRespCells = ttest(...
                squeeze(mean(longTCAV(lateWinFr,:,ind),1)),...
                squeeze(mean(longTCAV(basewin,:,ind),1)),...
                'dim',2,'tail','right','alpha',cellGroupsAlpha);
            lateSuppCells = ttest(...
                squeeze(mean(longTCAV(lateWinFr,:,ind),1)),...
                squeeze(mean(longTCAV(basewin,:,ind),1)),...
                'dim',2,'tail','left','alpha',cellGroupsAlpha);
        else
            lateRespCells = ttest(...
                squeeze(mean(longTCAV(lateWinFr,:,:),1)),...
                squeeze(mean(longTCAV(basewin,:,:),1)),...
                'dim',2,'tail','right','alpha',cellGroupsAlpha);
            lateSuppCells = ttest(...
                squeeze(mean(longTCAV(lateWinFr,:,:),1)),...
                squeeze(mean(longTCAV(basewin,:,:),1)),...
                'dim',2,'tail','left','alpha',cellGroupsAlpha);
        end

        lateCycRespCutoffPass = mean(mean(lateTCAV(respwin,:,:),3),1) > minRespThreshold_decode;

        eaTarRespCells = sum(cell2mat(cellfun(@(x) ...
            ttest(squeeze(mean(x(respwin_target,:,:),1)),...
            squeeze(mean(x(basewin_0_target,:,:),1)),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha),...
            targetDataExpt(iexp).av(visualTrials).tc,'unif',0)),2) > 0;


        allTargetTC = [];
        for itar = 1:length(targetDataExpt(iexp).av(visualTrials).tc)
            allTargetTC = cat(3,allTargetTC,...
                targetDataExpt(iexp).av(visualTrials).tc{itar});
        end
        allTarRespCells = ttest(squeeze(mean(allTargetTC(respwin_target,:,:),1)),...
            squeeze(mean(allTargetTC(basewin_0_target,:,:),1)),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha);

        eaTarRespCutoffPass = sum(cell2mat(cellfun(@(x) mean(mean(x(respwin_target,:,:),3),1) - mean(mean(x(basewin_0_target,:,:),3),1),...
            targetDataExpt(iexp).av(visualTrials).tc,'unif',0)') > minRespThreshold_decode,1) > 0 ...
            | (mean(mean(allTargetTC(respwin_target,:,:),3),1) - mean(mean(allTargetTC(basewin_0_target,:,:),3),1)) > minRespThreshold_decode; 

        respCellsExpt(iexp).exptName = antiDataExpt(iexp).exptName;
        respCellsExpt(iexp).firstRespCells = firstRespCells;
        respCellsExpt(iexp).lateRespCells = lateRespCells;
        respCellsExpt(iexp).lateSuppCells = lateSuppCells;
        respCellsExpt(iexp).lateCycRespCells = lateCycRespCells;
        respCellsExpt(iexp).targetRespCells = eaTarRespCells | allTarRespCells;
    %     respCellsExpt(iexp).decodeAnalysisCells = ...
    %         (firstRespCells & mean(mean(firstTCAV(respwin,:,:),3),1)' > minRespThreshold)...
    %         | ((eaTarRespCells | allTarRespCells) & eaTarRespCutoffPass');
        respCellsExpt(iexp).decodeAnalysisCells = ...
            (lateCycRespCells & lateCycRespCutoffPass')...
            | ((eaTarRespCells | allTarRespCells) & eaTarRespCutoffPass');
        if shortCycExptInd(iexp)
            isShortCycExpt = cat(1,isShortCycExpt,true(length(firstRespCells),1));
        else
            isShortCycExpt = cat(1,isShortCycExpt,false(length(firstRespCells),1));
        end

    end
    
if doDecoding
    if strcmp(ds,'FSAV_attentionV1') | strcmp(ds,'FSAV_attentionV1_noAttn')
%         dcInfo = struct;
        decodeAnalysis = struct;
        decodeAnalysis.av(visualTrials).name = 'Visual';
        decodeAnalysis.av(auditoryTrials).name = 'Auditory';
        for iexp = 1:nexp
            rng(0) % cells randomized and trials randomized
%             cellInd = respCellsExpt(iexp).decodeAnalysisCells & ...
%                 oriTuningExpt(iexp).tuningReliability' <= tuningReliabilityThresh_decode;
            cellInd = respCellsExpt(iexp).lateCycRespCells|...
                respCellsExpt(iexp).targetRespCells;
            decodeAnalysis(iexp).nCells = sum(cellInd);            
            
%             decodeAnalysis(iexp).nCellsSelected = sum(cellInd);
            decodeAnalysis(iexp).cellInd = cellInd;
            respOther = cell(1,2);
            trOutOther = cell(1,2);
            trSampleIndOther = cell(1,2);
            detectGLMOther = cell(1,2);
            targetGLMOther = cell(1,2);
            if isfield(decodeDataExpt(iexp).av(visualTrials),'catchOutcome')
                data = cat(2,...
                    decodeDataExpt(iexp).av(visualTrials).resp(cellInd,:),...
                    decodeDataExpt(iexp).av(auditoryTrials).resp(cellInd,:),...
                    decodeDataExpt(iexp).av(visualTrials).catchResp(cellInd,:))';               
                
                [coeffAllCells,scoresAllCells,latentAllCells] = pca(data);                
                respAllCells_av_withInv = zscore(scoresAllCells);
                nvis = size(decodeDataExpt(iexp).av(visualTrials).resp,2);
                naud = size(decodeDataExpt(iexp).av(auditoryTrials).resp,2);
                ninv = size(decodeDataExpt(iexp).av(visualTrials).catchResp,2);
                respAllCells_vis = respAllCells_av_withInv(1:nvis,:);
                respAllCells_aud = respAllCells_av_withInv((nvis+1):(nvis+naud),:);
                catchResp = respAllCells_av_withInv((end-ninv+1):end,:);
                
                data4Shuff = cell(1,2);
                data4Shuff{visualTrials} = data(1:nvis,:);
                data4Shuff{auditoryTrials} = data((nvis+1):(nvis+naud),:);
                                
                stim_eaTr = cat(2,...
                    decodeDataExpt(iexp).av(visualTrials).stim,...
                    decodeDataExpt(iexp).av(auditoryTrials).stim,...
                    decodeDataExpt(iexp).av(visualTrials).catchStim(...
                    1:length(decodeDataExpt(iexp).av(visualTrials).catchOutcome)));
                stims_all = unique(stim_eaTr);
                mod_eaTr = cat(2,...
                    ones(1,length(decodeDataExpt(iexp).av(visualTrials).stim)),...
                    2.*ones(1,length(decodeDataExpt(iexp).av(auditoryTrials).stim)),...
                    3.*ones(1,length(decodeDataExpt(iexp).av(visualTrials).catchOutcome)));
                data_resid = nan(size(data));
                for istim = 1:length(stims_all)
                    for imod = 1:3
                        ind = stim_eaTr == stims_all(istim) & mod_eaTr == imod;
                        data_resid(ind,:) = data(ind,:) - mean(data(ind,:));
                    end
                end
                
                [coeffAllCells_resid,scoresAllCells_resid] = pca(data_resid);                
                respAllCells_av_withInv_resid = zscore(scoresAllCells_resid);
                respAllCells_vis_resid = respAllCells_av_withInv_resid(1:nvis,:);
                respAllCells_aud_resid = respAllCells_av_withInv_resid((nvis+1):(nvis+naud),:);
                catchResp_resid = respAllCells_av_withInv_resid((end-ninv+1):end,:);    
%                 [coeffAllCells_vis,scoresAllCells_vis,latentAllCells_vis] = pca(respAllCells_vis(:,cellInd));
%                 [coeffAllCells_aud,scoresAllCells_aud,latentAllCells_aud] = pca(respAllCells_aud(:,cellInd));
%                 [coeffAllCells_catch,scoresAllCells_catch,latentAllCells_catch] = pca(catchResp(:,cellInd));

            else
                data = cat(2,...
                    decodeDataExpt(iexp).av(visualTrials).resp(cellInd,:),...
                    decodeDataExpt(iexp).av(auditoryTrials).resp(cellInd,:))';
                
                [coeffAllCells,scoresAllCells,latentAllCells] = pca(data);
                respAllCells_av_withInv = zscore(scoresAllCells);
                nvis = size(decodeDataExpt(iexp).av(visualTrials).resp,2);
                naud = size(decodeDataExpt(iexp).av(auditoryTrials).resp,2);
                respAllCells_vis = respAllCells_av_withInv(1:nvis,:);
                respAllCells_aud = respAllCells_av_withInv((nvis+1):end,:);

                
                data4Shuff = cell(1,2);
                data4Shuff{visualTrials} = data(1:nvis,:);
                data4Shuff{auditoryTrials} = data((nvis+1):(nvis+naud),:);
%                 [coeffAllCells_vis,scoresAllCells_vis,latentAllCells_vis] = pca(respAllCells_vis(:,cellInd));
%                 [coeffAllCells_aud,scoresAllCells_aud,latentAllCells_aud] = pca(respAllCells_aud(:,cellInd));

                stim_eaTr = cat(2,...
                    decodeDataExpt(iexp).av(visualTrials).stim,...
                    decodeDataExpt(iexp).av(auditoryTrials).stim);
                stims_all = unique(stim_eaTr);
                mod_eaTr = cat(2,...
                    ones(1,length(decodeDataExpt(iexp).av(visualTrials).stim)),...
                    2.*ones(1,length(decodeDataExpt(iexp).av(auditoryTrials).stim)));
                data_resid = nan(size(data));
                for istim = 1:length(stims_all)
                    for imod = 1:2
                        ind = stim_eaTr == stims_all(istim) & mod_eaTr == imod;
                        data_resid(ind,:) = data(ind,:) - mean(data(ind,:),1);
                    end
                end
                
                [coeffAllCells_resid,scoresAllCells_resid] = pca(data_resid);                
                respAllCells_av_withInv_resid = zscore(scoresAllCells_resid);
                respAllCells_vis_resid = respAllCells_av_withInv_resid(1:nvis,:);
                respAllCells_aud_resid = respAllCells_av_withInv_resid((nvis+1):(nvis+naud),:); 
            end
            rng(0)
            for iav = 1:2
                trOut = decodeDataExpt(iexp).av(iav).outcome;

                if iav == 1
                    respAllCells = respAllCells_vis;
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
                elseif iav == 2
                    respAllCells = respAllCells_aud;
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,ampBins);
                end
                nStimPerBin = histcounts(trStimID);
                minBinN = min(nStimPerBin(nStimPerBin >= minTrN_mdl));
                respStimSort = cell(1,nStimBins);
                trOutStimSort = cell(1,nStimBins);
                for istim = 1:nStimBins
                    ind = find(trStimID == istim);
                    if length(ind) >= minTrN_mdl
                        if istim == 1
                            matchTrialsInd = [];
                            if sum(nStimPerBin >= minTrN_mdl) == 2
                                n = minBinN;
                            elseif minBinN == nStimPerBin(istim)
                                error('not enough FA/CR trials')
                            else
                                n = (nStimBins-1).*minBinN;
                                if n > length(ind)
                                    error('not enough FA/CR trials')
                                end
                            end
                            indSample = randsample(ind,n);
                            matchTrialsInd = cat(2,matchTrialsInd,indSample);
                        else
                            indSample = randsample(ind,minBinN);
                            matchTrialsInd = cat(2,matchTrialsInd,indSample);
                        end
                        respStimSort{istim} = respAllCells(indSample,:);
                        trOutStimSort{istim} = trOut(indSample);
                    end
                end
                nMatchedTrials = cumsum(cellfun(@length,trOutStimSort));
                for istim = 1:nStimBins
                    if istim == 1
                        stimSortInd = cell(1,nStimBins);
                        stimSortInd{istim} = 1:nMatchedTrials;
                    else
                        stimSortInd{istim} = ...
                            (nMatchedTrials(istim-1)+1):nMatchedTrials(istim);
                    end
                end
                trSampleIndOther{iav} = matchTrialsInd;
                decodeAnalysis(iexp).av(iav).respAllCells = respAllCells(matchTrialsInd,:);
                decodeAnalysis(iexp).av(iav).trOut = trOut(matchTrialsInd);
                
                respStimSort_av{iav} = respStimSort;
                trOutStimSort_av{iav} = trOutStimSort;
                stimSortInd_av{iav} = stimSortInd;
                matchTrialsInd_av{iav} = matchTrialsInd;
            end
            
            data4Shuff_matched = cellfun(@(x,y) x(y,:),data4Shuff,matchTrialsInd_av,'unif',0);
            trOut_matched = cellfun(@(x,y) x(y),...
                {decodeDataExpt(iexp).av(visualTrials).outcome,...
                decodeDataExpt(iexp).av(auditoryTrials).outcome},matchTrialsInd_av,'unif',0);
            [detectTrInd,targetTrInd] = cellfun(@(x) getStimAndBehaviorYs(x),...
                trOut_matched,'unif',0);
            dataShuff_detect = cellfun(@(x,y) shuffWithinCats(x,y),data4Shuff_matched,...
                detectTrInd,'unif',0);
            dataShuff_target = cellfun(@(x,y) shuffWithinCats(x,y),data4Shuff_matched,...
                targetTrInd,'unif',0);
            
            n = cellfun(@(x) size(x,1),data4Shuff_matched);
            d = cell2mat(cellfun(@(x) cell2mat(x'),dataShuff_detect,'unif',0)');
            [~,scoresAllCells] = pca(d);
            pcaShuff_zscore = zscore(scoresAllCells);
            pcaShuff_detect_matched = {pcaShuff_zscore(1:n(1),:),pcaShuff_zscore((n(1)+1):end,:)};
            d = cell2mat(cellfun(@(x) cell2mat(x'),dataShuff_target,'unif',0)');
            [~,scoresAllCells] = pca(d);
            pcaShuff_zscore = zscore(scoresAllCells);
            pcaShuff_target_matched = {pcaShuff_zscore(1:n(1),:),pcaShuff_zscore((n(1)+1):end,:)};
            
            nTrials = min(cellfun(@length,matchTrialsInd_av));
            if round(nTrials.*fracPCsUsed) > maxPCs
                if round(nTrials.*fracPCsUsed) > sum(cellInd)
                    nPCs = sum(cellInd);
                else
                    nPCs = maxPCs;
                end
            elseif sum(cellInd) < round(nTrials.*fracPCsUsed)
                nPCs = sum(cellInd);
            else
                nPCs = round(nTrials.*fracPCsUsed);
            end
            decodeAnalysis(iexp).nPCs = nPCs;
            decodeAnalysis(iexp).minTrials = nTrials;
            
            for iav = 1:2
                trOut = decodeDataExpt(iexp).av(iav).outcome;
                respStimSort = respStimSort_av{iav};
                trOutStimSort = trOutStimSort_av{iav};
                stimSortInd = stimSortInd_av{iav};
                matchTrialsInd = matchTrialsInd_av{iav};

                if iav == 1
                    respAllCells = respAllCells_vis;
%                     trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
                elseif iav == 2
                    respAllCells = respAllCells_aud;
%                     trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,ampBins);
                end
                
                emptyInd = cellfun(@isempty,respStimSort);
                respStimSort(~emptyInd) = cellfun(@(x) x(:,1:nPCs),...
                    respStimSort(~emptyInd),'unif',0);
                
                resp = respAllCells(matchTrialsInd,1:nPCs);
%                 resp = respAllCells(matchTrialsInd,cellInd);
                [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOut(matchTrialsInd));

                detectCorr = corr(detectTrInd,resp);
                targetCorr = corr(targetTrInd,resp);

                C = eye(size(resp,2));
                p=1;
                [~,~,detectGLM] = glmfit(resp*C,detectTrInd,'binomial');
                [~,~,targetGLM] = glmfit(resp*C,targetTrInd,'binomial');
                                
                detectWeight = detectGLM.beta(2:end);
                targetWeight = targetGLM.beta(2:end);

                dv_detect = mean(detectTrInd);
                dv_target = mean(targetTrInd);
                
                fprintf('Expt %s, hold-out analysis\n',num2str(iexp))
                pctCorrectDetect_train = getPctCorr_trainData(detectGLM,resp,detectTrInd,dv_detect);
                pctCorrectDetect_ho = getPctCorr_hoData(resp,detectTrInd,dv_detect);

                pctCorrectTarget_train = getPctCorr_trainData(targetGLM,resp,targetTrInd,dv_target);
                pctCorrectTarget_ho = getPctCorr_hoData(resp,targetTrInd,dv_target);
 
                pctCorrDetect_xStim_train = nan(1,nStimBins);
                pctCorrDetect_xStim_ho = nan(1,nStimBins);
                pctCorrTarget_xStim_train = nan(1,nStimBins);
                pctCorrTarget_xStim_ho = nan(1,nStimBins);
                for istim = 1:nStimBins
                        if isempty(trOutStimSort{istim})
                            continue
                        end
                        [detectStimInd, targetStimInd] = getStimAndBehaviorYs(...
                            trOutStimSort{istim});
                        pctCorrDetect_xStim_train(istim) = getPctCorr_trainData(...
                            detectGLM,respStimSort{istim},detectStimInd,dv_detect);
                        pctCorrDetect_xStim_ho(istim) = getPctCorr_hoData_subGroup(...
                            resp,detectTrInd,stimSortInd{istim},dv_detect);
                        pctCorrTarget_xStim_train(istim) = getPctCorr_trainData(...
                            targetGLM,respStimSort{istim},targetStimInd,dv_target);
                        pctCorrTarget_xStim_ho(istim) = getPctCorr_hoData_subGroup(...
                            resp,targetTrInd,stimSortInd{istim},dv_target);
                end

                fprintf('Expt %s, starting resp win analysis\n',num2str(iexp))
                nwins = size(decodeDataExpt(iexp).av(iav).movWinResp,3);
                pctCorrectDetect_train_respwin = nan(1,nwins);
                pctCorrectDetect_ho_respwin = nan(1,nwins);
                pctCorrectTarget_train_respwin = nan(1,nwins);
                pctCorrectTarget_ho_respwin = nan(1,nwins);
                for iwin = 1:nwins
%                     rAllCells = zscore(decodeDataExpt(iexp).av(iav).movWinResp(:,:,iwin)');
%                     r = rAllCells(matchTrialsInd,cellInd);
                    rAllCells = zscore(decodeDataExpt(iexp).av(iav).movWinResp(:,:,iwin)');
                    [~,sAllCells] = pca(rAllCells(:,cellInd));
                    r = sAllCells(matchTrialsInd,1:nPCs);
                    [~,~,detectGLM_temp] = glmfit(r*C,detectTrInd,'binomial');
                    [~,~,targetGLM_temp] = glmfit(r*C,targetTrInd,'binomial');

                    pctCorrectDetect_train_respwin(iwin) = getPctCorr_trainData(...
                        detectGLM_temp,r,detectTrInd,dv_detect);
                    pctCorrectDetect_ho_respwin(iwin) = getPctCorr_hoData(r,detectTrInd,dv_detect);
                    pctCorrectTarget_train_respwin(iwin) = getPctCorr_trainData(...
                        targetGLM_temp,r,targetTrInd,dv_target);
                    pctCorrectTarget_ho_respwin(iwin) = getPctCorr_hoData(r,targetTrInd,dv_target);
                end

                respOther{iav} = resp;
                trOutOther{iav} = trOut(matchTrialsInd);
                detectGLMOther{iav} = detectGLM;
                targetGLMOther{iav} = targetGLM;  
                
                % train the model on the PCs from the residuals
                if iav == 1
                    respAllCells_resid = respAllCells_vis_resid;
                elseif iav == 2
                    respAllCells_resid = respAllCells_aud_resid;
                end
                resp_resid = respAllCells_resid(matchTrialsInd,1:nPCs);
                C = eye(size(resp_resid,2));
                [~,~,detectGLM_resid] = glmfit(resp_resid*C,detectTrInd,'binomial');
                [~,~,targetGLM_resid] = glmfit(resp_resid*C,targetTrInd,'binomial');
                 
                pctCorrectDetect_train_resid = getPctCorr_trainData(detectGLM_resid,resp_resid,detectTrInd,dv_detect);
                pctCorrectDetect_ho_resid = getPctCorr_hoData(resp_resid,detectTrInd,dv_detect);

                pctCorrectTarget_train_resid = getPctCorr_trainData(targetGLM_resid,resp_resid,targetTrInd,dv_target);
                pctCorrectTarget_ho_resid = getPctCorr_hoData(resp_resid,targetTrInd,dv_target);
 
                detectWeight_cells_resid = coeffAllCells_resid(:,1:nPCs)*detectGLM_resid.beta(2:(nPCs+1));
                targetWeight_cells_resid = coeffAllCells_resid(:,1:nPCs)*targetGLM_resid.beta(2:(nPCs+1));
                
                % train the model with stim or choices as predictors
                resp_withStim = cat(2,resp,targetTrInd);
                resp_withChoice = cat(2,resp,detectTrInd);
                C = eye(size(resp_withStim,2));
                [~,~,detectGLM_withStim] = glmfit(resp_withStim*C,detectTrInd,'binomial');
                pctCorrectDetect_ho_withStim = getPctCorr_hoData(resp_withStim,detectTrInd,dv_detect);
                C = eye(size(resp_withChoice,2));
                [~,~,targetGLM_withChoice] = glmfit(resp_withChoice*C,targetTrInd,'binomial');
                pctCorrectTarget_ho_withChoice = getPctCorr_hoData(resp_withChoice,targetTrInd,dv_target);
                
                detectWeight_cells_withStim = coeffAllCells(:,1:nPCs)*detectGLM_withStim.beta(2:(nPCs+1));
                targetWeight_cells_withChoice = coeffAllCells(:,1:nPCs)*targetGLM_withChoice.beta(2:(nPCs+1));
                modeledStimWeight = detectGLM_withStim.beta(nPCs+2);
                modeledChoiceWeight = targetGLM_withChoice.beta(nPCs+2);
                
%                 pctCorrectDetect_ho_withStim_pcsOnly = getPctCorr_hoData_choosePCs(...
%                     resp_withStim,detectTrInd,dv_detect,1:nPCs);
%                 pctCorrectTarget_ho_withChoice_pcsOnly = getPctCorr_hoData_choosePCs(...
%                     resp_withChoice,targetTrInd,dv_target,1:nPCs);
                
                
                % train the model with fixed stim or choice predictors
                weight_choice = glmfit(detectTrInd,targetTrInd,'binomial');
                weight_stim = glmfit(targetTrInd,detectTrInd,'binomial');
                
                dv_choiceEaTrial = glmval(weight_choice,detectTrInd,'logit');
                dv_stimEaTrial = glmval(weight_stim,targetTrInd,'logit');
                
                resp_withFixedChoice = cat(2,resp,dv_choiceEaTrial);
                resp_withFixedStim = cat(2,resp,dv_stimEaTrial);
                
                C = eye(size(resp_withFixedChoice,2));
                [~,~,targetGLM_withFixedChoice] = glmfit(resp_withFixedChoice*C,targetTrInd,'binomial');
                pctCorrectTarget_ho_withFixedChoice = getPctCorr_hoData(resp_withFixedChoice,targetTrInd,dv_target);
                C = eye(size(resp_withFixedStim,2));
                [~,~,detectGLM_withFixedStim] = glmfit(resp_withFixedStim*C,detectTrInd,'binomial');
                pctCorrectDetect_ho_withFixedStim = getPctCorr_hoData(resp_withFixedStim,detectTrInd,dv_detect);
                
                detectGLMOther_withStim{iav} = detectGLM_withFixedStim;
                targetGLMOther_withChoice{iav} = targetGLM_withFixedChoice; 
                
                pctCorrectDetect_ho_withFixedStim_pcsOnly = getPctCorr_hoData_choosePCs(...
                    resp_withFixedStim,detectTrInd,dv_detect,1:nPCs);
                pctCorrectTarget_ho_withFixedChoice_pcsOnly = getPctCorr_hoData_choosePCs(...
                    resp_withFixedChoice,targetTrInd,dv_target,1:nPCs);
                
                
                
%                 C = eye(size(dv_stimEaTrial,2));
%                 [~,~,detectGLM_stimOnly] = glmfit(dv_stimEaTrial*C,detectTrInd,'binomial');
                pctCorrectDetect_stimOnly = getPctCorr_hoData(dv_stimEaTrial,detectTrInd,dv_detect);
%                 C = eye(size(dv_choiceEaTrial,2));
%                 [~,~,targetGLM_choiceOnly] = glmfit(dv_choiceEaTrial*C,detectTrInd,'binomial');
                pctCorrectTarget_choiceOnly = getPctCorr_hoData(dv_choiceEaTrial,targetTrInd,dv_target);
                
                % bootstrap trials used to generate stim:choice correlation
                % distribution
                choiceModel_fixedStim_struct = struct;
                [choiceModel_fixedStim_struct.stimChoiceCorr_boot,...
                    choiceModel_fixedStim_struct.trialsUsed_boot,...
                    choiceModel_fixedStim_struct.choiceModelWeights_boot,...
                    choiceModel_fixedStim_struct.pctCorr_choice_fixedStim_boot,...
                    choiceModel_fixedStim_struct.pctCorr_choice_fixedStim_pcsOnly_boot,...
                    choiceModel_fixedStim_struct.pctCorr_choice_fixedStimOnly_boot,...
                    choiceModel_fixedStim_struct.pctCorr_pcs,...
                    choiceModel_fixedStim_struct.weights_pcs] = ...
                    getSubSampledModel_boot(...
                    targetTrInd,detectTrInd,resp_withFixedStim,100,round(0.5*length(detectTrInd)));
                choiceModel_fixedStim_struct.choiceModelWeights_boot = ...
                    cellfun(@(x,y) cat(1,coeffAllCells(:,1:nPCs)*x(2:(end-1)),x(end)),...
                    choiceModel_fixedStim_struct.choiceModelWeights_boot,'unif',0);
                choiceModel_fixedStim_struct.weights_pcs = ...
                    cellfun(@(x,y) coeffAllCells(:,1:nPCs)*x(2:(end)),...
                    choiceModel_fixedStim_struct.weights_pcs,'unif',0);
                
                stimModel_fixedChoice_struct = struct;
                [stimModel_fixedChoice_struct.stimChoiceCorr_boot,...
                    stimModel_fixedChoice_struct.trialsUsed_boot,...
                    stimModel_fixedChoice_struct.stimModelWeights_boot,...
                    stimModel_fixedChoice_struct.pctCorr_stim_fixedChoice_boot,...
                    stimModel_fixedChoice_struct.pctCorr_stim_fixedChoice_pcsOnly_boot,...
                    stimModel_fixedChoice_struct.pctCorr_stim_fixedChoiceOnly_boot,...
                    stimModel_fixedChoice_struct.pctCorr_pcs,...
                    stimModel_fixedChoice_struct.weights_pcs] = ...
                    getSubSampledModel_boot(...
                    detectTrInd,targetTrInd,resp_withFixedChoice,100,round(0.5*length(targetTrInd)));
                stimModel_fixedChoice_struct.stimModelWeights_boot = ...
                    cellfun(@(x,y) cat(1,coeffAllCells(:,1:nPCs)*x(2:(end-1)),x(end)),...
                    stimModel_fixedChoice_struct.stimModelWeights_boot,'unif',0);
                stimModel_fixedChoice_struct.weights_pcs = ...
                    cellfun(@(x,y) coeffAllCells(:,1:nPCs)*x(2:(end)),...
                    stimModel_fixedChoice_struct.weights_pcs,'unif',0);
                
                % sub-sample choices (with replacement) to get a 50% choice 
                % dv for both targets and distractors, thus de-correlating
                trialSubsetIndex = getDV50_choiceAndStim(...
                    targetTrInd,detectTrInd,length(targetTrInd),minTrN);
                if isnan(trialSubsetIndex)
                    pctCorrectDetect_ho_matchAll = [];

                    pctCorrectTarget_ho_matchAll = [];

                    detectWeight_cells_matchAll = [];
                    targetWeight_cells_matchAll = [];
                    n_repTrials_matchAll = nan(1,4);
                    f_repTrials_matchAll = nan(1,2);
                else
                    C = eye(size(resp(trialSubsetIndex,:),2));
                    p=1;
                    [~,~,detectGLM_matchAll] = glmfit(resp(trialSubsetIndex,:)*C,detectTrInd(trialSubsetIndex),'binomial');
                    [~,~,targetGLM_matchAll] = glmfit(resp(trialSubsetIndex,:)*C,targetTrInd(trialSubsetIndex),'binomial');

                    dv_detect = mean(detectTrInd(trialSubsetIndex));
                    dv_target = mean(targetTrInd(trialSubsetIndex));

                    pctCorrectDetect_train_matchAll = getPctCorr_trainData(detectGLM,resp(trialSubsetIndex,:),detectTrInd(trialSubsetIndex),dv_detect);
                    pctCorrectDetect_ho_matchAll = getPctCorr_hoData(resp(trialSubsetIndex,:),detectTrInd(trialSubsetIndex),dv_detect);

                    pctCorrectTarget_train_matchAll = getPctCorr_trainData(targetGLM,resp(trialSubsetIndex,:),targetTrInd(trialSubsetIndex),dv_target);
                    pctCorrectTarget_ho_matchAll = getPctCorr_hoData(resp(trialSubsetIndex,:),targetTrInd(trialSubsetIndex),dv_target);

                    detectWeight_cells_matchAll = coeffAllCells(:,1:nPCs)*detectGLM_matchAll.beta(2:(nPCs+1));
                    targetWeight_cells_matchAll = coeffAllCells(:,1:nPCs)*targetGLM_matchAll.beta(2:(nPCs+1));
                    
                    [n_repTrials_matchAll,n_trTypeID] = getRepeatedTrialsN(...
                        trialSubsetIndex,targetTrInd,detectTrInd);
                    trStimID = discretize(...
                        decodeDataExpt(iexp).av(iav).stim(matchTrialsInd),oriBins);
                    [f_repTrials_matchAll,f_trTypeID] = getFractionEasyTrials(...
                        trialSubsetIndex,targetTrInd,detectTrInd,trStimID);
                    
                    respOther_matchAll{iav} = resp(trialSubsetIndex,:);
                    trOutOther_matchAll{iav} = trOutOther{iav}(trialSubsetIndex);
                    detectGLMOther_matchAll{iav} = detectGLM_matchAll;
                    targetGLMOther_matchAll{iav} = targetGLM_matchAll;
                end                
                
                % train trial shuffle model
                detectTrInd_shuff = cat(1,zeros(length(detectTrInd)-sum(detectTrInd),1),...
                    ones(sum(detectTrInd),1));
                targetTrInd_shuff = cat(1,zeros(length(targetTrInd)-sum(targetTrInd),1),...
                    ones(sum(targetTrInd),1));
                resp_detect_shuff = pcaShuff_detect_matched{iav}(:,1:nPCs);
                resp_target_shuff = pcaShuff_target_matched{iav}(:,1:nPCs);
                
                detectCorr_shuff = corr(detectTrInd_shuff,resp_detect_shuff);
                targetCorr_shuff = corr(targetTrInd_shuff,resp_target_shuff);

                C = eye(size(resp_detect_shuff,2));
                p=1;
                [~,~,detectGLM_shuff] = glmfit(resp_detect_shuff*C,detectTrInd,'binomial');
                [~,~,targetGLM_shuff] = glmfit(resp_target_shuff*C,targetTrInd,'binomial');
                                
                detectWeight_shuff = detectGLM_shuff.beta(2:end);
                targetWeight_shuff = targetGLM_shuff.beta(2:end);

                dv_detect_shuff = mean(detectTrInd_shuff);
                dv_target_shuff = mean(targetTrInd_shuff);
                
                fprintf('Expt %s, trial shuffle analysis\n',num2str(iexp))
                pctCorrectDetect_train_shuff = getPctCorr_trainData(detectGLM_shuff,...
                    resp_detect_shuff,detectTrInd_shuff,dv_detect_shuff);
                pctCorrectDetect_ho_shuff = getPctCorr_hoData(resp_detect_shuff,...
                    detectTrInd_shuff,dv_detect_shuff);
                pctCorrectTarget_train_shuff = getPctCorr_trainData(targetGLM_shuff,...
                    resp_target_shuff,targetTrInd_shuff,dv_target_shuff);
                pctCorrectTarget_ho_shuff = getPctCorr_hoData(resp_target_shuff,...
                    targetTrInd_shuff,dv_target_shuff);
                                
                respOther_shuff{iav} = {resp_detect_shuff,resp_target_shuff};
                trOutInd_shuff{iav} = {detectTrInd_shuff,targetTrInd_shuff};
                detectGLMOther_shuff{iav} = detectGLM_shuff;
                targetGLMOther_shuff{iav} = targetGLM_shuff;  
                
                
                % train detect model with distractor or target choices only
                if iav == 1
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
                elseif iav == 2
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,ampBins);
                end
%                 trOut_matched = trOut(matchTrialsInd);
                % distractors
                distInd = strcmp(trOut,'fa')|strcmp(trOut,'cr');
                trOut_dist = trOut(distInd);
                resp_dist = respAllCells(distInd,1:nPCs);
                [detectTrInd_dist, ~] = getStimAndBehaviorYs(trOut_dist);
                dv_detect_dist = mean(detectTrInd_dist);
                C = eye(size(resp_dist,2));
                p=1;
                [~,~,detectGLM_dist] = glmfit(resp_dist*C,detectTrInd_dist,'binomial');
                pctCorrectDetect_ho_dist = getPctCorr_hoData(resp_dist,detectTrInd_dist,dv_detect_dist);
                respOther_dist{iav} = resp_dist;
                detectGLMOther_dist{iav} = detectGLM_dist;
                trOutOther_dist{iav} = trOut_dist;
                detectWeight_cells_distOnly = coeffAllCells(:,1:nPCs)*detectGLM_dist.beta(2:end);
                detectWeight_pcs_distOnly = detectGLM_dist.beta(2:end);
                
                % targets
                tarInd = strcmp(trOut,'h')|strcmp(trOut,'m');
                trOut_tar = trOut(tarInd);
                resp_tar = respAllCells(tarInd,1:nPCs);
                
                [detectTrInd_tar, ~] = getStimAndBehaviorYs(trOut_tar);
                dv_detect_tar = mean(detectTrInd_tar);
                C = eye(size(resp_tar,2));
                p=1;
                [~,~,detectGLM_tar] = glmfit(resp_tar*C,detectTrInd_tar,'binomial');
                pctCorrectDetect_ho_tar = getPctCorr_hoData(resp_tar,detectTrInd_tar,dv_detect_tar);
                respOther_tar{iav} = resp_tar;
                detectGLMOther_tar{iav} = detectGLM_tar;
                trOutOther_tar{iav} = trOut_tar;
                detectWeight_cells_tarOnly = coeffAllCells(:,1:nPCs)*detectGLM_tar.beta(2:end);
                
                
                detectWeight_cells = coeffAllCells(:,1:nPCs)*detectWeight;
                targetWeight_cells = coeffAllCells(:,1:nPCs)*targetWeight;

                decodeAnalysis(iexp).av(iav).dvDetect = dv_detect;
                decodeAnalysis(iexp).av(iav).dvTarget = dv_target;
                decodeAnalysis(iexp).av(iav).correlationDetect = detectCorr;
                decodeAnalysis(iexp).av(iav).correlationTarget = targetCorr;
                decodeAnalysis(iexp).av(iav).weightDetect = detectWeight_cells;
                decodeAnalysis(iexp).av(iav).weightTarget = targetWeight_cells;
                decodeAnalysis(iexp).av(iav).weightDetect_pcs = detectWeight;
                decodeAnalysis(iexp).av(iav).weightTarget_pcs = targetWeight;
                decodeAnalysis(iexp).av(iav).pctCorrectAllDetect_train = pctCorrectDetect_train;
                decodeAnalysis(iexp).av(iav).pctCorrectAllDetect_holdout = pctCorrectDetect_ho;
                decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_train = pctCorrectTarget_train;
                decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_holdout = pctCorrectTarget_ho;
                decodeAnalysis(iexp).av(iav).pctCorrectXStimDetect_train = pctCorrDetect_xStim_train;
                decodeAnalysis(iexp).av(iav).pctCorrectXStimDetect_holdout = pctCorrDetect_xStim_ho;
                decodeAnalysis(iexp).av(iav).pctCorrectXStimTarget_train = pctCorrTarget_xStim_train;
                decodeAnalysis(iexp).av(iav).pctCorrectXStimTarget_holdout = pctCorrTarget_xStim_ho;
                decodeAnalysis(iexp).av(iav).pctCorrectDetectMovRespWin_train = pctCorrectDetect_train_respwin;
                decodeAnalysis(iexp).av(iav).pctCorrectDetectMovRespWin_holdout = pctCorrectDetect_ho_respwin;
                decodeAnalysis(iexp).av(iav).pctCorrectTargetMovRespWin_train = pctCorrectTarget_train_respwin;
                decodeAnalysis(iexp).av(iav).pctCorrectTargetMovRespWin_holdout = pctCorrectTarget_ho_respwin;
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_distOnly_holdout = pctCorrectDetect_ho_dist;
                decodeAnalysis(iexp).av(iav).weightDetect_distOnly = detectWeight_cells_distOnly;
                decodeAnalysis(iexp).av(iav).weightDetect_pcs_distOnly = detectWeight_pcs_distOnly;
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_tarOnly_holdout = pctCorrectDetect_ho_tar;
                decodeAnalysis(iexp).av(iav).weightDetect_tarOnly = detectWeight_cells_tarOnly;
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_shuff_holdout = pctCorrectDetect_ho_shuff;
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_shuff_holdout = pctCorrectTarget_ho_shuff;
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_withStim_holdout = pctCorrectDetect_ho_withStim;
                decodeAnalysis(iexp).av(iav).weightDetect_withStim = detectWeight_cells_withStim;
                decodeAnalysis(iexp).av(iav).modeledStimWeight = modeledStimWeight;
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_withChoice_holdout = pctCorrectTarget_ho_withChoice;
                decodeAnalysis(iexp).av(iav).weightTarget_withChoice = targetWeight_cells_withChoice;
                decodeAnalysis(iexp).av(iav).modeledChoiceWeight = modeledChoiceWeight;
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_withFixedStim_holdout = pctCorrectDetect_ho_withFixedStim;
                decodeAnalysis(iexp).av(iav).weights_fixedStimGivenChoice_DT = [min(dv_stimEaTrial) max(dv_stimEaTrial)];
                decodeAnalysis(iexp).av(iav).weightDetect_withFixedStim = coeffAllCells(:,1:nPCs)*detectGLM_withFixedStim.beta(2:end-1);
                decodeAnalysis(iexp).av(iav).modeledStimWeight_fixed = detectGLM_withFixedStim.beta(end);
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_withFixedChoice_holdout = pctCorrectTarget_ho_withFixedChoice;
                decodeAnalysis(iexp).av(iav).weights_fixedChoiceGivenStim_NY = [min(dv_choiceEaTrial) max(dv_choiceEaTrial)];
                decodeAnalysis(iexp).av(iav).weightTarget_withFixedChoice = coeffAllCells(:,1:nPCs)*targetGLM_withFixedChoice.beta(2:end-1);
                decodeAnalysis(iexp).av(iav).modeledChoiceWeight_fixed = targetGLM_withFixedChoice.beta(end);
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_withFixedStim_holdout_pcsOnly = ...
                    pctCorrectDetect_ho_withFixedStim_pcsOnly;
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_withFixedChoice_holdout_pcsOnly = ...
                    pctCorrectTarget_ho_withFixedChoice_pcsOnly;
                
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_choiceOnly = pctCorrectTarget_choiceOnly;
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_stimOnly = pctCorrectDetect_stimOnly;
                
                decodeAnalysis(iexp).av(iav).choiceModel_fixedStim_bootStruct = choiceModel_fixedStim_struct;
                decodeAnalysis(iexp).av(iav).stimModel_fixedChoice_bootStruct = stimModel_fixedChoice_struct;
                
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_resid = pctCorrectTarget_ho_resid;
                decodeAnalysis(iexp).av(iav).weightTarget_resid = targetWeight_cells_resid;
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_resid = pctCorrectDetect_ho_resid;
                decodeAnalysis(iexp).av(iav).weightDetect_resid = detectWeight_cells_resid;
                
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_matchAll = pctCorrectTarget_ho_matchAll;
                decodeAnalysis(iexp).av(iav).weightTarget_matchAll = targetWeight_cells_matchAll;
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_matchAll = pctCorrectDetect_ho_matchAll;
                decodeAnalysis(iexp).av(iav).weightDetect_matchAll = detectWeight_cells_matchAll;
                decodeAnalysis(iexp).av(iav).nRepeatTrials_fa_cr_h_m = n_repTrials_matchAll;
                decodeAnalysis(iexp).av(iav).fractionEasy_h_m = f_repTrials_matchAll;
                
%                 decodeAnalysis(iexp).av(iav).validMatchedPctCorrectDetect_holdout = validCatchMatchDetect;
%                 decodeAnalysis(iexp).av(iav).validMatchedPctCorrectTarget_holdout = validCatchMatchTarget;
            end
            fprintf('Expt %s: Testing opposite model...\n',num2str(iexp))
            for iav = 1:2
                if iav == 1
                    otherAV = 2;
                    if isfield(decodeDataExpt(iexp).av(iav),'catchOutcome')
                        [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{otherAV});
                        [catchDetectInd,catchTargetInd] = getStimAndBehaviorYs(...
                            decodeDataExpt(iexp).av(iav).catchOutcome);
%                         catchTargetResp = catchResp(:,cellInd);
                        catchTargetResp = catchResp(:,1:nPCs);
                        nCatch = length(catchDetectInd);

                            [catchPctCorrectDetect,catchCorrectTrialDetect] = getPctCorr_trainData(...
                                detectGLMOther{iav},catchTargetResp,catchDetectInd,...
                                decodeAnalysis(iexp).av(iav).dvDetect);
                            [catchPctCorrectTarget,catchCorrectTrialTarget] = getPctCorr_trainData(...
                                targetGLMOther{iav},catchTargetResp,catchTargetInd,...
                                decodeAnalysis(iexp).av(iav).dvTarget);
                            [catchPctCorrectDetect_aud,catchCorrectTrialDetect_aud] = getPctCorr_trainData(...
                                detectGLMOther{auditoryTrials},catchTargetResp,...
                                catchDetectInd,...
                                decodeAnalysis(iexp).av(auditoryTrials).dvDetect);            
                            [catchPctCorrectTarget_aud,catchCorrectTrialTarget_aud] = getPctCorr_trainData(...
                                targetGLMOther{auditoryTrials},catchTargetResp,...
                                catchTargetInd,...
                                decodeAnalysis(iexp).av(auditoryTrials).dvTarget);
%                             catchPctCorrectDetect = getPctCorr_trainData(...
%                                 detectGLMOther{iav},catchRespBalanced,catchDetectIndBalanced,...
%                                 decodeAnalysis(iexp).av(iav).dvDetect);
%                             catchPctCorrectTarget = getPctCorr_trainData(...
%                                 targetGLMOther{iav},catchRespBalanced,catchTargetIndBalanced,...
%                                 decodeAnalysis(iexp).av(iav).dvTarget);
%                             catchPctCorrectDetect_audModel = getPctCorr_trainData(...
%                                 detectGLMOther{auditoryTrials},catchRespBalanced,...
%                                 catchDetectIndBalanced,...
%                                 decodeAnalysis(iexp).av(auditoryTrials).dvDetect);            
%                             catchPctCorrectTarget_audModel = getPctCorr_trainData(...
%                                 targetGLMOther{auditoryTrials},catchRespBalanced,...
%                                 catchTargetIndBalanced,...
%                                 decodeAnalysis(iexp).av(auditoryTrials).dvTarget);
%                         else
%                             catchPctCorrectDetect = nan;
%                             catchCorrectTrialDetect = nan(0,1);
%                             catchPctCorrectTarget = nan;
%                             catchCorrectTrialTarget = nan(0,1);
%                             catchPctCorrectDetect_audModel = nan;
%                             catchCorrectTrialDetect_audModel = nan(0,1);        
%                             catchPctCorrectTarget_audModel = nan;
%                             catchCorrectTrialTarget_audModel = nan(0,1);  
%                         end
                        

                        trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
                        matchedTrialsID = trStimID(trSampleIndOther{iav});
                        resp_val = respOther{iav};
%                         catchStimID = cat(2,ones(1,nCatch),...
%                             discretize(decodeDataExpt(iexp).av(iav).catchStim,oriBins));
                        catchStimID = discretize(decodeDataExpt(iexp).av(iav).catchStim,oriBins);
                        
                        [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{iav});
                        if ~any(ismember(matchedTrialsID,unique(catchStimID)))
                            catchMatchTrials = [];
                        else
                            catchMatchTrials = cell2mat(getMatchedValidTrialIndex(...
                                matchedTrialsID,catchStimID)); 
                        end
                        [validCatchMatchDetect,validCatchCorrectTrialDetect] = getPctCorr_hoData_subGroup(...
                            resp_val,detectTrInd,catchMatchTrials,decodeAnalysis(iexp).av(iav).dvDetect);  
                        [validCatchMatchTarget,validCatchCorrectTrialTarget] = getPctCorr_hoData_subGroup(...
                            resp_val,targetTrInd,catchMatchTrials,decodeAnalysis(iexp).av(iav).dvTarget);
                        
%                         trStimID = discretize(decodeDataExpt(iexp).av(auditoryTrials).stim,ampBins);
%                         matchedTrialsID = trStimID(trSampleIndOther{auditoryTrials});
%                         matchedTrialsID(matchedTrialsID > 1) = 2;
%                         catchStimID(catchStimID > 1) = 2;
%                         [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{auditoryTrials});
%                         catchMatchTrials = cell2mat(getMatchedValidTrialIndex(...
%                             matchedTrialsID,catchStimID));
                        [validCatchMatchDetect_aud,validCatchCorrectTrialDetect_aud] = getPctCorr_trainData(...
                            detectGLMOther{otherAV},resp_val(catchMatchTrials,:),...
                            detectTrInd(catchMatchTrials),decodeAnalysis(iexp).av(otherAV).dvDetect);
                        [validCatchMatchTarget_aud,validCatchCorrectTrialTarget_aud] = getPctCorr_trainData(...
                            targetGLMOther{otherAV},resp_val(catchMatchTrials,:),...
                            targetTrInd(catchMatchTrials),decodeAnalysis(iexp).av(otherAV).dvTarget); 
                       
                    else
                        validCatchMatchDetect = [];
                        validCatchCorrectTrialDetect = [];
                        validCatchMatchTarget = [];
                        validCatchCorrectTrialTarget = [];
                        validCatchMatchDetect_aud = [];
                        validCatchCorrectTrialDetect_aud = [];
                        validCatchMatchTarget_aud = [];
                        validCatchCorrectTrialTarget_aud = [];
                        catchPctCorrectDetect = [];
                        catchPctCorrectTarget = [];
                        catchPctCorrectDetect_aud = [];
                        catchPctCorrectTarget_aud = [];
                        catchCorrectTrialDetect = [];
                        catchCorrectTrialTarget = [];
                        catchCorrectTrialDetect_aud = [];
                        catchCorrectTrialTarget_aud = [];
                        
                    end
                elseif iav == 2
                    otherAV = 1;
                    validCatchMatchDetect = [];
                    validCatchMatchTarget = [];
                    validCatchMatchDetect_aud = [];
                    validCatchMatchTarget_aud = [];
                    catchPctCorrectDetect = [];
                    catchPctCorrectTarget = [];
                    catchPctCorrectDetect_aud = [];
                    catchPctCorrectTarget_aud = [];
                end
                resp = respOther{otherAV};
                [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{otherAV});

%                 dv_detect = mean(detectTrInd);
%                 dv_target = mean(targetTrInd);
                [detecttrind4model,targettrind4model] = getStimAndBehaviorYs(trOutOther{iav});
                dv_detect = mean(detecttrind4model);
                dv_target = mean(targettrind4model);
                
                pctCorrectDetect = getPctCorr_trainData(detectGLMOther{iav},resp,detectTrInd,dv_detect);
                pctCorrectTarget = getPctCorr_trainData(targetGLMOther{iav},resp,targetTrInd,dv_target);
                
                stimIndTarget = targetTrInd == 0;
                pctCorrectDetectxStim = nan(1,2);
                pctCorrectDetectxStim(1) = getPctCorr_trainData(...
                    detectGLMOther{iav},resp(stimIndTarget,:),...
                    detectTrInd(stimIndTarget),dv_detect);
                pctCorrectDetectxStim(2) = getPctCorr_trainData(...
                    detectGLMOther{iav},resp(~stimIndTarget,:),...
                    detectTrInd(~stimIndTarget),dv_detect);
                pctCorrectTargetxStim = nan(1,2);
                pctCorrectTargetxStim(1) = getPctCorr_trainData(...
                    targetGLMOther{iav},resp(stimIndTarget,:),...
                    targetTrInd(stimIndTarget),dv_target);
                pctCorrectTargetxStim(2) = getPctCorr_trainData(...
                    targetGLMOther{iav},resp(~stimIndTarget,:),...
                    targetTrInd(~stimIndTarget),dv_target);
                
                % test other model with fixedstim/choice predictors
                pctCorrectDetect_withFixedStim = getPctCorr_weightsOnly(...
                    detectGLMOther_withStim{iav}.beta(1:(nPCs+1)),resp,...
                    detectTrInd,dv_detect);
                pctCorrectTarget_withFixedChoice = getPctCorr_weightsOnly(...
                    targetGLMOther_withChoice{iav}.beta(1:(nPCs+1)),resp,...
                    targetTrInd,dv_target);
                
                % test other model with balanced/uncorreleted choice
                % (_matchAll)
                if isempty(decodeAnalysis(iexp).av(visualTrials).pctCorrectDetect_matchAll)||...
                        isempty(decodeAnalysis(iexp).av(auditoryTrials).pctCorrectDetect_matchAll)
                    
                    pctCorrectDetect_matchAll = [];
                    pctCorrectTarget_matchAll = [];
                
                else
                    resp = respOther_matchAll{otherAV};
                    [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther_matchAll{otherAV});

                    [detecttrind4model,targettrind4model] = getStimAndBehaviorYs(trOutOther{iav});
                    dv_detect = mean(detecttrind4model);
                    dv_target = mean(targettrind4model);

                    
                    pctCorrectDetect_matchAll = getPctCorr_trainData(detectGLMOther_matchAll{iav},resp,detectTrInd,dv_detect);
                    pctCorrectTarget_matchAll = getPctCorr_trainData(targetGLMOther_matchAll{iav},resp,targetTrInd,dv_target);
                
                end
                
                % test other distractor only detect model
                resp = respOther_dist{otherAV};
                [detectTrInd, ~] = getStimAndBehaviorYs(trOutOther_dist{otherAV});
                dv_detect = mean(getStimAndBehaviorYs(trOutOther_dist{iav}));
                pctCorrectDetect_dist = getPctCorr_trainData(detectGLMOther_dist{iav},resp,detectTrInd,dv_detect);
                % test other target only detect model
                resp = respOther_tar{otherAV};
                [detectTrInd, ~] = getStimAndBehaviorYs(trOutOther_tar{otherAV});
                dv_detect = mean(getStimAndBehaviorYs(trOutOther_tar{iav}));
                pctCorrectDetect_tar = getPctCorr_trainData(detectGLMOther_tar{iav},resp,detectTrInd,dv_detect);
                
                % test other shuffled model
                [resp_detect,resp_target] = deal(respOther_shuff{otherAV}{:});
                [detectTrInd,targetTrInd] = deal(trOutInd_shuff{otherAV}{:});
                dv_detect = mean(trOutInd_shuff{iav}{1});
                dv_target = mean(trOutInd_shuff{iav}{2});    
                pctCorrectDetect_shuff = getPctCorr_trainData(...
                    detectGLMOther_shuff{iav},resp_detect,detectTrInd,dv_detect);
                pctCorrectTarget_shuff = getPctCorr_trainData(...
                    targetGLMOther_shuff{iav},resp_target,targetTrInd,dv_target);
                
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_otherAV = pctCorrectDetect;
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_otherAV = pctCorrectTarget;
                decodeAnalysis(iexp).av(iav).pctCorrectDetectxStim_otherAV = pctCorrectDetectxStim;
                decodeAnalysis(iexp).av(iav).pctCorrectTargetxStim_otherAV = pctCorrectTargetxStim;
                decodeAnalysis(iexp).av(iav).validMatchedPctCorrectDetect_holdout = validCatchMatchDetect;
                decodeAnalysis(iexp).av(iav).validMatchedPctCorrectTarget_holdout = validCatchMatchTarget;
                decodeAnalysis(iexp).av(iav).invalidPctCorrectDetect = catchPctCorrectDetect;
                decodeAnalysis(iexp).av(iav).invalidPctCorrectTarget = catchPctCorrectTarget;
                decodeAnalysis(iexp).av(iav).invalidPctCorrectDetect_testAudModel = catchPctCorrectDetect_aud;
                decodeAnalysis(iexp).av(iav).invalidPctCorrectTarget_testAudModel = catchPctCorrectTarget_aud;
                decodeAnalysis(iexp).av(iav).validMatchedPctCorrectDetect_testAudModel = validCatchMatchDetect_aud;
                decodeAnalysis(iexp).av(iav).validMatchedPctCorrectTarget_testAudModel = validCatchMatchTarget_aud;
                
                decodeAnalysis(iexp).av(iav).validMatchedCorrectTrialsDetect = validCatchCorrectTrialDetect;
                decodeAnalysis(iexp).av(iav).validMatchedCorrectTrialsTarget = validCatchCorrectTrialTarget;
                decodeAnalysis(iexp).av(iav).invalidCorrectTrialsDetect = catchCorrectTrialDetect;
                decodeAnalysis(iexp).av(iav).invalidCorrectTrialsTarget = catchCorrectTrialTarget;
                decodeAnalysis(iexp).av(iav).invalidCorrectTrialsDetect_testAud = catchCorrectTrialDetect_aud;
                decodeAnalysis(iexp).av(iav).invalidCorrectTrialsTarget_testAud = catchCorrectTrialTarget_aud;
                decodeAnalysis(iexp).av(iav).validMatchedCorrectTrialsDetect_testAud = validCatchCorrectTrialDetect_aud;
                decodeAnalysis(iexp).av(iav).validMatchedCorrectTrialsTarget_testAud = validCatchCorrectTrialTarget_aud;
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_distOnly_otherAV = pctCorrectDetect_dist;
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_tarOnly_otherAV = pctCorrectDetect_tar;
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_shuff_otherAV = pctCorrectDetect_shuff;
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_shuff_otherAV = pctCorrectTarget_shuff;
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_withFixedStim_otherAV = ...
                    pctCorrectDetect_withFixedStim;
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_withFixedChoice_otherAV = ...
                    pctCorrectTarget_withFixedChoice;
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_matchAll_otherAV = ...
                    pctCorrectDetect_matchAll;
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_matchAll_otherAV = ...
                    pctCorrectTarget_matchAll;
                
            end
            rng(0)
            randInd = cellfun(@(x) randsample(length(x),floor(length(x)/2)),...
                trOutOther,'unif',0);
            otherRandInd = cellfun(@(x,y) setdiff(1:length(x),y),trOutOther,randInd,'unif',0);
            nt = cellfun(@length,randInd);
            nt_other = cellfun(@length,otherRandInd);
%             nt = cellfun(@length,trOutOther);
            respAV = cat(1,respOther{visualTrials}(randInd{visualTrials},:),...
                respOther{auditoryTrials}(randInd{auditoryTrials},:));
%             respAV = cat(1,respOther{visualTrials},...
%                 respOther{auditoryTrials});
            [detectTrInd, targetTrInd] = getStimAndBehaviorYs(cat(2,...
                trOutOther{visualTrials}(randInd{visualTrials}),trOutOther{auditoryTrials}(randInd{auditoryTrials})));
            [detectTrInd_other, targetTrInd_other] = getStimAndBehaviorYs(cat(2,...
                trOutOther{visualTrials}(otherRandInd{visualTrials}),trOutOther{auditoryTrials}(otherRandInd{auditoryTrials})));
%             [detectTrInd, targetTrInd] = getStimAndBehaviorYs(cat(2,...
%                 trOutOther{visualTrials},trOutOther{auditoryTrials}));
            dv_detect = mean(detectTrInd);
            dv_target = mean(targetTrInd);
            C = eye(size(respAV,2));
            p=1;
            [~,~,detectGLM] = glmfit(respAV*C,detectTrInd,'binomial');
            [~,~,targetGLM] = glmfit(respAV*C,targetTrInd,'binomial');
            pctCorrectDetect_combo_vis = cmbPctCorrect(nt(visualTrials),...
                getPctCorr_hoData_subGroup(respAV,detectTrInd,1:nt(visualTrials),dv_detect),...
                nt_other(visualTrials),...
                getPctCorr_trainData(detectGLM,respOther{visualTrials}(otherRandInd{visualTrials},:),...
                    detectTrInd_other(1:nt_other),dv_detect));
            pctCorrectTarget_combo_vis = cmbPctCorrect(nt(visualTrials),...
                getPctCorr_hoData_subGroup(respAV,targetTrInd,1:nt(visualTrials),dv_target),...
                nt_other(visualTrials),...
                getPctCorr_trainData(targetGLM,respOther{visualTrials}(otherRandInd{visualTrials},:),...
                    targetTrInd_other(1:nt_other),dv_target));
            pctCorrectDetect_combo_aud = cmbPctCorrect(nt(auditoryTrials),...
                getPctCorr_hoData_subGroup(respAV,detectTrInd,(nt(visualTrials)+1):sum(nt),dv_detect),...
                nt_other(auditoryTrials),...
                getPctCorr_trainData(detectGLM,respOther{auditoryTrials}(otherRandInd{auditoryTrials},:),...
                    detectTrInd_other((nt(visualTrials)+1):sum(nt)),dv_detect));
            pctCorrectTarget_combo_aud = cmbPctCorrect(nt(auditoryTrials),...
                getPctCorr_hoData_subGroup(respAV,targetTrInd,(nt(visualTrials)+1):sum(nt),dv_target),...
                nt_other(auditoryTrials),...
                getPctCorr_trainData(targetGLM,respOther{auditoryTrials}(otherRandInd{auditoryTrials},:),...
                    targetTrInd_other((nt(visualTrials)+1):sum(nt)),dv_target));
                
%             pctCorrectDetect_combo_vis = getPctCorr_hoData_subGroup(respAV,detectTrInd,1:nt(visualTrials),dv_detect);
%             pctCorrectTarget_combo_vis = getPctCorr_hoData_subGroup(respAV,targetTrInd,1:nt(visualTrials),dv_target);
%             pctCorrectDetect_combo_aud = getPctCorr_hoData_subGroup(respAV,detectTrInd,(nt(visualTrials)+1):sum(nt),dv_detect);
%             pctCorrectTarget_combo_aud = getPctCorr_hoData_subGroup(respAV,targetTrInd,(nt(visualTrials)+1):sum(nt),dv_target);
            
            [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{visualTrials});
            pctCorrectDetect_vis = getPctCorr_hoData_subGroup(...
                respOther{visualTrials},...
                detectTrInd,randInd{visualTrials},mean(detectTrInd));
            pctCorrectTarget_vis = getPctCorr_hoData_subGroup(...
                respOther{visualTrials},...
                targetTrInd,randInd{visualTrials},mean(targetTrInd));
            [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{auditoryTrials});
            pctCorrectDetect_aud = getPctCorr_hoData_subGroup(...
                respOther{auditoryTrials},...
                detectTrInd,randInd{auditoryTrials},mean(detectTrInd));
            pctCorrectTarget_aud = getPctCorr_hoData_subGroup(...
                respOther{auditoryTrials},...
                targetTrInd,randInd{auditoryTrials},mean(targetTrInd));
            
            
            decodeAnalysis(iexp).comboTrainWeightDetect = coeffAllCells(:,1:nPCs)*detectGLM.beta(2:end);
            decodeAnalysis(iexp).comboTrainWeightTarget = coeffAllCells(:,1:nPCs)*targetGLM.beta(2:end);
            decodeAnalysis(iexp).av(visualTrials).pctCorrectDetect_comboTrain = pctCorrectDetect_combo_vis;
            decodeAnalysis(iexp).av(visualTrials).pctCorrectTarget_comboTrain = pctCorrectTarget_combo_vis;
            decodeAnalysis(iexp).av(auditoryTrials).pctCorrectDetect_comboTrain = pctCorrectDetect_combo_aud;
            decodeAnalysis(iexp).av(auditoryTrials).pctCorrectTarget_comboTrain = pctCorrectTarget_combo_aud;
            decodeAnalysis(iexp).av(visualTrials).pctCorrectDetect_visTrainComboMatch = pctCorrectDetect_vis;
            decodeAnalysis(iexp).av(visualTrials).pctCorrectTarget_visTrainComboMatch = pctCorrectTarget_vis;
            decodeAnalysis(iexp).av(auditoryTrials).pctCorrectDetect_audTrainComboMatch = pctCorrectDetect_aud;
            decodeAnalysis(iexp).av(auditoryTrials).pctCorrectTarget_audTrainComboMatch = pctCorrectTarget_aud;
        end

    nPCsSelected = cell2mat({decodeAnalysis.nPCs});    
    fprintf('%s +/- %s (%s-%s) pcs selected\n',num2str(mean(nPCsSelected)),...
        num2str(ste(nPCsSelected,2)),num2str(min(nPCsSelected)),num2str(max(nPCsSelected)))
    elseif strcmp(ds,'FSAV_V1_audControl')
        
    else
        decodeAnalysis = struct;
        decodeAnalysis.av(visualTrials).name = 'Visual';
        decodeAnalysis.av(auditoryTrials).name = 'Auditory';
        for iexp = 1:nexp
            rng(0) % cells randomized and trials randomized
            cellInd = respCellsExpt(iexp).lateCycRespCells|...
                respCellsExpt(iexp).targetRespCells;
            decodeAnalysis(iexp).nCells = sum(cellInd);            
            
            decodeAnalysis(iexp).cellInd = cellInd;
            respOther = cell(1,2);
            trOutOther = cell(1,2);
            trSampleIndOther = cell(1,2);
            targetGLMOther = cell(1,2);
            if isfield(decodeDataExpt(iexp).av(visualTrials),'catchOutcome')
                [coeffAllCells,scoresAllCells,latentAllCells] = pca(cat(2,...
                    decodeDataExpt(iexp).av(visualTrials).resp(cellInd,:),...
                    decodeDataExpt(iexp).av(auditoryTrials).resp(cellInd,:),...
                    decodeDataExpt(iexp).av(visualTrials).catchResp(cellInd,:))');                
                respAllCells_av_withInv = zscore(scoresAllCells);
                nvis = size(decodeDataExpt(iexp).av(visualTrials).resp,2);
                naud = size(decodeDataExpt(iexp).av(auditoryTrials).resp,2);
                ninv = size(decodeDataExpt(iexp).av(visualTrials).catchResp,2);
                respAllCells_vis = respAllCells_av_withInv(1:nvis,:);
                respAllCells_aud = respAllCells_av_withInv((nvis+1):(nvis+naud),:);
                catchResp = respAllCells_av_withInv((end-ninv+1):end,:);
            else
                [coeffAllCells,scoresAllCells,latentAllCells] = pca(cat(2,...
                    decodeDataExpt(iexp).av(visualTrials).resp(cellInd,:),...
                    decodeDataExpt(iexp).av(auditoryTrials).resp(cellInd,:))');
                respAllCells_av_withInv = zscore(scoresAllCells);
                nvis = size(decodeDataExpt(iexp).av(visualTrials).resp,2);
                naud = size(decodeDataExpt(iexp).av(auditoryTrials).resp,2);
                respAllCells_vis = respAllCells_av_withInv(1:nvis,:);
                respAllCells_aud = respAllCells_av_withInv((nvis+1):end,:);

            end
            
            for iav = 1:2
                trOut = decodeDataExpt(iexp).av(iav).outcome;

                if iav == 1
                    respAllCells = respAllCells_vis;
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
                elseif iav == 2
                    respAllCells = respAllCells_aud;
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,ampBins);
                end
                nStimPerBin = histcounts(trStimID);
                minBinN = min(nStimPerBin(nStimPerBin >= minTrN_mdl));
                respStimSort = cell(1,nStimBins);
                trOutStimSort = cell(1,nStimBins);
                for istim = 1:nStimBins
                    ind = find(trStimID == istim);
                    if length(ind) >= minTrN_mdl
                        if istim == 1
                            matchTrialsInd = [];
                            if sum(nStimPerBin >= minTrN_mdl) == 2
                                n = minBinN;
                            elseif minBinN == nStimPerBin(istim)
                                error('not enough FA/CR trials')
                            else
                                n = (nStimBins-1).*minBinN;
                                if n > length(ind)
                                    error('not enough FA/CR trials')
                                end
                            end
                            indSample = randsample(ind,n);
                            matchTrialsInd = cat(2,matchTrialsInd,indSample);
                        else
                            indSample = randsample(ind,minBinN);
                            matchTrialsInd = cat(2,matchTrialsInd,indSample);
                        end
                        respStimSort{istim} = respAllCells(indSample,:);
                        trOutStimSort{istim} = trOut(indSample);
                    end
                end
                nMatchedTrials = cumsum(cellfun(@length,trOutStimSort));
                for istim = 1:nStimBins
                    if istim == 1
                        stimSortInd = cell(1,nStimBins);
                        stimSortInd{istim} = 1:nMatchedTrials;
                    else
                        stimSortInd{istim} = ...
                            (nMatchedTrials(istim-1)+1):nMatchedTrials(istim);
                    end
                end
                trSampleIndOther{iav} = matchTrialsInd;
                decodeAnalysis(iexp).av(iav).respAllCells = respAllCells(matchTrialsInd,:);
                decodeAnalysis(iexp).av(iav).trOut = trOut(matchTrialsInd);
                
                respStimSort_av{iav} = respStimSort;
                trOutStimSort_av{iav} = trOutStimSort;
                stimSortInd_av{iav} = stimSortInd;
                matchTrialsInd_av{iav} = matchTrialsInd;
            end
                        
            nTrials = min(cellfun(@length,matchTrialsInd_av));
            if sum(cellInd) < round(nTrials.*fracPCsUsed)
                nPCs = sum(cellInd);
            elseif round(nTrials.*fracPCsUsed) > maxPCs
                nPCs = maxPCs;
            else
                nPCs = round(nTrials.*fracPCsUsed);
            end
            decodeAnalysis(iexp).nPCs = nPCs;
            decodeAnalysis(iexp).minTrials = nTrials;
            
            for iav = 1:2
                trOut = decodeDataExpt(iexp).av(iav).outcome;
                respStimSort = respStimSort_av{iav};
                trOutStimSort = trOutStimSort_av{iav};
                stimSortInd = stimSortInd_av{iav};
                matchTrialsInd = matchTrialsInd_av{iav};

                if iav == 1
                    respAllCells = respAllCells_vis;
                elseif iav == 2
                    respAllCells = respAllCells_aud;
                end
                
                emptyInd = cellfun(@isempty,respStimSort);
                respStimSort(~emptyInd) = cellfun(@(x) x(:,1:nPCs),...
                    respStimSort(~emptyInd),'unif',0);
                
                resp = respAllCells(matchTrialsInd,1:nPCs);
%                 resp = respAllCells(matchTrialsInd,cellInd);
                [~, targetTrInd] = getStimAndBehaviorYs(trOut(matchTrialsInd));

                targetCorr = corr(targetTrInd,resp);

                C = eye(size(resp,2));
                p=1;
                [~,~,targetGLM] = glmfit(resp*C,targetTrInd,'binomial');
                                
                targetWeight = targetGLM.beta(2:end);

                dv_target = mean(targetTrInd);
                
                fprintf('Expt %s, hold-out analysis\n',num2str(iexp))

                pctCorrectTarget_train = getPctCorr_trainData(targetGLM,resp,targetTrInd,dv_target);
                pctCorrectTarget_ho = getPctCorr_hoData(resp,targetTrInd,dv_target);
 
                pctCorrTarget_xStim_train = nan(1,nStimBins);
                pctCorrTarget_xStim_ho = nan(1,nStimBins);
                for istim = 1:nStimBins
                    if nStimPerBin(istim) >= minTrN_mdl
                        if isempty(trOutStimSort{istim})
                            continue
                        end
                        [~, targetStimInd] = getStimAndBehaviorYs(...
                            trOutStimSort{istim});
                        pctCorrTarget_xStim_train(istim) = getPctCorr_trainData(...
                            targetGLM,respStimSort{istim},targetStimInd,dv_target);
                        pctCorrTarget_xStim_ho(istim) = getPctCorr_hoData_subGroup(...
                            resp,targetTrInd,stimSortInd{istim},dv_target);
                    end
                end

                respOther{iav} = resp;
                trOutOther{iav} = trOut(matchTrialsInd);
                targetGLMOther{iav} = targetGLM;        
                
                targetWeight_cells = coeffAllCells(:,1:nPCs)*targetWeight;

                decodeAnalysis(iexp).av(iav).dvTarget = dv_target;
                decodeAnalysis(iexp).av(iav).correlationTarget = targetCorr;
                decodeAnalysis(iexp).av(iav).weightTarget = targetWeight_cells;
                decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_train = pctCorrectTarget_train;
                decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_holdout = pctCorrectTarget_ho;
                decodeAnalysis(iexp).av(iav).pctCorrectXStimTarget_train = pctCorrTarget_xStim_train;
                decodeAnalysis(iexp).av(iav).pctCorrectXStimTarget_holdout = pctCorrTarget_xStim_ho;
                
            end
            fprintf('Expt %s: Testing opposite model...\n',num2str(iexp))
            for iav = 1:2
                if iav == 1
                    otherAV = 2;
                    if isfield(decodeDataExpt(iexp).av(iav),'catchOutcome')
                        [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{otherAV});
                        [~,catchTargetInd] = getStimAndBehaviorYs(...
                            decodeDataExpt(iexp).av(iav).catchOutcome);
%                         catchTargetResp = catchResp(:,cellInd);
                        catchTargetResp = catchResp(:,1:nPCs);
                        nCatch = length(catchTargetInd);
                        catchDistRespAll = respOther{otherAV}(targetTrInd==0,1:nPCs);
                        targetTrInd_distonly = targetTrInd(targetTrInd==0);
                        catchDistInd = randsample(sum(targetTrInd==0),nCatch);
                        catchRespBalanced = cat(1,catchDistRespAll(catchDistInd,:),...
                            catchTargetResp);
                        
                        catchTargetIndBalanced = cat(1,targetTrInd_distonly(catchDistInd),...
                            catchTargetInd);

                        catchPctCorrectTarget = getPctCorr_trainData(...
                            targetGLMOther{iav},catchRespBalanced,catchTargetIndBalanced,...
                            decodeAnalysis(iexp).av(iav).dvTarget);        
                        catchPctCorrectTarget_audModel = getPctCorr_trainData(...
                            targetGLMOther{auditoryTrials},catchRespBalanced,...
                            catchTargetIndBalanced,...
                            decodeAnalysis(iexp).av(auditoryTrials).dvTarget);
                        

                        trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
                        matchedTrialsID = trStimID(trSampleIndOther{iav});
                        catchStimID = cat(2,ones(1,nCatch),...
                            discretize(decodeDataExpt(iexp).av(iav).catchStim,oriBins));
                        
                        [~, targetTrInd] = getStimAndBehaviorYs(trOutOther{iav});
                        catchMatchTrials = cell2mat(getMatchedValidTrialIndex(...
                            matchedTrialsID,catchStimID)); 
                        validCatchMatchTarget = getPctCorr_hoData_subGroup(...
                            respOther{iav},targetTrInd,catchMatchTrials,decodeAnalysis(iexp).av(iav).dvTarget);
                        
                        trStimID = discretize(decodeDataExpt(iexp).av(auditoryTrials).stim,ampBins);
                        matchedTrialsID = trStimID(trSampleIndOther{auditoryTrials});
                        matchedTrialsID(matchedTrialsID > 1) = 2;
                        catchStimID(catchStimID > 1) = 2;
                        [~, targetTrInd] = getStimAndBehaviorYs(trOutOther{auditoryTrials});
                        catchMatchTrials = cell2mat(getMatchedValidTrialIndex(...
                            matchedTrialsID,catchStimID));
                        validCatchMatchTarget_aud = getPctCorr_hoData_subGroup(...
                            respOther{auditoryTrials},targetTrInd,catchMatchTrials,...
                            decodeAnalysis(iexp).av(auditoryTrials).dvTarget);  
                       
                    else
                        validCatchMatchTarget = [];
                        validCatchMatchTarget_aud = [];
                        catchPctCorrectTarget = [];
                        catchPctCorrectTarget_audModel = [];
                        
                    end
                elseif iav == 2
                    otherAV = 1;
                    validCatchMatchTarget = [];
                    validCatchMatchTarget_aud = [];
                    catchPctCorrectTarget = [];
                    catchPctCorrectTarget_audModel = [];
                end
                resp = respOther{otherAV};
                [~, targetTrInd] = getStimAndBehaviorYs(trOutOther{otherAV});

                dv_target = mean(targetTrInd);

                pctCorrectTarget = getPctCorr_trainData(targetGLMOther{iav},resp,targetTrInd,dv_target);
                
                stimIndTarget = targetTrInd == 0;
               
                pctCorrectTargetxStim = nan(1,2);
                pctCorrectTargetxStim(1) = getPctCorr_trainData(...
                    targetGLMOther{iav},resp(stimIndTarget,:),...
                    targetTrInd(stimIndTarget),dv_target);
                pctCorrectTargetxStim(2) = getPctCorr_trainData(...
                    targetGLMOther{iav},resp(~stimIndTarget,:),...
                    targetTrInd(~stimIndTarget),dv_target);
                
                
                
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_otherAV = pctCorrectTarget;
                decodeAnalysis(iexp).av(iav).pctCorrectTargetxStim_otherAV = pctCorrectTargetxStim;
                decodeAnalysis(iexp).av(iav).validMatchedPctCorrectTarget_holdout = validCatchMatchTarget;
                decodeAnalysis(iexp).av(iav).invalidPctCorrectTarget = catchPctCorrectTarget;
                decodeAnalysis(iexp).av(iav).invalidPctCorrectTarget_testAudModel = catchPctCorrectTarget_audModel;
                decodeAnalysis(iexp).av(iav).validMatchedPctCorrectTarget_testAudModel = validCatchMatchTarget_aud;
                
            end
            rng(0)
            randInd = cellfun(@(x) randsample(length(x),floor(length(x)/2)),...
                trOutOther,'unif',0);
            otherRandInd = cellfun(@(x,y) setdiff(1:length(x),y),trOutOther,randInd,'unif',0);
            nt = cellfun(@length,randInd);
            nt_other = cellfun(@length,otherRandInd);
%             nt = cellfun(@length,trOutOther);
            respAV = cat(1,respOther{visualTrials}(randInd{visualTrials},:),...
                respOther{auditoryTrials}(randInd{auditoryTrials},:));
%             respAV = cat(1,respOther{visualTrials},...
%                 respOther{auditoryTrials});
            [~, targetTrInd] = getStimAndBehaviorYs(cat(2,...
                trOutOther{visualTrials}(randInd{visualTrials}),trOutOther{auditoryTrials}(randInd{auditoryTrials})));
            [~, targetTrInd_other] = getStimAndBehaviorYs(cat(2,...
                trOutOther{visualTrials}(otherRandInd{visualTrials}),trOutOther{auditoryTrials}(otherRandInd{auditoryTrials})));
            dv_target = mean(targetTrInd);
            C = eye(size(respAV,2));
            p=1;
            [~,~,targetGLM] = glmfit(respAV*C,targetTrInd,'binomial');
            
            pctCorrectTarget_combo_vis = cmbPctCorrect(nt(visualTrials),...
                getPctCorr_hoData_subGroup(respAV,targetTrInd,1:nt(visualTrials),dv_target),...
                nt_other(visualTrials),...
                getPctCorr_trainData(targetGLM,respOther{visualTrials}(otherRandInd{visualTrials},:),...
                    targetTrInd_other(1:nt_other),dv_target));
           
            pctCorrectTarget_combo_aud = cmbPctCorrect(nt(auditoryTrials),...
                getPctCorr_hoData_subGroup(respAV,targetTrInd,(nt(visualTrials)+1):sum(nt),dv_target),...
                nt_other(auditoryTrials),...
                getPctCorr_trainData(targetGLM,respOther{auditoryTrials}(otherRandInd{auditoryTrials},:),...
                    targetTrInd_other((nt(visualTrials)+1):sum(nt)),dv_target));
                
            
            [~, targetTrInd] = getStimAndBehaviorYs(trOutOther{visualTrials});
            
            pctCorrectTarget_vis = getPctCorr_hoData_subGroup(...
                respOther{visualTrials},...
                targetTrInd,randInd{visualTrials},mean(targetTrInd));
            [~, targetTrInd] = getStimAndBehaviorYs(trOutOther{auditoryTrials});
            
            pctCorrectTarget_aud = getPctCorr_hoData_subGroup(...
                respOther{auditoryTrials},...
                targetTrInd,randInd{auditoryTrials},mean(targetTrInd));
            
            
            decodeAnalysis(iexp).comboTrainWeightTarget = coeffAllCells(:,1:nPCs)*targetGLM.beta(2:end);
            decodeAnalysis(iexp).av(visualTrials).pctCorrectTarget_comboTrain = pctCorrectTarget_combo_vis;
            decodeAnalysis(iexp).av(auditoryTrials).pctCorrectTarget_comboTrain = pctCorrectTarget_combo_aud;
            decodeAnalysis(iexp).av(visualTrials).pctCorrectTarget_visTrainComboMatch = pctCorrectTarget_vis;
            decodeAnalysis(iexp).av(auditoryTrials).pctCorrectTarget_audTrainComboMatch = pctCorrectTarget_aud;
        end

    nPCsSelected = cell2mat({decodeAnalysis.nPCs});    
    fprintf('%s +/- %s (%s-%s) pcs selected\n',num2str(mean(nPCsSelected)),...
        num2str(ste(nPCsSelected,2)),num2str(min(nPCsSelected)),num2str(max(nPCsSelected)))
    end
end

    antiAnalysis = struct;
    antiAnalysis.longTC = cell(1,3);
    antiAnalysis.longTCErr = cell(1,3);
    antiAnalysis.lateCycTC = cell(1,3);
    antiAnalysis.lateCycTCErr = cell(1,3);
    antiAnalysis.lateCycTC_shuff = cell(1,3);
    antiAnalysis.firstCycTC = cell(1,3);
    antiAnalysis.firstCycTCErr = cell(1,3);
    antiAnalysis.lateCycSI = [];
    antiAnalysis.lateCycSI_shuff = [];    
    antiAnalysis.lateCycAVauROC = []; 
    antiAnalysis.lateCycAVauROC_test = [];
    antiAnalysis.lateWinSI = [];
    antiAnalysis.lateCycAV95CITest = [];
    antiAnalysis.lateWinAV95CITest = [];
    antiAnalysis.lateCycAVShuffTest = [];
    antiAnalysis.lateWinAVShuffTest = [];
    antiAnalysis.lateCycFracSI = [];
    antiAnalysis.eaCycTC = cell(2,nCycles);
    antiAnalysis.eaCycTest = cell(1,nCycles);
    targetAnalysis = struct;
    targetAnalysis.tc = cell(3,2);
    targetAnalysis.targets = cell(1,2);
    targetAnalysis.targetsAuROC = [];
    targetAnalysis.targetsAuROCTest = [];
    targetAnalysis.trOutTC = cell(length(trOutType),2);
    firstStimAuROC = [];
    lateStimAuROC = [];
    firstStimAuROCTest = [];
    lateStimAuROCTest = [];
    firstStimAuROC_aud = [];
    lateStimAuROC_aud = [];
    firstStimAuROCTest_aud = [];
    lateStimAuROCTest_aud = [];
    taskTuningTest = [];
    taskTuningPref = [];
    nModCellsPerExpt = nan(1,nexp);
    siPerExpt = cell(1,nexp);
    siPerExpt_latewin = cell(6,nexp);
    antiAnalysis.eaCycSI = cell(1,nCycles);
    antiAnalysis.eaCycAVauroc = cell(1,nCycles);
    antiAnalysis.eaCycResp = cell(2,nCycles);
    antiAnalysis.earlyCycTC = cell(3,1);
    antiAnalysis.earlyCycTCErr = cell(3,1);
    antiAnalysis.earlyTest = [];
    antiAnalysis.gm = cell(3,1);
    antiAnalysis.rsc = cell(2,1);
    for iexp = 1:nexp
        fprintf('Anti Analysis Expt %s\n',num2str(iexp))
        longTC_vis = antiDataExpt(iexp).longTC{visualTrials};
        longTC_aud = antiDataExpt(iexp).longTC{auditoryTrials};

        nCells = size(longTC_vis,2);

        firstCycResp_vis = antiDataExpt(iexp).cycTC{visualTrials,1};
        firstCycResp_aud = antiDataExpt(iexp).cycTC{auditoryTrials,1};

        cycTC_vis = antiDataExpt(iexp).cycTC(visualTrials,:);
        cycTC_aud = antiDataExpt(iexp).cycTC(auditoryTrials,:);

        lateCycTC_vis = [];
        lateCycTC_aud = [];
        lateCycles = lateCyclesMin:size(antiDataExpt(iexp).cycTC,2);
        for icyc = 1:length(lateCycles)
            lateCycTC_vis = cat(3,lateCycTC_vis,...
                cycTC_vis{lateCycles(icyc)} - mean(cycTC_vis{lateCycles(icyc)}(basewin_0,:,:),1));
            lateCycTC_aud = cat(3,lateCycTC_aud,...
                cycTC_aud{lateCycles(icyc)} - mean(cycTC_aud{lateCycles(icyc)}(basewin_0,:,:),1));
        end
        
        earlyCycTC_vis = [];
        earlyCycTC_aud = [];
        for icyc = 1:2
            earlyCycTC_vis = cat(3,earlyCycTC_vis,...
                cycTC_vis{icyc} - mean(cycTC_vis{icyc}(basewin_0,:,:),1));
            earlyCycTC_aud = cat(3,earlyCycTC_aud,...
                cycTC_aud{icyc} - mean(cycTC_aud{icyc}(basewin_0,:,:),1));
        end
        earlyTCAV = cat(3,earlyCycTC_vis,earlyCycTC_aud);
        earlyTCTest = ttest(squeeze(mean(earlyTCAV(respwin,:,:))),...
            squeeze(mean(earlyTCAV(basewin_0,:,:))),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha)';
        
        cycAV = cell(1,length(cycTC_vis));
        for icyc = 1:length(cycTC_vis)
            cycAV{icyc} = cat(3,cycTC_vis{icyc},cycTC_aud{icyc});
        end
        
        lastCycResp = cell(2,1);
        lastCycTC = cell(2,1);
        lastCycAV = [];
        for icyc = nCycles:length(cycTC_vis)
            for iav = 1:2
                lastCycResp{iav} = cat(2,lastCycResp{iav}, cell2mat(cellfun(@(x) ...
                    squeeze(mean(x(respwin,:,:))) - squeeze(mean(x(basewin_0,:,:))),...
                    antiDataExpt(iexp).cycTC(iav,icyc),'unif',0)));
                lastCycTC{iav} = cat(3,lastCycTC{iav}, cell2mat(cellfun(@(x) ...
                    x - mean(x(basewin_0,:,:),1),...
                    antiDataExpt(iexp).cycTC(iav,icyc),'unif',0)));
            end
            lastCycAV = cat(3,lastCycAV,cycAV{icyc});
        end
        cycResp = cat(2,cellfun(@(x) ...
            squeeze(mean(x(respwin,:,:))) - squeeze(mean(x(basewin_0,:,:))),...
            antiDataExpt(iexp).cycTC(:,1:nCycles-1),'unif',0),lastCycResp);
        cycBL = cellfun(@(x) x-mean(x(basewin_0,:,:)),antiDataExpt(iexp).cycTC,'unif',0);
        cycTest = cellfun(@(x) ttest(squeeze(mean(x(respwin,:,:))),...
            squeeze(mean(x(basewin_0,:,:))),'dim',2,'tail','right','alpha',0.05/(nCycles-1))',...
            [cycAV(:,1:nCycles-1) {lastCycAV}],'unif',0);
        
        avAuroc = cell(1,nCycles);
        for icyc = 1:nCycles
                avAuroc{icyc} = nan(1,nCells);
            for icell = 1:nCells
                avAuroc{icyc}(icell) = roc_gh(...
                    cycResp{auditoryTrials,icyc}(icell,:),...
                    cycResp{visualTrials,icyc}(icell,:));
            end
        end        
%         lateCycAVTest = ttest2(squeeze(mean(lateCycTC_vis(respwin,:,:),1))',...
%             squeeze(mean(lateCycTC_aud(respwin,:,:),1))','dim',1,'alpha',0.05);
%         lateWinAVTest = ttest2(squeeze(mean(longTC_vis(lateWinFr,:,:),1))',...
%             squeeze(mean(longTC_aud(lateWinFr,:,:),1))','dim',1,'alpha',0.05);
        rng(0)
        [lateCycSI95CITest,lateCycAVShuffTest] = testSelectivityIndex(...
            squeeze(mean(lateCycTC_vis(respwin,:,:),1))',...
            squeeze(mean(lateCycTC_aud(respwin,:,:),1))');
        [lateWinSI95CITest,lateWinAVShuffTest] = testSelectivityIndex(...
            squeeze(mean(longTC_vis(lateWinFr,:,:),1))',...
            squeeze(mean(longTC_aud(lateWinFr,:,:),1))');
        nModCellsPerExpt(iexp) = sum(lateCycSI95CITest & respCellsExpt(iexp).lateCycRespCells);
        
        rng(0)
        lateCycSI = getSelectivityIndex(squeeze(mean(lateCycTC_vis(respwin,:,:),1))',...
            squeeze(mean(lateCycTC_aud(respwin,:,:),1))');
        lateWinSI = getSelectivityIndex(squeeze(mean(longTC_vis(lateWinFr,:,:),1))',...
            squeeze(mean(longTC_aud(lateWinFr,:,:),1))');
        [lateCycAVauROC, lateCycAVauROC_test] = getAVauROC(squeeze(mean(lateCycTC_aud(respwin,:,:),1))',...
            squeeze(mean(lateCycTC_vis(respwin,:,:),1))');
        if size(lateCycTC_vis,3) ~= size(lateCycTC_aud,3)
            if size(lateCycTC_vis,3) > size(lateCycTC_aud,3)
                randInd = randsample(size(lateCycTC_vis,3),size(lateCycTC_aud,3));
                lateCycRespAll = cat(1,squeeze(mean(lateCycTC_vis(respwin,:,randInd),1))',...
                    squeeze(mean(lateCycTC_aud(respwin,:,:),1))');
            else
                randInd = randsample(size(lateCycTC_aud,3),size(lateCycTC_vis,3));
                lateCycRespAll = cat(1,squeeze(mean(lateCycTC_vis(respwin,:,:),1))',...
                    squeeze(mean(lateCycTC_aud(respwin,:,randInd),1))');
            end
        else
            lateCycRespAll = cat(1,squeeze(mean(lateCycTC_vis(respwin,:,:),1))',...
                squeeze(mean(lateCycTC_aud(respwin,:,randInd),1))');                
        end
                
        [shuffResp_vis, shuffResp_aud] = twoRandGroups(lateCycRespAll);
        lateCycSI_shuff = getSelectivityIndex(shuffResp_vis,shuffResp_aud);
        v=mean(mean(lateCycTC_vis(respwin,:,:),3),1);
        v(v<0) = 0;
        a=mean(mean(lateCycTC_aud(respwin,:,:),3),1);
        a(a<0) = 0;
        lateCycFracSI = (v-a)./v;
        lateCycFracSI_aud = (a-v)./a;
        lateCycFracSI(a > v) = -lateCycFracSI_aud(a > v);

        eaCycSI = cellfun(@(x,y) getSelectivityIndex(x',y'),...
            cycResp(visualTrials,:),cycResp(auditoryTrials,:),'unif',0);
        
        firstCycTC_vis = antiDataExpt(iexp).cycTC{visualTrials,1};
        firstCycTC_aud = antiDataExpt(iexp).cycTC{auditoryTrials,1};

        firstVisResp = squeeze(mean(firstCycTC_vis(respwin,:,:),1));
        lateVisResp = squeeze(mean(lateCycTC_vis(respwin,:,:),1));

        tarVisResp = cellfun(@(x) ...
            squeeze(mean(x(respwin_target,:,:),1) - mean(x(basewin_0_target,:,:),1)),...
            cat(2,targetDataExpt(iexp).av(visualTrials).tc,...
            {cat(3,targetDataExpt(iexp).av(visualTrials).tc{1},...
            targetDataExpt(iexp).av(visualTrials).tc{2})}),'unif',0);
        auroc_first = nan(nCells,3);
        auroc_first_test = nan(nCells,3);
        auroc_late = nan(nCells,3);
        auroc_late_test = nan(nCells,3);
        auroc_target = nan(nCells,1);
        auroc_target_test = nan(nCells,1);
        for icell = 1:nCells
            auroc_target(icell) = roc_gh(tarVisResp{1}(icell,:),tarVisResp{2}(icell,:));
            auroc_target_test(icell) = ranksum(...
                tarVisResp{1}(icell,:),tarVisResp{2}(icell,:)) < 0.05;

            auroc_first(icell,:) = cellfun(@(x) ...
                roc_gh(firstVisResp(icell,:),x(icell,:)),tarVisResp);
            auroc_first_test(icell,:) = cellfun(@(x) ...
                ranksum(firstVisResp(icell,:),x(icell,:)),tarVisResp);
            auroc_late(icell,:) = cellfun(@(x) ...
                roc_gh(lateVisResp(icell,:),x(icell,:)),tarVisResp);
            auroc_late_test(icell,:) = cellfun(@(x) ...
                ranksum(lateVisResp(icell,:),x(icell,:)),tarVisResp);
        end

        firstAudResp = squeeze(mean(firstCycTC_aud(respwin,:,:),1));
        lateAudResp = squeeze(mean(lateCycTC_aud(respwin,:,:),1));
        tarAudResp = cell2mat(cellfun(@(x) ...
            squeeze(mean(x(respwin_target,:,:),1) - mean(x(basewin_0_target,:,:),1)),...
            {cat(3,targetDataExpt(iexp).av(auditoryTrials).tc{1},...
            targetDataExpt(iexp).av(auditoryTrials).tc{2})},'unif',0));

        auroc_first_aud = nan(nCells,1);
        auroc_first_aud_test = nan(nCells,1);
        auroc_late_aud = nan(nCells,1);
        auroc_late_aud_test = nan(nCells,1);
        for icell = 1:nCells        
            auroc_first_aud(icell,:) =  ...
                roc_gh(firstAudResp(icell,:),tarAudResp(icell,:));
            auroc_first_aud_test(icell,:) =  ...
                ranksum(firstAudResp(icell,:),tarAudResp(icell,:));
            auroc_late_aud(icell,:) = ...
                roc_gh(firstAudResp(icell,:),tarAudResp(icell,:));
            auroc_late_aud_test(icell,:) =  ...
                ranksum(firstAudResp(icell,:),tarAudResp(icell,:));
        end

        stimResp4anova = cat(2,firstVisResp,tarVisResp{1},tarVisResp{2});
        stimID4anova = cat(2,ones(1,size(firstVisResp,2)),...
            ones(1,size(tarVisResp{1},2))*2,ones(1,size(tarVisResp{2},2))*3);
        taskStimAnova = nan(1,nCells);
        for icell = 1:nCells
            taskStimAnova(icell) = anova1(stimResp4anova(icell,:),stimID4anova,'off');
        end
        [~,taskTuningID] = max(cat(2,mean(firstVisResp,2),...
            mean(tarVisResp{1},2),mean(tarVisResp{2},2)),[],2);

        antiAnalysis.longTC{visualTrials} = cat(2,antiAnalysis.longTC{visualTrials},...
            mean(longTC_vis,3));
        antiAnalysis.longTC{auditoryTrials} = cat(2,antiAnalysis.longTC{auditoryTrials},...
            mean(longTC_aud,3));
        antiAnalysis.longTC{allTrialsInd} = cat(2,antiAnalysis.longTC{allTrialsInd},...
            mean(cat(3,longTC_vis,longTC_aud),3));
        antiAnalysis.longTCErr{allTrialsInd} = cat(2,antiAnalysis.longTCErr{allTrialsInd},...
            ste(cat(3,longTC_vis,longTC_aud),3));
        antiAnalysis.lateCycTC{visualTrials} = cat(2,antiAnalysis.lateCycTC{visualTrials},...
            mean(lateCycTC_vis,3));
        antiAnalysis.lateCycTC{auditoryTrials} = cat(2,antiAnalysis.lateCycTC{auditoryTrials},...
            mean(lateCycTC_aud,3));
        antiAnalysis.lateCycTC{allTrialsInd} = cat(2,antiAnalysis.lateCycTC{allTrialsInd},...
            mean(cat(3,lateCycTC_vis,lateCycTC_aud),3));
        antiAnalysis.lateCycTCErr{visualTrials} = cat(2,antiAnalysis.lateCycTCErr{visualTrials},...
            ste(lateCycTC_vis,3));
        antiAnalysis.lateCycTCErr{auditoryTrials} = cat(2,antiAnalysis.lateCycTCErr{auditoryTrials},...
            ste(lateCycTC_aud,3));
        antiAnalysis.lateCycTCErr{allTrialsInd} = cat(2,antiAnalysis.lateCycTCErr{allTrialsInd},...
            ste(cat(3,lateCycTC_vis,lateCycTC_aud),3));
        antiAnalysis.lateCycSI = cat(2,antiAnalysis.lateCycSI,lateCycSI);
        antiAnalysis.lateCycSI_shuff = cat(2,antiAnalysis.lateCycSI_shuff,lateCycSI_shuff);
        antiAnalysis.lateCycAVauROC = cat(2,antiAnalysis.lateCycAVauROC,lateCycAVauROC);
        antiAnalysis.lateCycAVauROC_test = cat(2,antiAnalysis.lateCycAVauROC_test,lateCycAVauROC_test);
        antiAnalysis.lateWinSI = cat(2,antiAnalysis.lateWinSI,lateWinSI);
        antiAnalysis.lateCycAV95CITest = cat(2,antiAnalysis.lateCycAV95CITest,lateCycSI95CITest');
        antiAnalysis.lateWinAV95CITest = cat(2,antiAnalysis.lateWinAV95CITest,lateWinSI95CITest');
        antiAnalysis.lateCycAVShuffTest = cat(2,antiAnalysis.lateCycAVShuffTest,lateCycAVShuffTest');
        antiAnalysis.lateWinAVShuffTest = cat(2,antiAnalysis.lateWinAVShuffTest,lateWinAVShuffTest');
        antiAnalysis.lateCycFracSI = cat(2,antiAnalysis.lateCycFracSI,lateCycFracSI);
        antiAnalysis.firstCycTC{visualTrials} = cat(2,antiAnalysis.firstCycTC{visualTrials},...
            mean(firstCycTC_vis,3));
        antiAnalysis.firstCycTC{auditoryTrials} = cat(2,antiAnalysis.firstCycTC{auditoryTrials},...
            mean(firstCycTC_aud,3));
        antiAnalysis.firstCycTC{allTrialsInd} = cat(2,antiAnalysis.firstCycTC{allTrialsInd},...
            mean(cat(3,firstCycTC_vis,firstCycTC_aud),3));
        antiAnalysis.firstCycTCErr{visualTrials} = cat(2,antiAnalysis.firstCycTCErr{visualTrials},...
            ste(firstCycTC_vis,3));
        antiAnalysis.firstCycTCErr{auditoryTrials} = cat(2,antiAnalysis.firstCycTCErr{auditoryTrials},...
            ste(firstCycTC_aud,3));
        antiAnalysis.firstCycTCErr{allTrialsInd} = cat(2,antiAnalysis.firstCycTCErr{allTrialsInd},...
            ste(cat(3,firstCycTC_vis,firstCycTC_aud),3));
        antiAnalysis.eaCycSI = cellfun(@(x,y) cat(2,x,y),antiAnalysis.eaCycSI,eaCycSI,'unif',0);
        antiAnalysis.eaCycResp = cellfun(@(x,y) cat(1,x,y),antiAnalysis.eaCycResp,...
            cellfun(@(x)mean(x,2),cycResp,'unif',0),'unif',0);
        antiAnalysis.eaCycAVauroc = cellfun(@(x,y) cat(2,x,y),antiAnalysis.eaCycAVauroc,avAuroc,'unif',0);
        antiAnalysis.eaCycTC = cellfun(@(x,y) cat(2,x,mean(y,3)),antiAnalysis.eaCycTC,...
            [antiDataExpt(iexp).cycTC(:,1:nCycles-1) lastCycTC],'unif',0);
        antiAnalysis.eaCycTest = cellfun(@(x,y) cat(2,x,y), antiAnalysis.eaCycTest,cycTest,'unif',0);
        antiAnalysis.earlyCycTC = cellfun(@(x,y) cat(2,x,mean(y,3)),antiAnalysis.earlyCycTC,...
            {earlyCycTC_vis;earlyCycTC_aud;earlyTCAV},'unif',0);
        antiAnalysis.earlyCycTCErr = cellfun(@(x,y) cat(2,x,ste(y,3)),antiAnalysis.earlyCycTCErr,...
            {earlyCycTC_vis;earlyCycTC_aud;earlyTCAV},'unif',0);
        antiAnalysis.earlyTest = cat(2,antiAnalysis.earlyTest,earlyTCTest);
        %for noise correlations, build a matrix that is all cells x all
        %cells, with nans in place for cell pairs from different
        %experiments.
        if iexp==1
            antiAnalysis.rsc = respCellsExpt(iexp).rsc;
            antiAnalysis.gm = respCellsExpt(iexp).gm;
        else
            antiAnalysis.rsc = cellfun(@(x,y) ...
                cat(2,cat(1,x,nan(size(y,1),size(x,2))),cat(1,nan(size(x,1),size(y,2)),y)),...
                antiAnalysis.rsc,respCellsExpt(iexp).rsc,'unif',0);
            antiAnalysis.gm = cellfun(@(x,y) ...
                cat(2,cat(1,x,nan(size(y,1),size(x,2))),cat(1,nan(size(x,1),size(y,2)),y)),...
                antiAnalysis.gm,respCellsExpt(iexp).gm,'unif',0);
        end
        
        targetAnalysis.tc{1,visualTrials} = cat(2,targetAnalysis.tc{1,visualTrials},...
            mean(targetDataExpt(iexp).av(visualTrials).tc{1},3));
        targetAnalysis.tc{2,visualTrials} = cat(2,targetAnalysis.tc{2,visualTrials},...
            mean(targetDataExpt(iexp).av(visualTrials).tc{2},3));
        targetAnalysis.tc{3,visualTrials} = cat(2,targetAnalysis.tc{3,visualTrials},...
            mean(cat(3,targetDataExpt(iexp).av(visualTrials).tc{1},...
            targetDataExpt(iexp).av(visualTrials).tc{2}),3));
        targetAnalysis.tc{1,auditoryTrials} = cat(2,targetAnalysis.tc{1,auditoryTrials},...
            mean(targetDataExpt(iexp).av(auditoryTrials).tc{1},3));
        targetAnalysis.tc{2,auditoryTrials} = cat(2,targetAnalysis.tc{2,auditoryTrials},...
            mean(targetDataExpt(iexp).av(auditoryTrials).tc{2},3));
        targetAnalysis.tc{3,auditoryTrials} = cat(2,targetAnalysis.tc{3,auditoryTrials},...
            mean(cat(3,targetDataExpt(iexp).av(auditoryTrials).tc{1},...
            targetDataExpt(iexp).av(auditoryTrials).tc{2}),3));
        if strcmp(ds,'FSAV_attentionV1') | strcmp(ds,'FSAV_attentionV1')
            targetAnalysis.trOutTC(:,visualTrials) = cellfun(@(x,y) cat(2,x,mean(y(1:(nBaselineFr*2),:,:),3)),...
                targetAnalysis.trOutTC(:,visualTrials),...
                targetDataExpt(iexp).av(visualTrials).trOutTC','unif',0);
            targetAnalysis.trOutTC(:,auditoryTrials) = cellfun(@(x,y) cat(2,x,mean(y(1:(nBaselineFr*2),:,:),3)),...
                targetAnalysis.trOutTC(:,auditoryTrials),...
                targetDataExpt(iexp).av(auditoryTrials).trOutTC','unif',0);
        end
        targetAnalysis.targets{visualTrials} = cat(2,targetAnalysis.targets{visualTrials},...
            cat(1,ones(1,nCells).*mean(targetDataExpt(iexp).av(visualTrials).targets{1}),...
            ones(1,nCells).*mean(targetDataExpt(iexp).av(visualTrials).targets{2})));
        targetAnalysis.targetsAuROC = cat(1,targetAnalysis.targetsAuROC,auroc_target);
        targetAnalysis.targetsAuROCTest = cat(1,targetAnalysis.targetsAuROCTest,auroc_target_test);
        taskTuningTest = cat(1,taskTuningTest,taskStimAnova' < 0.05);
        taskTuningPref = cat(1,taskTuningPref,taskTuningID);

        firstStimAuROC = cat(1,firstStimAuROC,auroc_first);
        lateStimAuROC = cat(1,lateStimAuROC,auroc_late);
        firstStimAuROCTest = cat(1,firstStimAuROCTest,auroc_first_test);
        lateStimAuROCTest = cat(1,lateStimAuROCTest,auroc_late_test);

        firstStimAuROC_aud = cat(1,firstStimAuROC_aud,auroc_first_aud);
        lateStimAuROC_aud = cat(1,lateStimAuROC_aud,auroc_late_aud);
        firstStimAuROCTest_aud = cat(1,firstStimAuROCTest_aud,auroc_first_aud_test);
        lateStimAuROCTest_aud = cat(1,lateStimAuROCTest_aud,auroc_late_aud_test);

        siPerExpt{iexp} = lateCycSI;
        siPerExpt_latewin{iexp} = lateWinSI;
    %     ind = taskStimAnova' > 0.05 & ttest(firstVisResp,tarVisResp{3},'dim',2);
    %     taskTuningID_temp = taskTuningID;
    %     taskTuningID_temp(ind) = 4;
    %     taskTuningPref_temp =cat(1,taskTuningPref_temp,taskTuningID_temp);
    end
    % antiAnalysis.adapt = cellfun(@(x,y) ...
    %     mean(x(respwin,:),1)./mean(y(respwin,:),1),...
    %     antiAnalysis.lateCycTC,antiAnalysis.firstCycTC,'unif',0);
    
    pctModCellsExpt = nModCellsPerExpt./cellfun(@sum,{respCellsExpt.lateCycRespCells});
    fprintf('%s+/-%s mod cells, rew mice\n',num2str(mean(pctModCellsExpt(logical(rewExptInd)))),...
        num2str(ste(pctModCellsExpt(logical(rewExptInd)),2)))
    fprintf('%s+/-%s mod cells, no rew mice\n',num2str(mean(pctModCellsExpt(logical(~rewExptInd)))),...
        num2str(ste(pctModCellsExpt(logical(~rewExptInd)),2)))
    

    cellInfo = struct;
    cellInfo.firstRespCells = logical(cell2mat({respCellsExpt.firstRespCells}'));
    cellInfo.lateRespCells = logical(cell2mat({respCellsExpt.lateRespCells}'));
    cellInfo.lateSuppCells = logical(cell2mat({respCellsExpt.lateSuppCells}'));
    cellInfo.lateCycRespCells = logical(cell2mat({respCellsExpt.lateCycRespCells}'));
    cellInfo.minRespCells = (mean(antiAnalysis.firstCycTC{visualTrials}(respwin,:),1) > ...
        minRespThreshold & mean(antiAnalysis.firstCycTC{auditoryTrials}(respwin,:),1) > ...
        minRespThreshold)';
    cellInfo.targetRespCells = logical(cell2mat({respCellsExpt.targetRespCells}'));
    cellInfo.isShortCycExpt = isShortCycExpt;
    cellInfo.isTuned = logical(cell2mat({oriTuningExpt.isTuned}))';
    cellInfo.tuningReliability = cell2mat({oriTuningExpt.tuningReliability})';
    cellInfo.oriResp = cell2mat({oriTuningExpt.oriResp}');
    cellInfo.oriRespErr = cell2mat({oriTuningExpt.oriRespErr}');
    cellInfo.oriFit = cell2mat({oriTuningExpt.fit})';
    cellInfo.oriPref = cell2mat({oriTuningExpt.oriPref})';
    cellInfo.hwhm = hwhmFromOriFit(cellInfo.oriFit(:,1:180)',1:180)';
    cellInfo.taskTuningPref = taskTuningPref;
    % cellInfo.taskTuningPref = taskTuningPref_temp;
    cellInfo.taskTuningTest = taskTuningTest;

    auroc_first = nan(length(cellInfo.firstRespCells),1);
    auroc_late = nan(length(cellInfo.firstRespCells),1);
    auroc_first_test = nan(length(cellInfo.firstRespCells),1);
    auroc_late_test = nan(length(cellInfo.firstRespCells),1);
    auroc_late_HT = lateStimAuROC(:,1);
    auroc_late_HT_test = lateStimAuROCTest(:,1) < 0.05;
    for icell = 1:length(cellInfo.firstRespCells)
        if ~cellInfo.targetRespCells(icell)
            auroc_first(icell) = firstStimAuROC(icell,3);
            auroc_late(icell) = lateStimAuROC(icell,3);
            auroc_first_test(icell) = firstStimAuROCTest(icell,3) < 0.05;
            auroc_late_test(icell) = lateStimAuROCTest(icell,3) < 0.05;
        elseif ~targetAnalysis.targetsAuROCTest(icell)
            auroc_first(icell) = firstStimAuROC(icell,3);
            auroc_late(icell) = lateStimAuROC(icell,3);
            auroc_first_test(icell) = firstStimAuROCTest(icell,3) < 0.05;
            auroc_late_test(icell) = lateStimAuROCTest(icell,3) < 0.05;
        elseif targetAnalysis.targetsAuROCTest(icell) && targetAnalysis.targetsAuROC(icell) < 0.5
            auroc_first(icell) = firstStimAuROC(icell,1);
            auroc_late(icell) = lateStimAuROC(icell,1);
            auroc_first_test(icell) = firstStimAuROCTest(icell,1) < 0.05;
            auroc_late_test(icell) = lateStimAuROCTest(icell,1) < 0.05;
        elseif targetAnalysis.targetsAuROCTest(icell) && targetAnalysis.targetsAuROC(icell) >= 0.5
            auroc_first(icell) = firstStimAuROC(icell,2);
            auroc_late(icell) = lateStimAuROC(icell,2);
            auroc_first_test(icell) = firstStimAuROCTest(icell,2) < 0.05;
            auroc_late_test(icell) = lateStimAuROCTest(icell,2) < 0.05;
        end
        
    end
    cellInfo.firstStimAuROC = auroc_first;
    cellInfo.lateStimAuROC = auroc_late;
    cellInfo.firstStimAuROCTest = auroc_first_test;
    cellInfo.lateStimAuROCTest = auroc_late_test;
    cellInfo.audFirstStimAuROC = firstStimAuROC_aud;
    cellInfo.audLateStimAuROC = lateStimAuROC_aud;
    cellInfo.audFirstStimAuROCTest = firstStimAuROCTest_aud < 0.05;
    cellInfo.audLateStimAuROCTest = lateStimAuROCTest_aud < 0.05;
    if strcmp(ds,'FSAV_attentionV1')|strcmp(ds,'FSAV_attentionV1_noAttn')
        
        cellInfo.targetComboWeight = [];
        cellInfo.detectComboWeight = [];
        cellInfo.targetComboWeight = [];
        cellInfo.detectComboWeight = [];
        cellInfo.av(visualTrials).targetWeight = [];
        cellInfo.av(visualTrials).detectWeight = [];
        cellInfo.av(auditoryTrials).targetWeight = [];
        cellInfo.av(auditoryTrials).detectWeight = [];        
        for iexp = 1:nexp
            comboWeight_stim = nan(1,length(decodeAnalysis(iexp).cellInd));
            comboWeight_ch = nan(1,length(decodeAnalysis(iexp).cellInd));
            comboWeight_stim(decodeAnalysis(iexp).cellInd) = decodeAnalysis(iexp).comboTrainWeightTarget;
            comboWeight_ch(decodeAnalysis(iexp).cellInd) = decodeAnalysis(iexp).comboTrainWeightDetect;
            cellInfo.targetComboWeight = cat(2,cellInfo.targetComboWeight,comboWeight_stim);
            cellInfo.detectComboWeight = cat(2,cellInfo.detectComboWeight,comboWeight_ch);
            for iav = 1:2
                w_stim = nan(1,length(decodeAnalysis(iexp).cellInd));
                w_ch = nan(1,length(decodeAnalysis(iexp).cellInd));
                w_stim(decodeAnalysis(iexp).cellInd) = decodeAnalysis(iexp).av(iav).weightTarget;
                w_ch(decodeAnalysis(iexp).cellInd) = decodeAnalysis(iexp).av(iav).weightDetect;
                cellInfo.av(iav).targetWeight = cat(2,cellInfo.av(iav).targetWeight,w_stim);
                cellInfo.av(iav).detectWeight = cat(2,cellInfo.av(iav).detectWeight,w_ch);
            end
        end
        save([fnout 'imgAnalysisData'],...
            'decodeAnalysis','antiAnalysis','targetAnalysis','cellInfo','respCellsExpt','antiDataExpt','siPerExpt','oriTuningExpt')
    elseif strcmp(ds,'FSAV_V1_audControl')
        cellInfo.targetWeight = nan(size(cellInfo.firstRespCells));
        cellInfo.dcModelCells = logical(cell2mat({decodeAnalysis.cellInd}'));
    else
%         cellInfo.targetWeight = targetWeightBoot;
%         cellInfo.dcModelCells = logical(cell2mat({decodeAnalysis.cellInd}'));
    end
    save([fnout 'imgAnalysisData'],...
        'decodeAnalysis','antiAnalysis','targetAnalysis','cellInfo','respCellsExpt','antiDataExpt','siPerExpt','oriTuningExpt')
end
%% plotting params
respTCLim = [-0.005 0.05];
cycTCLim = [-0.005 0.015];
eaCycTCLim = [-0.005 0.05];
cycTCLim_minRespCells = [-0.005 0.025];
scatLim_win = [-0.2 0.6];
scatLim_cyc = [-0.035 0.085];
hmLim = [-0.1 0.1];
exCellTCLim = [-0.02 0.15];
oriRespLim = [-0.05 0.15];
siLim = [-10 10];
siOriLim = [-4 4];
oriBarLim_win = [0 0.08];
oriBarLim_resp = [0 0.04];
oriLim_taskResp = [-0.005 0.035];
oriNLim = [0 120];
oriTCLim = [-0.005 0.08];
targetTCLim = [-0.015 0.08];
outTCLim = [-0.005 0.04];
firstTCLim = [-0.005 0.04];
% adaptLim = [0 1];
suppTCLim = [-0.05 0.005];
suppScatLim_win = [-0.2 0.1];
suppScatLim_cyc = [-0.015 0.015];
cellRespTCLim = [-0.05 0.15];
exCellRespTCLim = [-0.1 0.1];
siBinLim = [-4 4];
stimRespLim = [-0.01 0.05];
avName = {'Vis','Aud'};

tcStartFrame = 26;
cycTCEndTimeMs = 350;
cycTCEndFr = 45;
ttLabel_long = 0:500:2500;
ttLabel_cyc = -200:100:cycTCEndTimeMs;
ttLabel_target = -1000:250:900;
preTargetStimLabel = -700:350:0;
nFr_long = size(antiAnalysis.longTC{1,1},1);
tt_longTC = ((tcStartFrame:nFr_long)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
ttLabelFr_long = ((ttLabel_long./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_cyc = ((ttLabel_cyc./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_target = ((ttLabel_target./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr_target)+1);

nFr_cyc = size(antiAnalysis.lateCycTC{1,1},1);
tt_cycTC = ((tcStartFrame:nFr_cyc)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
tt_targetTC = ((1:nFr_cyc)-(nBaselineFr+nVisDelayFr_target)).*(1000/frameRateHz);

lateWinTT = ([lateWinFr(1) lateWinFr(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
respWinTT = ([respwin(1) respwin(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
respWinTT_target = (...
    [respwin_target(1) respwin_target(end)] - (nBaselineFr+nVisDelayFr_target))...
    .*(1000/frameRateHz);
baseWinTT = (...
    [basewin_0(1) basewin_0(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);

movWinLabelFr = 30:(30+nMovWin-1);
movWinLabelMs = (movWinLabelFr - (nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);

weightLim = [-3 4.2];
binnedWeightLim = [-0.4 0.4];
weightLimSum = [-0.8 0.8];
siLimSum = [-0.5 2.5];


% siBins = [-10,-1,1,10];
siBins = [-10,0,10];
nSIBins = length(siBins)-1;
SI = antiAnalysis.lateCycSI;
[n,~,siBinID] = histcounts(SI,siBins);
% siBinID = nan(1,length(SI));
% siBinID(~antiAnalysis.lateCycAV95CITest) = 2;
% siBinID(antiAnalysis.lateCycAV95CITest & SI < 0) = 1;
% siBinID(antiAnalysis.lateCycAV95CITest & SI > 0) = 3;
% nSIBins = 3;

aurocBinID = nan(1,length(SI));
aurocBinID(~cellInfo.lateStimAuROCTest) = 3;
aurocBinID(cellInfo.lateStimAuROCTest & cellInfo.lateStimAuROC < 0.5) = 1;
aurocBinID(cellInfo.lateStimAuROCTest & cellInfo.lateStimAuROC > 0.5) = 2;


aurocBinID_aud = nan(1,length(SI));
aurocBinID_aud(~cellInfo.audLateStimAuROCTest) = 3;
aurocBinID_aud(cellInfo.audLateStimAuROCTest & cellInfo.audLateStimAuROC < 0.5) = 1;
aurocBinID_aud(cellInfo.audLateStimAuROCTest & cellInfo.audLateStimAuROC > 0.5) = 2;

Dr = cellfun(@(x) mean(x(respwin,:),1),antiAnalysis.lateCycTC,'unif',0);
Tr_vis = cellfun(@(x) mean(x(respwin_target,:),1),targetAnalysis.tc(:,visualTrials),'unif',0);
Tr_aud = cellfun(@(x) mean(x(respwin_target,:),1),targetAnalysis.tc(:,auditoryTrials),'unif',0);

nMice = length(mice);

%% plot anticipation analysis (Figure 2)

setFigParams4Print('landscape')

% heatmaps
lateWinRespAll = mean(antiAnalysis.longTC{allTrialsInd}(lateWinFr,:),1);
[~,lateWinSortInd] = sort(lateWinRespAll);
hm = flipud(antiAnalysis.longTC{allTrialsInd}(:,lateWinSortInd)');

figure
colormap(brewermap([],'*RdBu'));
% subplot 121
imagesc(hm(:,tcStartFrame:end))
hold on
% if strcmp(ds,'FSAV_attentionV1')
%     exCellInd = [exampleCell_1,exampleCell_2,exampleCell_3];
%     exCellMat = zeros(1,length(cellInfo.firstRespCells));
%     exCellMat(exCellInd) = 1;
%     exCellSortInd = find(flip(exCellMat(lateWinSortInd)));
%     hline(exCellSortInd,'k-')
% end
figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabelFr_long)
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title('All Trials, All Cells')
print([fnout 'heatmapAllTrialsAllCells'],'-dpdf','-fillpage')


if strcmp(ds,'FSAV_attentionV1')
    exCellInd = false(length(cellInfo.firstRespCells),1);
    exCellInd([exampleCell_1,exampleCell_3]) = true;
    exCellInd = flipud(exCellInd(lateWinSortInd));
elseif strcmp(ds,'FSAV_attentionV1_noAttn')
    exCellInd = false(length(cellInfo.firstRespCells),1);
    exCellInd(exampleCell_2) = true;
    exCellInd = flipud(exCellInd(lateWinSortInd));
else
    exCellInd = false(length(cellInfo.firstRespCells),1);
end

figure
colormap gray
% subplot 122
ind = cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells;
indBW = ind; indBW(ind) = 0; indBW(~ind) = 1;
imagesc(flipud(indBW(lateWinSortInd)))
figYAxis([],'Cell #',[])
figAxForm
colorbar
title('Anti. Resp. Cells (First, Late, Late Win.)')
print([fnout 'cellIDforHeatmap'],'-dpdf','-fillpage')

figure
colormap gray
subplot 121
ind = antiAnalysis.earlyTest'==1 |cellInfo.lateCycRespCells==1;
stimRespCellsBW = flipud(imcomplement(ind(lateWinSortInd)));
imagesc(stimRespCellsBW)
figYAxis([],'Cell #',[])
figAxForm
title('Resp. to any distractor')
subplot 122
ind = cellInfo.lateCycRespCells==1;
lateStimRespCellsBW = flipud(imcomplement(ind(lateWinSortInd)));
imagesc(lateStimRespCellsBW)
figYAxis([],'Cell #',[])
figAxForm
title('Resp. to late distractor')
print([fnout 'cellIDforHeatmap_cycResp'],'-dpdf','-fillpage')

save([fnout '_heatmaps'],'stimRespCellsBW','lateStimRespCellsBW','hm','exCellInd')

% time-courses and quantification scatters
setFigParams4Print('portrait')
figure
subplot 311
ind = cellInfo.lateCycRespCells|antiAnalysis.earlyTest'| cellInfo.lateRespCells;
% ind = cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells;
for iav = 1:2
    y = mean(antiAnalysis.longTC{iav}...
        ((tcStartFrame:end),(ind) & ...
        cellInfo.isShortCycExpt),2);
    yerr = ste(antiAnalysis.longTC{iav}...
        ((tcStartFrame:end),(ind) & cellInfo.isShortCycExpt),2);
    hold on
    shadedErrorBar_chooseColor(tt_longTC,y,yerr,cueColor{iav});
end
figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
figYAxis([],'dF/F',respTCLim)  
vline(lateWinTT,'k--')
hline(0,'k:')
figAxForm([],0)
title(sprintf('First Stim Responsive Cells (%s/%s)',...
    num2str(sum((...
    ind)...
    & cellInfo.isShortCycExpt)),...
    num2str(sum(cellInfo.isShortCycExpt))))

subplot 323
x = mean(antiAnalysis.longTC{visualTrials}(lateWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),1);
xerr = ste(x,2);
y = mean(antiAnalysis.longTC{auditoryTrials}(lateWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),1);
yerr = ste(y,2);
plot(x,y,'.')
hold on
errorbar(mean(x),mean(y),yerr,yerr,xerr,xerr,'.')
plot(scatLim_win,scatLim_win,'k--')
plot(scatLim_win,[0 0],'k:')
plot([0 0],scatLim_win,'k:')
[~,p] = ttest(x,y);
figXAxis([],'Visual (dF/F)',scatLim_win)
figYAxis([],'Auditory (dF/F)',scatLim_win)
figAxForm
title(sprintf('Late Window, All Resp. Cells (%s/%s), p = %s',...
    num2str(sum(cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),...
    num2str(length(cellInfo.firstRespCells)),...
    num2str(round(p,2,'significant'))))


subplot 325
for iav = 1:2
    y = mean(antiAnalysis.lateCycTC{iav}...
        ((tcStartFrame:end),(cellInfo.lateCycRespCells)),2);
    yerr = ste(antiAnalysis.lateCycTC{iav}...
        ((tcStartFrame:end),(cellInfo.lateCycRespCells)),2);
    hold on
    shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
end
figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
figYAxis([],'dF/F',cycTCLim)  
vline(lateWinTT,'k--')
hline(0,'k:')
vline(respWinTT,'k--')
figAxForm
title(sprintf('Late Cycle Resp. Cells (%s/%s)',...
    num2str(sum(cellInfo.lateCycRespCells)),...
    num2str(length(cellInfo.firstRespCells))))

subplot 326
x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,...
    (cellInfo.lateCycRespCells)),1);
xerr = ste(x,2);
y = mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,...
    (cellInfo.lateCycRespCells)),1);
yerr = ste(y,2);
plot(x,y,'.')
hold on
errorbar(mean(x),mean(y),yerr,yerr,xerr,xerr,'.')
plot(scatLim_cyc,scatLim_cyc,'k--')
plot(scatLim_cyc,[0 0],'k:')
plot([0 0],scatLim_cyc,'k:')
[~,p] = ttest(x,y);
figXAxis([],'Visual (dF/F)',scatLim_cyc)
figYAxis([],'Auditory (dF/F)',scatLim_cyc)
figAxForm
title(sprintf('Late Cycles, Late Cyc. Resp. Cells (%s/%s), p = %s',...
    num2str(sum(cellInfo.lateCycRespCells)),...
    num2str(length(cellInfo.firstRespCells)),...
    num2str(round(p,2,'significant'))))

print([fnout 'tcLongTrialsAndLateCycWithQuant_LateCycRespCells'],'-dpdf','-fillpage')

% selectivity
figure
subplot 121
y = antiAnalysis.lateCycSI(cellInfo.lateCycRespCells);
h = cdfplot(y);
hold on
vline(mean(y),'k-')
figXAxis([],'Selectivity Index',siLim)
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
title(sprintf('Late Cyc. Resp. Cells, mean = %s, ste = %s',...
    num2str(round(mean(y),2,'significant')),...
    num2str(round(ste(y,2),2,'significant'))))
subplot 122
y = antiAnalysis.lateCycSI_shuff(cellInfo.lateCycRespCells);
h = cdfplot(y);
hold on;
vline(mean(y),'k-')
figXAxis([],'Selectivity Index of Shuffle',siLim)
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
[~,p] = ttest(antiAnalysis.lateCycSI(cellInfo.lateCycRespCells), ...
    antiAnalysis.lateCycSI_shuff(cellInfo.lateCycRespCells));
title(sprintf('mean=%s, ste=%s, ttest p=%s',...
    num2str(round(mean(y),2,'significant')),...
    num2str(round(ste(y,2),2,'significant')),...
    num2str(round(p,2,'significant'))))

print([fnout 'selectivityCDF'],'-dpdf')

% print([fnout 'VAauROCCDF'],'-dpdf')

% selectivity for each mouse
meanSIPerExpt = cellfun(@(x,y) mean(x(y==1)),siPerExpt,{respCellsExpt.lateCycRespCells});
nCellsPerExpt = cellfun(@(x) sum(x==1),{respCellsExpt.lateCycRespCells});
nCellsNorm = nCellsPerExpt/30;
nCellsNorm(nCellsPerExpt >=30) = 1;
setFigParams4Print('landscape')
figure
if strcmp(ds,'FSAV_attentionV1')
    colors = brewermap(nMice,'Set2');
elseif strcmp(ds,'FSAV_attentionV1_noAttn')
    colors = brewermap(nMice,'Dark2');
else
    colors = brewermap(nMice,'Set1');
end
subplot 221
siPerMouse = nan(1,nMice);
siPerMouse_err = nan(1,nMice);
for im = 1:nMice
    hold on
    ind = strcmp({respCellsExpt.mouse},mice(im));
    ind = ind & nCellsNorm ~= 0;
    if sum(ind) ~= 0
        h=scatter(ones(1,sum(ind))*im,meanSIPerExpt(ind),1000.*(nCellsNorm(ind)),'Marker','.');
        h.MarkerFaceColor = colors(im,:);
        h.MarkerEdgeColor = colors(im,:);
        y = cellfun(@(x,y) x(y==1),siPerExpt(ind),{respCellsExpt(ind).lateCycRespCells},'unif',0);
        siPerMouse(im) = nanmean(cell2mat(y));
        siPerMouse_err(im) = ste(cell2mat(y),2);
        if length(cell2mat(y)) < minCellN
            continue
        end
        h=errorbar(im,mean(cell2mat(y)),ste(cell2mat(y),2),'.','MarkerSize',20);
        h.Color = colors(im,:);
    end
end
figXAxis([],'Mouse #',[0 nMice+1],1:nMice,mice)
figYAxis([],'SI',[-2.5 2.5])
figAxForm
title({'Mean SI; single points are expts';...
    'error points are pooled across expt';...
    'size is n cells (max size is 30 cells)'})
hline(0,'k:')
hline(mean(antiAnalysis.lateCycSI(cellInfo.lateCycRespCells)),'r-')

subplot 223
for im = 1:nMice
    hold on
    h=errorbar(im,siPerMouse(im),siPerMouse_err(im),'.','MarkerSize',20);
    h.Color = colors(im,:);
end
figXAxis([],'Mouse #',[0 nMice+1],1:nMice,mice)
figYAxis([],'SI',[-2.5 2.5])
figAxForm
hline(0,'k:')
hline(mean(antiAnalysis.lateCycSI(cellInfo.lateCycRespCells)),'r-')

if strcmp(ds,'FSAV_attentionV1')||strcmp(ds,'FSAV_attentionV1_noAttn')
    load(fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs','bxStats.mat'))
    subplot 222
    ind = ismember(bxStats.mouseNames,mice);
    [~, bxDataSortInd] = sort(bxStats.mouseNames(ind));
    hrDiffAll = bxStats.av(visualTrials).matchedHRDiff(ind);
    hrDiff = hrDiffAll(bxDataSortInd);
    for im = 1:nMice
        hold on
        h=errorbar(hrDiff(im),siPerMouse(im),siPerMouse_err(im),'.','MarkerSize',20);
        h.Color = colors(im,:);
    end
    figXAxis([],'Val - Inv Vis. HR (Training Data)',[-0.2 0.5])
    figYAxis([],'SI',[-2.5 2.5])
    figAxForm
    title('Mean SI pooled across expt for late resp neurons')
    hline(0,'k:')
    vline(0,'k:')
    hline(mean(antiAnalysis.lateCycSI(cellInfo.lateCycRespCells)),'r-')
    
    subplot 224
    ind = ismember(bxStats.mouseNames,mice);
    [~, bxDataSortInd] = sort(bxStats.mouseNames(ind));
    hrDiffAll = bxStats.av(auditoryTrials).matchedHRDiff(ind);
    hrDiff = hrDiffAll(bxDataSortInd);
    for im = 1:nMice
        hold on
        h=errorbar(hrDiff(im),siPerMouse(im),siPerMouse_err(im),'.','MarkerSize',20);
        h.Color = colors(im,:);
    end
    figXAxis([],'Val - Inv Aud. HR (Training Data)',[-0.2 0.5])
    figYAxis([],'SI',[-2.5 2.5])
    figAxForm
    title('Mean SI pooled across expt for late resp neurons')
    hline(0,'k:')
    vline(0,'k:')
    hline(mean(antiAnalysis.lateCycSI(cellInfo.lateCycRespCells)),'r-')
end
print([fnout 'perExptSI'],'-dpdf','-fillpage')

%selectivity for each cycle
figure
subplot 221
d = cellfun(@(x) x(cellInfo.lateCycRespCells),antiAnalysis.eaCycSI,'unif',0);
p = anova1(cell2mat(d')',[],'off');
y = cell2mat(cellfun(@(x) mean(x),d,'unif',0));
yerr = cell2mat(cellfun(@(x) ste(x,2),d,'unif',0));
hold on
h=errorbar(1:nCycles,y,yerr,'k.','MarkerSize',10);
figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,...
    cat(2,cellfun(@num2str,num2cell(1:nCycles-1),'unif',0),[num2str(nCycles) '+']))
figYAxis([],'SI',[-1.5 1.5])  
hline(0,'k:')
figAxForm
title({'mean across all late resp cells';sprintf('One-way ANOVA, p=%s',...
    num2str(round(p,2,'significant')))})
subplot 222
siEaMs = nan(nCycles,nMice);
L = [];
for im = 1:nMice
    hold on
    ind = strcmp({respCellsExpt.mouse},mice(im));
    ind = ind & nCellsNorm ~= 0;
    sidata = cell(1,nCycles);
    if sum(ind) ~= 0
        for iexp = 1:nexp
            if ind(iexp)
                sidata = cellfun(@(z,a) cat(2,z,a),sidata,cellfun(@(x) x(respCellsExpt(iexp).lateCycRespCells==1), ...
                    antiDataExpt(iexp).eaCycSI(1:nCycles),'unif',0),'unif',0);
            end
        end
        h = shadedErrorBar_chooseColor(1:nCycles,...
            cellfun(@(x) mean(x),sidata(1:nCycles)),...
            cellfun(@(x) ste(x,2),sidata(1:nCycles)),colors(im,:));
%         h.Color = colors(im,:);
        L(im) = h.mainLine;
    end
    siEaMs(:,im) = cellfun(@(x) mean(x),sidata(1:nCycles));
end
figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,...
    cat(2,cellfun(@num2str,num2cell(1:nCycles-1),'unif',0),[num2str(nCycles) '+']))
figYAxis([],'SI',[-3 3])  
hline(0,'k:')
figAxForm
legend(L,mice)

subplot 223
ind = cellInfo.lateCycRespCells|antiAnalysis.earlyTest';
d = cellfun(@(x,y) x(ind==1),antiAnalysis.eaCycSI,'unif',0);
% d = cellfun(@(x,y) x(y==1),antiAnalysis.eaCycSI,antiAnalysis.eaCycTest,'unif',0);
% p = anova1(cell2mat(d')',[],'off');
y = cell2mat(cellfun(@(x) mean(x),d,'unif',0));
yerr = cell2mat(cellfun(@(x) ste(x,2),d,'unif',0));
hold on
h=errorbar(1:nCycles,y,yerr,'k.','MarkerSize',10);
figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,...
    cat(2,cellfun(@num2str,num2cell(1:nCycles-1),'unif',0),[num2str(nCycles) '+']))
figYAxis([],'SI',[-1.5 1.5])  
hline(0,'k:')
figAxForm
title(sprintf('early or late resp cells, n=%s',num2str(sum(ind))))

subplot 224
h=errorbar(1:nCycles,mean(siEaMs,2),ste(siEaMs,2),'k.','MarkerSize',10);
figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,...
    cat(2,cellfun(@num2str,num2cell(1:nCycles-1),'unif',0),[num2str(nCycles) '+']))
figYAxis([],'SI',[-1.5 1.5])  
hline(0,'k:')
figAxForm
title('mean across mice')
print([fnout 'perStimAttnEffect'],'-dpdf','-fillpage')

figure
suptitle('Cells Resp. to Ea. Stim')
% ind = cellInfo.lateCycRespCells;
% ind = sum(cell2mat(antiAnalysis.eaCycTest'),1)'>0;

ind = cellInfo.lateCycRespCells|antiAnalysis.earlyTest';
for icyc = 1:nCycles
    subplot(2,4,icyc)
    for iav = 1:2
%         ind = antiAnalysis.eaCycTest{icyc}==1;
        y = mean(antiAnalysis.eaCycTC{iav,icyc}(tcStartFrame:end,ind),2)-...
            mean(mean(antiAnalysis.eaCycTC{iav,icyc}(basewin_0,ind),2));
        yerr = ste(antiAnalysis.eaCycTC{iav}(tcStartFrame:end,ind),2);
        hold on
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
    end
    figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',eaCycTCLim)  
    hline(0,'k:')
    vline(respWinTT,'k--')
    figAxForm
    title(sprintf('n=%s/%s',...
        num2str(sum(ind)),...
        num2str(length(ind))))
end
print([fnout 'perStimCycTC'],'-dpdf','-fillpage')

figure
ind = antiAnalysis.earlyTest'==1 &~cellInfo.lateCycRespCells==1;
subplot 121
for iav = 1:2
%         ind = antiAnalysis.eaCycTest{icyc}==1;
    y = mean(antiAnalysis.earlyCycTC{iav}(tcStartFrame:end,ind),2);
    yerr = ste(antiAnalysis.earlyCycTC{iav}(tcStartFrame:end,ind),2);
    hold on
    shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
end
[~,p] = ttest(mean(antiAnalysis.earlyCycTC{visualTrials}(respwin,ind),1),...
    mean(antiAnalysis.earlyCycTC{auditoryTrials}(respwin,ind),1));
figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
figYAxis([],'dF/F',eaCycTCLim)  
hline(0,'k:')
vline(respWinTT,'k--')
figAxForm
title(sprintf('Early Stim Resp.,n=%s/%s, p=%s',...
    num2str(sum(ind)),...
    num2str(length(ind)),num2str(round(p,2,'significant'))))
subplot 122
ind = cellInfo.lateCycRespCells==1;
for iav = 1:2
%         ind = antiAnalysis.eaCycTest{icyc}==1;
    y = mean(antiAnalysis.lateCycTC{iav}(tcStartFrame:end,ind),2);
    yerr = ste(antiAnalysis.lateCycTC{iav}(tcStartFrame:end,ind),2);
    hold on
    shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
end
[~,p] = ttest(mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1),...
    mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,ind),1));
figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
figYAxis([],'dF/F',eaCycTCLim)  
hline(0,'k:')
vline(respWinTT,'k--')
figAxForm
title(sprintf('Late Stim Resp., n=%s/%s, p=%s',...
    num2str(sum(ind)),...
    num2str(length(ind)),num2str(round(p,2,'significant'))))
print([fnout 'earlyLateStimCycTC'],'-dpdf','-fillpage')

%% heatmap for all behavior experiments
if strcmp(ds,'FSAV_attentionV1')
    hm_noAttn = load(fullfile(rc.ashley, 'Manuscripts',...
        'Attention V1','Matlab Figs',...
        ['attentionV1_noAttn_' attnAnalysisDate '__heatmaps.mat']));
    hm_attn = load(fullfile(rc.ashley, 'Manuscripts',...
        'Attention V1','Matlab Figs',...
        ['attentionV1_' attnAnalysisDate '__heatmaps.mat']));
    
    hmAll = cat(1,hm_attn.hm,hm_noAttn.hm);
    [~,sortInd] = sort(mean(hmAll(:,lateWinFr),2));
    stimRespInd = cat(1,hm_attn.stimRespCellsBW,hm_noAttn.stimRespCellsBW);
    lateRespInd = cat(1,hm_attn.lateStimRespCellsBW,hm_noAttn.lateStimRespCellsBW);
    exCellInd = cat(1,hm_attn.exCellInd,hm_noAttn.exCellInd);
    
    figure
    colormap(brewermap([],'*RdBu'));
    imagesc(flipud(hmAll(sortInd,tcStartFrame:end)))
    hold on
    figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabelFr_long)
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    caxis(hmLim)
    hline(find(flipud(exCellInd(sortInd))),'k-')
    title(sprintf('Total Cells=%s',num2str(length(sortInd))))
    print([fnout 'heatmapAllTrialsAllCells_allBx'],'-dpdf','-fillpage')
    
    figure
    colormap gray
    subplot 121
    imagesc(flipud(stimRespInd(sortInd)))
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    title('Resp. to early or late distractor')
    subplot 122
    imagesc(flipud(lateRespInd(sortInd)))
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    title('Resp. to late distractor')
    print([fnout 'cellIDforHeatmap_allBx'],'-dpdf','-fillpage')
end
%% create a structure with stats from analyses

if strcmp(ds,'FSAV_attentionV1')
    load(fullfile(rc.ashley, 'Manuscripts',...
        'Attention V1','Matlab Figs',...
        ['attentionV1_noAttn_' attnAnalysisDate '_imgStats.mat']));
    imgStats_noAttn = imgStats;
%     imgStats_naive = load(['V1_100ms_naive_' analysisDate '_imgStats.mat']);
end

imgStats = struct;
% imgStats.nTrialsPerExpt = nTrialsPerExpt;
% fprintf('%s +/- %s (%s-%s) trials per experiment\n',num2str(mean(nTrialsPerExpt)),...
%     num2str(ste(nTrialsPerExpt,2)), num2str(min(nTrialsPerExpt)),...
%     num2str(max(nTrialsPerExpt)))
imgStats.firstWinTimeMs = round(([firstWinFr(1) firstWinFr(end)]...
    - (nBaselineFr+nVisDelayFr)).*(1000/frameRateHz),4,'significant');
imgStats.lateWinTimeMs = round(lateWinTT,4,'significant');

%cell numbers
imgStats.nCells.totalCells = length(cellInfo.firstRespCells);
imgStats.nCells.allResp = sum(cellInfo.firstRespCells|cellInfo.lateCycRespCells|cellInfo.lateRespCells);
imgStats.nCells.firstResp = sum(cellInfo.firstRespCells);
imgStats.nCells.lateCycResp = sum(cellInfo.lateCycRespCells);
imgStats.nCells.lateOnlyResp = sum(...
    ~cellInfo.firstRespCells & (cellInfo.lateCycRespCells|cellInfo.lateRespCells));
% imgStats.nCells.lateTCResp = sum(cellInfo.lateRespCells);
imgStats.nCells.lateSupp = sum(cellInfo.lateSuppCells & ...
    ~cellInfo.lateCycRespCells & ~cellInfo.lateRespCells & ~cellInfo.firstRespCells);
imgStats.nCells.targetResp = sum(cellInfo.targetRespCells);
imgStats.nCells

x = mean(antiAnalysis.longTC{visualTrials}(firstWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),1);
y = mean(antiAnalysis.longTC{auditoryTrials}(firstWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),1);
imgStats.av(visualTrials).firstWinResp = mean(x);
imgStats.av(visualTrials).firstWinRespErr = ste(x,2);
imgStats.av(auditoryTrials).firstWinResp = mean(y);
imgStats.av(auditoryTrials).firstWinRespErr = ste(y,2);
[~,imgStats.firstWinTest] = ttest(x,y);
x = mean(antiAnalysis.longTC{visualTrials}(lateWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),1);
y = mean(antiAnalysis.longTC{auditoryTrials}(lateWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),1);
imgStats.av(visualTrials).lateWinResp = mean(x);
imgStats.av(visualTrials).lateWinRespErr = ste(x,2);
imgStats.av(auditoryTrials).lateWinResp = mean(y);
imgStats.av(auditoryTrials).lateWinRespErr = ste(y,2);
[~,imgStats.lateWinTest] = ttest(x,y);

% fprintf('Visual - Mean/Err First Win: %s/%s\n',...
%     num2str(round(imgStats.av(visualTrials).firstWinResp,2,'significant')),...
%     num2str(round(imgStats.av(visualTrials).firstWinRespErr,2,'significant')))
% fprintf('Auditory - Mean/Err First Win: %s/%s\n',...
%     num2str(round(imgStats.av(auditoryTrials).firstWinResp,2,'significant')),...
%     num2str(round(imgStats.av(auditoryTrials).firstWinRespErr,2,'significant')))
fprintf('First Win Test, p=%s\n',...
    num2str(round(imgStats.firstWinTest,2,'significant')))
% fprintf('Visual - Mean/Err Late Win: %s/%s\n',...
%     num2str(round(imgStats.av(visualTrials).lateWinResp,2,'significant')),...
%     num2str(round(imgStats.av(visualTrials).lateWinRespErr,2,'significant')))
% fprintf('Auditory - Mean/Err Late Win: %s/%s\n',...
%     num2str(round(imgStats.av(auditoryTrials).lateWinResp,2,'significant')),...
%     num2str(round(imgStats.av(auditoryTrials).lateWinRespErr,2,'significant')))
fprintf('Late Win Test, p=%s\n',...
    num2str(round(imgStats.lateWinTest,2,'significant')))

x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,...
    (cellInfo.lateCycRespCells)),1);
y = mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,...
    (cellInfo.lateCycRespCells)),1);
imgStats.av(visualTrials).lateCycResp = mean(x);
imgStats.av(visualTrials).lateCycRespErr = ste(x,2);
imgStats.av(auditoryTrials).lateCycResp = mean(y);
imgStats.av(auditoryTrials).lateCycRespErr = ste(y,2);
[~,imgStats.lateCycTest] = ttest(x,y);

% fprintf('Visual - Mean/Err Late Cyc: %s/%s\n',...
%     num2str(round(imgStats.av(visualTrials).lateCycResp,2,'significant')),...
%     num2str(round(imgStats.av(visualTrials).lateCycRespErr,2,'significant')))
% fprintf('Auditory - Mean/Err Late Cyc: %s/%s\n',...
%     num2str(round(imgStats.av(auditoryTrials).lateCycResp,2,'significant')),...
%     num2str(round(imgStats.av(auditoryTrials).lateCycRespErr,2,'significant')))
fprintf('Late Cyc Test, p=%s\n',...
    num2str(round(imgStats.lateCycTest,2,'significant')))

imgStats.respCellsSI = antiAnalysis.lateCycSI(cellInfo.lateCycRespCells);
[~,imgStats.respCellsSITest] = ttest(antiAnalysis.lateCycSI(...
    cellInfo.lateCycRespCells));
% imgStats.nAVModCells = sum(antiAnalysis.lateCycAVTest(cellInfo.lateCycRespCells));
% fprintf('Mean/Err SI: %s/%s\n',...
%     num2str(round(mean(imgStats.respCellsSI,2),2,'significant')),...
%     num2str(round(ste(imgStats.respCellsSI,2),2,'significant')))
fprintf('SI Test, p=%s\n',...
    num2str(round(imgStats.respCellsSITest,2,'significant')))

imgStats.siPerMouse = siPerMouse;
imgStats.mouseID = mice;

% save([fnout 'imgStats'],'imgStats')

if strcmp(ds,'FSAV_attentionV1')  
    siAllBxMice = cat(2,imgStats.siPerMouse,imgStats_noAttn.siPerMouse);
    allMiceID = cat(2,imgStats.mouseID,imgStats_noAttn.mouseID);
    [~,siSortInd] = sort(allMiceID);
    [~,bxSortInd] = sort(bxStats.mouseNames);
    [bxSIcorr,p] = corr(siAllBxMice(siSortInd)',bxStats.allMatchedHRDiff(bxSortInd));
    imgStats.bxVsSICorr = [bxSIcorr,p];
    fprintf('SI vs. Val-Inv HR, r=%s, p=%s\n',...
        num2str(round(imgStats.bxVsSICorr(1),2,'significant')),...
        num2str(round(imgStats.bxVsSICorr(2),2,'significant')))
    [visBxSIcorr,p] = corr(siAllBxMice(siSortInd)',...
        bxStats.av(visualTrials).matchedHRDiff(bxSortInd));
    imgStats.visBxVsSICorr = [visBxSIcorr,p];
    fprintf('SI vs. Vis Only Val-Inv HR, r=%s, p=%s\n',...
        num2str(round(imgStats.visBxVsSICorr(1),2,'significant')),...
        num2str(round(imgStats.visBxVsSICorr(2),2,'significant')))
    
    [~,imgStats.respCellsSITestVsNoAttn] = kstest2(imgStats.respCellsSI, ...
        imgStats_noAttn.respCellsSI,'tail','smaller');
    fprintf('Kolmogorov-Smirnov Test SI, attention vs. no attention, p = %s\n',...
        num2str(round(imgStats.respCellsSITestVsNaive,2,'significant')))
    
    [~,imgStats.respCellsSITestVsNaive] = kstest2(imgStats.respCellsSI, ...
        imgStats_naive.respCellsSI,'tail','smaller');
    fprintf('Kolmogorov-Smirnov Test SI, attention vs. naive, p = %s\n',...
        num2str(round(imgStats.respCellsSITestVsNaive,2,'significant')))
end
save([fnout 'imgStats'],'imgStats')

%%
% example cells time-courses
if strcmp(ds,'FSAV_attentionV1')
    setFigParams4Print('landscape')
    exCellInd = [exampleCell_1,exampleCell_3];
    figure
    for icell = 1:length(exCellInd)
        exampleCell = exCellInd(icell);
        if icell == 1
            plotOffset = 0;
        else
            plotOffset = plotOffset+4;
        end
        subplot(3,4,1+plotOffset)
        y = antiAnalysis.longTC{allTrialsInd}(tcStartFrame:end,exampleCell);
        yerr = antiAnalysis.longTCErr{allTrialsInd}(tcStartFrame:end,exampleCell);
        shadedErrorBar_chooseColor(tt_longTC,y,yerr,[0 0 0]);
        hold on
        figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
        if icell == 3
            figYAxis([],'dF/F',exCellTCLim-0.06)
        else
            figYAxis([],'dF/F',exCellTCLim)  
        end
        vline(lateWinTT,'k--')
        hline(0,'k:')
        figAxForm([],0)
        title(sprintf('Example Cell #%s',num2str(exampleCell)))
        
        subplot(3,4,2+plotOffset)
        y = antiAnalysis.earlyCycTC{allTrialsInd}(tcStartFrame:end,exampleCell);
        yerr = antiAnalysis.earlyCycTCErr{allTrialsInd}(tcStartFrame:end,exampleCell);
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,[0 0 0]);
        hold on
        figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
        if icell == 3
            figYAxis([],'dF/F',exCellTCLim-0.06)
        else
            figYAxis([],'dF/F',exCellTCLim)  
        end
        vline(lateWinTT,'k--')
        hline(0,'k:')
        figAxForm([],0)
        title(sprintf('Example Cell #%s',num2str(exampleCell)))
        
        subplot(3,4,3+plotOffset)
        y = antiAnalysis.lateCycTC{allTrialsInd}(tcStartFrame:end,exampleCell);
        yerr = antiAnalysis.lateCycTCErr{allTrialsInd}(tcStartFrame:end,exampleCell);
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,[0 0 0]);
        hold on
        figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
        if icell == 3
            figYAxis([],'dF/F',exCellTCLim-0.06)
        else
            figYAxis([],'dF/F',exCellTCLim)  
        end
        vline(lateWinTT,'k--')
        hline(0,'k:')
        figAxForm([],0)
        title(sprintf('Example Cell #%s',num2str(exampleCell)))

        subplot(3,4,4+plotOffset)
        y = cellInfo.oriResp(exampleCell,:);
        yerr = cellInfo.oriRespErr(exampleCell,:);
        errorbar(orientations,y,yerr,'.')
        hold on
        x = 0:180;
        y = cellInfo.oriFit(exampleCell,:);
        plot(x,y,'-')
        figXAxis([],'Orienation (deg)',[-10 190])
        figYAxis([],'dF/F',oriRespLim)
        figAxForm
        title(sprintf('Passive Ori. Tuning, Cell #%s',num2str(exampleCell)));
    end
    print([fnout 'exampleCellsTCWithTuning'],'-dpdf','-fillpage')
elseif strcmp(ds,'FSAV_attentionV1_noAttn')
    setFigParams4Print('landscape')
    exCellInd = [exampleCell_2];
    figure
    for icell = 1:length(exCellInd)
        exampleCell = exCellInd(icell);
        if icell == 1
            plotOffset = 0;
        else
            plotOffset = plotOffset+4;
        end
        subplot(3,4,1+plotOffset)
        y = antiAnalysis.longTC{allTrialsInd}(tcStartFrame:end,exampleCell);
        yerr = antiAnalysis.longTCErr{allTrialsInd}(tcStartFrame:end,exampleCell);
        shadedErrorBar_chooseColor(tt_longTC,y,yerr,[0 0 0]);
        hold on
        figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
        if icell == 3
            figYAxis([],'dF/F',exCellTCLim-0.06)
        else
            figYAxis([],'dF/F',exCellTCLim)  
        end
        vline(lateWinTT,'k--')
        hline(0,'k:')
        figAxForm([],0)
        title(sprintf('Example Cell #%s',num2str(exampleCell)))
        
        subplot(3,4,2+plotOffset)
        y = antiAnalysis.earlyCycTC{allTrialsInd}(tcStartFrame:end,exampleCell);
        yerr = antiAnalysis.earlyCycTCErr{allTrialsInd}(tcStartFrame:end,exampleCell);
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,[0 0 0]);
        hold on
        figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
        if icell == 3
            figYAxis([],'dF/F',exCellTCLim-0.06)
        else
            figYAxis([],'dF/F',exCellTCLim)  
        end
        vline(lateWinTT,'k--')
        hline(0,'k:')
        figAxForm([],0)
        title(sprintf('Example Cell #%s',num2str(exampleCell)))
        
        subplot(3,4,3+plotOffset)
        y = antiAnalysis.lateCycTC{allTrialsInd}(tcStartFrame:end,exampleCell);
        yerr = antiAnalysis.lateCycTCErr{allTrialsInd}(tcStartFrame:end,exampleCell);
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,[0 0 0]);
        hold on
        figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
        if icell == 3
            figYAxis([],'dF/F',exCellTCLim-0.06)
        else
            figYAxis([],'dF/F',exCellTCLim)  
        end
        vline(lateWinTT,'k--')
        hline(0,'k:')
        figAxForm([],0)
        title(sprintf('Example Cell #%s',num2str(exampleCell)))

        subplot(3,4,4+plotOffset)
        y = cellInfo.oriResp(exampleCell,:);
        yerr = cellInfo.oriRespErr(exampleCell,:);
        errorbar(orientations,y,yerr,'.')
        hold on
        x = 0:180;
        y = cellInfo.oriFit(exampleCell,:);
        plot(x,y,'-')
        figXAxis([],'Orienation (deg)',[-10 190])
        figYAxis([],'dF/F',oriRespLim)
        figAxForm
        title(sprintf('Passive Ori. Tuning, Cell #%s',num2str(exampleCell)));
    end
    print([fnout 'exampleCellsTCWithTuning'],'-dpdf','-fillpage')
end
%% plot orientation analysis (Figure 3)
if strcmp(ds,'FSAV_attentionV1')
    minCellN_SIFRmatch = 34;
elseif strcmp(ds,'FSAV_attentionV1_noAttn')
    minCellN_SIFRmatch = 12;
elseif strcmp(ds,'FSAV_V1_audControl')|strcmp(ds,'FSAV_V1_naive_GCaMP6m')
    minCellN_SIFRmatch = 1;
else
    minCellN_SIFRmatch = 1;
end
% adaptAnalysis = struct;
% ind = cellInfo.firstRespCells & cellInfo.minRespCells;
% adaptAnalysis.firstRespCells = cellfun(@(x) mean(x(ind)),antiAnalysis.adapt);
% adaptAnalysis.firstRespCellsErr = cellfun(@(x) ste(x(ind),2),antiAnalysis.adapt);
% [~,adaptAnalysis.testAV] = ttest(antiAnalysis.adapt{visualTrials}(ind),...
%     antiAnalysis.adapt{auditoryTrials}(ind),'tail','right');
% ind2 = ind & cellInfo.isTuned;
% adaptAnalysis.oriGroupsTest = anova1(antiAnalysis.adapt{3}(ind2),...
%     cellInfo.oriPref(ind2),'off');

% get first, late win, and late cycle response for passive ori tuning groups
oriGroups = struct;
oriGroups.n = nan(1,nOri);
oriGroups.nShortCyc = nan(1,nOri);
oriGroups.nTarOrDist = nan(1,nOri);
oriGroups.nTar = nan(1,nOri);
oriGroups.lateCycSI = nan(1,nOri);
oriGroups.lateCycSIErr = nan(1,nOri);
% oriGroups.lateCycAVauROC = nan(1,nOri);
% oriGroups.lateCycAVauROCErr = nan(1,nOri);
% oriGroups.lateCycAVauROCData = cell(1,nOri);
oriGroups.lateCycRespAll = nan(1,nOri);
oriGroups.lateCycRespAllErr = nan(1,nOri);
oriGroups.lateCycRespAllData = cell(1,nOri);
oriGroups.firstResp = nan(2,nOri);
oriGroups.firstRespErr = nan(2,nOri);
oriGroups.lateWin = nan(2,nOri);
oriGroups.lateWinErr = nan(2,nOri);
oriGroups.lateCycResp = nan(2,nOri);
oriGroups.lateCycRespErr = nan(2,nOri);
oriGroups.lateCycRespDiff = nan(1,nOri);
oriGroups.lateCycRespDiffErr = nan(1,nOri);
oriGroups.firstTC = cell(2,nOri);
oriGroups.firstTCErr = cell(2,nOri);
oriGroups.longTC = cell(2,nOri);
oriGroups.longTCErr = cell(2,nOri);
oriGroups.cycTC = cell(2,nOri);
oriGroups.cycTCErr = cell(2,nOri);
oriGroups.lateCycSITestEaOri = nan(1,nOri);
oriGroups.lateCycTest = nan(1,nOri);
% oriGroups.lateWinTest = nan(1,nOri);
oriGroups.targetTC = cell(1,nOri);
oriGroups.targetTCErr = cell(1,nOri);
oriGroups.targetTuningResp = cell(1,nOri);
oriGroups.targetTuningRespErr = cell(1,nOri);
oriGroups.targetTuningStim = cell(1,nOri);
oriGroups.targetTuningStimErr = cell(1,nOri);
% oriGroups.firstStimRespForTargetAnalysis = nan(1,nOri);
% oriGroups.firstStimRespErrForTargetAnalysis = nan(1,nOri);
oriGroups.matchSIxFR = nan(1,nOri);
oriGroups.matchSIxFRErr = nan(1,nOri);
oriGroups.matchedFR = nan(1,nOri);
oriGroups.matchedFRErr = nan(1,nOri);
oriGroups.matchSIxFRTestData = cell(1,nOri);
oriGroups.taskTuningData = cell(nOri,3);
% adaptAnalysis.oriGroupsN = nan(1,nOri);
% adaptAnalysis.oriGroups = nan(3,nOri);
% adaptAnalysis.oriGroupsErr = nan(3,nOri);
for iori = 1:nOri
    ind = cellInfo.isTuned & cellInfo.oriPref == iori &...
        (cellInfo.lateCycRespCells);
%     ind2 = cellInfo.isTuned & cellInfo.oriPref == iori & ...
%         cellInfo.firstRespCells & cellInfo.minRespCells;
    ind_tarAndDist = cellInfo.isTuned & cellInfo.oriPref == iori &...
        (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
    ind_tar = cellInfo.isTuned & cellInfo.oriPref == iori & cellInfo.targetRespCells;
    oriGroups.n(iori) = sum(ind);
    oriGroups.nShortCyc(iori) = sum(ind & cellInfo.isShortCycExpt);
    oriGroups.nTarOrDist(iori) = sum(ind_tarAndDist);
    oriGroups.nTar(iori) = sum(ind_tar);
%     adaptAnalysis.oriGroupsN(iori) = sum(ind2);
    oriGroups.lateCycSI(iori) = mean(antiAnalysis.lateCycSI(ind));
    oriGroups.lateCycSIErr(iori) = ste(antiAnalysis.lateCycSI(ind),2);
%     oriGroups.lateCycAVauROC(iori) = mean(antiAnalysis.lateCycAVauROC(ind),2);
%     oriGroups.lateCycAVauROCErr(iori) = ste(antiAnalysis.lateCycAVauROC(ind),2);
%     oriGroups.lateCycAVauROCData{iori} = antiAnalysis.lateCycAVauROC(ind);
    [~,oriGroups.lateCycSITestEaOri(iori)] = ttest(antiAnalysis.lateCycSI(ind));
    [~,oriGroups.lateWinTest(iori)] = ttest(mean(antiAnalysis.longTC{visualTrials}...
        (lateWinFr,ind & cellInfo.isShortCycExpt),1),...
        mean(antiAnalysis.longTC{auditoryTrials}...
        (lateWinFr,ind & cellInfo.isShortCycExpt),1));
    [~,oriGroups.lateCycTest(iori)] = ttest(mean(antiAnalysis.lateCycTC{visualTrials}...
        (respwin,ind),1),mean(antiAnalysis.lateCycTC{auditoryTrials}...
        (respwin,ind),1));
    oriGroups.lateCycRespDiff(iori) = mean(mean(antiAnalysis.lateCycTC{visualTrials}...
        (respwin,ind),1) - mean(antiAnalysis.lateCycTC{auditoryTrials}...
        (respwin,ind),1),2);
    oriGroups.lateCycRespDiffErr(iori) = ste(mean(antiAnalysis.lateCycTC{visualTrials}...
        (respwin,ind),1) - mean(antiAnalysis.lateCycTC{auditoryTrials}...
        (respwin,ind),1),2);
    
    oriGroups.lateCycRespAll(iori) = mean(mean(...
            antiAnalysis.lateCycTC{3}(respwin,ind),1));
    oriGroups.lateCycRespAllErr(iori) = ste(mean(...
            antiAnalysis.lateCycTC{3}(respwin,ind),1),2);
    oriGroups.lateCycRespAllData{iori} = mean(...
            antiAnalysis.lateCycTC{3}(respwin,ind),1);
    
    oriGroups.targetTC{iori} = cell2mat(cellfun(@(x) mean(x(:,ind_tar),2),...
        targetAnalysis.tc(1:2,visualTrials),'unif',0)');
    oriGroups.targetTCErr{iori} = cell2mat(cellfun(@(x) ste(x(:,ind_tar),2),...
        targetAnalysis.tc(1:2,visualTrials),'unif',0)');
    
%     oriGroups.targetTuningResp{iori} = cat(1,mean(mean(...
%             antiAnalysis.lateCycTC{visualTrials}(respwin_target,ind_tarAndDist),1)),cellfun(@(x) ...
%         mean(mean(x(respwin_target,ind_tarAndDist),1) - mean(x(basewin_0,ind_tarAndDist),1)),...
%         targetAnalysis.tc(1:2,visualTrials)));
%     oriGroups.targetTuningRespErr{iori} = cat(1,ste(mean(...
%             antiAnalysis.lateCycTC{visualTrials}(respwin_target,ind_tarAndDist),1),2),cellfun(@(x) ...
%         ste(mean(x(respwin_target,ind_tarAndDist),1) - mean(x(basewin_0,ind_tarAndDist),1),2),...
%         targetAnalysis.tc(1:2,visualTrials)));
%     oriGroups.firstStimRespForTargetAnalysis(iori) = mean(mean(...
%             antiAnalysis.firstCycTC{visualTrials}(respwin,ind_tarAndDist),1));
%     oriGroups.firstStimRespErrForTargetAnalysis(iori) = ste(mean(...
%             antiAnalysis.firstCycTC{visualTrials}(respwin,ind_tarAndDist),1),2);
%     oriGroups.targetTuningStim{iori} = cat(1,0,...
%         mean(targetAnalysis.targets{visualTrials}(:,ind_tarAndDist),2));
%     oriGroups.targetTuningStimErr{iori} = cat(1,0,...
%         ste(targetAnalysis.targets{visualTrials}(:,ind_tarAndDist),2));
    
    oriGroups.taskTuningData{iori,1} = mean(...
            antiAnalysis.lateCycTC{visualTrials}(respwin_target,ind),1)';
    oriGroups.taskTuningData(iori,2:3) = cellfun(@(x) ...
        mean(x(respwin_target,ind),1)' - mean(x(basewin_0,ind),1)',...
        targetAnalysis.tc(1:2,visualTrials),'unif',0);
    oriGroups.targetTuningResp{iori} = cat(1,mean(mean(...
            antiAnalysis.lateCycTC{visualTrials}(respwin_target,ind),1)),cellfun(@(x) ...
        mean(mean(x(respwin_target,ind),1) - mean(x(basewin_0,ind),1)),...
        targetAnalysis.tc(1:2,visualTrials)));
    oriGroups.targetTuningRespErr{iori} = cat(1,ste(mean(...
            antiAnalysis.lateCycTC{visualTrials}(respwin_target,ind),1),2),cellfun(@(x) ...
        ste(mean(x(respwin_target,ind),1) - mean(x(basewin_0,ind),1),2),...
        targetAnalysis.tc(1:2,visualTrials)));
    oriGroups.firstStimRespForTargetAnalysis(iori) = mean(mean(...
            antiAnalysis.firstCycTC{visualTrials}(respwin,ind),1));
    oriGroups.firstStimRespErrForTargetAnalysis(iori) = ste(mean(...
            antiAnalysis.firstCycTC{visualTrials}(respwin,ind),1),2);
    oriGroups.targetTuningStim{iori} = cat(1,0,...
        mean(targetAnalysis.targets{visualTrials}(:,ind),2));
    oriGroups.targetTuningStimErr{iori} = cat(1,0,...
        ste(targetAnalysis.targets{visualTrials}(:,ind),2));
    
    lateCycResp_grp = +Dr{3}(ind);
    lateCycSI_grp = antiAnalysis.lateCycSI(ind);
    [sortFR,sortFRInd] = sort(lateCycResp_grp,'descend');
    sortSI = lateCycSI_grp(sortFRInd);
    oriGroups.matchSIxFR(iori) = mean(sortSI(1:minCellN_SIFRmatch),2);
    oriGroups.matchSIxFRErr(iori) = ste(sortSI(1:minCellN_SIFRmatch),2);
    oriGroups.matchedFR(iori) = mean(sortFR(1:minCellN_SIFRmatch),2);
    oriGroups.matchedFRErr(iori) = ste(sortFR(1:minCellN_SIFRmatch),2);
    oriGroups.matchSIxFRTestData(iori) = {sortSI(1:minCellN_SIFRmatch)};
            
    for iav = 1:2
        oriGroups.firstResp(iav,iori) = mean(mean(...
            antiAnalysis.firstCycTC{iav}(respwin,ind),1));
        oriGroups.firstRespErr(iav,iori) = ste(mean(...
            antiAnalysis.firstCycTC{iav}(respwin,ind),1),2);
        
        oriGroups.firstTC(iav,iori) = {mean(...
            antiAnalysis.firstCycTC{iav}(:,ind),2)};
        oriGroups.firstTCErr(iav,iori) = {ste(...
            antiAnalysis.firstCycTC{iav}(:,ind),2)};
        
        oriGroups.lateWin(iav,iori) = mean(mean(...
            antiAnalysis.longTC{iav}(lateWinFr,ind),1));
        oriGroups.lateWinErr(iav,iori) = ste(mean(...
            antiAnalysis.longTC{iav}(lateWinFr,ind),1),2);
        
        oriGroups.lateCycResp(iav,iori) = mean(mean(...
            antiAnalysis.lateCycTC{iav}(respwin,ind),1));
        oriGroups.lateCycRespErr(iav,iori) = ste(mean(...
            antiAnalysis.lateCycTC{iav}(respwin,ind),1),2);
        
        oriGroups.longTC{iav,iori} = mean(...
            antiAnalysis.longTC{iav}(:,ind & cellInfo.isShortCycExpt),2);
        oriGroups.longTCErr{iav,iori} = ste(...
            antiAnalysis.longTC{iav}(:,ind & cellInfo.isShortCycExpt),2);
        
        oriGroups.cycTC{iav,iori} = mean(...
            antiAnalysis.lateCycTC{iav}(:,ind),2);
        oriGroups.cycTCErr{iav,iori} = ste(...
            antiAnalysis.lateCycTC{iav}(:,ind),2);
%         adaptAnalysis.oriGroups(iav,iori) = mean(...
%             antiAnalysis.adapt{iav}(ind2),2);
%         adaptAnalysis.oriGroupsErr(iav,iori) = ste(...
%             antiAnalysis.adapt{iav}(ind2),2);
%         if iav == 2
%             adaptAnalysis.oriGroups(3,iori) = mean(...
%                 antiAnalysis.adapt{3}(ind2),2);
%             adaptAnalysis.oriGroupsErr(3,iori) = ste(...
%                 antiAnalysis.adapt{3}(ind2),2);
%         end
        
        
    end
end
ind = cellInfo.isTuned &...
    (cellInfo.lateCycRespCells);
lateCycRespDiffData = ...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1) - ...
    mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,ind),1);
[oriGroups.lateCycRespTest,~,stats] = anova1(lateCycRespDiffData,cellInfo.oriPref(ind),'off');
oriGroups.lateCycRespTestPostHoc = multcompare(stats,[],'off');
% oriGroups.lateCycSITest = anova1(antiAnalysis.lateCycSI(ind),cellInfo.oriPref(ind),'off');
[lateCycSITest,~,stats] = anova1(antiAnalysis.lateCycSI(ind),cellInfo.oriPref(ind),'off');
% lmTest = fitlm(antiAnalysis.lateCycSI(ind),cellInfo.oriPref(ind));
oriGroups.lateCycSITest = lateCycSITest;
[~,oriGroups.FRmatchedSITestEaOri] = cellfun(@ttest,oriGroups.matchSIxFRTestData);

% [lateCycAVauROCTest,~,stats] = anova1(antiAnalysis.lateCycAVauROC(ind),cellInfo.oriPref(ind),'off');
% lmTest = fitlm(antiAnalysis.lateCycSI(ind),cellInfo.oriPref(ind));
% oriGroups.lateCycAVauROCTest = lateCycAVauROCTest;
% [~,oriGroups.lateCycAVauROCTestEaOri] = cellfun(@(x) ttest(x,0.5),oriGroups.lateCycAVauROCData);

% lateCycRespAV = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1);
% [oriGroups.lateCycRespAllTest,~,stats] = anova1(lateCycRespAV,cellInfo.oriPref(ind),'off');
% oriGroups.lateCycRespTestPostHoc = multcompare(stats);

% ind = cellInfo.isTuned & (cellInfo.firstRespCells|cellInfo.lateCycRespCells | cellInfo.targetRespCells);
% popTuningData = reshape(cat(1,mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1),...
%     cell2mat(cellfun(@(x) mean(x(respwin_target,ind),1) - mean(x(basewin_0,ind),1),...
%         targetAnalysis.tc(1:2,visualTrials),'unif',0))),[3*sum(ind),1]);
% popTuningGrp1 = repmat([1:3]',[sum(ind),1]);
% popTuningGrp2 = reshape(repmat(cellInfo.oriPref(ind)',[3,1]),[3*sum(ind),1]);
% oriGroups.popTuningTest = anovan(popTuningData,{popTuningGrp1,popTuningGrp2});

% taskTuningTest = twoWayAnovaUnmatched(oriGroups.taskTuningData,{'Ori';'Stim'});
% [p,~,stats] = anova1(Tr_vis{1}(cellInfo.lateCycRespCells),...
%     cellInfo.oriPref(cellInfo.lateCycRespCells),'off');
% tbl = multcompare(stats,[],'off');

% grpID = [];
% for igrp = 1:4
%     grpID = cat(2,grpID,ones(1,length(oriGroups.matchSIxFRTestData{igrp})).*igrp);
% end
% [oriGroups.matchSIxFRTest,~,stats] = anova1(cell2mat(...
%     oriGroups.matchSIxFRTestData),grpID,'off');

figure
suptitle('Tuned, Late Distractor Responsive Neurons')
for iori = 1:nOri
    subplot(2,2,iori)
    for iav = 1:2
        y = oriGroups.cycTC{iav,iori}(tcStartFrame:end);
        yerr = oriGroups.cycTCErr{iav,iori}(tcStartFrame:end);
        hold on
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
    end
    figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',cycTCLim)  
    vline(respWinTT,'k--')
    hline(0,'k:')
    figAxForm
    title(sprintf('Pref %s, n=%s, p=%s',num2str(orientations(iori)),...
        num2str(oriGroups.n(iori)),num2str(round(oriGroups.lateCycTest(iori),2,'significant'))))
    
    
end

print([fnout 'tuningAnalysis_oriGroupsLateCycTC'],'-dpdf','-fillpage')    


figure
suptitle('First Cyc; Tuned, First or Late Distractor Responsive Neurons')
for iori = 1:nOri
    subplot(2,2,iori)
    for iav = 1:2
        y = oriGroups.firstTC{iav,iori}(tcStartFrame:end);
        yerr = oriGroups.firstTCErr{iav,iori}(tcStartFrame:end);
        hold on
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
    end
    figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',firstTCLim)  
    vline(respWinTT,'k--')
    hline(0,'k:')
    figAxForm
    title(sprintf('Pref %s, n=%s, p=%s',num2str(orientations(iori)),...
        num2str(oriGroups.n(iori)),num2str(round(oriGroups.lateCycTest(iori),2,'significant'))))
    
end
print([fnout 'tuningAnalysis_oriGroupsFirstTC'],'-dpdf','-fillpage')  

setFigParams4Print('portrait')
figure;
suptitle('Tuned, First or Late Distractor Responsive Neurons, Short Cyc Expts')
for iori = 1:nOri
    subplot(4,2,iori)
    for iav = 1:2
        y = oriGroups.longTC{iav,iori}(tcStartFrame:end);
        yerr = oriGroups.longTCErr{iav,iori}(tcStartFrame:end);
        hold on
        shadedErrorBar_chooseColor(tt_longTC,y,yerr,cueColor{iav});
    end
    figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
    figYAxis([],'dF/F',oriTCLim)  
    vline(lateWinTT,'k--')
    hline(0,'k:')
    figAxForm([],0)
    title(sprintf('Pref %s, n=%s, p=%s',num2str(orientations(iori)),...
        num2str(oriGroups.nShortCyc(iori)),num2str(round(oriGroups.lateWinTest(iori),2,'significant'))))
    
%     subplot(4,2,iori)
%     for iav = 1:2
%         y = oriGroups.cycTC{iav,iori}(tcStartFrame:end);
%         yerr = oriGroups.cycTCErr{iav,iori}(tcStartFrame:end);
%         hold on
%         shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
%     end
%     figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
%     figYAxis([],'dF/F',oriTCLim)  
%     vline(respWinTT,'k--')
%     hline(0,'k:')
%     figAxForm
%     title(sprintf('Pref %s, n=%s, p=%s',num2str(orientations(iori)),...
%         num2str(oriGroups.n(iori)),num2str(round(oriGroups.lateCycTest(iori),2,'significant'))))
    
    subplot(4,2,iori+4)
    for itar = 1:2
        y = oriGroups.targetTC{iori}(tcStartFrame:end,itar);
        yerr = oriGroups.targetTCErr{iori}(tcStartFrame:end,itar);
        hold on
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,hiLoColor{itar});
    end
    figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',targetTCLim)  
    hline(0,'k:')
    vline(respWinTT,'k--')
    figAxForm
    title(sprintf('Pref %s, n=%s',num2str(orientations(iori)),...
        num2str(oriGroups.nTar(iori))))
    
end

print([fnout 'tuningAnalysis_oriGroupsTC'],'-dpdf','-fillpage')    

setFigParams4Print('portrait')
figure;
suptitle(sprintf('Tuned, Late Distractor-Responsive Neurons, %s/%s',...
    num2str(sum(cellInfo.isTuned & (cellInfo.lateCycRespCells))),...
    num2str(length(cellInfo.isTuned))))
xsub = [-0.25 +0.25];
subplot 231
for iav = 1:2
    hold on
    errorbar((1:nOri)+xsub(iav),oriGroups.firstResp(iav,:),...
        oriGroups.firstRespErr(iav,:),'.')
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_resp)
figAxForm
title('First Stim')

subplot 232
for iav = 1:2
    hold on
    errorbar((1:nOri)+xsub(iav),oriGroups.lateWin(iav,:),...
        oriGroups.lateWinErr(iav,:),'.')
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_win)
figAxForm
title('Late Window')

subplot 233
for iav = 1:2
    hold on
    errorbar((1:nOri)+xsub(iav),oriGroups.lateCycResp(iav,:),...
        oriGroups.lateCycRespErr(iav,:),'.')
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_resp)
figAxForm
title(sprintf('Late Stim Resp, One-Way ANOVA p=%s',num2str(round(...
    oriGroups.lateCycRespTest,2,'significant'))))

subplot 234
h = bar(oriGroups.n,'EdgeColor','none','BarWidth',0.5);
hold on
for iori = 1:nOri
    text(iori,oriGroups.n(iori)+1,num2str(oriGroups.n(iori)))
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriNLim)
figAxForm

subplot 235
h = bar(oriGroups.lateCycSI,'EdgeColor','none','BarWidth',0.5);
hold on
errorbar(1:nOri,oriGroups.lateCycSI,oriGroups.lateCycSIErr,'.')
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'Selectivity',siOriLim)
figAxForm
title(sprintf('Late Stim Resp, One-Way ANOVA p=%s',num2str(round(...
    oriGroups.lateCycSITest,2,'significant'))))

print([fnout 'tuningAnalysis_lateDistRespCells'],'-dpdf','-fillpage')

% population task tuning curves of ori pref groups
setFigParams4Print('landscape')    
figure
suptitle(sprintf('Tuned and Distractor or Target Resp. Neurons %s/%s',...
    num2str(sum(cellInfo.isTuned & ...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.targetRespCells))),...
    num2str(length(cellInfo.isTuned))))
for iori = 1:nOri
    subplot(2,2,iori)
    x = oriGroups.targetTuningStim{iori};
    xerr = oriGroups.targetTuningStimErr{iori};
    y = oriGroups.targetTuningResp{iori};
    yerr = oriGroups.targetTuningRespErr{iori};
    errorbar(0,oriGroups.firstStimRespForTargetAnalysis(iori),...
        oriGroups.firstStimRespErrForTargetAnalysis(iori),'.')
    hold on
    errorbar(x,y,yerr,yerr,xerr,xerr,'.')
    figXAxis([],'Task Ori. (deg)',[-10 100],x,oriBins([1,3,4]))
    figYAxis([],'dF/F',oriLim_taskResp)  
    vline(lateWinTT,'k--')
    hline(0,'k:')
    figAxForm([],0)
    title(sprintf('Pref %s, n=%s',num2str(orientations(iori)),...
        num2str(oriGroups.n(iori))))
end

print([fnout 'tuningAnalysis_oriGroupsTaskTuning'],'-dpdf','-fillpage')   

% orientation groups binned by SI
ind = cellInfo.lateCycRespCells & cellInfo.isTuned;
frOriInBin = nan(nSIBins,4);
frSIInOriBin = nan(4,nSIBins);
siInOriBinChiTest = nan(4,1);
% binnedSIxOri = nan(4,2);
for ibin = 1:nSIBins
    binInd = siBinID' == ibin & ind;
    frOriInBin(ibin,:) = histcounts(cellInfo.oriPref(binInd),4,...
        'Normalization','probability');
%     frOriInBin(ibin,:) = histcounts(cellInfo.oriPref(binInd),4);
end
for iori = 1:4
    binInd = ind & cellInfo.oriPref == iori;
    frSIInOriBin(iori,:) = histcounts(siBinID(binInd),1:(nSIBins+1),...
        'Normalization','probability');
%     binnedSIxOri{iori,1} = SI(binInd & SI' < 0);
    n = sum(binInd);
    n2 = sum(binInd' & siBinID == 1);
    [~,siInOriBinChiTest(iori)] = prop_test([n2, round(n/2)],[n, n], 0.05);
end
binnedHT = cell(4,nSIBins);
binnedET = cell(4,nSIBins);
binnedT = cell(4,nSIBins);
binnedD = cell(4,nSIBins);
binnedSI = cell(4,nSIBins);
for ibin = 1:nSIBins
    for iori = 1:4
        binInd = ind & siBinID' == ibin & cellInfo.oriPref == iori;
        binnedHT{iori,ibin} = Tr_vis{1}(binInd);
        binnedET{iori,ibin} = Tr_vis{2}(binInd);
        binnedT{iori,ibin} = Tr_vis{3}(binInd);
        binnedD{iori,ibin} = Dr{3}(binInd);
        binnedSI{iori,ibin} = SI(binInd);
    end
end

% matName = {'Stim';'SI'};
oriSIBinTest = nan(1,nSIBins);
for ibin = 1:nSIBins
    binInd = ind & siBinID' == ibin;
    oriSIBinTest(ibin)= anova1(SI(binInd),cellInfo.oriPref(binInd),'off');
end



figure
subplot 221
bar(frOriInBin,'stacked')
figXAxis([],'',[0 4],1:2,{'-SI','+SI'})
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
legend(cellfun(@num2str,num2cell(orientations),'unif',0),...
    'location','northeastoutside')
subplot 222
bar(frSIInOriBin,'stacked')
oriLabel = cellfun(@num2str,num2cell(orientations),'unif',0);
figXAxis([],'',[0 5],1:4,oriLabel)
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
legend({'-SI','+SI'},...
    'location','northeastoutside')
subplot 223
for ibin = 1:nSIBins
    y = cellfun(@mean,binnedSI(:,ibin));
    yerr = cellfun(@(x) ste(x,2),binnedSI(:,ibin));
    hold on
    errorbar(1:4,y,yerr,'.')
end
figXAxis([],'',[0 5],1:4,oriLabel)
figYAxis([],'SI',siBinLim)
figAxForm
legend({'-SI','+SI'},...
    'location','northeastoutside')

print([fnout 'tuningAnalysis_oriGroupsSIBinning'],'-dpdf','-fillpage') 

%% alternate selectivity index
resp_v = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,:),1);
resp_a = mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,:),1);
si_alt = (resp_v-resp_a)./(resp_v+resp_a);

ind = cellInfo.lateCycRespCells;
y = antiAnalysis.lateCycSI(ind);
x = si_alt(ind);

figure
subplot 231
plot(x,y,'.')
figXAxis([],'(V-A)./(V+A)',[-7 12])
figYAxis([],'(V-A)./stdpool(V,A)',[-7 12])
figAxForm
hline(0,'k:')
vline(0,'k:')
title(sprintf('r=%s',sigfigString(corr(x',y'))))
subplot 234
plot(x,y,'.')
figXAxis([],'(V-A)./(V+A)',[-2 2])
figYAxis([],'(V-A)./stdpool(V,A)',[-4 4])
figAxForm
hline(0,'k:')
vline(0,'k:')
title('zoomed in')
subplot 232
h1 = cdfplot(x);
vline(mean(x),'b')
figYAxis([],'fraction of cells',[0 1])
figXAxis([],'selectivity index',[-2 2])
figAxForm
title('(V-A)./(V+A)')
subplot 235
h2 = cdfplot(y);
vline(mean(y),'r')
figXAxis([],'selectivity index',siLim)
figYAxis([],'fraction of cells',[0 1])
figAxForm
title('(V-A)./stdpool(V,A)')

si_alt_ori = cell(1,nOri);
si_norm_ori = cell(1,nOri);
for iori = 1:nOri
    ind = cellInfo.lateCycRespCells & cellInfo.isTuned & cellInfo.oriPref == iori;
    si_alt_ori{iori} = si_alt(ind);
    si_norm_ori{iori} = antiAnalysis.lateCycSI(ind);
end

subplot 233
y = cellfun(@mean,si_alt_ori);
yerr = cellfun(@(x) ste(x,2),si_alt_ori);
errorbar(1:nOri,y,yerr,'.')
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'Selectivity',[])
figAxForm
title('SI = (V-A)./(V+A)')
subplot 236
y = cellfun(@mean,si_norm_ori);
yerr = cellfun(@(x) ste(x,2),si_norm_ori);
errorbar(1:nOri,y,yerr,'.')
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'Selectivity',[])
figAxForm
title('SI = (V-A)./stdpool(V,A)')
print([fnout 'alternateSelectivityAnalysis'],'-dpdf','-fillpage')
%% SI x tuning reliability
bins_rb = [0:5:30,100];
figure
subplot 121
ind = cellInfo.lateCycRespCells;
x = cellInfo.tuningReliability(ind);
y = antiAnalysis.lateCycSI(ind);
plot(x,y,'.')
[n,~,rb_binID] = histcounts(x,bins_rb);
rb = nan(1,length(bins_rb));
rb_ste = nan(1,length(bins_rb));
si_bin = nan(1,length(bins_rb));
si_bin_ste = nan(1,length(bins_rb));
for ibin = 1:length(bins_rb)-1
    ind = rb_binID == ibin;
    rb(ibin) = mean(x(ind),1);
    rb_ste(ibin) = ste(x(ind),1);
    si_bin(ibin) = mean(y(ind),2);
    si_bin_ste(ibin) = ste(y(ind),2);
end
hold on
errorbar(rb,si_bin,si_bin_ste,si_bin_ste,rb_ste,rb_ste,'.-')
figXAxis([],'Tuning Reliability',[0 90])
figYAxis([],'SI',siLim)
figAxForm

subplot 122
x2 = [mean(x(x<=tuningReliabilityThresh),1),mean(x(x>tuningReliabilityThresh),1)];
x2_err = [ste(x(x<=tuningReliabilityThresh),1),ste(x(x>tuningReliabilityThresh),1)];
y2 = [mean(y(x<=tuningReliabilityThresh),2),mean(y(x>tuningReliabilityThresh),2)];
y2_err = [ste(y(x<=tuningReliabilityThresh),2),ste(y(x>tuningReliabilityThresh),2)];
errorbar(x2,y2,y2_err,y2_err,x2_err,x2_err,'.-')
figXAxis([],'Tuning Reliability',[0 90])
figYAxis([],'SI',siLim)
figAxForm
print([fnout 'tuningAnalysis_tuningReliabilityBinning'],'-dpdf') 

%% SI x dF/F
figure
for iplot = 1:2
    if iplot == 1
        ind = cellInfo.lateCycRespCells;
    elseif iplot == 2
        ind = cellInfo.lateCycRespCells & cellInfo.isTuned;
    else
        error('index not added')
    end
    subplot(2,2,iplot)
    x = mean(antiAnalysis.lateCycTC{3}(respwin,ind));
    x(x<0) = nan;
    x = log(x);
    y = antiAnalysis.lateCycSI(ind);
    mdl = fitlm(x,y);
    yfit = predict(mdl,sort(x)');
    [corrmat,pmat] = corrcoef(x(~isnan(x)),y(~isnan(x)));
    hold on
%     plot(x,y,'.');
    plot(sort(x),yfit,'-')
    if iplot==1
        for iori = 1:5
            if iori==5
                ind_ori = ind & ~cellInfo.isTuned;
            else
                ind_ori = ind & cellInfo.isTuned & cellInfo.oriPref==iori;
            end
            x = mean(antiAnalysis.lateCycTC{3}(respwin,ind_ori));
            x(x<0) = nan;
            x = log(x);
            y = antiAnalysis.lateCycSI(ind_ori);
            plot(x,y,'.');
        end
        legend(cat(2,{'fit'},strsplit(num2str(orientations)),{'untuned'}),...
            'location','northeastoutside')
            
    else
        plot(x,y,'.');
    end
    % ax = gca;
    % ax.XScale = 'log';
    figXAxis([],'Late Stim Resp, (log(dF/F))',[])
    figYAxis([],'SI',siLim)
    figAxForm
    title(sprintf('tuned, late resp. cells n=%s, r=%s, p=%s',num2str(sum(ind)),...
        num2str(round(corrmat(1,2),2,'significant')),...
        num2str(round(pmat(1,2),2,'significant'))))
    if iplot==1
        subplot(2,2,iplot+2)
        errorbar(1:4,oriGroups.matchSIxFR,oriGroups.matchSIxFRErr,'.')
        figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
        figYAxis([],'Selectivity',siOriLim)
        figAxForm
        hline(0,'k:')
        title(sprintf('SI matched for FR, top %s neurons in ea. grp.',...
            num2str(minCellN_SIFRmatch)))
    end
end

print([fnout 'tuningAnalysis_oriGroupsMatchedSIxFR'],'-dpdf','-fillpage')  

%% noise correlations
load(fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs','bxStats.mat'))
ind = ismember(bxStats.mouseNames,mice);
[~, bxDataSortInd] = sort(bxStats.mouseNames(ind));
hrDiffAll = bxStats.av(visualTrials).matchedHRDiff(ind);
hrDiff = hrDiffAll(bxDataSortInd);
% example expt
if strcmp(ds,'FSAV_attentionV1')
    iexp = 5;
else
    iexp = 8;
end
ind = respCellsExpt(iexp).lateCycRespCells==1;
rsc = respCellsExpt(iexp).rsc{3}(ind,ind);
rsc(eye(sum(ind))==1) = nan;
[~,sortInd] = sort(nanmean(rsc,1));

figure
colormap(brewermap([],'*RdBu'))
imagesc(rsc(sortInd,sortInd))
colorbar
clim([-0.5 0.5])
figXAxis([],'Cell #',[])
figYAxis([],'Cell #',[])
figAxForm
title(sprintf('Noise Corr. Across All trials, sorted by mean noise corr (%s)',...
    antiDataExpt(iexp).exptName))
print([fnout 'corrAnalysis_exExpt'],'-dpdf','-fillpage') 

figure
ind = cellInfo.lateCycRespCells;
suptitle(sprintf('Late Resp. Neurons, n=%s',num2str(sum(ind))))
colormap(brewermap([],'*RdBu'));
for iav = 1:2
    subplot(4,2,iav)
    imagesc(antiAnalysis.rsc{iav}(ind,ind));
    caxis([-1 1])
    title(sprintf('%s Trials',avName{iav}))
    colorbar
    figXAxis([],'Cell #',[])
    figYAxis([],'Cell #',[])
    figAxForm
end

% rsc_ind = cellfun(@(x) tril(x(ind,ind),antiAnalysis.rsc,'unif',0);
% cellfun(@(x,y) splitapply(@mean,x,findgroups(discretize(
% gm = antiAnalysis.gm{3}(ind);

subplot 423
for iav = 1:2    
    rsc = antiAnalysis.rsc{iav}(ind,ind);
    halfInd = tril(true(sum(ind)),-1) & ~isnan(rsc);
    rsc = rsc(halfInd);
    gm = antiAnalysis.gm{iav}(ind,ind);
    gm = gm(halfInd);
    
    binEdges = [min(gm),quantile(gm,3),max(gm)];
    [~,~,binID] = histcounts(gm,binEdges);
    gm_binned = splitapply(@mean,gm,binID);
    rsc_binned = splitapply(@mean,rsc,binID);
    gm_binned_err = splitapply(@(x) ste(x,1),gm,binID);
    rsc_binned_err = splitapply(@(x) ste(x,1),rsc,binID);
%     subplot(2,2,iav+2)
%     plot(gm,rsc,'.')
    hold on
    errorbar(gm_binned,rsc_binned,rsc_binned_err,rsc_binned_err,...
        gm_binned_err,gm_binned_err,'.-')
    
end
figXAxis([],'Geometric Mean dF/F (All Trials)',[0 0.02])
%     ax = gca;
%     ax.XScale = 'log';
figYAxis([],'Mean Noise Corr.',[0 0.2])
figAxForm
legend(avName,'location','northeastoutside')

subplot 426
rsc = antiAnalysis.rsc{visualTrials}(ind,ind);
halfInd = tril(true(sum(ind)),-1) & ~isnan(rsc);
x = rsc(halfInd);
rsc = antiAnalysis.rsc{auditoryTrials}(ind,ind);
halfInd = tril(true(sum(ind)),-1) & ~isnan(rsc);
y = rsc(halfInd);
plot(x,y,'.')
hold on
figXAxis([],'Visual Noise Corr.',[-1 1])
figYAxis([],'Auditory Noise Corr.',[-1 1])
figAxForm
plot([-1 1],[-1 1],'k:')
[~,p] = ttest(x,y);
errorbar(mean(x),mean(y),ste(y,1),ste(y,1),ste(x,1),ste(x,1),'.','MarkerSize',20)
title(sprintf('Vis:%s, Aud:%s, p=%s',num2str(round(mean(x),2,'significant')),...
    num2str(round(mean(y),2,'significant')),num2str(round(p,2,'significant'))))

subplot 425
rsc = antiAnalysis.rsc{3}(ind,ind);
halfInd = tril(true(sum(ind)),-1) & ~isnan(rsc);
x = rsc(halfInd);
cdfplot(x)
figXAxis([],'Noise Corr (All Trials)',[-.2 .4])
figYAxis([],'Fraction of Pairs',[0 1])
figAxForm
title(sprintf('All Trials, mean=%s,ste=%s',num2str(round(mean(x),2,'significant')),...
    num2str(round(ste(x,1),2,'significant'))))

rsc_ori = nan(2,nOri);
rsc_ori_err = nan(2,nOri);
for iav = 1:2
    for iori = 1:nOri
        ind_ori = cellInfo.lateCycRespCells & cellInfo.oriPref == iori & cellInfo.isTuned;
        rsc = antiAnalysis.rsc{iav}(ind_ori,ind_ori);
        halfInd = tril(true(sum(ind_ori)),-1) & ~isnan(rsc);
        x = rsc(halfInd);
        
        rsc_ori(iav,iori) = mean(x);
        rsc_ori_err(iav,iori) = ste(x,1);
        
    end
end

subplot 427
hold on
for iav = 1:2
    errorbar(1:nOri,rsc_ori(iav,:),rsc_ori_err(iav,:),'.','MarkerSize',20)
end
figXAxis([],'Ori Pref',[0 nOri+1],1:nOri,cellfun(@num2str,num2cell(orientations),'unif',0))
figYAxis([],'Mean Noise Corr.',[0 0.2])
figAxForm

subplot 428
ind = cellInfo.lateCycRespCells & cellInfo.isTuned;
oriPref = cellInfo.oriPref(ind);
samePrefInd = (oriPref-oriPref')==0;
hold on
for iav = 1:2    
    rsc = antiAnalysis.rsc{iav}(ind,ind);
    halfInd = tril(true(sum(ind)),-1) & ~isnan(rsc);
    
    x = rsc(halfInd & samePrefInd);
    y = rsc(halfInd & ~samePrefInd);
    errorbar(1:2,[mean(x),mean(y)],[ste(x,1),ste(y,1)],'.-','MarkerSize',20)
end
figXAxis([],'Ori Prefs of Pair',[0 3],1:2,{'Same','Diff'})
figYAxis([],'Mean Noise Corr',[0 0.1])
figAxForm
title('tuned neurons')

subplot 424
ind = cellInfo.lateCycRespCells & cellInfo.isTuned;
oriPref = cellInfo.oriPref(ind);
samePrefInd = (oriPref-oriPref')==0;
hold on
for iav = 1:2    
    rsc = antiAnalysis.rsc{iav}(ind,ind);
    halfInd = tril(true(sum(ind)),-1) & ~isnan(rsc);
    rsc = rsc(halfInd & samePrefInd);
    gm = antiAnalysis.gm{iav}(ind,ind);
    gm = gm(halfInd & samePrefInd);
    
    binEdges = [min(gm),quantile(gm,3),max(gm)];
    [~,~,binID] = histcounts(gm,binEdges);
    gm_binned = splitapply(@mean,gm,binID);
    rsc_binned = splitapply(@mean,rsc,binID);
    gm_binned_err = splitapply(@(x) ste(x,1),gm,binID);
    rsc_binned_err = splitapply(@(x) ste(x,1),rsc,binID);
    hold on
    errorbar(gm_binned,rsc_binned,rsc_binned_err,rsc_binned_err,...
        gm_binned_err,gm_binned_err,'.-')
    
    rsc = antiAnalysis.rsc{iav}(ind,ind);
    halfInd = tril(true(sum(ind)),-1) & ~isnan(rsc);
    rsc = rsc(halfInd & ~samePrefInd);
    gm = antiAnalysis.gm{iav}(ind,ind);
    gm = gm(halfInd & ~samePrefInd);
    
    binEdges = [min(gm),quantile(gm,3),max(gm)];
    [~,~,binID] = histcounts(gm,binEdges);
    gm_binned = splitapply(@mean,gm,binID);
    rsc_binned = splitapply(@mean,rsc,binID);
    gm_binned_err = splitapply(@(x) ste(x,1),gm,binID);
    rsc_binned_err = splitapply(@(x) ste(x,1),rsc,binID);
    hold on
    errorbar(gm_binned,rsc_binned,rsc_binned_err,rsc_binned_err,...
        gm_binned_err,gm_binned_err,'.-')
end
figXAxis([],'Geometric Mean dF/F (All Trials)',[0 0.02])
%     ax = gca;
%     ax.XScale = 'log';
figYAxis([],'Mean Noise Corr.',[0 0.2])
figAxForm
legend({'Vis-Same','Vis-Dif','Aud-Same','Aud-Diff'},'location','northeastoutside')

print([fnout 'corrAnalysis'],'-dpdf','-fillpage')

figure
subplot 221
if strcmp(ds,'FSAV_attentionV1')
    colors = brewermap(nMice,'Set2');
elseif strcmp(ds,'FSAV_attentionV1_noAttn')
    colors = brewermap(nMice,'Dark2');
else
    colors = brewermap(nMice,'Set1');
end
rsc_exp = nan(2,nexp);
for iexp = 1:nexp
    msInd = strcmp(mice,{respCellsExpt(iexp).mouse});
    ind = respCellsExpt(iexp).lateCycRespCells==1;
    rsc = respCellsExpt(iexp).rsc{visualTrials}(ind,ind);
    halfInd = tril(true(sum(ind)),-1);
    x = rsc(halfInd);
    rsc = respCellsExpt(iexp).rsc{auditoryTrials}(ind,ind);
    halfInd = tril(true(sum(ind)),-1);
    y = rsc(halfInd);
    hold on
    h = plot(mean(x),mean(y),'.','MarkerSize',20);
    h.Color = colors(msInd,:);
    rsc_exp(visualTrials,iexp) = mean(x);
    rsc_exp(auditoryTrials,iexp) = mean(y);
end
plot([-1 1],[-1 1],'k:')
figXAxis([],'Visual Noise Corr.',[-0.1 0.25])
figYAxis([],'Auditory Noise Corr.',[-0.1 0.25])
figAxForm
rsc_eaMs = nan(2,nMice);
for im = 1:nMice
    for iav = 1:2
        msInd = find(strcmp({respCellsExpt.mouse},mice(im)));
        rsc = [];
        for i = 1:length(msInd)
            ind = respCellsExpt(msInd(i)).lateCycRespCells==1;
            rsc_temp = respCellsExpt(msInd(i)).rsc{iav}(ind,ind);
            halfInd = tril(true(sum(ind)),-1) & ~isnan(rsc_temp);
            rsc = cat(1,rsc,rsc_temp(halfInd));
        end
        rsc_eaMs(iav,im) = mean(rsc);
    end
end
subplot 222
hold on
for im = 1:nMice
    h = plot(hrDiff(im),rsc_eaMs(visualTrials,im) ...
        - rsc_eaMs(auditoryTrials,im),'.','MarkerSize',20);
    h.Color = colors(im,:);
end
figXAxis([],'HR Diff (Val-Inv)',[-0.2 0.5])
figYAxis([],'Noise Corr. (Vis-Aud)',[-0.05 0.02])
figAxForm
subplot 223
hold on
for im = 1:nMice
    h = plot(hrDiff(im),rsc_eaMs(visualTrials,im),'.','MarkerSize',20);
    h.Color = colors(im,:);
end
figXAxis([],'HR Diff (Val-Inv)',[-0.2 0.5])
figYAxis([],'Vis Noise Corr.',[-.1 0.2])
figAxForm
subplot 224
hold on
for im = 1:nMice
    h = plot(hrDiff(im),rsc_eaMs(auditoryTrials,im),'.','MarkerSize',20);
    h.Color = colors(im,:);
end
figXAxis([],'HR Diff (Val-Inv)',[-0.2 0.5])
figYAxis([],'Aud Noise Corr.',[-.1 0.2])
figAxForm

print([fnout 'corrAnalysisXmouse'],'-dpdf')

%%
imgStats.nCells.oriGroups = oriGroups.n;
fprintf('n 0/45/90/135 pref = %s\n', num2str(imgStats.nCells.oriGroups))

imgStats.lateCycRespAllxOri = oriGroups.lateCycRespAll;
imgStats.lateCycRespAllxOriErr = oriGroups.lateCycRespAllErr;
% fprintf('Late Cyc Resp (All Trials) 0/45/90/135 pref = %s\n', ...
%     num2str(round(imgStats.lateCycRespAllxOri,2,'significant')))
% fprintf('Late Cyc Resp Err (All Trials) 0/45/90/135 pref = %s\n', ...
%     num2str(round(imgStats.lateCycRespAllxOriErr,2,'significant')))

imgStats.lateCycRespDiffxOri = oriGroups.lateCycRespDiff;
imgStats.lateCycRespDiffxOriErr = oriGroups.lateCycRespDiffErr;
imgStats.lateCycRespTestxOri = oriGroups.lateCycTest;
imgStats.lateCycSITestxOri = oriGroups.lateCycSITestEaOri;
imgStats.lateCycSImatchedFRTestxOri = oriGroups.FRmatchedSITestEaOri;
% fprintf('Late Cyc Resp (Vis-Aud Trials) 0/45/90/135 pref = %s\n', ...
%     num2str(round(imgStats.lateCycRespDiffxOri,2,'significant')))
% fprintf('Late Cyc Resp Err (Vis-Aud Trials) 0/45/90/135 pref = %s\n', ...
%     num2str(round(imgStats.lateCycRespDiffxOriErr,2,'significant')))
fprintf('Late Cyc Resp AV Test x Ori, 0/45/90/135 pref, p=%s\n',num2str(round(...
    imgStats.lateCycRespTestxOri,2,'significant')))
fprintf('Late Cyc SI Test x Ori, 0/45/90/135 pref, a=%s, p=%s\n',...
    num2str(round(0.05/(nOri-1),2,'significant')),num2str(round(...
    imgStats.lateCycSITestxOri,2,'significant')))
fprintf('Late Cyc FR Matched SI Test x Ori, 0/45/90/135 pref, a=%s, p=%s\n',...
    num2str(round(0.05/(nOri-1),2,'significant')),num2str(round(...
    imgStats.lateCycSImatchedFRTestxOri,2,'significant')))
imgStats.lateCycSIxOriANOVA = oriGroups.lateCycSITest;
fprintf('SI x Ori 1-Way ANOVA, p=%s\n',num2str(round(...
    imgStats.lateCycSIxOriANOVA,2,'significant')))
imgStats.lateCycSIxOriMatchFRTest = oriGroups.matchSIxFRTest;
fprintf('SI x Ori Match FR, ANOVA p=%s\n',num2str(round(...
    imgStats.lateCycSIxOriMatchFRTest,2,'significant')))

%% suppressed cells analysis (Supplemental Figure 5)
setFigParams4Print('portrait')
figure
subplot 311
ind = cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells);
for iav = 1:2
    y = mean(antiAnalysis.longTC{iav}...
        ((tcStartFrame:end),(ind) & ...
        cellInfo.isShortCycExpt),2);
    yerr = ste(antiAnalysis.longTC{iav}...
        ((tcStartFrame:end),(ind) & cellInfo.isShortCycExpt),2);
    hold on
    shadedErrorBar_chooseColor(tt_longTC,y,yerr,cueColor{iav});
end
figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
figYAxis([],'dF/F',suppTCLim)  
vline(lateWinTT,'k--')
hline(0,'k:')
figAxForm([],0)
title(sprintf('Suppressed Only Cells (%s/%s)',...
    num2str(sum((...
    ind)...
    & cellInfo.isShortCycExpt)),...
    num2str(sum(cellInfo.isShortCycExpt))))

subplot 323
x = mean(antiAnalysis.longTC{visualTrials}(lateWinFr,...
    (ind)),1);
xerr = ste(x,2);
y = mean(antiAnalysis.longTC{auditoryTrials}(lateWinFr,...
    (ind)),1);
yerr = ste(y,2);
plot(x,y,'.')
hold on
errorbar(mean(x),mean(y),yerr,yerr,xerr,xerr,'.')
plot(suppScatLim_win,suppScatLim_win,'k--')
plot(suppScatLim_win,[0 0],'k:')
plot([0 0],suppScatLim_win,'k:')
[~,p] = ttest(x,y);
figXAxis([],'Visual (dF/F)',suppScatLim_win)
figYAxis([],'Auditory (dF/F)',suppScatLim_win)
figAxForm
title(sprintf('Late Window, All Resp. Cells (%s/%s), p = %s',...
    num2str(sum(ind)),...
    num2str(length(cellInfo.firstRespCells)),...
    num2str(round(p,2,'significant'))))

subplot 324
y = antiAnalysis.lateWinSI(cellInfo.lateSuppCells & ~cellInfo.lateCycRespCells);
[~,p]=ttest(y);
h = cdfplot(y);
hold on;
vline(mean(y),'k-')
figXAxis([],'Selectivity Index',siLim)
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
title(sprintf('Late Cyc. Resp. Cells, mean = %s+/-%s, p=%s',...
    num2str(round(mean(y),2,'significant')),...
    num2str(round(ste(y,2),2,'significant')),...
    num2str(round(p,2,'significant'))))

subplot 325
ind1 = cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells) & ...
    cellInfo.lateStimAuROCTest & cellInfo.lateStimAuROC > 0.5;
ind2 = cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells) & ...
    cellInfo.lateStimAuROCTest & cellInfo.lateStimAuROC < 0.5;
ind3 = cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells) & ...
    ~cellInfo.lateStimAuROCTest;

y = antiAnalysis.lateWinSI(ind1);
yerr = ste(y,2);
errorbar(1,mean(y),yerr,'.')
[~,p]=ttest(y);
text(1,0.5,num2str(round(p,2,'significant')))
hold on
y = antiAnalysis.lateWinSI(ind2);
yerr = ste(y,2);
errorbar(2,mean(y),yerr,'.')
[~,p]=ttest(y);
text(2,0.5,num2str(round(p,2,'significant')))
y = antiAnalysis.lateWinSI(ind3);
yerr = ste(y,2);
errorbar(3,mean(y),yerr,'.')
[~,p]=ttest(y);
text(3,0.5,num2str(round(p,2,'significant')))
figXAxis([],'Vis. Task Stim. Pref.',[0 4],1:3,{'T','D','NP'})
figYAxis([],'Selectivity (Late Win)',[-1 0.5])
figAxForm
hline(0,'k:')
legend({num2str(sum(ind1)),num2str(sum(ind2)),num2str(sum(ind3))})

subplot 326
ind1 = cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells) & ...
    cellInfo.audLateStimAuROCTest & cellInfo.audLateStimAuROC > 0.5;
ind2 = cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells) & ...
    cellInfo.audLateStimAuROCTest & cellInfo.audLateStimAuROC < 0.5;
ind3 = cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells) & ...
    ~cellInfo.audLateStimAuROCTest;

y = antiAnalysis.lateWinSI(ind1);
yerr = ste(y,2);
errorbar(1,mean(y),yerr,'.')
[~,p]=ttest(y);
text(1,0.5,num2str(round(p,2,'significant')))
hold on
y = antiAnalysis.lateWinSI(ind2);
yerr = ste(y,2);
errorbar(2,mean(y),yerr,'.')
[~,p]=ttest(y);
text(2,0.5,num2str(round(p,2,'significant')))
y = antiAnalysis.lateWinSI(ind3);
yerr = ste(y,2);
errorbar(3,mean(y),yerr,'.')
[~,p]=ttest(y);
text(3,0.5,num2str(round(p,2,'significant')))
figXAxis([],'Aud. Task Stim. Pref.',[0 4],1:3,{'T','D','NP'})
figYAxis([],'Selectivity (Late Win)',[-1 0.5])
figAxForm
hline(0,'k:')
legend({num2str(sum(ind1)),num2str(sum(ind2)),num2str(sum(ind3))})

print([fnout 'tcLongTrialsAndLateCycWithQuant_SuppCells'],'-dpdf','-fillpage')
%% Supplemental Figure - Fraction Mod Cells
ind = antiAnalysis.lateCycAV95CITest == 1;
pctLateAVMod_UD = ([sum(ind' & antiAnalysis.lateCycSI' <= 0 & cellInfo.lateCycRespCells),...
    sum(ind' & antiAnalysis.lateCycSI' > 0 & cellInfo.lateCycRespCells)])...
    ./sum(cellInfo.lateCycRespCells);
pctSuppAVMod_UD = ([sum(ind' & antiAnalysis.lateWinSI' <= 0 &  cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells)),...
    sum(ind' & antiAnalysis.lateWinSI' > 0 & cellInfo.lateSuppCells & ~cellInfo.lateCycRespCells)])...
    ./sum(cellInfo.lateSuppCells & ~cellInfo.lateCycRespCells);

figure
subplot 211
bar(cat(1,pctLateAVMod_UD,pctSuppAVMod_UD),'stacked');
figXAxis([],'Cell Type',[0 3],1:2,{'Late Resp.';'Late Supp.'})
figYAxis([],'Fraction V-A Modulated Cells',[0 1])
figAxForm
legend({'V<A';'V>A'},'location','northeastoutside')
title('Significance tested by 95% CI of bootstrap V-A S.I.')

ind = antiAnalysis.lateCycAVShuffTest == 1;
pctLateAVMod_UD = ([sum(ind' & antiAnalysis.lateCycSI' <= 0 & cellInfo.lateCycRespCells),...
    sum(ind' & antiAnalysis.lateCycSI' > 0 & cellInfo.lateCycRespCells)])...
    ./sum(cellInfo.lateCycRespCells);
pctSuppAVMod_UD = ([sum(ind' & antiAnalysis.lateWinSI' <= 0 &  cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells)),...
    sum(ind' & antiAnalysis.lateWinSI' > 0 & cellInfo.lateSuppCells & ~cellInfo.lateCycRespCells)])...
    ./sum(cellInfo.lateSuppCells & ~cellInfo.lateCycRespCells);

subplot 212
bar(cat(1,pctLateAVMod_UD,pctSuppAVMod_UD),'stacked');
figXAxis([],'Cell Type',[0 3],1:2,{'Late Resp.';'Late Supp.'})
figYAxis([],'Fraction V-A Modulated Cells',[0 1])
figAxForm
legend({'V<A';'V>A'},'location','northeastoutside')
title('Significant SI - Shuffled Trial ID')

print([fnout 'fractionVAModCells'],'-dpdf')

%%
ind = cellInfo.lateSuppCells & ~cellInfo.lateCycRespCells;
for iav = 1:2
    imgStats.av(iav).lateWinSupp = mean(mean(antiAnalysis.longTC{iav}(lateWinFr,...
        (ind)),1));
    imgStats.av(iav).lateWinSuppErr = ste(mean(antiAnalysis.longTC{iav}(lateWinFr,...
        (ind)),1),2);
    fprintf('%s Suppr. Late Win mean/err: %s/%s\n',avName{iav},...
        num2str(round(imgStats.av(iav).lateWinSupp,2,'significant')),...
        num2str(round(imgStats.av(iav).lateWinSuppErr,2,'significant')))
    
    imgStats.av(iav).lateCycSupp = mean(mean(antiAnalysis.lateCycTC{iav}(respwin,...
        (ind)),1));
    imgStats.av(iav).lateCycSuppErr = ste(mean(antiAnalysis.lateCycTC{iav}(respwin,...
        (ind)),1),2);
    fprintf('%s Suppr. Late Cyc Resp mean/err: %s/%s\n',avName{iav},...
        num2str(round(imgStats.av(iav).lateCycSupp,2,'significant')),...
        num2str(round(imgStats.av(iav).lateCycSuppErr,2,'significant')))
    
%     x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,...
%         (cellInfo.lateCycRespCells)),1);
%     y = mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,...
%         (cellInfo.lateCycRespCells)),1);
end
[~,imgStats.lateWinTestSuppCells] = ttest(...
    mean(antiAnalysis.longTC{visualTrials}(lateWinFr,ind),1),...
    mean(antiAnalysis.longTC{auditoryTrials}(lateWinFr,ind),1));
fprintf('Late Win Test, Suppr. Cells, p=%s\n',...
    num2str(round(imgStats.lateWinTestSuppCells,2,'significant')))
[~,imgStats.lateCycTestSuppCells] = ttest(...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1),...
    mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,ind),1));
fprintf('Late Cyc Resp Test, Suppr. Cells, p=%s\n',...
    num2str(round(imgStats.lateCycTestSuppCells,2,'significant')))

ind = antiAnalysis.lateCycAV95CITest == 1;
imgStats.modCells.nLateResp = [sum(ind' & antiAnalysis.lateCycSI' <= 0 & cellInfo.lateCycRespCells),...
    sum(ind' & antiAnalysis.lateCycSI' > 0 & cellInfo.lateCycRespCells)];
imgStats.modCells.nSupp = [sum(ind' & antiAnalysis.lateWinSI' <= 0 &  cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells)),...
    sum(ind' & antiAnalysis.lateWinSI' > 0 & cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells))];

fprintf('Fraction modulated late resp cells: %s\n',...
    num2str(round(sum(imgStats.modCells.nLateResp)./imgStats.nCells.lateCycResp,3,'significant')))
fprintf('Fraction modulated late supp cells: %s\n',...
    num2str(round(sum(imgStats.modCells.nSupp)./imgStats.nCells.lateSupp,3,'significant')))
fprintf('N +SI Mod. Late Resp Cells: %s/%s\n',num2str(imgStats.modCells.nLateResp(2)),...
    num2str(sum(imgStats.modCells.nLateResp)))
fprintf('N +SI Mod. Late Supp Cells: %s/%s\n',num2str(imgStats.modCells.nSupp(2)),...
    num2str(sum(imgStats.modCells.nSupp)))

save([fnout 'imgStats'],'imgStats')

%% trial outcome target and distractor time-courses

figure
ind = cellInfo.lateCycRespCells==1;
for iav = 1:2
    subplot(3,2,iav)
    L = [];
    r = nan(sum(ind),2);
    for iout = 1:2
        outType = iout+2;
        x = targetAnalysis.trOutTC{outType,iav}(:,ind);
        hold on
        h = shadedErrorBar_chooseColor(tt_cycTC,...
            mean(x(tcStartFrame:end,:)-mean(x(basewin_0,:),1),2)',...
            ste(x(tcStartFrame:end,:),2),trOutColor{iout});
        L(iout) = h.mainLine;
        r(:,iout) = mean(x(respwin,:),1)-mean(x(basewin_0,:),1);
    end
    [~,p] = ttest(r(:,1),r(:,2));
    figXAxis([],'Time (ms)',[ttLabel_cyc(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',[-0.005 0.06]) 
    hline(0,'k:')
    vline(respWinTT,'k--')
    figAxForm
    legend(L,trOutType(iout),'location','northwest')
    title(sprintf('%s Dist Resp Cells, n=%s,p=%s',avName{iav},...
        num2str(sum(ind)),num2str(round(p,2,'significant'))))
    if iav == 1
        vline(visRTwindow(1),'b--')
    else
        vline(audRTwindow(1),'b--')
    end
end
for iav = 1:2
    subplot(3,2,iav+2)
    if iav == 1
        ind = cellInfo.lateCycRespCells | cellInfo.targetRespCells;
    else
        ind = cellInfo.lateCycRespCells | cellInfo.targetRespCells;
    end
    L = [];
    r = nan(sum(ind),2);
    for iout = 1:2
        x = targetAnalysis.trOutTC{iout,iav}(:,ind);
        hold on
        h = shadedErrorBar_chooseColor(tt_targetTC(tcStartFrame:end),...
            nanmean(x(tcStartFrame:end,:)-mean(x(basewin_0_target,:),1),2),...
            ste(x(tcStartFrame:end,:),2),trOutColor{iout});
        L(iout) = h.mainLine;
        r(:,iout) = mean(x(respwin_target,:),1)-mean(x(basewin_0_target,:),1);
    end
    [~,p] = ttest(r(:,1),r(:,2));
    figXAxis([],'Time (ms)',[ttLabel_cyc(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',[-0.005 0.06]) 
    hline(0,'k:')
    vline(respWinTT_target,'k--')
    figAxForm
    legend(L,trOutType(1:2),'location','northwest')
    if iav == 1
        title(sprintf('%s Target or Dist Resp Cells, n=%s,p=%s',avName{iav},...
            num2str(sum(ind)),num2str(round(p,2,'significant'))))
        vline(visRTwindow(1),'b--')
    else
        title(sprintf('%s Dist Resp Cells, n=%s,p=%s',avName{iav},...
            num2str(sum(ind)),num2str(round(p,2,'significant'))))
        vline(audRTwindow(1),'b--')
    end
    
    subplot(3,2,iav+4)
    hold on
    ind = cellInfo.lateCycRespCells | cellInfo.targetRespCells;
    L = [];
    r = nan(sum(ind),2);
    for iout = 1:2
        x = targetAnalysis.trOutTC{iout+4,iav}(:,ind);
        h = shadedErrorBar_chooseColor(tt_targetTC(tcStartFrame:end),...
            nanmean(x(tcStartFrame:end,:)-mean(x(basewin_0_target,:),1),2),...
            ste(x(tcStartFrame:end,:),2),trOutColor{iout});
        L(iout) = h.mainLine;
        r(:,iout) = mean(x(respwin_target,:),1)-mean(x(basewin_0_target,:),1);
    end
    [~,p] = ttest(r(:,1),r(:,2));
    figXAxis([],'Time (ms)',[ttLabel_cyc(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',[-0.005 0.06]) 
    hline(0,'k:')
    vline(respWinTT_target,'k--')
    figAxForm
    legend(L,trOutType(5:6),'location','northwest')
    title(sprintf('%s Target or Dist Resp Cells, n=%s,p=%s',avName{iav},...
        num2str(sum(ind)),num2str(round(p,2,'significant'))))
    vline(visRTwindow(1),'b--')
end

print([fnout 'choiceTC'],'-dpdf','-fillpage')

figure
suptitle('Target Response')
for iav = 1:2
    if iav == 1
        ind = cellInfo.lateCycRespCells | cellInfo.targetRespCells;
    else
        ind = cellInfo.lateCycRespCells | cellInfo.targetRespCells;
    end
    subplot(2,2,iav+1)
    L = [];
    hold on
    for itar = 1:2
        if iav == 2 & itar == 1
            x = targetAnalysis.tc{3,iav}(:,ind);
            h = shadedErrorBar_chooseColor(tt_targetTC(tcStartFrame:end),...
                nanmean(x(tcStartFrame:end,:)-mean(x(basewin_0_target,:),1),2),...
                ste(x(tcStartFrame:end,:),2),hiLoColor{itar});
            L = cat(2,L,h.mainLine);
        elseif iav == 2 & itar == 2
            continue
        else
            x = targetAnalysis.tc{itar,iav}(:,ind);
            h = shadedErrorBar_chooseColor(tt_targetTC(tcStartFrame:end),...
                nanmean(x(tcStartFrame:end,:)-mean(x(basewin_0_target,:),1),2),...
                ste(x(tcStartFrame:end,:),2),hiLoColor{itar});
            L = cat(2,L,h.mainLine);
        end
    end
    x = circshift(antiAnalysis.lateCycTC{iav}(:,ind),-1);
    h = shadedErrorBar_chooseColor(tt_targetTC(tcStartFrame:end),...
        nanmean(x(tcStartFrame:end,:)-mean(x(basewin_0_target,:),1),2),...
        ste(x(tcStartFrame:end,:),2),[0.75 0.75 0.75]);
    L = cat(2,L,h.mainLine);
    figXAxis([],'Time (ms)',[ttLabel_cyc(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',[-0.005 0.06]) 
    hline(0,'k:')
    vline(respWinTT_target,'k--')
    figAxForm
    if iav == 1
        legend(L,{'HT','ET','D'},'location','northwest')
        title(sprintf('%s, Target or Dist Resp Cells, n=%s,',avName{iav},...
            num2str(sum(ind))))
        vline(visRTwindow(1),'b--')
        
        ind = cellInfo.lateCycRespCells;
        subplot(2,2,iav)
        L = [];
        for itar = 1:2
            hold on
            x = targetAnalysis.tc{itar,iav}(:,ind);
            h = shadedErrorBar_chooseColor(tt_targetTC(tcStartFrame:end),...
                nanmean(x(tcStartFrame:end,:)-mean(x(basewin_0_target,:),1),2),...
                ste(x(tcStartFrame:end,:),2),hiLoColor{itar});
            L = cat(2,L,h.mainLine);
        end
        x = circshift(antiAnalysis.lateCycTC{iav}(:,ind),-1);
        h = shadedErrorBar_chooseColor(tt_targetTC(tcStartFrame:end),...
            nanmean(x(tcStartFrame:end,:)-mean(x(basewin_0_target,:),1),2),...
            ste(x(tcStartFrame:end,:),2),[0.75 0.75 0.75]);
        L = cat(2,L,h.mainLine);
        figXAxis([],'Time (ms)',[ttLabel_cyc(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
        figYAxis([],'dF/F',[-0.005 0.06]) 
        hline(0,'k:')
        vline(respWinTT_target,'k--')
        figAxForm
        legend(L,{'HT','ET','D'},'location','northwest')
        title(sprintf('%s, Dist Resp Cells, n=%s,',avName{iav},...
            num2str(sum(ind))))
        vline(visRTwindow(1),'b--')
        
        subplot(2,2,4)
        ind = cellInfo.lateCycRespCells| cellInfo.targetRespCells;
        L = [];
        hold on
        x = targetAnalysis.tc{3,iav}(:,ind);
        h = shadedErrorBar_chooseColor(tt_targetTC(tcStartFrame:end),...
            nanmean(x(tcStartFrame:end,:)-mean(x(basewin_0_target,:),1),2),...
            ste(x(tcStartFrame:end,:),2),hiLoColor{itar});
        L = cat(2,L,h.mainLine);
        x = circshift(antiAnalysis.lateCycTC{iav}(:,ind),-1);
        h = shadedErrorBar_chooseColor(tt_targetTC(tcStartFrame:end),...
            nanmean(x(tcStartFrame:end,:)-mean(x(basewin_0_target,:),1),2),...
            ste(x(tcStartFrame:end,:),2),[0.75 0.75 0.75]);
        L = cat(2,L,h.mainLine);
        figXAxis([],'Time (ms)',[ttLabel_cyc(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
        figYAxis([],'dF/F',[-0.005 0.06]) 
        hline(0,'k:')
        vline(respWinTT_target,'k--')
        figAxForm
        legend(L,{'T','D'},'location','northwest')
        title(sprintf('%s, Dist Resp Cells, n=%s,',avName{iav},...
            num2str(sum(ind))))
        vline(visRTwindow(1),'b--')
    else
        legend(L,{'T','D'},'location','northwest')
        title(sprintf('%s, Dist Resp Cells, n=%s',avName{iav},...
            num2str(sum(ind))))
        vline(audRTwindow(1),'b--')
    end
end
print([fnout 'targetTC'],'-dpdf','-fillpage')
%% decoding analysis figures

% setFigParams4Print('landscape')
% figure
% colormap(brewermap([],'*RdBu'));
% 
% targetRespAll = mean(targetAnalysis.tc{allTrialsInd,visualTrials}...
%     (respwin_target,:),1);
% [~,targetSortInd] = sort(targetRespAll);
% subplot 131
% % lateCycRespAll = mean(antiAnalysis.lateCycTC{allTrialsInd}(respwin,:),1);
% % [~,lateCycSortInd] = sort(lateCycRespAll);
% hm = flipud(antiAnalysis.lateCycTC{visualTrials}...
%     (tcStartFrame:cycTCEndFr,targetSortInd)');
% imagesc(hm)
% caxis(hmLim)
% colorbar
% figXAxis([],'Time (ms)',[],ttLabelFr_cyc,ttLabel_cyc)
% figYAxis([],'Cell # (Target Resp Sorted)',[])
% figAxForm
% title('Late Dist. Resp, All Cells')
% 
% subplot 132
% % targetRespAll = mean(targetAnalysis.tc{allTrialsInd,visualTrials}...
% %     (respwin_target,:),1);
% hm = flipud(targetAnalysis.tc{1,visualTrials}...
%     ((tcStartFrame-1):(cycTCEndFr-1),targetSortInd)');
% imagesc(hm)
% caxis(hmLim)
% colorbar
% figXAxis([],'Time (ms)',[],ttLabelFr_cyc,ttLabel_cyc)
% figYAxis([],'Cell # (Target Resp Sorted)',[])
% figAxForm
% title('Hard Target Resp, All Cells')
% 
% subplot 133
% hm = flipud(targetAnalysis.tc{2,visualTrials}...
%     ((tcStartFrame-1):(cycTCEndFr-1),targetSortInd)');
% imagesc(hm)
% caxis(hmLim)
% colorbar
% figXAxis([],'Time (ms)',[],ttLabelFr_cyc-1,ttLabel_cyc)
% figYAxis([],'Cell # (Target Resp Sorted)',[])
% figAxForm
% title('Easy Target Resp, All Cells')
% 
% print([fnout 'heatmapVisTrialsAllCells_dist&target'],'-dpdf','-fillpage')
% 
% figure
% colormap gray
% suptitle('Responsive Cells IDs, sorted by target response')
% 
% subplot 131
% ind = cellInfo.lateCycRespCells;
% indBW = ind; indBW(ind) = 0; indBW(~ind) = 1;
% imagesc(flipud(indBW(targetSortInd)))
% title(sprintf('Late Cyc. Resp. Cells n=%s',num2str(sum(ind))))
% figYAxis([],'Cell # (Target Resp Sorted)',[])
% figAxForm
% 
% subplot 132
% ind = cellInfo.targetRespCells;
% indBW = ind; indBW(ind) = 0; indBW(~ind) = 1;
% imagesc(flipud(indBW(targetSortInd)))
% title(sprintf('Target Resp. Cells n=%s',num2str(sum(ind))))
% figYAxis([],'Cell # (Target Resp Sorted)',[])
% figAxForm
% 
% print([fnout 'cellIDforHeatmap_dist&target'],'-dpdf','-fillpage')
% 
% if strcmp(ds,'FSAV_attentionV1') | strcmp(ds,'FSAV_attentionV1_noAttn')
%    
%     dcModel = struct;
%     dcModel(1).cellTypeName = {'Tar.';'Dist.';'N.P.'};
%     dcModel(1).name = 'Stimulus';
%     dcModel(2).name = 'Choice';
%     cellTypeIDAllCells = sum(cat(2,(cellInfo.lateCycRespCells & ~cellInfo.targetRespCells).*1,...
%         (cellInfo.lateCycRespCells & cellInfo.targetRespCells).*2,...
%         (~cellInfo.lateCycRespCells & cellInfo.targetRespCells).*3),2);
%     dcModel(1).cellInd = [];
% %     dcModel(1).cellTypeID = cellTypeIDAllCells(cellInfo.dcModelCells);
%     for imod = 1:2
%         dcModel(imod).av(1).name = 'Visual';
%         dcModel(imod).av(2).name = 'Auditory';
%         dcModel(imod).comboWeight = [];
%         for iav = 1:2
%             dcModel(2).av(iav).pctYes = nan(3,nexp);
%             dcModel(imod).av(iav).dv = nan(1,nexp);
%             dcModel(imod).av(iav).correlation = [];
%             dcModel(imod).av(iav).weight = [];
%             dcModel(imod).av(iav).pctCorrect = nan(4,nexp);
%             dcModel(imod).av(iav).testPerformance = nan(4,nexp);
%             dcModel(imod).av(iav).pctCorrect_movRespWin = nan(nMovWin,nexp);
%             dcModel(imod).av(iav).pctCorrect_otherAV = nan(1,nexp);
%         end
%     end
%     for iav = 1:2
%         dcModel(1).comboWeight = cat(1,dcModel(1).comboWeight,...
%                 decodeAnalysis(iexp).comboTrainWeightTarget);
%         dcModel(2).comboWeight = cat(1,dcModel(2).comboWeight,...
%                 decodeAnalysis(iexp).comboTrainWeightDetect);
%         for iexp = 1:nexp
%             if iav == 1
%                 dcModel(1).cellInd = cat(1,dcModel(1).cellInd,...
%                     decodeAnalysis(iexp).cellInd);
%             end
%             if decodeAnalysis(iexp).av(visualTrials).pctCorrectAllTarget_holdout < pctCorrThresh || ...
%                     decodeAnalysis(iexp).av(visualTrials).pctCorrectAllDetect_holdout < pctCorrThresh
%                 continue
%             end
%             d = decodeAnalysis(iexp).av(iav);
%             if isempty(d.weightTarget)
%                 d.weightTarget = [];
%                 d.weightDetect = [];
%             end
%             dcModel(1).av(iav).dv(iexp) = d.dvTarget;
%             dcModel(2).av(iav).dv(iexp) = d.dvDetect;
% %             dcModel(2).av(iav).pctYes(:,iexp) = d.pctYes;
%             dcModel(1).av(iav).correlation = cat(2, ...
%                 dcModel(1).av(iav).correlation,d.correlationTarget);
%             dcModel(1).av(iav).weight = cat(1, ...
%                 dcModel(1).av(iav).weight,d.weightTarget);
%             dcModel(1).av(iav).pctCorrect(:,iexp) = cat(2,...
%                 d.pctCorrectXStimTarget_holdout,d.pctCorrectAllTarget_holdout);
%             dcModel(1).av(iav).testPerformance(:,iexp) = cat(2,...
%                 d.pctCorrectXStimTarget_train,d.pctCorrectAllTarget_train) - ...
%                 cat(2,...
%                 d.pctCorrectXStimTarget_holdout,d.pctCorrectAllTarget_holdout);
%             dcModel(1).av(iav).pctCorrect_movRespWin(:,iexp) = d.pctCorrectTargetMovRespWin_holdout;
%             dcModel(1).av(iav).pctCorrect_otherAV(iexp) = d.pctCorrectTarget_otherAV;
%             dcModel(2).av(iav).correlation = cat(2, ...
%                 dcModel(2).av(iav).correlation,d.correlationDetect);
%             dcModel(2).av(iav).weight = cat(1, ...
%                 dcModel(2).av(iav).weight,d.weightDetect);
%             dcModel(2).av(iav).pctCorrect(:,iexp) = cat(2,...
%                 d.pctCorrectXStimDetect_holdout,d.pctCorrectAllDetect_holdout);
%             dcModel(2).av(iav).testPerformance(:,iexp) = cat(2,...
%                 d.pctCorrectXStimDetect_train,d.pctCorrectAllDetect_train) - ...
%                 cat(2,...
%                 d.pctCorrectXStimDetect_holdout,d.pctCorrectAllDetect_holdout);
%             dcModel(2).av(iav).pctCorrect_movRespWin(:,iexp) = d.pctCorrectDetectMovRespWin_holdout;
%             dcModel(2).av(iav).pctCorrect_otherAV(iexp) = d.pctCorrectDetect_otherAV;
%             
%         end
%     end
%     stimLabel = {'0';'HT';'ET';'All'};
%     mdlXStimAlpha = 0.05./2;
% 
% % response matrix (Figure 4 Schematic)
% if strcmp(ds,'FSAV_attentionV1') 
%     iexp=1;
%     for iav = 1:2
%     dc = decodeDataExpt(iexp).av(iav);
%     ind = decodeAnalysis(iexp).cellInd;
% %     trialInd = dcTrialsUsed{iexp,iav};
% %     trOut = dc.outcome(trialInd);
%     trOut = dc.outcome;
% 
% %     r = dc.tc(:,ind,trialInd);
%     r = dc.tc(:,ind,:);
%     [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOut);
%     wDetect = decodeAnalysis(iexp).av(iav).weightDetect;
%     wTarget = decodeAnalysis(iexp).av(iav).weightTarget;
%     [~,detectWeightSort] = sort(wDetect);
%     [~,targetWeightSort] = sort(wTarget);
% 
%     rDetect = {mean(r(:,:,detectTrInd==0),3);mean(r(:,:,detectTrInd==1),3)};
%     rTarget = {mean(r(:,:,targetTrInd==0),3);mean(r(:,:,targetTrInd==1),3)};
%     rDetectErr = {ste(r(:,:,detectTrInd==0),3);ste(r(:,:,detectTrInd==1),3)};
%     rTargetErr = {ste(r(:,:,targetTrInd==0),3);ste(r(:,:,targetTrInd==1),3)};
% 
%     dfig_n = figure;
%     suptitle(sprintf('%s, %s, No Resp, Detect Weight Sorted',...
%         decodeDataExpt(iexp).exptName,avName{iav}))
%     dfig_y = figure;
%     suptitle(sprintf('%s, %s, Yes Resp, Detect Weight Sorted',...
%         decodeDataExpt(iexp).exptName,avName{iav}))
%     tfig_d = figure;
%     suptitle(sprintf('%s, %s, Dist. Resp, Target Weight Sorted',...
%         decodeDataExpt(iexp).exptName,avName{iav}))
%     tfig_t = figure;
%     suptitle(sprintf('%s, %s, Tar. Resp, Target Weight Sorted',...
%         decodeDataExpt(iexp).exptName,avName{iav}))
%     for icell = 1:sum(ind)
%         cInd = detectWeightSort(icell);
%         figure(dfig_n);
%         subplot(5,2,icell)
% %         plot(tt_cycTC,rDetect{1}(tcStartFrame:nFr_cyc,cInd),'k-')
%         shadedErrorBar_chooseColor(tt_cycTC,rDetect{1}(tcStartFrame:nFr_cyc,cInd),...
%             rDetectErr{1}(tcStartFrame:nFr_cyc,cInd),[0 0 0]);
%         figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],[0,200],[0,200])
%         figYAxis([],'dF/F',cellRespTCLim)  
%         hline(0,'r:')
%         vline(respWinTT,'k--')
%         figAxForm
%         title(sprintf('N Resp., Cell #%s, W=%s',num2str(cInd),...
%             num2str(round(wDetect(cInd),2,'significant'))))
%         figure(dfig_y);
%         subplot(5,2,icell)
% %         plot(tt_cycTC,rDetect{2}(tcStartFrame:nFr_cyc,cInd),'k-')
%         shadedErrorBar_chooseColor(tt_cycTC,rDetect{2}(tcStartFrame:nFr_cyc,cInd),...
%             rDetectErr{2}(tcStartFrame:nFr_cyc,cInd),[0 0 0]);
%         figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],[0,200],[0,200])
%         figYAxis([],'dF/F',cellRespTCLim)  
%         hline(0,'r:')
%         vline(respWinTT,'k--')
%         figAxForm
%         title(sprintf('Y Resp., Cell #%s, W=%s',num2str(cInd),...
%             num2str(round(wDetect(cInd),2,'significant'))))
%     end
%     figure(dfig_n);
%     print(...
%         [fnout 'detectModel_' avName{iav}(1:3) '_noResp_' decodeDataExpt(iexp).exptName],...
%         '-dpdf','-fillpage')
%     figure(dfig_y);
%     print(...
%         [fnout 'detectModel_' avName{iav}(1:3) '_yesResp_' decodeDataExpt(iexp).exptName],...
%         '-dpdf','-fillpage')
%     for icell = 1:sum(ind)
%         figure(tfig_d);
%         subplot(5,2,icell)
%         cInd = targetWeightSort(icell);
% %         plot(tt_cycTC,rTarget{1}(tcStartFrame:nFr_cyc,cInd),'k-')
%         shadedErrorBar_chooseColor(tt_cycTC,rTarget{1}(tcStartFrame:nFr_cyc,cInd),...
%             rTargetErr{1}(tcStartFrame:nFr_cyc,cInd),[0 0 0]);
%         figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],[0,200],[0,200])
%         figYAxis([],'dF/F',cellRespTCLim)  
%         hline(0,'r:')
%         vline(respWinTT,'k--')
%         figAxForm
%         title(sprintf('D Resp., Cell #%s, W=%s',num2str(cInd),...
%             num2str(round(wTarget(cInd),2,'significant'))))
%         figure(tfig_t);
%         subplot(5,2,icell)
% %         plot(tt_cycTC,rTarget{2}(tcStartFrame:nFr_cyc,cInd),'k-')
%         shadedErrorBar_chooseColor(tt_cycTC,rTarget{2}(tcStartFrame:nFr_cyc,cInd),...
%             rTargetErr{2}(tcStartFrame:nFr_cyc,cInd),[0 0 0]);
%         figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],[0,200],[0,200])
%         figYAxis([],'dF/F',cellRespTCLim)  
%         hline(0,'r:')
%         vline(respWinTT,'k--')
%         figAxForm
%         title(sprintf('T Resp., Cell #%s, W=%s',num2str(cInd),...
%             num2str(round(wTarget(cInd),2,'significant'))))
%     end
%     figure(tfig_d);
%     print(...
%         [fnout 'targetModel_' avName{iav}(1:3) '_distResp_' decodeDataExpt(iexp).exptName],...
%         '-dpdf','-fillpage')
%     figure(tfig_t);
%     print(...
%         [fnout 'targetModel_' avName{iav}(1:3) '_tarResp_' decodeDataExpt(iexp).exptName],...
%         '-dpdf','-fillpage')
% end
% end
% % within modality decoding in V1  (Figure 4)
%     figure;
%     for imod = 1:2
%         if imod == 1
%             iplot = 0;
%         else
%             iplot = 2;
%         end
%         for iav = 1:2
%             y = dcModel(imod).av(iav).pctCorrect;
%             yerr = ste(dcModel(imod).av(iav).pctCorrect,2);
%             stimInd = sum(~isnan(y),2) > minTrN;
%             y = y(stimInd,:);
%             yerr = yerr(stimInd,:);
%             x = 1:sum(stimInd);
%             [mdlTest,mdlP] = ttest(y,0.5,'dim',2,'alpha',mdlXStimAlpha);
%             subplot(2,2,iplot+iav)
%             for iexp = 1:nexp
%     %             ind = ~isnan(y(:,iexp));
%                 hold on
%                 plot(1:length(x),y(:,iexp),'k.-');
%             end
%             for itest = 1:length(x)
%                 if mdlTest(itest) == 1 
%                     text(x(itest),1,sprintf('*%s',...
%                         num2str(round(mdlP(itest),2,'significant'))))
%                 else
%                     text(x(itest),1,sprintf('%s',...
%                         num2str(round(mdlP(itest),2,'significant'))))
%                 end
%             end
%             errorbar(1:length(x),nanmean(y,2),yerr,'.-','MarkerSize',15)
%             figXAxis([],'Stim',[0 length(x)+1],x,stimLabel(stimInd))
%             hline(0.5,'k:')
%             figYAxis([],'Fraction Correct',[0 1])
%             figAxForm
%             title(sprintf('%s Model, Train %s, Test %s',dcModel(imod).name,...
%                 dcModel(imod).av(iav).name,dcModel(imod).av(iav).name))
%         end
%     end
% 
%     print([fnout 'decodeModelPctCorrectXStim'],'-dpdf','-fillpage')
% 
%     figure;
%     suptitle('Performace of trained model on hold-out trials')
%     for imod = 1:2
%         if imod == 1
%             iplot = 0;
%         else
%             iplot = 2;
%         end
%         for iav = 1:2
%             y = dcModel(imod).av(iav).testPerformance;
%             yerr = ste(dcModel(imod).av(iav).testPerformance,2);
%             stimInd = sum(~isnan(y),2) > minTrN;
%             y = y(stimInd,:);
%             yerr = yerr(stimInd,:);
%             x = 1:sum(stimInd);
%             [mdlTest,mdlP] = ttest(y,0,'dim',2,'alpha',mdlXStimAlpha);
%             subplot(2,2,iplot+iav)
%             for iexp = 1:nexp
%                 ind = ~isnan(y(:,iexp));
%                 hold on
%                 plot(x(ind),y(ind,iexp),'k.-');
%             end
%             for itest = 1:length(x)
%                 if mdlTest(itest) == 1
%                     text(x(itest),0.5,sprintf('*%s',...
%                         num2str(round(mdlP(itest),2,'significant'))))
%                 else
%                     text(x(itest),0.5,sprintf('%s',...
%                         num2str(round(mdlP(itest),2,'significant'))))
%                 end
%             end
%             errorbar(x,nanmean(y,2),yerr,'.-','MarkerSize',15)
%             figXAxis([],'Stim',[0 length(x)+1],1:length(x),stimLabel(stimInd))
%             hline(0,'k:')
%             figYAxis([],'Train - HO Frac. Correct',[-0.5 0.5])
%             figAxForm
%             title(sprintf('%s Model, Train %s, Test %s',dcModel(imod).name,...
%                 dcModel(imod).av(iav).name,dcModel(imod).av(iav).name))
%         end
%     end
% 
%     print([fnout 'decodeModelPerformanceXStim'],'-dpdf','-fillpage')
% 
% % timing of decoding (Suppl. Fig. 8)
% nwins=15;
%     figure
%     suptitle('Percent Correct of responses analyzed at different time windows around stimulus; red is used analysis window; blue is min RT')
%     winInd = [1:3,5:nwins];
%     x = movWinLabelMs(winInd);
%     for imod = 1:2
%         if imod == 1
%             iplot = 0;
%         else
%             iplot = 2;
%         end
%         for iav = 1:2
%             y = dcModel(imod).av(iav).pctCorrect_movRespWin(winInd,:);
%             yerr = ste(dcModel(imod).av(iav).pctCorrect_movRespWin(winInd,:),2);
%             targetTimes = x(x>0);
%             firstTimes = nan(1,nexp);
%             for i = 1:nexp
%                 ind = find(y(x>0,i)>pctCorrThresh,1);
%                 if ~isnan(ind)
%                     firstTimes(i) = targetTimes(ind);
%                 end
%             end
%             subplot(2,2,iplot+iav)
%             for iexp = 1:nexp
%                 ind = ~isnan(y(:,iexp));
%                 hold on
%                 plot(x(ind),y(ind,iexp),'k.-');
%             end
%             errorbar(x,nanmean(y,2),yerr,'.-','MarkerSize',15)
%             errorbar(nanmean(firstTimes),0.9,[],[],ste(firstTimes,2),ste(firstTimes,2),'.')
%             text(0,0.95,[num2str(round(nanmean(firstTimes),2,'significant'))...
%                 '+/-' num2str(round(ste(firstTimes,2),2,'significant'))])
%             figXAxis([],'Resp. Win. Center (Relative to Stim)',[x(1) x(end)],x(1:2:length(x)),round(x(1:2:length(x))))
%             hline(0.5,'k:')
%             figYAxis([],'Frac. Correct',[0 1])
%             figAxForm
%             title(sprintf('%s Model, Train %s, Test %s',dcModel(imod).name,...
%                 dcModel(imod).av(iav).name,dcModel(imod).av(iav).name))
%             vline(respWinTT,'r--')
%             vline(baseWinTT,'k--')
%             vline(0,'k:')
%             if iav == 1
%                 vline(visRTwindow(1),'b--')
%             else
%                 vline(audRTwindow(1),'b--')
%             end
%         end
%     end
% 
%     print([fnout 'decodeModelPerformanceXRespWindow'],'-dpdf','-fillpage')
% 
% % within modality weights    
% for iav = 1:2
%         figure
%         suptitle(sprintf('%s,Task correlation and weights of cells used in model',...
%             dcModel(1).av(iav).name))
% 
%         weightInd = dcModel(imod).av(iav).weight > 0;
%         binnedTargetWeight = nan(1,2);
%         binnedTargetWeightErr = nan(1,2);
%         binnedDetectWeight = nan(1,2);
%         binnedDetectWeightErr = nan(1,2);
%         binnedDetectWeightTest = nan(1,2);
%         binnedTargetWeight(1) = mean(dcModel(1).av(iav).weight(~weightInd));
%         binnedTargetWeight(2) = mean(dcModel(1).av(iav).weight(weightInd));
%         binnedTargetWeightErr(1) = ste(dcModel(1).av(iav).weight(~weightInd),1);
%         binnedTargetWeightErr(2) = ste(dcModel(1).av(iav).weight(weightInd),1);
%         binnedDetectWeight(1) = mean(dcModel(2).av(iav).weight(~weightInd));
%         binnedDetectWeight(2) = mean(dcModel(2).av(iav).weight(weightInd));
%         binnedDetectWeightErr(1) = ste(dcModel(2).av(iav).weight(~weightInd),1);
%         binnedDetectWeightErr(2) = ste(dcModel(2).av(iav).weight(weightInd),1);
%         [~,binnedDetectWeightTest(1)] = ttest(dcModel(2).av(iav).weight(~weightInd));
%         [~,binnedDetectWeightTest(2)] = ttest(dcModel(2).av(iav).weight(weightInd));
% 
% %         for imod = 1:2
% %             x = dcModel(imod).av(iav).correlation;
% %             y = dcModel(imod).av(iav).weight;
% % 
% %             subplot(2,2,imod)
% %             plot(x,y,'.')
% %             figXAxis([],'Correlation',[-1 1])
% %             figYAxis([],'Weight',weightLim)
% %             figAxForm
% %             vline(0,'k:')
% %             hline(0,'k:')
% %             title(sprintf('%s Model',dcModel(imod).name))
% %         end
% 
%         subplot(2,2,3)
%         x = dcModel(1).av(iav).weight;
%         y = dcModel(2).av(iav).weight;
%         [corrmat,pmat] = corrcoef(x,y);
%         plot(x,y,'.')
%         figXAxis([],'Target Weight',[-2 3])
%         figYAxis([],'Detect Weight',[-2 3])
%         figAxForm
%         hline(0,'k:')
%         vline(0,'k:')
%         title(sprintf('Corr=%s, p=%s',num2str(round(corrmat(1,2),2,'significant')),...
%             num2str(round(pmat(1,2),2,'significant'))))
% 
%         subplot(2,2,4)
%         hold on
%         for i = 1:2
%             if i == 1
%                 ind = ~weightInd;
%             elseif i == 2
%                 ind = weightInd;
%             end
%             y = dcModel(2).av(iav).weight(ind);
%             x = binnedTargetWeight(i).*(ones(sum(ind),1));
%             plot(x,y,'.')
%             text(binnedTargetWeight(i),weightLim(end),num2str(round(...
%                 binnedDetectWeightTest(i),2,'significant')))
%         end
%         errorbar(binnedTargetWeight,binnedDetectWeight,binnedDetectWeightErr,'.-')
%         figXAxis([],'Target Weight',binnedWeightLim)
%         figYAxis([],'Detect Weight',weightLim)
%         figAxForm
%         hline(0,'k:')
%         vline(0,'k:')
% 
%         dcModel(1).av(iav).binnedWeight = binnedTargetWeight;
%         dcModel(1).av(iav).binnedWeightErr = binnedTargetWeightErr;
%         dcModel(2).av(iav).binnedWeight = binnedDetectWeight;
%         dcModel(2).av(iav).binnedWeightErr = binnedDetectWeightErr;
%         dcModel(2).av(iav).binnedWeightTest = binnedDetectWeightTest;
% 
%         print([fnout 'decodeModelWeights_' dcModel(1).av(iav).name],'-dpdf','-fillpage')
% end
% 
% % orthogonal coding (Figure 5)
%      
%     WxAttnFig = figure;
%     suptitle('Target Model Attention - Late Resp. Cells')
%     WFig = figure;
%     suptitle('Target Model Weight x vis AuROC Group - Late Resp. Cells')
%     weightDataFor2way = cell(1,2);
%     weightIndFor2way = cell(1,2);
%     for iav = 1:2
%         figure(WxAttnFig)
%         subplot(2,2,iav)
%     %     ind = cellInfo.lateCycRespCells;
%         ind1 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & cellInfo.lateStimAuROCTest...
%             & cellInfo.lateStimAuROC > 0.5;
%         ind2 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & cellInfo.lateStimAuROCTest...
%             & cellInfo.lateStimAuROC < 0.5;
%         ind3 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & ~cellInfo.lateStimAuROCTest;
%         n = {[dcModel(1).cellTypeName{1} ';n=' num2str(sum(ind1))];...
%             [dcModel(1).cellTypeName{2} ';n=' num2str(sum(ind2))];...
%             [dcModel(1).cellTypeName{3} ';n=' num2str(sum(ind3))]};
% 
%         x = antiAnalysis.lateCycSI(ind1);
%         y = cellInfo.targetWeight{iav}(ind1);    
%         plot(x,y,'.')
%         hold on
%         x = antiAnalysis.lateCycSI(ind2);
%         y = cellInfo.targetWeight{iav}(ind2);    
%         plot(x,y,'.')
%         x = antiAnalysis.lateCycSI(ind3);
%         y = cellInfo.targetWeight{iav}(ind3);    
%         plot(x,y,'.')
%         [corrmat,pmat] = corrcoef(antiAnalysis.lateCycSI(cellInfo.lateCycRespCells),...
%             cellInfo.targetWeight{iav}(cellInfo.lateCycRespCells));
%         r = corrmat(1,2);
%         p = pmat(1,2);
% 
%         figXAxis([],'Selectivity',siLim)
%         figYAxis([],'Weight',weightLim)
%         figAxForm
%         vline(0,'k:')
%         hline(0,'k:')
%         title(sprintf('Train %s, R=%s, p=%s',dcModel(1).av(iav).name,...
%             num2str(round(r,2,'significant')),num2str(round(p,2,'significant'))))
% 
%         subplot(2,2,iav+2)
% 
%         x = mean(antiAnalysis.lateCycSI(ind1));
%         xerr = ste(antiAnalysis.lateCycSI(ind1),2);
%         y = nanmean(cellInfo.targetWeight{iav}(ind1));
%         yerr = ste(cellInfo.targetWeight{iav}(ind1),2);    
%         errorbar(x,y,yerr,yerr,xerr,xerr,'.')
%         hold on
%         x = mean(antiAnalysis.lateCycSI(ind2));
%         xerr = ste(antiAnalysis.lateCycSI(ind2),2);
%         y = nanmean(cellInfo.targetWeight{iav}(ind2));
%         yerr = ste(cellInfo.targetWeight{iav}(ind2),2);    
%         errorbar(x,y,yerr,yerr,xerr,xerr,'.')
%         x = mean(antiAnalysis.lateCycSI(ind3));
%         xerr = ste(antiAnalysis.lateCycSI(ind3),2);
%         y = nanmean(cellInfo.targetWeight{iav}(ind3));
%         yerr = ste(cellInfo.targetWeight{iav}(ind3),2);    
%         errorbar(x,y,yerr,yerr,xerr,xerr,'.')
% 
%         figXAxis([],'Selectivity',siLimSum)
%         figYAxis([],'Weight',weightLimSum)
%         figAxForm
%         vline(0,'k:')
%         hline(0,'k:')
%         title(sprintf('Train %s',dcModel(1).av(iav).name))
%         L = legend(n,'location','southwest');
%         title(L,'Neuron Pref. in Vis Trials')
%         
%         figure(WFig)
%         pt = [];
%         subplot(1,2,iav)
%         y = nanmean(cellInfo.targetWeight{iav}(ind1));
%         yerr = ste(cellInfo.targetWeight{iav}(ind1),2);    
%         errorbar(1,y,yerr,yerr,'.')
%         [~,pt(1)] = ttest(cellInfo.targetWeight{iav}(ind1));
%         y = nanmean(cellInfo.targetWeight{iav}(ind2));
%         yerr = ste(cellInfo.targetWeight{iav}(ind2),2);   
%         hold on
%         errorbar(2,y,yerr,yerr,'.')
%         [~,pt(2)] = ttest(cellInfo.targetWeight{iav}(ind2));
%         y = nanmean(cellInfo.targetWeight{iav}(ind3));
%         yerr = ste(cellInfo.targetWeight{iav}(ind3),2);    
%         errorbar(3,y,yerr,yerr,'.')
%         [~,pt(3)] = ttest(cellInfo.targetWeight{iav}(ind3));
%         x = cat(2,ones(1,sum(ind1)),ones(1,sum(ind2)).*2,ones(1,sum(ind3)).*3);
%         y = cat(2,cellInfo.targetWeight{iav}(ind1),...
%             cellInfo.targetWeight{iav}(ind2),cellInfo.targetWeight{iav}(ind3));
%         [p,~,stats] = anova1(y,x,'off');
%         posthoc = multcompare(stats,[],'off');
%         fprintf('%s Target Weight x Vis auROC Model\n',avName{iav})
%         disp(posthoc(:,[1,2,end]))
%         fprintf('%s Target Weight Vs. Zero x Vis auROC\n',avName{iav})
%         disp(round(pt,2,'significant'))
%         figXAxis([],'Vis. auROC Group',[0 4],1:3,{'T','D','NP'})
%         figYAxis([],'Weight',weightLimSum)
%         figAxForm
%         hline(0,'k:')
%         title(sprintf('Train %s, p=%s',dcModel(1).av(iav).name,...
%             num2str(round(p,2,'significant'))))
%         L = legend(n,'location','southwest');
%         title(L,'Neuron Pref. in Vis Trials, p=%s')
%         
%         weightDataFor2way{iav} = y;
%         weightIndFor2way{iav} = x;
%     end
% %     [p,~,stats] = anovan(cell2mat(weightDataFor2way),...
% %         {cell2mat(weightIndFor2way) cat(2,ones(1,length(weightDataFor2way{1})),...
% %         ones(1,length(weightDataFor2way{1}))*2)},'model','interaction');
% %     multcompare(stats)
%     
%     figure(WFig)
%     print([fnout 'targetModelWeightxGroup'],'-dpdf')
%     figure(WxAttnFig)
%     print([fnout 'targetModelWeightXAttn'],'-dpdf','-fillpage')
% 
%     WxAttnFig = figure;
%     suptitle('Detect Model Attention - Late Resp. Cells')
%     WFig = figure;
%     suptitle('Detect Model Weight x vis AuROC Group - Late Resp. Cells')
%     for iav = 1:2
%         figure(WxAttnFig)
%         subplot(2,2,iav)
%     %     ind = cellInfo.lateCycRespCells;
%         ind1 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & cellInfo.lateStimAuROCTest...
%             & cellInfo.lateStimAuROC > 0.5;
%         ind2 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & cellInfo.lateStimAuROCTest...
%             & cellInfo.lateStimAuROC < 0.5;
%         ind3 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & ~cellInfo.lateStimAuROCTest;
%         n = {[dcModel(1).cellTypeName{1} ';n=' num2str(sum(ind1))];...
%             [dcModel(1).cellTypeName{2} ';n=' num2str(sum(ind2))];...
%             [dcModel(1).cellTypeName{3} ';n=' num2str(sum(ind3))]};
% 
%         x = antiAnalysis.lateCycSI(ind1);
%         y = cellInfo.detectWeight{iav}(ind1);    
%         plot(x,y,'.')
%         hold on
%         x = antiAnalysis.lateCycSI(ind2);
%         y = cellInfo.detectWeight{iav}(ind2);    
%         plot(x,y,'.')
%         x = antiAnalysis.lateCycSI(ind3);
%         y = cellInfo.detectWeight{iav}(ind3);    
%         plot(x,y,'.')
%         [corrmat,pmat] = corrcoef(antiAnalysis.lateCycSI(cellInfo.lateCycRespCells),...
%             cellInfo.detectWeight{iav}(cellInfo.lateCycRespCells));
%         r = corrmat(1,2);
%         p = pmat(1,2);
% 
%         figXAxis([],'Selectivity',siLim)
%         figYAxis([],'Weight',weightLim)
%         figAxForm
%         vline(0,'k:')
%         hline(0,'k:')
%         title(sprintf('Train %s, R=%s, p=%s',dcModel(2).av(iav).name,...
%             num2str(round(r,2,'significant')),num2str(round(p,2,'significant'))))
% 
%         subplot(2,2,iav+2)
% 
%         x = mean(antiAnalysis.lateCycSI(ind1));
%         xerr = ste(antiAnalysis.lateCycSI(ind1),2);
%         y = nanmean(cellInfo.detectWeight{iav}(ind1));
%         yerr = ste(cellInfo.detectWeight{iav}(ind1),2);    
%         errorbar(x,y,yerr,yerr,xerr,xerr,'.')
%         hold on
%         x = mean(antiAnalysis.lateCycSI(ind2));
%         xerr = ste(antiAnalysis.lateCycSI(ind2),2);
%         y = nanmean(cellInfo.detectWeight{iav}(ind2));
%         yerr = ste(cellInfo.detectWeight{iav}(ind2),2);    
%         errorbar(x,y,yerr,yerr,xerr,xerr,'.')
%         x = mean(antiAnalysis.lateCycSI(ind3));
%         xerr = ste(antiAnalysis.lateCycSI(ind3),2);
%         y = nanmean(cellInfo.detectWeight{iav}(ind3));
%         yerr = ste(cellInfo.detectWeight{iav}(ind3),2);    
%         errorbar(x,y,yerr,yerr,xerr,xerr,'.')
% 
%         figXAxis([],'Selectivity',siLimSum)
%         figYAxis([],'Weight',weightLimSum)
%         figAxForm
%         vline(0,'k:')
%         hline(0,'k:')
%         title(sprintf('Train %s',dcModel(1).av(iav).name))
%         L = legend(n,'location','southwest');
%         title(L,'Neuron Pref. in Vis Trials')
%         
%         figure(WFig)
%         subplot(1,2,iav)
%         y = nanmean(cellInfo.detectWeight{iav}(ind1));
%         yerr = ste(cellInfo.detectWeight{iav}(ind1),2);    
%         errorbar(1,y,yerr,yerr,'.')
%         [~,pt(1)] = ttest(cellInfo.targetWeight{iav}(ind1));
%         y = nanmean(cellInfo.detectWeight{iav}(ind2));
%         yerr = ste(cellInfo.detectWeight{iav}(ind2),2);   
%         hold on
%         errorbar(2,y,yerr,yerr,'.')
%         [~,pt(2)] = ttest(cellInfo.targetWeight{iav}(ind2));
%         y = nanmean(cellInfo.detectWeight{iav}(ind3));
%         yerr = ste(cellInfo.detectWeight{iav}(ind3),2);    
%         errorbar(3,y,yerr,yerr,'.')
%         [~,pt(3)] = ttest(cellInfo.targetWeight{iav}(ind3));
%         x = cat(2,ones(1,sum(ind1)),ones(1,sum(ind2)).*2,ones(1,sum(ind3)).*3);
%         y = cat(2,cellInfo.detectWeight{iav}(ind1),...
%             cellInfo.detectWeight{iav}(ind2),cellInfo.detectWeight{iav}(ind3));
%         [p,~,stats] = anova1(y,x,'off');
%         posthoc = multcompare(stats,[],'off');
%         fprintf('%s Detect Weight x Vis auROC Model\n',avName{iav})
%         disp(posthoc(:,[1,2,end]))
%         fprintf('%s Target Weight Vs. Zero x Vis auROC\n',avName{iav})
%         disp(round(pt,2,'significant'))
%         figXAxis([],'Vis. auROC Group',[0 4],1:3,{'T','D','NP'})
%         figYAxis([],'Weight',weightLimSum)
%         figAxForm
%         hline(0,'k:')
%         title(sprintf('Train %s,p=%s',dcModel(1).av(iav).name,...
%             num2str(round(p,2,'significant'))))
%         L = legend(n,'location','southwest');
%         title(L,'Neuron Pref. in Vis Trials')
%     end
%     
%     figure(WFig)
%     print([fnout 'detectModelWeightxGroup'],'-dpdf')
%     figure(WxAttnFig)
%     print([fnout 'detectModelWeightXAttn'],'-dpdf','-fillpage')
% 
%     WxAttnFig = figure;
%     suptitle('Target Model Attention - Late Resp. Cells')
%     WFig = figure;
%     suptitle('Target Model Weight x aud AuROC Group - Late Resp. Cells')
%     for iav = 1:2
%         figure(WxAttnFig)
%         subplot(2,2,iav)
%     %     ind = (cellInfo.lateCycRespCells|cellInfo.targetRespCells);
%         ind1 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & cellInfo.audLateStimAuROCTest...
%             & cellInfo.audLateStimAuROC > 0.5;
%         ind2 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & cellInfo.audLateStimAuROCTest...
%             & cellInfo.audLateStimAuROC < 0.5;
%         ind3 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & ~cellInfo.audLateStimAuROCTest;
%         n = {[dcModel(1).cellTypeName{1} ';n=' num2str(sum(ind1))];...
%             [dcModel(1).cellTypeName{2} ';n=' num2str(sum(ind2))];...
%             [dcModel(1).cellTypeName{3} ';n=' num2str(sum(ind3))]};
% 
%         x = antiAnalysis.lateCycSI(ind1);
%         y = cellInfo.targetWeight{iav}(ind1);    
%         plot(x,y,'.')
%         hold on
%         x = antiAnalysis.lateCycSI(ind2);
%         y = cellInfo.targetWeight{iav}(ind2);    
%         plot(x,y,'.')
%         x = antiAnalysis.lateCycSI(ind3);
%         y = cellInfo.targetWeight{iav}(ind3);    
%         plot(x,y,'.')
% 
%         figXAxis([],'Selectivity',siLim)
%         figYAxis([],'Weight',weightLim)
%         figAxForm
%         vline(0,'k:')
%         hline(0,'k:')
%         title(sprintf('Train %s',dcModel(1).av(iav).name))
% 
%         subplot(2,2,iav+2)
% 
%         x = mean(antiAnalysis.lateCycSI(ind1));
%         xerr = ste(antiAnalysis.lateCycSI(ind1),2);
%         y = nanmean(cellInfo.targetWeight{iav}(ind1));
%         yerr = ste(cellInfo.targetWeight{iav}(ind1),2);    
%         errorbar(x,y,yerr,yerr,xerr,xerr,'.')
%         hold on
%         x = mean(antiAnalysis.lateCycSI(ind2));
%         xerr = ste(antiAnalysis.lateCycSI(ind2),2);
%         y = nanmean(cellInfo.targetWeight{iav}(ind2));
%         yerr = ste(cellInfo.targetWeight{iav}(ind2),2);    
%         errorbar(x,y,yerr,yerr,xerr,xerr,'.')
%         x = mean(antiAnalysis.lateCycSI(ind3));
%         xerr = ste(antiAnalysis.lateCycSI(ind3),2);
%         y = nanmean(cellInfo.targetWeight{iav}(ind3));
%         yerr = ste(cellInfo.targetWeight{iav}(ind3),2);    
%         errorbar(x,y,yerr,yerr,xerr,xerr,'.')
% 
%         figure(WFig)
%         subplot(1,2,iav)
%         y = nanmean(cellInfo.targetWeight{iav}(ind1));
%         yerr = ste(cellInfo.targetWeight{iav}(ind1),2);    
%         errorbar(1,y,yerr,yerr,'.')
%         [~,pt(1)] = ttest(cellInfo.targetWeight{iav}(ind1));
%         y = nanmean(cellInfo.targetWeight{iav}(ind2));
%         yerr = ste(cellInfo.targetWeight{iav}(ind2),2);   
%         hold on
%         errorbar(2,y,yerr,yerr,'.')
%         [~,pt(2)] = ttest(cellInfo.targetWeight{iav}(ind2));
%         y = nanmean(cellInfo.targetWeight{iav}(ind3));
%         yerr = ste(cellInfo.targetWeight{iav}(ind3),2);    
%         errorbar(3,y,yerr,yerr,'.')
%         [~,pt(3)] = ttest(cellInfo.targetWeight{iav}(ind3));
%         x = cat(2,ones(1,sum(ind1)),ones(1,sum(ind2)).*2,ones(1,sum(ind3)).*3);
%         y = cat(2,cellInfo.targetWeight{iav}(ind1),...
%             cellInfo.targetWeight{iav}(ind2),cellInfo.targetWeight{iav}(ind3));
%         [p,~,stats] = anova1(y,x,'off');
%         posthoc = multcompare(stats,[],'off');
%         fprintf('%s Target Weight x Aud auROC Model',avName{iav})
%         disp(posthoc(:,[1,2,end]))
%         fprintf('%s Target Weight Vs. Zero x Aud auROC\n',avName{iav})
%         disp(round(pt,2,'significant'))
%         figXAxis([],'auROC Group',[0 4],1:3,{'T','D','NP'})
%         figYAxis([],'Weight',weightLimSum)
%         figAxForm
%         hline(0,'k:')
%         title(sprintf('Train %s,p=%s',dcModel(1).av(iav).name,...
%             num2str(round(p,2,'significant'))))
%         L = legend(n,'location','southwest');
%         title(L,'Neuron Pref. in Aud Trials')
%     end
% 
%     figure(WFig)
%     print([fnout 'targetModelAudSortWeightXGroup'],'-dpdf')
%     figure(WxAttnFig)
%     print([fnout 'targetModelAudSortWeightXAttn'],'-dpdf','-fillpage')
% 
%     WxAttnFig = figure;
%     suptitle('Detect Model Attention - Late Resp. Cells')
%     WFig = figure;
%     suptitle('Detect Model Weight x aud AuROC Group - Late Resp. Cells')
%     for iav = 1:2
%         figure(WxAttnFig)
%         subplot(2,2,iav)
%     %     ind = (cellInfo.lateCycRespCells|cellInfo.targetRespCells);
%         ind1 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & cellInfo.audLateStimAuROCTest...
%             & cellInfo.audLateStimAuROC > 0.5;
%         ind2 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & cellInfo.audLateStimAuROCTest...
%             & cellInfo.audLateStimAuROC < 0.5;
%         ind3 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & ~cellInfo.audLateStimAuROCTest;
%         n = {[dcModel(1).cellTypeName{1} ';n=' num2str(sum(ind1))];...
%             [dcModel(1).cellTypeName{2} ';n=' num2str(sum(ind2))];...
%             [dcModel(1).cellTypeName{3} ';n=' num2str(sum(ind3))]};
% 
%         x = antiAnalysis.lateCycSI(ind1);
%         y = cellInfo.detectWeight{iav}(ind1);    
%         plot(x,y,'.')
%         hold on
%         x = antiAnalysis.lateCycSI(ind2);
%         y = cellInfo.detectWeight{iav}(ind2);    
%         plot(x,y,'.')
%         x = antiAnalysis.lateCycSI(ind3);
%         y = cellInfo.detectWeight{iav}(ind3);    
%         plot(x,y,'.')
% 
%         figXAxis([],'Selectivity',siLim)
%         figYAxis([],'Weight',weightLim)
%         figAxForm
%         vline(0,'k:')
%         hline(0,'k:')
%         title(sprintf('Train %s',dcModel(1).av(iav).name))
% 
%         subplot(2,2,iav+2)
% 
%         x = mean(antiAnalysis.lateCycSI(ind1));
%         xerr = ste(antiAnalysis.lateCycSI(ind1),2);
%         y = nanmean(cellInfo.detectWeight{iav}(ind1));
%         yerr = ste(cellInfo.detectWeight{iav}(ind1),2);    
%         errorbar(x,y,yerr,yerr,xerr,xerr,'.')
%         hold on
%         x = mean(antiAnalysis.lateCycSI(ind2));
%         xerr = ste(antiAnalysis.lateCycSI(ind2),2);
%         y = nanmean(cellInfo.detectWeight{iav}(ind2));
%         yerr = ste(cellInfo.detectWeight{iav}(ind2),2);    
%         errorbar(x,y,yerr,yerr,xerr,xerr,'.')
%         x = mean(antiAnalysis.lateCycSI(ind3));
%         xerr = ste(antiAnalysis.lateCycSI(ind3),2);
%         y = nanmean(cellInfo.detectWeight{iav}(ind3));
%         yerr = ste(cellInfo.detectWeight{iav}(ind3),2);    
%         errorbar(x,y,yerr,yerr,xerr,xerr,'.')
% 
%         figXAxis([],'Selectivity',siLimSum)
%         figYAxis([],'Weight',weightLimSum)
%         figAxForm
%         vline(0,'k:')
%         hline(0,'k:')
%         title(sprintf('Train %s',dcModel(1).av(iav).name))
%         L = legend(n,'location','southwest');
%         title(L,'Neuron Pref. in Vis Trials')
%         
%         figure(WFig)
%         subplot(1,2,iav)
%         y = nanmean(cellInfo.detectWeight{iav}(ind1));
%         yerr = ste(cellInfo.detectWeight{iav}(ind1),2);    
%         errorbar(1,y,yerr,yerr,'.')
%         [~,pt(1)] = ttest(cellInfo.targetWeight{iav}(ind1));
%         y = nanmean(cellInfo.detectWeight{iav}(ind2));
%         yerr = ste(cellInfo.detectWeight{iav}(ind2),2);   
%         hold on
%         errorbar(2,y,yerr,yerr,'.')
%         [~,pt(2)] = ttest(cellInfo.targetWeight{iav}(ind2));
%         y = nanmean(cellInfo.detectWeight{iav}(ind3));
%         yerr = ste(cellInfo.detectWeight{iav}(ind3),2);    
%         errorbar(3,y,yerr,yerr,'.')
%         [~,pt(3)] = ttest(cellInfo.targetWeight{iav}(ind3));
%         x = cat(2,ones(1,sum(ind1)),ones(1,sum(ind2)).*2,ones(1,sum(ind3)).*3);
%         y = cat(2,cellInfo.detectWeight{iav}(ind1),...
%             cellInfo.detectWeight{iav}(ind2),cellInfo.detectWeight{iav}(ind3));
%         [p,~,stats] = anova1(y,x,'off');
%         posthoc = multcompare(stats,[],'off');
%         fprintf('%s Target Weight x Vis auROC Model',avName{iav})
%         disp(posthoc(:,[1,2,end]))
%         fprintf('%s Target Weight Vs. Zero x Aud auROC\n',avName{iav})
%         disp(round(pt,2,'significant'))
%         figXAxis([],'auROC Group',[0 4],1:3,{'T','D','NP'})
%         figYAxis([],'Weight',weightLimSum)
%         figAxForm
%         hline(0,'k:')
%         title(sprintf('Train %s,p=%s',dcModel(1).av(iav).name,...
%             num2str(round(p,2,'significant'))))
%         L = legend(n,'location','southwest');
%         title(L,'Neuron Pref. in Aud Trials')
%     end
% 
%     figure(WFig)
%     print([fnout 'detectModelAudSortWeightXGroup'],'-dpdf')
%     figure(WxAttnFig)
%     print([fnout 'detectModelAudSortWeightXAttn'],'-dpdf','-fillpage')
% 
%     figure
%     for imod = 1:2
%         for iav = 1:2
%             if imod == 1
%                 subplot(2,2,iav)
%             elseif imod == 2
%                 subplot(2,2,iav+2)
%             end
%             if iav == visualTrials
%                 otherAV = 'Auditory';
%             else
%                 otherAV = 'Visual';
%             end
%             x = 1:2;
%             y = [dcModel(imod).av(iav).pctCorrect(end,:)', dcModel(imod).av(iav).pctCorrect_otherAV'];
%             hold on
%             plot(x,y,'k')
%             errorbar(x,mean(y,1),ste(y,1),'.')
%             [~,p] = ttest(y,0.5);
%             [~,p2] = ttest(y(:,1),y(:,2));
% 
%             figXAxis([],'Test Trials',[0 3],x,{dcModel(imod).av(iav).name,otherAV})
%             figYAxis([],'Fraction Correct',[0 1])
%             figAxForm
%             hline(0.5,'k:')
%             title(sprintf('%s Model, Train %s',dcModel(imod).name,dcModel(imod).av(iav).name))
%             for i = 1:2
%                 if p(i) < 0.05/7
%                     text(i,1,sprintf('*%s',num2str(round(p(i),2,'significant'))))
%                 else
%                     text(i,1,num2str(round(p(i),2,'significant')))
%                 end
%             end
%             if p2 < 0.05/4
%                 text(1.5,0,sprintf('*%s',num2str(round(p2,2,'significant'))))
%             else
%                 text(1.5,0,sprintf('%s',num2str(round(p2,2,'significant'))))
%             end
%         end
%     end
% 
%     print([fnout 'decodeModelPctCorrectOtherAVModel'],'-dpdf','-fillpage')
% 
%     % compare visual and auditory weights
%     figure
%     subplot 221
% %     ind = cellInfo.lateCycRespCells | cellInfo.targetRespCells;
%     x = dcModel(1).av(visualTrials).weight;
%     y = dcModel(1).av(auditoryTrials).weight;
%     hold on
%     plot(x,y,'.')
%     [corrmat,pmat] = corrcoef(x,y);
%     c = corrmat(1,2);
%     p = pmat(1,2);
%     mdl = fitlm(x,y);
%     rsq = mdl.Rsquared.Ordinary;
%     yfit = predict(mdl,x);
%     plot(x,yfit,'-')
%     [xMax,yInd] = max(x);
%     plot([0 xMax],[0 y(yInd)],'-')
%     figXAxis([],'Visual Weight',weightLim)
%     figYAxis([],'Auditory Weight',weightLim)
%     figAxForm
%     title(sprintf('Target Model, Rsq = %s, corr = %s, p = %s',num2str(round(rsq,2,'significant')),...
%         num2str(round(c,2,'significant')),num2str(round(p,2,'significant')))) 
%     
%     oriEdges = [0:15:90];
%     subplot 222
%     targetWeightAVangle = rad2deg(abs(atan(y./x)));
%     histogram(targetWeightAVangle,oriEdges,'Normalization','probability')
%     figXAxis([],'Angle of A/V Weight',[-10 100],0:15:90)
%     figYAxis([],'Fraction of Cells',[0 .3])
%     figAxForm
%     title(sprintf('Target Model,exCell Angle=%s',num2str(round(targetWeightAVangle(yInd),2,'significant'))))
%     
%     
%     subplot 223
%     x = dcModel(2).av(visualTrials).weight;
%     y = dcModel(2).av(auditoryTrials).weight;
%     hold on
%     plot(x,y,'.')
%     [corrmat,pmat] = corrcoef(x,y);
%     p = pmat(1,2);
%     c = corrmat(1,2);
%     mdl = fitlm(x,y);
%     rsq = mdl.Rsquared.Ordinary;
%     yfit = predict(mdl,x);
%     plot(x,yfit,'-')
%     [yMax,xInd] = max(y);
%     plot([0 x(xInd)],[0 yMax],'-')
%     figXAxis([],'Visual Weight',weightLim)
%     figYAxis([],'Auditory Weight',weightLim)
%     figAxForm
%     title(sprintf('Detect Model, Rsq = %s, corr = %s, p = %s',num2str(round(rsq,2,'significant')),...
%         num2str(round(c,2,'significant')),num2str(round(p,2,'significant')))) 
% 
%     subplot 224
%     detectWeightAVangle = rad2deg(abs(atan(y./x)));
%     histogram(detectWeightAVangle,oriEdges,'Normalization','probability')
%     figXAxis([],'Angle of A/V Weight',[-10 100],0:15:90)
%     figYAxis([],'Fraction of Cells',[0 .3])
%     figAxForm
%     title(sprintf('Detect Model,exCell Angle=%s',num2str(round(detectWeightAVangle(yInd),2,'significant'))))
%     
%     print([fnout 'decodeModelCompareAVWeights'],'-dpdf','-fillpage')
%     
%     figure
%     suptitle('Visual Valid and Invalid Trials Matched for Orientation')
%     subplot 121
%     y = [];
%     for iexp = 1:nexp
%         if ~isempty(decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectTarget)
%             y = cat(1,y,[decodeAnalysis(iexp).av(visualTrials).validMatchedPctCorrectTarget_holdout,...
%                 decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectTarget, ...
%                 decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectTarget_testAudModel]);
%         end
%     end
%     plot(1:3,y,'k-')
%     hold on
%     errorbar(1:3,nanmean(y,1),ste(y,1),'.')
%     [~,p]= ttest(y,0.5);
%     for i = 1:3
%         text(i,1,num2str(round(p(i),2,'significant')))
%     end
% %     [~,p] = ttest(y(:,1),y(:,2));
%     [p,~,stats] = anova1(y,[],'off');
%     anovaCmpTarget = multcompare(stats,'display','off');
% %     text(1.5,.1,num2str(round(p,2,'significant')))
%     figXAxis([],'Train-Test-Cue',[0 4],1:3,{'Vis-Vis-Val','Vis-Vis-Inv','Aud-Vis-Inv'})
%     figYAxis([],'Fraction Correct',[0 1])
%     figAxForm
%     hline(0.5,'k:')
%     title(sprintf('Target Model, p=%s',num2str(round(p,2,'significant'))))
%     subplot 122
%     y = [];
%     for iexp = 1:nexp
%         if ~isempty(decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectDetect)
%             y = cat(1,y,[decodeAnalysis(iexp).av(visualTrials).validMatchedPctCorrectDetect_holdout,...
%                 decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectDetect, ...
%                 decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectDetect_testAudModel]);
%         end
%     end
%     plot(1:3,y,'k-')
%     hold on
%     errorbar(1:3,nanmean(y,1),ste(y,1),'.')
%     [~,p]= ttest(y,0.5);
%     for i = 1:3
%         text(i,1,num2str(round(p(i),2,'significant')))
%     end
% %     [~,p] = ttest(y(:,1),y(:,2));
%     [p,~,stats] = anova1(y,[],'off');
%     anovaCmpDetect = multcompare(stats,'display','off');
% %     text(1.5,.1,num2str(round(p,2,'significant')))
%     figXAxis([],'Train-Test-Cue',[0 4],1:3,{'Vis-Vis-Val','Vis-Vis-Inv','Aud-Vis-Inv'})
%     figYAxis([],'Fraction Correct',[0 1])
%     figAxForm
%     hline(0.5,'k:')
%     title(sprintf('Detect Model, p=%s',num2str(round(p,2,'significant'))))  
%     print([fnout 'decodeModelCompareInvPctCorrect'],'-dpdf')
%     
% % test distractor only model between modalities
%    figure;
%     suptitle('Detect Model')
%     distOnlyData = cell(1,2);
%     for iexp = 1:nexp
%         for iav = 1:2
%             subplot(1,2,iav)
%             hold on
%             d = decodeAnalysis(iexp).av(iav);
%             y = [d.pctCorrectAllDetect_holdout,...
%                 d.pctCorrectDetect_otherAV,...
%                 d.pctCorrectDetect_distOnly_holdout,...
%                 d.pctCorrectDetect_distOnly_otherAV];
%             h = plot(1:4,y,'-');
%             h.Color = [0.5 0.5 0.5];
%             distOnlyData{iav}(iexp,:) = y(3:4);
%         end
%     end
%     for iav = 1:2
%         [~,p] = ttest(distOnlyData{iav}(:,1),distOnlyData{iav}(:,2))
%         subplot(1,2,iav)
%         if iav == 1
%             figXAxis([],'Model-Test',[0 5],1:4,...
%                 {'VisTD-VisTD','VisTD-AudTD','VisD-VisD','VisD-AudD'})
%         else
%             figXAxis([],'Model-Test',[0 5],1:4,...
%                 {'AudTD-AudTD','AudTD-VisTD','AudD-AudD','AudD-VisD'})
%         end
%         title(sprintf('compare D-withinAV:D-otherAV, p=%s',num2str(round(p,2,'significant'))))
%         xtickangle(30)
%         figYAxis([],'Fraction Correct',[0 1])
%         figAxForm
%         hline(0.5,'k:')
%     end
%     print([fnout 'decodeModelCompareDistOnlyModelBwModalities'],'-dpdf') 
% % test combo model (Suppl. Fig 9)
%     avName = {'Visual','Auditory'};
%     setFigParams4Print('portrait')
%     figure
%     suptitle('p val in title is within mod vs. combo')
% for iav = 1:2
%     if iav == 1
%         iplot = 1;
%         otherAV = 2;
%     else
%         iplot = 3;
%         otherAV = 1;
%     end
%     subplot(2,2,iplot)
%     y = [];
%     for iexp = 1:nexp
% %             if iav == 1
%             y = cat(1,y,[decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_holdout,...
%                 decodeAnalysis(iexp).av(iav).pctCorrectTarget_comboTrain]);
% %             else
% %                 y = cat(1,y,[decodeAnalysis(iexp).av(iav).pctCorrectTarget_audTrainComboMatch,...
% %                     decodeAnalysis(iexp).av(iav).pctCorrectTarget_comboTrain]);
% %             end
%     end
%     plot(1:2,y,'k-')
%     hold on
%     errorbar(1:2,mean(y,1),ste(y,1),'.')
%     [~,p]= ttest(y,0.5);
%     for i = 1:2
%         text(i,1,num2str(round(p(i),2,'significant')))
%     end
%     figXAxis([],'Train Target',[0 3],1:2,{avName{iav};'Vis+Aud'})
%     figYAxis([],'Fraction Correct',[0 1])
%     figAxForm
%     hline(0.5,'k:')
%     [~,p]= ttest(y(:,1),y(:,2));
%     title(sprintf('Test %s, Target Model, p=%s',avName{iav},num2str(round(p,2,'significant'))))
%     subplot(2,2,iplot+1)
%     y = [];
%     for iexp = 1:nexp
% %             if iav == 1
%             y = cat(1,y,[decodeAnalysis(iexp).av(iav).pctCorrectAllDetect_holdout,...
%                 decodeAnalysis(iexp).av(iav).pctCorrectTarget_comboTrain]);
% %             else
% %                 y = cat(1,y,[decodeAnalysis(iexp).av(iav).pctCorrectDetect_audTrainComboMatch,...
% %                     decodeAnalysis(iexp).av(iav).pctCorrectTarget_comboTrain]);
% %             end
%     end
%     plot(1:2,y,'k-')
%     hold on
%     errorbar(1:2,mean(y,1),ste(y,1),'.')
%     [~,p]= ttest(y,0.5);
%     for i = 1:2
%         text(i,1,num2str(round(p(i),2,'significant')))
%     end
%     figXAxis([],'Train Detect',[0 3],1:2,{avName{iav};'Vis+Aud'})
%     figYAxis([],'Fraction Correct',[0 1])
%     figAxForm
%     hline(0.5,'k:')
%     [~,p]= ttest(y(:,1),y(:,2));
%     title(sprintf('Test %s, Detect Model, p=%s',avName{iav},num2str(round(p,2,'significant'))))
% end
%     print([fnout 'decodeModelComboTrainPctCorrect'],'-dpdf','-fillpage')
%     
%     setFigParams4Print('landscape')
%     figure;
%     suptitle('Task-responsive neurons, bootstrapped weights')
%     ind = cellInfo.targetRespCells | cellInfo.lateCycRespCells;
%     subplot 241
%     x = cellInfo.detectComboWeight(ind);
%     y = cellInfo.detectWeight{visualTrials}(ind);
%     plot(x,y,'.','MarkerSize',10)
%     [c,p] = corrcoef(x,y);
%     hold on
%     figXAxis([],'Combo Weight',weightLim)
%     figYAxis([],'Vis. Only Weight',weightLim)
%     figAxForm
%     hline(0,'k:')
%     vline(0,'k:')
%     title(sprintf('Visual Detect: R=%s,p=%s',...
%         num2str(round(c(1,2),2,'significant')),...
%         num2str(round(p(1,2),2,'significant'))))
%     subplot 242
%     x = cellInfo.detectComboWeight(ind);
%     y = cellInfo.detectWeight{auditoryTrials}(ind);
%     plot(x,y,'.','MarkerSize',10)
%     [c,p] = corrcoef(x,y);
%     hold on
%     figXAxis([],'Combo Weight',weightLim)
%     figYAxis([],'Aud. Only Weight',weightLim)
%     figAxForm
%     hline(0,'k:')
%     vline(0,'k:')
%     title(sprintf('Auditory Detect: R=%s,p=%s',...
%         num2str(round(c(1,2),2,'significant')),...
%         num2str(round(p(1,2),2,'significant'))))
%     subplot 243
%     x = cellInfo.detectWeight{visualTrials}(ind);
%     y = cellInfo.detectWeight{auditoryTrials}(ind);
%     plot(x,y,'.','MarkerSize',10)
%     [c,p] = corrcoef(x,y);
%     hold on
%     figXAxis([],'Vis. Only Weight',weightLim)
%     figYAxis([],'Aud. Only Weight',weightLim)
%     figAxForm
%     hline(0,'k:')
%     vline(0,'k:')
%     title(sprintf('Single Detect: R=%s,p=%s',...
%         num2str(round(c(1,2),2,'significant')),...
%         num2str(round(p(1,2),2,'significant'))))
%     subplot 244
%     y = abs(cat(1,cellInfo.detectWeight{visualTrials}(ind),...
%         cellInfo.detectWeight{auditoryTrials}(ind),...
%         cellInfo.detectComboWeight(ind)));
%     yerr = ste(y,2);
%     errorbar(1:3,mean(y,2),yerr,'.')
%     [p,~,stats] = anova1(y',[],'off');
%     c = multcompare(stats,[],'off');
%     disp(c(:,[1,2,6]))
%     figXAxis([],'Train',[0 4],1:3,{'Visual','Auditory','Vis+Aud'})
%     figYAxis([],'Weight Magnitude',[0 0.5])
%     figAxForm
%     title(sprintf('Detect Model, p=%s',num2str(round(p,2,'significant'))))
%     subplot 245
%     x = cellInfo.targetComboWeight(ind);
%     y = cellInfo.targetWeight{visualTrials}(ind);
%     plot(x,y,'.','MarkerSize',10)
%     [c,p] = corrcoef(x,y);
%     hold on
%     figXAxis([],'Combo Weight',weightLim)
%     figYAxis([],'Vis. Only Weight',weightLim)
%     figAxForm
%     hline(0,'k:')
%     vline(0,'k:')
%     title(sprintf('Visual Target: R=%s,p=%s',...
%         num2str(round(c(1,2),2,'significant')),...
%         num2str(round(p(1,2),2,'significant'))))
%     subplot 246
%     x = cellInfo.targetComboWeight(ind);
%     y = cellInfo.targetWeight{auditoryTrials}(ind);
%     plot(x,y,'.','MarkerSize',10)
%     [c,p] = corrcoef(x,y);
%     hold on
%     figXAxis([],'Combo Weight',weightLim)
%     figYAxis([],'Aud. Only Weight',weightLim)
%     figAxForm
%     hline(0,'k:')
%     vline(0,'k:')
%     title(sprintf('Auditory Target: R=%s,p=%s',...
%         num2str(round(c(1,2),2,'significant')),...
%         num2str(round(p(1,2),2,'significant'))))
%     subplot 247
%     x = cellInfo.targetWeight{visualTrials}(ind);
%     y = cellInfo.targetWeight{auditoryTrials}(ind);
%     plot(x,y,'.','MarkerSize',10)
%     [c,p] = corrcoef(x,y);
%     hold on
%     figXAxis([],'Vis. Only Weight',weightLim)
%     figYAxis([],'Aud. Only Weight',weightLim)
%     figAxForm
%     hline(0,'k:')
%     vline(0,'k:')
%     title(sprintf('Single Target: R=%s,p=%s',...
%         num2str(round(c(1,2),2,'significant')),...
%         num2str(round(p(1,2),2,'significant'))))
%     subplot 248
%     y = abs(cat(1,cellInfo.targetWeight{visualTrials}(ind),...
%         cellInfo.targetWeight{auditoryTrials}(ind),...
%         cellInfo.targetComboWeight(ind)));
%     yerr = ste(y,2);
%     errorbar(1:3,mean(y,2),yerr,'.')
%     [p,~,stats] = anova1(y',[],'off');
%     c = multcompare(stats,[],'off');
%     disp(c(:,[1,2,6]))
%     figXAxis([],'Train',[0 4],1:3,{'Visual','Auditory','Vis+Aud'})
%     figYAxis([],'Weight Magnitude',[0 0.5])
%     figAxForm
%     title(sprintf('Target Model, p=%s',num2str(round(p,2,'significant'))))
%     print([fnout 'decodeModelComboTrainWeights'],'-dpdf','-fillpage')
% 
% %% stats for modeling
% avName = {'Vis','Aud'};
% for iav = 1:2
%     nav = avName{iav};
%     if iav == 2
%         stimID = [1,3];
%     else
%         stimID = 1:3;
%     end
%     if iav == 1
%         [p,~,stats] = anova1(dcModel(1).av(iav).pctCorrect(stimID,:)',[],'off');
%         tbl = multcompare(stats,[],'off');
%         tbl = tbl(:,[1:2,end]);
%     else
%         y = dcModel(1).av(iav).pctCorrect(stimID,:)';
%         [~,p] = ttest(y(:,1),y(:,2));
%         tbl = [];
%     end
%     imgStats.av(iav).targetPerformanceXStimAnova = p;
%     imgStats.av(iav).targetPerformanceXStimPostHocTests = tbl;
%     fprintf('%s Target Performance x Stim ANOVA, p=%s\n',nav,num2str(round(...
%         imgStats.av(iav).targetPerformanceXStimAnova,2,'significant')))
%     disp(imgStats.av(iav).targetPerformanceXStimPostHocTests)
%     
%     if iav == 1
%         [p,~,stats] = anova1(dcModel(2).av(iav).pctCorrect(stimID,:)',[],'off');
%         tbl = multcompare(stats,[],'off');
%         tbl = tbl(:,[1:2,end]);
%     else
%         y = dcModel(2).av(iav).pctCorrect(stimID,:)';
%         [~,p] = ttest(y(:,1),y(:,2));
%         tbl = [];
%     end
%     imgStats.av(iav).detectPerformanceXStimAnova = p;
%     imgStats.av(iav).detectPerformanceXStimPostHocTests = tbl;
%     fprintf('%s Detect Performance x Stim ANOVA, p=%s\n',nav,num2str(round(...
%         imgStats.av(iav).detectPerformanceXStimAnova,2,'significant')))
%     disp(imgStats.av(iav).detectPerformanceXStimPostHocTests)
% end
% 
% y = [];
% z = [];
% for iexp = 1:nexp
%     if ~isempty(decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectTarget)
%         y = cat(1,y,[decodeAnalysis(iexp).av(visualTrials).validMatchedPctCorrectTarget_holdout,...
%             decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectTarget]);
%     z = cat(1,z,decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectTarget_testAudModel);
%     end
% end
% [~,p] = ttest(y(:,1),y(:,2));
% imgStats.visInvVsVisValTargetTest = p;
% [~,p] = ttest(z,0.5);
% imgStats.visInvTestAudTarget = p;
% fprintf('Paired t-Test Vis. Val vs. Vis. Inv. in Vis Target Model,p=%s\n',...
%     num2str(round(imgStats.visInvVsVisValTargetTest,2,'significant')))
% fprintf('t-Test Vis. Inv. vs. 0.5 in Aud Target Model,p=%s\n',...
%     num2str(round(imgStats.visInvTestAudTarget,2,'significant')))
% 
% y = [];
% z = [];
% a = [];
% for iexp = 1:nexp
%     if ~isempty(decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectDetect)
%         y = cat(1,y,[decodeAnalysis(iexp).av(visualTrials).validMatchedPctCorrectDetect_holdout,...
%             decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectDetect]);
%         z = cat(1,z,decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectDetect_testAudModel);
%         a = cat(1,a,decodeAnalysis(iexp).av(visualTrials).validMatchedPctCorrectDetect_testAudModel);
%     end
% end
% [~,p] = ttest(y(:,1),y(:,2));
% imgStats.visInvVsVisValDetectTest = p;
% [~,p] = ttest(z,y(:,1));
% imgStats.visInvTestAudVsVisInvTestVisDetect = p;
% [~,p] = ttest(y(:,1),z);
% imgStats.visInvTestAudVsVisValTestVisDetect = p;
% fprintf('Paired t-Test Vis. Val vs. Vis. Inv. in Vis Detect Model,p=%s\n',...
%     num2str(round(imgStats.visInvVsVisValDetectTest,2,'significant')))
% fprintf('Paired t-Test Vis. Val. in Vis Detect vs. Vis. Inv. in Aud Detect Model,p=%s\n',...
%     num2str(round(imgStats.visInvTestAudVsVisValTestVisDetect,2,'significant')))
% 
% save([fnout 'imgStats'],'imgStats')
% end
% %% Naive decoding figure
% 
% if ~strcmp(ds,'FSAV_attentionV1')
%     
%     dcModel = struct;
%     dcModel(1).cellTypeName = {'Tar.';'Dist.';'N.P.'};
%     dcModel(1).name = 'Target';
% %     cellTypeIDAllCells = sum(cat(2,(cellInfo.lateCycRespCells & ~cellInfo.targetRespCells).*1,...
% %         (cellInfo.lateCycRespCells & cellInfo.targetRespCells).*2,...
% %         (~cellInfo.lateCycRespCells & cellInfo.targetRespCells).*3),2);
% %     dcModel(1).cellTypeID = cellTypeIDAllCells(cellInfo.dcModelCells);
%     dcModel.av(1).name = 'Visual';
%     dcModel.av(2).name = 'Auditory';
%     for iav = 1:2
%         dcModel.av(iav).dv = nan(1,nexp);
%         dcModel.av(iav).correlation = [];
%         dcModel.av(iav).weight = [];
%         dcModel.av(iav).pctCorrect = nan(4,nexp);
%         dcModel.av(iav).testPerformance = nan(4,nexp);
%         dcModel.av(iav).pctCorrect_otherAV = nan(1,nexp);
%     end
%     for iav = 1:2
%         for iexp = 1:nexp
%             d = decodeAnalysis(iexp).av(iav);
%             dcModel(1).av(iav).dv(iexp) = d.dvTarget;
% 
%             dcModel(1).av(iav).correlation = cat(2, ...
%                 dcModel(1).av(iav).correlation,d.correlationTarget);
%             dcModel(1).av(iav).weight = cat(1, ...
%                 dcModel(1).av(iav).weight,d.weightTarget);
%             dcModel(1).av(iav).pctCorrect(:,iexp) = cat(2,...
%                 d.pctCorrectXStimTarget_holdout,d.pctCorrectAllTarget_holdout);
%             dcModel(1).av(iav).testPerformance(:,iexp) = cat(2,...
%                 d.pctCorrectXStimTarget_train,d.pctCorrectAllTarget_train) - ...
%                 cat(2,...
%                 d.pctCorrectXStimTarget_holdout,d.pctCorrectAllTarget_holdout);
%             dcModel(1).av(iav).pctCorrect_otherAV(iexp) = d.pctCorrectTarget_otherAV;
%         end
%     end
% 
%     stimLabel = {'0';'HT';'ET';'All'};
%     mdlXStimAlpha = 0.05./3;
%     
%         iexp=1;
%     for iav = 1:2
%         dc = decodeDataExpt(iexp).av(iav);
%         ind = decodeAnalysis(iexp).cellInd;
%     %     trialInd = dcTrialsUsed{iexp,iav};
%     %     trOut = dc.outcome(trialInd);
%         trOut = dc.outcome;
% 
%     %     r = dc.tc(:,ind,trialInd);
%         r = dc.tc(:,ind,:);
%         [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOut);
%         wTarget = decodeAnalysis(iexp).av(iav).weightTarget;
%         [~,targetWeightSort] = sort(wTarget);
% 
%         rTarget = {mean(r(:,:,targetTrInd==0),3);mean(r(:,:,targetTrInd==1),3)};
%         rTargetErr = {ste(r(:,:,targetTrInd==0),3);ste(r(:,:,targetTrInd==1),3)};
% 
%         tfig_d = figure;
%         suptitle(sprintf('%s, %s, Dist. Resp, Target Weight Sorted',...
%             decodeDataExpt(iexp).exptName,avName{iav}))
%         tfig_t = figure;
%         suptitle(sprintf('%s, %s, Tar. Resp, Target Weight Sorted',...
%             decodeDataExpt(iexp).exptName,avName{iav}))
%         for icell = 1:sum(ind)
%             figure(tfig_d);
%             subplot(5,3,icell)
%             cInd = targetWeightSort(icell);
%     %         plot(tt_cycTC,rTarget{1}(tcStartFrame:nFr_cyc,cInd),'k-')
%             shadedErrorBar_chooseColor(tt_cycTC,rTarget{1}(tcStartFrame:nFr_cyc,cInd),...
%                 rTargetErr{1}(tcStartFrame:nFr_cyc,cInd),[0 0 0]);
%             figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],[0,200],[0,200])
%             figYAxis([],'dF/F',cellRespTCLim)  
%             hline(0,'r:')
%             vline(respWinTT,'k--')
%             figAxForm
%             title(sprintf('D Resp., Cell #%s, W=%s',num2str(cInd),...
%                 num2str(round(wTarget(cInd),2,'significant'))))
%             figure(tfig_t);
%             subplot(5,3,icell)
%     %         plot(tt_cycTC,rTarget{2}(tcStartFrame:nFr_cyc,cInd),'k-')
%             shadedErrorBar_chooseColor(tt_cycTC,rTarget{2}(tcStartFrame:nFr_cyc,cInd),...
%                 rTargetErr{2}(tcStartFrame:nFr_cyc,cInd),[0 0 0]);
%             figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],[0,200],[0,200])
%             figYAxis([],'dF/F',cellRespTCLim)  
%             hline(0,'r:')
%             vline(respWinTT,'k--')
%             figAxForm
%             title(sprintf('T Resp., Cell #%s, W=%s',num2str(cInd),...
%                 num2str(round(wTarget(cInd),2,'significant'))))
%         end
%         figure(tfig_d);
%         print(...
%             [fnout 'targetModel_' avName{iav}(1:3) '_distResp_' decodeDataExpt(iexp).exptName],...
%             '-dpdf','-fillpage')
%         figure(tfig_t);
%         print(...
%             [fnout 'targetModel_' avName{iav}(1:3) '_tarResp_' decodeDataExpt(iexp).exptName],...
%             '-dpdf','-fillpage')
%     end
% 
%     
%     figure
%     for iav = 1:2
%         y = dcModel.av(iav).pctCorrect;
%         yerr = ste(dcModel.av(iav).pctCorrect,2);
%         stimInd = sum(~isnan(y),2) > minTrN;
%         y = y(stimInd,:);
% %         yerr = yerr(stimInd,:);
%         x = 1:sum(stimInd);
%         [mdlTest,mdlP] = ttest(y,0.5,'dim',2,'tail','right','alpha',mdlXStimAlpha);
%         subplot(1,2,iav)
%         for iexp = 1:nexp
%             ind = ~isnan(y(:,iexp));
%             hold on
%             plot(x(ind),y(ind,iexp),'k.-');
%         end
%         for itest = 1:length(x)
%             if mdlTest(itest) == 1
%                 text(x(itest),1,sprintf('*%s',...
%                     num2str(round(mdlP(itest),2,'significant'))))
%             else
%                 text(x(itest),1,sprintf('%s',...
%                     num2str(round(mdlP(itest),2,'significant'))))
%             end
%         end
%         errorbar(x,nanmean(y,2),yerr,'.-','MarkerSize',15)
%         figXAxis([],'Stim',[0 length(x)+1],1:length(x),stimLabel(stimInd))
%         hline(0.5,'k:')
%         figYAxis([],'Train - HO Frac. Correct',[0 1])
%         figAxForm
%         title(sprintf('%s Model, Train %s, Test %s',dcModel.name,...
%             dcModel.av(iav).name,dcModel.av(iav).name))
%     end
%     
%     print([fnout 'decodeModelPctCorrectXStim'],'-dpdf','-fillpage')
%     
%     figure
%     ind = cellInfo.lateCycRespCells & ~isnan(cellInfo.targetWeight{visualTrials})';
%     subplot 221
%     x = antiAnalysis.lateCycSI(ind)';
%     y = cellInfo.targetWeight{visualTrials}(ind)';
%     hold on
%     plot(x,y,'.')
%     [corrmat,pmat] = corrcoef(x,y);
%     c = corrmat(1,2);
%     p = pmat(1,2);
%     mdl = fitlm(x,y);
%     rsq = mdl.Rsquared.Ordinary;
%     yfit = predict(mdl,x);
%     plot(x,yfit,'-')
%     figXAxis([],'Selectivity',siLim)
%     figYAxis([],'Visual Weight',weightLim)
%     figAxForm
%     title(sprintf('Visual Target Model, Rsq = %s, corr = %s, p = %s',num2str(round(rsq,2,'significant')),...
%         num2str(round(c,2,'significant')),num2str(round(p,2,'significant'))))
%     subplot 222
%     y = cellInfo.targetWeight{auditoryTrials}(ind)';
%     hold on
%     plot(x,y,'.')
%     [corrmat,pmat] = corrcoef(x,y);
%     c = corrmat(1,2);
%     p = pmat(1,2);
%     mdl = fitlm(x,y);
%     rsq = mdl.Rsquared.Ordinary;
%     yfit = predict(mdl,x);
%     plot(x,yfit,'-')
%     figXAxis([],'Selectivity',siLim)
%     figYAxis([],'Auditory Weight',weightLim)
%     figAxForm
%     title(sprintf('Auditory Target Model, Rsq = %s, corr = %s, p = %s',num2str(round(rsq,2,'significant')),...
%         num2str(round(c,2,'significant')),num2str(round(p,2,'significant'))))
%     subplot 223
%     ind = cellInfo.lateCycRespCells & ~isnan(cellInfo.targetWeight{visualTrials})';
%     x = cellInfo.targetWeight{visualTrials}(ind)';
%     y = cellInfo.targetWeight{auditoryTrials}(ind)';
%     hold on
%     plot(x,y,'.')
%     [corrmat,pmat] = corrcoef(x,y);
%     c = corrmat(1,2);
%     p = pmat(1,2);
%     mdl = fitlm(x,y);
%     rsq = mdl.Rsquared.Ordinary;
%     yfit = predict(mdl,x);
%     plot(x,yfit,'-')
%     figXAxis([],'Visual Weight',weightLim)
%     figYAxis([],'Auditory Weight',weightLim)
%     figAxForm
%     title(sprintf('Target Model, Rsq = %s, corr = %s, p = %s',num2str(round(rsq,2,'significant')),...
%         num2str(round(c,2,'significant')),num2str(round(p,2,'significant'))))
%     
%     print([fnout 'decodeModelWeights'],'-dpdf','-fillpage')
% end