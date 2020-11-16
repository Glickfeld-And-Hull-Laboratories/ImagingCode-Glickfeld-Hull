clear all
close all
ds = 'FSAV_V1_naive_GCaMP6m'; % 'FSAV_V1_naive_GCaMP6m'  'FSAV_attentionV1'   'FSAV_attentionV1_noAttn'
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
        [titleStr '_rvs_' datestr(now,'yymmdd') '_']); 
elseif cellsOrDendrites == 2
    load(fullfile(rc.caOutputDir,ds,...
        [mouse_str '_trOutcomeStruct_d' ds(5:end) '.mat']));
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
                        
            if strcmp(ds,'FSAV_attentionV1')|strcmp(ds,'FSAV_attentionV1_noAttn')
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
            if length(d.av(1).align) == 5 &  isfield(de,'catchOutcome')
                iav = 1;
                de = d.av(iav).align(5);
                outs = de.outcome;
                outs(strcmp(outs,'failure')) = {'fa'};
                outs(strcmp(outs,'success')) = {'h'};
                outs(strcmp(outs,'ignore')) = {'m'};
                outcomeEaTrial = outs;
                if isfield(de,'catchOutcome')
                    catch_outs = de.catchOutcome;
                    catch_outs(cellfun(@isempty,catch_outs)) = {nan};
                    catch_outs(strcmp(catch_outs,'FA')) = {'h'};
                    catch_outs(strcmp(catch_outs,'CR')) = {'m'};
                    catchOutcomeEaTrial = catch_outs;
                    catchOriEaTr = de.catchOri;
                    catchOriEaTr(cellfun(@isnan,catch_outs)) = nan;
                end                
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
                
                if isfield(de,'catchOutcome')
                    decodeDataExpt(exptN).catch.cycRespEaTrial = cycRespEaTrial;
                    decodeDataExpt(exptN).catch.catchOutcomeEaTrial = catchOutcomeEaTrial;
                    decodeDataExpt(exptN).catch.catchStimEaTrial = catchOriEaTr;
                    decodeDataExpt(exptN).catch.outcomeEaTrial = outcomeEaTrial;
                    decodeDataExpt(exptN).catch.stimEaTrial = stimEaTrial;
                    decodeDataExpt(exptN).catch.catchCycN = cCyc;
                end
                
            end
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
                elseif strcmp(ds,'FSAV_V1_naive_GCaMP6m')
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
                if isfield(d.av(visualTrials).align(5),'catchOutcome')
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
            [~,tuningReliabilitySortInd] = sort(oriTuningExpt(iexp).tuningReliability');
            responsiveCells = respCellsExpt(iexp).lateCycRespCells|...
                respCellsExpt(iexp).targetRespCells;
            tunedCells = tuningReliabilitySortInd(responsiveCells);
            cellInd = false(length(responsiveCells),1);
%             cellInd = respCellsExpt(iexp).decodeAnalysisCells & ...
%                 oriTuningExpt(iexp).tuningReliability' <= tuningReliabilityThresh_decode;
%             cellInd = respCellsExpt(iexp).lateCycRespCells|respCellsExpt(iexp).targetRespCells;
            decodeAnalysis(iexp).nCells = length(tunedCells);            
            
            respOther = cell(1,2);
            trOutOther = cell(1,2);
            trSampleIndOther = cell(1,2);
            detectGLMOther = cell(1,2);
            targetGLMOther = cell(1,2);
            respOther_dist = cell(1,2);
            detectGLMOther_dist = cell(1,2);
            trOutOther_dist = cell(1,2);
            respOther_pcs = cell(1,2);
            trOutOther_pcs = cell(1,2);
            trSampleIndOther_pcs = cell(1,2);
            detectGLMOther_pcs = cell(1,2);
            targetGLMOther_pcs = cell(1,2);
            respOther_pcs_dist = cell(1,2);
            detectGLMOther_pcs_dist = cell(1,2);
            trOutOther_pcs_dist = cell(1,2);
            
            rng(0)
            for iav = 1:2
                trOut = decodeDataExpt(iexp).av(iav).outcome;

                if iav == 1
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
                elseif iav == 2
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,ampBins);
                end
                nStimPerBin = histcounts(trStimID);
                minBinN = min(nStimPerBin(nStimPerBin >= minTrN_mdl));
                trials_StimSort = cell(1,nStimBins);
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
                        trials_StimSort{istim} = indSample;
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
                decodeAnalysis(iexp).av(iav).trOut = trOut(matchTrialsInd);
                
                trials_StimSort_av{iav} = trials_StimSort;
                trOutStimSort_av{iav} = trOutStimSort;
                stimSortInd_av{iav} = stimSortInd;
                matchTrialsInd_av{iav} = matchTrialsInd;
            end
            
            rng(0)
            nTrials = min(cellfun(@length,matchTrialsInd_av));
            if round(nTrials.*fracCellsUsed) > maxCellN
                nCells = maxCellN;
            elseif length(cellInd) < round(nTrials.*fracCellsUsed)
                nCells = length(cellInd);
            elseif round(nTrials.*fracCellsUsed) < minCellN
                nCells = minCellN;
            else
                nCells = round(nTrials.*fracCellsUsed);
            end
            cellInd(tuningReliabilitySortInd(1:nCells)) = true;
            
            if nCells ~= sum(cellInd)
                ind = randsample(find(cellInd),nCells);
                cellInd = false(length(cellInd),1);
                cellInd(ind) = true;
            end
            decodeAnalysis(iexp).cellInd = cellInd;
            decodeAnalysis(iexp).nCellsUsed = nCells;
            
            nPCs = nCells;
            
            decodeAnalysis(iexp).nPCs = nPCs;
            decodeAnalysis(iexp).minTrials = nTrials;   
            
            data = cat(2,...
                decodeDataExpt(iexp).av(visualTrials).resp(cellInd,matchTrialsInd_av{visualTrials}),...
                decodeDataExpt(iexp).av(auditoryTrials).resp(cellInd,matchTrialsInd_av{auditoryTrials}))';
            nvis = length(matchTrialsInd_av{visualTrials});
            naud = length(matchTrialsInd_av{auditoryTrials});
            
            data_z = zscore(data);
            respAllCells_vis = data_z(1:nvis,:);
            respAllCells_aud = data_z((nvis+1):end,:);

            [coeffAllCells,scoresAllCells] = pca(data);
            respAllCells_pcs_av = zscore(scoresAllCells);
            respAllCells_pcs_vis = respAllCells_pcs_av(1:nvis,1:nPCs);
            respAllCells_pcs_aud = respAllCells_pcs_av((nvis+1):end,1:nPCs);
            
            
            data_dist = cat(2,...
                decodeDataExpt(iexp).av(visualTrials).resp(cellInd,...
                strcmp(decodeDataExpt(iexp).av(visualTrials).outcome,'fa')|...
                strcmp(decodeDataExpt(iexp).av(visualTrials).outcome,'cr')),...
                decodeDataExpt(iexp).av(auditoryTrials).resp(cellInd,...
                strcmp(decodeDataExpt(iexp).av(auditoryTrials).outcome,'fa')|...
                strcmp(decodeDataExpt(iexp).av(auditoryTrials).outcome,'cr')))';
            nvis = sum(strcmp(decodeDataExpt(iexp).av(visualTrials).outcome,'fa')|...
                strcmp(decodeDataExpt(iexp).av(visualTrials).outcome,'cr'));
            naud = sum(strcmp(decodeDataExpt(iexp).av(auditoryTrials).outcome,'fa')|...
                strcmp(decodeDataExpt(iexp).av(auditoryTrials).outcome,'cr'));
            data_dist_z = zscore(data_dist);
            respAllCells_dist_vis = data_dist_z(1:nvis,:);
            respAllCells_dist_aud = data_dist_z((nvis+1):end,:);
            
            [coeffAllCells_dist,scoresAllCells_dist] = pca(data_dist);
            respAllCells_dist_pcs_av = zscore(scoresAllCells_dist);
            respAllCells_dist_pcs_vis = respAllCells_dist_pcs_av(1:nvis,1:nPCs);
            respAllCells_dist_pcs_aud = respAllCells_dist_pcs_av((nvis+1):end,1:nPCs);
            
            respStimSort_av = cell(1,2);
            for iav = 1:2
                if iav == 1
                    respAllCells = respAllCells_vis;
                else
                    respAllCells = respAllCells_aud;
                end
                decodeAnalysis(iexp).av(iav).respAllCells = respAllCells;
                    n = [0,cumsum(cellfun(@length,trials_StimSort_av{iav}))];
                    respStimSort_av{iav} = cell(1,3);
                for i = 1:3
                    if ~isempty(trials_StimSort_av{iav}{i})
                        ind = (n(i)+1):(n(i+1));
                        respStimSort_av{iav}{i} = respAllCells(ind,:);
                    end
                end
            end
          
            
            for iav = 1:2
                trOut = decodeDataExpt(iexp).av(iav).outcome;
                respStimSort = respStimSort_av{iav};
                trOutStimSort = trOutStimSort_av{iav};
                stimSortInd = stimSortInd_av{iav};
                matchTrialsInd = matchTrialsInd_av{iav};

                if iav == 1
                    resp = respAllCells_vis;
                    resp_pcs = respAllCells_pcs_vis;
                elseif iav == 2
                    resp = respAllCells_aud;
                    resp_pcs = respAllCells_pcs_aud;
                end
                
                emptyInd = cellfun(@isempty,respStimSort);
                respStimSort(~emptyInd) = cellfun(@(x) x(:,1:nPCs),...
                    respStimSort(~emptyInd),'unif',0);
                
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
                 
                % test model with ea stim
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
                
                % response window analysis
                fprintf('Expt %s, starting resp win analysis\n',num2str(iexp))
                nwins = size(decodeDataExpt(iexp).av(iav).movWinResp,3);
                pctCorrectDetect_train_respwin = nan(1,nwins);
                pctCorrectDetect_ho_respwin = nan(1,nwins);
                pctCorrectTarget_train_respwin = nan(1,nwins);
                pctCorrectTarget_ho_respwin = nan(1,nwins);
                for iwin = 1:nwins
%                     rAllCells = zscore(decodeDataExpt(iexp).av(iav).movWinResp(:,:,iwin)');
%                     r = rAllCells(matchTrialsInd,cellInd);
                    rAllCells = decodeDataExpt(iexp).av(iav).movWinResp(cellInd,matchTrialsInd,iwin)';
%                     [~,sAllCells] = pca(rAllCells(:,cellInd));
                    r = zscore(rAllCells);
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
                    
                % model with pcs                
                C = eye(size(resp_pcs,2));
                p=1;
                [~,~,detectGLM] = glmfit(resp_pcs*C,detectTrInd,'binomial');
                [~,~,targetGLM] = glmfit(resp_pcs*C,targetTrInd,'binomial');
                                
                detectWeight_pcs = detectGLM.beta(2:end);
                targetWeight_pcs = targetGLM.beta(2:end);
                
                fprintf('Expt %s, hold-out analysis\n',num2str(iexp))
                pctCorrectDetect_pcs_ho = getPctCorr_hoData(resp_pcs,detectTrInd,dv_detect);

                pctCorrectTarget_pcs_ho = getPctCorr_hoData(resp_pcs,targetTrInd,dv_target);
 
                detectWeight_pc2neur = coeffAllCells(:,1:nPCs)*detectGLM.beta(2:end);
                targetWeight_pc2neur = coeffAllCells(:,1:nPCs)*targetGLM.beta(2:end);
                
                respOther_pcs{iav} = resp_pcs;
                detectGLMOther_pcs{iav} = detectGLM;
                targetGLMOther_pcs{iav} = targetGLM; 
                
%                 % train detect model with distractor choices only
                fprintf('Expt %s, distractor only model\n',num2str(iexp))
                % distractors
%                 distInd = strcmp(trOut_matched,'fa')|strcmp(trOut_matched,'cr');
%                 resp_dist = resp(distInd,:);
%                 trOut_matched = trOut(matchTrialsInd);
%                 trOut_dist = trOut_matched(distInd);
                if iav == 1
                    resp_dist = respAllCells_dist_vis;
                    trOut_dist = trOut(strcmp(decodeDataExpt(iexp).av(visualTrials).outcome,'fa')|...
                        strcmp(decodeDataExpt(iexp).av(visualTrials).outcome,'cr'));
                else
                    resp_dist = respAllCells_dist_aud;
                    trOut_dist = trOut(strcmp(decodeDataExpt(iexp).av(auditoryTrials).outcome,'fa')|...
                        strcmp(decodeDataExpt(iexp).av(auditoryTrials).outcome,'cr'));
                end
                [detectTrInd_dist, ~] = getStimAndBehaviorYs(trOut_dist);
                dv_detect_dist = mean(detectTrInd_dist);
                C = eye(size(resp_dist,2));
                p=1;
                [~,~,detectGLM_dist] = glmfit(resp_dist*C,detectTrInd_dist,'binomial');
                pctCorrectDetect_ho_dist = getPctCorr_hoData(resp_dist,detectTrInd_dist,dv_detect_dist);
                respOther_dist{iav} = resp_dist;
                detectGLMOther_dist{iav} = detectGLM_dist;
                trOutOther_dist{iav} = trOut_dist;
                detectWeight_cells_distOnly = detectGLM_dist.beta(2:end);
                
                % distractor model with pcs
                if iav == 1
                    resp_dist = respAllCells_dist_pcs_vis;
                else
                    resp_dist = respAllCells_dist_pcs_aud;
                end
                [~,~,detectGLM_dist] = glmfit(resp_dist*C,detectTrInd_dist,'binomial');
                pctCorrectDetect_pcs_ho_dist = getPctCorr_hoData(resp_dist,detectTrInd_dist,dv_detect_dist);
                respOther_pcs_dist{iav} = resp_dist;
                detectGLMOther_pcs_dist{iav} = detectGLM_dist;
                detectWeight_pcs_distOnly = detectGLM_dist.beta(2:end);
                
                detectWeight_distOnly_pc2neur = coeffAllCells(:,1:nPCs)*detectGLM.beta(2:end);
                
                decodeAnalysis(iexp).av(iav).dvDetect = dv_detect;
                decodeAnalysis(iexp).av(iav).dvTarget = dv_target;
                decodeAnalysis(iexp).av(iav).correlationDetect = detectCorr;
                decodeAnalysis(iexp).av(iav).correlationTarget = targetCorr;
                decodeAnalysis(iexp).av(iav).weightDetect = detectWeight;
                decodeAnalysis(iexp).av(iav).weightTarget = targetWeight;
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
                
                decodeAnalysis(iexp).av(iav).weightDetect_pcs = detectWeight_pcs;
                decodeAnalysis(iexp).av(iav).weightTarget_pcs = targetWeight_pcs;
                decodeAnalysis(iexp).av(iav).pctCorrectAllDetect_pcs_holdout = pctCorrectDetect_pcs_ho;
                decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_pcs_holdout = pctCorrectTarget_pcs_ho;
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_pcs_distOnly_holdout = pctCorrectDetect_pcs_ho_dist;
                decodeAnalysis(iexp).av(iav).weightDetect_pcs_distOnly = detectWeight_pcs_distOnly;
                
            end            
            fprintf('Expt %s: Testing opposite model...\n',num2str(iexp))
            for iav = 1:2
                if iav == 1
                    otherAV = 2;
                elseif iav == 2
                    otherAV = 1;
                end
                resp = respOther{otherAV};
                [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{otherAV});

                [detecttrind4model,targettrind4model] = getStimAndBehaviorYs(trOutOther{iav});
                dv_detect = mean(detecttrind4model);
                dv_target = mean(targettrind4model);
                
                pctCorrectDetect = getPctCorr_trainData(detectGLMOther{iav},resp,detectTrInd,dv_detect);
                pctCorrectTarget = getPctCorr_trainData(targetGLMOther{iav},resp,targetTrInd,dv_target);
                                
                distInd = targetTrInd == 0;
                pctCorrectDetect_test_zero = getPctCorr_trainData(detectGLMOther{iav},resp(distInd,:),detectTrInd(distInd),dv_detect);
                
                % test other pc models
                resp = respOther_pcs{otherAV};
                pctCorrectDetect_pcs = getPctCorr_trainData(detectGLMOther_pcs{iav},resp,detectTrInd,dv_detect);
                
                
                % test other distractor only detect model
                resp = respOther_dist{otherAV};
                [detectTrInd, ~] = getStimAndBehaviorYs(trOutOther_dist{otherAV});
                dv_detect = mean(getStimAndBehaviorYs(trOutOther_dist{iav}));
                pctCorrectDetect_dist = getPctCorr_trainData(detectGLMOther_dist{iav},resp,detectTrInd,dv_detect);
                
                resp = respOther_pcs_dist{otherAV};
                pctCorrectDetect_pcs_dist = getPctCorr_trainData(detectGLMOther_pcs_dist{iav},resp,detectTrInd,dv_detect);
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_otherAV = pctCorrectDetect;
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_otherAV = pctCorrectTarget;
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_onlyZero_otherAV = pctCorrectDetect_test_zero;
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_pcs_otherAV = pctCorrectDetect_pcs;
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_distOnly_otherAV = pctCorrectDetect_dist;
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_pcs_distOnly_otherAV = pctCorrectDetect_pcs_dist;
                
            end
            
        end

    elseif strcmp(ds,'FSAV_V1_naive_GCaMP6m')
        
        decodeAnalysis = struct;
        decodeAnalysis.av(visualTrials).name = 'Visual';
        decodeAnalysis.av(auditoryTrials).name = 'Auditory';
        for iexp = 1:nexp
            rng(0) % cells randomized and trials randomized
            [~,tuningReliabilitySortInd] = sort(oriTuningExpt(iexp).tuningReliability');
            responsiveCells = respCellsExpt(iexp).lateCycRespCells|...
                respCellsExpt(iexp).targetRespCells;
            tunedCells = tuningReliabilitySortInd(responsiveCells);
            cellInd = false(length(responsiveCells),1);
%             cellInd = respCellsExpt(iexp).decodeAnalysisCells & ...
%                 oriTuningExpt(iexp).tuningReliability' <= tuningReliabilityThresh_decode;
%             cellInd = respCellsExpt(iexp).lateCycRespCells|respCellsExpt(iexp).targetRespCells;
            decodeAnalysis(iexp).nCells = length(tunedCells);            
            
            respOther = cell(1,2);
            trOutOther = cell(1,2);
            trSampleIndOther = cell(1,2);
            targetGLMOther = cell(1,2);
            respOther_pcs = cell(1,2);
            trOutOther_pcs = cell(1,2);
            trSampleIndOther_pcs = cell(1,2);
            targetGLMOther_pcs = cell(1,2);
            respOther_pcs_dist = cell(1,2);
            trOutOther_pcs_dist = cell(1,2);
            
            rng(0)
            for iav = 1:2
                trOut = decodeDataExpt(iexp).av(iav).outcome;

                if iav == 1
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
                elseif iav == 2
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,ampBins);
                end
                nStimPerBin = histcounts(trStimID);
                minBinN = min(nStimPerBin(nStimPerBin >= minTrN_mdl));
                trials_StimSort = cell(1,nStimBins);
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
                        trials_StimSort{istim} = indSample;
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
                decodeAnalysis(iexp).av(iav).trOut = trOut(matchTrialsInd);
                
                trials_StimSort_av{iav} = trials_StimSort;
                trOutStimSort_av{iav} = trOutStimSort;
                stimSortInd_av{iav} = stimSortInd;
                matchTrialsInd_av{iav} = matchTrialsInd;
            end
            
            rng(0)
            nTrials = min(cellfun(@length,matchTrialsInd_av));
            if length(cellInd) < round(nTrials.*fracCellsUsed)
                nCells = length(cellInd);
            elseif round(nTrials.*fracCellsUsed) > maxCellN
                nCells = maxCellN;
            elseif round(nTrials.*fracCellsUsed) < minCellN
                nCells = minCellN;
            else
                nCells = round(nTrials.*fracCellsUsed);
            end
            cellInd(tuningReliabilitySortInd(1:nCells)) = true;
            
            if nCells ~= sum(cellInd)
                ind = randsample(find(cellInd),nCells);
                cellInd = false(length(cellInd),1);
                cellInd(ind) = true;
            end
            decodeAnalysis(iexp).cellInd = cellInd;
            decodeAnalysis(iexp).nCellsUsed = nCells;
            
            nPCs = nCells;
            
            decodeAnalysis(iexp).nPCs = nPCs;
            decodeAnalysis(iexp).minTrials = nTrials;   
            
            data = cat(2,...
                decodeDataExpt(iexp).av(visualTrials).resp(cellInd,matchTrialsInd_av{visualTrials}),...
                decodeDataExpt(iexp).av(auditoryTrials).resp(cellInd,matchTrialsInd_av{auditoryTrials}))';
            nvis = length(matchTrialsInd_av{visualTrials});
            naud = length(matchTrialsInd_av{auditoryTrials});
            
            data_z = zscore(data);
            respAllCells_vis = data_z(1:nvis,:);
            respAllCells_aud = data_z((nvis+1):end,:);

            [coeffAllCells,scoresAllCells] = pca(data);
            respAllCells_pcs_av = zscore(scoresAllCells);
            respAllCells_pcs_vis = respAllCells_pcs_av(1:nvis,1:nPCs);
            respAllCells_pcs_aud = respAllCells_pcs_av((nvis+1):end,1:nPCs);
            
            
            data_dist = cat(2,...
                decodeDataExpt(iexp).av(visualTrials).resp(cellInd,...
                strcmp(decodeDataExpt(iexp).av(visualTrials).outcome,'fa')|...
                strcmp(decodeDataExpt(iexp).av(visualTrials).outcome,'cr')),...
                decodeDataExpt(iexp).av(auditoryTrials).resp(cellInd,...
                strcmp(decodeDataExpt(iexp).av(auditoryTrials).outcome,'fa')|...
                strcmp(decodeDataExpt(iexp).av(auditoryTrials).outcome,'cr')))';
            nvis = sum(strcmp(decodeDataExpt(iexp).av(visualTrials).outcome,'fa')|...
                strcmp(decodeDataExpt(iexp).av(visualTrials).outcome,'cr'));
            naud = sum(strcmp(decodeDataExpt(iexp).av(auditoryTrials).outcome,'fa')|...
                strcmp(decodeDataExpt(iexp).av(auditoryTrials).outcome,'cr'));
            data_dist_z = zscore(data_dist);
            respAllCells_dist_vis = data_dist_z(1:nvis,:);
            respAllCells_dist_aud = data_dist_z((nvis+1):end,:);
            
            [coeffAllCells_dist,scoresAllCells_dist] = pca(data_dist);
            respAllCells_dist_pcs_av = zscore(scoresAllCells_dist);
            respAllCells_dist_pcs_vis = respAllCells_dist_pcs_av(1:nvis,1:nPCs);
            respAllCells_dist_pcs_aud = respAllCells_dist_pcs_av((nvis+1):end,1:nPCs);
            
            respStimSort_av = cell(1,2);
            for iav = 1:2
                if iav == 1
                    respAllCells = respAllCells_vis;
                else
                    respAllCells = respAllCells_aud;
                end
                decodeAnalysis(iexp).av(iav).respAllCells = respAllCells;
                    n = [0,cumsum(cellfun(@length,trials_StimSort_av{iav}))];
                    respStimSort_av{iav} = cell(1,3);
                for i = 1:3
                    if ~isempty(trials_StimSort_av{iav}{i})
                        ind = (n(i)+1):(n(i+1));
                        respStimSort_av{iav}{i} = respAllCells(ind,:);
                    end
                end
            end
          
            
            for iav = 1:2
                trOut = decodeDataExpt(iexp).av(iav).outcome;
                respStimSort = respStimSort_av{iav};
                trOutStimSort = trOutStimSort_av{iav};
                stimSortInd = stimSortInd_av{iav};
                matchTrialsInd = matchTrialsInd_av{iav};

                if iav == 1
                    resp = respAllCells_vis;
                    resp_pcs = respAllCells_pcs_vis;
                elseif iav == 2
                    resp = respAllCells_aud;
                    resp_pcs = respAllCells_pcs_aud;
                end
                
                emptyInd = cellfun(@isempty,respStimSort);
                respStimSort(~emptyInd) = cellfun(@(x) x(:,1:nPCs),...
                    respStimSort(~emptyInd),'unif',0);
                
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
                 
                % test model with ea stim
                pctCorrTarget_xStim_train = nan(1,nStimBins);
                pctCorrTarget_xStim_ho = nan(1,nStimBins);
                for istim = 1:nStimBins
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
                
                % response window analysis
                fprintf('Expt %s, starting resp win analysis\n',num2str(iexp))
                nwins = size(decodeDataExpt(iexp).av(iav).movWinResp,3);
                pctCorrectTarget_train_respwin = nan(1,nwins);
                pctCorrectTarget_ho_respwin = nan(1,nwins);
                for iwin = 1:nwins
%                     rAllCells = zscore(decodeDataExpt(iexp).av(iav).movWinResp(:,:,iwin)');
%                     r = rAllCells(matchTrialsInd,cellInd);
                    rAllCells = decodeDataExpt(iexp).av(iav).movWinResp(cellInd,matchTrialsInd,iwin)';
%                     [~,sAllCells] = pca(rAllCells(:,cellInd));
                    r = zscore(rAllCells);
                    [~,~,targetGLM_temp] = glmfit(r*C,targetTrInd,'binomial');
                    pctCorrectTarget_train_respwin(iwin) = getPctCorr_trainData(...
                        targetGLM_temp,r,targetTrInd,dv_target);
                    pctCorrectTarget_ho_respwin(iwin) = getPctCorr_hoData(r,targetTrInd,dv_target);
                end

                respOther{iav} = resp;
                trOutOther{iav} = trOut(matchTrialsInd);
                targetGLMOther{iav} = targetGLM;  
                    
                % model with pcs                
                C = eye(size(resp_pcs,2));
                p=1;
                [~,~,targetGLM] = glmfit(resp_pcs*C,targetTrInd,'binomial');
                                
                targetWeight_pcs = targetGLM.beta(2:end);
                
                fprintf('Expt %s, hold-out analysis\n',num2str(iexp))
               
                pctCorrectTarget_pcs_ho = getPctCorr_hoData(resp_pcs,targetTrInd,dv_target);
 
                targetWeight_pc2neur = coeffAllCells(:,1:nPCs)*targetGLM.beta(2:end);
                
                respOther_pcs{iav} = resp_pcs;
                targetGLMOther_pcs{iav} = targetGLM; 
                
%                 % train detect model with distractor choices only
                fprintf('Expt %s, distractor only model\n',num2str(iexp))
                decodeAnalysis(iexp).av(iav).dvTarget = dv_target;
                decodeAnalysis(iexp).av(iav).correlationTarget = targetCorr;
                decodeAnalysis(iexp).av(iav).weightTarget = targetWeight;
                decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_train = pctCorrectTarget_train;
                decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_holdout = pctCorrectTarget_ho;
                decodeAnalysis(iexp).av(iav).pctCorrectXStimTarget_train = pctCorrTarget_xStim_train;
                decodeAnalysis(iexp).av(iav).pctCorrectXStimTarget_holdout = pctCorrTarget_xStim_ho;
                decodeAnalysis(iexp).av(iav).pctCorrectTargetMovRespWin_train = pctCorrectTarget_train_respwin;
                decodeAnalysis(iexp).av(iav).pctCorrectTargetMovRespWin_holdout = pctCorrectTarget_ho_respwin;
                
                decodeAnalysis(iexp).av(iav).weightTarget_pcs = targetWeight_pcs;
                decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_pcs_holdout = pctCorrectTarget_pcs_ho;
                
                
            end            
            fprintf('Expt %s: Testing opposite model...\n',num2str(iexp))
            for iav = 1:2
                if iav == 1
                    otherAV = 2;
                elseif iav == 2
                    otherAV = 1;
                end
                resp = respOther{otherAV};
                [~, targetTrInd] = getStimAndBehaviorYs(trOutOther{otherAV});

                [~,targettrind4model] = getStimAndBehaviorYs(trOutOther{iav});
                dv_target = mean(targettrind4model);
                
                pctCorrectTarget = getPctCorr_trainData(targetGLMOther{iav},resp,targetTrInd,dv_target);
                
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_otherAV = pctCorrectTarget;
                
                
            end
            
        end
        
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
    
    save([fnout 'imgAnalysisData'],...
        'decodeAnalysis','antiAnalysis','targetAnalysis','cellInfo','respCellsExpt','antiDataExpt','siPerExpt','oriTuningExpt')
end