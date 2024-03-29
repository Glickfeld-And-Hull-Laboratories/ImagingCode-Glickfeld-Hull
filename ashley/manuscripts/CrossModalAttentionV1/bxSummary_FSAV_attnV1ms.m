% run eaMsBxSummary_attnV1ms and check bxParams_FSAV_attnV1ms before starting
%% load data
clear all
close all
% useRandSeed = true;
rc = behavConstsAV;
exptSummaryDir = fullfile(rc.ashley,...
    'Manuscripts','Attention V1','Mouse Info.xlsx');
exptSummaryInfo = readtable(exptSummaryDir);
fnout = fullfile(rc.ashley,'Manuscripts','Attention V1','Matlab Figs');
% 
% ms2analyze = cellfun(@num2str,num2cell(exptSummaryInfo.SubjectNumber),...
%     'unif',0);
ms2analyze = exptSummaryInfo.SubjectNumber';
nMice = length(ms2analyze);

exampleMouse = '668';
exMsInd = find(strcmp(ms2analyze,exampleMouse));

bxParams_FSAV_attnV1ms

doPlot = true;

% if useRandSeed
%     load(fullfile(fnout,'bxStats'))
% end
    rng(0);

%% compile data

trainRew = strcmp(exptSummaryInfo.TrainingType,'Reward');

msSumStruct = [];
msExptStruct = cell(1,nMice);
nTrialsPerExpt = cell(1,nMice);
nInvPerExpt = cell(1,nMice);
matchedHRall = cell(1,nMice);
matchedHR = cell(2,nMice);
for im = 1:nMice
    mouseName = ms2analyze{im};
    fn = fullfile(rc.ashleyAnalysis,mouseName,'behavior');
    load(fullfile(fn,[mouseName,'bxSummary_dataAnalyzed_attnV1ms']))
    msSumStruct = cat(1,msSumStruct,msCmlvData);
%     if im == exMsInd
%         exMsExptInfo = msExptAnalyzed;
%     end
    load(fullfile(fn,[mouseName,'bxSummary_data.mat']))
    msExpt = msExptInfo(exptInd);
    nTrialsPerExpt{im} = cellfun(@length,{msExpt.hit});
    nInvPerExpt{im} = cellfun(@(x,y,z,a)...
        sum((x > 0 | y > 0) & (z > 0 | a > 0)),...
        {msExpt.tInvVisTargets},{msExpt.tInvAudTargets},...
        {msExpt.invHit},{msExpt.invMiss});
    
    nExpt = sum(exptInd);
    matchedHRall{im} = nan(2,nExpt);
    matchedHR{visualTrials,im} = nan(2,nExpt);
    matchedHR{auditoryTrials,im} = nan(2,nExpt);
    for iexp = 1:nExpt
        [~,visHRall,visInvHRall] = allTrialsAttnTest(msExpt(iexp).tVisTargets,...
            msExpt(iexp).tInvVisTargets,msExpt(iexp).hit,msExpt(iexp).miss,...
            msExpt(iexp).invHit,msExpt(iexp).invMiss);
        [~,audHRall,audInvHRall] = allTrialsAttnTest(msExpt(iexp).tAudTargets,...
            msExpt(iexp).tInvAudTargets,msExpt(iexp).hit,msExpt(iexp).miss,...
            msExpt(iexp).invHit,msExpt(iexp).invMiss);
        tInvVis = msExpt(iexp).tInvVisTargets;
        tInvVis(isnan(tInvVis)) = 0;
        tInvAud = msExpt(iexp).tInvAudTargets;
        tInvAud(isnan(tInvAud)) = 0;
        [~,allHRall,allInvHRall] = allTrialsAttnTest(...
            msExpt(iexp).tVisTargets+msExpt(iexp).tAudTargets,...
            double(tInvVis)+double(tInvAud),...
            msExpt(iexp).hit,msExpt(iexp).miss,...
            msExpt(iexp).invHit,msExpt(iexp).invMiss);
        matchedHRall{im}(valid,iexp) = allHRall;
        matchedHRall{im}(invalid,iexp) = allInvHRall;
        matchedHR{visualTrials,im}(valid,iexp) = visHRall;
        matchedHR{visualTrials,im}(invalid,iexp) = visInvHRall;
        matchedHR{auditoryTrials,im}(valid,iexp) = audHRall;
        matchedHR{auditoryTrials,im}(invalid,iexp) = audInvHRall;
    end
end

msHR = struct;
allData = struct;
rewSortMsData = struct;
rewSortStart = zeros(1,2);
HR = cell(2,2);
HR(:) = {nan(nMice,nBins)};
targets = cell(2,2);
targets(:) = {nan(nMice,nBins)};
for im = 1:nMice
    msCmlvData = msSumStruct(im);
    if im == 1
        allData.cue(valid).av(visualTrials).targets = msCmlvData.tVisTargets;
        allData.cue(valid).hit = msCmlvData.hit;
        allData.cue(valid).miss = msCmlvData.miss;
        allData.cue(valid).av(visualTrials).nFAandDist = msCmlvData.visNFAandDistractors;
        allData.cue(invalid).av(visualTrials).targets = msCmlvData.tInvVisTargets;
        allData.cue(invalid).hit = msCmlvData.invHit;
        allData.cue(invalid).miss = msCmlvData.invMiss;
        allData.cue(valid).av(auditoryTrials).targets = msCmlvData.tAudTargets;
        allData.cue(valid).av(auditoryTrials).nFAandDist = msCmlvData.audNFAandDistractors;
        allData.cue(invalid).av(auditoryTrials).targets = msCmlvData.tInvAudTargets;
    else
        allData.cue(valid).av(visualTrials).targets = cat(2,...
           allData.cue(valid).av(visualTrials).targets, msCmlvData.tVisTargets);
        allData.cue(valid).hit = cat(2,...
           allData.cue(valid).hit,msCmlvData.hit);
        allData.cue(valid).miss =cat(2,...
            allData.cue(valid).miss,msCmlvData.miss);
        allData.cue(valid).av(visualTrials).nFAandDist = ...
                    allData.cue(valid).av(visualTrials).nFAandDist +...
                    msCmlvData.visNFAandDistractors;

        allData.cue(invalid).av(visualTrials).targets = cat(2,...
           allData.cue(invalid).av(visualTrials).targets,msCmlvData.tInvVisTargets);
        allData.cue(invalid).hit = cat(2,...
           allData.cue(invalid).hit,msCmlvData.invHit);
        allData.cue(invalid).miss = cat(2,...
           allData.cue(invalid).miss,msCmlvData.invMiss);
        allData.cue(valid).av(auditoryTrials).targets = cat(2,...
           allData.cue(valid).av(auditoryTrials).targets,msCmlvData.tAudTargets);
        allData.cue(valid).av(auditoryTrials).nFAandDist = ...
                    allData.cue(valid).av(auditoryTrials).nFAandDist +...
                    msCmlvData.audNFAandDistractors;
        allData.cue(invalid).av(auditoryTrials).targets = cat(2,...
           allData.cue(invalid).av(auditoryTrials).targets,msCmlvData.tInvAudTargets);
    end
    
%     visTargets = unique(msCmlvData.tVisTargets);
%     visTargets = visTargets(2:end);
%     audTargets = unique(msCmlvData.tAudTargets);
%     audTargets = audTargets(2:end);

%     visBinEdges = exp(linspace(log(min(visTargets)-1),log(max(visTargets)),nBins+1));
%     audBinEdges = exp(linspace(...
%         log(min(audTargets(audTargets > 0.00001))-...
%         (0.5*min(audTargets(audTargets > 0.00001)))),...
%         log(max(audTargets)),nBins+1));
    [~,~,visBinInd] = histcounts(msCmlvData.tVisTargets,visBinEdges);
    [~,~,audBinInd] = histcounts(msCmlvData.tAudTargets,audBinEdges);
    [~,~,invVisBinInd] = histcounts(msCmlvData.tInvVisTargets,visBinEdges);
    [~,~,invAudBinInd] = histcounts(msCmlvData.tInvAudTargets,audBinEdges);
        
    [visRTBinned,~,visRTtargets,~,visRTanovaP] = getBinnedRT(msCmlvData.tVisTargets,visBinEdges,...
        nBins,minTrN_ms,msCmlvData.valRT,msCmlvData.hit);
    [audRTBinned,~,audRTtargets,~,audRTanovaP] = getBinnedRT(msCmlvData.tAudTargets,audBinEdges,...
        nBins,minTrN_ms,msCmlvData.valRT,msCmlvData.hit);
    close all
    
    nVisHits = nan(1,nBins);
    nVisMisses = nan(1,nBins);
    nAudHits = nan(1,nBins);
    nAudMisses = nan(1,nBins);
    visTargetsBinned = nan(1,nBins);
    audTargetsBinned = nan(1,nBins);
    visTargetsSte = nan(1,nBins);
    audTargetsSte = nan(1,nBins);
    nInvVisHits = nan(1,nBins);
    nInvVisMisses = nan(1,nBins);
    nInvAudHits = nan(1,nBins);
    nInvAudMisses = nan(1,nBins);
    invVisTargetsBinned = nan(1,nBins);
    invAudTargetsBinned = nan(1,nBins);
    invVisTargetsSte = nan(1,nBins);
    invAudTargetsSte = nan(1,nBins);   
%     visRTBinned = nan(1,nBins);
%     audRTBinned = nan(1,nBins);
%     visRTBinnedSte = nan(1,nBins);
%     audRTBinnedSte = nan(1,nBins);
   
    for ibin = 1:nBins
        ind = (msCmlvData.hit | msCmlvData.miss) & visBinInd == ibin;
        if sum(ind) > minTrN_ms
            nVisHits(ibin) = sum(msCmlvData.hit & visBinInd == ibin);
            nVisMisses(ibin) = sum(msCmlvData.miss & visBinInd == ibin);
            visTargetsBinned(ibin) = mean(msCmlvData.tVisTargets(ind));
            visTargetsSte(ibin) = ste(msCmlvData.tVisTargets(ind),2);
%             visRTBinned(ibin) = mean(msCmlvData.valRT(ind & msCmlvData.hit));
%             visRTBinnedSte(ibin) = ste(msCmlvData.valRT(ind & msCmlvData.hit),2);
        end

        ind = (msCmlvData.hit | msCmlvData.miss) & audBinInd == ibin;
        if sum(ind) > minTrN_ms
            nAudHits(ibin) = sum(msCmlvData.hit & audBinInd == ibin);
            nAudMisses(ibin) = sum(msCmlvData.miss & audBinInd == ibin);
            audTargetsBinned(ibin) = mean(msCmlvData.tAudTargets(ind));
            audTargetsSte(ibin) = ste(msCmlvData.tAudTargets(ind),2);
%             audRTBinned(ibin) = mean(msCmlvData.valRT(ind & msCmlvData.hit));
%             audRTBinnedSte(ibin) = ste(msCmlvData.valRT(ind & msCmlvData.hit),2);
        end

        ind = (msCmlvData.invHit | msCmlvData.invMiss) & invVisBinInd == ibin;
        if sum(ind) > minTrN_ms
            nInvVisHits(ibin) = sum(invVisBinInd == ibin & msCmlvData.invHit);
            nInvVisMisses(ibin) = sum(invVisBinInd == ibin & msCmlvData.invMiss);
            invVisTargetsBinned(ibin) = mean(msCmlvData.tInvVisTargets(ind));
            invVisTargetsSte(ibin) = ste(msCmlvData.tInvVisTargets(ind),2);
        end

        ind = (msCmlvData.invHit | msCmlvData.invMiss) & invAudBinInd == ibin;
        if sum(ind) > minTrN_ms  
            nInvAudHits(ibin) = sum(invAudBinInd == ibin & msCmlvData.invHit);
            nInvAudMisses(ibin) = sum(invAudBinInd == ibin & msCmlvData.invMiss);
            invAudTargetsBinned(ibin) = mean(msCmlvData.tInvAudTargets(ind));
            invAudTargetsSte(ibin) = ste(msCmlvData.tInvAudTargets(ind),2);
        end
    end
    visInd = ~isnan(nVisHits);
    invVisInd = ~isnan(nInvVisHits); 
    audInd = ~isnan(nAudHits);
    invAudInd = ~isnan(nInvAudHits); 
    
    [visAttnTest,visHRall,visInvHRall] = allTrialsAttnTest(msCmlvData.tVisTargets,...
        msCmlvData.tInvVisTargets,msCmlvData.hit,msCmlvData.miss,...
        msCmlvData.invHit,msCmlvData.invMiss);
    [audAttnTest,audHRall,audInvHRall] = allTrialsAttnTest(msCmlvData.tAudTargets,...
        msCmlvData.tInvAudTargets,msCmlvData.hit,msCmlvData.miss,...
        msCmlvData.invHit,msCmlvData.invMiss);
    
    [allAttnTest,allHRall,allInvHRall] = allTrialsAttnTest(...
        msCmlvData.tVisTargets+msCmlvData.tAudTargets,...
        msCmlvData.tInvVisTargets+msCmlvData.tInvAudTargets,...
        msCmlvData.hit,msCmlvData.miss,...
        msCmlvData.invHit,msCmlvData.invMiss);
    if allAttnTest < attnTestAlpha
        allAttnPower = sampsizepwr('p',allHRall,allInvHRall,sampleSizePower,[],'tail','left');
    else
        allAttnPower = nan;
    end
    allAttnN = sum((msCmlvData.tInvVisTargets+msCmlvData.tInvAudTargets) > 0);
        
    [visHR,visHR95ci] = binofit(nVisHits(visInd),nVisHits(visInd)+nVisMisses(visInd));
    [audHR,audHR95ci] = binofit(nAudHits(audInd),nAudHits(audInd)+nAudMisses(audInd));
    if sum(nInvVisHits(invVisInd)+nInvVisMisses(invVisInd)) > 0
        [invVisHR,invVisHR95ci] = binofit(nInvVisHits(invVisInd),...
            nInvVisHits(invVisInd)+nInvVisMisses(invVisInd));
    else
        invVisHR = nan;
        invVisHR95ci = nan(1,2);
    end
    if sum(nInvAudHits(invAudInd)+nInvAudMisses(invAudInd)) > 0
        [invAudHR,invAudHR95ci] = binofit(nInvAudHits(invAudInd),...
            nInvAudHits(invAudInd)+nInvAudMisses(invAudInd));
    else
        invAudHR = nan;
        invAudHR95ci = nan(1,2);
    end
    
    nTrials = nVisHits(visInd)+nVisMisses(visInd);
    msVisFit = weibullFitLG(visTargetsBinned(visInd), visHR, 0,0, {'nTrials',nTrials});

    maxI = max(visTargetsBinned(visInd));
    minI = min(visTargetsBinned(visInd));
    msVisXGrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);

    nTrials = nAudHits(audInd)+nAudMisses(audInd);
    try
        msAudFit = weibullFitLG(audTargetsBinned(audInd), audHR, 0,0, {'nTrials',nTrials});
    catch
        msAudFit = nan;
    end
    maxI = max(audTargetsBinned(audInd));
    minI = min(audTargetsBinned(audInd));
    msAudXGrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);
    
    FAR_vis = msCmlvData.visNFAandDistractors(1)./msCmlvData.visNFAandDistractors(2);
    FAR_aud = msCmlvData.audNFAandDistractors(1)./msCmlvData.audNFAandDistractors(2);
    
    visHRWithFA = cat(2,FAR_vis,visHR);
    visTargetsWithFA = cat(2,0,visTargetsBinned(visInd));
    nVisTrialsWithFA = cat(2,msCmlvData.visNFAandDistractors(2),...
        nVisHits(visInd)+nVisMisses(visInd));
    msVisFitWithFA = weibullFitLG(visTargetsWithFA, visHRWithFA, ...
        0,0, {'nTrials',nVisTrialsWithFA});

    audHRWithFA = cat(2,FAR_aud,audHR);
    audTargetsWithFA = cat(2,0,audTargetsBinned(audInd));
    nAudTrialsWithFA = cat(2,msCmlvData.audNFAandDistractors(2),...
        nAudHits(audInd)+nAudMisses(audInd));
    msAudFitWithFA = weibullFitLG(audTargetsWithFA, audHRWithFA, ...
        0,0, {'nTrials',nAudTrialsWithFA});
    if allAttnTest < attnTestAlpha
        if trainRew(im)
            rewSortInd = 1;
        else
            rewSortInd = 2;     
        end
        if rewSortStart(rewSortInd) == 1
            rewSortMsData(rewSortInd).av(visualTrials).cue(valid).targets = cat(2,msCmlvData.tVisTargets,...
                rewSortMsData(rewSortInd).av(visualTrials).cue(valid).targets);
            rewSortMsData(rewSortInd).hit = cat(2,msCmlvData.hit,...
                rewSortMsData(rewSortInd).hit);
            rewSortMsData(rewSortInd).miss = cat(2,msCmlvData.miss,...
                rewSortMsData(rewSortInd).miss);

            rewSortMsData(rewSortInd).av(visualTrials).cue(invalid).targets = cat(2,msCmlvData.tInvVisTargets,...
                rewSortMsData(rewSortInd).av(visualTrials).cue(invalid).targets);
            rewSortMsData(rewSortInd).invHit = cat(2,msCmlvData.invHit,...
                rewSortMsData(rewSortInd).invHit);
            rewSortMsData(rewSortInd).invMiss = cat(2,msCmlvData.invMiss,...
                rewSortMsData(rewSortInd).invMiss);

            rewSortMsData(rewSortInd).av(auditoryTrials).cue(valid).targets = cat(2,msCmlvData.tAudTargets,...
                rewSortMsData(rewSortInd).av(auditoryTrials).cue(valid).targets);

            rewSortMsData(rewSortInd).av(auditoryTrials).cue(invalid).targets = cat(2,msCmlvData.tInvAudTargets,...
                rewSortMsData(rewSortInd).av(auditoryTrials).cue(invalid).targets);
        else
            rewSortMsData(rewSortInd).av(visualTrials).cue(valid).targets = msCmlvData.tVisTargets;
            rewSortMsData(rewSortInd).hit = msCmlvData.hit;
            rewSortMsData(rewSortInd).miss = msCmlvData.miss;

            rewSortMsData(rewSortInd).av(visualTrials).cue(invalid).targets = msCmlvData.tInvVisTargets;
            rewSortMsData(rewSortInd).invHit = msCmlvData.invHit;
            rewSortMsData(rewSortInd).invMiss = msCmlvData.invMiss;

            rewSortMsData(rewSortInd).av(auditoryTrials).cue(valid).targets = msCmlvData.tAudTargets;

            rewSortMsData(rewSortInd).av(auditoryTrials).cue(invalid).targets = msCmlvData.tInvAudTargets;
        end
        if trainRew(im)
            rewSortStart(1) = 1;
        else   
            rewSortStart(2) = 1;
        end
    end
    
    f = msVisFitWithFA.modelFun(msVisFitWithFA.coefEsts, msVisXGrid);
    oriAtHighThresh = msVisXGrid(find(f > highThreshold,1));
    if isempty(oriAtHighThresh)
        oriAtHighThresh = msVisXGrid(find(f > 0.85,1));
    end
    f = msAudFitWithFA.modelFun(msAudFitWithFA.coefEsts, msAudXGrid);
    ampAtHighThresh = msAudXGrid(find(f > highThreshold,1));
    lowVisAttnTest = belowThreshAttnTest(oriAtHighThresh,msCmlvData.tVisTargets,...
        msCmlvData.tInvVisTargets,msCmlvData.hit,msCmlvData.miss,...
        msCmlvData.invHit,msCmlvData.invMiss);
    lowAudAttnTest = belowThreshAttnTest(ampAtHighThresh,msCmlvData.tAudTargets,...
        msCmlvData.tInvAudTargets,msCmlvData.hit,msCmlvData.miss,...
        msCmlvData.invHit,msCmlvData.invMiss);
    highVisAttnTest = aboveThreshAttnTest(oriAtHighThresh,msCmlvData.tVisTargets,...
        msCmlvData.tInvVisTargets,msCmlvData.hit,msCmlvData.miss,...
        msCmlvData.invHit,msCmlvData.invMiss);
    highAudAttnTest = aboveThreshAttnTest(ampAtHighThresh,msCmlvData.tAudTargets,...
        msCmlvData.tInvAudTargets,msCmlvData.hit,msCmlvData.miss,...
        msCmlvData.invHit,msCmlvData.invMiss);
    
    ind1 = msCmlvData.tVisTargets > 0 & msCmlvData.tVisTargets < oriAtHighThresh;
    ind2 = msCmlvData.tAudTargets > 0 & msCmlvData.tAudTargets < ampAtHighThresh;
    ind3 = msCmlvData.tInvVisTargets > 0 & msCmlvData.tInvVisTargets < oriAtHighThresh;
    ind4 = msCmlvData.tInvAudTargets > 0 & msCmlvData.tInvAudTargets < ampAtHighThresh;
    [loThreshAllAttnTest,loThreshHRall,loThreshInvHRall] = allTrialsAttnTest(...
        cat(2,msCmlvData.tVisTargets(ind1),msCmlvData.tAudTargets(ind2)),...
        cat(2,msCmlvData.tInvVisTargets(ind3),msCmlvData.tInvAudTargets(ind4)),...
        cat(2,msCmlvData.hit(ind1),msCmlvData.hit(ind2)),...
        cat(2,msCmlvData.miss(ind1),msCmlvData.miss(ind2)),...
        cat(2,msCmlvData.invHit(ind3),msCmlvData.invHit(ind4)),...
        cat(2,msCmlvData.invMiss(ind3),msCmlvData.invMiss(ind4)));
    ind1 = msCmlvData.tVisTargets > oriAtHighThresh;
    ind2 = msCmlvData.tAudTargets > ampAtHighThresh;
    ind3 = msCmlvData.tInvVisTargets > oriAtHighThresh;
    ind4 = msCmlvData.tInvAudTargets > ampAtHighThresh;
    [hiThreshAllAttnTest,hiThreshHRall,hiThreshInvHRall] = allTrialsAttnTest(...
        cat(2,msCmlvData.tVisTargets(ind1),msCmlvData.tAudTargets(ind2)),...
        cat(2,msCmlvData.tInvVisTargets(ind3),msCmlvData.tInvAudTargets(ind4)),...
        cat(2,msCmlvData.hit(ind1),msCmlvData.hit(ind2)),...
        cat(2,msCmlvData.miss(ind1),msCmlvData.miss(ind2)),...
        cat(2,msCmlvData.invHit(ind3),msCmlvData.invHit(ind4)),...
        cat(2,msCmlvData.invMiss(ind3),msCmlvData.invMiss(ind4)));
    
    
    
    valRT_vis_loHi = nan(1,2);
    valRT_aud_loHi = nan(1,2);
    for ibin = 1:2
        if ibin == 1
            visRTInd = msCmlvData.tVisTargets > 0 &...
                msCmlvData.tVisTargets < oriAtHighThresh;
            audRTInd = msCmlvData.tAudTargets > 0 &...
                msCmlvData.tAudTargets < ampAtHighThresh;
        else
            visRTInd = msCmlvData.tVisTargets > oriAtHighThresh;
            audRTInd = msCmlvData.tAudTargets > ampAtHighThresh;
        end
        
        valRT_vis_loHi(ibin) = mean(msCmlvData.valRT(visRTInd & msCmlvData.hit));
        valRT_aud_loHi(ibin) = mean(msCmlvData.valRT(audRTInd & msCmlvData.hit));
    end
    
    invInd = msCmlvData.tInvVisTargets > 0 & msCmlvData.invHit;
    matches = cell2mat(getMatchedValidTrialIndex(msCmlvData.tVisTargets,...
        msCmlvData.tInvVisTargets(invInd)));
    valMatchHits = ismember(matches,find(msCmlvData.hit));
    valInd = matches(valMatchHits);
    invRT_vis = mean(msCmlvData.invRT(invInd));
    valRT_vis = mean(msCmlvData.valRT(valInd));

    invInd = msCmlvData.tInvAudTargets > 0 & msCmlvData.invHit;
    matches = cell2mat(getMatchedValidTrialIndex(msCmlvData.tAudTargets,...
        msCmlvData.tInvAudTargets(invInd)));
    valMatchHits = ismember(matches,find(msCmlvData.hit));
    valInd = matches(valMatchHits);
    invRT_aud = mean(msCmlvData.invRT(invInd));
    valRT_aud = mean(msCmlvData.valRT(valInd));
      
    [valHR_highThreshold_vis, invHR_highThreshold_vis] = ...
        getMatchedHighThresholdHR(visTargetsBinned,msVisFitWithFA,highThreshold,...
        msCmlvData.tVisTargets,msCmlvData.tInvVisTargets,...
        msCmlvData.hit,msCmlvData.miss,msCmlvData.invHit,msCmlvData.invMiss);
    
    [valHR_highThreshold_aud, invHR_highThreshold_aud] = ...
        getMatchedHighThresholdHR(audTargetsBinned,msAudFitWithFA,highThreshold,...
        msCmlvData.tAudTargets,msCmlvData.tInvAudTargets,...
        msCmlvData.hit,msCmlvData.miss,msCmlvData.invHit,msCmlvData.invMiss);
    
    msHR(im).matchedHRall = [allHRall,allInvHRall];
    msHR(im).attnTestAll = allAttnTest;
    msHR(im).attnTestPowerTest = allAttnN >= allAttnPower;
    msHR(im).matchedHRhiLo = [loThreshHRall,loThreshInvHRall;...
        hiThreshHRall,hiThreshInvHRall];
    msHR(im).attnTestHiLo = [loThreshAllAttnTest,hiThreshAllAttnTest];
    
    msHR(im).av(visualTrials).cue(valid).HR = visHR.*100;
    msHR(im).av(visualTrials).cue(valid).HR95CI = visHR95ci.*100;
    msHR(im).av(visualTrials).cue(valid).targets = visTargetsBinned(visInd);
    msHR(im).av(visualTrials).cue(valid).targetsErr = visTargetsSte(visInd);
    msHR(im).av(visualTrials).cue(valid).fit = msVisFit;
    msHR(im).av(visualTrials).cue(valid).fitGrid = msVisXGrid;
    msHR(im).av(visualTrials).cue(valid).hiLoHR = valHR_highThreshold_vis;
    msHR(im).av(visualTrials).cue(valid).RT = valRT_vis;
    msHR(im).av(visualTrials).cue(valid).RTloHi = valRT_vis_loHi;
    msHR(im).av(visualTrials).cue(valid).RTbinned = visRTBinned;
    msHR(im).av(visualTrials).cue(valid).RTanovaTest = visRTanovaP;
    msHR(im).av(visualTrials).cue(valid).RTtargets = visRTtargets;
    msHR(im).av(visualTrials).attnTest = visAttnTest;   
    msHR(im).av(visualTrials).belowThreshAttnTest = lowVisAttnTest; 
    msHR(im).av(visualTrials).aboveThreshAttnTest = highVisAttnTest;   
    msHR(im).av(visualTrials).matchedHRall = [visHRall,visInvHRall];
    msHR(im).av(visualTrials).FAR = FAR_vis;
    msHR(im).av(visualTrials).cue(valid).fitWithFAR = msVisFitWithFA;
%     msHR(im).av(visualTrials).threshAttnTest = lowVisAttnTest;    
    
    HR{visualTrials,valid}(im,visInd) = ...
        msHR(im).av(visualTrials).cue(valid).HR;
    targets{visualTrials,valid}(im,visInd) = ...
        msHR(im).av(visualTrials).cue(valid).targets;
    
    msHR(im).av(auditoryTrials).cue(valid).HR = audHR.*100;
    msHR(im).av(auditoryTrials).cue(valid).HR95CI = audHR95ci.*100;
    msHR(im).av(auditoryTrials).cue(valid).targets = audTargetsBinned(audInd);
    msHR(im).av(auditoryTrials).cue(valid).targetsErr = audTargetsSte(audInd);
    msHR(im).av(auditoryTrials).cue(valid).fit = msAudFit;
    msHR(im).av(auditoryTrials).cue(valid).fitGrid = msAudXGrid;
    msHR(im).av(auditoryTrials).cue(valid).hiLoHR = valHR_highThreshold_aud;
    msHR(im).av(auditoryTrials).cue(valid).RT = valRT_aud;
    msHR(im).av(auditoryTrials).cue(valid).RTloHi = valRT_aud_loHi;
    msHR(im).av(auditoryTrials).cue(valid).RTbinned = audRTBinned;
    msHR(im).av(auditoryTrials).cue(valid).RTanovaTest = audRTanovaP;
    msHR(im).av(auditoryTrials).cue(valid).RTtargets = audRTtargets;
    msHR(im).av(auditoryTrials).attnTest = audAttnTest;
    msHR(im).av(auditoryTrials).belowThreshAttnTest = lowAudAttnTest; 
    msHR(im).av(auditoryTrials).aboveThreshAttnTest = highAudAttnTest; 
    msHR(im).av(auditoryTrials).matchedHRall = [audHRall,audInvHRall];
    msHR(im).av(auditoryTrials).FAR = FAR_aud;
    msHR(im).av(auditoryTrials).cue(valid).fitWithFAR = msAudFitWithFA;
    
    
    HR{auditoryTrials,valid}(im,audInd) = ...
        msHR(im).av(auditoryTrials).cue(valid).HR;
    targets{auditoryTrials,valid}(im,audInd) = ...
        msHR(im).av(auditoryTrials).cue(valid).targets;
    
    msHR(im).av(visualTrials).cue(invalid).HR = invVisHR.*100;
    msHR(im).av(visualTrials).cue(invalid).HR95CI = invVisHR95ci.*100;
    msHR(im).av(visualTrials).cue(invalid).targets = invVisTargetsBinned(invVisInd);
        msHR(im).av(visualTrials).cue(invalid).targetsErr = invVisTargetsSte(invVisInd);
    msHR(im).av(visualTrials).cue(invalid).hiLoHR = invHR_highThreshold_vis;
    msHR(im).av(visualTrials).cue(invalid).RT = invRT_vis;
    
    HR{visualTrials,invalid}(im,invVisInd) = ...
        msHR(im).av(visualTrials).cue(invalid).HR;
    targets{visualTrials,invalid}(im,invVisInd) = ...
        msHR(im).av(visualTrials).cue(invalid).targets;
    
    msHR(im).av(auditoryTrials).cue(invalid).HR = invAudHR.*100;
    msHR(im).av(auditoryTrials).cue(invalid).HR95CI = invAudHR95ci.*100;
    msHR(im).av(auditoryTrials).cue(invalid).targets = invAudTargetsBinned(invAudInd);
    msHR(im).av(auditoryTrials).cue(invalid).targetsErr = invAudTargetsSte(invAudInd);
    msHR(im).av(auditoryTrials).cue(invalid).hiLoHR = invHR_highThreshold_aud;
    msHR(im).av(auditoryTrials).cue(invalid).RT = invRT_aud;
    
    HR{auditoryTrials,invalid}(im,invAudInd) = ...
        msHR(im).av(auditoryTrials).cue(invalid).HR;
    targets{auditoryTrials,invalid}(im,invAudInd) = ...
        msHR(im).av(auditoryTrials).cue(invalid).targets;    
    
    msHR(im).name = ms2analyze{im};
    msHR(im).nTrialsRange = [];
end

visTargets = unique(allData.cue(valid).av(visualTrials).targets);
visTargets = visTargets(2:end);
audTargets = unique(allData.cue(valid).av(auditoryTrials).targets);
audTargets = audTargets(2:end);
visBinEdges = exp(linspace(log(min(visTargets)-1),log(max(visTargets)),nBins+1));
audBinEdges = exp(linspace(...
    log(min(audTargets(audTargets > 0.00001))-...
    (0.5*min(audTargets(audTargets > 0.00001)))),...
    log(max(audTargets)),nBins+1));
binEdges = {visBinEdges,audBinEdges};
allHRstruct = struct;
for icue = 1:2
    for iav = 1:2
        if icue == valid
            [allHRstruct.cue(icue).av(iav).HR,allHRstruct.cue(icue).av(iav).HRci,...
                allHRstruct.cue(icue).av(iav).targets,allHRstruct.cue(icue).av(iav).targetsErr,...
                allHRstruct.cue(icue).av(iav).fit] = ...
                getBinnedHRandFit(allData.cue(icue).av(iav).targets,...
                binEdges{iav},nBins,50,allData.cue(icue).hit,allData.cue(icue).miss,...
                allData.cue(icue).av(iav).nFAandDist(1)./allData.cue(icue).av(iav).nFAandDist(2),...
                allData.cue(icue).av(iav).nFAandDist(2));
            
            [allHRstruct.cue(icue).av(iav).attnTest,...
                allHRstruct.cue(valid).av(iav).allHRstruct, allHRstruct.cue(invalid).av(iav).allHRstruct] = ...
                allTrialsAttnTest(allData.cue(valid).av(iav).targets,...
                allData.cue(invalid).av(iav).targets,allData.cue(valid).hit,...
                allData.cue(valid).miss,allData.cue(invalid).hit,...
                allData.cue(invalid).miss);
        else
            [allHRstruct.cue(icue).av(iav).HR,allHRstruct.cue(icue).av(iav).HRci,...
                allHRstruct.cue(icue).av(iav).targets,allHRstruct.cue(icue).av(iav).targetsErr] = ...
                getBinnedHRandFit(allData.cue(icue).av(iav).targets,...
                binEdges{iav},nBins,50,allData.cue(icue).hit,allData.cue(icue).miss);
        end
    end
end

rewSortMsHR = struct;
rewSortMsHR(1).name = 'Training Rewarded';
rewSortMsHR(2).name = 'Training Not Rewarded';
for im = 1:2
    hit = rewSortMsData(im).hit;
    miss = rewSortMsData(im).miss;
    invHit = rewSortMsData(im).invHit;
    invMiss = rewSortMsData(im).invMiss;
    
    tVisTargets = rewSortMsData(im).av(visualTrials).cue(valid).targets;
    tAudTargets = rewSortMsData(im).av(auditoryTrials).cue(valid).targets;
    visTargets = unique(tVisTargets);
    visTargets = visTargets(2:end);
    audTargets = unique(tAudTargets);
    audTargets = audTargets(2:end);
    
    tInvVisTargets = rewSortMsData(im).av(visualTrials).cue(invalid).targets;
    tInvAudTargets = rewSortMsData(im).av(auditoryTrials).cue(invalid).targets;
    
%     visBinEdges = exp(linspace(log(min(visTargets)-1),log(max(visTargets)),nBins+1));
%     audBinEdges = exp(linspace(...
%         log(min(audTargets(audTargets > 0.00001))-...
%         (0.5*min(audTargets(audTargets > 0.00001)))),...
%         log(max(audTargets)),nBins+1));
    [~,~,visBinInd] = histcounts(tVisTargets,visBinEdges);
    [~,~,audBinInd] = histcounts(tAudTargets,audBinEdges);
    [~,~,invVisBinInd] = histcounts(tInvVisTargets,visBinEdges);
    [~,~,invAudBinInd] = histcounts(tInvAudTargets,audBinEdges);
        
    nVisHits = nan(1,nBins);
    nVisMisses = nan(1,nBins);
    nAudHits = nan(1,nBins);
    nAudMisses = nan(1,nBins);
    visTargetsBinned = nan(1,nBins);
    audTargetsBinned = nan(1,nBins);
    visTargetsSte = nan(1,nBins);
    audTargetsSte = nan(1,nBins);
    nInvVisHits = nan(1,nBins);
    nInvVisMisses = nan(1,nBins);
    nInvAudHits = nan(1,nBins);
    nInvAudMisses = nan(1,nBins);
    invVisTargetsBinned = nan(1,nBins);
    invAudTargetsBinned = nan(1,nBins);
    invVisTargetsSte = nan(1,nBins);
    invAudTargetsSte = nan(1,nBins); 
%     visRTBinned = nan(1,nBins);
%     audRTBinned = nan(1,nBins);
%     visRTBinnedSte = nan(1,nBins);
%     audRTBinnedSte = nan(1,nBins);
    for ibin = 1:nBins
        ind = (hit | miss) & visBinInd == ibin;
        if sum(ind) > minTrN_ms
            nVisHits(ibin) = sum(hit & visBinInd == ibin);
            nVisMisses(ibin) = sum(miss & visBinInd == ibin);
            visTargetsBinned(ibin) = mean(tVisTargets(ind));
            visTargetsSte(ibin) = ste(tVisTargets(ind),2);
%             visRTBinned(ibin) = mean(msCmlvData.valRT(ind & hit));
%             visRTBinnedSte(ibin) = ste(msCmlvData.valRT(ind & hit),2);
        end

        ind = (hit | miss) & audBinInd == ibin;
        if sum(ind) > minTrN_ms
            nAudHits(ibin) = sum(hit & audBinInd == ibin);
            nAudMisses(ibin) = sum(miss & audBinInd == ibin);
            audTargetsBinned(ibin) = mean(tAudTargets(ind));
            audTargetsSte(ibin) = ste(tAudTargets(ind),2);
%             audRTBinned(ibin) = mean(msCmlvData.valRT(ind & hit));
%             audRTBinned(ibin) = ste(msCmlvData.valRT(ind & hit),2);
        end

        ind = (invHit | invMiss) & invVisBinInd == ibin;
        if sum(ind) > minTrN_ms
            nInvVisHits(ibin) = sum(invVisBinInd == ibin & invHit);
            nInvVisMisses(ibin) = sum(invVisBinInd == ibin & invMiss);
            invVisTargetsBinned(ibin) = mean(tInvVisTargets(ind));
            invVisTargetsSte(ibin) = ste(tInvVisTargets(ind),2);
        end

        ind = (invHit | invMiss) & invAudBinInd == ibin;
        if sum(ind) > minTrN_ms  
            nInvAudHits(ibin) = sum(invAudBinInd == ibin & invHit);
            nInvAudMisses(ibin) = sum(invAudBinInd == ibin & invMiss);
            invAudTargetsBinned(ibin) = mean(tInvAudTargets(ind));
            invAudTargetsSte(ibin) = ste(tInvAudTargets(ind),2);
        end
    end
    visInd = ~isnan(nVisHits);
    invVisInd = ~isnan(nInvVisHits); 
    audInd = ~isnan(nAudHits);
    invAudInd = ~isnan(nInvAudHits); 
    
    [visHR,visHR95ci] = binofit(nVisHits(visInd),nVisHits(visInd)+nVisMisses(visInd));
    [audHR,audHR95ci] = binofit(nAudHits(audInd),nAudHits(audInd)+nAudMisses(audInd));
    if sum(nInvVisHits(invVisInd)+nInvVisMisses(invVisInd)) > 0
        [invVisHR,invVisHR95ci] = binofit(nInvVisHits(invVisInd),...
            nInvVisHits(invVisInd)+nInvVisMisses(invVisInd));
    else
        invVisHR = nan;
        invVisHR95ci = nan(1,2);
    end
    if sum(nInvAudHits(invAudInd)+nInvAudMisses(invAudInd)) > 0
        [invAudHR,invAudHR95ci] = binofit(nInvAudHits(invAudInd),...
            nInvAudHits(invAudInd)+nInvAudMisses(invAudInd));
    else
        invAudHR = nan;
        invAudHR95ci = nan(1,2);
    end
    FAR_vis = msCmlvData.visNFAandDistractors(1)./msCmlvData.visNFAandDistractors(2);
    FAR_aud = msCmlvData.audNFAandDistractors(1)./msCmlvData.audNFAandDistractors(2);
    
    visHRWithFA = cat(2,FAR_vis,visHR);
    visTargetsWithFA = cat(2,0,visTargetsBinned(visInd));
    nVisTrialsWithFA = cat(2,msCmlvData.visNFAandDistractors(2),...
        nVisHits(visInd)+nVisMisses(visInd));
    msVisFitWithFA = weibullFitLG(visTargetsWithFA, visHRWithFA, ...
        0,0, {'nTrials',nVisTrialsWithFA});
    
    audHRWithFA = cat(2,FAR_aud,audHR);
    audTargetsWithFA = cat(2,0,audTargetsBinned(audInd));
    nAudTrialsWithFA = cat(2,msCmlvData.audNFAandDistractors(2),...
        nAudHits(audInd)+nAudMisses(audInd));
    msAudFitWithFA = weibullFitLG(audTargetsWithFA, audHRWithFA, ...
        0,0, {'nTrials',nAudTrialsWithFA});
    
    [valHR_highThreshold_vis, invHR_highThreshold_vis,...
        valCI_highThreshold_vis, invCI_highThreshold_vis,...
        nVal_highThreshold_vis, nInv_highThreshold_vis] = ...
        getMatchedHighThresholdHR(visTargetsBinned,msVisFitWithFA,highThreshold,...
        tVisTargets,tInvVisTargets,...
        hit,miss,invHit,invMiss);
    
    [valHR_highThreshold_aud, invHR_highThreshold_aud,...
        valCI_highThreshold_aud, invCI_highThreshold_aud,...
        nVal_highThreshold_aud, nInv_highThreshold_aud] = ...
        getMatchedHighThresholdHR(audTargetsBinned,msAudFitWithFA,highThreshold,...
        tAudTargets,tInvAudTargets,...
        hit,miss,invHit,invMiss);
    
    invAllHR_vis = sum(tInvVisTargets > 0 & invHit)./...
        sum(tInvVisTargets > 0 &(invHit | invMiss));
    valMatchInd = ismember(tVisTargets,unique(tInvVisTargets(tInvVisTargets > 0)));
    valAllHR_vis = sum(valMatchInd & hit)./...
        sum(valMatchInd &(hit | miss));
    nMatchedVis = [sum(valMatchInd &(hit | miss)),...
        sum(tInvVisTargets > 0 &(invHit | invMiss))];
    
    invAllHR_aud = sum(tInvAudTargets > 0 & invHit)./...
        sum(tInvAudTargets > 0 &(invHit | invMiss));
    valMatchInd = ismember(tAudTargets,unique(tInvAudTargets(tInvAudTargets > 0)));
    valAllHR_aud = sum(valMatchInd & hit)./...
        sum(valMatchInd &(hit | miss));
    nMatchedAud = [sum(valMatchInd &(hit | miss)),...
        sum(tInvAudTargets > 0 &(invHit | invMiss))];
    
    invAllHR_all = sum((tInvVisTargets > 0|tInvAudTargets > 0) & invHit)./...
        sum((tInvVisTargets > 0|tInvAudTargets > 0) &(invHit | invMiss));
%     valMatchInd = ismember(tVisTargets,unique(tInvVisTargets(tInvVisTargets > 0)));
    invInd = (tInvVisTargets+tInvAudTargets)>0;
    valMatchInd = cell2mat(getMatchedValidTrialIndex(tVisTargets+tAudTargets,...
        tInvVisTargets(invInd)+tInvAudTargets(invInd)));
    valAllHR_all = sum(hit(valMatchInd))./...
        sum(hit(valMatchInd) | miss(valMatchInd));
    nMatchedAll = [sum(hit(valMatchInd) | miss(valMatchInd)),...
        sum((tInvVisTargets+tInvAudTargets) > 0 &(invHit | invMiss))];
    
    
    rewSortMsHR(im).valInvAllTrialsHR = [valAllHR_all invAllHR_all];
    rewSortMsHR(im).matchedTrialN = nMatchedAll;
    rewSortMsHR(im).av(visualTrials).cue(valid).HR = visHR.*100;
    rewSortMsHR(im).av(visualTrials).cue(valid).HR95ci = visHR95ci.*100;
    rewSortMsHR(im).av(visualTrials).cue(valid).hiLoHR = valHR_highThreshold_vis.*100;
    rewSortMsHR(im).av(visualTrials).cue(valid).hiLoHR95ci = valCI_highThreshold_vis.*100;
    rewSortMsHR(im).av(visualTrials).cue(valid).nHiLo = nVal_highThreshold_vis;
    rewSortMsHR(im).av(visualTrials).cue(valid).targets = visTargetsBinned(visInd);
    rewSortMsHR(im).av(visualTrials).cue(valid).targetsErr = visTargetsSte(visInd);
    rewSortMsHR(im).av(visualTrials).fit = msVisFitWithFA;    
    rewSortMsHR(im).av(visualTrials).cue(invalid).HR = invVisHR.*100;
    rewSortMsHR(im).av(visualTrials).cue(invalid).HR95ci = invVisHR95ci.*100;
    rewSortMsHR(im).av(visualTrials).cue(invalid).hiLoHR = invHR_highThreshold_vis.*100;
    rewSortMsHR(im).av(visualTrials).cue(invalid).hiLoHR95ci = invCI_highThreshold_vis.*100;
    rewSortMsHR(im).av(visualTrials).cue(invalid).nHiLo = nInv_highThreshold_vis;
    rewSortMsHR(im).av(visualTrials).cue(invalid).targets = invVisTargetsBinned(~isnan(invVisTargetsBinned));
    rewSortMsHR(im).av(visualTrials).cue(invalid).targetsErr = invVisTargetsSte(~isnan(invVisTargetsBinned));
    rewSortMsHR(im).av(visualTrials).valInvAllTrialsHR = [valAllHR_vis, invAllHR_vis];
    rewSortMsHR(im).av(visualTrials).matchedTrialN = nMatchedVis;
    
    rewSortMsHR(im).av(auditoryTrials).cue(valid).HR = audHR.*100;
    rewSortMsHR(im).av(auditoryTrials).cue(valid).HR95ci = audHR95ci.*100;
    rewSortMsHR(im).av(auditoryTrials).cue(valid).hiLoHR = valHR_highThreshold_aud.*100;
    rewSortMsHR(im).av(auditoryTrials).cue(valid).hiLoHR95ci = valCI_highThreshold_aud.*100;
    rewSortMsHR(im).av(auditoryTrials).cue(valid).nHiLo = nVal_highThreshold_aud;
    rewSortMsHR(im).av(auditoryTrials).cue(valid).targets = audTargetsBinned(audInd);
    rewSortMsHR(im).av(auditoryTrials).cue(valid).targetsErr = audTargetsSte(audInd);
    rewSortMsHR(im).av(auditoryTrials).fit = msAudFitWithFA;    
    rewSortMsHR(im).av(auditoryTrials).cue(invalid).HR = invAudHR.*100;
    rewSortMsHR(im).av(auditoryTrials).cue(invalid).HR95ci = invAudHR95ci.*100;
    rewSortMsHR(im).av(auditoryTrials).cue(invalid).hiLoHR = invHR_highThreshold_aud.*100;
    rewSortMsHR(im).av(auditoryTrials).cue(invalid).hiLoHR95ci = invCI_highThreshold_aud.*100;
    rewSortMsHR(im).av(auditoryTrials).cue(invalid).nHiLo = nInv_highThreshold_aud;
    rewSortMsHR(im).av(auditoryTrials).cue(invalid).targets = invAudTargetsBinned(~isnan(invAudTargetsBinned));
    rewSortMsHR(im).av(auditoryTrials).cue(invalid).targetsErr = invAudTargetsSte(~isnan(invAudTargetsSte));
    rewSortMsHR(im).av(auditoryTrials).valInvAllTrialsHR = [valAllHR_aud, invAllHR_aud];
    rewSortMsHR(im).av(auditoryTrials).matchedTrialN = nMatchedAud;
end


%% Summary Data Structure
visTargets = unique(allData.cue(valid).av(visualTrials).targets);
visTargets = visTargets(2:end);
audTargets = allData.cue(valid).av(auditoryTrials).targets;
audTargets(audTargets < 0.0000001) = 0; % extremely small values are effectively zero, there are only a handful of trials like this
audTargets = unique(audTargets);
audTargets = audTargets(2:end);
[~,sessionAttnTest] = cellfun(@(x) ttest(x(valid,:),x(invalid,:),'tail','right'),...
    matchedHRall,'unif',0);
sessionAttnPower = cellfun(@(x) sampsizepwr('t',...
    [mean(x(valid,:),2), std(x(invalid,:),[],2)],...
    mean(x(invalid,:),2)),matchedHRall);
sessionAttnPowerTest = cellfun(@(x) size(x,2),matchedHRall) >= sessionAttnPower;

for im = 1:nMice
    if im == 1
        visLR = [];
        audLR = [];
        visFAR = [];
        audFAR = [];
        
        visAttnP = [];
        audAttnP = [];
        
        allAttnP = [];
        allAttnPower = [];
        allAttnHiLoP = [];
        allVisAudHR = [];
        allLoHiHRDiff = [];
        
        visHR = [];
        audHR = [];
        
        visRT = [];
        audRT = [];
%         valVisRTLoHiDiff = [];
%         valAudRTLoHiDiff = [];
        visRT_loHi = nan(nMice,2);
        audRT_loHi = nan(nMice,2);        
    end
    visLR = cat(2,visLR,1 - msHR(im).av(visualTrials).cue(valid).HR(end)./100);
    audLR = cat(2,audLR,1 - msHR(im).av(auditoryTrials).cue(valid).HR(end)./100);
    visFAR = cat(2,visFAR,msHR(im).av(visualTrials).FAR);
    audFAR = cat(2,audFAR,msHR(im).av(auditoryTrials).FAR);
    visAttnP = cat(2,visAttnP,msHR(im).av(visualTrials).attnTest);
    audAttnP = cat(2,audAttnP,msHR(im).av(auditoryTrials).attnTest);
    
    allAttnP = cat(2,allAttnP,msHR(im).attnTestAll);
    allAttnPower = cat(2,allAttnPower,msHR(im).attnTestPowerTest);
    allAttnHiLoP = cat(1,allAttnHiLoP,msHR(im).attnTestHiLo);
    allVisAudHR = cat(1,allVisAudHR,msHR(im).matchedHRall);
    allLoHiHRDiff = cat(2,allLoHiHRDiff,...
        msHR(im).matchedHRhiLo(:,valid) - msHR(im).matchedHRhiLo(:,invalid));
    
    visHR = cat(1,visHR,msHR(im).av(visualTrials).matchedHRall);
    audHR = cat(1,audHR,msHR(im).av(auditoryTrials).matchedHRall);
    
    visRT = cat(1,visRT,cat(2,msHR(im).av(visualTrials).cue(valid).RT, ...
            msHR(im).av(visualTrials).cue(invalid).RT));
    audRT = cat(1,audRT,cat(2,msHR(im).av(auditoryTrials).cue(valid).RT, ...
            msHR(im).av(auditoryTrials).cue(invalid).RT));
    for i = 1:2
        visRT_loHi(im,i) = msHR(im).av(visualTrials).cue(valid).RTloHi(i);
        audRT_loHi(im,i) = msHR(im).av(auditoryTrials).cue(valid).RTloHi(i);
    end
end

bxStats = struct;
% bxStats.randGeneratorSeed = rng;
bxStats.nTrialsPerSessionRange = [min(cell2mat(nTrialsPerExpt)),...
    max(cell2mat(nTrialsPerExpt))];
bxStats.nTrialsPerMouse = cellfun(@length,{msSumStruct.tVisTargets});
bxStats.nSessionsPerMouse = cellfun(@length,nTrialsPerExpt);

pctInvPerSession = cell2mat(cellfun(@(x,y) x./y,nInvPerExpt,nTrialsPerExpt,'unif',0));
bxStats.pctInv = mean(pctInvPerSession(pctInvPerSession > 0));
bxStats.pctInvErr = ste(pctInvPerSession(pctInvPerSession > 0),2);

bxStats.allAttnBinomialTest = allAttnP;
bxStats.allAttnPowerTest = allAttnPower;
bxStats.sessionAttnTTest = cell2mat(sessionAttnTest);
bxStats.sessionPowerTest = sessionAttnPowerTest;
bxStats.attnMiceInd = bxStats.allAttnBinomialTest < attnTestAlpha & bxStats.allAttnPowerTest;

fprintf('%s +/- %s (%s-%s) sessions per mouse\n', num2str(mean(bxStats.nSessionsPerMouse)),...
    num2str(ste(bxStats.nSessionsPerMouse,2)),num2str(min(bxStats.nSessionsPerMouse)),...
    num2str(max(bxStats.nSessionsPerMouse)))
fprintf('%s +/- %s (%s-%s) sessions per attn mouse\n', ...
    num2str(mean(bxStats.nSessionsPerMouse(bxStats.attnMiceInd))),...
    num2str(ste(bxStats.nSessionsPerMouse(bxStats.attnMiceInd),2)),...
    num2str(min(bxStats.nSessionsPerMouse(bxStats.attnMiceInd))),...
    num2str(max(bxStats.nSessionsPerMouse(bxStats.attnMiceInd))))
fprintf('%s +/- %s (%s-%s) sessions per mouse\n', num2str(mean(bxStats.nTrialsPerMouse)),...
    num2str(ste(bxStats.nTrialsPerMouse,2)),num2str(min(bxStats.nTrialsPerMouse)),...
    num2str(max(bxStats.nTrialsPerMouse)))

bxStats.av(visualTrials).targets = visTargets;
bxStats.av(auditoryTrials).targets = audTargets;

bxStats.av(visualTrials).lapseRate = mean(visLR);
bxStats.av(auditoryTrials).lapseRate = mean(audLR);
bxStats.av(visualTrials).lapseRateErr = ste(visLR,2);
bxStats.av(auditoryTrials).lapseRateErr = ste(audLR,2);

bxStats.av(visualTrials).falseAlarmRate = mean(visFAR);
bxStats.av(auditoryTrials).falseAlarmRate = mean(audFAR);
bxStats.av(visualTrials).falseAlarmRateErr = ste(visFAR,2);
bxStats.av(auditoryTrials).falseAlarmRateErr = ste(audFAR,2);

[~,hiLoAttnTest] = ttest(allLoHiHRDiff(1,bxStats.attnMiceInd),...
    allLoHiHRDiff(2,bxStats.attnMiceInd));
bxStats.allMatchedHRDiff = allVisAudHR(:,valid) - allVisAudHR(:,invalid);
bxStats.allMatchedLoHiHRDiff = allLoHiHRDiff;
bxStats.loHiAttnTest = hiLoAttnTest;

fprintf('Mean HR Diff/Err: %s/%s\n',...
    num2str(mean(bxStats.allMatchedHRDiff(bxStats.attnMiceInd),1)),...
    num2str(ste(bxStats.allMatchedHRDiff(bxStats.attnMiceInd),1)))
fprintf('HR Diff Range: %s-%s\n',...
    num2str(min(bxStats.allMatchedHRDiff(bxStats.attnMiceInd))),...
    num2str(max(bxStats.allMatchedHRDiff(bxStats.attnMiceInd))))
fprintf('Mean HR Hard Trials Diff/Err: %s/%s\n',...
    num2str(mean(bxStats.allMatchedLoHiHRDiff(1,bxStats.attnMiceInd),2)),...
    num2str(ste(bxStats.allMatchedLoHiHRDiff(1,bxStats.attnMiceInd),2)))
% fprintf('HR Hard Trials Diff Range: %s-%s\n',...
%     num2str(min(bxStats.allMatchedLoHiHRDiff(1,bxStats.attnMiceInd))),...
%     num2str(max(bxStats.allMatchedLoHiHRDiff(1,bxStats.attnMiceInd))))
fprintf('Mean HR Easy Trials Diff/Err: %s/%s\n',...
    num2str(mean(bxStats.allMatchedLoHiHRDiff(2,bxStats.attnMiceInd),2)),...
    num2str(ste(bxStats.allMatchedLoHiHRDiff(2,bxStats.attnMiceInd),2)))
% fprintf('HR Easy Trials Diff Range: %s-%s\n',...
%     num2str(min(bxStats.allMatchedLoHiHRDiff(2,bxStats.attnMiceInd))),...
%     num2str(max(bxStats.allMatchedLoHiHRDiff(2,bxStats.attnMiceInd))))

bxStats.av(visualTrials).matchedHRDiff = visHR(:,valid) - visHR(:,invalid);
bxStats.av(auditoryTrials).matchedHRDiff = audHR(:,valid) - audHR(:,invalid);

fprintf('Mean Visual HR Diff/Err: %s/%s\n',...
    num2str(mean(bxStats.av(visualTrials).matchedHRDiff(bxStats.attnMiceInd),1)),...
    num2str(ste(bxStats.av(visualTrials).matchedHRDiff(bxStats.attnMiceInd),1)))
fprintf('Mean Visual HR Diff Range: %s/%s\n',...
    num2str(min(bxStats.av(visualTrials).matchedHRDiff(bxStats.attnMiceInd))),...
    num2str(max(bxStats.av(visualTrials).matchedHRDiff(bxStats.attnMiceInd))))
fprintf('Mean Auditory HR Diff/Err: %s/%s\n',...
    num2str(mean(bxStats.av(auditoryTrials).matchedHRDiff(bxStats.attnMiceInd),1)),...
    num2str(ste(bxStats.av(auditoryTrials).matchedHRDiff(bxStats.attnMiceInd),1)))
fprintf('Mean Auditory HR Diff Range: %s/%s\n',...
    num2str(min(bxStats.av(auditoryTrials).matchedHRDiff(bxStats.attnMiceInd))),...
    num2str(max(bxStats.av(auditoryTrials).matchedHRDiff(bxStats.attnMiceInd))))

bxStats.av(visualTrials).RTdiff = visRT(:,valid) - visRT(:,invalid);
[~,bxStats.av(visualTrials).RTdiffTest] = ttest(visRT(bxStats.attnMiceInd,valid),...
    visRT(bxStats.attnMiceInd,invalid),'tail','right');
bxStats.av(auditoryTrials).RTdiff = audRT(:,valid) - audRT(:,invalid);
[~,bxStats.av(auditoryTrials).RTdiffTest] = ttest(audRT(bxStats.attnMiceInd,valid),...
    audRT(bxStats.attnMiceInd,invalid),'tail','right');
bxStats.RTdiffAV = visRT(:,valid) - audRT(:,valid);
[~,bxStats.RTdiffAVTest] = ttest(visRT(bxStats.attnMiceInd,valid),...
    audRT(bxStats.attnMiceInd,valid),'tail','right');

fprintf('Mean Visual RT Diff Mean/Err: %s/%s; p=%s\n',...
    num2str(mean(bxStats.av(visualTrials).RTdiff(bxStats.attnMiceInd),1)),...
    num2str(ste(bxStats.av(visualTrials).RTdiff(bxStats.attnMiceInd),1)),...
    num2str(round(bxStats.av(visualTrials).RTdiffTest,2,'significant')))
fprintf('Mean Auditory RT Diff Mean/Err: %s/%s; p=%s\n',...
    num2str(mean(bxStats.av(auditoryTrials).RTdiff(bxStats.attnMiceInd),1)),...
    num2str(ste(bxStats.av(auditoryTrials).RTdiff(bxStats.attnMiceInd),1)),...
    num2str(round(bxStats.av(auditoryTrials).RTdiffTest,2,'significant')))
fprintf('Mean Vis-Aud RT Diff Mean/Err: %s/%s; p=%s\n',...
    num2str(mean(bxStats.RTdiffAV(bxStats.attnMiceInd),1)),...
    num2str(ste(bxStats.RTdiffAV(bxStats.attnMiceInd),1)),...
    num2str(round(bxStats.RTdiffAVTest,2,'significant')))

[~,bxStats.RTAVTest_loHi(1)] = ttest(visRT_loHi(bxStats.attnMiceInd,1),...
    audRT_loHi(bxStats.attnMiceInd,1));
[~,bxStats.RTAVTest_loHi(2)] = ttest(visRT_loHi(bxStats.attnMiceInd,2),...
    audRT_loHi(bxStats.attnMiceInd,2));

bxStats.av(visualTrials).attnTest = visAttnP;
bxStats.av(auditoryTrials).attnTest = audAttnP;

normHR_rew = nan(1,sum(bxStats.attnMiceInd));
normHR_norew = nan(1,sum(bxStats.attnMiceInd));
for im = 1:nMice
    if im == 1
        attnMiceID = 1;
    else
        attnMiceID = attnMiceID+1;        
    end
    if bxStats.attnMiceInd(im)
        valInvHRAll = msHR(im).matchedHRall./msHR(im).matchedHRall(1);
        if trainRew(im)
            normHR_rew(attnMiceID) = valInvHRAll(2);
        else
            normHR_norew(attnMiceID) = valInvHRAll(2);
        end
    end
end

% normHR_rew = rewSortMsHR(1).valInvAllTrialsHR(2)./rewSortMsHR(1).valInvAllTrialsHR(1);
% normHR_norew = rewSortMsHR(2).valInvAllTrialsHR(2)./rewSortMsHR(2).valInvAllTrialsHR(1);
% [~,p] = prop_test(...
%     [round(rewSortMsHR(1).matchedTrialN(2)*normHR_rew),round(rewSortMsHR(2).matchedTrialN(2)*normHR_norew)],...
%     [rewSortMsHR(1).matchedTrialN(2),rewSortMsHR(2).matchedTrialN(2)],0);

% bxStats.trainingRewVsNorewTest = p;
fprintf('Norm Inv. HR, training rew: %s+/-%s\n', num2str(nanmean(normHR_rew)), num2str(ste(normHR_rew,2)))
fprintf('Norm Inv. HR, training no rew: %s+/-%s\n', num2str(nanmean(normHR_norew)), num2str(ste(normHR_norew,2)))
% fprintf('Training Rew. vs No Rew., Chi-Square Test on normalized catch trial HR, p=%s\n',...
%     num2str(round(p,2,'significant')))

[~,sessionAttnTest] = cellfun(@(x) ttest(x(1,:),x(2,:),'tail','right'),matchedHRall);
bxStats.sessionAttnTest = sessionAttnTest;
disp(sessionAttnTest)

save(fullfile(fnout,'bxStats'),'bxStats')

%%

HR_pct_lim = [0 100];
HR_pct_label = [0:20:100];

HR_lim = [0 1];
HR_label = [0:0.2:1];

FAR_lim = [0 0.12];
FAR_label = 0:0.02:0.12;

RT_lim = [200 400];
RT_label = 200:50:400;

msColors = brewermap(nMice,'Dark2');

visLevels_lim = [min(targets{visualTrials,valid}(:))-1 110];
visLevels_label = [11.25 22.5 45 90];
audLevels_lim = [min(targets{auditoryTrials,valid}(:))-...
    (0.5*min(audTargets(audTargets > 0.00001))) 1.1];
audLevels_label = [0 0.001 0.01 0.1 1];

attnMiceInd = bxStats.attnMiceInd;
avLabel = {'Visual','Auditory'};
%% plot HR summary
setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',12)
fitsOris = cell(1,2);
fitsOris(:) = {nan(100,nMice)};
fitsHR = cell(1,2);
fitsHR(:) = {nan(100,nMice)};
HRfig = figure;
for im = 1:nMice
    for iav = 1:2
        x = msHR(im).av(iav).cue(valid).targets;
        xerr = msHR(im).av(iav).cue(valid).targetsErr;
        y = msHR(im).av(iav).cue(valid).HR;
        ylerr = y - (msHR(im).av(iav).cue(valid).HR95CI(:,1)');
        yuerr = (msHR(im).av(iav).cue(valid).HR95CI(:,2)') - y;
%         f = msHR(im).av(iav).cue(valid).fit;
        f = msHR(im).av(iav).cue(valid).fitWithFAR;
        fitX = msHR(im).av(iav).cue(valid).fitGrid;
        fitY = f.modelFun(f.coefEsts, fitX);
        fitsOri{iav}(:,im) = fitX;
        fitsHR{iav}(:,im) = fitY;
        subplot(4,2,iav)
        hold on
        h = plot(fitX,fitY,'-');
        h.Color = cueColor{valid};
        if ~attnMiceInd(im)
            h.LineStyle = ':';
        end
    end
end 

hiLoHR = cell(2,2);
hiLoHR(:) = {nan(nMice,2)};
allHR = cell(1,2);
allHR(:) = {nan(nMice,2)};
for im = 1:nMice
    for iav = 1:2
        for ithresh = 1:2
            x = 1:2;
            y = [msHR(im).av(iav).cue(1).hiLoHR(ithresh),...
                msHR(im).av(iav).cue(2).hiLoHR(ithresh)];
            hiLoHR{iav,ithresh}(im,:) = y;            
            if ithresh == 1
                subplot(4,2,iav+2)
                hold on
                title(sprintf('Targets < %s Threshold',num2str(highThreshold)))
                attnTestP = msHR(im).av(iav).belowThreshAttnTest;
            elseif ithresh == 2
                subplot(4,2,iav+4)
                hold on
                title(sprintf('Targets > %s Threshold',num2str(highThreshold)))
                attnTestP = msHR(im).av(iav).aboveThreshAttnTest;
            end
            h = plot(x,y,'-');
            if ~attnMiceInd(im)
                h.LineStyle = ':';
            end
            h.Color = hiLoColor{ithresh};
            if im == nMice
                y = mean(hiLoHR{iav,ithresh}(attnMiceInd,:),1);
                yerr = ste(hiLoHR{iav,ithresh}(attnMiceInd,:),1);
                h = errorbar(x,y,yerr,'.');
                h.Color = hiLoColor{ithresh};
                figXAxis([],'',[0 3],x,{'Vaild';'Invalid'})
                figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
                figAxForm
            end

        end

        x = 1:2;
        y = msHR(im).av(iav).matchedHRall;
%         p = msHR(im).av(iav).attnTest;
        allHR{iav}(im,:) = y;
        subplot(4,2,iav+6)
        hold on
        h = plot(x,y,'-');
            h.Color = 'k';
        if attnMiceInd(im)
            h.LineStyle = '-';
        else
            h.LineStyle = ':';
        end
        if im == nMice
            y = mean(allHR{iav}(attnMiceInd,:),1);
            yerr = ste(allHR{iav}(attnMiceInd,:),1);
            h = errorbar(x,y,yerr,'.');
            h.Color = 'k';
            figXAxis([],'',[0 3],x,{'Vaild';'Invalid'})
            figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
            figAxForm
            title(sprintf('All Trials, alpha = %s',num2str(round(attnTestAlpha,2,'significant'))))
        end
    end
end
subplot(4,2,visualTrials)
ax = gca;
ax.XScale = 'log';
figXAxis([],'Orientation Change (deg)',visLevels_lim,visLevels_label,visLevels_label)
figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
figAxForm
title('Visual Trials')
hold on
hline(highThreshold,'k:')

subplot(4,2,auditoryTrials)
ax = gca;
ax.XScale = 'log';
figXAxis([],'Tone Volume (?)',audLevels_lim,audLevels_label,audLevels_label)
figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
figAxForm
title('Auditory Trials')
hold on
hline(highThreshold,'k:')
print(fullfile(fnout,'bxSummary_allMice'),'-dpdf','-fillpage')

%matched hit rate across all visual and auditory trials
HRallTrials = [];
hiLoHRallTrials = cell(1,2);
figure
suptitle({'Attn Test across all (matched) visual and auditory trials';...
    'Solid lines are mice that have signif. diference across all trials types';...
    'Points are average across "attention" mice'})
subplot 131
for im = 1:nMice
    subplot 131
    x = 1:2;
    y = msHR(im).matchedHRall;
    HRallTrials(im,:) = y;  
    hold on
    h = plot(x,y,'-');
    if ~attnMiceInd(im)
        h.LineStyle = ':';
    end
    h.Color = 'k';
    title('All Trials')
    if im == nMice
        y = mean(HRallTrials(attnMiceInd,:),1);
        yerr = ste(HRallTrials(attnMiceInd,:),1);
        h = errorbar(x,y,yerr,'.');
        h.Color = hiLoColor{ithresh};
        figXAxis([],'',[0 3],x,{'Vaild';'Invalid'})
        figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
        figAxForm
    end
    for ithresh = 1:2
        subplot(1,3,1+ithresh)
        x = 1:2;
        y = msHR(im).matchedHRhiLo(ithresh,:);
        hiLoHRallTrials{ithresh}(im,:) = y;  
        attnTestP = msHR(im).attnTestHiLo(ithresh);          
        if ithresh == 1
            hold on
            title(sprintf('Targets < %s Threshold',num2str(highThreshold)))
        elseif ithresh == 2
            hold on
            title(sprintf('Targets > %s Threshold',num2str(highThreshold)))
        end
        h = plot(x,y,'-');
        if ~attnMiceInd(im)
            h.LineStyle = ':';
        end
        h.Color = hiLoColor{ithresh};
        if im == nMice
            y = mean(hiLoHRallTrials{ithresh}(attnMiceInd,:),1);
            yerr = ste(hiLoHRallTrials{ithresh}(attnMiceInd,:),1);
            h = errorbar(x,y,yerr,'.');
            h.Color = hiLoColor{ithresh};
            figXAxis([],'',[0 3],x,{'Vaild';'Invalid'})
            figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
            figAxForm
        end
    end
end
print(fullfile(fnout,'matchedHR_allTrials_allMice'),'-dpdf','-fillpage')


%false alarm and lapse rate
FAR_LR_RT_fig = figure;
subplot 521
FAR = nan(nMice,2);
for im = 1:nMice
    x = 1:2;
    y = cat(2,msHR(im).av(visualTrials).FAR, msHR(im).av(auditoryTrials).FAR);
    hold on
    h = plot(x,y,'k-');
    if attnMiceInd(im)
        h.LineStyle = '-';
    else
        h.LineStyle = ':';
    end
    FAR(im,:) = y;
end
y = mean(FAR(attnMiceInd,:),1);
yerr = ste(FAR(attnMiceInd,:),1);
errorbar(x,y,yerr,'k.')
figXAxis([],'',[0 3],x,{'Vis';'Aud'})
figYAxis([],'FA Rate',FAR_lim,FAR_label,FAR_label)
figAxForm
[~,p] = ttest(FAR(attnMiceInd,1),FAR(attnMiceInd,2));
title(sprintf('p = %s',num2str(round(p,2,'significant'))))
subplot 522
LR = nan(nMice,2);
for im = 1:nMice
    x = 1:2;
    y = cat(2,1 - msHR(im).av(visualTrials).cue(valid).HR(end)./100,...
        1 - msHR(im).av(auditoryTrials).cue(valid).HR(end)./100);
    LR(im,:) = y;
    hold on
    h = plot(x,y,'k-');
    if attnMiceInd(im)
        h.LineStyle = '-';
    else
        h.LineStyle = ':';
    end
end
y = mean(LR(attnMiceInd,:),1);
yerr = ste(LR(attnMiceInd,:),1);
errorbar(x,y,yerr,'k.')
figXAxis([],'',[0 3],x,{'Vis';'Aud'})
figYAxis([],'Lapse Rate',FAR_lim,FAR_label,FAR_label)
figAxForm
[~,p] = ttest(LR(attnMiceInd,1),LR(attnMiceInd,2));
title(sprintf('p = %s',num2str(round(p,2,'significant'))))

%reaction time
RTanovaTestAV = nan(1,2);
for iav = 1:2
    subplot(5,2,iav+2)
    RT = nan(sum(attnMiceInd),2);
    for im = 1:nMice
        if ~attnMiceInd(im)
            continue
        end
        x = 1:2;
        y = cat(2,msHR(im).av(iav).cue(valid).RT, ...
            msHR(im).av(iav).cue(invalid).RT);
        hold on
        h = plot(x,y,'k-');
%         h.Color = msColors(im,:);
        if attnMiceInd(im)
            h.LineStyle = '-';
        else
            h.LineStyle = ':';
        end
        RT(im,:) = y;
    end
    y = mean(squeeze(RT(:,:)),1);
    yerr = ste(squeeze(RT(:,:)),1);
    h = errorbar(x,y,yerr,'k.');
    figXAxis([],'',[0 3],x,{'Val';'Inv'})
    figYAxis([],'Reaction Time (ms)',RT_lim,RT_label,RT_label)
    figAxForm
    [~,p] = ttest(RT(:,1),RT(:,2),'tail','right');
    title(sprintf('%s, p=%s',avLabel{iav},num2str(p)))
    subplot(5,2,iav+4)
    RT = nan(sum(attnMiceInd),nBins);
    RTtargets = nan(nMice,nBins);
    for im = 1:nMice
        if ~attnMiceInd(im)
            continue
        end
        x = msHR(im).av(iav).cue(valid).RTtargets;
        y = msHR(im).av(iav).cue(valid).RTbinned;
        hold on
        h = plot(x,y,'k.-');
%         h.Color = msColors(im,:);
        RT(im,:) = y;
        RTtargets(im,:) = x;
        if attnMiceInd(im)
            h.LineStyle = '-';
        else
            h.LineStyle = ':';
        end
        if msHR(im).av(iav).cue(valid).RTanovaTest > RTanovaAlpha
            h.Color = [0.5 0.5 0.5];
        end
    end
    p = anova1(RT,[],'off');
    RTanovaTestAV(iav) = p;
    figure(FAR_LR_RT_fig)
    subplot(5,2,iav+4)
    if iav == 1
        title(sprintf('Visual Trials, ANOVA p=%s',num2str(round(p,2,'significant'))))
        ax = gca;
        ax.XScale = 'log';
        figXAxis([],'Orientation Change (deg)',...
            visLevels_lim,visLevels_label,visLevels_label)
    else
        title(sprintf('Auditory Trials, ANOVA p=%s',num2str(round(p,2,'significant'))))
        ax = gca;
        ax.XScale = 'log';
        figXAxis([],'Tone Volume (?)',...
            audLevels_lim,audLevels_label,audLevels_label)
%         legend(ms2analyze(attnMiceInd))
    end
    figYAxis([],'Reaction Time (ms)',RT_lim,RT_label,RT_label)
    figAxForm
    subplot(5,2,iav+6)
    y = nanmean(RT,1);
    yerr = ste(RT,1);
    x = nanmean(RTtargets,1);
    xerr = ste(RTtargets,1);
    errorbar(x,y,yerr,yerr,xerr,xerr,'k.')
    if iav == 1
        title('Visual Trials')
        ax = gca;
        ax.XScale = 'log';
        figXAxis([],'Orientation Change (deg)',...
            visLevels_lim,visLevels_label,visLevels_label)
    else
        title('Auditory Trials')
        ax = gca;
        ax.XScale = 'log';
        figXAxis([],'Tone Volume (?)',...
            audLevels_lim,audLevels_label,audLevels_label)
%         legend(ms2analyze(attnMiceInd))
    end
    figYAxis([],'Reaction Time (ms)',RT_lim,RT_label,RT_label)
    figAxForm
    
end
% subplot 529
% RT = nan(sum(attnMiceInd),2);
% for im = 1:nMice
%     if ~attnMiceInd(im)
%         continue
%     end
%     for iav = 1:2
%         RT(im,iav) = msHR(im).av(iav).cue(valid).RT;
%     end
%     hold on
%     plot(1:2,RT(im,:),'k-');
% end
% y = mean(RT,1);
% yerr = ste(RT,1);
% errorbar(1:2,y,yerr,'k')
% figXAxis([],'',[0 3],1:2,{'Vis';'Aud'})
% figYAxis([],'Reaction Time (ms)',RT_lim,RT_label,RT_label)
% figAxForm
for irt = 1:2
    subplot(5,2,8+irt)
    y = cat(2,visRT_loHi(bxStats.attnMiceInd,irt),audRT_loHi(bxStats.attnMiceInd,irt));
    yerr = ste(y,1);
    plot(1:2,y,'k-')
    hold on
    errorbar(1:2,mean(y,1),yerr,'.')
    figXAxis([],'',[0 3],1:2,{'Vis';'Aud'})
    figYAxis([],'Reaction Time (ms)',RT_lim,RT_label,RT_label)
    figAxForm
    if irt == 1
        title(sprintf('Below %s Thresh., p=%s',num2str(highThreshold),...
            num2str(round(bxStats.RTAVTest_loHi(irt),2,'significant'))))
    else
        title(sprintf('Above %s Thresh., p=%s',num2str(highThreshold),...
            num2str(round(bxStats.RTAVTest_loHi(irt),2,'significant'))))
    end
end

print(fullfile(fnout,'FAR_LR_RT_allMice'),'-dpdf','-fillpage')

% reward vs. no reward
avName = {'Visual';'Auditory'};
if doPlot
    setFigParams4Print('portrait')
    set(0,'defaultAxesFontSize',10)
    figure
    suptitle('Compare Training Types')
    for im = 1:2
        for iav = 1:2
            iplot = iav + ((im-1)*2);
    %         iplot2 = iplot + 2;
            subplot(5,2,iplot)
            for icue = 1:2
                x = rewSortMsHR(im).av(iav).cue(icue).targets;
                xerr = rewSortMsHR(im).av(iav).cue(icue).targetsErr;
                y = rewSortMsHR(im).av(iav).cue(icue).HR;
                yerrl = y' - rewSortMsHR(im).av(iav).cue(icue).HR95ci(:,1);
                yerru = rewSortMsHR(im).av(iav).cue(icue).HR95ci(:,2) - y';

                hold on
                h = errorbar(x,y,yerrl,yerru,xerr,xerr,'o-');
                h.Color = cueColor{icue};
                h.LineStyle = rewardedLine{im};

                if icue == 1
                    f = rewSortMsHR(im).av(iav).fit;
                    minI = min(x);
                    maxI = max(x);
                    fitX = logspace(log10(minI*0.1),log10(maxI*1.5),100);
                    fitY = f.modelFun(f.coefEsts, fitX).*100;
                    h = plot(fitX,fitY,'-');
                    h.Color = cueColor{valid};
                    plot([f.thresh f.thresh],[0 f.threshY.*100],'r-')
                    text(f.thresh.*1.1,f.threshY.*100.*0.5,...
                        num2str(round(f.thresh,2,'significant')))
                    highThreshX = fitX(find(fitY >= highThreshold.*100,1));
                    plot([highThreshX highThreshX],[0 highThreshold.*100],'b-')
                    text(highThreshX.*1.1,highThreshold.*100.*0.5,...
                        num2str(round(highThreshX,2,'significant')))
                end
            end
            if iav == 1
                ax = gca;
                ax.XScale = 'log';
                figXAxis([],'Orientation Change (deg)',visLevels_lim,visLevels_label,visLevels_label)
                figYAxis([],'Hit Rate (%)',HR_pct_lim,HR_pct_label,HR_pct_label)
                figAxForm
            else
                ax = gca;
                ax.XScale = 'log';
                figXAxis([],'Tone Volume (?)',audLevels_lim,audLevels_label,audLevels_label)
                figYAxis([],'Hit Rate (%)',HR_pct_lim,HR_pct_label,HR_pct_label)
                figAxForm
            end
            hline(highThreshold.*100,'k--')
            title(sprintf('%s Trials,%s',avName{iav},rewSortMsHR(im).name))
        end
    end

    for im = 1:2
        for iav = 1:2
            iplot = iav+4;
            iplot2 = iplot + 2;
            iplot3 = iplot + 4;

            subplot(5,2,iplot)
            normHR = cat(2,rewSortMsHR(im).av(iav).cue(1).hiLoHR(1)./rewSortMsHR(im).av(iav).cue(1).hiLoHR(1),...
                rewSortMsHR(im).av(iav).cue(2).hiLoHR(1)./rewSortMsHR(im).av(iav).cue(1).hiLoHR(1));
            hold on
            h = plot(1:2,normHR,'o-');
            h.Color = cueColor{1};
            h.LineStyle = rewardedLine{im};
            text(3,normHR(2),num2str(round(normHR(2),2,'significant')))
            figXAxis([],'',[0 4],1:3,{'Vaild';'Invalid';'norm HR'})
            figYAxis([],'Norm. Hit Rate',[0 1],HR_pct_label./100,HR_pct_label./100)
            figAxForm
            if im == 2
                normHR_ms1 = cat(2,rewSortMsHR(1).av(iav).cue(1).hiLoHR(1)./rewSortMsHR(1).av(iav).cue(1).hiLoHR(1),...
                    rewSortMsHR(1).av(iav).cue(2).hiLoHR(1)./rewSortMsHR(1).av(iav).cue(1).hiLoHR(1));
                n1 = rewSortMsHR(1).av(iav).cue(2).nHiLo(1);
                n2 = rewSortMsHR(2).av(iav).cue(2).nHiLo(1);
                nHit1 = round(n1.*normHR_ms1(2),2,'significant');
                nHit2 = round(n2.*normHR(2),2,'significant');
                [~,p] = prop_test([nHit1 nHit2],[n1 n2],0);
                 title({sprintf('%s, Below %s% Threshold',avName{iav},...
                    num2str(highThreshold.*100));...
                    sprintf('Inv HR p = %s (chi-square test)',num2str(round(p,2,'significant')))})   
            end

            subplot(5,2,iplot2)
            normHR = cat(2,rewSortMsHR(im).av(iav).cue(1).hiLoHR(2)./rewSortMsHR(im).av(iav).cue(1).hiLoHR(2),...
                rewSortMsHR(im).av(iav).cue(2).hiLoHR(2)./rewSortMsHR(im).av(iav).cue(1).hiLoHR(2));
            hold on
            h = plot(1:2,normHR,'o-');
            h.Color = cueColor{1};
            h.LineStyle = rewardedLine{im};
            text(3,normHR(2),num2str(round(normHR(2),2,'significant')))
            figXAxis([],'',[0 4],1:3,{'Vaild';'Invalid';'norm HR'})
            figYAxis([],'Norm. Hit Rate',[0 1],HR_pct_label./100,HR_pct_label./100)
            figAxForm
            if im == 2
                normHR_ms1 = cat(2,rewSortMsHR(1).av(iav).cue(1).hiLoHR(2)./rewSortMsHR(1).av(iav).cue(1).hiLoHR(2),...
                    rewSortMsHR(1).av(iav).cue(2).hiLoHR(2)./rewSortMsHR(1).av(iav).cue(1).hiLoHR(2));
                n1 = rewSortMsHR(1).av(iav).cue(2).nHiLo(2);
                n2 = rewSortMsHR(2).av(iav).cue(2).nHiLo(2);
                nHit1 = round(n1.*normHR_ms1(2),2,'significant');
                nHit2 = round(n2.*normHR(2),2,'significant');
                [~,p] = prop_test([nHit1 nHit2],[n1 n2],0);
                 title({sprintf('%s, Above %s% Threshold',avName{iav},...
                    num2str(highThreshold.*100));...
                    sprintf('Inv HR p = %s (chi-square test)',num2str(round(p,2,'significant')))})   
            end

            subplot(5,2,iplot3)
            normHR = rewSortMsHR(im).av(iav).valInvAllTrialsHR./...
                rewSortMsHR(im).av(iav).valInvAllTrialsHR(1);
            hold on
            h = plot(1:2,normHR,'o-');
            h.Color = cueColor{1};
            h.LineStyle = rewardedLine{im};
            text(3,normHR(2),num2str(round(normHR(2),2,'significant')))
            figXAxis([],'',[0 4],1:3,{'Vaild';'Invalid';'norm HR'})
            figYAxis([],'Norm. Hit Rate',[0 1],HR_pct_label./100,HR_pct_label./100)
            figAxForm
            if im == 2
                normHR_ms1 = rewSortMsHR(1).av(iav).valInvAllTrialsHR./...
                    rewSortMsHR(1).av(iav).valInvAllTrialsHR(1);
                n1 = rewSortMsHR(1).av(iav).matchedTrialN(2);
                n2 = rewSortMsHR(2).av(iav).matchedTrialN(2);
                nHit1 = round(n1.*normHR_ms1(2),2,'significant');
                nHit2 = round(n2.*normHR(2),2,'significant');
                [~,p] = prop_test([nHit1 nHit2],[n1 n2],0);
                 title({sprintf('%s, All Trials',avName{iav},...
                    num2str(highThreshold.*100));...
                    sprintf('Inv HR p = %s (chi-square test)',num2str(round(p,2,'significant')))})   
            end
        end
    end
    print(fullfile(fnout,'compareTrainingTypes'),'-dpdf','-fillpage')
end

%% add some stats to structure
bxStats.av(visualTrials).RTanova = RTanovaTestAV(visualTrials);
bxStats.av(auditoryTrials).RTanova = RTanovaTestAV(auditoryTrials);

save(fullfile(fnout,'bxStats'),'bxStats')
%% example mouse
sessionAttnFig = figure;
[nSessRows,nSessCols] = optimizeSubplotDim(nMice);
for im = 1:nMice
    mouseName = ms2analyze{im};
    fn = fullfile(rc.ashleyAnalysis,mouseName,'behavior');
    load(fullfile(fn,[mouseName,'bxSummary_dataAnalyzed_attnV1ms']))
    msExptInfo = msExptAnalyzed;

    nexp = size(msExptInfo,2);
    exExptHR = cell(2,2);
    exExptHR(:) = {nan(2,nexp)};
    exExptHighThresh = cell(2,2);
    exExptHighThresh(:) = {nan(2,nexp)};
    for i = 1:nexp 
        for iav = 1:2
            for icue = 1:2
                exExptHR{iav,icue}(:,i) = msExptInfo(i).av(iav).cue(icue).highThreshHR;
            end
        end
    end

    if doPlot
        exMsFig = figure;
        suptitle(ms2analyze(im))
        for iav = 1:2
            for icue = 1:2
                x = msHR(im).av(iav).cue(icue).targets;
                xerr = msHR(im).av(iav).cue(icue).targetsErr;
                y = msHR(im).av(iav).cue(icue).HR;
                ylerr = y -  msHR(im).av(iav).cue(icue).HR95CI(:,1)';
                yuerr = msHR(im).av(iav).cue(icue).HR95CI(:,2)' - y;

                subplot(4,2,iav)
                hold on
                h = errorbar(x,y,ylerr,yuerr,xerr,xerr,'.');
                h.Color = cueColor{icue};
                h.LineWidth = 1;
                h.MarkerFaceColor = [1 1 1];

                if icue == valid
    %                 f = msHR(im).av(iav).cue(icue).fit;
                    f = msHR(im).av(iav).cue(icue).fitWithFAR;
                    x = msHR(im).av(iav).cue(icue).fitGrid;
                    y = f.modelFun(f.coefEsts, x).*100;
                    subplot(4,2,iav)
                    hold on
                    h = plot(x,y,'-');
                    h.LineWidth = 1;
                    h.Color = cueColor{icue};
                end
            end
        end
        subplot(4,2,visualTrials)
        ax = gca;
        ax.XScale = 'log';
        figXAxis([],'Orientation Change (deg)',visLevels_lim,visLevels_label,visLevels_label)
        figYAxis([],'Hit Rate (%)',HR_pct_lim,HR_pct_label,HR_pct_label)
        figAxForm
        title('Visual Trials - Combined Expt')
           hold on
           hline(highThreshold*100,'k:')

        subplot(4,2,auditoryTrials)
        ax = gca;
        ax.XScale = 'log';
        figXAxis([],'Tone Volume (?)',audLevels_lim,audLevels_label,audLevels_label)
        figYAxis([],'Hit Rate (%)',HR_pct_lim,HR_pct_label,HR_pct_label)
        figAxForm
        title('Auditory Trials - Combined Expt')
       hold on
       hline(highThreshold*100,'k:')
        for iav = 1:2
            for ithresh = 1:2
                subplot(4,2,iav+2)
                hold on
                x = exExptHR{iav,1}(ithresh,:);
                y = exExptHR{iav,2}(ithresh,:);
                ind = ~isnan(y);
                h = plot(x,y,'o');
                h.Color = hiLoColor{ithresh};
                h.MarkerFaceColor = hiLoColor{ithresh};
                h = errorbar(nanmean(x(ind)),nanmean(y(ind)),...
                    ste(y(ind),2),ste(y(ind),2),ste(x(ind),2),ste(x(ind),2),'o');
                h.Color = hiLoColor{ithresh};
                h.MarkerFaceColor = [1 1 1];
                h.LineWidth = 1;
            end
            plot(HR_lim,HR_lim,'k--')
            figXAxis([],'Valid Hit Rate (%)',HR_lim,HR_label,HR_label)
            figYAxis([],'Invalid Hit Rate (%)',HR_lim,HR_label,HR_label)
            figAxForm
            if iav == 1
                title('Visual Trials - Each Expt')
            elseif iav == 2
                title('Auditory Trials - Each Expt')
            end
            
            subplot(4,2,iav+4)
            x = cell2mat(cellfun(@(x) x(valid,:),...
                matchedHR(iav,im),'unif',0));
            y = cell2mat(cellfun(@(x) x(invalid,:),...
            matchedHR(iav,im),'unif',0));h = plot(x,y,'o');
            ind = ~isnan(y);
            h.Color = hiLoColor{ithresh};
            h = plot(x,y,'k.');
            hold on
            h = errorbar(nanmean(x(ind)),nanmean(y(ind)),...
                ste(y(ind),2),ste(y(ind),2),ste(x(ind),2),ste(x(ind),2),'k.');
            h.LineWidth = 1;
            plot(HR_lim,HR_lim,'k--')
            figXAxis([],'Valid Hit Rate (%)',HR_lim,HR_label,HR_label)
            figYAxis([],'Invalid Hit Rate (%)',HR_lim,HR_label,HR_label)
            figAxForm
            if iav == 1
                title('Visual Trials - Each Expt')
            elseif iav == 2
                title('Auditory Trials - Each Expt')
            end
        end
        subplot(4,2,7)
        x = matchedHRall{im}(valid,:);
        y = matchedHRall{im}(invalid,:);
        ind = ~isnan(y);
        h.Color = hiLoColor{ithresh};
        h = plot(x,y,'k.');
        hold on
        h = errorbar(nanmean(x(ind)),nanmean(y(ind)),...
            ste(y(ind),2),ste(y(ind),2),ste(x(ind),2),ste(x(ind),2),'k.');
        h.LineWidth = 1;
        plot(HR_lim,HR_lim,'k--')
        figXAxis([],'Valid Hit Rate (%)',HR_lim,HR_label,HR_label)
        figYAxis([],'Invalid Hit Rate (%)',HR_lim,HR_label,HR_label)
        figAxForm
        title('All Matched Trials')
        print(fullfile(fnout,['bxSummary_exampleMouse_' ms2analyze{im}]),'-dpdf','-fillpage')
        
        figure(sessionAttnFig)
        subplot(nSessRows,nSessCols,im)
        h = plot(x,y,'k.');
        [~,p] = ttest(x,y,'tail','right');
        hold on
        h = errorbar(nanmean(x(ind)),nanmean(y(ind)),...
            ste(y(ind),2),ste(y(ind),2),ste(x(ind),2),ste(x(ind),2),'k.');
        h.LineWidth = 1;
        plot(HR_lim,HR_lim,'k--')
        figXAxis([],'Valid Hit Rate (%)',HR_lim,HR_label,HR_label)
        figYAxis([],'Invalid Hit Rate (%)',HR_lim,HR_label,HR_label)
        figAxForm
        title(sprintf('%s, p=%s',ms2analyze{im},num2str(round(p,2,'significant'))))
    end
end
figure(sessionAttnFig)
print(fullfile(fnout,'bxSummary_sessionAttn_allMice'),'-dpdf','-fillpage')