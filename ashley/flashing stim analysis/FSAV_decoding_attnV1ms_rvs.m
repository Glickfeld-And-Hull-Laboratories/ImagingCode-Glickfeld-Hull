clear all
close all
ds = 'FSAV_attentionV1_noAttn';
analysisDate = '201016'; %200527
titleStr = ds(6:end);
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms
fn = fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
    [titleStr '_rvs_' datestr(now,'yymmdd') '_']); 
load([fn(1:end-7) analysisDate '_imgAnalysisData'])
dc_noAttn = decodeAnalysis;
siPerExpt_noAttn = siPerExpt;
respCellsExpt_noAttn = respCellsExpt;
oriTuning_noAttn = oriTuningExpt;
mouseName_noAttn = {respCellsExpt.mouse};

ds = 'FSAV_V1_naive_GCaMP6m';
analysisDate = '201019';
titleStr = ds(6:end);
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms
fn = fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
    [titleStr '_rvs_' datestr(now,'yymmdd') '_']); 
load([fn(1:end-7) analysisDate '_imgAnalysisData'])
dc_naive = decodeAnalysis;
siPerExpt_naive = siPerExpt;
respCellsExpt_naive = respCellsExpt;
oriTuning_naive = oriTuningExpt;
mouseName_naive = {respCellsExpt.mouse};

ds = 'FSAV_attentionV1';
analysisDate = '201016'; %200527
titleStr = ds(6:end);
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms
fn = fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
    [titleStr '_rvs_' datestr(now,'yymmdd') '_']); 
load([fn(1:end-7) analysisDate '_imgAnalysisData'])
dc_attn = decodeAnalysis;
siPerExpt_attn = siPerExpt;
respCellsExpt_attn = respCellsExpt;
oriTuning_attn = oriTuningExpt;
mouseName_attn = {respCellsExpt.mouse};

fnout = fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
    [titleStr '_' datestr(now,'yymmdd') '_decoding_']);

%%
load(fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs','bxStats.mat'))
    
%% indexing and naming variables
avName = {'Vis','Aud'};
targetName = {'D','HT','ET'};
modelName = {'Stimulus','Choice'};

nexp_attn = size(dc_attn,2);
nexp_noAttn = size(dc_noAttn,2);
nexp_naive = size(dc_naive,2);

stimModInd = 1;
choiceModInd = 2;
%% performance all trials
%model performance - attention mice
mdl_attn = struct;
for iav = 1:2
    mdl_attn(iav).name = avName{iav};
    mdl_attn(iav).mdl(stimModInd).name = 'Stimulus';
    mdl_attn(iav).mdl(choiceModInd).name = 'Choice';
    mdl_attn(iav).mdl(stimModInd).pctCorr_all = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorr_all = nan(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).pctCorr_all(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectAllTarget_holdout;
        mdl_attn(iav).mdl(choiceModInd).pctCorr_all(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectAllDetect_holdout;
    end
end
goodPerfInd_attn = mdl_attn(visualTrials).mdl(stimModInd).pctCorr_all...
    > pctCorrThresh;
%model performance - no attention mice
mdl_noAttn = struct;
for iav = 1:2
    mdl_noAttn(iav).name = avName{iav};
    mdl_noAttn(iav).mdl(stimModInd).name = 'Stimulus';
    mdl_noAttn(iav).mdl(choiceModInd).name = 'Choice';
    mdl_noAttn(iav).mdl(stimModInd).pctCorr_all = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all = nan(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(stimModInd).pctCorr_all(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectAllTarget_holdout;
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectAllDetect_holdout;
    end
end
goodPerfInd_noAttn = mdl_noAttn(visualTrials).mdl(stimModInd).pctCorr_all...
    > pctCorrThresh;

%model performance - naive mice
mdl_naive = struct;
for iav = 1:2
    mdl_naive(iav).name = avName{iav};
    mdl_naive(iav).mdl(stimModInd).name = 'Stimulus';
    mdl_naive(iav).mdl(stimModInd).pctCorr_all = nan(1,nexp_naive);
    for iexp = 1:nexp_naive
        mdl_naive(iav).mdl(stimModInd).pctCorr_all(iexp) = ...
            dc_naive(iexp).av(iav).pctCorrectAllTarget_holdout;
    end
end
goodPerfInd_naive = mdl_naive(visualTrials).mdl(stimModInd).pctCorr_all...
    > pctCorrThresh;

mdlPerf_fig = figure;

suptitle(sprintf(...
    'All Trials Performance (held-out data), pct corr must be > %s in vis stim model',...
    num2str(pctCorrThresh)))
for iav = 1:2
    for imod = 1:2
        if iav == 1
            iplot = imod;
        else
            iplot = imod+2;
        end
        subplot(2,2,iplot)
        y1 = mdl_attn(iav).mdl(imod).pctCorr_all(goodPerfInd_attn);
        plot(ones(1,length(y1)),y1,'k.','MarkerSize',10)
        hold on
        y2 = mdl_noAttn(iav).mdl(imod).pctCorr_all(goodPerfInd_noAttn);
        plot(ones(1,length(y2)).*2,y2,'k.','MarkerSize',10)
        if imod == 1
            x = 1:3;
            y3 = mdl_naive(iav).mdl(imod).pctCorr_all(goodPerfInd_naive);
            plot(ones(1,length(y3))*3,y3,'k.','MarkerSize',10);
            y_all = [mean(y1),mean(y2),mean(y3)];
            y_all_ste = [ste(y1,2),ste(y2,2),ste(y3,2)];
            figXAxis([],'',[0 4],1:3,{'Attn','No Attn','Naive'})
        else
            x = 1:2;
            y_all = [mean(y1),mean(y2)];
            y_all_ste = [ste(y1,2),ste(y2,2)];
            figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
        end
        errorbar(x,y_all,y_all_ste,'.','MarkerSize',20)
        figYAxis([],'% Correct',[0 1],0:0.2:1)
        figAxForm([],0)
        hline(pctCorrThresh,'k:')
        hline(0.5,'k-')
        [~,p] = ttest2(y1,y2);
        title(sprintf('%s %s Model, p=%s',mdl_attn(iav).name,...
            mdl_attn(iav).mdl(imod).name,num2str(round(p,2,'significant'))))
    end
end

print([fnout 'pctCorrect'],'-dpdf','-fillpage')

fprintf('Attn mice - vis stim model successful expt: %spct %s/%s expt\n',...
    sigfigString(sum(goodPerfInd_attn)/length(goodPerfInd_attn)*100),...
    num2str(sum(goodPerfInd_attn)),num2str(length(goodPerfInd_attn)))
fprintf('No Attn mice - vis stim model successful expt: %spct %s/%s expt\n',...
    sigfigString(sum(goodPerfInd_noAttn)/length(goodPerfInd_noAttn)*100),...
    num2str(sum(goodPerfInd_noAttn)),num2str(length(goodPerfInd_noAttn)))
fprintf('Naive mice - vis stim model successful expt: %spct %s/%s expt\n',...
    sigfigString(sum(goodPerfInd_naive)/length(goodPerfInd_naive)*100),...
    num2str(sum(goodPerfInd_naive)),num2str(length(goodPerfInd_naive)))
%% compare stim and choice performance
setFigParams4Print('portrait')
figure
suptitle(sprintf(...
    'pct corr must be > %s in vis stim model',...
    num2str(pctCorrThresh)))
for iav = 1:2
    subplot(2,2,iav)
    hold on
    y1 = mdl_attn(iav).mdl(stimModInd).pctCorr_all(goodPerfInd_attn);
    y2 = mdl_attn(iav).mdl(choiceModInd).pctCorr_all(goodPerfInd_attn);
    plot(1:2,[y1',y2'],'k-')
    errorbar(1:2,mean([y1',y2']),ste([y1',y2'],1),'.','MarkerSize',20)
    figXAxis([],'Model',[0 3],1:2,{'Stim','Choice'})
    figYAxis([],'Fraction Correct',[0 1],0:0.2:1)
    figAxForm([],0)
    hline(pctCorrThresh,'k:')
    hline(0.5,'k-')
    [~,p] = ttest(y1,y2);
    title(sprintf('Attn, %s Models, p=%s',avName{iav},sigfigString(p)))
    
    subplot(2,2,iav+2)
    hold on
    y1 = mdl_noAttn(iav).mdl(stimModInd).pctCorr_all(goodPerfInd_noAttn);
    y2 = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all(goodPerfInd_noAttn);
    plot(1:2,[y1',y2'],'k-')
    errorbar(1:2,mean([y1',y2']),ste([y1',y2'],1),'.','MarkerSize',20)
    figXAxis([],'Model',[0 3],1:2,{'Stim','Choice'})
    figYAxis([],'Fraction Correct',[0 1],0:0.2:1)
    figAxForm([],0)
    hline(pctCorrThresh,'k:')
    hline(0.5,'k-')
    [~,p] = ttest(y1,y2);
    title(sprintf('No Attn, %s Models, p=%s',avName{iav},sigfigString(p)))
end    
print([fnout 'pctCorrect_stimVsChoice'],'-dpdf','-fillpage') 
%% performance across stim types
for iav = 1:2
    mdl_attn(iav).mdl(stimModInd).pctCorr_xstim = nan(3,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorr_xstim = nan(3,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).pctCorr_xstim(:,iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectXStimTarget_holdout;
        mdl_attn(iav).mdl(choiceModInd).pctCorr_xstim(:,iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectXStimDetect_holdout;
    end
    mdl_noAttn(iav).mdl(stimModInd).pctCorr_xstim = nan(3,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_xstim = nan(3,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(stimModInd).pctCorr_xstim(:,iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectXStimTarget_holdout;
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_xstim(:,iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectXStimDetect_holdout;
    end
    
    mdl_naive(iav).mdl(stimModInd).pctCorr_xstim = nan(3,nexp_naive);
    for iexp = 1:nexp_naive
        mdl_naive(iav).mdl(stimModInd).pctCorr_xstim(:,iexp) = ...
            dc_naive(iexp).av(iav).pctCorrectXStimTarget_holdout;
    end    
end

mdlPerfxStim_fig = figure;
suptitle(sprintf(...
    'Ea. Stim Performance, pct corr must be > %s within that model type',...
    num2str(pctCorrThresh)))
x=1:3;
for iav = 1:2
    for imod = 1:2
        if iav == 1
            iplot = imod;
        else
            iplot = imod+2;
        end
        subplot(3,4,iplot)
        hold on
        y = mdl_attn(iav).mdl(imod).pctCorr_xstim;
        good_ind = mdl_attn(iav).mdl(imod).pctCorr_all > pctCorrThresh;
        ind = sum(isnan(y),2)<minTrN;
        for i = 1:size(y,2)
            if ~good_ind(i)
                continue
            end
            nanInd = ~isnan(y(:,i)) & ind;
            plot(x(nanInd),y(nanInd,i),'k-')
        end
        errorbar(x(ind),nanmean(y(ind,good_ind),2),ste(y(ind,good_ind),2),'.-','MarkerSize',20)
        figXAxis([],'Test',[0 4],x,targetName)
        figYAxis([],'% Correct',[0 1],0:0.2:1)
        figAxForm
        hline(pctCorrThresh,'k:')
        hline(0.5,'k-')
        [~,pChance] = ttest(y',0.5);
        anovaTest = anova1(y',[],'off');
        title({'Attn Mice';sprintf('%s %s Model, anovaP=%s',mdl_attn(iav).name,...
            mdl_attn(iav).mdl(imod).name,num2str(round(anovaTest,2,'significant')))})
        for i = 1:3
            text(i,0.2,num2str(round(pChance(i),2,'significant')))
        end
        text(0,0.2,'pChnc=')
        
        if iav == 1
            iplot = imod+4;
        else
            iplot = imod+6;
        end
        subplot(3,4,iplot)
        hold on
        y = mdl_noAttn(iav).mdl(imod).pctCorr_xstim;
        good_ind = mdl_noAttn(iav).mdl(imod).pctCorr_all > pctCorrThresh;
        ind = sum(isnan(y),2)<minTrN;
        for i = 1:size(y,2)
            if ~good_ind(i)
                continue
            end
            nanInd = ~isnan(y(:,i)) & ind;
            plot(x(nanInd),y(nanInd,i),'k-')
        end
        errorbar(x(ind),nanmean(y(ind,:),2),ste(y(ind,:),2),'.-','MarkerSize',20)
        figXAxis([],'Test',[0 4],x,targetName)
        figYAxis([],'% Correct',[0 1],0:0.2:1)
        figAxForm
        hline(pctCorrThresh,'k:')
        hline(0.5,'k-')
        [~,pChance] = ttest(y',0.5);
        anovaTest = anova1(y',[],'off');
        title({'No Attn Mice';sprintf('%s %s Model, anovaP=%s',mdl_noAttn(iav).name,...
            mdl_noAttn(iav).mdl(imod).name,num2str(round(anovaTest,2,'significant')))})
        for i = 1:3
            text(i,0.2,num2str(round(pChance(i),2,'significant')))
        end
        text(0,0.2,'pChnc=')
        if imod == 1
            if iav == 1
                iplot = imod+8;
            else
                iplot = imod+10;
            end
            subplot(3,4,iplot)
            hold on
            y = mdl_naive(iav).mdl(imod).pctCorr_xstim;
            good_ind = mdl_naive(iav).mdl(imod).pctCorr_all > pctCorrThresh;
            ind = sum(isnan(y),2)<minTrN;
            for i = 1:size(y,2)
                if ~good_ind(i)
                    continue
                end
                nanInd = ~isnan(y(:,i)) & ind;
                plot(x(nanInd),y(nanInd,i),'k-')
            end
            errorbar(x(ind),nanmean(y(ind,:),2),ste(y(ind,:),2),'.-','MarkerSize',20)
            figXAxis([],'Test',[0 4],x,targetName)
            figYAxis([],'% Correct',[0 1],0:0.2:1)
            figAxForm
            hline(pctCorrThresh,'k:')
            hline(0.5,'k-')
            [~,pChance] = ttest(y',0.5);
            anovaTest = anova1(y',[],'off');
            title({'Naive Mice';sprintf('%s %s Model, anovaP=%s',mdl_naive(iav).name,...
                mdl_naive(iav).mdl(imod).name,num2str(round(anovaTest,2,'significant')))})
            for i = 1:3
                text(i,0.2,num2str(round(pChance(i),2,'significant')))
            end
            text(0,0.2,'pChnc=')
        end
    end 
end

print([fnout 'pctCorrectXStim'],'-dpdf','-fillpage')
%% cross-modality tested model performance
for iav = 1:2
    mdl_attn(iav).mdl(stimModInd).pctCorr_otherAV = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorr_otherAV = nan(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).pctCorr_otherAV(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectTarget_otherAV;
        mdl_attn(iav).mdl(choiceModInd).pctCorr_otherAV(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_otherAV;
    end
    mdl_noAttn(iav).mdl(stimModInd).pctCorr_otherAV = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_otherAV = nan(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(stimModInd).pctCorr_otherAV(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectTarget_otherAV;
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_otherAV(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_otherAV;
    end
    
    mdl_naive(iav).mdl(stimModInd).pctCorr_otherAV = nan(1,nexp_naive);
    for iexp = 1:nexp_naive
        mdl_naive(iav).mdl(stimModInd).pctCorr_otherAV(iexp) = ...
            dc_naive(iexp).av(iav).pctCorrectTarget_otherAV;
    end
end

xModalMdlPerf_fig = figure;
suptitle(sprintf(...
    'Cross-modality model test, pct corr must be > %s in vis stim model',...
    num2str(pctCorrThresh)))
xModel_compAttn_fig = figure;
for iav = 1:2
    for imod = 1:2
        if iav == 1
            iplot = imod;
            otherAVInd = 2;
        else
            iplot = imod+2;
            otherAVInd = 1;
        end
        figure(xModalMdlPerf_fig)
        subplot(2,4,iplot)
        ind = mdl_attn(iav).mdl(imod).pctCorr_all...
            > pctCorrThresh;
        y1 = mdl_attn(iav).mdl(imod).pctCorr_all(ind);
        y2 = mdl_attn(iav).mdl(imod).pctCorr_otherAV(ind);
        diff_attn = y1-y2;
        plot(1:2,[y1',y2'],'k-')
        hold on
        errorbar(1:2,[mean(y1),mean(y2)],[ste(y1,2),ste(y2,2)],...
            '.','MarkerSize',20)
        figXAxis([],'Test',[0 3],1:2,...
            {mdl_attn(iav).name,mdl_attn(otherAVInd).name})
        figYAxis([],'% Correct',[0 1],0:0.2:1)
        figAxForm([],0)
        hline(pctCorrThresh,'k:')
        hline(0.5,'k-')
        [~,p] = ttest2(y1,y2);
        title({'Attn Mice';sprintf('Train %s, %s Model, p=%s',mdl_attn(iav).name,...
            mdl_attn(iav).mdl(imod).name,num2str(round(p,2,'significant')))})
        [~,p_chance] = ttest([y1',y2'],0.5);
        text([0:2],ones(1,3)*.2,{'pCh',num2str(round(p_chance(1),2,'significant')),...
            num2str(round(p_chance(2),2,'significant'))})
        
        if iav == 1
            iplot = imod+4;
        else
            iplot = imod+6;
        end
        subplot(2,4,iplot)
        ind = mdl_noAttn(iav).mdl(imod).pctCorr_all...
            > pctCorrThresh;
        y1 = mdl_noAttn(iav).mdl(imod).pctCorr_all(ind);
        y2 = mdl_noAttn(iav).mdl(imod).pctCorr_otherAV(ind);
        diff_noAttn = y1-y2;
        plot(1:2,[y1',y2'],'k-')
        hold on
        errorbar(1:2,[mean(y1),mean(y2)],[ste(y1,2),ste(y2,2)],...
            '.','MarkerSize',20)
        figXAxis([],'Test',[0 3],1:2,...
            {mdl_noAttn(iav).name,mdl_noAttn(otherAVInd).name})
        figYAxis([],'% Correct',[0 1],0:0.2:1)
        figAxForm([],0)
        hline(pctCorrThresh,'k:')
        hline(0.5,'k-')
        [~,p] = ttest2(y1,y2);
        title({'No Attn Mice';sprintf('Train %s, %s Model, p=%s',mdl_noAttn(iav).name,...
            mdl_noAttn(iav).mdl(imod).name,num2str(round(p,2,'significant')))})
        [~,p_chance] = ttest([y1',y2'],0.5);
        text([0:2],ones(1,3)*.2,{'pCh',num2str(round(p_chance(1),2,'significant')),...
            num2str(round(p_chance(2),2,'significant'))})
        
        figure(xModel_compAttn_fig)
        if iav == 1
            iplot = imod;
        else
            iplot = imod+2;
        end
        subplot(2,2,iplot)
        hold on
        plot(ones(1,length(diff_attn)),diff_attn,'k.')
        errorbar(1,mean(diff_attn),ste(diff_attn,2),'.')
        plot(ones(1,length(diff_noAttn)).*2,diff_noAttn,'k.')
        errorbar(2,mean(diff_noAttn),ste(diff_noAttn,2),'.')
        [~,p] = ttest2(diff_attn,diff_noAttn);
        title(sprintf('Train %s, Test %s, %s Model,p=%s',avName{iav},avName{otherAVInd},...
            modelName{imod},sigfigString(p)))
        figXAxis([],'Test',[0 3],1:2,...
            {'Attn','No Attn'})
        figYAxis([],'% Correct',[-1 1])
        figAxForm
        hline(0,'k-')
    end 
end
figure(xModalMdlPerf_fig)
print([fnout 'pctCorrectOtherAVModel'],'-dpdf','-fillpage')

figure
for iav = 1:2
    if iav == 1
        otherAVInd = 2;
    else
        otherAVInd = 1;
    end
    subplot(1,2,iav)
    ind = mdl_naive(iav).mdl(stimModInd).pctCorr_all...
            > pctCorrThresh;
    y1 = mdl_naive(iav).mdl(stimModInd).pctCorr_all(ind);
    y2 = mdl_naive(iav).mdl(stimModInd).pctCorr_otherAV(ind);
    plot(1:2,[y1',y2'],'k-')
    hold on
    errorbar(1:2,[mean(y1),mean(y2)],[ste(y1,2),ste(y2,2)],...
        '.','MarkerSize',20)
    figXAxis([],'Test',[0 3],1:2,...
        {mdl_naive(iav).name,mdl_naive(otherAVInd).name})
    figYAxis([],'% Correct',[0 1],0:0.2:1)
    figAxForm([],0)
    hline(pctCorrThresh,'k:')
    hline(0.5,'k-')
    [~,p] = ttest2(y1,y2);
    title({'Naive Mice';sprintf('Train %s, %s Model, p=%s',mdl_naive(iav).name,...
        mdl_naive(iav).mdl(stimModInd).name,num2str(round(p,2,'significant')))})
    [~,p_chance] = ttest([y1',y2'],0.5);
    text([0:2],ones(1,3)*.2,{'pCh',num2str(round(p_chance(1),2,'significant')),...
        num2str(round(p_chance(2),2,'significant'))})
end

print([fnout 'pctCorrectOtherAVModel_naive'],'-dpdf','-fillpage')

%% weights
weightEdges = [-2.8:0.2:4.8];
nShuff = 1000;

for iav = 1:2
    mdl_attn(iav).mdl(stimModInd).weights = cell(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).weights = cell(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).weights{iexp} = ...
            dc_attn(iexp).av(iav).weightTarget;
        mdl_attn(iav).mdl(choiceModInd).weights{iexp} = ...
            dc_attn(iexp).av(iav).weightDetect;
    end
   
    mdl_noAttn(iav).mdl(stimModInd).weights = cell(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).weights = cell(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(stimModInd).weights{iexp} = ...
            dc_noAttn(iexp).av(iav).weightTarget;
        mdl_noAttn(iav).mdl(choiceModInd).weights{iexp} = ...
            dc_noAttn(iexp).av(iav).weightDetect;
    end
end
goodPerfInd_av_attn = mdl_attn(visualTrials).mdl(stimModInd).pctCorr_all...
    > pctCorrThresh |...
    mdl_attn(auditoryTrials).mdl(stimModInd).pctCorr_all > pctCorrThresh |...
    mdl_attn(visualTrials).mdl(choiceModInd).pctCorr_all > pctCorrThresh |...
    mdl_attn(auditoryTrials).mdl(choiceModInd).pctCorr_all > pctCorrThresh;
goodPerfInd_av_noAttn = mdl_noAttn(visualTrials).mdl(stimModInd).pctCorr_all...
    > pctCorrThresh |...
    mdl_noAttn(auditoryTrials).mdl(stimModInd).pctCorr_all > pctCorrThresh |...
    mdl_noAttn(visualTrials).mdl(choiceModInd).pctCorr_all > pctCorrThresh |...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).pctCorr_all > pctCorrThresh;

for iav = 1:2
    for imod = 1:2
        rng(0)
        mdl_attn(iav).mdl(imod).shuffWeights = cellfun(...
            @(x) x(randperm(length(x))),mdl_attn(iav).mdl(imod).weights,'unif',0);      
        mdl_noAttn(iav).mdl(imod).shuffWeights = cellfun(...
            @(x) x(randperm(length(x))),mdl_noAttn(iav).mdl(imod).weights,'unif',0);
        if iav == 1
            for ishuff = 1:nShuff
                shuffWeights_vis = cellfun(@(x) x(randperm(length(x))),...
                    mdl_attn(visualTrials).mdl(imod).weights(goodPerfInd_av_attn),'unif',0)';
                shuffWeights_aud = cellfun(@(x) x(randperm(length(x))),...
                    mdl_attn(auditoryTrials).mdl(imod).weights(goodPerfInd_av_attn),'unif',0)';
                if ishuff == 1
                    mdl_attn(visualTrials).mdl(imod).weightAVangle_shuff = nan(sum(cellfun(@length,shuffWeights_vis)),nShuff);
                end
                mdl_attn(visualTrials).mdl(imod).weightAVangle_shuff(:,ishuff) = cell2mat(cellfun(@(x,y) ...
                    rad2deg(abs(atan(x./y))),shuffWeights_aud,...
                    shuffWeights_vis,'unif',0));
%                 
                shuffWeights_vis = cellfun(@(x) x(randperm(length(x))),...
                    mdl_noAttn(visualTrials).mdl(imod).weights(goodPerfInd_av_noAttn),'unif',0)';
                shuffWeights_aud = cellfun(@(x) x(randperm(length(x))),...
                    mdl_noAttn(auditoryTrials).mdl(imod).weights(goodPerfInd_av_noAttn),'unif',0)';
                if ishuff == 1
                    mdl_noAttn(visualTrials).mdl(imod).weightAVangle_shuff = nan(sum(cellfun(@length,shuffWeights_vis)),nShuff);
                end
                mdl_noAttn(visualTrials).mdl(imod).weightAVangle_shuff(:,ishuff) = cell2mat(cellfun(@(x,y) ...
                    rad2deg(abs(atan(x./y))),shuffWeights_aud,...
                    shuffWeights_vis,'unif',0));
            end
        end
    end    
end

for iav = 1:2
    for imod = 1:2
        rng(0)
        if iav == 1
            for ishuff = 1:nShuff
                shuffWeights_vis = cellfun(@(x) x(randperm(length(x))),...
                    mdl_attn(visualTrials).mdl(imod).weights(goodPerfInd_av_attn),'unif',0)';
                shuffWeights_aud = cellfun(@(x) x(randperm(length(x))),...
                    mdl_attn(auditoryTrials).mdl(imod).weights(goodPerfInd_av_attn),'unif',0)';
                if ishuff == 1
                    mdl_attn(visualTrials).mdl(imod).weightAVangle_noAbs_shuff = nan(sum(cellfun(@length,shuffWeights_vis)),nShuff);
                end
                mdl_attn(visualTrials).mdl(imod).weightAVangle_noAbs_shuff(:,ishuff) = cell2mat(cellfun(@(x,y) ...
                    rad2deg(atan(x./y)),shuffWeights_aud,...
                    shuffWeights_vis,'unif',0));
%                 
                shuffWeights_vis = cellfun(@(x) x(randperm(length(x))),...
                    mdl_noAttn(visualTrials).mdl(imod).weights(goodPerfInd_av_noAttn),'unif',0)';
                shuffWeights_aud = cellfun(@(x) x(randperm(length(x))),...
                    mdl_noAttn(auditoryTrials).mdl(imod).weights(goodPerfInd_av_noAttn),'unif',0)';
                if ishuff == 1
                    mdl_noAttn(visualTrials).mdl(imod).weightAVangle_noAbs_shuff = nan(sum(cellfun(@length,shuffWeights_vis)),nShuff);
                end
                mdl_noAttn(visualTrials).mdl(imod).weightAVangle_noAbs_shuff(:,ishuff) = cell2mat(cellfun(@(x,y) ...
                    rad2deg(atan(x./y)),shuffWeights_aud,...
                    shuffWeights_vis,'unif',0));
                
            end
        end
    end    
end

weights_AVcorr_attn = cell(1,2);
weights_AVcorr_noAttn = cell(1,2);
for imod = 1:2
    weights_AVcorr_attn{imod} = cellfun(@(x,y) corr(x,y),...
        mdl_attn(visualTrials).mdl(imod).weights(goodPerfInd_av_attn),...
        mdl_attn(auditoryTrials).mdl(imod).weights(goodPerfInd_av_attn));
    weights_AVcorr_noAttn{imod} = cellfun(@(x,y) corr(x,y),...
        mdl_noAttn(visualTrials).mdl(imod).weights(goodPerfInd_av_noAttn),...
        mdl_noAttn(auditoryTrials).mdl(imod).weights(goodPerfInd_av_noAttn));
end
[~,p_wAVCorr]=cellfun(@(x,y) ttest2(x,y),weights_AVcorr_attn,weights_AVcorr_noAttn);

figure
for imod = 1:2
    subplot(1,2,imod)
    y1 = cell2mat(mdl_attn(visualTrials).mdl(imod).weights(goodPerfInd_av_attn)');
    y2 = cell2mat(mdl_attn(auditoryTrials).mdl(imod).weights(goodPerfInd_av_attn)');
    n_attn = length(y1);
%     [r_attn,p_attn,rCIl_attn,rCIu_attn] = corrcoef(y1,y2);
    [r_attn,p_attn] = corr(y1,y2);
    y1 = cell2mat(mdl_noAttn(visualTrials).mdl(imod).weights(goodPerfInd_av_noAttn)');
    y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(imod).weights(goodPerfInd_av_noAttn)');
    n_noAttn = length(y1);
    [r_noAttn,p_noAttn] = corr(y1,y2);
    bar(1,r_attn)
    text(1,0.1,sigfigString(p_attn))
    hold on
    bar(2,r_noAttn)
    text(2,0.1,sigfigString(p_noAttn))
    figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
    figYAxis([],'Vis-Aud Corr',[0 1])
    figAxForm
    title(sprintf('%s Model, p=%s',modelName{imod},...
        sigfigString(compare_correlation_coefficients(...
        r_attn,n_attn,r_noAttn,n_noAttn))))
end
print([fnout 'avWeightCorr'],'-dpdf')

figure
for iav = 1:2
    for imod = 1:2
        if iav == 1
            iplot = 0;
        else
            iplot = 2;
        end
        subplot(4,2,imod+iplot)
        y = cell2mat(mdl_attn(iav).mdl(imod).weights(goodPerfInd_av_attn)');
        histogram(y,weightEdges,'Normalization','probability')
        figXAxis([],'Weight',[-1.8 3.2])
        figYAxis([],'Fraction of Cells',[0 .7])
        figAxForm
        title({'Attn Mice';sprintf('%s %s Model',...
            avName{iav},modelName{imod})})
        vline(0,'k:')
        subplot(4,2,imod+iplot+4)
        y = cell2mat(mdl_noAttn(iav).mdl(imod).weights(goodPerfInd_av_noAttn)');
        histogram(y,weightEdges,'Normalization','probability')
        figXAxis([],'Weight',[weightEdges(1) weightEdges(end)])
        figYAxis([],'Fraction of Cells',[0 .7])
        figAxForm
        title({'No Attn Mice';sprintf('%s %s Model',...
            avName{iav},modelName{imod})})
        vline(0,'k:')
    end
end
print([fnout 'weightHist'],'-dpdf','-fillpage')

figure
for iav = 1:2
    for imod = 1:2
        if iav == 1
            iplot = 0;
        else
            iplot = 2;
        end
        subplot(2,2,imod+iplot)
        y1 = cell2mat(mdl_attn(iav).mdl(imod).weights(goodPerfInd_av_attn)');
        y2 = cell2mat(mdl_noAttn(iav).mdl(imod).weights(goodPerfInd_av_noAttn)');
        cdfplot(y1)
        hold on
        cdfplot(y2)
        figXAxis([],'Weight',[-1 2])
        figYAxis([],'Fraction of Cells',[0 1])
        figAxForm
        title(sprintf('%s %s Model',...
            avName{iav},modelName{imod}))
        vline(0,'k:')
        legend({'Attn','No Attn'},'location','southeast')
    end
end
print([fnout 'weightCDF'],'-dpdf','-fillpage')

mdlWeightsSub_fig = figure;
mdlWeights_fig = figure;
suptitle({'Weights - AV angles without abs val';...
    sprintf('pct corr must be > %s in vis & aud stim models',...
    num2str(pctCorrThresh))})
oriEdges_noAbs = [-90:15:90];
for imod = 1:2
    figure(mdlWeights_fig)
    subplot(3,4,imod)
    y1 = cell2mat(mdl_attn(visualTrials).mdl(imod).weights(goodPerfInd_av_attn)');
    y2 = cell2mat(mdl_attn(auditoryTrials).mdl(imod).weights(goodPerfInd_av_attn)');
    plot(y1,y2,'k.','MarkerSize',10)
    hold on
    mdl_fitline = fitlm(y1,y2);
    rsq = mdl_fitline.Rsquared.Ordinary;
    yfit = predict(mdl_fitline,y1);
    plot(y1,yfit,'-')
    figXAxis([],'Visual Weight',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Auditory Weight',[weightEdges(1) weightEdges(end)])
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    r = corr(y1,y2);
    title({'Attn Mice';sprintf('%s Model, r=%s',...
        mdl_attn(iav).mdl(imod).name,num2str(round(r,2,'significant')))})
    
    subplot(3,4,imod+4)
    weightAVangle = rad2deg(atan(y2./y1));
    histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
    figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
    figYAxis([],'Fraction of Cells',[0 .4])
    figAxForm
    title(sprintf('Attn Mice %s Model',...
        mdl_attn(iav).mdl(imod).name))
        
    subplot(3,4,imod+8)
    avAngleShuff = mdl_attn(visualTrials).mdl(imod).weightAVangle_noAbs_shuff;
    h = nan(length(oriEdges_noAbs)-1,nShuff);
    for ishuff = 1:nShuff
        h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges_noAbs,'Normalization','Probability');
    end    
    h_sort = sort(h,2);
    y = mean(h_sort,2);
    yerrl = y-h_sort(:,round(nShuff*0.025));
    yerru = h_sort(:,round(nShuff*0.975))-y;
    shadedErrorBar_chooseColor(oriEdges_noAbs(1:end-1),y,[yerrl';yerru'],[0 0 0]);
    hold on
    plot(oriEdges_noAbs(1:end-1),histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability'),'.-')
    controlSub_attn = histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability') - y';
    figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
    figYAxis([],'Fraction of Cells',[0 .4])
    figAxForm
    title(sprintf('Attn Mice-Shuff, %s Model',...
        mdl_attn(iav).mdl(imod).name))
    
    subplot(3,4,imod+2)
    y1 = cell2mat(mdl_noAttn(visualTrials).mdl(imod).weights(goodPerfInd_av_noAttn)');
    y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(imod).weights(goodPerfInd_av_noAttn)');
    plot(y1,y2,'k.','MarkerSize',10)
    hold on
    mdl_fitline = fitlm(y1,y2);
    rsq = mdl_fitline.Rsquared.Ordinary;
    yfit = predict(mdl_fitline,y1);
    plot(y1,yfit,'-')
    figXAxis([],'Visual Weight',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Auditory Weight',[weightEdges(1) weightEdges(end)])
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    r = corr(y1,y2);
    title({'No Attn Mice';sprintf('%s Model, r=%s',...
        mdl_noAttn(iav).mdl(imod).name,num2str(round(r,2,'significant')))})
    
    subplot(3,4,imod+6)
    weightAVangle = rad2deg(atan(y2./y1));
    histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
    figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
    figYAxis([],'Fraction of Cells',[0 .4])
    figAxForm
    title(sprintf('No Attn Mice, %s Model',...
        mdl_attn(iav).mdl(imod).name))
       
    subplot(3,4,imod+10)
    avAngleShuff = mdl_noAttn(visualTrials).mdl(imod).weightAVangle_noAbs_shuff;
    h = nan(length(oriEdges_noAbs)-1,nShuff);
    for ishuff = 1:nShuff
        h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges_noAbs,'Normalization','Probability');
    end    
    h_sort = sort(h,2);
    y = mean(h_sort,2);
    yerrl = y-h_sort(:,round(nShuff*0.025));
    yerru = h_sort(:,round(nShuff*0.975))-y;
    shadedErrorBar_chooseColor(oriEdges_noAbs(1:end-1),y,[yerrl';yerru'],[0 0 0]);
    hold on
    plot(oriEdges_noAbs(1:end-1),histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability'),'.-')
    controlSub_noAttn = histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability') - y';
    figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
    figYAxis([],'Fraction of Cells',[0 .4])
    figAxForm
    title(sprintf('Attn Mice-Shuff, %s Model',...
        mdl_noAttn(iav).mdl(imod).name))
    
    figure(mdlWeightsSub_fig)
    subplot(1,2,imod)
    plot(oriEdges_noAbs(1:end-1),controlSub_attn,'.-')
    hold on
%     shadedErrorBar_chooseColor(oriEdges(1:end-1),controlSub_noAttn,[yerrl';yerru'],[0.5 0.5 0.5]);
    plot(oriEdges_noAbs(1:end-1),controlSub_noAttn,'.-')
    figXAxis([],'Angle of Norm. A/V Weight',[-100 100],-90:15:90)
    figYAxis([],'Control Subtr. Fraction of Cells',[-0.15 0.15])
    figAxForm
    hline(0,'k:')
    title(sprintf('%s Model',modelName{imod}))
       
end 
figure(mdlWeights_fig)
print([fnout 'modelWeights_noAbs'],'-dpdf','-fillpage')
figure(mdlWeightsSub_fig)
print([fnout 'modelWeights_noAbs_ciSub'],'-dpdf')

%% stim and choice weight correlation with bootstrap quantification
nboot = 1000;
figure
for imod = 1:2
    subplot(2,4,imod)
    y1 = cell2mat(mdl_attn(visualTrials).mdl(imod).weights(goodPerfInd_av_attn)');
    y2 = cell2mat(mdl_attn(auditoryTrials).mdl(imod).weights(goodPerfInd_av_attn)');
    plot(y1,y2,'k.','MarkerSize',10)
    hold on
    mdl_fitline = fitlm(y1,y2);
    yfit = predict(mdl_fitline,y1);
    plot(y1,yfit,'-')
    figXAxis([],'Visual Weight',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Auditory Weight',[weightEdges(1) weightEdges(end)])
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    title({'Attn Mice';sprintf('%s Model',...
        mdl_attn(iav).mdl(imod).name)})
    
    r_attn = nan(1,nboot);
    for iboot = 1:nboot
        n = length(y1);
        ind = randsample(1:n,n,1);
        y1_samp = y1(ind);
        y2_samp = y2(ind);
        r_attn(iboot) = corr(y1_samp,y2_samp);
    end
    
    subplot(2,4,imod+2)
    y1 = cell2mat(mdl_noAttn(visualTrials).mdl(imod).weights(goodPerfInd_av_noAttn)');
    y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(imod).weights(goodPerfInd_av_noAttn)');
    plot(y1,y2,'k.','MarkerSize',10)
    hold on
    mdl_fitline = fitlm(y1,y2);
    yfit = predict(mdl_fitline,y1);
    plot(y1,yfit,'-')
    figXAxis([],'Visual Weight',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Auditory Weight',[weightEdges(1) weightEdges(end)])
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    title({'Attn Mice';sprintf('%s Model',...
        mdl_noAttn(iav).mdl(imod).name)})
    
    r_noAttn = nan(1,nboot);
    for iboot = 1:nboot
        n = length(y1);
        ind = randsample(1:n,n,1);
        y1_samp = y1(ind);
        y2_samp = y2(ind);
        r_noAttn(iboot) = corr(y1_samp,y2_samp);
    end
    
    subplot(2,4,4+imod); hold on
    bar(1:2,[mean(r_attn),mean(r_noAttn)])
    [ylerr,yuerr] = ciFromBoot(r_attn',95);
    errorbar(1,mean(r_attn),ylerr,yuerr,'.','MarkerSize',20)
    [ylerr,yuerr] = ciFromBoot(r_noAttn',95);
    errorbar(2,mean(r_noAttn),ylerr,yuerr,'.','MarkerSize',20)
    figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
    xtickangle(-45)
    figYAxis([],'Stim:Choice Weight Corr',[0 1])
    figAxForm
    title(sprintf('%s Model',mdl_noAttn(iav).mdl(imod).name))
end
print([fnout 'modelWeights_bootCorr'],'-dpdf')
%% weights from matched model
weightEdges = [-2.8:0.2:4.8];
nShuff = 1000;

for iav = 1:2
    mdl_attn(iav).mdl(stimModInd).weights_matchAll = cell(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).weights_matchAll = cell(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).weights_matchAll{iexp} = ...
            dc_attn(iexp).av(iav).weightTarget_matchAll;
        mdl_attn(iav).mdl(choiceModInd).weights_matchAll{iexp} = ...
            dc_attn(iexp).av(iav).weightDetect_matchAll;
    end
   
    mdl_noAttn(iav).mdl(stimModInd).weights_matchAll = cell(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).weights_matchAll = cell(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(stimModInd).weights_matchAll{iexp} = ...
            dc_noAttn(iexp).av(iav).weightTarget_matchAll;
        mdl_noAttn(iav).mdl(choiceModInd).weights_matchAll{iexp} = ...
            dc_noAttn(iexp).av(iav).weightDetect_matchAll;
    end
end

for iav = 1:2
    for imod = 1:2
        rng(0)
        mdl_attn(iav).mdl(imod).shuffWeights_matchAll = cellfun(...
            @(x) x(randperm(length(x))),mdl_attn(iav).mdl(imod).weights_matchAll,'unif',0);      
        mdl_noAttn(iav).mdl(imod).shuffWeights_matchAll = cellfun(...
            @(x) x(randperm(length(x))),mdl_noAttn(iav).mdl(imod).weights_matchAll,'unif',0);
        if iav == 1
            for ishuff = 1:nShuff
                ind = cellfun(@(x) ~isempty(x),mdl_attn(visualTrials).mdl(imod).weights_matchAll) & ...
                    cellfun(@(x) ~isempty(x),mdl_attn(auditoryTrials).mdl(imod).weights_matchAll);
                shuffWeights_vis = cellfun(@(x) x(randperm(length(x))),...
                    mdl_attn(visualTrials).mdl(imod).weights_matchAll(goodPerfInd_av_attn & ind),'unif',0)';
                shuffWeights_aud = cellfun(@(x) x(randperm(length(x))),...
                    mdl_attn(auditoryTrials).mdl(imod).weights_matchAll(goodPerfInd_av_attn & ind),'unif',0)';
                if ishuff == 1
                    mdl_attn(visualTrials).mdl(imod).weightAVangle_shuff_matchAll = nan(sum(cellfun(@length,shuffWeights_vis)),nShuff);
                end
                mdl_attn(visualTrials).mdl(imod).weightAVangle_shuff_matchAll(:,ishuff) = cell2mat(cellfun(@(x,y) ...
                    rad2deg(abs(atan(x./y))),shuffWeights_aud,...
                    shuffWeights_vis,'unif',0));
%               
                ind = cellfun(@(x) ~isempty(x),mdl_noAttn(visualTrials).mdl(imod).weights_matchAll) & ...
                    cellfun(@(x) ~isempty(x),mdl_noAttn(auditoryTrials).mdl(imod).weights_matchAll);  
                shuffWeights_vis = cellfun(@(x) x(randperm(length(x))),...
                    mdl_noAttn(visualTrials).mdl(imod).weights_matchAll(goodPerfInd_av_noAttn & ind),'unif',0)';
                shuffWeights_aud = cellfun(@(x) x(randperm(length(x))),...
                    mdl_noAttn(auditoryTrials).mdl(imod).weights_matchAll(goodPerfInd_av_noAttn & ind),'unif',0)';
                if ishuff == 1
                    mdl_noAttn(visualTrials).mdl(imod).weightAVangle_shuff_matchAll = nan(sum(cellfun(@length,shuffWeights_vis)),nShuff);
                end
                mdl_noAttn(visualTrials).mdl(imod).weightAVangle_shuff_matchAll(:,ishuff) = cell2mat(cellfun(@(x,y) ...
                    rad2deg(abs(atan(x./y))),shuffWeights_aud,...
                    shuffWeights_vis,'unif',0));
            end
        end
    end    
end

for iav = 1:2
    for imod = 1:2
        rng(0)
        if iav == 1
            for ishuff = 1:nShuff
                ind = cellfun(@(x) ~isempty(x),mdl_attn(visualTrials).mdl(imod).weights_matchAll) & ...
                    cellfun(@(x) ~isempty(x),mdl_attn(auditoryTrials).mdl(imod).weights_matchAll);
                shuffWeights_vis = cellfun(@(x) x(randperm(length(x))),...
                    mdl_attn(visualTrials).mdl(imod).weights_matchAll(goodPerfInd_av_attn & ind),'unif',0)';
                shuffWeights_aud = cellfun(@(x) x(randperm(length(x))),...
                    mdl_attn(auditoryTrials).mdl(imod).weights_matchAll(goodPerfInd_av_attn & ind),'unif',0)';
                if ishuff == 1
                    mdl_attn(visualTrials).mdl(imod).weightAVangle_noAbs_shuff_matchAll = nan(sum(cellfun(@length,shuffWeights_vis)),nShuff);
                end
                mdl_attn(visualTrials).mdl(imod).weightAVangle_noAbs_shuff_matchAll(:,ishuff) = cell2mat(cellfun(@(x,y) ...
                    rad2deg(atan(x./y)),shuffWeights_aud,...
                    shuffWeights_vis,'unif',0));
%                 
                ind = cellfun(@(x) ~isempty(x),mdl_noAttn(visualTrials).mdl(imod).weights_matchAll) & ...
                    cellfun(@(x) ~isempty(x),mdl_noAttn(auditoryTrials).mdl(imod).weights_matchAll);
                shuffWeights_vis = cellfun(@(x) x(randperm(length(x))),...
                    mdl_noAttn(visualTrials).mdl(imod).weights_matchAll(goodPerfInd_av_noAttn & ind),'unif',0)';
                shuffWeights_aud = cellfun(@(x) x(randperm(length(x))),...
                    mdl_noAttn(auditoryTrials).mdl(imod).weights_matchAll(goodPerfInd_av_noAttn & ind),'unif',0)';
                if ishuff == 1
                    mdl_noAttn(visualTrials).mdl(imod).weightAVangle_noAbs_shuff_matchAll = nan(sum(cellfun(@length,shuffWeights_vis)),nShuff);
                end
                mdl_noAttn(visualTrials).mdl(imod).weightAVangle_noAbs_shuff_matchAll(:,ishuff) = cell2mat(cellfun(@(x,y) ...
                    rad2deg(atan(x./y)),shuffWeights_aud,...
                    shuffWeights_vis,'unif',0));
                
            end
        end
    end    
end

weights_AVcorr_attn = cell(1,2);
weights_AVcorr_noAttn = cell(1,2);
for imod = 1:2
    weights_AVcorr_attn{imod} = cellfun(@(x,y) corr(x,y),...
        mdl_attn(visualTrials).mdl(imod).weights(goodPerfInd_av_attn),...
        mdl_attn(auditoryTrials).mdl(imod).weights(goodPerfInd_av_attn));
    weights_AVcorr_noAttn{imod} = cellfun(@(x,y) corr(x,y),...
        mdl_noAttn(visualTrials).mdl(imod).weights(goodPerfInd_av_noAttn),...
        mdl_noAttn(auditoryTrials).mdl(imod).weights(goodPerfInd_av_noAttn));
end
[~,p_wAVCorr]=cellfun(@(x,y) ttest2(x,y),weights_AVcorr_attn,weights_AVcorr_noAttn);

figure
for imod = 1:2
    subplot(1,2,imod)
    ind = cellfun(@(x) ~isempty(x),mdl_attn(visualTrials).mdl(imod).weights_matchAll) & ...
        cellfun(@(x) ~isempty(x),mdl_attn(auditoryTrials).mdl(imod).weights_matchAll);        
    y1 = cell2mat(mdl_attn(visualTrials).mdl(imod).weights_matchAll(goodPerfInd_av_attn & ind)');
    y2 = cell2mat(mdl_attn(auditoryTrials).mdl(imod).weights_matchAll(goodPerfInd_av_attn & ind)');
    n_attn = length(y1);
%     [r_attn,p_attn,rCIl_attn,rCIu_attn] = corrcoef(y1,y2);
    [r_attn,p_attn] = corr(y1,y2);
    ind = cellfun(@(x) ~isempty(x),mdl_noAttn(visualTrials).mdl(imod).weights_matchAll) & ...
        cellfun(@(x) ~isempty(x),mdl_noAttn(auditoryTrials).mdl(imod).weights_matchAll);
    y1 = cell2mat(mdl_noAttn(visualTrials).mdl(imod).weights_matchAll(goodPerfInd_av_noAttn & ind)');
    y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(imod).weights_matchAll(goodPerfInd_av_noAttn & ind)');
    n_noAttn = length(y1);
    [r_noAttn,p_noAttn] = corr(y1,y2);
    bar(1,r_attn)
    text(1,0.1,sigfigString(p_attn))
    hold on
    bar(2,r_noAttn)
    text(2,0.1,sigfigString(p_noAttn))
    figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
    figYAxis([],'Vis-Aud Corr',[0 1])
    figAxForm
    title(sprintf('%s Model, p=%s',modelName{imod},...
        sigfigString(compare_correlation_coefficients(...
        r_attn,n_attn,r_noAttn,n_noAttn))))
end
print([fnout 'avWeightCorr_matchAll'],'-dpdf')

figure
for iav = 1:2
    for imod = 1:2
        if iav == 1
            iplot = 0;
        else
            iplot = 2;
        end
        subplot(4,2,imod+iplot)
        ind = cellfun(@(x) ~isempty(x),mdl_attn(visualTrials).mdl(imod).weights_matchAll) & ...
            cellfun(@(x) ~isempty(x),mdl_attn(auditoryTrials).mdl(imod).weights_matchAll);        
        y = cell2mat(mdl_attn(iav).mdl(imod).weights_matchAll(goodPerfInd_av_attn & ind)');
        histogram(y,weightEdges,'Normalization','probability')
        figXAxis([],'Weight',[-1.8 3.2])
        figYAxis([],'Fraction of Cells',[0 .7])
        figAxForm
        title({'Attn Mice';sprintf('%s %s Model',...
            avName{iav},modelName{imod})})
        vline(0,'k:')
        subplot(4,2,imod+iplot+4)
        ind = cellfun(@(x) ~isempty(x),mdl_noAttn(visualTrials).mdl(imod).weights_matchAll) & ...
            cellfun(@(x) ~isempty(x),mdl_noAttn(auditoryTrials).mdl(imod).weights_matchAll);        
        y = cell2mat(mdl_noAttn(iav).mdl(imod).weights_matchAll(goodPerfInd_av_noAttn & ind)');
        histogram(y,weightEdges,'Normalization','probability')
        figXAxis([],'Weight',[weightEdges(1) weightEdges(end)])
        figYAxis([],'Fraction of Cells',[0 .7])
        figAxForm
        title({'No Attn Mice';sprintf('%s %s Model',...
            avName{iav},modelName{imod})})
        vline(0,'k:')
    end
end
print([fnout 'weightHist_matchAll'],'-dpdf','-fillpage')

figure
for iav = 1:2
    for imod = 1:2
        if iav == 1
            iplot = 0;
        else
            iplot = 2;
        end
        subplot(2,2,imod+iplot)
        ind = cellfun(@(x) ~isempty(x),mdl_attn(visualTrials).mdl(imod).weights_matchAll) & ...
            cellfun(@(x) ~isempty(x),mdl_attn(auditoryTrials).mdl(imod).weights_matchAll);        
        y1 = cell2mat(mdl_attn(iav).mdl(imod).weights_matchAll(goodPerfInd_av_attn & ind)');
        ind = cellfun(@(x) ~isempty(x),mdl_noAttn(visualTrials).mdl(imod).weights_matchAll) & ...
            cellfun(@(x) ~isempty(x),mdl_noAttn(auditoryTrials).mdl(imod).weights_matchAll);        
        y2 = cell2mat(mdl_noAttn(iav).mdl(imod).weights_matchAll(goodPerfInd_av_noAttn & ind)');
        cdfplot(y1)
        hold on
        cdfplot(y2)
        figXAxis([],'Weight',[-1 2])
        figYAxis([],'Fraction of Cells',[0 1])
        figAxForm
        title(sprintf('%s %s Model',...
            avName{iav},modelName{imod}))
        vline(0,'k:')
        legend({'Attn','No Attn'},'location','southeast')
    end
end
print([fnout 'weightCDF_matchAll'],'-dpdf','-fillpage')

mdlWeightsSub_fig = figure;
mdlWeights_fig = figure;
suptitle({'Weights - AV angles without abs val';...
    sprintf('pct corr must be > %s in vis & aud stim models',...
    num2str(pctCorrThresh))})
oriEdges_noAbs = [-90:15:90];
for imod = 1:2
    figure(mdlWeights_fig)
    subplot(3,4,imod)
    ind = cellfun(@(x) ~isempty(x),mdl_attn(visualTrials).mdl(imod).weights_matchAll) & ...
        cellfun(@(x) ~isempty(x),mdl_attn(auditoryTrials).mdl(imod).weights_matchAll);        
    y1 = cell2mat(mdl_attn(visualTrials).mdl(imod).weights_matchAll(goodPerfInd_av_attn & ind)');
    y2 = cell2mat(mdl_attn(auditoryTrials).mdl(imod).weights_matchAll(goodPerfInd_av_attn & ind)');
    plot(y1,y2,'k.','MarkerSize',10)
    hold on
    mdl_fitline = fitlm(y1,y2);
    rsq = mdl_fitline.Rsquared.Ordinary;
    yfit = predict(mdl_fitline,y1);
    plot(y1,yfit,'-')
    figXAxis([],'Visual Weight',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Auditory Weight',[weightEdges(1) weightEdges(end)])
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    r = corr(y1,y2);
    title({'Attn Mice';sprintf('%s Model, r=%s',...
        mdl_attn(iav).mdl(imod).name,num2str(round(r,2,'significant')))})
    
    subplot(3,4,imod+4)
    weightAVangle = rad2deg(atan(y2./y1));
    histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
    figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
    figYAxis([],'Fraction of Cells',[0 .4])
    figAxForm
    title(sprintf('Attn Mice %s Model',...
        mdl_attn(iav).mdl(imod).name))
        
    subplot(3,4,imod+8)
    avAngleShuff = mdl_attn(visualTrials).mdl(imod).weightAVangle_noAbs_shuff_matchAll;
    h = nan(length(oriEdges_noAbs)-1,nShuff);
    for ishuff = 1:nShuff
        h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges_noAbs,'Normalization','Probability');
    end    
    h_sort = sort(h,2);
    y = mean(h_sort,2);
    yerrl = y-h_sort(:,round(nShuff*0.025));
    yerru = h_sort(:,round(nShuff*0.975))-y;
    shadedErrorBar_chooseColor(oriEdges_noAbs(1:end-1),y,[yerrl';yerru'],[0 0 0]);
    hold on
    plot(oriEdges_noAbs(1:end-1),histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability'),'.-')
    controlSub_attn = histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability') - y';
    figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
    figYAxis([],'Fraction of Cells',[0 .4])
    figAxForm
    title(sprintf('Attn Mice-Shuff, %s Model',...
        mdl_attn(iav).mdl(imod).name))
    
    subplot(3,4,imod+2)
    ind = cellfun(@(x) ~isempty(x),mdl_noAttn(visualTrials).mdl(imod).weights_matchAll) & ...
        cellfun(@(x) ~isempty(x),mdl_noAttn(auditoryTrials).mdl(imod).weights_matchAll);        
    y1 = cell2mat(mdl_noAttn(visualTrials).mdl(imod).weights_matchAll(goodPerfInd_av_noAttn & ind)');
    y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(imod).weights_matchAll(goodPerfInd_av_noAttn & ind)');
    plot(y1,y2,'k.','MarkerSize',10)
    hold on
    mdl_fitline = fitlm(y1,y2);
    rsq = mdl_fitline.Rsquared.Ordinary;
    yfit = predict(mdl_fitline,y1);
    plot(y1,yfit,'-')
    figXAxis([],'Visual Weight',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Auditory Weight',[weightEdges(1) weightEdges(end)])
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    r = corr(y1,y2);
    title({'No Attn Mice';sprintf('%s Model, r=%s',...
        mdl_noAttn(iav).mdl(imod).name,num2str(round(r,2,'significant')))})
    
    subplot(3,4,imod+6)
    weightAVangle = rad2deg(atan(y2./y1));
    histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
    figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
    figYAxis([],'Fraction of Cells',[0 .4])
    figAxForm
    title(sprintf('No Attn Mice, %s Model',...
        mdl_attn(iav).mdl(imod).name))
       
    subplot(3,4,imod+10)
    avAngleShuff = mdl_noAttn(visualTrials).mdl(imod).weightAVangle_noAbs_shuff_matchAll;
    h = nan(length(oriEdges_noAbs)-1,nShuff);
    for ishuff = 1:nShuff
        h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges_noAbs,'Normalization','Probability');
    end    
    h_sort = sort(h,2);
    y = mean(h_sort,2);
    yerrl = y-h_sort(:,round(nShuff*0.025));
    yerru = h_sort(:,round(nShuff*0.975))-y;
    shadedErrorBar_chooseColor(oriEdges_noAbs(1:end-1),y,[yerrl';yerru'],[0 0 0]);
    hold on
    plot(oriEdges_noAbs(1:end-1),histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability'),'.-')
    controlSub_noAttn = histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability') - y';
    figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
    figYAxis([],'Fraction of Cells',[0 .4])
    figAxForm
    title(sprintf('Attn Mice-Shuff, %s Model',...
        mdl_noAttn(iav).mdl(imod).name))
    
    figure(mdlWeightsSub_fig)
    subplot(1,2,imod)
    plot(oriEdges_noAbs(1:end-1),controlSub_attn,'.-')
    hold on
%     shadedErrorBar_chooseColor(oriEdges(1:end-1),controlSub_noAttn,[yerrl';yerru'],[0.5 0.5 0.5]);
    plot(oriEdges_noAbs(1:end-1),controlSub_noAttn,'.-')
    figXAxis([],'Angle of Norm. A/V Weight',[-100 100],-90:15:90)
    figYAxis([],'Control Subtr. Fraction of Cells',[-0.15 0.15])
    figAxForm
    hline(0,'k:')
    title(sprintf('%s Model',modelName{imod}))
       
end 
figure(mdlWeights_fig)
print([fnout 'modelWeights_noAbs_matchAll'],'-dpdf','-fillpage')
figure(mdlWeightsSub_fig)
print([fnout 'modelWeights_noAbs_ciSub_matchAll'],'-dpdf')

%% example experiment PCs
exExpt = 11;
npcs = dc_attn(exExpt).nPCs;
stimExPCs_fig = figure;
colormap(brewermap([],'Blues'))
choiceExPCs_fig = figure;
colormap(brewermap([],'Reds'))
weightsExPCs_fig = figure;
colormap gray
pcNum = cellfun(@num2str,num2cell(1:npcs),'unif',0)';
for iav = 1:2
    pcs = dc_attn(exExpt).av(iav).respAllCells(:,1:npcs);
    pc_wChoice = dc_attn(exExpt).av(iav).weightDetect_pcs(1:npcs);
    pc_wStim = dc_attn(exExpt).av(iav).weightTarget_pcs(1:npcs);
    trOut = dc_attn(exExpt).av(iav).trOut;
    [choiceTrOut, stimTrOut] = getStimAndBehaviorYs(trOut);
    choiceResp = nan(npcs,2);
    stimResp = nan(npcs,2);
    for i = 1:2
        choiceResp(:,i) = mean(pcs(choiceTrOut==(i-1),:),1);
        stimResp(:,i) = mean(pcs(stimTrOut==(i-1),:),1);
    end
    if iav == 1
%         [~,sortInd_choice] = sort(pc_wChoice);
        [~,sortInd_stim] = sort(pc_wStim);
    end
    
    pcLabel_stim = cellfun(@(x,y) [x '-' y],...
        cellfun(@num2str,num2cell(round(pc_wStim(sortInd_stim),2,'significant')),'unif',0),...
        pcNum(sortInd_stim),'unif',0);
    pcLabel_choice = cellfun(@(x,y) [x '-' y],...
        cellfun(@num2str,num2cell(round(pc_wChoice(sortInd_stim),2,'significant')),'unif',0),...
        pcNum(sortInd_stim),'unif',0);
    
    figure(stimExPCs_fig)
    subplot(1,2,iav)
    imagesc(stimResp(sortInd_stim,:))
	clim([-0.6 1.2])
    colorbar
    figXAxis([],'Stimulus',[],1:2,{'Dist','Tar'})
    figYAxis([],'PC #',[],1:npcs,pcLabel_stim)
    figAxForm([],0)
    title(avName{iav})
    
    figure(choiceExPCs_fig)
    subplot(1,2,iav)
    imagesc(choiceResp(sortInd_stim,:))
    clim([-0.6 1.2])
    colorbar
    figXAxis([],'Choice',[],1:2,{'No','Yes'})
    figYAxis([],'PC #',[],1:npcs,pcLabel_choice)
    figAxForm([],0)
    title(avName{iav})
    
    figure(weightsExPCs_fig)
    subplot(2,2,iav)
    imagesc(pc_wStim(sortInd_stim))
    clim([-1 3])
    colorbar
    figXAxis([],'',[])
    figYAxis([],'Stimulus Weight',[])
    figAxForm([],0)
    title(avName{iav})
    subplot(2,2,iav+2)
    imagesc(pc_wChoice(sortInd_stim))
    clim([-1 3])
    colorbar
    figXAxis([],'',[])
    figYAxis([],'Choice Weight',[])
    figAxForm([],0)
    title(avName{iav})
end
figure(stimExPCs_fig)
print([fnout 'exExpt_pcHM_stim'],'-dpdf')
figure(choiceExPCs_fig)
print([fnout 'exExpt_pcHM_choice'],'-dpdf')
figure(weightsExPCs_fig)
print([fnout 'exExpt_pcHM_weights'],'-dpdf')

% %% catch trial performance
% catchExptInd_attn = true(1,nexp_attn);
% for imod = 1:2
%     mdl_attn(visualTrials).mdl(imod).pctCorrect_inv = nan(1,nexp_attn);
%     mdl_attn(visualTrials).mdl(imod).pctCorrect_valMatched = nan(1,nexp_attn);
%     mdl_attn(visualTrials).mdl(imod).pctCorrect_invAudTest = nan(1,nexp_attn);
%     mdl_attn(visualTrials).mdl(imod).pctCorrect_valMatchedAudTest = nan(1,nexp_attn);
%     mdl_attn(visualTrials).mdl(imod).catchTrialResult_inv = cell(1,nexp_attn);
%     mdl_attn(visualTrials).mdl(imod).catchTrialResult_valMatched = cell(1,nexp_attn);
%     mdl_attn(visualTrials).mdl(imod).catchTrialResult_invAudTest = cell(1,nexp_attn);
%     mdl_attn(visualTrials).mdl(imod).catchTrialResult_valMatchedAudTest = cell(1,nexp_attn);
%     
%     for iexp = 1:nexp_attn
%         if isempty(dc_attn(iexp).av(visualTrials).invalidPctCorrectDetect)
%             catchExptInd_attn(iexp) = false;
%             continue
%         end
%         if imod == 1
%             modType = 'Target';
%         else
%             modType = 'Detect';
%         end
%         mdl_attn(visualTrials).mdl(imod).pctCorrect_inv(iexp) = ...
%             eval(['dc_attn(iexp).av(visualTrials).invalidPctCorrect' modType]);
%         mdl_attn(visualTrials).mdl(imod).pctCorrect_valMatched(iexp) = ...
%             eval(['dc_attn(iexp).av(visualTrials).validMatchedPctCorrect'...
%             modType '_holdout']);            
%         mdl_attn(visualTrials).mdl(imod).pctCorrect_invAudTest(iexp) = ...
%             eval(['dc_attn(iexp).av(visualTrials).invalidPctCorrect'...
%             modType '_testAudModel']);  
%         mdl_attn(visualTrials).mdl(imod).pctCorrect_valMatchedAudTest(iexp) = ...
%             eval(['dc_attn(iexp).av(visualTrials).validMatchedPctCorrect'...
%             modType '_testAudModel']);  
%         
%         mdl_attn(visualTrials).mdl(imod).catchTrialResult_inv{iexp} = ...
%             eval(['dc_attn(iexp).av(visualTrials).invalidCorrectTrials' modType])';
%         mdl_attn(visualTrials).mdl(imod).catchTrialResult_valMatched{iexp} = ...
%             eval(['dc_attn(iexp).av(visualTrials).validMatchedCorrectTrials' modType])';            
%         mdl_attn(visualTrials).mdl(imod).catchTrialResult_invAudTest{iexp} = ...
%             eval(['dc_attn(iexp).av(visualTrials).invalidCorrectTrials'...
%             modType '_testAud'])';  
%         mdl_attn(visualTrials).mdl(imod).catchTrialResult_valMatchedAudTest{iexp} = ...
%             eval(['dc_attn(iexp).av(visualTrials).validMatchedCorrectTrials'...
%             modType '_testAud'])';  
%     end
% end
% 
% catchExptInd_noAttn = true(1,nexp_noAttn);
% for imod = 1:2
%     mdl_noAttn(visualTrials).mdl(imod).pctCorrect_inv = nan(1,nexp_noAttn);
%     mdl_noAttn(visualTrials).mdl(imod).pctCorrect_valMatched = nan(1,nexp_noAttn);
%     mdl_noAttn(visualTrials).mdl(imod).pctCorrect_invAudTest = nan(1,nexp_noAttn);
%     mdl_noAttn(visualTrials).mdl(imod).pctCorrect_valMatchedAudTest = nan(1,nexp_noAttn);
%     mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_inv = cell(1,nexp_noAttn);
%     mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_valMatched = cell(1,nexp_noAttn);
%     mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_invAudTest = cell(1,nexp_noAttn);
%     mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_valMatchedAudTest = cell(1,nexp_noAttn);
%     
%     for iexp = 1:nexp_noAttn
%         if isempty(dc_noAttn(iexp).av(visualTrials).invalidPctCorrectDetect)
%             catchExptInd_noAttn(iexp) = false;
%             continue
%         end
%         if imod == 1
%             modType = 'Target';
%         else
%             modType = 'Detect';
%         end
%         mdl_noAttn(visualTrials).mdl(imod).pctCorrect_inv(iexp) = ...
%             eval(['dc_noAttn(iexp).av(visualTrials).invalidPctCorrect' modType]);
%         mdl_noAttn(visualTrials).mdl(imod).pctCorrect_valMatched(iexp) = ...
%             eval(['dc_noAttn(iexp).av(visualTrials).validMatchedPctCorrect'...
%             modType '_holdout']);            
%         mdl_noAttn(visualTrials).mdl(imod).pctCorrect_invAudTest(iexp) = ...
%             eval(['dc_noAttn(iexp).av(visualTrials).invalidPctCorrect'...
%             modType '_testAudModel']);  
%         mdl_noAttn(visualTrials).mdl(imod).pctCorrect_valMatchedAudTest(iexp) = ...
%             eval(['dc_noAttn(iexp).av(visualTrials).validMatchedPctCorrect'...
%             modType '_testAudModel']);  
%         
%         mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_inv{iexp} = ...
%             eval(['dc_noAttn(iexp).av(visualTrials).invalidCorrectTrials' modType])';
%         mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_valMatched{iexp} = ...
%             eval(['dc_noAttn(iexp).av(visualTrials).validMatchedCorrectTrials' modType])';            
%         mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_invAudTest{iexp} = ...
%             eval(['dc_noAttn(iexp).av(visualTrials).invalidCorrectTrials'...
%             modType '_testAud'])';  
%         mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_valMatchedAudTest{iexp} = ...
%             eval(['dc_noAttn(iexp).av(visualTrials).validMatchedCorrectTrials'...
%             modType '_testAud'])';  
%     end
% end
% 
% figure
% for imod = 1:2
%     ind = mdl_attn(visualTrials).mdl(imod).pctCorrect_valMatched > pctCorrThresh ...
%         & catchExptInd_attn;
%     y = cat(1,mdl_attn(visualTrials).mdl(imod).pctCorrect_valMatched(ind),...
%         mdl_attn(visualTrials).mdl(imod).pctCorrect_inv(ind),...
%         mdl_attn(visualTrials).mdl(imod).pctCorrect_invAudTest(ind));
%     
%     subplot(2,2,imod)
%     plot(1:3,y','k.-','MarkerSize',10)
%     hold on
%     errorbar(1:3,mean(y',1),ste(y',1),'.','MarkerSize',20)
%     figXAxis([],'Train-Test',[0 4],1:3,{'Vis-ValVis','Vis-InvVis','Aud-InvVis'})
%     figYAxis([],'Pct Correct',[0 1])
%     figAxForm
%     title(sprintf('Attn, %s Model',modelName{imod}))
%     hline(0.5,'k:')
%     
%     ind = mdl_noAttn(visualTrials).mdl(imod).pctCorrect_valMatched > pctCorrThresh ...
%         & catchExptInd_noAttn;
%     y = cat(1,mdl_noAttn(visualTrials).mdl(imod).pctCorrect_valMatched(ind),...
%         mdl_noAttn(visualTrials).mdl(imod).pctCorrect_inv(ind),...
%         mdl_noAttn(visualTrials).mdl(imod).pctCorrect_invAudTest(ind));
%     
%     subplot(2,2,imod+2)
%     plot(1:3,y','k.-','MarkerSize',10)
%     hold on
%     errorbar(1:3,mean(y',1),ste(y',1),'.','MarkerSize',20)
%     figXAxis([],'Train-Test',[0 4],1:3,{'Vis-ValVis','Vis-InvVis','Aud-InvVis'})
%     figYAxis([],'Pct Correct',[0 1])
%     figAxForm
%     title(sprintf('No Attn, %s Model',modelName{imod}))
%     hline(0.5,'k:')
% end
% print([fnout 'invalidTrialPerf'],'-dpdf','-fillpage')
% 
% figure
% for imod = 1:2
%     subplot(2,2,imod)
%     ind = cellfun(@(x) ~isempty(x),mdl_attn(visualTrials).mdl(imod).catchTrialResult_inv) & ...
%         mdl_attn(visualTrials).mdl(imod).pctCorr_all > pctCorrThresh;
%     invTrialResult = cell2mat(mdl_attn(visualTrials).mdl(imod).catchTrialResult_inv(ind));
%     [invHR,invCI] = binofit(sum(invTrialResult),length(invTrialResult));
%     valTrialResult = cell2mat(mdl_attn(visualTrials).mdl(imod).catchTrialResult_valMatched(ind));
%     [valHR,valCI] = binofit(sum(valTrialResult),length(valTrialResult));
%     invTrialResult_audModel = cell2mat(mdl_attn(visualTrials).mdl(imod).catchTrialResult_invAudTest(ind));
%     [invHR_audModel,invCI_audModel] = binofit(sum(invTrialResult_audModel),length(invTrialResult_audModel));
%     valTrialResult_audModel = cell2mat(mdl_attn(visualTrials).mdl(imod).catchTrialResult_valMatchedAudTest(ind));
%     [valHR_audModel,valCI_audModel] = binofit(sum(valTrialResult_audModel),length(valTrialResult_audModel));
%     
%     y = [valHR,invHR,valHR_audModel,invHR_audModel];
%     y_lerr = y - [valCI(1),invCI(1),valCI_audModel(1),invCI_audModel(1)];
%     y_uerr = [valCI(2),invCI(2),valCI_audModel(2),invCI_audModel(2)]-y;
%     errorbar(1:4,y, y_lerr,y_uerr,'.','MarkerSize',20)
%     figXAxis([],'Train-Test',[0 5],1:4,{'Vis-Val','Vis-Inv','Aud-Val','Aud-Inv'})
%     figYAxis([],'Pct Correct',[0 1])
%     figAxForm
%     hline(0.5,'k:')
%     title(sprintf('Attn; %s Model, Vis Trials',modelName{imod}))
%     
%     subplot(2,2,imod+2)
%     ind = cellfun(@(x) ~isempty(x),mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_inv) & ...
%         mdl_noAttn(visualTrials).mdl(imod).pctCorr_all > pctCorrThresh;
%     invTrialResult = cell2mat(mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_inv(ind));
%     [invHR,invCI] = binofit(sum(invTrialResult),length(invTrialResult));
%     valTrialResult = cell2mat(mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_valMatched(ind));
%     [valHR,valCI] = binofit(sum(valTrialResult),length(valTrialResult));
%     invTrialResult_audModel = cell2mat(mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_invAudTest(ind));
%     [invHR_audModel,invCI_audModel] = binofit(sum(invTrialResult_audModel),length(invTrialResult_audModel));
%     valTrialResult_audModel = cell2mat(mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_valMatchedAudTest(ind));
%     [valHR_audModel,valCI_audModel] = binofit(sum(valTrialResult_audModel),length(valTrialResult_audModel));
%     
%     y = [valHR,invHR,valHR_audModel,invHR_audModel];
%     y_lerr = y - [valCI(1),invCI(1),valCI_audModel(1),invCI_audModel(1)];
%     y_uerr = [valCI(2),invCI(2),valCI_audModel(2),invCI_audModel(2)]-y;
%     errorbar(1:4,y, y_lerr,y_uerr,'.','MarkerSize',20)
%     figXAxis([],'Train-Test',[0 5],1:4,{'Vis-Val','Vis-Inv','Aud-Val','Aud-Inv'})
%     figYAxis([],'Pct Correct',[0 1])
%     figAxForm
%     hline(0.5,'k:')
%     title(sprintf('No Attn; %s Model, Vis Trials',modelName{imod}))
% end
% print([fnout 'invalidTrialPerf_combExpt'],'-dpdf','-fillpage')
% 
% figure
% ind = goodPerfInd_attn & catchExptInd_attn;
% x = mdl_attn(visualTrials).mdl(stimModInd).pctCorrect_valMatched(ind) -...
%     mdl_attn(visualTrials).mdl(stimModInd).pctCorrect_inv(ind);
% y = mdl_attn(visualTrials).mdl(choiceModInd).pctCorrect_valMatched(ind) -...
%     mdl_attn(visualTrials).mdl(choiceModInd).pctCorrect_inv(ind);
% subplot(2,2,1)
% plot(x,y,'.','MarkerSize',20)
% hold on
% plot([-0.5 0.7],[-0.5 0.7],'k:')
% figXAxis([],'Stim Model, Val-Inv',[-0.2 0.7])
% figYAxis([],'Choice Model, Val-Inv',[-0.2 0.7])
% figAxForm
% vline(0,'k:')
% hline(0,'k:')
% title('No Attn')
% 
% valInvDiff_stim_eaMouse = nan(1,length(mice_attn));
% valInvDiff_choice_eaMouse = nan(1,length(mice_attn));
% for im = 1:length(mice_attn)
%     ind = goodPerfInd_attn & catchExptInd_attn & ...
%         strcmp(mouseName_attn,mice_attn(im));
%     x = mdl_attn(visualTrials).mdl(stimModInd).pctCorrect_valMatched(ind) -...
%         mdl_attn(visualTrials).mdl(stimModInd).pctCorrect_inv(ind);
%     y = mdl_attn(visualTrials).mdl(choiceModInd).pctCorrect_valMatched(ind) -...
%         mdl_attn(visualTrials).mdl(choiceModInd).pctCorrect_inv(ind);
%     valInvDiff_stim_eaMouse(im) = mean(x);
%     valInvDiff_choice_eaMouse(im) = mean(y);
% end
% 
% subplot(2,2,3)
% plot(hrDiff_attn,valInvDiff_stim_eaMouse,'.','MarkerSize',20)
% hold on
% figXAxis([],'HR Val-Inv',[-0 0.5])
% figYAxis([],'Stim Model, Val-Inv',[-0.2 0.7])
% figAxForm
% vline(0,'k:')
% hline(0,'k:')
% subplot(2,2,4)
% plot(hrDiff_attn,valInvDiff_choice_eaMouse,'.','MarkerSize',20)
% hold on
% figXAxis([],'HR Val-Inv',[-0 0.5])
% figYAxis([],'Stim Model, Val-Inv',[-0.2 0.7])
% figAxForm
% vline(0,'k:')
% hline(0,'k:')
% 
% ind = goodPerfInd_noAttn & catchExptInd_noAttn;
% x = mdl_noAttn(visualTrials).mdl(stimModInd).pctCorrect_valMatched(ind) -...
%     mdl_noAttn(visualTrials).mdl(stimModInd).pctCorrect_inv(ind);
% y = mdl_noAttn(visualTrials).mdl(choiceModInd).pctCorrect_valMatched(ind) -...
%     mdl_noAttn(visualTrials).mdl(choiceModInd).pctCorrect_inv(ind);
% subplot(2,2,2)
% plot(x,y,'.','MarkerSize',20)
% hold on
% plot([-0.5 0.7],[-0.5 0.7],'k:')
% figXAxis([],'Stim Model, Val-Inv',[-0.2 0.7])
% figYAxis([],'Choice Model, Val-Inv',[-0.2 0.7])
% figAxForm
% vline(0,'k:')
% hline(0,'k:')
% title('No Attn')
% 
% valInvDiff_stim_eaMouse = nan(1,length(mice_noAttn));
% valInvDiff_choice_eaMouse = nan(1,length(mice_noAttn));
% for im = 1:length(mice_noAttn)
%     ind = goodPerfInd_noAttn & catchExptInd_noAttn & ...
%         strcmp(mouseName_noAttn,mice_noAttn(im));
%     x = mdl_noAttn(visualTrials).mdl(stimModInd).pctCorrect_valMatched(ind) -...
%         mdl_noAttn(visualTrials).mdl(stimModInd).pctCorrect_inv(ind);
%     y = mdl_noAttn(visualTrials).mdl(choiceModInd).pctCorrect_valMatched(ind) -...
%         mdl_noAttn(visualTrials).mdl(choiceModInd).pctCorrect_inv(ind);
%     valInvDiff_stim_eaMouse(im) = mean(x);
%     valInvDiff_choice_eaMouse(im) = mean(y);
% end
% 
% subplot(2,2,3)
% plot(hrDiff_noAttn,valInvDiff_stim_eaMouse,'.','MarkerSize',20)
% hold on
% figXAxis([],'HR Val-Inv',[-0 0.5])
% figYAxis([],'Stim Model, Val-Inv',[-0.2 0.7])
% figAxForm
% vline(0,'k:')
% hline(0,'k:')
% legend({'Attn','No Attn'})
% subplot(2,2,4)
% plot(hrDiff_noAttn,valInvDiff_choice_eaMouse,'.','MarkerSize',20)
% hold on
% figXAxis([],'HR Val-Inv',[-0 0.5])
% figYAxis([],'Stim Model, Val-Inv',[-0.2 0.7])
% figAxForm
% vline(0,'k:')
% hline(0,'k:')
% legend({'Attn','No Attn'})
% print([fnout 'invalidTrialPerfxBehavior'],'-dpdf','-fillpage')

%% distractor only choice model performance
for iav = 1:2
    mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly_otherAV = nan(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_distOnly_holdout;
        mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly_otherAV(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_distOnly_otherAV;
    end
    
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly_otherAV = nan(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_distOnly_holdout;
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly_otherAV(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_distOnly_otherAV;
    end
end

figure
suptitle('Distractor Only Trained')
for iav = 1:2
    if iav == 1
        otherAV = 2;
    else
        otherAV = 1;
    end
    subplot(2,2,iav)
    y = cat(1,mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly,...
        mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly_otherAV);
    ind = y(1,:)>pctCorrThresh;
    plot(1:2,y(:,ind),'k-')
    hold on
    errorbar(1:2,mean(y(:,ind),2),ste(y(:,ind),2),'.','MarkerSize',20);
    [~,p] = ttest(y(1,:),y(2,:));
    [~,p_chance] = ttest(y');
    figXAxis([],'Test',[0 3],1:2,{avName{iav},avName{otherAV}})
    figYAxis([],'Pct Correct',[0 1])
    figAxForm
    hline(0.5,'k-')
    hline(pctCorrThresh,'k:')
    title(sprintf('Attn, %s Choice Model, p=%s',avName{iav},num2str(round(p,2,'significant'))))    
    text([0:2],ones(1,3)*.2,{'pCh',num2str(round(p_chance(1),2,'significant')),...
        num2str(round(p_chance(2),2,'significant'))})
    
    subplot(2,2,iav+2)
    y = cat(1,mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly,...
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly_otherAV);
    ind = y(1,:)>pctCorrThresh;
    plot(1:2,y(:,ind),'k-')
    hold on
    errorbar(1:2,mean(y(:,ind),2),ste(y(:,ind),2),'.','MarkerSize',20);
    [~,p] = ttest(y(1,:),y(2,:));
    [~,p_chance] = ttest(y');
    figXAxis([],'Test',[0 3],1:2,{avName{iav},avName{otherAV}})
    figYAxis([],'Pct Correct',[0 1])
    figAxForm
    hline(0.5,'k-')
    hline(pctCorrThresh,'k:')
    title(sprintf('No Attn, %s Choice Model, p=%s',avName{iav},num2str(round(p,2,'significant'))))    
    text([0:2],ones(1,3)*.2,{'pCh',num2str(round(p_chance(1),2,'significant')),...
        num2str(round(p_chance(2),2,'significant'))})
end
print([fnout 'pctCorrectOtherAVModel_distOnly'],'-dpdf','-fillpage')

figure
suptitle('Distractor Only Trained')
for iav = 1:2
    subplot(1,2,iav)
    y1 = mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly(goodPerfInd_attn);
    plot(ones(1,sum(goodPerfInd_attn)),y1,'k.','MarkerSize',10)
    hold on
    errorbar(1,mean(y1),ste(y1,2),'.','MarkerSize',20)
    y2 = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly(goodPerfInd_noAttn);
    plot(ones(1,sum(goodPerfInd_noAttn)).*2,y2,'k.','MarkerSize',10)
    hold on
    errorbar(2,mean(y2),ste(y2,2),'.','MarkerSize',20)
    [~,p] = ttest2(y1,y2);
    figXAxis([],'Test',[0 3],1:2,{'Attn','No Attn'})
    figYAxis([],'Pct Correct',[0 1])
    figAxForm
    hline(0.5,'k-')
    hline(pctCorrThresh,'k:')
    title(sprintf('%s Choice Model',avName{iav}))     
end
print([fnout 'pctCorrect_distOnly'],'-dpdf','-fillpage')    

%% distractor only model weights

for iav = 1:2
    mdl_attn(iav).mdl(choiceModInd).weights_distOnly = cell(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(choiceModInd).weights_distOnly{iexp} = ...
            dc_attn(iexp).av(iav).weightDetect_distOnly;
    end
    mdl_attn(iav).mdl(choiceModInd).normWeights_distOnly = ...
        cellfun(@(x) x./norm(x), ...
        mdl_attn(iav).mdl(choiceModInd).weights_distOnly,'unif',0);
    
    mdl_noAttn(iav).mdl(choiceModInd).weights_distOnly = cell(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(choiceModInd).weights_distOnly{iexp} = ...
            dc_noAttn(iexp).av(iav).weightDetect_distOnly;
    end
    mdl_noAttn(iav).mdl(choiceModInd).normWeights_distOnly = ...
        cellfun(@(x) x./norm(x), ...
        mdl_noAttn(iav).mdl(choiceModInd).weights_distOnly,'unif',0);
end

weightTooHigh_distOnly_attn = cellfun(@(x) max(abs(x)),...
    mdl_attn(visualTrials).mdl(choiceModInd).weights_distOnly)> 30 |...
    cellfun(@(x) max(abs(x)),...
    mdl_attn(auditoryTrials).mdl(choiceModInd).weights_distOnly) > 30;
weightTooHigh_distOnly_noAttn = cellfun(@(x) max(abs(x)),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).weights_distOnly)> 30 |...
    cellfun(@(x) max(abs(x)),...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_distOnly) > 30;
goodPerfInd_av_distOnly_attn = (mdl_attn(visualTrials).mdl(choiceModInd).pctCorr_distOnly > pctCorrThresh |...
    mdl_attn(auditoryTrials).mdl(choiceModInd).pctCorr_distOnly > pctCorrThresh)...
    & ~weightTooHigh_distOnly_attn;
goodPerfInd_av_distOnly_noAttn = (mdl_noAttn(visualTrials).mdl(choiceModInd).pctCorr_distOnly > pctCorrThresh |...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).pctCorr_distOnly > pctCorrThresh)...
    & ~weightTooHigh_distOnly_noAttn;

for iav = 1:2
    rng(0)
    mdl_attn(iav).mdl(choiceModInd).shuffWeights_distOnly = cellfun(...
        @(x) x(randperm(length(x))),mdl_attn(iav).mdl(choiceModInd).weights_distOnly,'unif',0);      
    mdl_noAttn(iav).mdl(choiceModInd).shuffWeights_distOnly = cellfun(...
        @(x) x(randperm(length(x))),mdl_noAttn(iav).mdl(choiceModInd).weights_distOnly,'unif',0);  
    if iav == 1
        for ishuff = 1:nShuff
%                 
            % shuffle normalized weights
            shuffWeights_vis = cellfun(@(x) x(randperm(length(x))),...
                mdl_attn(visualTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_attn),'unif',0)';
            shuffWeights_aud = cellfun(@(x) x(randperm(length(x))),...
                mdl_attn(auditoryTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_attn),'unif',0)';
            if ishuff == 1
                mdl_attn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_distOnly = ...
                    nan(sum(cellfun(@length,shuffWeights_vis)),nShuff);
            end
            mdl_attn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_distOnly(:,ishuff) = ...
                cell2mat(cellfun(@(x,y) ...
                rad2deg(atan(x./y)),shuffWeights_aud,...
                shuffWeights_vis,'unif',0));

            shuffWeights_vis = cellfun(@(x) x(randperm(length(x))),...
                mdl_noAttn(visualTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_noAttn),'unif',0)';
            shuffWeights_aud = cellfun(@(x) x(randperm(length(x))),...
                mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_noAttn),'unif',0)';
            if ishuff == 1
                mdl_noAttn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_distOnly = ...
                    nan(sum(cellfun(@length,shuffWeights_vis)),nShuff);
            end
            mdl_noAttn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_distOnly(:,ishuff) = ...
                cell2mat(cellfun(@(x,y) ...
                rad2deg(atan(x./y)),shuffWeights_aud,...
                shuffWeights_vis,'unif',0));
        end
    end
end


figure
y1 = cell2mat(mdl_attn(visualTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_attn)');
y2 = cell2mat(mdl_attn(auditoryTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_attn)');
n_attn = length(y1);
%     [r_attn,p_attn,rCIl_attn,rCIu_attn] = corrcoef(y1,y2);
[r_attn,p_attn] = corr(y1,y2);
y1 = cell2mat(mdl_noAttn(visualTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_noAttn)');
y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_noAttn)');
n_noAttn = length(y1);
[r_noAttn,p_noAttn] = corr(y1,y2);
bar(1,r_attn)
text(1,0.1,sigfigString(p_attn))
hold on
bar(2,r_noAttn)
text(2,0.1,sigfigString(p_noAttn))
figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
figYAxis([],'Vis-Aud Corr',[-0.2 1])
figAxForm
title(sprintf('%s Model, Distractors Only, p=%s',modelName{choiceModInd},...
    sigfigString(compare_correlation_coefficients(...
    r_attn,n_attn,r_noAttn,n_noAttn))))
print([fnout 'avWeightCorr_distOnly'],'-dpdf')

figure
suptitle('Distractor Only Model')
for iav = 1:2
    if iav == 1
        iplot = 0;
    else
        iplot = 2;
    end
    subplot(1,4,iav)
    y = cell2mat(mdl_attn(iav).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_attn)');
    histogram(y,weightEdges,'Normalization','probability')
    figXAxis([],'Weight',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Fraction of Cells',[0 .7])
    figAxForm
    title({'Attn Mice';sprintf('%s %s Model',...
        avName{iav},modelName{choiceModInd})})
    vline(0,'k:')
    subplot(1,4,iav+2)
    y = cell2mat(mdl_noAttn(iav).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_noAttn)');
    histogram(y,weightEdges,'Normalization','probability')
    figXAxis([],'Weight',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Fraction of Cells',[0 .7])
    figAxForm
    title({'No Attn Mice';sprintf('%s %s Model',...
        avName{iav},modelName{choiceModInd})})
    vline(0,'k:')
end
print([fnout 'weightHist_distOnly'],'-dpdf','-fillpage')

figure
suptitle('Distractor Only Model')
for iav = 1:2
    subplot(1,2,iav)
    y1 = cell2mat(mdl_attn(iav).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_attn)');
    y2 = cell2mat(mdl_noAttn(iav).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_noAttn)');
    cdfplot(y1)
    hold on
    cdfplot(y2)
    figXAxis([],'Weight',[-1 1])
    figYAxis([],'Fraction of Cells',[0 1])
    figAxForm
    title({'Attn Mice';sprintf('%s %s Model',...
        avName{iav},modelName{choiceModInd})})
    vline(0,'k:')
    legend({'Attn','No Attn'},'location','southeast')
end
print([fnout 'weightCDF_distOnly'],'-dpdf','-fillpage')

mdlWeightsSub_fig = figure;
mdlWeights_fig = figure;
suptitle({'Distractor Only Model,Weights';...
    sprintf('pct corr must be > %s in vis & aud stim models',...
    num2str(pctCorrThresh))})

figure(mdlWeights_fig)
subplot(3,2,1)
y1 = cell2mat(mdl_attn(visualTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_attn)');
y2 = cell2mat(mdl_attn(auditoryTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_attn)');
plot(y1,y2,'k.','MarkerSize',10)
hold on
mdl_fitline = fitlm(y1,y2);
rsq = mdl_fitline.Rsquared.Ordinary;
yfit = predict(mdl_fitline,y1);
plot(y1,yfit,'-')
figXAxis([],'Visual Weight',[weightEdges(1) weightEdges(end)])
figYAxis([],'Auditory Weight',[weightEdges(1) weightEdges(end)])
figAxForm
hline(0,'k:')
vline(0,'k:')
r = corr(y1,y2);
title({'Attn Mice';sprintf('%s Model, r=%s',...
    modelName{choiceModInd},num2str(round(r,2,'significant')))})

subplot(3,2,3)
weightAVangle = rad2deg(atan(y2./y1));
histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
s = skewness(weightAVangle);
title({'Attn Mice';sprintf('%s Model, skew=%s,ks=%s',...
    modelName{choiceModInd},num2str(round(s,2,'significant')))})

subplot(3,2,5)
%     y = mean(mdl_attn(visualTrials).mdl(choiceModInd).normWeightAVangle_shuff,2);
%     histogram(y,oriEdges,'Normalization','probability')
avAngleShuff = mdl_attn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_distOnly;
h = nan(length(oriEdges_noAbs)-1,nShuff);
for ishuff = 1:nShuff
    h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges_noAbs,'Normalization','Probability');
end    
h_sort = sort(h,2);
y = mean(h_sort,2);
yerrl = y-h_sort(:,round(nShuff*0.025));
yerru = h_sort(:,round(nShuff*0.975))-y;
shadedErrorBar_chooseColor(oriEdges_noAbs(1:end-1),y,[yerrl';yerru'],[0 0 0]);
hold on
plot(oriEdges_noAbs(1:end-1),histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability'),'.-')
controlSub_attn = histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability') - y';
figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
s = skewness(y);
[~,p] = kstest2(h_sort(:),weightAVangle);
title({'Attn Mice-Shuff';sprintf('%s Model, skew=%s,p=%s',...
    modelName{choiceModInd},num2str(round(s,2,'significant')),...
    sigfigString(p))})

figure(mdlWeightsSub_fig)
plot(oriEdges_noAbs(1:end-1),controlSub_attn,'.-')
%     shadedErrorBar_chooseColor(oriEdges(1:end-1),controlSub_attn,[yerrl';yerru'],[0 0 0]);

figure(mdlWeights_fig)
subplot(3,2,2)
y1 = cell2mat(mdl_noAttn(visualTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_noAttn)');
y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_noAttn)');
plot(y1,y2,'k.','MarkerSize',10)
hold on
mdl_fitline = fitlm(y1,y2);
rsq = mdl_fitline.Rsquared.Ordinary;
yfit = predict(mdl_fitline,y1);
plot(y1,yfit,'-')
figXAxis([],'Visual Weight',[weightEdges(1) weightEdges(end)])
figYAxis([],'Auditory Weight',[weightEdges(1) weightEdges(end)])
figAxForm
hline(0,'k:')
vline(0,'k:')
r = corr(y1,y2);
title({'No Attn Mice';sprintf('%s Model, r=%s',...
    modelName{choiceModInd},num2str(round(r,2,'significant')))})

subplot(3,2,4)
weightAVangle = rad2deg(atan(y2./y1));
histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
s = skewness(weightAVangle);
title({'No Attn Mice';sprintf('%s Model, skew=%s',...
    modelName{choiceModInd},num2str(round(s,2,'significant')))})

subplot(3,2,6)
%     y = mean(mdl_noAttn(visualTrials).mdl(choiceModInd).normWeightAVangle_shuff,2);
%     histogram(y,oriEdges,'Normalization','probability')
avAngleShuff = mdl_noAttn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_distOnly;
h = nan(length(oriEdges_noAbs)-1,nShuff);
for ishuff = 1:nShuff
    h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges_noAbs,'Normalization','Probability');
end    
h_sort = sort(h,2);
y = mean(h_sort,2);
yerrl = y-h_sort(:,round(nShuff*0.025));
yerru = h_sort(:,round(nShuff*0.975))-y;
shadedErrorBar_chooseColor(oriEdges_noAbs(1:end-1),y,[yerrl';yerru'],[0 0 0]);
hold on
plot(oriEdges_noAbs(1:end-1),histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability'),'.-')
controlSub_noAttn = histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability') - y';
figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
s = skewness(y);
title({'No Attn Mice-Shuff';sprintf('%s Model, skew=%s',...
    modelName{choiceModInd},num2str(round(s,2,'significant')))})

figure(mdlWeightsSub_fig)
hold on
plot(oriEdges_noAbs(1:end-1),controlSub_noAttn,'.-')
figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
figYAxis([],'Control Subtr. Fraction of Cells',[-0.15 0.15])
figAxForm
hline(0,'k:')
title(sprintf('%s Model, Distractors Only',modelName{choiceModInd}))

figure(mdlWeights_fig)
print([fnout 'modelWeights_distOnly'],'-dpdf','-fillpage')
figure(mdlWeightsSub_fig)
print([fnout 'modelWeights_ciSub_distOnly'],'-dpdf')

%% compare PC weights
weight_lim_pcs = [-5 6];
for iav = 1:2
    mdl_attn(iav).mdl(choiceModInd).weights_pcs = cell(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(choiceModInd).weights_pcs{iexp} = ...
            dc_attn(iexp).av(iav).weightDetect_pcs;
    end
    
    mdl_noAttn(iav).mdl(choiceModInd).weights_pcs = cell(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(choiceModInd).weights_pcs{iexp} = ...
            dc_noAttn(iexp).av(iav).weightDetect_pcs;
    end
end
for iav = 1:2
    mdl_attn(iav).mdl(choiceModInd).weights_pcs_distOnly = cell(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(choiceModInd).weights_pcs_distOnly{iexp} = ...
            dc_attn(iexp).av(iav).weightDetect_pcs_distOnly;
    end
    
    mdl_noAttn(iav).mdl(choiceModInd).weights_pcs_distOnly = cell(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(choiceModInd).weights_pcs_distOnly{iexp} = ...
            dc_noAttn(iexp).av(iav).weightDetect_pcs_distOnly;
    end
end

% *** which experiments: goodPerfInd_av_distOnly_attn   goodPerfInd_attn
nboot = 1000;

figure
suptitle('PC weights')
subplot 231; hold on
y1 = cell2mat(mdl_attn(visualTrials).mdl(choiceModInd).weights_pcs(goodPerfInd_attn)');
y2 = cell2mat(mdl_attn(auditoryTrials).mdl(choiceModInd).weights_pcs(goodPerfInd_attn)');
plot(y1,y2,'k.','MarkerSize',10)
hold on
mdl_fitline = fitlm(y1,y2);
rsq = mdl_fitline.Rsquared.Ordinary;
yfit = predict(mdl_fitline,y1);
plot(y1,yfit,'-')
figXAxis([],'Visual Weight',weight_lim_pcs)
figYAxis([],'Auditory Weight',weight_lim_pcs)
figAxForm
hline(0,'k:')
vline(0,'k:')
r = corr(y1,y2);
title({'Attn Mice';sprintf('%s Model, r=%s',...
    modelName{choiceModInd},num2str(round(r,2,'significant')))})
r_attn = nan(1,nboot);
for iboot = 1:nboot
    n = length(y1);
    ind = randsample(1:n,n,1);
    y1_samp = y1(ind);
    y2_samp = y2(ind);
    r_attn(iboot) = corr(y1_samp,y2_samp);
end
subplot 232; hold on
y1 = cell2mat(mdl_noAttn(visualTrials).mdl(choiceModInd).weights_pcs(goodPerfInd_noAttn)');
y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_pcs(goodPerfInd_noAttn)');
plot(y1,y2,'k.','MarkerSize',10)
hold on
mdl_fitline = fitlm(y1,y2);
rsq = mdl_fitline.Rsquared.Ordinary;
yfit = predict(mdl_fitline,y1);
plot(y1,yfit,'-')
figXAxis([],'Visual Weight',weight_lim_pcs)
figYAxis([],'Auditory Weight',weight_lim_pcs)
figAxForm
hline(0,'k:')
vline(0,'k:')
r = corr(y1,y2);
title({'No Attn Mice';sprintf('%s Model, r=%s',...
    modelName{choiceModInd},num2str(round(r,2,'significant')))})
r_noAttn = nan(1,nboot);
for iboot = 1:nboot
    n = length(y1);
    ind = randsample(1:n,n,1);
    y1_samp = y1(ind);
    y2_samp = y2(ind);
    r_noAttn(iboot) = corr(y1_samp,y2_samp);
end
subplot 233; hold on
bar(1:2,[mean(r_attn),mean(r_noAttn)])
[ylerr,yuerr] = ciFromBoot(r_attn',95);
errorbar(1,mean(r_attn),ylerr,yuerr,'.','MarkerSize',20)
[ylerr,yuerr] = ciFromBoot(r_noAttn',95);
errorbar(2,mean(r_noAttn),ylerr,yuerr,'.','MarkerSize',20)
figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
xtickangle(-45)
figYAxis([],'Vis:Aud Weight Corr',[0 1])
figAxForm
title(sprintf('%s Model',mdl_noAttn(iav).mdl(imod).name))

subplot 234; hold on
y1 = cell2mat(mdl_attn(visualTrials).mdl(choiceModInd).weights_pcs_distOnly(goodPerfInd_av_distOnly_attn)');
y2 = cell2mat(mdl_attn(auditoryTrials).mdl(choiceModInd).weights_pcs_distOnly(goodPerfInd_av_distOnly_attn)');
plot(y1,y2,'k.','MarkerSize',10)
hold on
mdl_fitline = fitlm(y1,y2);
rsq = mdl_fitline.Rsquared.Ordinary;
yfit = predict(mdl_fitline,y1);
plot(y1,yfit,'-')
figXAxis([],'Visual Weight',weight_lim_pcs)
figYAxis([],'Auditory Weight',weight_lim_pcs)
figAxForm
hline(0,'k:')
vline(0,'k:')
r = corr(y1,y2);
title({'Attn Mice';sprintf('Dist Only %s Model, r=%s',...
    modelName{choiceModInd},num2str(round(r,2,'significant')))})
r_attn = nan(1,nboot);
for iboot = 1:nboot
    n = length(y1);
    ind = randsample(1:n,n,1);
    y1_samp = y1(ind);
    y2_samp = y2(ind);
    r_attn(iboot) = corr(y1_samp,y2_samp);
end
subplot 235; hold on
y1 = cell2mat(mdl_noAttn(visualTrials).mdl(choiceModInd).weights_pcs_distOnly(goodPerfInd_av_distOnly_noAttn)');
y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_pcs_distOnly(goodPerfInd_av_distOnly_noAttn)');
plot(y1,y2,'k.','MarkerSize',10)
hold on
mdl_fitline = fitlm(y1,y2);
rsq = mdl_fitline.Rsquared.Ordinary;
yfit = predict(mdl_fitline,y1);
plot(y1,yfit,'-')
figXAxis([],'Visual Weight',weight_lim_pcs)
figYAxis([],'Auditory Weight',weight_lim_pcs)
figAxForm
hline(0,'k:')
vline(0,'k:')
r = corr(y1,y2);
title({'No Attn Mice';sprintf('Dist Only %s Model, r=%s',...
    modelName{choiceModInd},num2str(round(r,2,'significant')))})
r_noAttn = nan(1,nboot);
for iboot = 1:nboot
    n = length(y1);
    ind = randsample(1:n,n,1);
    y1_samp = y1(ind);
    y2_samp = y2(ind);
    r_noAttn(iboot) = corr(y1_samp,y2_samp);
end
subplot 236; hold on
bar(1:2,[mean(r_attn),mean(r_noAttn)])
[ylerr,yuerr] = ciFromBoot(r_attn',95);
errorbar(1,mean(r_attn),ylerr,yuerr,'.','MarkerSize',20)
[ylerr,yuerr] = ciFromBoot(r_noAttn',95);
errorbar(2,mean(r_noAttn),ylerr,yuerr,'.','MarkerSize',20)
figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
xtickangle(-45)
figYAxis([],'Vis:Aud Weight Corr',[0 1])
figAxForm
title(sprintf('Dist Only %s Model',mdl_noAttn(iav).mdl(imod).name))
print([fnout 'modelWeights_pcs_all_distOnly'],'-dpdf','-fillpage')

%% how does pc weight change when moving from inclusive model to distractor only model?
figure
suptitle('PC weights')
for iav = 1:2
    if iav == 1
        iplot = 0;
    else
        iplot = 3;
    end
    subplot(2,3,1+iplot); hold on
    y1 = cell2mat(mdl_attn(iav).mdl(choiceModInd).weights_pcs(goodPerfInd_av_distOnly_attn)');
    y2 = cell2mat(mdl_attn(iav).mdl(choiceModInd).weights_pcs_distOnly(goodPerfInd_av_distOnly_attn)');
    plot(y1,y2,'k.','MarkerSize',10)
    hold on
    mdl_fitline = fitlm(y1,y2);
    rsq = mdl_fitline.Rsquared.Ordinary;
    yfit = predict(mdl_fitline,y1);
    plot(y1,yfit,'-')
    figXAxis([],'All Stim Weight',weight_lim_pcs)
    figYAxis([],'Dist Onlly Weight',weight_lim_pcs)
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    r = corr(y1,y2);
    title({'Attn Mice';sprintf('%s %s Model, r=%s',...
        avName{iav},modelName{choiceModInd},num2str(round(r,2,'significant')))})
    r_attn = nan(1,nboot);
    for iboot = 1:nboot
        n = length(y1);
        ind = randsample(1:n,n,1);
        y1_samp = y1(ind);
        y2_samp = y2(ind);
        r_attn(iboot) = corr(y1_samp,y2_samp);
    end

    subplot(2,3,2+iplot); hold on
    y1 = cell2mat(mdl_noAttn(iav).mdl(choiceModInd).weights_pcs(goodPerfInd_av_distOnly_noAttn)');
    y2 = cell2mat(mdl_noAttn(iav).mdl(choiceModInd).weights_pcs_distOnly(goodPerfInd_av_distOnly_noAttn)');
    plot(y1,y2,'k.','MarkerSize',10)
    hold on
    mdl_fitline = fitlm(y1,y2);
    rsq = mdl_fitline.Rsquared.Ordinary;
    yfit = predict(mdl_fitline,y1);
    plot(y1,yfit,'-')
    figXAxis([],'All Stim Weight',weight_lim_pcs)
    figYAxis([],'Dist Onlly Weight',weight_lim_pcs)
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    r = corr(y1,y2);
    title({'No Attn Mice';sprintf('%s %s Model, r=%s',...
        avName{iav},modelName{choiceModInd},num2str(round(r,2,'significant')))})
    r_noAttn = nan(1,nboot);
    for iboot = 1:nboot
        n = length(y1);
        ind = randsample(1:n,n,1);
        y1_samp = y1(ind);
        y2_samp = y2(ind);
        r_noAttn(iboot) = corr(y1_samp,y2_samp);
    end

    subplot(2,3,3+iplot); hold on
    bar(1:2,[mean(r_attn),mean(r_noAttn)])
    [ylerr,yuerr] = ciFromBoot(r_attn',95);
    errorbar(1,mean(r_attn),ylerr,yuerr,'.','MarkerSize',20)
    [ylerr,yuerr] = ciFromBoot(r_noAttn',95);
    errorbar(2,mean(r_noAttn),ylerr,yuerr,'.','MarkerSize',20)
    figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
    xtickangle(-45)
    figYAxis([],'Weight Corr',[0 1])
    figAxForm
    title(sprintf('%s Model',mdl_noAttn(iav).mdl(imod).name))
end
%% target only choice model performance

for iav = 1:2
    mdl_attn(iav).mdl(choiceModInd).pctCorr_tarOnly = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorr_tarOnly_otherAV = nan(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(choiceModInd).pctCorr_tarOnly(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_tarOnly_holdout;
        mdl_attn(iav).mdl(choiceModInd).pctCorr_tarOnly_otherAV(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_tarOnly_otherAV;
    end
    
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_tarOnly = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_tarOnly_otherAV = nan(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_tarOnly(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_tarOnly_holdout;
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_tarOnly_otherAV(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_tarOnly_otherAV;
    end
end
figure
suptitle('Target Only Trained')
for iav = 1:2
    if iav == 1
        otherAV = 2;
    else
        otherAV = 1;
    end
    subplot(2,2,iav)
    y = cat(1,mdl_attn(iav).mdl(choiceModInd).pctCorr_tarOnly,...
        mdl_attn(iav).mdl(choiceModInd).pctCorr_tarOnly_otherAV);
    ind = y(1,:)>pctCorrThresh;
    plot(1:2,y(:,ind),'k-')
    hold on
    errorbar(1:2,mean(y(:,ind),2),ste(y(:,ind),2),'.','MarkerSize',20);
    [~,p] = ttest(y(1,:),y(2,:));
    [~,p_chance] = ttest(y');
    figXAxis([],'Test',[0 3],1:2,{avName{iav},avName{otherAV}})
    figYAxis([],'Pct Correct',[0 1])
    figAxForm
    hline(0.5,'k-')
    hline(pctCorrThresh,'k:')
    title(sprintf('Attn, %s Choice Model, p=%s',avName{iav},num2str(round(p,2,'significant'))))    
    text([0:2],ones(1,3)*.2,{'pCh',num2str(round(p_chance(1),2,'significant')),...
        num2str(round(p_chance(2),2,'significant'))})
    
    subplot(2,2,iav+2)
    y = cat(1,mdl_noAttn(iav).mdl(choiceModInd).pctCorr_tarOnly,...
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_tarOnly_otherAV);
    ind = y(1,:)>pctCorrThresh;
    plot(1:2,y(:,ind),'k-')
    hold on
    errorbar(1:2,mean(y(:,ind),2),ste(y(:,ind),2),'.','MarkerSize',20);
    [~,p] = ttest(y(1,:),y(2,:));
    [~,p_chance] = ttest(y');
    figXAxis([],'Test',[0 3],1:2,{avName{iav},avName{otherAV}})
    figYAxis([],'Pct Correct',[0 1])
    figAxForm
    hline(0.5,'k-')
    hline(pctCorrThresh,'k:')
    title(sprintf('No Attn, %s Choice Model, p=%s',avName{iav},num2str(round(p,2,'significant'))))    
    text([0:2],ones(1,3)*.2,{'pCh',num2str(round(p_chance(1),2,'significant')),...
        num2str(round(p_chance(2),2,'significant'))})
end
print([fnout 'pctCorrectOtherAVModel_tarOnly'],'-dpdf','-fillpage')

%% target only model weights

for iav = 1:2
    mdl_attn(iav).mdl(choiceModInd).weights_tarOnly = cell(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(choiceModInd).weights_tarOnly{iexp} = ...
            dc_attn(iexp).av(iav).weightDetect_tarOnly;
    end
    mdl_attn(iav).mdl(choiceModInd).normWeights_tarOnly = ...
        cellfun(@(x) x./norm(x), ...
        mdl_attn(iav).mdl(choiceModInd).weights_tarOnly,'unif',0);
    
    mdl_noAttn(iav).mdl(choiceModInd).weights_tarOnly = cell(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(choiceModInd).weights_tarOnly{iexp} = ...
            dc_noAttn(iexp).av(iav).weightDetect_tarOnly;
    end
    mdl_noAttn(iav).mdl(choiceModInd).normWeights_tarOnly = ...
        cellfun(@(x) x./norm(x), ...
        mdl_noAttn(iav).mdl(choiceModInd).weights_tarOnly,'unif',0);
end

weightTooHigh_distOnly_attn = cellfun(@(x) max(abs(x)),...
    mdl_attn(visualTrials).mdl(choiceModInd).weights_tarOnly)> 30 |...
    cellfun(@(x) max(abs(x)),...
    mdl_attn(auditoryTrials).mdl(choiceModInd).weights_tarOnly) > 30;
weightTooHigh_distOnly_noAttn = cellfun(@(x) max(abs(x)),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).weights_tarOnly)> 30 |...
    cellfun(@(x) max(abs(x)),...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_tarOnly) > 30;
goodPerfInd_av_tarOnly_attn = (mdl_attn(visualTrials).mdl(choiceModInd).pctCorr_tarOnly > pctCorrThresh |...
    mdl_attn(auditoryTrials).mdl(choiceModInd).pctCorr_tarOnly > pctCorrThresh)...
    & ~weightTooHigh_distOnly_attn;
goodPerfInd_av_tarOnly_noAttn = (mdl_noAttn(visualTrials).mdl(choiceModInd).pctCorr_tarOnly > pctCorrThresh |...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).pctCorr_tarOnly > pctCorrThresh)...
    & ~weightTooHigh_distOnly_noAttn;

for iav = 1:2
    rng(0)
    mdl_attn(iav).mdl(choiceModInd).shuffWeights_tarOnly = cellfun(...
        @(x) x(randperm(length(x))),mdl_attn(iav).mdl(choiceModInd).weights_tarOnly,'unif',0);      
    mdl_noAttn(iav).mdl(choiceModInd).shuffWeights_tarOnly = cellfun(...
        @(x) x(randperm(length(x))),mdl_noAttn(iav).mdl(choiceModInd).weights_tarOnly,'unif',0);  
    if iav == 1
        for ishuff = 1:nShuff
%                 
            % shuffle normalized weights
            shuffWeights_vis = cellfun(@(x) x(randperm(length(x))),...
                mdl_attn(visualTrials).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_attn),'unif',0)';
            shuffWeights_aud = cellfun(@(x) x(randperm(length(x))),...
                mdl_attn(auditoryTrials).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_attn),'unif',0)';
            if ishuff == 1
                mdl_attn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_tarOnly = ...
                    nan(sum(cellfun(@length,shuffWeights_vis)),nShuff);
            end
            mdl_attn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_tarOnly(:,ishuff) = ...
                cell2mat(cellfun(@(x,y) ...
                rad2deg(atan(x./y)),shuffWeights_aud,...
                shuffWeights_vis,'unif',0));

            shuffWeights_vis = cellfun(@(x) x(randperm(length(x))),...
                mdl_noAttn(visualTrials).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_noAttn),'unif',0)';
            shuffWeights_aud = cellfun(@(x) x(randperm(length(x))),...
                mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_noAttn),'unif',0)';
            if ishuff == 1
                mdl_noAttn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_tarOnly = ...
                    nan(sum(cellfun(@length,shuffWeights_vis)),nShuff);
            end
            mdl_noAttn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_tarOnly(:,ishuff) = ...
                cell2mat(cellfun(@(x,y) ...
                rad2deg(atan(x./y)),shuffWeights_aud,...
                shuffWeights_vis,'unif',0));
        end
    end
end

figure
y1 = cell2mat(mdl_attn(visualTrials).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_attn)');
y2 = cell2mat(mdl_attn(auditoryTrials).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_attn)');
n_attn = length(y1);
%     [r_attn,p_attn,rCIl_attn,rCIu_attn] = corrcoef(y1,y2);
[r_attn,p_attn] = corr(y1,y2);
y1 = cell2mat(mdl_noAttn(visualTrials).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_noAttn)');
y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_noAttn)');
n_noAttn = length(y1);
[r_noAttn,p_noAttn] = corr(y1,y2);
bar(1,r_attn)
text(1,0.1,sigfigString(p_attn))
hold on
bar(2,r_noAttn)
text(2,0.1,sigfigString(p_noAttn))
figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
figYAxis([],'Vis-Aud Corr',[-0.2 1])
figAxForm
title(sprintf('%s Model, Targets Only, p=%s',modelName{choiceModInd},...
    sigfigString(compare_correlation_coefficients(...
    r_attn,n_attn,r_noAttn,n_noAttn))))
print([fnout 'avWeightCorr_tarOnly'],'-dpdf')

figure
suptitle('Targets Only Model')
for iav = 1:2
    if iav == 1
        iplot = 0;
    else
        iplot = 2;
    end
    subplot(1,4,iav)
    y = cell2mat(mdl_attn(iav).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_attn)');
    histogram(y,weightEdges,'Normalization','probability')
    figXAxis([],'Weight',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Fraction of Cells',[0 .7])
    figAxForm
    title({'Attn Mice';sprintf('%s %s Model',...
        avName{iav},modelName{choiceModInd})})
    vline(0,'k:')
    subplot(1,4,iav+2)
    y = cell2mat(mdl_noAttn(iav).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_noAttn)');
    histogram(y,weightEdges,'Normalization','probability')
    figXAxis([],'Weight',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Fraction of Cells',[0 .7])
    figAxForm
    title({'No Attn Mice';sprintf('%s %s Model',...
        avName{iav},modelName{choiceModInd})})
    vline(0,'k:')
end
print([fnout 'weightHist_tarOnly'],'-dpdf','-fillpage')

figure
suptitle('Targets Only Model')
for iav = 1:2
    subplot(1,2,iav)
    y1 = cell2mat(mdl_attn(iav).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_attn)');
    y2 = cell2mat(mdl_noAttn(iav).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_noAttn)');
    cdfplot(y1)
    hold on
    cdfplot(y2)
    figXAxis([],'Weight',[-1 1])
    figYAxis([],'Fraction of Cells',[0 1])
    figAxForm
    title({'Attn Mice';sprintf('%s %s Model',...
        avName{iav},modelName{choiceModInd})})
    vline(0,'k:')
    legend({'Attn','No Attn'},'location','southeast')
end
print([fnout 'weightCDF_tarOnly'],'-dpdf','-fillpage')

mdlWeightsSub_fig = figure;
mdlWeights_fig = figure;
suptitle({'Targets Only Model,Weights';...
    sprintf('pct corr must be > %s in vis & aud stim models',...
    num2str(pctCorrThresh))})

figure(mdlWeights_fig)
subplot(3,2,1)
y1 = cell2mat(mdl_attn(visualTrials).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_attn)');
y2 = cell2mat(mdl_attn(auditoryTrials).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_attn)');
plot(y1,y2,'k.','MarkerSize',10)
hold on
mdl_fitline = fitlm(y1,y2);
rsq = mdl_fitline.Rsquared.Ordinary;
yfit = predict(mdl_fitline,y1);
plot(y1,yfit,'-')
figXAxis([],'Visual Weight',[weightEdges(1) weightEdges(end)])
figYAxis([],'Auditory Weight',[weightEdges(1) weightEdges(end)])
figAxForm
hline(0,'k:')
vline(0,'k:')
r = corr(y1,y2);
title({'Attn Mice';sprintf('%s Model, r=%s',...
    modelName{choiceModInd},num2str(round(r,2,'significant')))})

subplot(3,2,3)
weightAVangle = rad2deg(atan(y2./y1));
histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
s = skewness(weightAVangle);
title({'Attn Mice';sprintf('%s Model, skew=%s,ks=%s',...
    modelName{choiceModInd},num2str(round(s,2,'significant')))})

subplot(3,2,5)
%     y = mean(mdl_attn(visualTrials).mdl(choiceModInd).normWeightAVangle_shuff,2);
%     histogram(y,oriEdges,'Normalization','probability')
avAngleShuff = mdl_attn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_tarOnly;
h = nan(length(oriEdges_noAbs)-1,nShuff);
for ishuff = 1:nShuff
    h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges_noAbs,'Normalization','Probability');
end    
h_sort = sort(h,2);
y = mean(h_sort,2);
yerrl = y-h_sort(:,round(nShuff*0.025));
yerru = h_sort(:,round(nShuff*0.975))-y;
shadedErrorBar_chooseColor(oriEdges_noAbs(1:end-1),y,[yerrl';yerru'],[0 0 0]);
hold on
plot(oriEdges_noAbs(1:end-1),histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability'),'.-')
controlSub_attn = histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability') - y';
figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
s = skewness(y);
[~,p] = kstest2(h_sort(:),weightAVangle);
title({'Attn Mice-Shuff';sprintf('%s Model, skew=%s,p=%s',...
    modelName{choiceModInd},num2str(round(s,2,'significant')),...
    sigfigString(p))})

figure(mdlWeightsSub_fig)
plot(oriEdges_noAbs(1:end-1),controlSub_attn,'.-')
%     shadedErrorBar_chooseColor(oriEdges(1:end-1),controlSub_attn,[yerrl';yerru'],[0 0 0]);

figure(mdlWeights_fig)
subplot(3,2,2)
y1 = cell2mat(mdl_noAttn(visualTrials).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_noAttn)');
y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_noAttn)');
plot(y1,y2,'k.','MarkerSize',10)
hold on
mdl_fitline = fitlm(y1,y2);
rsq = mdl_fitline.Rsquared.Ordinary;
yfit = predict(mdl_fitline,y1);
plot(y1,yfit,'-')
figXAxis([],'Visual Weight',[weightEdges(1) weightEdges(end)])
figYAxis([],'Auditory Weight',[weightEdges(1) weightEdges(end)])
figAxForm
hline(0,'k:')
vline(0,'k:')
r = corr(y1,y2);
title({'No Attn Mice';sprintf('%s Model, r=%s',...
    modelName{choiceModInd},num2str(round(r,2,'significant')))})

subplot(3,2,4)
weightAVangle = rad2deg(atan(y2./y1));
histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
s = skewness(weightAVangle);
title({'No Attn Mice';sprintf('%s Model, skew=%s',...
    modelName{choiceModInd},num2str(round(s,2,'significant')))})

subplot(3,2,6)
%     y = mean(mdl_noAttn(visualTrials).mdl(choiceModInd).normWeightAVangle_shuff,2);
%     histogram(y,oriEdges,'Normalization','probability')
avAngleShuff = mdl_noAttn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_tarOnly;
h = nan(length(oriEdges_noAbs)-1,nShuff);
for ishuff = 1:nShuff
    h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges_noAbs,'Normalization','Probability');
end    
h_sort = sort(h,2);
y = mean(h_sort,2);
yerrl = y-h_sort(:,round(nShuff*0.025));
yerru = h_sort(:,round(nShuff*0.975))-y;
shadedErrorBar_chooseColor(oriEdges_noAbs(1:end-1),y,[yerrl';yerru'],[0 0 0]);
hold on
plot(oriEdges_noAbs(1:end-1),histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability'),'.-')
controlSub_noAttn = histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability') - y';
figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
s = skewness(y);
title({'No Attn Mice-Shuff';sprintf('%s Model, skew=%s',...
    modelName{choiceModInd},num2str(round(s,2,'significant')))})

figure(mdlWeightsSub_fig)
hold on
plot(oriEdges_noAbs(1:end-1),controlSub_noAttn,'.-')
figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
figYAxis([],'Control Subtr. Fraction of Cells',[-0.15 0.15])
figAxForm
hline(0,'k:')
title(sprintf('%s Model, Targets Only',modelName{choiceModInd}))

figure(mdlWeights_fig)
print([fnout 'modelWeights_tarOnly'],'-dpdf','-fillpage')
figure(mdlWeightsSub_fig)
print([fnout 'modelWeights_ciSub_tarOnly'],'-dpdf')
%% weights verses selectivity
mdl_attn(visualTrials).si = cellfun(@(x,y) x(y == 1),...
    siPerExpt_attn,{dc_attn.cellInd},'unif',0);
mdl_noAttn(visualTrials).si = cellfun(@(x,y) x(y == 1),...
    siPerExpt_noAttn,{dc_noAttn.cellInd},'unif',0);
% inclusive model
figure
for imod = 1:2
    y1 = cell2mat(mdl_attn(visualTrials).mdl(imod).weights(goodPerfInd_av_attn)');
    y2 = cell2mat(mdl_attn(auditoryTrials).mdl(imod).weights(goodPerfInd_av_attn)');
    weightAVangle = rad2deg(atan(y2./y1));
    si = cell2mat(mdl_attn(visualTrials).si(goodPerfInd_av_attn));
    
    subplot(4,2,imod)
    plot(si,weightAVangle,'k.')
    figXAxis([],'V-A selectivity',[-10 10])
    figYAxis([],'Angle of A/V Weight',[-90 90],-90:30:90)
    figAxForm
    title(sprintf('Attn, %s Model, All Inclusive',modelName{imod}))
    
    y1 = cell2mat(mdl_noAttn(visualTrials).mdl(imod).weights(goodPerfInd_av_noAttn)');
    y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(imod).weights(goodPerfInd_av_noAttn)');
    weightAVangle = rad2deg(atan(y2./y1));
    si = cell2mat(mdl_noAttn(visualTrials).si(goodPerfInd_av_noAttn));
    
    subplot(4,2,imod+4)
    plot(si,weightAVangle,'k.')
    figXAxis([],'V-A selectivity',[-10 10])
    figYAxis([],'Angle of A/V Weight',[-90 90],-90:30:90)
    figAxForm
    title(sprintf('No Attn, %s Model, All Inclusive',modelName{imod}))
end

% distractor only choice model
y1 = cell2mat(mdl_attn(visualTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_attn)');
y2 = cell2mat(mdl_attn(auditoryTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_attn)');
weightAVangle = rad2deg(atan(y2./y1));
si = cell2mat(mdl_attn(visualTrials).si(goodPerfInd_av_distOnly_attn));

subplot(4,2,3)
plot(si,weightAVangle,'k.')
figXAxis([],'V-A selectivity',[-10 10])
figYAxis([],'Angle of Norm. A/V Weight',[-0 90],0:15:90)
figAxForm
title(sprintf('Attn, %s Model, Distractors Only',modelName{imod}))

y1 = cell2mat(mdl_noAttn(visualTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_noAttn)');
y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_noAttn)');
weightAVangle = rad2deg(atan(y2./y1));
si = cell2mat(mdl_noAttn(visualTrials).si(goodPerfInd_av_distOnly_noAttn));

subplot(4,2,7)
plot(si,weightAVangle,'k.')
figXAxis([],'V-A selectivity',[-10 10])
    figYAxis([],'Angle of A/V Weight',[-90 90],-90:30:90)
figAxForm
title(sprintf('No Attn, %s Model, Distractors Only',modelName{imod}))

% target only choice model
y1 = cell2mat(mdl_attn(visualTrials).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_attn)');
y2 = cell2mat(mdl_attn(auditoryTrials).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_attn)');
weightAVangle = rad2deg(atan(y2./y1));
si = cell2mat(mdl_attn(visualTrials).si(goodPerfInd_av_tarOnly_attn));

subplot(4,2,4)
plot(si,weightAVangle,'k.')
figXAxis([],'V-A selectivity',[-10 10])
    figYAxis([],'Angle of A/V Weight',[-90 90],-90:30:90)
figAxForm
title(sprintf('Attn, %s Model, Targets Only',modelName{imod}))

y1 = cell2mat(mdl_noAttn(visualTrials).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_noAttn)');
y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_tarOnly(goodPerfInd_av_tarOnly_noAttn)');
weightAVangle = rad2deg(atan(y2./y1));
si = cell2mat(mdl_noAttn(visualTrials).si(goodPerfInd_av_tarOnly_noAttn));

subplot(4,2,8)
plot(si,weightAVangle,'k.')
figXAxis([],'V-A selectivity',[-10 10])
    figYAxis([],'Angle of A/V Weight',[-90 90],-90:30:90)
figAxForm
title(sprintf('No Attn, %s Model, Targets Only',modelName{imod}))

print([fnout 'siXweights'],'-dpdf','-fillpage')
%% performance across response windows
nBaselineFr = 30;
nwins=15;
winInd = [1:2,6:nwins];
movWinLabelFr = 30:(30+nwins-1);
movWinLabelMs = (movWinLabelFr - (nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
respWinTT = ([respwin(1) respwin(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
baseWinTT = (...
    [basewin_0(1) basewin_0(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);

for iav = 1:2
    mdl_attn(iav).mdl(stimModInd).pctCorr_respWin = nan(nwins,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorr_respWin = nan(nwins,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).pctCorr_respWin(:,iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectTargetMovRespWin_holdout;        
        mdl_attn(iav).mdl(choiceModInd).pctCorr_respWin(:,iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetectMovRespWin_holdout;
    end
    mdl_noAttn(iav).mdl(stimModInd).pctCorr_respWin = nan(nwins,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_respWin = nan(nwins,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(stimModInd).pctCorr_respWin(:,iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectTargetMovRespWin_holdout;        
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_respWin(:,iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetectMovRespWin_holdout;
    end
end

figure
x = movWinLabelMs(winInd);
targetTimes = x(x>0);
for iav = 1:2
    for imod = 1:2
        if imod == 1
            iplot = 0;
        else
            iplot = 2;
        end
        subplot(2,2,iplot+iav)
        hold on
        y1 = mdl_attn(iav).mdl(imod).pctCorr_respWin(winInd,:);
        ind1 = mdl_attn(iav).mdl(imod).pctCorr_all > pctCorrThresh;
        y2 = mdl_noAttn(iav).mdl(imod).pctCorr_respWin(winInd,:);
        ind2 = mdl_noAttn(iav).mdl(imod).pctCorr_all > pctCorrThresh;
        y = cat(2,y1(:,ind1),y2(:,ind2));
        h = plot(x,y,'k-');
        errorbar(x,mean(y,2),ste(y,2),'.-')
        firstTimes = nan(1,nexp_attn+nexp_noAttn);
        for iexp = 1:nexp_attn
            firstTimes(iexp) = targetTimes(find(y1(x>0,iexp)>pctCorrThresh,1));
        end
        for iexp = 1:nexp_noAttn
            firstTimes(iexp+nexp_attn) = targetTimes(find(y2(x>0,iexp)>pctCorrThresh,1));
        end
        errorbar(mean(firstTimes(ind)),0.9,[],[],...
            ste(firstTimes(ind),2),ste(firstTimes(ind),2),'.')
        figXAxis([],'Resp. Win. Center (Relative to Stim)',[x(1) x(end)],x(1:2:length(x)),round(x(1:2:length(x))))
        hline(0.5,'k-')
        figYAxis([],'Frac. Correct',[0 1])
        vline(0,'k:')
        figAxForm
        title(sprintf('All Mice %s %s Model',avName{iav},modelName{imod}))
        vline(respWinTT,'r--')
        vline(baseWinTT,'k--')
        hline(pctCorrThresh,'k:')
        if iav == 1
            vline(visRTwindow(1),'b--')
        else
            vline(audRTwindow(1),'b--')
        end
        
    end
end

print([fnout 'pctCorrectXrespWin'],'-dpdf','-fillpage')   

%% stim and choice weight correlation
setFigParams4Print('landscape')
figure
for iav = 1:2
    ind = mdl_attn(iav).mdl(stimModInd).pctCorr_all > pctCorrThresh;
    subplot(2,4,iav)
    x = cell2mat(mdl_attn(iav).mdl(stimModInd).weights(ind)');
    y = cell2mat(mdl_attn(iav).mdl(choiceModInd).weights(ind)');
    hold on
    plot(x,y,'.','MarkerSize',10)
    
    mdl_fitline = fitlm(x,y);
    rsq = mdl_fitline.Rsquared.Ordinary;
    yfit = predict(mdl_fitline,x);
    plot(x,yfit,'-')
    [r,p] = corr(x,y);
    title(sprintf('Attn, %s Models,r=%s,p=%s',avName{iav},...
        num2str(round(r,2,'significant')),num2str(round(p,2,'significant'))))    
    figXAxis([],'Stimulus Weight',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Choice Weight',[weightEdges(1) weightEdges(end)])
    figAxForm
    
    ind = mdl_noAttn(iav).mdl(stimModInd).pctCorr_all > pctCorrThresh;
    subplot(2,4,iav+2)
    x = cell2mat(mdl_noAttn(iav).mdl(stimModInd).weights(ind)');
    y = cell2mat(mdl_noAttn(iav).mdl(choiceModInd).weights(ind)');
    hold on
    plot(x,y,'.','MarkerSize',10)
    
    mdl_fitline = fitlm(x,y);
    rsq = mdl_fitline.Rsquared.Ordinary;
    yfit = predict(mdl_fitline,x);
    plot(x,yfit,'-')
    [r,p] = corr(x,y);
    title(sprintf('No Attn, %s Models,r=%s,p=%s',avName{iav},...
        num2str(round(r,2,'significant')),num2str(round(p,2,'significant'))))    
    figXAxis([],'Stimulus Weight',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Choice Weight',[weightEdges(1) weightEdges(end)])
    figAxForm
end
subplot 245
ind = goodPerfInd_av_attn;
x = cellfun(@corr,mdl_attn(visualTrials).mdl(stimModInd).weights(ind),...
    mdl_attn(visualTrials).mdl(choiceModInd).weights(ind));
y = cellfun(@corr,mdl_attn(auditoryTrials).mdl(stimModInd).weights(ind),...
    mdl_attn(auditoryTrials).mdl(choiceModInd).weights(ind));
plot(x,y,'.','MarkerSize',20)
hold on
plot(0:1,0:1,'k--')
figXAxis([],'Vis Stim W:Choice W Corr',[0 1])
figYAxis([],'Aud Stim W:Choice W Corr',[0 1])
figAxForm
title('Attn Mice')

subplot 246
plot(1:2,[x',y'],'k-','MarkerSize',20)
hold on
errorbar(1,mean(x),ste(x,2),'.','MarkerSize',20)
% plot(ones(1,sum(ind)).*2,y(ind),'.','MarkerSize',20)
errorbar(2,mean(y),ste(y,2),'.','MarkerSize',20)
figXAxis([],'',[0 3],1:2,avName)
figYAxis([],'Stim:Choice Weight Corr',[0 1])
figAxForm

subplot 247
ind = goodPerfInd_av_noAttn;
x = cellfun(@corr,mdl_noAttn(visualTrials).mdl(stimModInd).weights,...
    mdl_noAttn(visualTrials).mdl(choiceModInd).weights);
y = cellfun(@corr,mdl_noAttn(auditoryTrials).mdl(stimModInd).weights,...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights);
plot(x,y,'.','MarkerSize',20)
hold on
plot(0:1,0:1,'k--')
figXAxis([],'Vis Stim W:Choice W Corr',[0 1])
figYAxis([],'Aud Stim W:Choice W Corr',[0 1])
figAxForm
title('No Attn Mice')

subplot 248
plot(1:2,[x(ind)',y(ind)'],'k-','MarkerSize',20)
hold on
errorbar(1,mean(x(ind)),ste(x(ind),2),'.','MarkerSize',20)
errorbar(2,mean(y(ind)),ste(y(ind),2),'.','MarkerSize',20)
figXAxis([],'',[0 3],1:2,avName)
figYAxis([],'Stim:Choice Weight Corr',[0 1])
figAxForm

print([fnout 'stimVsChoiceWeights'],'-dpdf','-fillpage')

%% bootstrapped model correlation analysis
setFigParams4Print('portrait')
msColors_attn = cat(1,brewermap(9,'Set1'),brewermap(nexp_attn-9,'Set2'));
msColors_noAttn = cat(1,brewermap(9,'Set1'),brewermap(nexp_noAttn-9,'Set2'));
fig_choice = figure;
suptitle('Attn  Choice Model')
fig_stim = figure;
suptitle('Attn Stimulus Model')

for iexp = 1:nexp_attn
    for iav = 1:2
        figure(fig_choice)
        bootModel = dc_attn(iexp).av(iav).choiceModel_fixedStim_bootStruct;
        stimW = cellfun(@(x) x(end),bootModel.choiceModelWeights_boot);
        ind = abs(stimW) > 20;
        stimW(ind) = nan;
        c = bootModel.stimChoiceCorr_boot;
        c(ind) = nan;
        p_pcsOnly = bootModel.pctCorr_choice_fixedStim_pcsOnly_boot;
        p_pcsOnly(ind) = nan;
        p_stimOnly = bootModel.pctCorr_choice_fixedStimOnly_boot;
        p_stimOnly(ind) = nan;
        p_all = bootModel.pctCorr_choice_fixedStim_boot;
        p_all(ind) = nan;

        subplot(4,2,iav)
        hold on
        h = plot(c,p_stimOnly,'.','MarkerSize',10);
        h.Color = msColors_attn(iexp,:);
        figXAxis([],'Stimulus:Choice Correlation',[0 1])
        figYAxis([],'Pct Correct',[0 1])
        figAxForm
        hline(pctCorrThresh,'k--')
        title(sprintf('%s - Stimulus Only Parameter',avName{iav}))

        subplot(4,2,iav+2)
        hold on
        h = plot(c,p_all,'.','MarkerSize',10);
        h.Color = msColors_attn(iexp,:);
        figXAxis([],'Stimulus:Choice Correlation',[0 1])
        figYAxis([],'Pct Correct',[0 1])
        figAxForm
        hline(pctCorrThresh,'k--')
        title(sprintf('%s - PCs+Stimulus',avName{iav}))

        subplot(4,2,iav+4)
        hold on
        h = plot(c,p_pcsOnly,'.','MarkerSize',10);
        h.Color = msColors_attn(iexp,:);
        figXAxis([],'Stimulus:Choice Correlation',[0 1])
        figYAxis([],'Pct Correct',[0 1])
        figAxForm
        hline(pctCorrThresh,'k--')
        title(sprintf('%s - (PCs+Stimulus)-Stimulus',avName{iav}))

        subplot(4,2,iav+6)
        hold on
        h = plot(c,stimW,'.','MarkerSize',10);
        h.Color = msColors_attn(iexp,:);
        figXAxis([],'Stimulus:Choice Correlation',[0 1])
        figYAxis([],'Stimulus Weight',[])
        figAxForm
        title('PCs+Stimulus')
        title(sprintf('%s - PCs+Stimulus',avName{iav}))
        
        figure(fig_stim)
        bootModel = dc_attn(iexp).av(iav).stimModel_fixedChoice_bootStruct;
        choiceW = cellfun(@(x) x(end),bootModel.stimModelWeights_boot);
        ind = abs(choiceW) > 20;
        choiceW(ind) = nan;
        c = bootModel.stimChoiceCorr_boot;
        c(ind) = nan;
        p_pcsOnly = bootModel.pctCorr_stim_fixedChoice_pcsOnly_boot;
        p_pcsOnly(ind) = nan;
        p_stimOnly = bootModel.pctCorr_stim_fixedChoiceOnly_boot;
        p_stimOnly(ind) = nan;
        p_all = bootModel.pctCorr_stim_fixedChoice_boot;
        p_all(ind) = nan;

        subplot(4,2,iav)
        hold on
        h = plot(c,p_stimOnly,'.','MarkerSize',10);
        h.Color = msColors_attn(iexp,:);
        figXAxis([],'Stimulus:stim Correlation',[0 1])
        figYAxis([],'Pct Correct',[0 1])
        figAxForm
        hline(pctCorrThresh,'k--')
        title(sprintf('%s - Stimulus Only Parameter',avName{iav}))

        subplot(4,2,iav+2)
        hold on
        h = plot(c,p_all,'.','MarkerSize',10);
        h.Color = msColors_attn(iexp,:);
        figXAxis([],'Stimulus:stim Correlation',[0 1])
        figYAxis([],'Pct Correct',[0 1])
        figAxForm
        hline(pctCorrThresh,'k--')
        title(sprintf('%s - PCs+Stimulus',avName{iav}))

        subplot(4,2,iav+4)
        hold on
        h = plot(c,p_pcsOnly,'.','MarkerSize',10);
        h.Color = msColors_attn(iexp,:);
        figXAxis([],'Stimulus:stim Correlation',[0 1])
        figYAxis([],'Pct Correct',[0 1])
        figAxForm
        hline(pctCorrThresh,'k--')
        title(sprintf('%s - (PCs+Stimulus)-Stimulus',avName{iav}))

        subplot(4,2,iav+6)
        hold on
        h = plot(c,choiceW,'.','MarkerSize',10);
        h.Color = msColors_attn(iexp,:);
        figXAxis([],'Stimulus:stim Correlation',[0 1])
        figYAxis([],'Stimulus Weight',[])
        figAxForm
        title('PCs+Stimulus')
        title(sprintf('%s - PCs+Stimulus',avName{iav}))
    end
end
figure(fig_choice)
print([fnout 'choiceModelBoot_attn'],'-dpdf','-fillpage')
figure(fig_stim)
print([fnout 'stimModelBoot_attn'],'-dpdf','-fillpage')
figure
suptitle('No Attn  Choice Model')
for iexp = 1:nexp_noAttn
    for iav = 1:2
    
        bootModel = dc_noAttn(iexp).av(iav).choiceModel_fixedStim_bootStruct;
        stimW = cellfun(@(x) x(end),bootModel.choiceModelWeights_boot);
        ind = abs(stimW) > 10;
        stimW(ind) = nan;
        c = bootModel.stimChoiceCorr_boot;
        c(ind) = nan;
        p_pcsOnly = bootModel.pctCorr_choice_fixedStim_pcsOnly_boot;
        p_pcsOnly(ind) = nan;
        p_stimOnly = bootModel.pctCorr_choice_fixedStimOnly_boot;
        p_stimOnly(ind) = nan;
        p_all = bootModel.pctCorr_choice_fixedStim_boot;
        p_all(ind) = nan;

        subplot(4,2,iav)
        hold on
        h = plot(c,p_stimOnly,'.','MarkerSize',10);
        h.Color = msColors_noAttn(iexp,:);
        figXAxis([],'Stimulus:Choice Correlation',[0 1])
        figYAxis([],'Pct Correct',[0 1])
        figAxForm
        hline(pctCorrThresh,'k--')
        title(sprintf('%s - Stimulus Only Parameter',avName{iav}))

        subplot(4,2,iav+2)
        hold on
        h = plot(c,p_all,'.','MarkerSize',10);
        h.Color = msColors_noAttn(iexp,:);
        figXAxis([],'Stimulus:Choice Correlation',[0 1])
        figYAxis([],'Pct Correct',[0 1])
        figAxForm
        hline(pctCorrThresh,'k--')
        title(sprintf('%s - PCs+Stimulus',avName{iav}))

        subplot(4,2,iav+4)
        hold on
        h = plot(c,p_pcsOnly,'.','MarkerSize',10);
        h.Color = msColors_noAttn(iexp,:);
        figXAxis([],'Stimulus:Choice Correlation',[0 1])
        figYAxis([],'Pct Correct',[0 1])
        figAxForm
        hline(pctCorrThresh,'k--')
        title(sprintf('%s - (PCs+Stimulus)-Stimulus',avName{iav}))

        subplot(4,2,iav+6)
        hold on
        h = plot(c,stimW,'.','MarkerSize',10);
        h.Color = msColors_noAttn(iexp,:);
        figXAxis([],'Stimulus:Choice Correlation',[0 1])
        figYAxis([],'Stimulus Weight',[])
        figAxForm
        title('PCs+Stimulus')
        title(sprintf('%s - PCs+Stimulus',avName{iav}))
    end
end
print([fnout 'choiceModelBoot_noAttn'],'-dpdf','-fillpage')

figure
suptitle(sprintf('Only bootstraps with pct correct pcs alone > %s',num2str(pctCorrThresh)))
for iexp = 1:nexp_attn
    for iav = 1:2
        bootModel = dc_attn(iexp).av(iav).choiceModel_fixedStim_bootStruct;        
        subplot(4,2,iav)
        hold on
        ind = bootModel.pctCorr_pcs > pctCorrThresh;
        y = bootModel.pctCorr_pcs(ind) - ...
            bootModel.pctCorr_choice_fixedStim_pcsOnly_boot(ind);
        x = bootModel.stimChoiceCorr_boot(ind);
        h = plot(x,y,'.','MarkerSize',10);
        h.Color = msColors_attn(iexp,:);
        figXAxis([],'Stimulus:Choice Correlation',[0 1])
        figYAxis([],'Pct Correct PCs - Pct Correct PCs from fixed model',[-0.25 1])
        figAxForm
        title(sprintf('Attn %s Choice Model', avName{iav}))
        
        
        bootModel = dc_attn(iexp).av(iav).stimModel_fixedChoice_bootStruct; 
        subplot(4,2,iav+2)
        hold on
        ind = bootModel.pctCorr_pcs > pctCorrThresh;
        y = bootModel.pctCorr_pcs(ind) - ...
            bootModel.pctCorr_stim_fixedChoice_pcsOnly_boot(ind);
        x = bootModel.stimChoiceCorr_boot(ind);
        h = plot(x,y,'.','MarkerSize',10);
        h.Color = msColors_attn(iexp,:);
        figXAxis([],'Stimulus:Choice Correlation',[0 1])
        figYAxis([],'Pct Correct PCs - Pct Correct PCs from fixed model',[-0.25 1])
        figAxForm
        title(sprintf('Attn %s Stimulus Model', avName{iav}))
        
    end    
end
for iexp = 1:nexp_noAttn
    for iav = 1:2
        
        bootModel = dc_noAttn(iexp).av(iav).choiceModel_fixedStim_bootStruct;        
        subplot(4,2,iav+4)
        hold on
        ind = bootModel.pctCorr_pcs > pctCorrThresh;
        y = bootModel.pctCorr_pcs(ind) - ...
            bootModel.pctCorr_choice_fixedStim_pcsOnly_boot(ind);
        x = bootModel.stimChoiceCorr_boot(ind);
        h = plot(x,y,'.','MarkerSize',10);
        h.Color = msColors_noAttn(iexp,:);
        figXAxis([],'Stimulus:Choice Correlation',[0 1])
        figYAxis([],'Pct Correct PCs - Pct Correct PCs from fixed model',[-0.25 1])
        figAxForm
        title(sprintf('No Attn %s Choice Model', avName{iav}))
        
        
        bootModel = dc_noAttn(iexp).av(iav).stimModel_fixedChoice_bootStruct; 
        subplot(4,2,iav+6)
        hold on
        ind = bootModel.pctCorr_pcs > pctCorrThresh;
        y = bootModel.pctCorr_pcs(ind) - ...
            bootModel.pctCorr_stim_fixedChoice_pcsOnly_boot(ind);
        x = bootModel.stimChoiceCorr_boot(ind);
        h = plot(x,y,'.','MarkerSize',10);
        h.Color = msColors_noAttn(iexp,:);
        figXAxis([],'Stimulus:Choice Correlation',[0 1])
        figYAxis([],'Pct Correct PCs - Pct Correct PCs from fixed model',[-0.25 1])
        figAxForm
        title(sprintf('No Attn %s Stimulus Model', avName{iav}))
    end
end
print([fnout 'modelBoot_noStimVsFixed'],'-dpdf','-fillpage')


corrBins = [-1 0:0.1:1];
corrColors = brewermap(11,'YlOrRd');
corrColors = corrColors(end-length(corrBins)+2:end,:);
figure
for iav = 1:2
    binW = cell(1,length(corrBins)-1);
    for iexp = 1:nexp_attn    
        bootModel = dc_attn(iexp).av(iav).choiceModel_fixedStim_bootStruct;
        [n,~,corrBinID] = histcounts(bootModel.stimChoiceCorr_boot,corrBins);
        stimW = cellfun(@(x) x(end),bootModel.choiceModelWeights_boot);
        
        for ibin = 1:length(corrBins)-1
            ind = corrBinID == ibin & stimW < 20;
            d = abs(cell2mat(cellfun(@(x) x(1:(end-1)),bootModel.choiceModelWeights_boot(ind),'unif',0)'));
            binW{ibin} = cat(1,binW{ibin},d);
        end
    end
    L_label = [];
    for ibin = 1:length(corrBins)-1        
        if ~isempty(binW{ibin})
            subplot(2,2,iav)
            hold on
            h = cdfplot(binW{ibin});
            h.Color = corrColors(ibin,:);
            L_label = cat(2,L_label,{[num2str(corrBins(ibin)) ' to ' num2str(corrBins(ibin+1))]});
        end
    end
    legend(L_label,'location','southeast')
    figXAxis([],'Neuron Weight Magnitude',[0 3])
    figYAxis([],'Fraction of Neurons',[0 1])
    figAxForm
    title(sprintf('Attn, %s Choice, PCs+Stim, Model',avName{iav}))
    
    binW = cell(1,length(corrBins)-1);
    for iexp = 1:nexp_noAttn    
        bootModel = dc_noAttn(iexp).av(iav).choiceModel_fixedStim_bootStruct;
        [n,~,corrBinID] = histcounts(bootModel.stimChoiceCorr_boot,corrBins);
        stimW = cellfun(@(x) x(end),bootModel.choiceModelWeights_boot);
        
        for ibin = 1:length(corrBins)-1
            ind = corrBinID == ibin & stimW < 20;
            d = abs(cell2mat(cellfun(@(x) x(1:(end-1)),bootModel.choiceModelWeights_boot(ind),'unif',0)'));
            binW{ibin} = cat(1,binW{ibin},d);
        end
    end
    L_label = [];
    for ibin = 1:length(corrBins)-1        
        if ~isempty(binW{ibin})
            subplot(2,2,iav+2)
            hold on
            h = cdfplot(binW{ibin});
            h.Color = corrColors(ibin,:);
            L_label = cat(2,L_label,{[num2str(corrBins(ibin)) ' to ' num2str(corrBins(ibin+1))]});
        end
    end
    legend(L_label,'location','southeast')
    figXAxis([],'Neuron Weight Magnitude',[0 3])
    figYAxis([],'Fraction of Neurons',[0 1])
    figAxForm
    title(sprintf('No Attn, %s Choice, PCs+Stim, Model',avName{iav}))
end
print([fnout 'choiceModelBoot_neuronWeights'],'-dpdf','-fillpage')

%% combo model performance
mdl_attn(visualTrials).mdl(stimModInd).comboWeight = cell(1,nexp_attn);
mdl_attn(visualTrials).mdl(choiceModInd).comboWeight = cell(1,nexp_attn);
mdl_noAttn(visualTrials).mdl(stimModInd).comboWeight = cell(1,nexp_noAttn);
mdl_noAttn(visualTrials).mdl(choiceModInd).comboWeight = cell(1,nexp_noAttn);
for iav = 1:2
    mdl_attn(iav).mdl(stimModInd).pctCorr_combo = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorr_combo = nan(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).pctCorr_combo(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectTarget_comboTrain;
        mdl_attn(iav).mdl(choiceModInd).pctCorr_combo(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_comboTrain;
        if iav==1
            mdl_attn(visualTrials).mdl(stimModInd).comboWeight{iexp} = ...
                dc_attn(iexp).comboTrainWeightTarget;
            mdl_attn(visualTrials).mdl(choiceModInd).comboWeight{iexp} = ...
                dc_attn(iexp).comboTrainWeightDetect;
        end
    end
    mdl_noAttn(iav).mdl(stimModInd).pctCorr_combo = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_combo = nan(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(stimModInd).pctCorr_combo(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectTarget_comboTrain;
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_combo(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_comboTrain;
        if iav==1
            mdl_noAttn(visualTrials).mdl(stimModInd).comboWeight{iexp} = ...
                dc_noAttn(iexp).comboTrainWeightTarget;
            mdl_noAttn(visualTrials).mdl(choiceModInd).comboWeight{iexp} = ...
                dc_noAttn(iexp).comboTrainWeightDetect;
        end
    end
end

for imod = 1:2
    mdl_attn(visualTrials).mdl(imod).normComboWeight = cellfun(@(x) ...
        x./max(abs(x)),mdl_attn(visualTrials).mdl(imod).comboWeight,'unif',0);
    mdl_noAttn(visualTrials).mdl(imod).normComboWeight = cellfun(@(x) ...
        x./max(abs(x)),mdl_noAttn(visualTrials).mdl(imod).comboWeight,'unif',0);
end

figure
for iav = 1:2
    for imod = 1:2
        if imod == 1
            iplot = 0;
        else
            iplot = 2;
        end
        subplot(4,2,iav+iplot)
        ind = mdl_attn(iav).mdl(imod).pctCorr_all > pctCorrThresh;
        y = cat(2,mdl_attn(iav).mdl(imod).pctCorr_all(ind)',...
            mdl_attn(iav).mdl(imod).pctCorr_combo(ind)');
        [~,p] = ttest(y(:,1),y(:,2));
        plot(1:2,y,'k-')
        hold on
        errorbar(1:2,mean(y,1),ste(y,1),'.')
        figXAxis([],'Train',[0 3],1:2,{avName{iav},'Vis+Aud'})
        figYAxis([],'Frac. Correct',[0 1])
        figAxForm
        hline(0.5,'k-')
        hline(pctCorrThresh,'k:')
        title(sprintf('Attn; Test %s, %s Model,p=%s',avName{iav},modelName{imod},...
            num2str(round(p,2,'significant'))))
        
        if imod == 1
            iplot = 4;
        else
            iplot = 6;
        end
        subplot(4,2,iav+iplot)
        ind = mdl_noAttn(iav).mdl(imod).pctCorr_all > pctCorrThresh;
        y = cat(2,mdl_noAttn(iav).mdl(imod).pctCorr_all(ind)',...
            mdl_noAttn(iav).mdl(imod).pctCorr_combo(ind)');
        [~,p] = ttest(y(:,1),y(:,2));
        plot(1:2,y,'k-')
        hold on
        errorbar(1:2,mean(y,1),ste(y,1),'.')
        figXAxis([],'Train',[0 3],1:2,{avName{iav},'Vis+Aud'})
        figYAxis([],'Frac. Correct',[0 1])
        figAxForm
        hline(0.5,'k-')
        hline(pctCorrThresh,'k:')
        title(sprintf('No Attn; Test %s, %s Model,p=%s',avName{iav},modelName{imod},...
            num2str(round(p,2,'significant'))))
    end
end
print([fnout 'comboModelPerf'],'-dpdf','-fillpage')


figure
for imod = 1:2
    subplot(2,2,imod)
    ind =  goodPerfInd_av_attn;
%     ind = mdl_attn(visualTrials).mdl(imod).pctCorr_all > pctCorrThresh & ...
%         mdl_attn(auditoryTrials).mdl(imod).pctCorr_all > pctCorrThresh;
%     maxWeight = cellfun(@(x,y,z) max(abs(cat(1,x,y,z))),...
%         mdl_attn(visualTrials).mdl(imod).weights,...
%         mdl_attn(auditoryTrials).mdl(imod).weights,...
%         mdl_attn(visualTrials).mdl(imod).comboWeight,'unif',0);
%     y = [cell2mat(cellfun(@(x,m) abs(x./m),...
%         mdl_attn(visualTrials).mdl(imod).weights(ind),maxWeight(ind),'unif',0)'),...
%         cell2mat(cellfun(@(x,m) abs(x./m),...
%         mdl_attn(auditoryTrials).mdl(imod).weights(ind),maxWeight(ind),'unif',0)'),...
%         cell2mat(cellfun(@(x,m) abs(x./m),...
%         mdl_attn(visualTrials).mdl(imod).comboWeight(ind),maxWeight(ind),'unif',0)')];
    y = abs([cell2mat(mdl_attn(visualTrials).mdl(imod).weights(ind)'),...
        cell2mat(mdl_attn(auditoryTrials).mdl(imod).weights(ind)'),...
        cell2mat(mdl_attn(visualTrials).mdl(imod).comboWeight(ind)')]);
    yNorm = y./max(y,[],2);
    [p,~,stats] = anova1(yNorm,[],'off');
    fprintf('Attn %s Model Max Weight = %s\n',modelName{imod},sigfigString(max(y(:))))
    tbl = multcompare(stats,'display','off');
    fprintf('Attn %s Model:\n',modelName{imod})
    disp(tbl(:,[1:2,6]))
    hold on
    errorbar(1:3,mean(yNorm,1),ste(yNorm,1),'.')
    figXAxis([],'',[0 4],1:3,cat(2,avName,{'Vis+Aud'}))
%     figYAxis([],'Abs Weight (norm across conditions)',[0 0.5])
    figAxForm
    title(sprintf('Attn %s Model,p=%s',modelName{imod},sigfigString(p)))

    subplot(2,2,imod+2)
        ind = goodPerfInd_av_noAttn;
%     ind = mdl_noAttn(visualTrials).mdl(imod).pctCorr_all > pctCorrThresh & ...
%         mdl_noAttn(auditoryTrials).mdl(imod).pctCorr_all > pctCorrThresh;
%     maxWeight = cellfun(@(x,y,z) max(abs(cat(1,x,y,z))),...
%         mdl_noAttn(visualTrials).mdl(imod).weights,...
%         mdl_noAttn(auditoryTrials).mdl(imod).weights,...
%         mdl_noAttn(visualTrials).mdl(imod).comboWeight,'unif',0);
%     y = [cell2mat(cellfun(@(x,m) abs(x./m),...
%         mdl_noAttn(visualTrials).mdl(imod).weights(ind),maxWeight(ind),'unif',0)'),...
%         cell2mat(cellfun(@(x,m) abs(x./m),...
%         mdl_noAttn(auditoryTrials).mdl(imod).weights(ind),maxWeight(ind),'unif',0)'),...
%         cell2mat(cellfun(@(x,m) abs(x./m),...
%         mdl_noAttn(visualTrials).mdl(imod).comboWeight(ind),maxWeight(ind),'unif',0)')];
    y = abs([cell2mat(mdl_noAttn(visualTrials).mdl(imod).weights(ind)'),...
        cell2mat(mdl_noAttn(auditoryTrials).mdl(imod).weights(ind)'),...
        cell2mat(mdl_noAttn(visualTrials).mdl(imod).comboWeight(ind)')]);
    yNorm = y./max(y,[],2);
    fprintf('No Attn %s Model Max Weight = %s\n',modelName{imod},sigfigString(max(y(:))))
    [p,~,stats] = anova1(yNorm,[],'off');
    tbl = multcompare(stats,'display','off');
    fprintf('No Attn %s Model:\n',modelName{imod})
    disp(tbl(:,[1:2,6]))
    hold on
    errorbar(1:3,mean(yNorm,1),ste(yNorm,1),'.')
    figXAxis([],'',[0 4],1:3,cat(2,avName,{'Vis+Aud'}))
%     figYAxis([],'Abs Weight (norm across conditions)',[0 0.5])
    figAxForm
    title(sprintf('No Attn %s Model,p=%s',modelName{imod},sigfigString(p)))
end

print([fnout 'comboModelWeights'],'-dpdf','-fillpage')

%% shuffle model performance
for iav = 1:2
    mdl_attn(iav).mdl(stimModInd).pctCorr_shuff = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorr_shuff = nan(1,nexp_attn);
    mdl_attn(iav).mdl(stimModInd).pctCorr_shuff_otherAV = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorr_shuff_otherAV = nan(1,nexp_attn);   
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).pctCorr_shuff(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectTarget_shuff_holdout;        
        mdl_attn(iav).mdl(choiceModInd).pctCorr_shuff(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_shuff_holdout;
        mdl_attn(iav).mdl(stimModInd).pctCorr_shuff_otherAV(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectTarget_shuff_otherAV;
        mdl_attn(iav).mdl(choiceModInd).pctCorr_shuff_otherAV(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_shuff_otherAV;
    end
    mdl_noAttn(iav).mdl(stimModInd).pctCorr_shuff = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_shuff = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(stimModInd).pctCorr_shuff_otherAV = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_shuff_otherAV = nan(1,nexp_noAttn);   
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(stimModInd).pctCorr_shuff(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectTarget_shuff_holdout;        
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_shuff(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_shuff_holdout;
        mdl_noAttn(iav).mdl(stimModInd).pctCorr_shuff_otherAV(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectTarget_shuff_otherAV;
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_shuff_otherAV(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_shuff_otherAV;
    end
end

attn_fig = figure;
suptitle('Attn Mice')
noAttn_fig = figure;
suptitle('No Attn Mice')
for iav = 1:2
    if iav == 1
        otherAV = 2;
    else 
        otherAV = 1;
    end
    for imod = 1:2
        if imod == 1
            iplot = 0;
        else
            iplot = 2;
        end
        figure(attn_fig)
        ind = mdl_attn(iav).mdl(imod).pctCorr_all >pctCorrThresh;
        pc_shuff = mdl_attn(iav).mdl(imod).pctCorr_shuff(ind);
        pc = mdl_attn(iav).mdl(imod).pctCorr_all(ind);
        pc_otherShuff = mdl_attn(iav).mdl(imod).pctCorr_shuff_otherAV(ind);
        
        subplot(4,2,iav+iplot)
        plot(1:2,[pc;pc_shuff],'k-')
        hold on
        errorbar(1:2,[mean(pc),mean(pc_shuff)],[ste(pc,2),ste(pc_shuff,2)],'.');
        [~,p] = ttest(pc,pc_shuff);
        figXAxis([],'trial correlations',[0 3],1:2,{'intact','shuff'})
        figYAxis([],'Fraction Correct',[0 1])
        figAxForm
        title(sprintf('%s %s Model, p=%s',avName{iav},modelName{imod},sigfigString(p)))
        
        subplot(4,2,iav+iplot+4)
        plot(1:2,[pc_shuff;pc_otherShuff],'k-')
        hold on
        errorbar(1:2,[mean(pc_shuff),mean(pc_otherShuff)],...
            [ste(pc_shuff,2),ste(pc_otherShuff,2)],'.');
        [~,p] = ttest(pc_shuff,pc_otherShuff);
        figXAxis([],'test data',[0 3],1:2,{avName{iav},avName{otherAV}})
        figYAxis([],'Fraction Correct',[0 1])
        figAxForm
        title(sprintf('Shuffle %s %s Model, p=%s',avName{iav},modelName{imod},sigfigString(p)))

        figure(noAttn_fig)
        ind = mdl_noAttn(iav).mdl(imod).pctCorr_all >pctCorrThresh;
        pc_shuff = mdl_noAttn(iav).mdl(imod).pctCorr_shuff(ind);
        pc = mdl_noAttn(iav).mdl(imod).pctCorr_all(ind);
        pc_otherShuff = mdl_noAttn(iav).mdl(imod).pctCorr_shuff_otherAV(ind);
        
        subplot(4,2,iav+iplot)
        plot(1:2,[pc;pc_shuff],'k-')
        hold on
        errorbar(1:2,[mean(pc),mean(pc_shuff)],[ste(pc,2),ste(pc_shuff,2)],'.');
        [~,p] = ttest(pc,pc_shuff);
        figXAxis([],'trial correlations',[0 3],1:2,{'intact','shuff'})
        figYAxis([],'Fraction Correct',[0 1])
        figAxForm
        title(sprintf('%s %s Model, p=%s',avName{iav},modelName{imod},sigfigString(p)))
        
        subplot(4,2,iav+iplot+4)
        plot(1:2,[pc_shuff;pc_otherShuff],'k-')
        hold on
        errorbar(1:2,[mean(pc_shuff),mean(pc_otherShuff)],...
            [ste(pc_shuff,2),ste(pc_otherShuff,2)],'.');
        [~,p] = ttest(pc_shuff,pc_otherShuff);
        figXAxis([],'test data',[0 3],1:2,{avName{iav},avName{otherAV}})
        figYAxis([],'Fraction Correct',[0 1])
        figAxForm
        title(sprintf('Shuffle %s %s Model, p=%s',avName{iav},modelName{imod},sigfigString(p)))
        
    end
end
figure(attn_fig)
print([fnout 'shufflePctCorrect_attn'],'-dpdf','-fillpage')
figure(noAttn_fig)
print([fnout 'shufflePctCorrect_noAttn'],'-dpdf','-fillpage')

%% weight analysis of stimulus-inclusive dc model
for iav = 1:2
    mdl_attn(iav).mdl(stimModInd).pctCorrect_withChoice = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorrect_withStim = nan(1,nexp_attn);
    mdl_attn(iav).mdl(stimModInd).modeledChoiceWeight = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).modeledStimWeight = nan(1,nexp_attn);
    mdl_attn(iav).mdl(stimModInd).weights_withChoice = cell(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).weights_withStim = cell(1,nexp_attn);
    mdl_attn(iav).mdl(stimModInd).weights_withFixedChoice = cell(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).weights_withFixedStim = cell(1,nexp_attn);
    
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).pctCorrect_withChoice(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectTarget_withChoice_holdout;
        mdl_attn(iav).mdl(choiceModInd).pctCorrect_withStim(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_withStim_holdout;
        mdl_attn(iav).mdl(stimModInd).modeledChoiceWeight(iexp) = ...
            dc_attn(iexp).av(iav).modeledChoiceWeight;
        mdl_attn(iav).mdl(choiceModInd).modeledStimWeight(iexp) = ...
            dc_attn(iexp).av(iav).modeledStimWeight;
        mdl_attn(iav).mdl(stimModInd).weights_withChoice{iexp} = ...
            dc_attn(iexp).av(iav).weightTarget_withChoice;
        mdl_attn(iav).mdl(choiceModInd).weights_withStim{iexp} = ...
            dc_attn(iexp).av(iav).weightDetect_withStim; 
        
        mdl_attn(iav).mdl(stimModInd).weights_withFixedChoice{iexp} = ...
            dc_attn(iexp).av(iav).weightTarget_withFixedChoice;
        mdl_attn(iav).mdl(choiceModInd).weights_withFixedStim{iexp} = ...
            dc_attn(iexp).av(iav).weightDetect_withFixedStim; 
    end
    
    mdl_noAttn(iav).mdl(stimModInd).pctCorrect_withChoice = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_withStim = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(stimModInd).modeledChoiceWeight = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).modeledStimWeight = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(stimModInd).weights_withChoice = cell(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).weights_withStim = cell(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(stimModInd).weights_withFixedChoice = cell(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).weights_withFixedStim = cell(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(stimModInd).pctCorrect_withChoice(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectTarget_withChoice_holdout;
        mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_withStim(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_withStim_holdout;
        mdl_noAttn(iav).mdl(stimModInd).modeledChoiceWeight(iexp) = ...
            dc_noAttn(iexp).av(iav).modeledChoiceWeight;
        mdl_noAttn(iav).mdl(choiceModInd).modeledStimWeight(iexp) = ...
            dc_noAttn(iexp).av(iav).modeledStimWeight;
        mdl_noAttn(iav).mdl(stimModInd).weights_withChoice{iexp} = ...
            dc_noAttn(iexp).av(iav).weightTarget_withChoice;
        mdl_noAttn(iav).mdl(choiceModInd).weights_withStim{iexp} = ...
            dc_noAttn(iexp).av(iav).weightDetect_withStim;     
        mdl_noAttn(iav).mdl(stimModInd).weights_withFixedChoice{iexp} = ...
            dc_noAttn(iexp).av(iav).weightTarget_withFixedChoice;
        mdl_noAttn(iav).mdl(choiceModInd).weights_withFixedStim{iexp} = ...
            dc_noAttn(iexp).av(iav).weightDetect_withFixedStim;    
    end
end
weightTooHigh_fixed_attn = cellfun(@(x) max(abs(x)),...
    mdl_attn(visualTrials).mdl(stimModInd).weights_withFixedChoice)> 30 |...
    cellfun(@(x) max(abs(x)),...
    mdl_attn(visualTrials).mdl(stimModInd).weights_withFixedChoice)> 30 |...
    cellfun(@(x) max(abs(x)),...
    mdl_attn(visualTrials).mdl(choiceModInd).weights_withFixedStim)> 30 |...
    cellfun(@(x) max(abs(x)),...
    mdl_attn(auditoryTrials).mdl(choiceModInd).weights_withFixedStim) > 30 | ...
    mdl_attn(visualTrials).mdl(stimModInd).modeledChoiceWeight > 30 | ...
    mdl_attn(visualTrials).mdl(choiceModInd).modeledStimWeight > 30 | ...
    mdl_attn(auditoryTrials).mdl(stimModInd).modeledChoiceWeight > 30 | ...
    mdl_attn(auditoryTrials).mdl(choiceModInd).modeledStimWeight > 30;
weightTooHigh_fixed_noAttn = cellfun(@(x) max(abs(x)),...
    mdl_noAttn(visualTrials).mdl(stimModInd).weights_withFixedChoice)> 30 |...
    cellfun(@(x) max(abs(x)),...
    mdl_noAttn(visualTrials).mdl(stimModInd).weights_withFixedChoice)> 30 |...
    cellfun(@(x) max(abs(x)),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).weights_withFixedStim)> 30 |...
    cellfun(@(x) max(abs(x)),...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_withFixedStim) > 30 | ...
    mdl_noAttn(visualTrials).mdl(stimModInd).modeledChoiceWeight > 30 | ...
    mdl_noAttn(visualTrials).mdl(choiceModInd).modeledStimWeight > 30 | ...
    mdl_noAttn(auditoryTrials).mdl(stimModInd).modeledChoiceWeight > 30 | ...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).modeledStimWeight > 30;
goodPerfInd_av_fixed_attn = goodPerfInd_av_attn & ~weightTooHigh_fixed_attn;
goodPerfInd_av_fixed_noAttn = goodPerfInd_av_noAttn & ~weightTooHigh_fixed_noAttn;

for iav = 1:2
    rng(0)
    if iav == 1
        for ishuff = 1:nShuff
%                 
            % shuffle weights
            shuffWeights_vis = cellfun(@(x) x(randperm(length(x))),...
                mdl_attn(visualTrials).mdl(choiceModInd).weights_withFixedStim(goodPerfInd_av_fixed_attn),'unif',0)';
            shuffWeights_aud = cellfun(@(x) x(randperm(length(x))),...
                mdl_attn(auditoryTrials).mdl(choiceModInd).weights_withFixedStim(goodPerfInd_av_fixed_attn),'unif',0)';
            if ishuff == 1
                mdl_attn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_withFixedStim = ...
                    nan(sum(cellfun(@length,shuffWeights_vis)),nShuff);
            end
            mdl_attn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_withFixedStim(:,ishuff) = ...
                cell2mat(cellfun(@(x,y) ...
                rad2deg(atan(x./y)),shuffWeights_aud,...
                shuffWeights_vis,'unif',0));
            
            shuffWeights_vis = cellfun(@(x) x(randperm(length(x))),...
                mdl_attn(visualTrials).mdl(stimModInd).weights_withFixedChoice(goodPerfInd_av_fixed_attn),'unif',0)';
            shuffWeights_aud = cellfun(@(x) x(randperm(length(x))),...
                mdl_attn(auditoryTrials).mdl(stimModInd).weights_withFixedChoice(goodPerfInd_av_fixed_attn),'unif',0)';
            if ishuff == 1
                mdl_attn(visualTrials).mdl(stimModInd).weightAVangle_shuff_withFixedChoice = ...
                    nan(sum(cellfun(@length,shuffWeights_vis)),nShuff);
            end
            mdl_attn(visualTrials).mdl(stimModInd).weightAVangle_shuff_withFixedChoice(:,ishuff) = ...
                cell2mat(cellfun(@(x,y) ...
                rad2deg(atan(x./y)),shuffWeights_aud,...
                shuffWeights_vis,'unif',0));

            shuffWeights_vis = cellfun(@(x) x(randperm(length(x))),...
                mdl_noAttn(visualTrials).mdl(choiceModInd).weights_withFixedStim(goodPerfInd_av_fixed_noAttn),'unif',0)';
            shuffWeights_aud = cellfun(@(x) x(randperm(length(x))),...
                mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_withFixedStim(goodPerfInd_av_fixed_noAttn),'unif',0)';
            if ishuff == 1
                mdl_noAttn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_withFixedStim = ...
                    nan(sum(cellfun(@length,shuffWeights_vis)),nShuff);
            end
            mdl_noAttn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_withFixedStim(:,ishuff) = ...
                cell2mat(cellfun(@(x,y) ...
                rad2deg(atan(x./y)),shuffWeights_aud,...
                shuffWeights_vis,'unif',0));
            
            shuffWeights_vis = cellfun(@(x) x(randperm(length(x))),...
                mdl_noAttn(visualTrials).mdl(stimModInd).weights_withFixedChoice(goodPerfInd_av_fixed_noAttn),'unif',0)';
            shuffWeights_aud = cellfun(@(x) x(randperm(length(x))),...
                mdl_noAttn(auditoryTrials).mdl(stimModInd).weights_withFixedChoice(goodPerfInd_av_fixed_noAttn),'unif',0)';
            if ishuff == 1
                mdl_noAttn(visualTrials).mdl(stimModInd).weightAVangle_shuff_withFixedChoice = ...
                    nan(sum(cellfun(@length,shuffWeights_vis)),nShuff);
            end
            mdl_noAttn(visualTrials).mdl(stimModInd).weightAVangle_shuff_withFixedChoice(:,ishuff) = ...
                cell2mat(cellfun(@(x,y) ...
                rad2deg(atan(x./y)),shuffWeights_aud,...
                shuffWeights_vis,'unif',0));
        end
    end
end

figure
suptitle('Attn Mice')
for iav = 1:2
    x = mdl_attn(iav).mdl(stimModInd).pctCorr_all(goodPerfInd_av_fixed_attn);
    y = mdl_attn(iav).mdl(stimModInd).pctCorrect_withChoice(goodPerfInd_av_fixed_attn);
    
    subplot(2,2,iav)
    hold on
    plot(x,y,'.','MarkerSize',20)
    figXAxis([],'pct corr. (pcs only)',[0.5 1])
    figYAxis([],'pct corr. (pcs + choice)',[0.5 1])
    figAxForm
    plot(0:1,0:1,'k--')
    title(sprintf('%s Stimulus Model',avName{iav}))
    
    x = mdl_attn(iav).mdl(choiceModInd).pctCorr_all(goodPerfInd_av_fixed_attn);
    y = mdl_attn(iav).mdl(choiceModInd).pctCorrect_withStim(goodPerfInd_av_fixed_attn);
    
    subplot(2,2,iav+2)
    hold on
    plot(x,y,'.','MarkerSize',20)
    figXAxis([],'pct corr. (pcs only)',[0.5 1])
    figYAxis([],'pct corr. (pcs + choice)',[0.5 1])
    figAxForm
    plot(0:1,0:1,'k--')
    title(sprintf('%s Choice Model',avName{iav}))
    
end
print([fnout 'pctCorrect_withVsWithoutStimCh_attn'],'-dpdf','-fillpage')
figure
suptitle('No Attn Mice')
for iav = 1:2
    x = mdl_noAttn(iav).mdl(stimModInd).pctCorr_all(goodPerfInd_av_fixed_noAttn);
    y = mdl_noAttn(iav).mdl(stimModInd).pctCorrect_withChoice(goodPerfInd_av_fixed_noAttn);
    
    subplot(2,2,iav)
    hold on
    plot(x,y,'.','MarkerSize',20)
    figXAxis([],'pct corr. (pcs only)',[0.5 1])
    figYAxis([],'pct corr. (pcs + choice)',[0.5 1])
    figAxForm
    plot(0:1,0:1,'k--')
    title(sprintf('%s Stimulus Model',avName{iav}))
    
    x = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all(goodPerfInd_av_fixed_noAttn);
    y = mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_withStim(goodPerfInd_av_fixed_noAttn);
    
    subplot(2,2,iav+2)
    hold on
    plot(x,y,'.','MarkerSize',20)
    figXAxis([],'pct corr. (pcs only)',[0.5 1])
    figYAxis([],'pct corr. (pcs + choice)',[0.5 1])
    figAxForm
    plot(0:1,0:1,'k--')
    title(sprintf('%s Choice Model',avName{iav}))
    
end
print([fnout 'pctCorrect_withVsWithoutStimCh_noAttn'],'-dpdf','-fillpage')

figure
for iav = 1:2
    x = mdl_attn(iav).mdl(stimModInd).modeledChoiceWeight(goodPerfInd_av_fixed_attn);
    y = mdl_noAttn(iav).mdl(stimModInd).modeledChoiceWeight(goodPerfInd_av_fixed_noAttn);
    
    subplot(2,2,iav)
    hold on
    plot(ones(length(x),1),x,'k.','MarkerSize',20)
    plot(ones(length(y),1).*2,y,'k.','MarkerSize',20)
    errorbar(1,mean(x),ste(x,2),'.')
    errorbar(2,mean(y),ste(y,2),'.')
    figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
    figYAxis([],'modeled choice weight',[-1 10])
    figAxForm
    title(sprintf('%s Stimulus Model',avName{iav}))
    
    x = mdl_attn(iav).mdl(choiceModInd).modeledStimWeight(goodPerfInd_av_fixed_attn);
    y = mdl_noAttn(iav).mdl(choiceModInd).modeledStimWeight(goodPerfInd_av_fixed_noAttn);
    
    subplot(2,2,iav+2)
    hold on
    plot(ones(length(x),1),x,'k.','MarkerSize',20)
    plot(ones(length(y),1).*2,y,'k.','MarkerSize',20)
    errorbar(1,mean(x),ste(x,2),'.')
    errorbar(2,mean(y),ste(y,2),'.')
    figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
    figYAxis([],'modeled stimulus weight',[-1 10])
    figAxForm
    title(sprintf('%s Choice Model',avName{iav}))
end
print([fnout 'compareStimChoiceWeights'],'-dpdf','-fillpage')

figure
x = mdl_attn(visualTrials).mdl(stimModInd).modeledChoiceWeight(goodPerfInd_av_fixed_attn);
y = mdl_attn(auditoryTrials).mdl(stimModInd).modeledChoiceWeight(goodPerfInd_av_fixed_attn);
subplot 221
hold on
plot(x,y,'.','MarkerSize',20)
plot(-20:20,-20:20,'k--')
figXAxis([],'visual model choice weight',[0 10])
figYAxis([],'auditory model choice weight',[0 10])
figAxForm
title('Attn Mice, Stimulus Model')
x = mdl_attn(visualTrials).mdl(choiceModInd).modeledStimWeight(goodPerfInd_av_fixed_attn);
y = mdl_attn(auditoryTrials).mdl(choiceModInd).modeledStimWeight(goodPerfInd_av_fixed_attn);
subplot 222
hold on
plot(x,y,'.','MarkerSize',20)
plot(-20:20,-20:20,'k--')
figXAxis([],'visual model stim weight',[0 10])
figYAxis([],'auditory model stim weight',[0 10])
figAxForm
title('Attn Mice, Choice Model')

x = mdl_noAttn(visualTrials).mdl(stimModInd).modeledChoiceWeight(goodPerfInd_av_fixed_noAttn);
y = mdl_noAttn(auditoryTrials).mdl(stimModInd).modeledChoiceWeight(goodPerfInd_av_fixed_noAttn);
subplot 223
hold on
plot(x,y,'.','MarkerSize',20)
plot(-20:20,-20:20,'k--')
figXAxis([],'visual model choice weight',[0 10])
figYAxis([],'auditory model choice weight',[0 10])
figAxForm
title('No Attn Mice, Stimulus Model')
x = mdl_noAttn(visualTrials).mdl(choiceModInd).modeledStimWeight(goodPerfInd_av_fixed_noAttn);
y = mdl_noAttn(auditoryTrials).mdl(choiceModInd).modeledStimWeight(goodPerfInd_av_fixed_noAttn);
subplot 224
hold on
plot(x,y,'.','MarkerSize',20)
plot(-20:20,-20:20,'k--')
figXAxis([],'visual model stim weight',[0 10])
figYAxis([],'auditory model stim weight',[0 10])
figAxForm
title('No Attn Mice, Choice Model')

print([fnout 'compareVisVsAudStimChoiceWeights'],'-dpdf','-fillpage')
%%
figure
suptitle('Stim or Choice Info Model')
for iav = 1:2
    subplot(2,2,iav)
    y1 = cell2mat(mdl_attn(iav).mdl(stimModInd).weights_withFixedChoice(goodPerfInd_av_fixed_attn)');
    y2 = cell2mat(mdl_noAttn(iav).mdl(stimModInd).weights_withFixedChoice(goodPerfInd_av_fixed_noAttn)');
    cdfplot(y1)
    hold on
    cdfplot(y2)
    figXAxis([],'Weight',[-1 1])
    figYAxis([],'Fraction of Cells',[0 1])
    figAxForm
    title(sprintf('%s %s Model',...
        avName{iav},modelName{stimModInd}))
    vline(0,'k:')
    legend({'Attn','No Attn'},'location','southeast')
    
    subplot(2,2,iav+2)
    y1 = cell2mat(mdl_attn(iav).mdl(choiceModInd).weights_withFixedStim(goodPerfInd_av_fixed_attn)');
    y2 = cell2mat(mdl_noAttn(iav).mdl(choiceModInd).weights_withFixedStim(goodPerfInd_av_fixed_noAttn)');
    cdfplot(y1)
    hold on
    cdfplot(y2)
    figXAxis([],'Weight',[-1 1])
    figYAxis([],'Fraction of Cells',[0 1])
    figAxForm
    title(sprintf('%s %s Model',...
        avName{iav},modelName{choiceModInd}))
    vline(0,'k:')
    legend({'Attn','No Attn'},'location','southeast')
end
print([fnout 'cdfStimChoiceWeights'],'-dpdf','-fillpage')

mdlWeightsSub_fig = figure;
mdlWeights_fig = figure;
suptitle({'Choice Model with Stim Info, Weights';...
    sprintf('pct corr must be > %s in vis & aud stim models',...
    num2str(pctCorrThresh))})

figure(mdlWeights_fig)
subplot(3,2,1)
y1 = cell2mat(mdl_attn(visualTrials).mdl(choiceModInd).weights_withFixedStim(goodPerfInd_av_fixed_attn)');
y2 = cell2mat(mdl_attn(auditoryTrials).mdl(choiceModInd).weights_withFixedStim(goodPerfInd_av_fixed_attn)');
plot(y1,y2,'k.','MarkerSize',10)
hold on
mdl_fitline = fitlm(y1,y2);
rsq = mdl_fitline.Rsquared.Ordinary;
yfit = predict(mdl_fitline,y1);
plot(y1,yfit,'-')
figXAxis([],'Visual Weight',[weightEdges(1) weightEdges(end)])
figYAxis([],'Auditory Weight',[weightEdges(1) weightEdges(end)])
figAxForm
hline(0,'k:')
vline(0,'k:')
r = corr(y1,y2);
title({'Attn Mice';sprintf('%s Model, r=%s',...
    modelName{choiceModInd},num2str(round(r,2,'significant')))})


subplot(3,2,3)
weightAVangle = rad2deg(atan(y2./y1));
histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
figXAxis([],'Angle of A/V Weight',[-100 100],oriEdges_noAbs)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
s = skewness(weightAVangle);
title({'Attn Mice';sprintf('%s Model, skew=%s,ks=%s',...
    modelName{choiceModInd},num2str(round(s,2,'significant')))})

subplot(3,2,5)
avAngleShuff = mdl_attn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_withFixedStim;
h = nan(length(oriEdges_noAbs)-1,nShuff);
for ishuff = 1:nShuff
    h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges_noAbs,'Normalization','Probability');
end    
h_sort = sort(h,2);
y = mean(h_sort,2);
yerrl = y-h_sort(:,round(nShuff*0.025));
yerru = h_sort(:,round(nShuff*0.975))-y;
shadedErrorBar_chooseColor(oriEdges_noAbs(1:end-1),y,[yerrl';yerru'],[0 0 0]);
hold on
plot(oriEdges_noAbs(1:end-1),histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability'),'.-')
controlSub_attn = histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability') - y';
figXAxis([],'Angle of A/V Weight',[-100 100],oriEdges_noAbs)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
s = skewness(y);
[~,p] = kstest2(h_sort(:),weightAVangle);
title({'Attn Mice-Shuff';sprintf('%s Model, skew=%s,p=%s',...
    modelName{choiceModInd},num2str(round(s,2,'significant')),...
    sigfigString(p))})

figure(mdlWeightsSub_fig); subplot 122
plot(oriEdges_noAbs(1:end-1),controlSub_attn,'.-')
%     shadedErrorBar_chooseColor(oriEdges(1:end-1),controlSub_attn,[yerrl';yerru'],[0 0 0]);

figure(mdlWeights_fig)
subplot(3,2,2)
y1 = cell2mat(mdl_noAttn(visualTrials).mdl(choiceModInd).weights_withFixedStim(goodPerfInd_av_fixed_noAttn)');
y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_withFixedStim(goodPerfInd_av_fixed_noAttn)');
plot(y1,y2,'k.','MarkerSize',10)
hold on
mdl_fitline = fitlm(y1,y2);
rsq = mdl_fitline.Rsquared.Ordinary;
yfit = predict(mdl_fitline,y1);
plot(y1,yfit,'-')
figXAxis([],'Visual Weight',[weightEdges(1) weightEdges(end)])
figYAxis([],'Auditory Weight',[weightEdges(1) weightEdges(end)])
figAxForm
hline(0,'k:')
vline(0,'k:')
r = corr(y1,y2);
title({'No Attn Mice';sprintf('%s Model, r=%s',...
    modelName{choiceModInd},num2str(round(r,2,'significant')))})

subplot(3,2,4)
weightAVangle = rad2deg(atan(y2./y1));
histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
figXAxis([],'Angle of A/V Weight',[-100 100],oriEdges_noAbs)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
s = skewness(weightAVangle);
title({'No Attn Mice';sprintf('%s Model, skew=%s',...
    modelName{choiceModInd},num2str(round(s,2,'significant')))})

subplot(3,2,6)
%     y = mean(mdl_noAttn(visualTrials).mdl(choiceModInd).normWeightAVangle_shuff,2);
%     histogram(y,oriEdges,'Normalization','probability')
avAngleShuff = mdl_noAttn(visualTrials).mdl(choiceModInd).weightAVangle_shuff_withFixedStim;
h = nan(length(oriEdges_noAbs)-1,nShuff);
for ishuff = 1:nShuff
    h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges_noAbs,'Normalization','Probability');
end    
h_sort = sort(h,2);
y = mean(h_sort,2);
yerrl = y-h_sort(:,round(nShuff*0.025));
yerru = h_sort(:,round(nShuff*0.975))-y;
shadedErrorBar_chooseColor(oriEdges_noAbs(1:end-1),y,[yerrl';yerru'],[0 0 0]);
hold on
plot(oriEdges_noAbs(1:end-1),histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability'),'.-')
controlSub_noAttn = histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability') - y';
figXAxis([],'Angle of A/V Weight',[-100 100],oriEdges_noAbs)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
s = skewness(y);
title({'No Attn Mice-Shuff';sprintf('%s Model, skew=%s',...
    modelName{choiceModInd},num2str(round(s,2,'significant')))})
print([fnout 'neuronWeights_stimInclModel'],'-dpdf','-fillpage')

figure(mdlWeightsSub_fig); subplot 122
hold on
%     shadedErrorBar_chooseColor(oriEdges(1:end-1),controlSub_noAttn,[yerrl';yerru'],[0.5 0.5 0.5]);
plot(oriEdges_noAbs(1:end-1),controlSub_noAttn,'.-')
figXAxis([],'Angle of A/V Weight',[-100 100],oriEdges_noAbs)
figYAxis([],'Control Subtr. Fraction of Cells',[-0.15 0.15])
figAxForm
hline(0,'k:')
title(sprintf('%s Model, Model Including Stim',modelName{choiceModInd}))

%%
mdlWeights_fig = figure;
suptitle({'Stim Model with Choice Info, Weights';...
    sprintf('pct corr must be > %s in vis & aud stim models',...
    num2str(pctCorrThresh))})

figure(mdlWeights_fig)
subplot(3,2,1)
y1 = cell2mat(mdl_attn(visualTrials).mdl(stimModInd).weights_withFixedChoice(goodPerfInd_av_fixed_attn)');
y2 = cell2mat(mdl_attn(auditoryTrials).mdl(stimModInd).weights_withFixedChoice(goodPerfInd_av_fixed_attn)');
plot(y1,y2,'k.','MarkerSize',10)
hold on
mdl_fitline = fitlm(y1,y2);
rsq = mdl_fitline.Rsquared.Ordinary;
yfit = predict(mdl_fitline,y1);
plot(y1,yfit,'-')
figXAxis([],'Visual Weight',[weightEdges(1) weightEdges(end)])
figYAxis([],'Auditory Weight',[weightEdges(1) weightEdges(end)])
figAxForm
hline(0,'k:')
vline(0,'k:')
r = corr(y1,y2);
title({'Attn Mice';sprintf('%s Model, r=%s',...
    modelName{stimModInd},num2str(round(r,2,'significant')))})


subplot(3,2,3)
weightAVangle = rad2deg(atan(y2./y1));
histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
figXAxis([],'Angle of A/V Weight',[-100 100],oriEdges_noAbs)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
s = skewness(weightAVangle);
title({'Attn Mice';sprintf('%s Model, skew=%s,ks=%s',...
    modelName{stimModInd},num2str(round(s,2,'significant')))})

subplot(3,2,5)
avAngleShuff = mdl_attn(visualTrials).mdl(stimModInd).weightAVangle_shuff_withFixedChoice;
h = nan(length(oriEdges_noAbs)-1,nShuff);
for ishuff = 1:nShuff
    h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges_noAbs,'Normalization','Probability');
end    
h_sort = sort(h,2);
y = mean(h_sort,2);
yerrl = y-h_sort(:,round(nShuff*0.025));
yerru = h_sort(:,round(nShuff*0.975))-y;
shadedErrorBar_chooseColor(oriEdges_noAbs(1:end-1),y,[yerrl';yerru'],[0 0 0]);
hold on
plot(oriEdges_noAbs(1:end-1),histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability'),'.-')
controlSub_attn = histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability') - y';
figXAxis([],'Angle of A/V Weight',[-100 100],oriEdges_noAbs)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
s = skewness(y);
[~,p] = kstest2(h_sort(:),weightAVangle);
title({'Attn Mice-Shuff';sprintf('%s Model, skew=%s,p=%s',...
    modelName{stimModInd},num2str(round(s,2,'significant')),...
    sigfigString(p))})

figure(mdlWeightsSub_fig); subplot 121
plot(oriEdges_noAbs(1:end-1),controlSub_attn,'.-')
%     shadedErrorBar_chooseColor(oriEdges(1:end-1),controlSub_attn,[yerrl';yerru'],[0 0 0]);

figure(mdlWeights_fig)
subplot(3,2,2)
y1 = cell2mat(mdl_noAttn(visualTrials).mdl(stimModInd).weights_withFixedChoice(goodPerfInd_av_fixed_noAttn)');
y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(stimModInd).weights_withFixedChoice(goodPerfInd_av_fixed_noAttn)');
plot(y1,y2,'k.','MarkerSize',10)
hold on
mdl_fitline = fitlm(y1,y2);
rsq = mdl_fitline.Rsquared.Ordinary;
yfit = predict(mdl_fitline,y1);
plot(y1,yfit,'-')
figXAxis([],'Visual Weight',[weightEdges(1) weightEdges(end)])
figYAxis([],'Auditory Weight',[weightEdges(1) weightEdges(end)])
figAxForm
hline(0,'k:')
vline(0,'k:')
r = corr(y1,y2);
title({'No Attn Mice';sprintf('%s Model, r=%s',...
    modelName{stimModInd},num2str(round(r,2,'significant')))})

subplot(3,2,4)
weightAVangle = rad2deg(atan(y2./y1));
histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
figXAxis([],'Angle of A/V Weight',[-100 100],oriEdges_noAbs)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
s = skewness(weightAVangle);
title({'No Attn Mice';sprintf('%s Model, skew=%s',...
    modelName{stimModInd},num2str(round(s,2,'significant')))})

subplot(3,2,6)
avAngleShuff = mdl_noAttn(visualTrials).mdl(stimModInd).weightAVangle_shuff_withFixedChoice;
h = nan(length(oriEdges_noAbs)-1,nShuff);
for ishuff = 1:nShuff
    h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges_noAbs,'Normalization','Probability');
end    
h_sort = sort(h,2);
y = mean(h_sort,2);
yerrl = y-h_sort(:,round(nShuff*0.025));
yerru = h_sort(:,round(nShuff*0.975))-y;
shadedErrorBar_chooseColor(oriEdges_noAbs(1:end-1),y,[yerrl';yerru'],[0 0 0]);
hold on
plot(oriEdges_noAbs(1:end-1),histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability'),'.-')
controlSub_noAttn = histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability') - y';
figXAxis([],'Angle of A/V Weight',[-100 100],oriEdges_noAbs)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
s = skewness(y);
title({'No Attn Mice-Shuff';sprintf('%s Model, skew=%s',...
    modelName{stimModInd},num2str(round(s,2,'significant')))})
print([fnout 'neuronWeights_choiceInclModel'],'-dpdf','-fillpage')

figure(mdlWeightsSub_fig); subplot 121
hold on
plot(oriEdges_noAbs(1:end-1),controlSub_noAttn,'.-')
figXAxis([],'Angle of A/V Weight',[-100 100],oriEdges_noAbs)
figYAxis([],'Control Subtr. Fraction of Cells',[-0.15 0.15])
figAxForm
hline(0,'k:')
title(sprintf('%s Model, Model Including Stim',modelName{stimModInd}))
print([fnout 'neuronWeightAVangle_stimChoiceInclModel_ciSub'],'-dpdf','-fillpage')

%% compare A and V weights

VsubAIndex_attn = {cell2mat(cellfun(@(x,y) (abs(x)-abs(y))./(abs(x)+abs(y)),...
    mdl_attn(visualTrials).mdl(stimModInd).weights_withFixedChoice,...
    mdl_attn(auditoryTrials).mdl(stimModInd).weights_withFixedChoice,'unif',0)');...
    cell2mat(cellfun(@(x,y) (abs(x)-abs(y))./(abs(x)+abs(y)),...
    mdl_attn(visualTrials).mdl(choiceModInd).weights_withFixedStim,...
    mdl_attn(auditoryTrials).mdl(choiceModInd).weights_withFixedStim,'unif',0)')};

VsubAIndex_noAttn = {cell2mat(cellfun(@(x,y) (abs(x)-abs(y))./(abs(x)+abs(y)),...
    mdl_noAttn(visualTrials).mdl(stimModInd).weights_withFixedChoice,...
    mdl_noAttn(auditoryTrials).mdl(stimModInd).weights_withFixedChoice,'unif',0)');...
    cell2mat(cellfun(@(x,y) (abs(x)-abs(y))./(abs(x)+abs(y)),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).weights_withFixedStim,...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_withFixedStim,'unif',0)')};

binEdges_vsuba = [-1:0.1:1];
figure
suptitle('all weights are absolute values')
for imod = 1:2
    subplot(2,2,imod)
    x = VsubAIndex_attn{imod};
    histogram(x,binEdges_vsuba,'Normalization','probability')
    figXAxis([],'(V-A)/(V+A)',[-1 1])
    figYAxis([],'Fraction of Neurons',[0 0.2])
    figAxForm([],0)
    title(sprintf('Attn %s Models',modelName{imod}))
    
    subplot(2,2,imod+2)
    x = VsubAIndex_noAttn{imod};
    histogram(x,binEdges_vsuba,'Normalization','probability')
    figXAxis([],'(V-A)/(V+A)',[-1 1])
    figYAxis([],'Fraction of Neurons',[0 0.2])
    figAxForm([],0)
    title(sprintf('No Attn %s Models',modelName{imod}))
end

print([fnout 'modelWeights_VSubAIndex_hist'],'-dpdf','-fillpage')

figure
suptitle('all weights are absolute values')
for imod = 1:2
    subplot(1,2,imod)
    hold on
    a = VsubAIndex_attn{imod};
    b = VsubAIndex_noAttn{imod};
    cdfplot(a)
    cdfplot(b)
    figXAxis([],'(V-A)/(V+A)',[-1 1])
    figYAxis([],'Fraction of Neurons',[0 1])
    figAxForm
    title(sprintf('%s Models',modelName{imod}))
    legend({'Attn','No Attn'},'location','northwest')
end

print([fnout 'modelWeights_VSubAIndex_cdf'],'-dpdf','-fillpage')

%% compare stim and choice parameter weights
weightLim = [weightEdges(1) weightEdges(end)];
figure
subplot 321
hold on
y1 = cell2mat(mdl_attn(visualTrials).mdl(choiceModInd).weights_withFixedStim(goodPerfInd_av_attn)');
y2 = cell2mat(mdl_attn(auditoryTrials).mdl(choiceModInd).weights_withFixedStim(goodPerfInd_av_attn)');
weightAVangle = rad2deg(atan(y2./y1));

ind1 = weightAVangle < -75 | weightAVangle > 75;
plot(y1(ind1),y2(ind1),'.','MarkerSize',10)
ind2 = weightAVangle > -15 & weightAVangle < 15;
plot(y1(ind2),y2(ind2),'.','MarkerSize',10)
ind3 = weightAVangle > 30 & weightAVangle < 60;
plot(y1(ind3),y2(ind3),'.','MarkerSize',10)
avAngle_groupCounts = [sum(ind1),sum(ind2),sum(ind3),sum(~(ind1|ind2|ind3))]./length(ind1);

plot(y1(~(ind1|ind2|ind3)),y2(~(ind1|ind2|ind3)),'k.','MarkerSize',10)
hold on
figXAxis([],'Visual Weight',weightLim)
figYAxis([],'Auditory Weight',weightLim)
figAxForm
hline(0,'k:')
vline(0,'k:')
r = corr(y1,y2);
title('Attn Choice Models')

subplot 325
h=bar(avAngle_groupCounts);
h.BarWidth = 1;
figXAxis([],'Cell Type',[0 5],1:4,{'pref aud','pref vis','vis=aud','other'})
ax = gca;
ax.XTickLabelRotation = -45;
figYAxis([],'Fraction of Cells in Group',[0 .5])
figAxForm
title('Attn Choice Models')

subplot 323
histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
%     y = histcounts(weightAVangle,oriEdges);
%     bar(oriEdges(1:end-1),y./sum(y),1);
figXAxis([],'Angle of A/V Weight',[-100 100],oriEdges_noAbs)
figYAxis([],'Fraction of Cells',[0 .3])
figAxForm
vline([-75,-15,15,30,60,75],'k--')
title('Attn Choice Models')

subplot 322
hold on
y1 = cell2mat(mdl_noAttn(visualTrials).mdl(choiceModInd).weights_withFixedStim(goodPerfInd_av_noAttn)');
y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_withFixedStim(goodPerfInd_av_noAttn)');
weightAVangle = rad2deg(atan(y2./y1));

ind1 = weightAVangle < -75 | weightAVangle > 75;
plot(y1(ind1),y2(ind1),'.','MarkerSize',10)
ind2 = weightAVangle > -15 & weightAVangle < 15;
plot(y1(ind2),y2(ind2),'.','MarkerSize',10)
ind3 = weightAVangle > 30 & weightAVangle < 60;
plot(y1(ind3),y2(ind3),'.','MarkerSize',10)
avAngle_groupCounts = [sum(ind1),sum(ind2),sum(ind3),sum(~(ind1|ind2|ind3))]./length(ind1);

plot(y1(~(ind1|ind2|ind3)),y2(~(ind1|ind2|ind3)),'k.','MarkerSize',10)
hold on
figXAxis([],'Visual Weight',weightLim)
figYAxis([],'Auditory Weight',weightLim)
figAxForm
hline(0,'k:')
vline(0,'k:')
r = corr(y1,y2);
title('No Attn Choice Models')


subplot 326
h=bar(avAngle_groupCounts);
h.BarWidth = 1;
figXAxis([],'Cell Type',[0 5],1:4,{'pref aud','pref vis','vis=aud','other'})
ax = gca;
ax.XTickLabelRotation = -45;
figYAxis([],'Fraction of Cells in Group',[0 .5])
figAxForm
title('No Attn Choice Models')

subplot 324
histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
%     y = histcounts(weightAVangle,oriEdges);
%     bar(oriEdges(1:end-1),y./sum(y),1);
figXAxis([],'Angle of A/V Weight',[-100 100],oriEdges_noAbs)
figYAxis([],'Fraction of Cells',[0 .3])
figAxForm
vline([-75,-15,15,30,60,75],'k--')
title('No Attn Choice Models')
print([fnout 'neuronWeightAVangle_stimIncl_ChoiceModel_angleCoded'],'-dpdf','-fillpage')

figure
subplot 321
hold on
y1 = cell2mat(mdl_attn(visualTrials).mdl(stimModInd).weights_withFixedChoice(goodPerfInd_av_attn)');
y2 = cell2mat(mdl_attn(auditoryTrials).mdl(stimModInd).weights_withFixedChoice(goodPerfInd_av_attn)');
weightAVangle = rad2deg(atan(y2./y1));

ind1 = weightAVangle < -75 | weightAVangle > 75;
plot(y1(ind1),y2(ind1),'.','MarkerSize',10)
ind2 = weightAVangle > -15 & weightAVangle < 15;
plot(y1(ind2),y2(ind2),'.','MarkerSize',10)
ind3 = weightAVangle > 30 & weightAVangle < 60;
plot(y1(ind3),y2(ind3),'.','MarkerSize',10)
avAngle_groupCounts = [sum(ind1),sum(ind2),sum(ind3),sum(~(ind1|ind2|ind3))]./length(ind1);

plot(y1(~(ind1|ind2|ind3)),y2(~(ind1|ind2|ind3)),'k.','MarkerSize',10)
hold on
figXAxis([],'Visual Weight',weightLim)
figYAxis([],'Auditory Weight',weightLim)
figAxForm
hline(0,'k:')
vline(0,'k:')
r = corr(y1,y2);
title('Attn Stimulus Models')

subplot 325
h=bar(avAngle_groupCounts);
h.BarWidth = 1;
figXAxis([],'Cell Type',[0 5],1:4,{'pref aud','pref vis','vis=aud','other'})
ax = gca;
ax.XTickLabelRotation = -45;
figYAxis([],'Fraction of Cells in Group',[0 .5])
figAxForm
title('Attn Stimulus Models')

subplot 323
histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
%     y = histcounts(weightAVangle,oriEdges);
%     bar(oriEdges(1:end-1),y./sum(y),1);
figXAxis([],'Angle of A/V Weight',[-100 100],oriEdges_noAbs)
figYAxis([],'Fraction of Cells',[0 .3])
figAxForm
vline([-75,-15,15,30,60,75],'k--')
title('Attn Stimulus Models')

subplot 322
hold on
y1 = cell2mat(mdl_noAttn(visualTrials).mdl(stimModInd).weights_withFixedChoice(goodPerfInd_av_noAttn)');
y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(stimModInd).weights_withFixedChoice(goodPerfInd_av_noAttn)');
weightAVangle = rad2deg(atan(y2./y1));

ind1 = weightAVangle < -75 | weightAVangle > 75;
plot(y1(ind1),y2(ind1),'.','MarkerSize',10)
ind2 = weightAVangle > -15 & weightAVangle < 15;
plot(y1(ind2),y2(ind2),'.','MarkerSize',10)
ind3 = weightAVangle > 30 & weightAVangle < 60;
plot(y1(ind3),y2(ind3),'.','MarkerSize',10)
avAngle_groupCounts = [sum(ind1),sum(ind2),sum(ind3),sum(~(ind1|ind2|ind3))]./length(ind1);

plot(y1(~(ind1|ind2|ind3)),y2(~(ind1|ind2|ind3)),'k.','MarkerSize',10)
hold on
figXAxis([],'Visual Weight',weightLim)
figYAxis([],'Auditory Weight',weightLim)
figAxForm
hline(0,'k:')
vline(0,'k:')
r = corr(y1,y2);
title('No Attn Stimulus Models')

subplot 326
h=bar(avAngle_groupCounts);
h.BarWidth = 1;
figXAxis([],'Cell Type',[0 5],1:4,{'pref aud','pref vis','vis=aud','other'})
ax = gca;
ax.XTickLabelRotation = -45;
figYAxis([],'Fraction of Cells in Group',[0 .5])
figAxForm
title('No Attn Stimulus Models')

subplot 324
histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
%     y = histcounts(weightAVangle,oriEdges);
%     bar(oriEdges(1:end-1),y./sum(y),1);
figXAxis([],'Angle of A/V Weight',[-100 100],oriEdges_noAbs)
figYAxis([],'Fraction of Cells',[0 .3])
figAxForm
vline([-75,-15,15,30,60,75],'k--')
title('No Attn Stimulus Models')
print([fnout 'neuronWeightAVangle_choiceIncl_StimModel_angleCoded'],'-dpdf','-fillpage')


%% test stim/choice model performances when stim or choice fixed in prediction
for iav = 1:2
    mdl_attn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorrect_fixedStimW = nan(1,nexp_attn);
    mdl_attn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW_pcsOnly = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorrect_fixedStimW_pcsOnly = nan(1,nexp_attn);
    mdl_attn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW_otherAV = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorrect_fixedStimW_otherAV = nan(1,nexp_attn);
    mdl_attn(iav).mdl(stimModInd).pctCorrect_choiceOnly = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorrect_stimOnly = nan(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectTarget_withFixedChoice_holdout;
        mdl_attn(iav).mdl(choiceModInd).pctCorrect_fixedStimW(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_withFixedStim_holdout;
        mdl_attn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW_pcsOnly(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectTarget_withFixedChoice_holdout_pcsOnly;
        mdl_attn(iav).mdl(choiceModInd).pctCorrect_fixedStimW_pcsOnly(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_withFixedStim_holdout_pcsOnly;
        mdl_attn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW_otherAV(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectTarget_withFixedChoice_otherAV;
        mdl_attn(iav).mdl(choiceModInd).pctCorrect_fixedStimW_otherAV(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_withFixedStim_otherAV;
        mdl_attn(iav).mdl(stimModInd).pctCorrect_choiceOnly(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectTarget_choiceOnly;        
        mdl_attn(iav).mdl(choiceModInd).pctCorrect_stimOnly(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_stimOnly;
    end
    
    mdl_noAttn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_fixedStimW = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW_pcsOnly = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_fixedStimW_pcsOnly = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW_otherAV = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_fixedStimW_otherAV = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(stimModInd).pctCorrect_choiceOnly = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_stimOnly = nan(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectTarget_withFixedChoice_holdout;
        mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_fixedStimW(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_withFixedStim_holdout;
        mdl_noAttn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW_pcsOnly(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectTarget_withFixedChoice_holdout_pcsOnly;
        mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_fixedStimW_pcsOnly(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_withFixedStim_holdout_pcsOnly;
        mdl_noAttn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW_otherAV(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectTarget_withFixedChoice_otherAV;
        mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_fixedStimW_otherAV(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_withFixedStim_otherAV;
        mdl_noAttn(iav).mdl(stimModInd).pctCorrect_choiceOnly(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectTarget_choiceOnly;        
        mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_stimOnly(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_stimOnly;
    end
end

setFigParams4Print('portrait')
pctCorr_fixedVsNorm = figure;
suptitle(sprintf('pct corr must be > %s in no stim/no choice model',...
    num2str(pctCorrThresh)))
pctCorr_fixedVsNorm_quant = figure;
for iav = 1:2
    figure(pctCorr_fixedVsNorm)
    subplot(4,2,iav)
    hold on
    ind = mdl_attn(iav).mdl(stimModInd).pctCorr_all > pctCorrThresh;
    y1 = mdl_attn(iav).mdl(stimModInd).pctCorr_all(ind);
    y2 = mdl_attn(iav).mdl(stimModInd).pctCorrect_withChoice(ind);
    y3 = mdl_attn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW_pcsOnly(ind);
    y4 = mdl_attn(iav).mdl(stimModInd).pctCorrect_choiceOnly(ind);
    plot(1:4,[y1',y2',y3',y4'],'k-')
    errorbar(1,mean(y1),ste(y1,2),'.','MarkerSize',20)
    errorbar(2,mean(y2),ste(y2,2),'.','MarkerSize',20)
    errorbar(3,mean(y3),ste(y3,2),'.','MarkerSize',20)
    errorbar(4,mean(y4),ste(y4,2),'.','MarkerSize',20)
    figXAxis([],'model',[0 5],1:4,{'PCs only','PCs+ch','(PCs+ch)-ch','ch only'})
    ax = gca;
    ax.XTickLabelRotation = -45;
    figYAxis([],'% Correct',[0 1],0:0.2:1)
    figAxForm
    hline(pctCorrThresh,'k:')
    hline(0.5,'k-')
    title(sprintf('Attn; %s Stimulus Model',avName{iav}))
    
    figure(pctCorr_fixedVsNorm_quant)
    subplot(2,2,iav)
    hold on
    y3 = (y2-y1)./y1*100;
    plot(ones(1,sum(ind)),y3,'k.','MarkerSize',10)
    errorbar(1,mean(y3),ste(y3,2),'.','MarkerSize',20)
    
    figure(pctCorr_fixedVsNorm)
    subplot(4,2,iav+2)
    hold on
    ind = mdl_attn(iav).mdl(choiceModInd).pctCorr_all > pctCorrThresh;
    y1 = mdl_attn(iav).mdl(choiceModInd).pctCorr_all(ind);
    y2 = mdl_attn(iav).mdl(choiceModInd).pctCorrect_withStim(ind);
    y3 = mdl_attn(iav).mdl(choiceModInd).pctCorrect_fixedStimW_pcsOnly(ind);
    y4 = mdl_attn(iav).mdl(choiceModInd).pctCorrect_stimOnly(ind);
    plot(1:4,[y1',y2',y3',y4'],'k-')
    errorbar(1,mean(y1),ste(y1,2),'.','MarkerSize',20)
    errorbar(2,mean(y2),ste(y2,2),'.','MarkerSize',20)
    errorbar(3,mean(y3),ste(y3,2),'.','MarkerSize',20)
    errorbar(4,mean(y4),ste(y4,2),'.','MarkerSize',20)
    figXAxis([],'model',[0 5],1:4,{'PCs only','PCs+s','(PCs+s)-s','s only'})
    ax = gca;
    ax.XTickLabelRotation = -45;
    figYAxis([],'% Correct',[0 1],0:0.2:1)
    figAxForm
    hline(pctCorrThresh,'k:')
    hline(0.5,'k-')
    title(sprintf('Attn; %s Choice Model',avName{iav}))

    figure(pctCorr_fixedVsNorm_quant)
    subplot(2,2,iav+2)
    hold on
    y3 = (y2-y1)./y1*100;
    plot(ones(1,sum(ind)),y3,'k.','MarkerSize',10)
    errorbar(1,mean(y3),ste(y3,2),'.','MarkerSize',20)
    
    figure(pctCorr_fixedVsNorm)
    subplot(4,2,iav+4)
    hold on
    ind = mdl_noAttn(iav).mdl(stimModInd).pctCorr_all > pctCorrThresh;
    y1 = mdl_noAttn(iav).mdl(stimModInd).pctCorr_all(ind);
    y2 = mdl_noAttn(iav).mdl(stimModInd).pctCorrect_withChoice(ind);
    y3 = mdl_noAttn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW_pcsOnly(ind);
    y4 = mdl_noAttn(iav).mdl(stimModInd).pctCorrect_choiceOnly(ind);
    plot(1:4,[y1',y2',y3',y4'],'k-')
    errorbar(1,mean(y1),ste(y1,2),'.','MarkerSize',20)
    errorbar(2,mean(y2),ste(y2,2),'.','MarkerSize',20)
    errorbar(3,mean(y3),ste(y3,2),'.','MarkerSize',20)
    errorbar(4,mean(y4),ste(y4,2),'.','MarkerSize',20)
    figXAxis([],'model',[0 5],1:4,{'PCs only','PCs+ch','(PCs+ch)-ch','ch only'})
    ax = gca;
    ax.XTickLabelRotation = -45;
    figYAxis([],'% Correct',[0 1],0:0.2:1)
    figAxForm
    hline(pctCorrThresh,'k:')
    hline(0.5,'k-')
    title(sprintf('No Attn; %s Stimulus Model',avName{iav}))
    
    figure(pctCorr_fixedVsNorm_quant)
    subplot(2,2,iav)
    hold on
    y3 = (y2-y1)./y1*100;
    plot(ones(1,sum(ind)).*2,y3,'k.','MarkerSize',10)
    errorbar(2,mean(y3),ste(y3,2),'.','MarkerSize',20)
    figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
    figYAxis([],'Pct Change from PC Only',[0 100])
    figAxForm([],0)
    title(sprintf('%s Stimulus Model',avName{iav}))
    
    figure(pctCorr_fixedVsNorm)
    subplot(4,2,iav+6)
    hold on
    ind = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all > pctCorrThresh;
    y1 = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all(ind);
    y2 = mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_withStim(ind);
    y3 = mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_fixedStimW_pcsOnly(ind);
    y4 = mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_stimOnly(ind);
    plot(1:4,[y1',y2',y3',y4'],'k-')
    errorbar(1,mean(y1),ste(y1,2),'.','MarkerSize',20)
    errorbar(2,mean(y2),ste(y2,2),'.','MarkerSize',20)
    errorbar(3,mean(y3),ste(y3,2),'.','MarkerSize',20)
    errorbar(4,mean(y4),ste(y4,2),'.','MarkerSize',20)
    figXAxis([],'model',[0 5],1:4,{'PCs only','PCs+s','(PCs+s)-s','s only'})
    ax = gca;
    ax.XTickLabelRotation = -45;
    figYAxis([],'% Correct',[0 1],0:0.2:1)
    figAxForm
    hline(pctCorrThresh,'k:')
    hline(0.5,'k-')
    title(sprintf('No Attn; %s Choice Model',avName{iav}))

    figure(pctCorr_fixedVsNorm_quant)
    subplot(2,2,iav+2)
    hold on
    y3 = (y2-y1)./y1*100;
    plot(ones(1,sum(ind)).*2,y3,'k.','MarkerSize',10)
    errorbar(2,mean(y3),ste(y3,2),'.','MarkerSize',20)
    figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
    figYAxis([],'Pct Change from PC Only',[0 100])
    figAxForm([],0)
    title(sprintf('%s Choice Model',avName{iav}))
end
figure(pctCorr_fixedVsNorm)
print([fnout 'pctCorrect_compareInclStimChoice'],'-dpdf','-fillpage')
figure(pctCorr_fixedVsNorm_quant)
print([fnout 'pctCorrect_compareInclStimChoice_quant'],'-dpdf')

%% test stim/choice models with opposite modalities
figure
suptitle('Stim or choice weight included for model, but removed for prediction')
for iav = 1:2
    if iav == 1
        otherAV = 2;
    else
        otherAV = 1;
    end
    subplot(4,2,iav)
    hold on
    ind = mdl_attn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW_pcsOnly > pctCorrThresh;
    y1 = mdl_attn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW_pcsOnly(ind);
    y2 = mdl_attn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW_otherAV(ind);
    plot(1:2,[y1',y2'],'k-')
    errorbar(1,mean(y1),ste(y1,2),'.','MarkerSize',20)
    errorbar(2,mean(y2),ste(y2,2),'.','MarkerSize',20)
    figXAxis([],'test',[0 3],1:2,{avName{iav},avName{otherAV}})
    figYAxis([],'% Correct',[0 1],0:0.2:1)
    figAxForm
    hline(pctCorrThresh,'k:')
    hline(0.5,'k-')
    [~,p] = ttest(y1,y2);
    title({sprintf('%s Stimulus Model, Ch. W Rmv',avName{iav});...
        sprintf('Attn, p=%s',sigfigString(p))})
    
    subplot(4,2,iav+2)
    hold on
    ind = mdl_attn(iav).mdl(choiceModInd).pctCorrect_fixedStimW_pcsOnly > pctCorrThresh;
    y1 = mdl_attn(iav).mdl(choiceModInd).pctCorrect_fixedStimW_pcsOnly(ind);
    y2 = mdl_attn(iav).mdl(choiceModInd).pctCorrect_fixedStimW_otherAV(ind);
    plot(1:2,[y1',y2'],'k-')
    errorbar(1,mean(y1),ste(y1,2),'.','MarkerSize',20)
    errorbar(2,mean(y2),ste(y2,2),'.','MarkerSize',20)
    figXAxis([],'test',[0 3],1:2,{avName{iav},avName{otherAV}})
    figYAxis([],'% Correct',[0 1],0:0.2:1)
    figAxForm
    hline(pctCorrThresh,'k:')
    hline(0.5,'k-')
    [~,p] = ttest(y1,y2);
    title({sprintf('%s Choice Model, Stim W Rmv',avName{iav});...
        sprintf('Attn, p=%s',sigfigString(p))})
    
    subplot(4,2,iav+4)
    hold on
    ind = mdl_noAttn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW_pcsOnly > pctCorrThresh;
    y1 = mdl_noAttn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW_pcsOnly(ind);
    y2 = mdl_noAttn(iav).mdl(stimModInd).pctCorrect_fixedChoiceW_otherAV(ind);
    plot(1:2,[y1',y2'],'k-')
    errorbar(1,mean(y1),ste(y1,2),'.','MarkerSize',20)
    errorbar(2,mean(y2),ste(y2,2),'.','MarkerSize',20)
    figXAxis([],'test',[0 3],1:2,{avName{iav},avName{otherAV}})
    figYAxis([],'% Correct',[0 1],0:0.2:1)
    figAxForm
    hline(pctCorrThresh,'k:')
    hline(0.5,'k-')
    [~,p] = ttest(y1,y2);
    title({sprintf('%s Stimulus Model, Ch. W Rmv',avName{iav});...
        sprintf('No Attn, p=%s',sigfigString(p))})
    
    subplot(4,2,iav+6)
    hold on
    ind = mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_fixedStimW_pcsOnly > pctCorrThresh;
    y1 = mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_fixedStimW_pcsOnly(ind);
    y2 = mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_fixedStimW_otherAV(ind);
    plot(1:2,[y1',y2'],'k-')
    errorbar(1,mean(y1),ste(y1,2),'.','MarkerSize',20)
    errorbar(2,mean(y2),ste(y2,2),'.','MarkerSize',20)
    figXAxis([],'test',[0 3],1:2,{avName{iav},avName{otherAV}})
    figYAxis([],'% Correct',[0 1],0:0.2:1)
    figAxForm
    hline(pctCorrThresh,'k:')
    hline(0.5,'k-')
    [~,p] = ttest(y1,y2);
    title({sprintf('%s Choice Model, Stim W Rmv',avName{iav});...
        sprintf('No Attn, p=%s',sigfigString(p))})
end
print([fnout 'pctCorrect_inclStimChoiceModelNotPred_testOtherAV'],'-dpdf','-fillpage')

%% weight analysis with fixed stim and choice paramerters
mdlWeightsSub_fig = figure;
mdlWeights_fig = figure;
suptitle({'Weights - from fixed stim/choice models';...
    sprintf('pct corr must be > %s in vis & aud stim models',...
    num2str(pctCorrThresh))})
oriEdges_noAbs = [-90:15:90];
for imod = 1:2
    if imod == 1
        w = 'weights_withFixedChoice';
    else
        w = 'weights_withFixedStim';
    end
    figure(mdlWeights_fig)
    subplot(3,4,imod)
    
    y1 = cell2mat(eval(['mdl_attn(visualTrials).mdl(imod).' w '(goodPerfInd_av_attn)'])');
    y2 = cell2mat(eval(['mdl_attn(auditoryTrials).mdl(imod).' w '(goodPerfInd_av_attn)'])');
    plot(y1,y2,'k.','MarkerSize',10)
    hold on
    mdl_fitline = fitlm(y1,y2);
    rsq = mdl_fitline.Rsquared.Ordinary;
    yfit = predict(mdl_fitline,y1);
    plot(y1,yfit,'-')
    figXAxis([],'Visual Weight',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Auditory Weight',[weightEdges(1) weightEdges(end)])
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    r = corr(y1,y2);
    title({'Attn Mice';sprintf('%s Model, r=%s',...
        mdl_attn(iav).mdl(imod).name,num2str(round(r,2,'significant')))})
    
    subplot(3,4,imod+4)
    weightAVangle = rad2deg(atan(y2./y1));
    histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
    figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
    figYAxis([],'Fraction of Cells',[0 .4])
    figAxForm
    title(sprintf('Attn Mice %s Model',...
        mdl_attn(iav).mdl(imod).name))
        
    subplot(3,4,imod+8)
    avAngleShuff = mdl_attn(visualTrials).mdl(imod).weightAVangle_noAbs_shuff;
    h = nan(length(oriEdges_noAbs)-1,nShuff);
    for ishuff = 1:nShuff
        h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges_noAbs,'Normalization','Probability');
    end    
    h_sort = sort(h,2);
    y = mean(h_sort,2);
    yerrl = y-h_sort(:,round(nShuff*0.025));
    yerru = h_sort(:,round(nShuff*0.975))-y;
    shadedErrorBar_chooseColor(oriEdges_noAbs(1:end-1),y,[yerrl';yerru'],[0 0 0]);
    hold on
    plot(oriEdges_noAbs(1:end-1),histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability'),'.-')
    controlSub_attn = histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability') - y';
    figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
    figYAxis([],'Fraction of Cells',[0 .4])
    figAxForm
    title(sprintf('Attn Mice-Shuff, %s Model',...
        mdl_attn(iav).mdl(imod).name))
    
    subplot(3,4,imod+2)
    y1 = cell2mat(eval(['mdl_noAttn(visualTrials).mdl(imod).' w '(goodPerfInd_av_noAttn)'])');
    y2 = cell2mat(eval(['mdl_noAttn(auditoryTrials).mdl(imod).' w '(goodPerfInd_av_noAttn)'])');
    plot(y1,y2,'k.','MarkerSize',10)
    hold on
    mdl_fitline = fitlm(y1,y2);
    rsq = mdl_fitline.Rsquared.Ordinary;
    yfit = predict(mdl_fitline,y1);
    plot(y1,yfit,'-')
    figXAxis([],'Visual Weight',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Auditory Weight',[weightEdges(1) weightEdges(end)])
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    r = corr(y1,y2);
    title({'No Attn Mice';sprintf('%s Model, r=%s',...
        mdl_noAttn(iav).mdl(imod).name,num2str(round(r,2,'significant')))})
    
    subplot(3,4,imod+6)
    weightAVangle = rad2deg(atan(y2./y1));
    histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
    figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
    figYAxis([],'Fraction of Cells',[0 .4])
    figAxForm
    title(sprintf('No Attn Mice, %s Model',...
        mdl_attn(iav).mdl(imod).name))
       
    subplot(3,4,imod+10)
    avAngleShuff = mdl_noAttn(visualTrials).mdl(imod).weightAVangle_noAbs_shuff;
    h = nan(length(oriEdges_noAbs)-1,nShuff);
    for ishuff = 1:nShuff
        h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges_noAbs,'Normalization','Probability');
    end    
    h_sort = sort(h,2);
    y = mean(h_sort,2);
    yerrl = y-h_sort(:,round(nShuff*0.025));
    yerru = h_sort(:,round(nShuff*0.975))-y;
    shadedErrorBar_chooseColor(oriEdges_noAbs(1:end-1),y,[yerrl';yerru'],[0 0 0]);
    hold on
    plot(oriEdges_noAbs(1:end-1),histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability'),'.-')
    controlSub_noAttn = histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability') - y';
    figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
    figYAxis([],'Fraction of Cells',[0 .4])
    figAxForm
    title(sprintf('Attn Mice-Shuff, %s Model',...
        mdl_noAttn(iav).mdl(imod).name))
    
    figure(mdlWeightsSub_fig)
    subplot(1,2,imod)
    plot(oriEdges_noAbs(1:end-1),controlSub_attn,'.-')
    hold on
%     shadedErrorBar_chooseColor(oriEdges(1:end-1),controlSub_noAttn,[yerrl';yerru'],[0.5 0.5 0.5]);
    plot(oriEdges_noAbs(1:end-1),controlSub_noAttn,'.-')
    figXAxis([],'Angle of Norm. A/V Weight',[-100 100],-90:15:90)
    figYAxis([],'Control Subtr. Fraction of Cells',[-0.15 0.15])
    figAxForm
    hline(0,'k:')
    title(sprintf('%s Model',modelName{imod}))
       
end 
figure(mdlWeights_fig)
print([fnout 'modelWeights_fixedStimChoice'],'-dpdf','-fillpage')
figure(mdlWeightsSub_fig)
print([fnout 'modelWeights_fixedStimChoice_ciSub'],'-dpdf')
%% compare neuron weights across model types

figure
suptitle('Choice Model')
for iav = 1:2
    subplot(2,2,iav)
    w1 = cell2mat(mdl_attn(iav).mdl(choiceModInd).weights(goodPerfInd_av_attn)');
    w2 = cell2mat(mdl_attn(iav).mdl(choiceModInd).weights_withStim(goodPerfInd_av_attn)');
    plot(w1,w2,'k.','MarkerSize',10)
    hold on
    plot(weightEdges,weightEdges,'k--')
    vline(0,'k--')
    hline(0,'k--')
    figXAxis([],'Weight - PCs only',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Weight - Fixed Stim',[weightEdges(1) weightEdges(end)])
    figAxForm
    title(sprintf('Attn %s Choice Model',avName{iav}))
    subplot(2,2,iav+2)
    w1 = cell2mat(mdl_noAttn(iav).mdl(choiceModInd).weights(goodPerfInd_av_noAttn)');
    w2 = cell2mat(mdl_noAttn(iav).mdl(choiceModInd).weights_withStim(goodPerfInd_av_noAttn)');
    plot(w1,w2,'k.','MarkerSize',10)
    hold on
    plot(weightEdges,weightEdges,'k--')
    vline(0,'k--')
    hline(0,'k--')
    figXAxis([],'Weight - PCs only',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Weight - Fixed Stim',[weightEdges(1) weightEdges(end)])
    figAxForm
    title(sprintf('No Attn %s Choice Model',avName{iav}))
end
print([fnout 'compareChoiceModelWeights_fixedVsPCs'],'-dpdf','-fillpage')

figure
suptitle('Stimulus Model')
for iav = 1:2
    subplot(2,2,iav)
    w1 = cell2mat(mdl_attn(iav).mdl(stimModInd).weights(goodPerfInd_av_attn)');
    w2 = cell2mat(mdl_attn(iav).mdl(stimModInd).weights_withChoice(goodPerfInd_av_attn)');
    plot(w1,w2,'k.','MarkerSize',10)
    hold on
    plot(weightEdges,weightEdges,'k--')
    vline(0,'k--')
    hline(0,'k--')
    figXAxis([],'Weight - PCs only',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Weight - Fixed Choice',[weightEdges(1) weightEdges(end)])
    figAxForm
    title(sprintf('Attn %s Stim Model',avName{iav}))
    subplot(2,2,iav+2)
    w1 = cell2mat(mdl_noAttn(iav).mdl(stimModInd).weights(goodPerfInd_av_noAttn)');
    w2 = cell2mat(mdl_noAttn(iav).mdl(stimModInd).weights_withChoice(goodPerfInd_av_noAttn)');
    plot(w1,w2,'k.','MarkerSize',10)
    hold on
    plot(weightEdges,weightEdges,'k--')
    vline(0,'k--')
    hline(0,'k--')
    figXAxis([],'Weight - PCs only',[weightEdges(1) weightEdges(end)])
    figYAxis([],'Weight - Fixed Choice',[weightEdges(1) weightEdges(end)])
    figAxForm
    title(sprintf('No Attn %s Stim Model',avName{iav}))
end
print([fnout 'compareStimModelWeights_fixedVsPCs'],'-dpdf','-fillpage')

%% compare fixed model difference and performance to mouse performance

msHR_attn = nan(2,nexp_attn);
for i = 1:nexp_attn
    for iav = 1:2
        trOut = dc_attn(i).av(iav).trOut;
        HR = sum(strcmp(trOut,'h'))/sum(strcmp(trOut,'h')|strcmp(trOut,'m'));
        msHR_attn(iav,i) = HR;
    end
end

msHR_noAttn = nan(2,nexp_noAttn);
for i = 1:nexp_noAttn
    for iav = 1:2
        trOut = dc_noAttn(i).av(iav).trOut;
        HR = sum(strcmp(trOut,'h'))/sum(strcmp(trOut,'h')|strcmp(trOut,'m'));
        msHR_noAttn(iav,i) = HR;
    end
end


figure
suptitle('Choice Models')
for iav = 1:2    
    ind = mdl_attn(iav).mdl(choiceModInd).pctCorr_all > pctCorrThresh;
    y1 = mdl_attn(iav).mdl(choiceModInd).pctCorr_all(ind);
    y2 = mdl_attn(iav).mdl(choiceModInd).pctCorrect_withStim(ind);
    y3 = (y2-y1)./y1;
    subplot(2,3,iav)
    plot(msHR_attn(ind),y3,'.','MarkerSize',20)
    figXAxis([],'Mouse Hit Rate',[0.5 1])
    figYAxis([], 'PCs to Ch. Inlc. Pct Change in Correct',[0 0.5])
    figAxForm
    title(sprintf('Attn %s Trials',avName{iav}))
    
    subplot(2,3,3)
    hold on
    plot(msHR_attn(ind),y3,'.','MarkerSize',20)
    if iav == 1
        yall = y3;
    else
        yall = cat(2,yall,y3);
        legend(avName,'location','northwest')      
        figXAxis([],'Mouse Hit Rate',[0.5 1])
        figYAxis([], 'PCs to Ch. Inlc. Pct Change in Correct',[0 0.5])
        figAxForm
        title(sprintf('Attn All Trials',avName{iav}))
    end
    
    figXAxis([],'Mouse Hit Rate',[0.5 1])
    figYAxis([], 'PCs to Ch. Inlc. Pct Change in Correct',[0 0.5])
    figAxForm
    title(sprintf('Attn %s Trials',avName{iav}))
    
    ind = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all > pctCorrThresh;
    y1 = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all(ind);
    y2 = mdl_noAttn(iav).mdl(choiceModInd).pctCorrect_withStim(ind);
    y3 = (y2-y1)./y1;
    subplot(2,3,iav+3)
    plot(msHR_noAttn(ind),y3,'.','MarkerSize',20)
    figXAxis([],'Mouse Hit Rate',[0.5 1])
    figYAxis([], 'PCs to Ch. Inlc. Pct Change in Correct',[0 0.5])
    figAxForm
    title(sprintf('No Attn %s Trials',avName{iav}))
    
    subplot(2,3,6)
    hold on
    plot(msHR_noAttn(ind),y3,'.','MarkerSize',20)
    if iav == 1
        yall = y3;
    else
        yall = cat(2,yall,y3);
        legend(avName,'location','northwest')        
        figXAxis([],'Mouse Hit Rate',[0.5 1])
        figYAxis([], 'PCs to Ch. Inlc. Pct Change in Correct',[0 0.5])
        figAxForm
        title(sprintf('No Attn All Trials',avName{iav}))
    end
end

print([fnout 'compareHRtoModelPctCorrect_fixedVsPCs'],'-dpdf','-fillpage')


%% naive weight analysis
for iav = 1:2
    mdl_naive(iav).mdl(stimModInd).weights = cell(1,nexp_naive);
    for iexp = 1:nexp_naive
        mdl_naive(iav).mdl(stimModInd).weights{iexp} = ...
            dc_naive(iexp).av(iav).weightTarget;
    end
end
goodPerfInd_av_naive = mdl_naive(visualTrials).mdl(stimModInd).pctCorr_all...
    > pctCorrThresh |...
    mdl_naive(auditoryTrials).mdl(stimModInd).pctCorr_all > pctCorrThresh;

for iav = 1:2
    rng(0)
    if iav == 1
        for ishuff = 1:nShuff
            shuffWeights_vis = cellfun(@(x) x(randperm(length(x))),...
                mdl_naive(visualTrials).mdl(stimModInd).weights(goodPerfInd_av_naive),'unif',0)';
            shuffWeights_aud = cellfun(@(x) x(randperm(length(x))),...
                mdl_naive(auditoryTrials).mdl(stimModInd).weights(goodPerfInd_av_naive),'unif',0)';
            if ishuff == 1
                mdl_naive(visualTrials).mdl(stimModInd).weightAVangle_noAbs_shuff = nan(sum(cellfun(@length,shuffWeights_vis)),nShuff);
            end
            mdl_naive(visualTrials).mdl(stimModInd).weightAVangle_noAbs_shuff(:,ishuff) = cell2mat(cellfun(@(x,y) ...
                rad2deg(atan(x./y)),shuffWeights_aud,...
                shuffWeights_vis,'unif',0));          
        end
    end
end

figure
subplot 221
y1 = cell2mat(mdl_naive(visualTrials).mdl(stimModInd).weights(goodPerfInd_av_naive)');
y2 = cell2mat(mdl_naive(auditoryTrials).mdl(stimModInd).weights(goodPerfInd_av_naive)');
plot(y1,y2,'k.','MarkerSize',10)
hold on
mdl_fitline = fitlm(y1,y2);
rsq = mdl_fitline.Rsquared.Ordinary;
yfit = predict(mdl_fitline,y1);
plot(y1,yfit,'-')
figXAxis([],'Visual Weight',[weightEdges(1) weightEdges(end)])
figYAxis([],'Auditory Weight',[weightEdges(1) weightEdges(end)])
figAxForm
hline(0,'k:')
vline(0,'k:')
r = corr(y1,y2);
title({'Naive Mice';sprintf('%s Model, r=%s',...
    mdl_naive(iav).mdl(stimModInd).name,num2str(round(r,2,'significant')))})

subplot 222
weightAVangle = rad2deg(atan(y2./y1));
histogram(weightAVangle,oriEdges_noAbs,'Normalization','probability')
figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
title(sprintf('Naive Mice %s Model',...
    mdl_naive(iav).mdl(stimModInd).name))

subplot 224
avAngleShuff = mdl_naive(visualTrials).mdl(stimModInd).weightAVangle_noAbs_shuff;
h = nan(length(oriEdges_noAbs)-1,nShuff);
for ishuff = 1:nShuff
    h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges_noAbs,'Normalization','Probability');
end    
h_sort = sort(h,2);
y = mean(h_sort,2);
yerrl = y-h_sort(:,round(nShuff*0.025));
yerru = h_sort(:,round(nShuff*0.975))-y;
shadedErrorBar_chooseColor(oriEdges_noAbs(1:end-1),y,[yerrl';yerru'],[0 0 0]);
hold on
plot(oriEdges_noAbs(1:end-1),histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability'),'.-')
controlSub_naive = histcounts(weightAVangle,oriEdges_noAbs,'Normalization','Probability') - y';
figXAxis([],'Angle of A/V Weight',[-100 100],-90:30:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
title(sprintf('Naive Mice-Shuff, %s Model',...
    mdl_naive(iav).mdl(stimModInd).name))

subplot 223
for iav = 1:2
    hold on
    y1 = cell2mat(mdl_naive(iav).mdl(stimModInd).weights(goodPerfInd_av_naive)');
    cdfplot(y1)
    figXAxis([],'Weight',[-1 2])
    figYAxis([],'Fraction of Cells',[0 1])
    figAxForm
    title(sprintf('Naive %s Stimulus Model',...
        avName{iav}))
    vline(0,'k:')
end
print([fnout 'weightAnalysis_naive'],'-dpdf','-fillpage')

%% weight analysis of cell-type groups
mdl_attn(visualTrials).mdl(choiceModInd).cellInd_lateResp = cell(1,nexp_attn);
mdl_attn(visualTrials).mdl(choiceModInd).cellInd_targetResp = cell(1,nexp_attn);
mdl_attn(visualTrials).mdl(choiceModInd).cellInd_tuned = cell(1,nexp_attn);
mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_lateResp = cell(1,nexp_noAttn);
mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_targetResp = cell(1,nexp_noAttn);
mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_tuned = cell(1,nexp_noAttn);
for iexp = 1:nexp_noAttn
    ind = respCellsExpt_attn(iexp).lateCycRespCells(dc_attn(iexp).cellInd);
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_lateResp{iexp} = ind;
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_targetResp{iexp} = ~ind;
    ind = oriTuning_attn(iexp).isTuned(dc_attn(iexp).cellInd);
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_tuned{iexp} = ind;
    ind = respCellsExpt_noAttn(iexp).lateCycRespCells(dc_noAttn(iexp).cellInd);
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_lateResp{iexp} = ind;
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_targetResp{iexp} = ~ind;
    ind = oriTuning_noAttn(iexp).isTuned(dc_noAttn(iexp).cellInd);
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_tuned{iexp} = ind;
end

setFigParams4Print('landscape')
figure
subplot 241
y1 = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_attn(visualTrials).mdl(choiceModInd).normWeights_distOnly(goodPerfInd_av_attn),...
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_attn),'unif',0)');
y2 = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_attn(auditoryTrials).mdl(choiceModInd).normWeights_distOnly(goodPerfInd_av_attn),...
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_attn),'unif',0)');
weightAVangle_attn = rad2deg(abs(atan(y2./y1)));
n_attn = length(y1);
%     [r_attn,p_attn,rCIl_attn,rCIu_attn] = corrcoef(y1,y2);
[r_attn,p_attn] = corr(y1,y2);
y1 = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).normWeights_distOnly(goodPerfInd_av_noAttn),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_noAttn),'unif',0)');
y2 = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).normWeights_distOnly(goodPerfInd_av_noAttn),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_noAttn),'unif',0)');n_attn = length(y1);
weightAVangle_noAttn = rad2deg(abs(atan(y2./y1)));
n_noAttn = length(y1);
[r_noAttn,p_noAttn] = corr(y1,y2);
bar(1,r_attn)
text(1,0.1,sigfigString(p_attn))
hold on
bar(2,r_noAttn)
text(2,0.1,sigfigString(p_noAttn))
figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
figYAxis([],'Vis-Aud Corr',[-0.1 1])
figAxForm
title(sprintf('%s Model, Distractors Only, late resp cells, p=%s',modelName{choiceModInd},...
    sigfigString(compare_correlation_coefficients(...
    r_attn,n_attn,r_noAttn,n_noAttn))))

subplot 242
cellInd = cell2mat(mdl_attn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_attn)')==1;
visW = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_attn(visualTrials).mdl(choiceModInd).normWeights_distOnly(goodPerfInd_av_attn),...
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_attn),'unif',0)');
audW = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_attn(auditoryTrials).mdl(choiceModInd).normWeights_distOnly(goodPerfInd_av_attn),...
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_attn),'unif',0)');
avAngleShuff = nan(sum(cellInd),nShuff);
for ishuff = 1:nShuff
%                 
    % shuffle normalized weights
    shuffWeights_vis = visW(randperm(sum(cellInd)));    
    shuffWeights_aud = audW(randperm(sum(cellInd)));
    
    avAngleShuff(:,ishuff) = rad2deg(abs(atan(shuffWeights_aud./shuffWeights_vis)));
end
h = nan(length(oriEdges)-1,nShuff);
for ishuff = 1:nShuff
    h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges,'Normalization','Probability');
end    
h_sort = sort(h,2);
y = mean(h_sort,2);
yerrl = y-h_sort(:,round(nShuff*0.025));
yerru = h_sort(:,round(nShuff*0.975))-y;
shadedErrorBar_chooseColor(oriEdges(1:end-1),y,[yerrl';yerru'],[0 0 0]);
hold on
plot(oriEdges(1:end-1),histcounts(weightAVangle_attn,oriEdges,'Normalization','Probability'),'.-')
controlSub_attn = histcounts(weightAVangle_attn,oriEdges,'Normalization','Probability') - y';
figXAxis([],'Angle of Norm. A/V Weight',[-10 100],0:15:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
title({'Attn Mice-Shuff';sprintf('%s Model, late resp cells',...
    modelName{choiceModInd})})

subplot 243
cellInd = cell2mat(mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_noAttn)')==1;
visW = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).normWeights_distOnly(goodPerfInd_av_noAttn),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_noAttn),'unif',0)');
audW = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).normWeights_distOnly(goodPerfInd_av_noAttn),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_noAttn),'unif',0)');
avAngleShuff = nan(sum(cellInd),nShuff);
for ishuff = 1:nShuff
%                 
    % shuffle normalized weights
    shuffWeights_vis = visW(randperm(sum(cellInd)));    
    shuffWeights_aud = audW(randperm(sum(cellInd)));
    
    avAngleShuff(:,ishuff) = rad2deg(abs(atan(shuffWeights_aud./shuffWeights_vis)));
end
h = nan(length(oriEdges)-1,nShuff);
for ishuff = 1:nShuff
    h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges,'Normalization','Probability');
end    
h_sort = sort(h,2);
y = mean(h_sort,2);
yerrl = y-h_sort(:,round(nShuff*0.025));
yerru = h_sort(:,round(nShuff*0.975))-y;
shadedErrorBar_chooseColor(oriEdges(1:end-1),y,[yerrl';yerru'],[0 0 0]);
hold on
plot(oriEdges(1:end-1),histcounts(weightAVangle_noAttn,oriEdges,'Normalization','Probability'),'.-')
controlSub_noAttn = histcounts(weightAVangle_noAttn,oriEdges,'Normalization','Probability') - y';
figXAxis([],'Angle of Norm. A/V Weight',[-10 100],0:15:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
title({'Attn Mice-Shuff';sprintf('%s Model, late resp cells',...
    modelName{choiceModInd})})

subplot 244
hold on
plot(oriEdges(1:end-1),controlSub_attn,'.-')
plot(oriEdges(1:end-1),controlSub_noAttn,'.-')
figXAxis([],'Angle of Norm. A/V Weight',[-10 100],0:15:90)
figYAxis([],'Control Subtr. Fraction of Cells',[-0.15 0.15])
figAxForm
hline(0,'k:')
title(sprintf('%s Model, Distractors Only, late resp cells',modelName{choiceModInd}))


subplot 245
y1 = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_attn(visualTrials).mdl(choiceModInd).normWeights_distOnly(goodPerfInd_av_attn),...
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_attn),'unif',0)');
y2 = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_attn(auditoryTrials).mdl(choiceModInd).normWeights_distOnly(goodPerfInd_av_attn),...
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_attn),'unif',0)');
weightAVangle_attn = rad2deg(abs(atan(y2./y1)));
n_attn = length(y1);
%     [r_attn,p_attn,rCIl_attn,rCIu_attn] = corrcoef(y1,y2);
[r_attn,p_attn] = corr(y1,y2);
y1 = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).normWeights_distOnly(goodPerfInd_av_noAttn),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_noAttn),'unif',0)');
y2 = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).normWeights_distOnly(goodPerfInd_av_noAttn),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_noAttn),'unif',0)');n_attn = length(y1);
weightAVangle_noAttn = rad2deg(abs(atan(y2./y1)));
n_noAttn = length(y1);
[r_noAttn,p_noAttn] = corr(y1,y2);
bar(1,r_attn)
text(1,0.1,sigfigString(p_attn))
hold on
bar(2,r_noAttn)
text(2,0.1,sigfigString(p_noAttn))
figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
figYAxis([],'Vis-Aud Corr',[-0.1 1])
figAxForm
title(sprintf('%s Model, Distractors Only, target resp cells, p=%s',modelName{choiceModInd},...
    sigfigString(compare_correlation_coefficients(...
    r_attn,n_attn,r_noAttn,n_noAttn))))

subplot 246
cellInd = cell2mat(mdl_attn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_attn)')==1;
visW = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_attn(visualTrials).mdl(choiceModInd).normWeights_distOnly(goodPerfInd_av_attn),...
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_attn),'unif',0)');
audW = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_attn(auditoryTrials).mdl(choiceModInd).normWeights_distOnly(goodPerfInd_av_attn),...
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_attn),'unif',0)');
avAngleShuff = nan(sum(cellInd),nShuff);
for ishuff = 1:nShuff
%                 
    % shuffle normalized weights
    shuffWeights_vis = visW(randperm(sum(cellInd)));    
    shuffWeights_aud = audW(randperm(sum(cellInd)));
    
    avAngleShuff(:,ishuff) = rad2deg(abs(atan(shuffWeights_aud./shuffWeights_vis)));
end
h = nan(length(oriEdges)-1,nShuff);
for ishuff = 1:nShuff
    h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges,'Normalization','Probability');
end    
h_sort = sort(h,2);
y = mean(h_sort,2);
yerrl = y-h_sort(:,round(nShuff*0.025));
yerru = h_sort(:,round(nShuff*0.975))-y;
shadedErrorBar_chooseColor(oriEdges(1:end-1),y,[yerrl';yerru'],[0 0 0]);
hold on
plot(oriEdges(1:end-1),histcounts(weightAVangle_attn,oriEdges,'Normalization','Probability'),'.-')
controlSub_attn = histcounts(weightAVangle_attn,oriEdges,'Normalization','Probability') - y';
figXAxis([],'Angle of Norm. A/V Weight',[-10 100],0:15:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
title({'Attn Mice-Shuff';sprintf('%s Model, target resp cells',...
    modelName{choiceModInd})})

subplot 247
cellInd = cell2mat(mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_noAttn)')==1;
visW = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).normWeights_distOnly(goodPerfInd_av_noAttn),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_noAttn),'unif',0)');
audW = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).normWeights_distOnly(goodPerfInd_av_noAttn),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_noAttn),'unif',0)');
avAngleShuff = nan(sum(cellInd),nShuff);
for ishuff = 1:nShuff
%                 
    % shuffle normalized weights
    shuffWeights_vis = visW(randperm(sum(cellInd)));    
    shuffWeights_aud = audW(randperm(sum(cellInd)));
    
    avAngleShuff(:,ishuff) = rad2deg(abs(atan(shuffWeights_aud./shuffWeights_vis)));
end
h = nan(length(oriEdges)-1,nShuff);
for ishuff = 1:nShuff
    h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges,'Normalization','Probability');
end    
h_sort = sort(h,2);
y = mean(h_sort,2);
yerrl = y-h_sort(:,round(nShuff*0.025));
yerru = h_sort(:,round(nShuff*0.975))-y;
shadedErrorBar_chooseColor(oriEdges(1:end-1),y,[yerrl';yerru'],[0 0 0]);
hold on
plot(oriEdges(1:end-1),histcounts(weightAVangle_noAttn,oriEdges,'Normalization','Probability'),'.-')
controlSub_noAttn = histcounts(weightAVangle_noAttn,oriEdges,'Normalization','Probability') - y';
figXAxis([],'Angle of Norm. A/V Weight',[-10 100],0:15:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
title({'Attn Mice-Shuff';sprintf('%s Model, target resp cells',...
    modelName{choiceModInd})})

subplot 248
hold on
plot(oriEdges(1:end-1),controlSub_attn,'.-')
plot(oriEdges(1:end-1),controlSub_noAttn,'.-')
figXAxis([],'Angle of Norm. A/V Weight',[-10 100],0:15:90)
figYAxis([],'Control Subtr. Fraction of Cells',[-0.15 0.15])
figAxForm
hline(0,'k:')
title(sprintf('%s Model, Distractors Only, target resp cells',modelName{choiceModInd}))

print([fnout 'weightAnalysis_distOnlyChoice_cellTypes'],'-dpdf','-fillpage')
%%
setFigParams4Print('landscape')
figure
subplot 241
y1 = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_attn(visualTrials).mdl(choiceModInd).normWeights_tarOnly(goodPerfInd_av_attn),...
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_attn),'unif',0)');
y2 = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_attn(auditoryTrials).mdl(choiceModInd).normWeights_tarOnly(goodPerfInd_av_attn),...
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_attn),'unif',0)');
weightAVangle_attn = rad2deg(abs(atan(y2./y1)));
n_attn = length(y1);
%     [r_attn,p_attn,rCIl_attn,rCIu_attn] = corrcoef(y1,y2);
[r_attn,p_attn] = corr(y1,y2);
y1 = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).normWeights_tarOnly(goodPerfInd_av_noAttn),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_noAttn),'unif',0)');
y2 = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).normWeights_tarOnly(goodPerfInd_av_noAttn),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_noAttn),'unif',0)');n_attn = length(y1);
weightAVangle_noAttn = rad2deg(abs(atan(y2./y1)));
n_noAttn = length(y1);
[r_noAttn,p_noAttn] = corr(y1,y2);
bar(1,r_attn)
text(1,0.1,sigfigString(p_attn))
hold on
bar(2,r_noAttn)
text(2,0.1,sigfigString(p_noAttn))
figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
figYAxis([],'Vis-Aud Corr',[-0.1 1])
figAxForm
title(sprintf('%s Model, Targets Only, late resp cells, p=%s',modelName{choiceModInd},...
    sigfigString(compare_correlation_coefficients(...
    r_attn,n_attn,r_noAttn,n_noAttn))))

subplot 242
cellInd = cell2mat(mdl_attn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_attn)')==1;
visW = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_attn(visualTrials).mdl(choiceModInd).normWeights_tarOnly(goodPerfInd_av_attn),...
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_attn),'unif',0)');
audW = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_attn(auditoryTrials).mdl(choiceModInd).normWeights_tarOnly(goodPerfInd_av_attn),...
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_attn),'unif',0)');
avAngleShuff = nan(sum(cellInd),nShuff);
for ishuff = 1:nShuff
%                 
    % shuffle normalized weights
    shuffWeights_vis = visW(randperm(sum(cellInd)));    
    shuffWeights_aud = audW(randperm(sum(cellInd)));
    
    avAngleShuff(:,ishuff) = rad2deg(abs(atan(shuffWeights_aud./shuffWeights_vis)));
end
h = nan(length(oriEdges)-1,nShuff);
for ishuff = 1:nShuff
    h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges,'Normalization','Probability');
end    
h_sort = sort(h,2);
y = mean(h_sort,2);
yerrl = y-h_sort(:,round(nShuff*0.025));
yerru = h_sort(:,round(nShuff*0.975))-y;
shadedErrorBar_chooseColor(oriEdges(1:end-1),y,[yerrl';yerru'],[0 0 0]);
hold on
plot(oriEdges(1:end-1),histcounts(weightAVangle_attn,oriEdges,'Normalization','Probability'),'.-')
controlSub_attn = histcounts(weightAVangle_attn,oriEdges,'Normalization','Probability') - y';
figXAxis([],'Angle of Norm. A/V Weight',[-10 100],0:15:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
title({'Attn Mice-Shuff';sprintf('%s Model, late resp cells',...
    modelName{choiceModInd})})

subplot 243
cellInd = cell2mat(mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_noAttn)')==1;
visW = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).normWeights_tarOnly(goodPerfInd_av_noAttn),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_noAttn),'unif',0)');
audW = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).normWeights_tarOnly(goodPerfInd_av_noAttn),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_lateResp(goodPerfInd_av_noAttn),'unif',0)');
avAngleShuff = nan(sum(cellInd),nShuff);
for ishuff = 1:nShuff
%                 
    % shuffle normalized weights
    shuffWeights_vis = visW(randperm(sum(cellInd)));    
    shuffWeights_aud = audW(randperm(sum(cellInd)));
    
    avAngleShuff(:,ishuff) = rad2deg(abs(atan(shuffWeights_aud./shuffWeights_vis)));
end
h = nan(length(oriEdges)-1,nShuff);
for ishuff = 1:nShuff
    h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges,'Normalization','Probability');
end    
h_sort = sort(h,2);
y = mean(h_sort,2);
yerrl = y-h_sort(:,round(nShuff*0.025));
yerru = h_sort(:,round(nShuff*0.975))-y;
shadedErrorBar_chooseColor(oriEdges(1:end-1),y,[yerrl';yerru'],[0 0 0]);
hold on
plot(oriEdges(1:end-1),histcounts(weightAVangle_noAttn,oriEdges,'Normalization','Probability'),'.-')
controlSub_noAttn = histcounts(weightAVangle_noAttn,oriEdges,'Normalization','Probability') - y';
figXAxis([],'Angle of Norm. A/V Weight',[-10 100],0:15:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
title({'Attn Mice-Shuff';sprintf('%s Model, late resp cells',...
    modelName{choiceModInd})})

subplot 244
hold on
plot(oriEdges(1:end-1),controlSub_attn,'.-')
plot(oriEdges(1:end-1),controlSub_noAttn,'.-')
figXAxis([],'Angle of Norm. A/V Weight',[-10 100],0:15:90)
figYAxis([],'Control Subtr. Fraction of Cells',[-0.15 0.15])
figAxForm
hline(0,'k:')
title(sprintf('%s Model, Targets Only, late resp cells',modelName{choiceModInd}))


subplot 245
y1 = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_attn(visualTrials).mdl(choiceModInd).normWeights_tarOnly(goodPerfInd_av_attn),...
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_attn),'unif',0)');
y2 = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_attn(auditoryTrials).mdl(choiceModInd).normWeights_tarOnly(goodPerfInd_av_attn),...
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_attn),'unif',0)');
weightAVangle_attn = rad2deg(abs(atan(y2./y1)));
n_attn = length(y1);
%     [r_attn,p_attn,rCIl_attn,rCIu_attn] = corrcoef(y1,y2);
[r_attn,p_attn] = corr(y1,y2);
y1 = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).normWeights_tarOnly(goodPerfInd_av_noAttn),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_noAttn),'unif',0)');
y2 = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).normWeights_tarOnly(goodPerfInd_av_noAttn),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_noAttn),'unif',0)');n_attn = length(y1);
weightAVangle_noAttn = rad2deg(abs(atan(y2./y1)));
n_noAttn = length(y1);
[r_noAttn,p_noAttn] = corr(y1,y2);
bar(1,r_attn)
text(1,0.1,sigfigString(p_attn))
hold on
bar(2,r_noAttn)
text(2,0.1,sigfigString(p_noAttn))
figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
figYAxis([],'Vis-Aud Corr',[-0.1 1])
figAxForm
title(sprintf('%s Model, Targets Only, target resp cells, p=%s',modelName{choiceModInd},...
    sigfigString(compare_correlation_coefficients(...
    r_attn,n_attn,r_noAttn,n_noAttn))))

subplot 246
cellInd = cell2mat(mdl_attn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_attn)')==1;
visW = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_attn(visualTrials).mdl(choiceModInd).normWeights_tarOnly(goodPerfInd_av_attn),...
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_attn),'unif',0)');
audW = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_attn(auditoryTrials).mdl(choiceModInd).normWeights_tarOnly(goodPerfInd_av_attn),...
    mdl_attn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_attn),'unif',0)');
avAngleShuff = nan(sum(cellInd),nShuff);
for ishuff = 1:nShuff
%                 
    % shuffle normalized weights
    shuffWeights_vis = visW(randperm(sum(cellInd)));    
    shuffWeights_aud = audW(randperm(sum(cellInd)));
    
    avAngleShuff(:,ishuff) = rad2deg(abs(atan(shuffWeights_aud./shuffWeights_vis)));
end
h = nan(length(oriEdges)-1,nShuff);
for ishuff = 1:nShuff
    h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges,'Normalization','Probability');
end    
h_sort = sort(h,2);
y = mean(h_sort,2);
yerrl = y-h_sort(:,round(nShuff*0.025));
yerru = h_sort(:,round(nShuff*0.975))-y;
shadedErrorBar_chooseColor(oriEdges(1:end-1),y,[yerrl';yerru'],[0 0 0]);
hold on
plot(oriEdges(1:end-1),histcounts(weightAVangle_attn,oriEdges,'Normalization','Probability'),'.-')
controlSub_attn = histcounts(weightAVangle_attn,oriEdges,'Normalization','Probability') - y';
figXAxis([],'Angle of Norm. A/V Weight',[-10 100],0:15:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
title({'Attn Mice-Shuff';sprintf('%s Model, target resp cells',...
    modelName{choiceModInd})})

subplot 247
cellInd = cell2mat(mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_noAttn)')==1;
visW = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).normWeights_tarOnly(goodPerfInd_av_noAttn),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_noAttn),'unif',0)');
audW = cell2mat(cellfun(@(x,y) x(y==1),...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).normWeights_tarOnly(goodPerfInd_av_noAttn),...
    mdl_noAttn(visualTrials).mdl(choiceModInd).cellInd_targetResp(goodPerfInd_av_noAttn),'unif',0)');
avAngleShuff = nan(sum(cellInd),nShuff);
for ishuff = 1:nShuff
%                 
    % shuffle normalized weights
    shuffWeights_vis = visW(randperm(sum(cellInd)));    
    shuffWeights_aud = audW(randperm(sum(cellInd)));
    
    avAngleShuff(:,ishuff) = rad2deg(abs(atan(shuffWeights_aud./shuffWeights_vis)));
end
h = nan(length(oriEdges)-1,nShuff);
for ishuff = 1:nShuff
    h(:,ishuff) = histcounts(avAngleShuff(:,ishuff),oriEdges,'Normalization','Probability');
end    
h_sort = sort(h,2);
y = mean(h_sort,2);
yerrl = y-h_sort(:,round(nShuff*0.025));
yerru = h_sort(:,round(nShuff*0.975))-y;
shadedErrorBar_chooseColor(oriEdges(1:end-1),y,[yerrl';yerru'],[0 0 0]);
hold on
plot(oriEdges(1:end-1),histcounts(weightAVangle_noAttn,oriEdges,'Normalization','Probability'),'.-')
controlSub_noAttn = histcounts(weightAVangle_noAttn,oriEdges,'Normalization','Probability') - y';
figXAxis([],'Angle of Norm. A/V Weight',[-10 100],0:15:90)
figYAxis([],'Fraction of Cells',[0 .4])
figAxForm
title({'Attn Mice-Shuff';sprintf('%s Model, target resp cells',...
    modelName{choiceModInd})})

subplot 248
hold on
plot(oriEdges(1:end-1),controlSub_attn,'.-')
plot(oriEdges(1:end-1),controlSub_noAttn,'.-')
figXAxis([],'Angle of Norm. A/V Weight',[-10 100],0:15:90)
figYAxis([],'Control Subtr. Fraction of Cells',[-0.15 0.15])
figAxForm
hline(0,'k:')
title(sprintf('%s Model, Targets Only, target resp cells',modelName{choiceModInd}))

print([fnout 'weightAnalysis_targetOnlyChoice_cellTypes'],'-dpdf','-fillpage')


