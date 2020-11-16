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
    [titleStr '_' datestr(now,'yymmdd') '_decoding_tunedNeurons_pcs_neurons_lessneurons']);
 
%% indexing and naming variables
avName = {'Vis','Aud'};
targetName = {'D','HT','ET'};
modelName = {'Stimulus','Choice'};

nexp_attn = size(dc_attn,2);
nexp_noAttn = size(dc_noAttn,2);

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
        x = 1:2;
        y_all = [mean(y1),mean(y2)];
        y_all_ste = [ste(y1,2),ste(y2,2)];
        figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})        
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

print([fnout 'pctCorrect_neurons'],'-dpdf','-fillpage')

fprintf('Attn mice - vis stim model successful expt: %spct %s/%s expt\n',...
    sigfigString(sum(goodPerfInd_attn)/length(goodPerfInd_attn)*100),...
    num2str(sum(goodPerfInd_attn)),num2str(length(goodPerfInd_attn)))
fprintf('No Attn mice - vis stim model successful expt: %spct %s/%s expt\n',...
    sigfigString(sum(goodPerfInd_noAttn)/length(goodPerfInd_noAttn)*100),...
    num2str(sum(goodPerfInd_noAttn)),num2str(length(goodPerfInd_noAttn)))

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
    title({'No Attn Mice';sprintf('%s Model',...
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

%% distractor only choice weight correlation with bootstrap quantification
weight_lim_dist = [-6.5 6.5];
nboot = 1000;
figure
suptitle('Distractor Only Model, Neurons')

    subplot 131
    y1 = cell2mat(mdl_attn(visualTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_attn)');
    y2 = cell2mat(mdl_attn(auditoryTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_attn)');
    plot(y1,y2,'k.','MarkerSize',10)
    hold on
    mdl_fitline = fitlm(y1,y2);
    yfit = predict(mdl_fitline,y1);
    plot(y1,yfit,'-')
    figXAxis([],'Visual Weight',[weight_lim_dist])
    figYAxis([],'Auditory Weight',[weight_lim_dist])
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    title({'Attn Mice';sprintf('%s Model',...
        mdl_attn(iav).mdl(choiceModInd).name)})
    
    r_attn = nan(1,nboot);
    for iboot = 1:nboot
        n = length(y1);
        ind = randsample(1:n,n,1);
        y1_samp = y1(ind);
        y2_samp = y2(ind);
        r_attn(iboot) = corr(y1_samp,y2_samp);
    end
    
    subplot 132
    y1 = cell2mat(mdl_noAttn(visualTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_noAttn)');
    y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_noAttn)');
    plot(y1,y2,'k.','MarkerSize',10)
    hold on
    mdl_fitline = fitlm(y1,y2);
    yfit = predict(mdl_fitline,y1);
    plot(y1,yfit,'-')
%     figXAxis([],'Visual Weight',[weight_lim_dist])
%     figYAxis([],'Auditory Weight',[weight_lim_dist])
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    title({'No Attn Mice';sprintf('%s Model',...
        mdl_noAttn(iav).mdl(imod).name)})
    
    r_noAttn = nan(1,nboot);
    for iboot = 1:nboot
        n = length(y1);
        ind = randsample(1:n,n,1);
        y1_samp = y1(ind);
        y2_samp = y2(ind);
        r_noAttn(iboot) = corr(y1_samp,y2_samp);
    end
    
    subplot 133; hold on
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

print([fnout 'modelWeights_distOnly_bootCorr'],'-dpdf')


%% compare PC weights
weight_lim_pcs = [-7 6];
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
weight_lim_neurons = [-2.5 2.5];
figure
suptitle('Neuron Weights')
for iav = 1:2
    if iav == 1
        iplot = 0;
    else
        iplot = 3;
    end
    subplot(2,3,1+iplot); hold on
    y1 = cell2mat(mdl_attn(iav).mdl(choiceModInd).weights(goodPerfInd_av_distOnly_attn)');
    y2 = cell2mat(mdl_attn(iav).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_attn)');
    plot(y1,y2,'k.','MarkerSize',10)
    hold on
    mdl_fitline = fitlm(y1,y2);
    rsq = mdl_fitline.Rsquared.Ordinary;
    yfit = predict(mdl_fitline,y1);
    plot(y1,yfit,'-')
    figXAxis([],'All Stim Weight',weight_lim_neurons)
    figYAxis([],'Dist Onlly Weight',weight_lim_neurons)
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
    y1 = cell2mat(mdl_noAttn(iav).mdl(choiceModInd).weights(goodPerfInd_av_distOnly_noAttn)');
    y2 = cell2mat(mdl_noAttn(iav).mdl(choiceModInd).weights_distOnly(goodPerfInd_av_distOnly_noAttn)');
    plot(y1,y2,'k.','MarkerSize',10)
    hold on
    mdl_fitline = fitlm(y1,y2);
    rsq = mdl_fitline.Rsquared.Ordinary;
    yfit = predict(mdl_fitline,y1);
    plot(y1,yfit,'-')
    figXAxis([],'All Stim Weight',weight_lim_neurons)
    figYAxis([],'Dist Onlly Weight',weight_lim_neurons)
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

print([fnout 'modelWeights_compareAll2Dist'],'-dpdf','-fillpage')
%% direct comparison all trials model and distractor only model performance

for iav = 1:2
    mdl_attn(iav).mdl(choiceModInd).pctCorr_all_testDist = nan(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(choiceModInd).pctCorr_all_testDist(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectXStimDetect_holdout(1);
    end
    
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all_testDist = nan(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all_testDist(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectXStimDetect_holdout(1);
    end
end
figure
suptitle('Test Same Modality Neurons Only')
x_all = 1:3;
for iav = 1:2
    subplot(2,2,iav); hold on
    y1 = mdl_attn(iav).mdl(choiceModInd).pctCorr_all;
    y2 = mdl_attn(iav).mdl(choiceModInd).pctCorr_all_testDist;
    y3 = mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly;
    
    ind1 = y1>pctCorrThresh;
    ind2 = y3>pctCorrThresh;
    
    for i = 1:nexp_attn
        if ind1(i) & ind2(i)
            y = cat(1,y1(i),y2(i),y3(i));
            x = x_all;
        elseif ind1(i) & ~ind2(i)
            y = cat(1,y1(i),y2(i));
            x = x_all(1:2);
        elseif ~ind1(i) & ind2(i)
            y = y3(i);
            x = x_all(3);
        end
        plot(x,y,'k.-','MarkerSize',10)
    end
    figXAxis([],'Model-Test',[0 4],x_all,{'All-All','All-Dist','Dist-Dist'})        
    figYAxis([],'% Correct',[0 1],0:0.2:1)
    figAxForm([],0)
    hline(pctCorrThresh,'k:')
    hline(0.5,'k-')
    title(sprintf('Attn %s %s Model',mdl_attn(iav).name,...
        mdl_attn(iav).mdl(choiceModInd).name))
    
    subplot(2,2,iav+2); hold on
    y1 = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all;
    y2 = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all_testDist;
    y3 = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly;
    
    ind1 = y1>pctCorrThresh;
    ind2 = y3>pctCorrThresh;
    
    for i = 1:nexp_noAttn
        if ind1(i) & ind2(i)
            y = cat(1,y1(i),y2(i),y3(i));
            x = x_all;
        elseif ind1(i) & ~ind2(i)
            y = cat(1,y1(i),y2(i));
            x = x_all(1:2);
        elseif ~ind1(i) & ind2(i)
            y = y3(i);
            x = x_all(3);
        end
        plot(x,y,'k.-','MarkerSize',10)
    end
    figXAxis([],'Model-Test',[0 4],x_all,{'All-All','All-Dist','Dist-Dist'})        
    figYAxis([],'% Correct',[0 1],0:0.2:1)
    figAxForm([],0)
    hline(pctCorrThresh,'k:')
    hline(0.5,'k-')
    title(sprintf('No Attn %s %s Model',mdl_noAttn(iav).name,...
        mdl_noAttn(iav).mdl(choiceModInd).name))
end
print([fnout 'modelPerformance_compareAll2Dist'],'-dpdf','-fillpage')

for iav = 1:2
    mdl_attn(iav).mdl(choiceModInd).pctCorr_all_otherAV_testDist = nan(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(choiceModInd).pctCorr_all_otherAV_testDist(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_onlyZero_otherAV;
    end
    
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all_otherAV_testDist = nan(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all_otherAV_testDist(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_onlyZero_otherAV;
    end
end

figure
suptitle('Test Opposite Modality (good experiments for same modality model')
x_all = 1:3;
for iav = 1:2
    subplot(2,2,iav); hold on
    y1 = mdl_attn(iav).mdl(choiceModInd).pctCorr_otherAV;
    y2 = mdl_attn(iav).mdl(choiceModInd).pctCorr_all_otherAV_testDist;
    y3 = mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly_otherAV;
    
    ind1 = mdl_attn(iav).mdl(choiceModInd).pctCorr_all>pctCorrThresh;
    ind2 = mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly>pctCorrThresh;
    
    for i = 1:nexp_attn
        if ind1(i) & ind2(i)
            y = cat(1,y1(i),y2(i),y3(i));
            x = x_all;
        elseif ind1(i) & ~ind2(i)
            y = cat(1,y1(i),y2(i));
            x = x_all(1:2);
        elseif ~ind1(i) & ind2(i)
            y = y3(i);
            x = x_all(3);
        end
        plot(x,y,'k.-','MarkerSize',10)
    end
    figXAxis([],'Model-Test',[0 4],x_all,{'All-All','All-Dist','Dist-Dist'})        
    figYAxis([],'% Correct',[0 1],0:0.2:1)
    figAxForm([],0)
    hline(pctCorrThresh,'k:')
    hline(0.5,'k-')
    title(sprintf('Attn %s %s Model',mdl_attn(iav).name,...
        mdl_attn(iav).mdl(choiceModInd).name))
    
    subplot(2,2,iav+2); hold on
    y1 = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_otherAV;
    y2 = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all_otherAV_testDist;
    y3 = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly_otherAV;
    
    ind1 = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all>pctCorrThresh;
    ind2 = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly>pctCorrThresh;
    
    for i = 1:nexp_noAttn
        if ind1(i) & ind2(i)
            y = cat(1,y1(i),y2(i),y3(i));
            x = x_all;
        elseif ind1(i) & ~ind2(i)
            y = cat(1,y1(i),y2(i));
            x = x_all(1:2);
        elseif ~ind1(i) & ind2(i)
            y = y3(i);
            x = x_all(3);
        end
        plot(x,y,'k.-','MarkerSize',10)
    end
    figXAxis([],'Model-Test',[0 4],x_all,{'All-All','All-Dist','Dist-Dist'})        
    figYAxis([],'% Correct',[0 1],0:0.2:1)
    figAxForm([],0)
    hline(pctCorrThresh,'k:')
    hline(0.5,'k-')
    title(sprintf('No Attn %s %s Model',mdl_noAttn(iav).name,...
        mdl_noAttn(iav).mdl(choiceModInd).name))
end
print([fnout 'modelPerformance_compareAll2Dist_otherAV'],'-dpdf','-fillpage')