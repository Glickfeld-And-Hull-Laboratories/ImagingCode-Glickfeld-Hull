clear all
close all
ds = 'FS_V1_naive_LG';
bxExpt = 1;
cellsOrDendrites = 1;
doLoadPreviousAnalysis = true;
%%
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms

eval(ds)
titleStr = ds(6:end);
mice = unique({expt.SubNum});
mouse_str = ['i' strjoin(mice,'_i')];
fnout = fullfile(rc.lindseyAnalysis,ds);
load(fullfile(fnout, 'trOutcomeStruct_cells.mat'));

%%
nBaselineFr = mouse(1).expt(1).info.preAlignFrames;
nFrames1s = frameRateHz;
nexp = size(expt,2);
nCycles = 7;
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

trOutType = {'h';'m';'fa';'cr'};
trOutTypeName = {'H-All';'H-HT';'H-ET';'M-All';'M-HT';'M-ET';'FA';'CR'};
%% pool experiment data

    antiDataExpt = struct;
    oriTuningExpt = struct;
    for imouse = 1:size(mouse,2)
        for iexp = 1:size(mouse(imouse).expt,2)
            if imouse == 1 && iexp == 1
                exptN = 1;
            else
                exptN = exptN+1;
            end
            fprintf([mouse(imouse).mouse_name ' ' mouse(imouse).expt(iexp).date '\n'])
            d = mouse(imouse).expt(iexp);
            exptID = strcmp({expt.SubNum},mouse(imouse).mouse_name) & strcmp({expt.date},d.date);

            cycLengthFr = d.info.cycTimeFrames;
            nCycLongTC = ceil(longTrialLengthFr./cycLengthFr);

            maxCycles = max(cat(2,d.align(alignStart).nCycles,...
                d.align(alignStart).nCycles));

            cycTC = cell(1,maxCycles);
%             longTC = [];
            dd = d.align(alignStart);
            if bxExpt
                hits = strcmp(dd.outcome,'success');
                misses = strcmp(dd.outcome,'ignore');
            else
                hits = true(1,length(dd.outcome));
                misses = false(1,length(dd.outcome));
            end
            for icyc = 1:maxCycles
                tc = dd.respTC(:,:,dd.nCycles >= icyc & (hits | misses));
                cycStartOffset = ((icyc-1).*cycLengthFr)+nBaselineFr;
                cycTC{icyc} = tc(...
                    (cycStartOffset-nBaselineFr+1):(cycStartOffset+nFrames1s),:,:);
            end
            longTC = dd.respTC(1:(longTrialLengthFr+nBaselineFr),:,...
                dd.nCycles >= nCycLongTC & (hits | misses));

            antiDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
            antiDataExpt(exptN).exptCycLengthFr = cycLengthFr;
            antiDataExpt(exptN).longTC = longTC;
            antiDataExpt(exptN).cycTC = cycTC;
            fprintf(['TC- n=' num2str(size(cycTC{1},2)) '\n'])
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
            fprintf(['Ori- n=' num2str(length(oriPref)) '\n'])
        end
    end

    shortCycExptInd = cell2mat({antiDataExpt.exptCycLengthFr}) == 11;
%     isShortCycExpt = [];
    respCellsExpt = struct;
    for iexp = 1:nexp
        firstTC = antiDataExpt(iexp).cycTC{1};

        firstRespCells = ttest(...
            squeeze(mean(firstTC(respwin,:,:),1)),...
            squeeze(mean(firstTC(basewin_0,:,:),1)),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha);

        respCellsExpt(iexp).exptName = antiDataExpt(iexp).exptName;
        respCellsExpt(iexp).firstRespCells = firstRespCells;
    end

    antiAnalysis = struct;
    antiAnalysis.longTC = [];
    antiAnalysis.longTCErr = [];
    antiAnalysis.cycTC = cell(1,nCycles);
    antiAnalysis.cycTCErr = cell(1,nCycles);
    for iexp = 1:nexp
        longTC = antiDataExpt(iexp).longTC;
        cycTC = antiDataExpt(iexp).cycTC(1:nCycles);

        antiAnalysis.longTC = cat(2,antiAnalysis.longTC,...
            mean(longTC,3));
        antiAnalysis.longTCErr = cat(2,antiAnalysis.longTCErr,...
            ste(longTC,3));
        antiAnalysis.cycTC = cellfun(@(x,y) cat(2,x,mean(y,3)),antiAnalysis.cycTC,...
            cycTC,'unif',0);
        antiAnalysis.cycTCErr = cellfun(@(x,y) cat(2,x,ste(y,3)),antiAnalysis.cycTCErr,...
            cycTC,'unif',0);
    end

    cellInfo = struct;
    cellInfo.firstRespCells = logical(cell2mat({respCellsExpt.firstRespCells}'));
    cellInfo.minRespCells = (mean(antiAnalysis.cycTC{1}(respwin,:),1) > ...
        minRespThreshold)';
    cellInfo.isTuned = logical(cell2mat({oriTuningExpt.isTuned}))';
    cellInfo.oriResp = cell2mat({oriTuningExpt.oriResp}');
    cellInfo.oriRespErr = cell2mat({oriTuningExpt.oriRespErr}');
    cellInfo.oriFit = cell2mat({oriTuningExpt.fit})';
    cellInfo.oriPref = cell2mat({oriTuningExpt.oriPref})';
    cellInfo.hwhm = hwhmFromOriFit(cellInfo.oriFit(:,1:180)',1:180)';

save(fullfile(fnout, [ds '_adaptAnalysis']),'antiAnalysis','cellInfo')
%% plotting params
respTCLim = [-0.005 0.05];
cycTCLim = [-0.005 0.05];
cycTCLim_minRespCells = [-0.005 0.025];
scatLim_win = [-0.2 0.6];
scatLim_cyc = [-0.035 0.085];
hmLim = [-0.1 0.1];
exCellTCLim = [-0.02 0.15];
oriRespLim = [-0.05 0.15];
siLim = [-10 10];
siOriLim = [-3 3];
oriBarLim_win = [0 0.08];
oriBarLim_resp = [0 0.04];
oriLim_taskResp = [-0.005 0.035];
oriNLim = [0 120];
oriTCLim = [-0.005 0.08];
targetTCLim = [-0.015 0.08];
outTCLim = [-0.005 0.04];
firstTCLim = [-0.005 0.04];
adaptLim = [0 1];
suppTCLim = [-0.05 0.005];
suppScatLim_win = [-0.2 0.1];
suppScatLim_cyc = [-0.015 0.015];

tcStartFrame = 26;
cycTCEndTimeMs = 350;
cycTCEndFr = 45;
% ttLabel_long = 0:500:2500;
ttLabel_cyc = -200:100:cycTCEndTimeMs;
ttLabel_target = -1000:250:900;
preTargetStimLabel = -700:350:0;
% nFr_long = size(antiAnalysis.longTC{1,1},1);
% tt_longTC = ((tcStartFrame:nFr_long)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
% ttLabelFr_long = ((ttLabel_long./1000)*frameRateHz)+...
%     ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_cyc = ((ttLabel_cyc./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
% ttLabelFr_target = ((ttLabel_target./1000)*frameRateHz)+...
%     ((nBaselineFr+nVisDelayFr_target)+1);

nFr_cyc = size(antiAnalysis.cycTC{1,1},1);
tt_cycTC = ((tcStartFrame:nFr_cyc)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
% tt_targetTC = ((1:nFr_cyc)-(nBaselineFr+nVisDelayFr_target)).*(1000/frameRateHz);

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

% lateCycRespAll = mean(antiAnalysis.lateCycTC{1}(respwin,:),1);

%%
ind = cellInfo.firstRespCells & cellInfo.minRespCells & cellInfo.isTuned;
cycResp = cellfun(@(x) mean(x(respwin,:),1) - mean(x(basewin_0,:),1),antiAnalysis.cycTC,'unif',0);

adaptResp = cellfun(@(x) x./cycResp{1},cycResp,'unif',0);

setFigParams4Print('landscape')
figure
suptitle([ds ': Tuned, Min. & First Stim Resp. Neurons'])
for icyc = 1:nCycles
    subplot(2,nCycles,icyc)
    bl = mean(antiAnalysis.cycTC{icyc}(basewin_0,ind),1);
    y = antiAnalysis.cycTC{icyc}(26:end,ind) - bl;
    yerr = ste(antiAnalysis.cycTC{icyc}(26:end,ind),2);
    shadedErrorBar_chooseColor(tt_cycTC,mean(y,2),yerr,[0 0 0]);
    figXAxis([],'Time from Stim (ms)',[tt_cycTC(1) 350])
    figYAxis([],'dF/F',cycTCLim)
    figAxForm
    hline(0,'k:')
    hold on
    vline(respWinTT,'k--');
    if icyc == 1
        title(sprintf('Stim #%s (%s/%s)',num2str(icyc),num2str(sum(ind)),...
            num2str(length(ind))))
    else
        title(sprintf('Stim #%s',num2str(icyc)))
    end
end

subplot 245
y = cellfun(@(x) mean(x(ind)),cycResp);
yerr = cellfun(@(x) ste(x(ind),2),cycResp);
errorbar(1:nCycles,y,yerr,'.')
figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'dF/F',cycTCLim)
figAxForm
title(sprintf('First Stim Resp. Cells (%s/%s)',num2str(sum(ind)),num2str(length(ind))))
subplot 246
y = cellfun(@(x) mean(x(ind)),adaptResp);
yerr = cellfun(@(x) ste(x(ind),2),adaptResp);
errorbar(1:nCycles,y,yerr,'.')
figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'Norm dF/F',[0 1])
figAxForm
title(sprintf('All Tuned Cells, n=%s', num2str(sum(ind))))
untunInd = cellInfo.firstRespCells & cellInfo.minRespCells & ~cellInfo.isTuned;
subplot 247
y = cellfun(@(x) mean(x(untunInd)),adaptResp);
yerr = cellfun(@(x) ste(x(untunInd),2),adaptResp);
errorbar(1:nCycles,y,yerr,'.')
figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'Norm dF/F',[0 1])
figAxForm
title(sprintf('All Untuned Cells, n=%s', num2str(sum(untunInd))))

oriColors = brewermap(5,'Blues');
oriColors = oriColors(2:end,:);
oriAdaptData = cell(nCycles,4);
oriAdaptGroup = cell(nCycles,4);
L = [];
for iori = 1:4
%     if iori ==5
%         indOri = ind;
%     else
        indOri = ind & cellInfo.oriPref == iori;
%     end
    
    subplot(2,4,8)
    y = cellfun(@(x) mean(x(indOri)),adaptResp);
    yerr = cellfun(@(x) ste(x(indOri),2),adaptResp);
    h = errorbar(1:nCycles,y,yerr,'.-');
    h.Color = oriColors(iori,:);
    L(iori) = h;
    hold on
    
    oriAdaptData(:,iori) = cellfun(@(x) x(indOri),adaptResp,'unif',0);
    oriAdaptGroup(:,iori) = repmat({ones(1,sum(indOri)).*iori},[nCycles 1]);
end
n = cellfun(@length,oriAdaptData(1,:));
leg = legend(L,cellfun(@(x,y) [x '-' y],...
    cellfun(@num2str,num2cell(orientations),'unif',0),...
    cellfun(@num2str,num2cell(n),'unif',0),'unif',0),...
    'location','northeast');
title(leg, 'Ori. Pref.')
figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'Norm dF/F',[0 1])
figAxForm

oriAdaptCycGroup = cell(nCycles,4);
for icyc = 1:nCycles
    oriAdaptCycGroup(icyc,:) = cellfun(@(x) ones(1,length(x))*icyc,oriAdaptData(icyc,:),'unif',0);
end
grp1 = cell2mat(oriAdaptGroup);
grp2 = cell2mat(oriAdaptCycGroup);
d = cell2mat(oriAdaptData);
[oriAdaptTest,~,stats] = anovan(d(:),{grp1(:);grp2(:)},'model','interaction');
title(sprintf('Two-way ANOVA: Ori p=%s, Stim # p=%s, OrixStim p=%s',...
    num2str(round(oriAdaptTest(1),2,'significant')),...
    num2str(round(oriAdaptTest(2),2,'significant')),...
    num2str(round(oriAdaptTest(3),2,'significant'))))

print(fullfile(fnout, [ds '_adaptation']),'-dpdf','-fillpage')
%% All resp cells
figure
suptitle([ds ': Min. & First Stim Resp. Neurons'])
ind = cellInfo.firstRespCells & cellInfo.minRespCells;
for icyc = 1:nCycles
    subplot(5,nCycles,icyc)
    bl = mean(antiAnalysis.cycTC{icyc}(basewin_0,ind),1);
    y = antiAnalysis.cycTC{icyc}(26:end,ind) - bl;
    yerr = ste(antiAnalysis.cycTC{icyc}(26:end,ind),2);
    shadedErrorBar_chooseColor(tt_cycTC,mean(y,2),yerr,[0 0 0]);
    figXAxis([],'Time from Stim (ms)',[tt_cycTC(1) 350])
    figYAxis([],'dF/F',cycTCLim)
    figAxForm
    hline(0,'k:')
    hold on
    vline(respWinTT,'k--');
    if icyc == 1
        title(sprintf('Stim #%s (%s/%s)',num2str(icyc),num2str(sum(ind)),...
            num2str(length(ind))))
    else
        title(sprintf('Stim #%s',num2str(icyc)))
    end
end

for iori = 1:4
    for icyc = 1:nCycles
        subplot(5,nCycles,icyc+iori*nCycles)
        indOri = ind & cellInfo.oriPref == iori;
        bl = mean(antiAnalysis.cycTC{icyc}(basewin_0,indOri),1);
        y = antiAnalysis.cycTC{icyc}(26:end,indOri) - bl;
        yerr = ste(antiAnalysis.cycTC{icyc}(26:end,indOri),2);
        shadedErrorBar_chooseColor(tt_cycTC,mean(y,2),yerr,[0 0 0]);
        figXAxis([],'Time from Stim (ms)',[tt_cycTC(1) 350])
        figYAxis([],'dF/F',cycTCLim)
        figAxForm
        hline(0,'k:')
        hold on
        vline(respWinTT,'k--');
        if icyc == 1
            title(sprintf('Stim #%s (%s/%s)',num2str(orientations(iori)),num2str(sum(indOri)),...
                num2str(length(indOri))))
        else
            title(sprintf('Stim #%s',num2str(icyc)))
        end
    end
end
print(fullfile(fnout, [ds '_adaptationByCycAndOri_all']),'-dpdf','-fillpage')

figure
subplot(1,2,1)
suptitle([ds ': Min. & First Stim Resp. Neurons'])
oriColors = brewermap(5,'Blues');
oriColors = oriColors(2:end,:);
oriAdaptData = cell(nCycles,4);
oriAdaptGroup = cell(nCycles,4);
L = [];
for iori = 1:4
    indOri = ind & cellInfo.oriPref == iori;
    y = cellfun(@(x) mean(x(indOri)),cycResp);
    yerr = cellfun(@(x) ste(x(indOri),2),cycResp);
    h = errorbar(1:nCycles,y,yerr,'.-');
    h.Color = oriColors(iori,:);
    L(iori) = h;
    hold on
    
    oriAdaptData(:,iori) = cellfun(@(x) x(indOri),cycResp,'unif',0);
    oriAdaptGroup(:,iori) = repmat({ones(1,sum(indOri)).*iori},[nCycles 1]);
end
n = cellfun(@length,oriAdaptData(1,:));
leg = legend(L,cellfun(@(x,y) [x '-' y],...
    cellfun(@num2str,num2cell(orientations),'unif',0),...
    cellfun(@num2str,num2cell(n),'unif',0),'unif',0),...
    'location','northeast');
title(leg, 'Ori. Pref.')
figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'dF/F',[0 0.05])
figAxForm

subplot(1,2,2)
oriColors = brewermap(5,'Blues');
oriColors = oriColors(2:end,:);
oriAdaptData = cell(nCycles,4);
oriAdaptGroup = cell(nCycles,4);
L = [];
for iori = 1:4
    indOri = ind & cellInfo.oriPref == iori;
    y = cellfun(@(x) mean(x(indOri)),adaptResp);
    yerr = cellfun(@(x) ste(x(indOri),2),adaptResp);
    h = errorbar(1:nCycles,y,yerr,'.-');
    h.Color = oriColors(iori,:);
    L(iori) = h;
    hold on
    
    oriAdaptData(:,iori) = cellfun(@(x) x(indOri),adaptResp,'unif',0);
    oriAdaptGroup(:,iori) = repmat({ones(1,sum(indOri)).*iori},[nCycles 1]);
end
n = cellfun(@length,oriAdaptData(1,:));
leg = legend(L,cellfun(@(x,y) [x '-' y],...
    cellfun(@num2str,num2cell(orientations),'unif',0),...
    cellfun(@num2str,num2cell(n),'unif',0),'unif',0),...
    'location','northeast');
title(leg, 'Ori. Pref.')
figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'Norm dF/F',[0 1])
figAxForm
print(fullfile(fnout, [ds '_adaptationByCycAndOri_AVG_all']),'-dpdf','-fillpage')

figure
for iori = 1:4
    for icyc = 1
        indOri = ind & cellInfo.oriPref == iori;
        bl = mean(antiAnalysis.cycTC{icyc}(basewin_0,indOri),1);
        y = antiAnalysis.cycTC{icyc}(26:end,indOri) - bl;
        plot(tt_cycTC,mean(y,2));
        figXAxis([],'Time from Stim (ms)',[tt_cycTC(1) 350])
        figYAxis([],'dF/F',cycTCLim)
        hold on
    end
end
title('Stim 1')
print(fullfile(fnout, [ds '_Stim1TCs_all']),'-dpdf','-fillpage')
%% Only tuned cells
figure
suptitle([ds ': Tuned, Min. & First Stim Resp. Neurons'])
ind = cellInfo.firstRespCells & cellInfo.minRespCells & cellInfo.isTuned;
for icyc = 1:nCycles
    subplot(5,nCycles,icyc)
    bl = mean(antiAnalysis.cycTC{icyc}(basewin_0,ind),1);
    y = antiAnalysis.cycTC{icyc}(26:end,ind) - bl;
    yerr = ste(antiAnalysis.cycTC{icyc}(26:end,ind),2);
    shadedErrorBar_chooseColor(tt_cycTC,mean(y,2),yerr,[0 0 0]);
    figXAxis([],'Time from Stim (ms)',[tt_cycTC(1) 350])
    figYAxis([],'dF/F',cycTCLim)
    figAxForm
    hline(0,'k:')
    hold on
    vline(respWinTT,'k--');
    if icyc == 1
        title(sprintf('Stim #%s (%s/%s)',num2str(icyc),num2str(sum(ind)),...
            num2str(length(ind))))
    else
        title(sprintf('Stim #%s',num2str(icyc)))
    end
end

for iori = 1:4
    for icyc = 1:nCycles
        subplot(5,nCycles,icyc+iori*nCycles)
        indOri = ind & cellInfo.oriPref == iori;
        bl = mean(antiAnalysis.cycTC{icyc}(basewin_0,indOri),1);
        y = antiAnalysis.cycTC{icyc}(26:end,indOri) - bl;
        yerr = ste(antiAnalysis.cycTC{icyc}(26:end,indOri),2);
        shadedErrorBar_chooseColor(tt_cycTC,mean(y,2),yerr,[0 0 0]);
        figXAxis([],'Time from Stim (ms)',[tt_cycTC(1) 350])
        figYAxis([],'dF/F',cycTCLim)
        figAxForm
        hline(0,'k:')
        hold on
        vline(respWinTT,'k--');
        if icyc == 1
            title(sprintf('Stim #%s (%s/%s)',num2str(orientations(iori)),num2str(sum(indOri)),...
                num2str(length(indOri))))
        else
            title(sprintf('Stim #%s',num2str(icyc)))
        end
    end
end
print(fullfile(fnout, [ds '_adaptationByCycAndOri_tuned']),'-dpdf','-fillpage')

figure
suptitle([ds ': Tuned, Min. & First Stim Resp. Neurons'])
subplot(1,2,1)
oriColors = brewermap(5,'Blues');
oriColors = oriColors(2:end,:);
oriAdaptData = cell(nCycles,4);
oriAdaptGroup = cell(nCycles,4);
L = [];
for iori = 1:4
    indOri = ind & cellInfo.oriPref == iori;
    y = cellfun(@(x) mean(x(indOri)),cycResp);
    yerr = cellfun(@(x) ste(x(indOri),2),cycResp);
    h = errorbar(1:nCycles,y,yerr,'.-');
    h.Color = oriColors(iori,:);
    L(iori) = h;
    hold on
    
    oriAdaptData(:,iori) = cellfun(@(x) x(indOri),cycResp,'unif',0);
    oriAdaptGroup(:,iori) = repmat({ones(1,sum(indOri)).*iori},[nCycles 1]);
end
n = cellfun(@length,oriAdaptData(1,:));
leg = legend(L,cellfun(@(x,y) [x '-' y],...
    cellfun(@num2str,num2cell(orientations),'unif',0),...
    cellfun(@num2str,num2cell(n),'unif',0),'unif',0),...
    'location','northeast');
title(leg, 'Ori. Pref.')
figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'dF/F',[0 0.05])
figAxForm

subplot(1,2,2)
oriColors = brewermap(5,'Blues');
oriColors = oriColors(2:end,:);
oriAdaptData = cell(nCycles,4);
oriAdaptGroup = cell(nCycles,4);
L = [];
for iori = 1:4
    indOri = ind & cellInfo.oriPref == iori;
    y = cellfun(@(x) mean(x(indOri)),adaptResp);
    yerr = cellfun(@(x) ste(x(indOri),2),adaptResp);
    h = errorbar(1:nCycles,y,yerr,'.-');
    h.Color = oriColors(iori,:);
    L(iori) = h;
    hold on
    
    oriAdaptData(:,iori) = cellfun(@(x) x(indOri),adaptResp,'unif',0);
    oriAdaptGroup(:,iori) = repmat({ones(1,sum(indOri)).*iori},[nCycles 1]);
end
n = cellfun(@length,oriAdaptData(1,:));
leg = legend(L,cellfun(@(x,y) [x '-' y],...
    cellfun(@num2str,num2cell(orientations),'unif',0),...
    cellfun(@num2str,num2cell(n),'unif',0),'unif',0),...
    'location','northeast');
title(leg, 'Ori. Pref.')
figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'Norm dF/F',[0 1])
figAxForm
print(fullfile(fnout, [ds '_adaptationByCycAndOri_AVG_tuned']),'-dpdf','-fillpage')
%% 
if ~bxExpt
    bxData = load(fullfile(rc.ashleyAnalysis, 'Expt summaries','adaptation','behavior',...
        ['attentionV1_adaptAnalysis']));
    bxCycResp = cellfun(@(x) mean(x(respwin,:),1) - mean(x(basewin_0,:),1),bxData.antiAnalysis.cycTC,'unif',0);
    bxAdaptResp = cellfun(@(x) x./bxCycResp{1},bxCycResp,'unif',0);
    bxInd = bxData.cellInfo.firstRespCells & bxData.cellInfo.minRespCells & bxData.cellInfo.isTuned;

    
    
    figure    
    subplot 121
    y = cellfun(@(x) mean(x(ind)),adaptResp);
    yerr = cellfun(@(x) ste(x(ind),2),adaptResp);
    errorbar(1:nCycles,y,yerr,'.')
    hold on
    y = cellfun(@(x) mean(x(bxInd)),bxAdaptResp);
    yerr = cellfun(@(x) ste(x(bxInd),2),bxAdaptResp);
    errorbar(1:nCycles,y,yerr,'.')
    hold on
    figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
    figYAxis([],'Norm dF/F',[0 1])
    figAxForm
    title(sprintf('Naive: n=%s, Behav: n=%s', num2str(sum(ind)),...
        num2str(sum(bxInd))))
    legend({'Naive';'Behavior'})
    
    subplot 122
    [~,x] = max(cellInfo.oriFit(ind,:),[],2);
    y = adaptResp{4}(ind);
    y(y<0) = 0;
    plot(x,y,'.','MarkerSize',10)
    hold on
    [~,x] = max(bxData.cellInfo.oriFit(bxInd,:),[],2);
    y = bxAdaptResp{4}(bxInd);
    y(y<0) = 0;
    plot(x,y,'.','MarkerSize',10)
    figXAxis([],'Ori. Pref.',[-1 182])
    figYAxis([],'Norm. dF/F',[0 1]);
    figAxForm
    
    print(fullfile(rc.ashleyAnalysis, 'Expt summaries',...
        'adaptation','adaptAnalysisBehaviorVsNaive'),'-dpdf','-fillpage')
end