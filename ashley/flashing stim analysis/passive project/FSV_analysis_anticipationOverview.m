
nMice = length(mice);
channelColor = {'green','red'};
avName = {'visual','auditory'};

%%

%% plotting params
respTCLim = [-0.005 0.05];
cycTCLim = [-0.005 0.015];
eaCycTCLim = [-0.005 0.05];
cycTCLim_minRespCells = [-0.005 0.025];
scatLim_win = [-0.2 0.6];
scatLim_cyc = [-0.035 0.085];
hmLim = [-0.2 0.2];
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

offset_longBL = 10;
tcStartFrame = 26;
tcStartFrame_longBL = tcStartFrame-offset_longBL;
tcEndFrame = 45;
cycTCEndTimeMs = 350;
cycTCEndFr = 45;
ttLabel_long = 0:500:2500;
ttLabel_cyc = -200:100:cycTCEndTimeMs;
ttLabel_target = -1000:250:900;
preTargetStimLabel = -700:350:0;
nFr_long = 118;
tt_longTC = ((tcStartFrame:nFr_long)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
tt_longTC_longBL = ((tcStartFrame_longBL:nFr_long)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
ttLabelFr_long = ((ttLabel_long./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_cyc = ((ttLabel_cyc./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_target = ((ttLabel_target./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr_target)+1);

nFr_cyc = 45;
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

%% grab SOM data
nexp = length(som_expt);
expt_somExptOnly = struct;
expt_sortID = nan(1,nexp);
for im = 1:nMice
    ne = size(mouse(im).expt,2);
    if any(strcmp({som_expt.SubNum},mouse(im).mouse_name))
        for iexp = 1:ne
            if isempty(fieldnames(expt_somExptOnly))
                exptN = 1;
            else
                exptN = exptN+1;
            end
            d = mouse(im).expt(iexp);
            expt_somExptOnly(exptN).expt_name = [mouse(im).mouse_name '-' d.date];
            ind = strcmp({expt.SubNum},mouse(im).mouse_name) & ...
                strcmp({expt.date},d.date);
            if expt(ind).indicator{2}(1:4) == 'flex'
                expt_somExptOnly(exptN).indicator = expt(ind).indicator{2}(6:end);
            else
                expt_somExptOnly(exptN).indicator = expt(ind).indicator{2};
            end
            expt_somExptOnly(exptN).isBx = logical(expt(ind).isBehav);
            expt_sortID(exptN) = find(ind);
            
            if ~isempty(expt(ind).passExpt)
                mouseInd = strcmp({mousePass.mouse_name},mouse(im).mouse_name);
                exptInd = strcmp({mousePass(mouseInd).expt.date},d.date);
                d_pass = mousePass(mouseInd).expt(exptInd);
            else
                d_pass = [];
            end
            for itag = 1:2
                if itag == 1
                    tag_name = expt(ind).greenChannelLabel;
                    if strcmp(tag_name,'SOM')
                        somIsGreen = true;
                    else
                        somIsGreen = false;
                    end
                else
                    tag_name = expt(ind).redChannelLabel;
                    if strcmp(tag_name,'SOM')
                        somIsGreen = false;
                    end
                end
                expt_somExptOnly(exptN).tag(itag).name = tag_name;
                if ~isempty(tag_name)
                    for iav = 1:2
                        if somIsGreen                            
                            expt_somExptOnly(exptN).tag(1).av(iav) = d.tag(2).av(iav);
                            if ~isempty(d_pass)
                                expt_somExptOnly(exptN).tag(1).av_pass(iav) = ...
                                    d_pass.tag(2).av(iav);
                            else
                                expt_somExptOnly(exptN).tag(1).av_pass = [];
                            end
                        else
                            expt_somExptOnly(exptN).tag(itag).av(iav) = d.tag(itag).av(iav);
                            if ~isempty(d_pass)
                                expt_somExptOnly(exptN).tag(itag).av_pass(iav) = ...
                                    d_pass.tag(itag).av(iav);
                                fprintf([num2str(exptN) '\n'])
                            else
                                expt_somExptOnly(exptN).tag(itag).av_pass = [];
                            end
                        end
                    end
                end

            end
        end
    end
end

indicators = unique({expt_somExptOnly.indicator});

% long TC and each stimulus response
cycLengthFr = 11;
resp_SOM = struct;
resp_Other = struct;
for iexp = 1:nexp
    if expt_somExptOnly(iexp).tag(1).name == 'SOM'
        somInd = 1;
        otherInd = 2;
    else
        somInd = 2;
        otherInd = 1;
    end
    % SOM neuron data
    resp_SOM(iexp).isBx = expt_somExptOnly(iexp).isBx;
    resp_SOM(iexp).channel = channelColor{somInd};
    resp_SOM(iexp).indicator = expt_somExptOnly(iexp).indicator;
    if ~isempty(expt_somExptOnly(iexp).tag(somInd).av(1).align(1).respTC)
        for iav = 1:2
            d = expt_somExptOnly(iexp).tag(somInd).av(iav).align(1);
            ind_outcome = strcmp(d.outcome,'success')|strcmp(d.outcome,'ignore');
            ind_long = d.nCycles >= max(lateCycles);
            resp_SOM(iexp).av(iav).long = d.respTC(:,:,ind_long) - ...
                mean(d.respTC(basewin,:,ind_long),1);
            cycTC = cell(1,maxCycles);
            lateCycTC = [];
            for icyc = 1:maxCycles
                tc = d.respTC(:,:,d.nCycles >= icyc & ind_outcome);
                cycStartOffset = ((icyc-1).*cycLengthFr)+nBaselineFr;
                cycTC{1,icyc} = tc(...
                    (cycStartOffset-nBaselineFr+1):(cycStartOffset+nFrames1s),:,:);
                if icyc >= min(lateCycles)
                    lateCycTC = cat(3,lateCycTC,cycTC{1,icyc});
                end
            end
            resp_SOM(iexp).av(iav).cycTC = cycTC;
            resp_SOM(iexp).av(iav).lateCycTC = lateCycTC;
            if isfield(expt_somExptOnly(iexp).tag(somInd),'av_pass')
                fprintf('.\n')
                if ~isempty(expt_somExptOnly(iexp).tag(somInd).av_pass)
                    d = expt_somExptOnly(iexp).tag(somInd).av_pass(iav).align(1);
                    ind_long = d.nCycles >= max(lateCycles);
                    resp_SOM(iexp).av_pass(iav).long = d.respTC(:,:,ind_long) - ...
                        mean(d.respTC(basewin,:,ind_long),1);
                end
            end
        end
    end
    % Other neuron data
    resp_Other(iexp).isBx = expt_somExptOnly(iexp).isBx;
    resp_Other(iexp).channel = channelColor{otherInd};
    resp_Other(iexp).indicator = expt_somExptOnly(iexp).indicator;
    if ~isempty(expt_somExptOnly(iexp).tag(otherInd).av)
        for iav = 1:2
            d = expt_somExptOnly(iexp).tag(otherInd).av(iav).align(1);
            ind_outcome = strcmp(d.outcome,'success')|strcmp(d.outcome,'ignore');
            ind_long = d.nCycles >= max(lateCycles);
            resp_Other(iexp).av(iav).long = d.respTC(:,:,ind_long) - ...
                mean(d.respTC(basewin,:,ind_long),1);
            cycTC = cell(1,maxCycles);
            lateCycTC = [];
            for icyc = 1:maxCycles
                tc = d.respTC(:,:,d.nCycles >= icyc & ind_outcome);
                cycStartOffset = ((icyc-1).*cycLengthFr)+nBaselineFr;
                cycTC{1,icyc} = tc(...
                    (cycStartOffset-nBaselineFr+1):(cycStartOffset+nFrames1s),:,:);
                if icyc >= min(lateCycles)
                    lateCycTC = cat(3,lateCycTC,cycTC{1,icyc});
                end
            end
            resp_Other(iexp).av(iav).cycTC = cycTC;
            resp_Other(iexp).av(iav).lateCycTC = lateCycTC;
            if isfield(expt_somExptOnly(iexp).tag(otherInd),'av_pass')
                if ~isempty(expt_somExptOnly(iexp).tag(otherInd).av_pass)
                    d = expt_somExptOnly(iexp).tag(otherInd).av_pass(iav).align(1);
                    ind_long = d.nCycles >= max(lateCycles);
                    resp_Other(iexp).av_pass(iav).long = d.respTC(:,:,ind_long) - ...
                        mean(d.respTC(basewin,:,ind_long),1);
                end
            end
        end
    end
end

cellInfo_SOM = struct;
cellInfo_Other = struct;
for iexp = 1:nexp
    d_av_long = cat(3,resp_SOM(iexp).av(visualTrials).long,...
        resp_SOM(iexp).av(auditoryTrials).long);
    d_av_first = cat(3,resp_SOM(iexp).av(visualTrials).cycTC{1},...
        resp_SOM(iexp).av(auditoryTrials).cycTC{1});
    d_av_late = cat(3,resp_SOM(iexp).av(visualTrials).lateCycTC,...
        resp_SOM(iexp).av(auditoryTrials).lateCycTC);
    cellInfo_SOM(iexp).responsive_first = ttest(squeeze(mean(d_av_first(respwin,:,:))),...
        squeeze(mean(d_av_first(basewin,:,:))),'dim',2,'tail','right');
    cellInfo_SOM(iexp).responsive_late = ttest(squeeze(mean(d_av_late(respwin,:,:))),...
        squeeze(mean(d_av_late(basewin,:,:))),'dim',2,'tail','right');
    cellInfo_SOM(iexp).responsive_lateWindow = ttest(squeeze(mean(...
        d_av_long((cycLengthFr*length(lateCycles)+1):end,:,:))),...
        squeeze(mean(d_av_long(basewin,:,:))),'dim',2,'tail','right');
    cellInfo_SOM(iexp).suppressed_lateWindow = ttest(squeeze(mean(...
        d_av_long((cycLengthFr*length(lateCycles)+1):end,:,:))),...
        squeeze(mean(d_av_long(basewin,:,:))),'dim',2,'tail','left');
    if ~isempty(resp_Other(iexp).av)
        d_av_long = cat(3,resp_Other(iexp).av(visualTrials).long,...
            resp_Other(iexp).av(auditoryTrials).long);
        d_av_first = cat(3,resp_Other(iexp).av(visualTrials).cycTC{1},...
            resp_Other(iexp).av(auditoryTrials).cycTC{1});
        d_av_late = cat(3,resp_Other(iexp).av(visualTrials).lateCycTC,...
            resp_Other(iexp).av(auditoryTrials).lateCycTC);
        cellInfo_Other(iexp).responsive_first = ttest(squeeze(mean(d_av_first(respwin,:,:))),...
            squeeze(mean(d_av_first(basewin,:,:))),'dim',2,'tail','right');
        cellInfo_Other(iexp).responsive_late = ttest(squeeze(mean(d_av_late(respwin,:,:))),...
            squeeze(mean(d_av_late(basewin,:,:))),'dim',2,'tail','right');
        cellInfo_Other(iexp).responsive_lateWindow = ttest(squeeze(mean(...
            d_av_long((cycLengthFr*length(lateCycles)+1):end,:,:))),...
            squeeze(mean(d_av_long(basewin,:,:))),'dim',2,'tail','right');
        cellInfo_Other(iexp).suppressed_lateWindow = ttest(squeeze(mean(...
            d_av_long((cycLengthFr*length(lateCycles)+1):end,:,:))),...
            squeeze(mean(d_av_long(basewin,:,:))),'dim',2,'tail','left');  
    else
        cellInfo_Other(iexp).responsive_first = [];
    end
end

%% plot heatmaps for naive and training data.
bxInd = cell2mat({resp_SOM.isBx}) == 1;

figure
colormap(brewermap([],'*RdBu'));

d_bx = [];
d_pass = [];
d_naive = [];
for iexp = 1:nexp
    if bxInd(iexp)
        ind = cellInfo_SOM(iexp).responsive_first | ...
            cellInfo_SOM(iexp).responsive_late | ...
            cellInfo_SOM(iexp).responsive_lateWindow | ...
            cellInfo_SOM(iexp).suppressed_lateWindow;
        d_bx = cat(2,d_bx,mean(resp_SOM(iexp).av(visualTrials).long(:,ind,:),3));
        if isfield(resp_SOM(iexp),'av_pass')
            if ~isempty(resp_SOM(iexp).av_pass)
                d_pass = cat(2,d_pass,...
                    mean(resp_SOM(iexp).av_pass(visualTrials).long(:,ind,:),3));
            end
        end
    else
        ind = cellInfo_SOM(iexp).responsive_first | ...
            cellInfo_SOM(iexp).responsive_late | ...
            cellInfo_SOM(iexp).responsive_lateWindow | ...
            cellInfo_SOM(iexp).suppressed_lateWindow;
        d_naive = cat(2,d_naive,mean(resp_SOM(iexp).av(visualTrials).long(:,ind,:),3));
    end
end

subplot 231
[~,sortInd] = sort(mean(d_bx(lateWinFr,:),1));
hm = flipud(d_bx(:,sortInd)');
imagesc(hm(:,tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr)))
hold on
figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title(sprintf('SOM+ - Behav, n=%s',num2str(size(hm,1))))

subplot 232
[~,sortInd] = sort(mean(d_pass(lateWinFr,:),1));
hm = flipud(d_pass(:,sortInd)');
imagesc(hm(:,tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr)))
hold on
figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title(sprintf('SOM+ - Passive, n=%s',num2str(size(hm,1))))

subplot 233
[~,sortInd] = sort(mean(d_naive(lateWinFr,:),1));
hm = flipud(d_naive(:,sortInd)');
imagesc(hm(:,tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr)))
hold on
figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title(sprintf('SOM+ - Naive, n=%s',num2str(size(hm,1))))

d_bx = [];
d_naive = [];
d_pass = [];
for iexp = 1:nexp
    if bxInd(iexp)
        ind = cellInfo_Other(iexp).responsive_first | ...
            cellInfo_Other(iexp).responsive_late | ...
            cellInfo_Other(iexp).responsive_lateWindow | ...
            cellInfo_Other(iexp).suppressed_lateWindow;
        if ~isempty(resp_Other(iexp).av)
            d_bx = cat(2,d_bx,mean(resp_Other(iexp).av(visualTrials).long(:,ind,:),3));
            if isfield(resp_Other(iexp),'av_pass')
                if ~isempty(resp_Other(iexp).av_pass)
                    d_pass = cat(2,d_pass,...
                        mean(resp_Other(iexp).av_pass(visualTrials).long(:,ind,:),3));
                end
            end
        end
    else
        ind = cellInfo_Other(iexp).responsive_first | ...
            cellInfo_Other(iexp).responsive_late | ...
            cellInfo_Other(iexp).responsive_lateWindow | ...
            cellInfo_Other(iexp).suppressed_lateWindow;
        if ~isempty(resp_Other(iexp).av)
            d_naive = cat(2,d_naive,...
                mean(resp_Other(iexp).av(visualTrials).long(:,ind,:),3));
        end
    end
end

subplot 234
[~,sortInd] = sort(mean(d_bx(lateWinFr,:),1));
hm = flipud(d_bx(:,sortInd)');
imagesc(hm(:,tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr)))
hold on
figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title(sprintf('Other - Behav, n=%s',num2str(size(hm,1))))

subplot 235
[~,sortInd] = sort(mean(d_pass(lateWinFr,:),1));
hm = flipud(d_pass(:,sortInd)');
imagesc(hm(:,tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr)))
hold on
figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title(sprintf('Other - Behav, n=%s',num2str(size(hm,1))))

subplot 236
[~,sortInd] = sort(mean(d_naive(lateWinFr,:),1));
hm = flipud(d_naive(:,sortInd)');
imagesc(hm(:,tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr)))
hold on
figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title(sprintf('Other - Behav, n=%s',num2str(size(hm,1))))

print([fnout 'suppressedHeatmaps'],'-dpdf','-fillpage')
%% plot TC of suppressed cells
bxInd = cell2mat({resp_SOM.isBx}) == 1;
engage_label = {'Behav','Passive','Naive'};
engage_colors = {[0 0.7 0],[0 0 1],[0 0 0]};
tcSOM = cell(1,3);
tcOther = cell(1,3);
nSOM = zeros(1,3);
nOther = zeros(1,3);
for iexp = 1:nexp
    indSOM = cellInfo_SOM(iexp).suppressed_lateWindow == 1;
    indOther = cellInfo_Other(iexp).suppressed_lateWindow == 1;
    if bxInd(iexp)
        tcSOM{1} = cat(2,tcSOM{1},...
            mean(resp_SOM(iexp).av(visualTrials).long(:,indSOM,:),3));
        nSOM(1) = nSOM(1)+length(indSOM);
        if ~isempty(resp_Other(iexp).av)
            tcOther{1} = cat(2,tcOther{1},...
                mean(resp_Other(iexp).av(visualTrials).long(:,indOther,:),3));
            nOther(1) = nOther(1)+length(indOther);
        end
        if isfield(resp_SOM(iexp),'av_pass')
            if ~isempty(resp_SOM(iexp).av_pass)
                tcSOM{2} = cat(2,tcSOM{2},...
                    mean(resp_SOM(iexp).av_pass(visualTrials).long(:,indSOM,:),3));
                nSOM(2) = nSOM(2)+length(indSOM);
                tcOther{2} = cat(2,tcOther{2},...
                    mean(resp_Other(iexp).av_pass(visualTrials).long(:,indOther,:),3));
                nOther(2) = nOther(2)+length(indOther);
            end
        end
    else
        tcSOM{3} = cat(2,tcSOM{3},...
            mean(resp_SOM(iexp).av(visualTrials).long(:,indSOM,:),3));
        nSOM(3) = nSOM(3)+length(indSOM);        
        if ~isempty(resp_Other(iexp).av)
            tcOther{3} = cat(2,tcOther{3},...
                mean(resp_Other(iexp).av(visualTrials).long(:,indOther,:),3));
            nOther(3) = nOther(3)+length(indOther);
        end
    end
end

figure
subplot 211
leg_lines = [];
leg_label = cell(1,3);
for i = 1:3
    y = mean(tcSOM{i}(tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr),:),2);
    yerr = ste(tcSOM{i}(tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr),:),2);
    hold on
    h = shadedErrorBar_chooseColor(tt_longTC,y,yerr,engage_colors{i});
    leg_lines(i) = h.mainLine;
    leg_label{i} = sprintf('%s: n=%s/%s',engage_label{i},...
        num2str(size(tcSOM{i},2)),num2str(nSOM(i)));
end
legend(leg_lines,leg_label,'location','southwest')
figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
figYAxis([],'dF/F',[])
figAxForm([],0)
title('SOM+')
hline(0,'k:')

subplot 212
leg_lines = [];
leg_label = cell(1,3);
for i = 1:3
    y = mean(tcOther{i}(tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr),:),2);
    yerr = ste(tcOther{i}(tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr),:),2);
    hold on
    h = shadedErrorBar_chooseColor(tt_longTC,y,yerr,engage_colors{i});
    leg_lines(i) = h.mainLine;
    leg_label{i} = sprintf('%s: n=%s/%s',engage_label{i},...
        num2str(size(tcOther{i},2)),num2str(nSOM(i)));
end
legend(leg_lines,leg_label,'location','southwest')
figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
figYAxis([],'dF/F',[])
figAxForm([],0)
title('Other')
hline(0,'k:')

print([fnout 'suppressedTC'],'-dpdf','-fillpage')

%% behaving versus passive - same cells

bxInd = cell2mat({resp_SOM.isBx}) == 1;
passInd = false(1,nexp);
for iexp = 1:nexp
    if isfield(resp_SOM(iexp),'av_pass')
        if ~isempty(resp_SOM(iexp).av_pass)
            passInd(iexp) = true;
        end
    end
end

% Each experiment
setFigParams4Print('landscape')
for iexp = 1:nexp
    if bxInd(iexp) & passInd(iexp)
        ind = cellInfo_SOM(iexp).responsive_first | ...
            cellInfo_SOM(iexp).responsive_late | ...
            cellInfo_SOM(iexp).responsive_lateWindow | ...
            cellInfo_SOM(iexp).suppressed_lateWindow;
        ind2 = cellInfo_SOM(iexp).suppressed_lateWindow(ind) == 1;
        d_bx = mean(resp_SOM(iexp).av(visualTrials).long(:,ind,:),3);
        d_pass = mean(resp_SOM(iexp).av_pass(visualTrials).long(:,ind,:),3);
        
        figure
        suptitle({[expt_somExptOnly(iexp).expt_name '-' resp_SOM(iexp).indicator];...
            'cells sorted same across behav & pass';
            sprintf('n=%s active cells of %s ',...
            num2str(sum(ind)),num2str(length(ind)))})
        colormap(brewermap([],'*RdBu'));
        
        subplot 231
        [~,sortInd] = sort(mean(d_bx(lateWinFr,:),1));
        hm = flipud(d_bx(:,sortInd)');
        imagesc(hm(:,tcStartFrame_longBL:((lateCycles(end).*cycLengthFr)+nBaselineFr)))
        hold on
        figXAxis([],'Time from Start (ms)',[],ttLabelFr_long+offset_longBL,ttLabel_long)
        figYAxis([],'Cell #',[])
        figAxForm
        colorbar
        caxis(hmLim)
        title(sprintf('SOM+ - Behav, n=%s',num2str(size(hm,1))))
        
        subplot 232
        hm = flipud(d_pass(:,sortInd)');
        imagesc(hm(:,tcStartFrame_longBL:((lateCycles(end).*cycLengthFr)+nBaselineFr)))
        hold on
        figXAxis([],'Time from Start (ms)',[],ttLabelFr_long+offset_longBL,ttLabel_long)
        figYAxis([],'Cell #',[])
        figAxForm
        colorbar
        caxis(hmLim)
        title(sprintf('SOM+ - Passive, n=%s',num2str(size(hm,1))))
        
        subplot 233
        if any(ind2)
            leg_lines = [];
            leg_label = cell(1,2);
            tc = {d_bx(:,ind2), d_pass(:,ind2)};
            for i = 1:2
                y = mean(tc{i}(tcStartFrame_longBL:((lateCycles(end).*cycLengthFr)+nBaselineFr),:),2);
                yerr = ste(tc{i}(tcStartFrame_longBL:((lateCycles(end).*cycLengthFr)+nBaselineFr),:),2);
                hold on
                h = shadedErrorBar_chooseColor(tt_longTC_longBL,y,yerr,engage_colors{i});
                leg_lines(i) = h.mainLine;
                leg_label{i} = sprintf('%s',engage_label{i});
            end
            legend(leg_lines,leg_label,'location','southwest')
            figXAxis([],'Time (ms)',[tt_longTC_longBL(1) tt_longTC_longBL(end)],ttLabel_long,ttLabel_long)
            figYAxis([],'dF/F',[])
            figAxForm([],0)
            title(sprintf('Supp Cells, n=%s',num2str(sum(ind2))))
            hline(0,'k:')
        else
            title('No Supp Cells')
        end
        
        ind = cellInfo_Other(iexp).responsive_first | ...
            cellInfo_Other(iexp).responsive_late | ...
            cellInfo_Other(iexp).responsive_lateWindow | ...
            cellInfo_Other(iexp).suppressed_lateWindow;
        ind2 = cellInfo_Other(iexp).suppressed_lateWindow(ind) == 1;
        d_bx = mean(resp_Other(iexp).av(visualTrials).long(:,ind,:),3);
        d_pass = mean(resp_Other(iexp).av_pass(visualTrials).long(:,ind,:),3);
        
        subplot 234
        [~,sortInd] = sort(mean(d_bx(lateWinFr,:),1));
        hm = flipud(d_bx(:,sortInd)');
        imagesc(hm(:,tcStartFrame_longBL:((lateCycles(end).*cycLengthFr)+nBaselineFr)))
        hold on
        figXAxis([],'Time from Start (ms)',[],ttLabelFr_long+offset_longBL,ttLabel_long)
        figYAxis([],'Cell #',[])
        figAxForm
        colorbar
        caxis(hmLim)
        title(sprintf('Other - Behav, n=%s',num2str(size(hm,1))))
        
        subplot 235
        hm = flipud(d_pass(:,sortInd)');
        imagesc(hm(:,tcStartFrame_longBL:((lateCycles(end).*cycLengthFr)+nBaselineFr)))
        hold on
        figXAxis([],'Time from Start (ms)',[],ttLabelFr_long+offset_longBL,ttLabel_long)
        figYAxis([],'Cell #',[])
        figAxForm
        colorbar
        caxis(hmLim)
        title(sprintf('Other - Passive, n=%s',num2str(size(hm,1))))
                
        subplot 236
        if any(ind2)
            leg_lines = [];
            leg_label = cell(1,2);
            tc = {d_bx(:,ind2); d_pass(:,ind2)};
            for i = 1:2
                y = mean(tc{i}(tcStartFrame_longBL:((lateCycles(end).*cycLengthFr)+nBaselineFr),:),2);
                yerr = ste(tc{i}(tcStartFrame_longBL:((lateCycles(end).*cycLengthFr)+nBaselineFr),:),2);
                hold on
                h = shadedErrorBar_chooseColor(tt_longTC_longBL,y,yerr,engage_colors{i});
                leg_lines(i) = h.mainLine;
                leg_label{i} = sprintf('%s',engage_label{i});
            end
            legend(leg_lines,leg_label,'location','southwest')
            figXAxis([],'Time (ms)',[tt_longTC_longBL(1) tt_longTC_longBL(end)],ttLabel_long,ttLabel_long)
            figYAxis([],'dF/F',[])
            figAxForm([],0)
            title(sprintf('Supp Cells, n=%s',num2str(sum(ind2))))
            hline(0,'k:')
        else
            title('No Supp Cells')
        end
        if iexp == 20
            print([fnout 'exExpt_hm_suppCellsTC'],'-dpdf','-fillpage')
        end
    end
end

% across experiments
colors_expt = brewermap(sum(bxInd&passInd),'Dark2');
figure
suptitle('late window response, suppressed cells')
for iexp = 1:nexp
    if iexp == 1
        exptN = 0;
    end
    if bxInd(iexp) & passInd(iexp)
        exptN = exptN+1;
        subplot 221
        hold on
        ind = cellInfo_SOM(iexp).suppressed_lateWindow == 1;
        x =  mean(mean(resp_SOM(iexp).av(visualTrials).long(lateWinFr,ind,:),3),1);
        y =  mean(mean(resp_SOM(iexp).av_pass(visualTrials).long(lateWinFr,ind,:),3),1);
        if any(ind)
            h = plot(x,y,'.','MarkerSize',10);
            h.Color = colors_expt(exptN,:);
            subplot 222
            hold on
            h = errorbar(1,mean(x-y),ste(x-y,2),'.','MarkerSize',20);
            h.Color = (colors_expt(exptN,:));
        end
        
        subplot 223
        hold on
        ind = cellInfo_Other(iexp).suppressed_lateWindow == 1;
        x =  mean(mean(resp_Other(iexp).av(visualTrials).long(lateWinFr,ind,:),3),1);
        y =  mean(mean(resp_Other(iexp).av_pass(visualTrials).long(lateWinFr,ind,:),3),1);
        if any(ind)
            h = plot(x,y,'.','MarkerSize',10);
            h.Color = colors_expt(exptN,:);
            subplot 224
            hold on
            h = errorbar(1,mean(x-y),ste(x-y,2),'.','MarkerSize',20);
            h.Color = (colors_expt(exptN,:));
        end
        
    end
end
plotInd = [1,3];
suppRespLim = [-0.15 0.15];
for i = 1:2
    subplot(2,2,plotInd(i))
    figXAxis([],'Behavior',suppRespLim)
    figYAxis([],'Passive',suppRespLim)
    figAxForm
    plot(suppRespLim,suppRespLim,'k--')
    if i == 1
        title('SOM+')
    else
        title('Other')
    end
end
plotInd = [2,4];
suppRespLim = [-0.15 0.15];
for i = 1:2
    subplot(2,2,plotInd(i))
    figXAxis([],'',[0 2],1,'Passive')
    ax = gca;
    ax.XTickLabelRotation = -45;
    figYAxis([],'Behavior-Passive',[-0.04 0.01])
    figAxForm
    if i == 1
        title('SOM+')
    else
        title('Other')
    end
end
  
print([fnout 'lateWin_BxVsPass_suppCells'],'-dpdf','-fillpage')
    