
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
exptIDs_mouse_date = cell(2,nexp);
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
            exptIDs_mouse_date{1,exptN} = mouse(im).mouse_name;
            exptIDs_mouse_date{2,exptN} = d.date;
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
%%

bxInd = cell2mat({resp_SOM.isBx}) == 1;
passInd = false(1,nexp);
for iexp = 1:nexp
    if isfield(resp_SOM(iexp),'av_pass')
        if ~isempty(resp_SOM(iexp).av_pass)
            passInd(iexp) = true;
        end
    end
end

%% grab FOV masks
respLim = [-0.05 0.05];
iexp=18;
masks = load(fullfile(rc.ashleyAnalysis,exptIDs_mouse_date{1,iexp},...
    'two-photon imaging',exptIDs_mouse_date{2,iexp},'data processing',...
    'final_mask_cells.mat'));

n_red = length(unique(masks.red_mask_cell(:)))-1;
n_green = length(unique(masks.green_mask_cell(:)))-1;

cellInd_Other = cellInfo_Other(iexp).responsive_first | cellInfo_Other(iexp).responsive_late | ...
    cellInfo_Other(iexp).responsive_lateWindow | cellInfo_Other(iexp).suppressed_lateWindow;
cellInd_SOM = cellInfo_SOM(iexp).responsive_first | cellInfo_SOM(iexp).responsive_late | ...
    cellInfo_SOM(iexp).responsive_lateWindow | cellInfo_SOM(iexp).suppressed_lateWindow;
[rows,cols] = size(masks.red_mask_cell);

respMask_SOM = zeros(rows*cols,1);
for i = 1:n_red
    if cellInd_SOM(i)
        ind = masks.red_mask_cell(:) == i;
        r = mean(mean(resp_SOM(iexp).av(visualTrials).long(lateWinFr,i,:),3),1);
        respMask_SOM(ind) = r;
    end
end
respMask_Other = zeros(rows*cols,1);
for i = 1:n_green
    if cellInd_Other(i)
        ind = masks.green_mask_cell(:) == i;
        r = mean(mean(resp_Other(iexp).av(visualTrials).long(lateWinFr,i,:),3),1);
        respMask_Other(ind) = r;
    end
end

img_SOM = reshape(respMask_SOM,[rows,cols]);
img_Other = reshape(respMask_Other,[rows,cols]);

x = img_SOM;
x(x ~=0) = 1;
y = bwlabel(x);

[coordX_SOM,coordY_SOM] = roiCentersXY(abs(img_SOM),false);
[coordX_Other,coordY_Other] = roiCentersXY(abs(img_Other),false);

setFigParams4Print('landscape')
figure
colormap(brewermap([],'*RdBu'))
subplot 221
imagesc(img_SOM)
colorbar
clim(respLim)
figAxForm([],0)
title('SOM+ Neurons')
subplot 222
imagesc(img_Other)
colorbar
clim(respLim)
figAxForm([],0)
title('Other Neurons')
subplot 223
imagesc(img_SOM+img_Other)
hold on
plot(coordX_SOM,coordY_SOM,'ko')
colorbar
clim(respLim)
figAxForm([],0)
subplot 224
ind = strcmp({expt.SubNum},exptIDs_mouse_date{1,iexp}) & ...
    strcmp({expt.date},exptIDs_mouse_date{2,iexp});
[~,~,sb_img_50um] = scalebarCalib(exptIDs_mouse_date{2,iexp},expt(ind).obj, img_SOM,1.7);
imagesc(sb_img_50um)
colorbar
clim(respLim)
figAxForm([],0)

print([fnout 'FOVheatmaps_exExpt'],'-dpdf','-fillpage')
%% plot neuron responses as function of distance to SOM neuron

d = getDistanceBetweenNeuronTypes(coordX_Other,coordY_Other,...
    coordX_SOM,coordY_SOM,img_SOM,exptIDs_mouse_date{2,iexp},...
    '16x', 1.7);
figure
subplot 121
ind = cellInfo_SOM(iexp).responsive_lateWindow(cellInd_SOM) == 1;
r = mean(mean(resp_Other(iexp).av(visualTrials).long(lateWinFr,cellInd_Other,:),3),1);
plot(mean(cell2mat(d(ind)'),1),r,'r.','MarkerSize',20)
figXAxis([],'Avg. Distance to Responsive SOM+ Neuron (um)',[])
figYAxis([],'late dF/F',[-0.25 0.25])
figAxForm
title('Responsive SOM+ neurons')

subplot 122
ind = cellInfo_SOM(iexp).suppressed_lateWindow(cellInd_SOM) == 1;
r = mean(mean(resp_Other(iexp).av(visualTrials).long(lateWinFr,cellInd_Other,:),3),1);
plot(mean(cell2mat(d(ind)'),1),r,'b.','MarkerSize',20)
figXAxis([],'Avg. Distance to Suppressed SOM+ Neurons (um)',[])
figYAxis([],'late dF/F',[-0.25 0.25])
figAxForm
title('Suppressed SOM+ neurons')

print([fnout 'respXDist2SOM_exExpt'],'-dpdf','-fillpage')