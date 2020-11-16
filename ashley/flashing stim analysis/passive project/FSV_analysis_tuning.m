
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
                            if isfield(d.tag(2),'oriTuning')
                                expt_somExptOnly(exptN).tag(1).oriTuning = d.tag(2).oriTuning;
                            else
                                expt_somExptOnly(exptN).tag(1).oriTuning = [];
                            end
                            if ~isempty(d_pass)
                                expt_somExptOnly(exptN).tag(1).av_pass(iav) = ...
                                    d_pass.tag(2).av(iav);
                            else
                                expt_somExptOnly(exptN).tag(1).av_pass = [];
                            end
                        else
                            expt_somExptOnly(exptN).tag(itag).av(iav) = d.tag(itag).av(iav);
                            if isfield(d.tag(itag),'oriTuning')
                                expt_somExptOnly(exptN).tag(itag).oriTuning = d.tag(itag).oriTuning;
                            else
                                expt_somExptOnly(exptN).tag(itag).oriTuning = [];
                            end
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

% stimulus response and outcomes
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
            clear r_all outcome stim cycTC
            % target response
            d = expt_somExptOnly(iexp).tag(somInd).av(iav).align(4);
            ind_outcome = strcmp(d.outcome,'success')|strcmp(d.outcome,'ignore');
            r_all = squeeze(mean(d.respTC(respwin,:,ind_outcome),1))...
                -squeeze(mean(d.respTC(basewin_0,:,ind_outcome),1));
            cycTC = squeeze(d.respTC(:,:,ind_outcome))...
                -mean(d.respTC(basewin_0,:,ind_outcome),1);            
            outcome = d.outcome(ind_outcome);
            outcome(strcmp(outcome,'success')) = {'h'};
            outcome(strcmp(outcome,'ignore')) = {'m'};
            if iav == 1
                stim = d.ori(ind_outcome);
            else
                stim = d.amp(ind_outcome);
            end
            
            % correct reject response
            d = expt_somExptOnly(iexp).tag(somInd).av(iav).align(3);
            r = squeeze(mean(d.respTC(respwin,:,:),1))...
                -squeeze(mean(d.respTC(basewin_0,:,:),1));
            r_all = cat(2,r_all,r);            
            cycTC = cat(3,cycTC,squeeze(d.respTC(:,:,:))...
                -mean(d.respTC(basewin_0,:,:),1));
            outcome = cat(2,outcome,repmat({'cr'},[1,length(d.outcome)]));
            stim = cat(2,stim,zeros(size(d.outcome)));
            
            if resp_SOM(iexp).isBx
                % false alarm response
                d = expt_somExptOnly(iexp).tag(somInd).av(iav).align(2);
                r = squeeze(mean(d.respTC(respwin,:,:),1))...
                    -squeeze(mean(d.respTC(basewin_0,:,:),1));
                r_all = cat(2,r_all,r);
                cycTC = cat(3,cycTC,squeeze(d.respTC(:,:,:))...
                    -mean(d.respTC(basewin_0,:,:),1));
                outcome = cat(2,outcome,repmat({'fa'},[1,length(d.outcome)]));
                stim = cat(2,stim,zeros(size(d.outcome)));
            end
            
            resp_SOM(iexp).av(iav).resp = r_all;
            resp_SOM(iexp).av(iav).cycTC = cycTC;
            resp_SOM(iexp).av(iav).outcome = outcome;
            resp_SOM(iexp).av(iav).stim = stim;
            
            if isfield(expt_somExptOnly(iexp).tag(somInd),'av_pass')
                fprintf('.\n')
                if ~isempty(expt_somExptOnly(iexp).tag(somInd).av_pass)
                    clear r_all outcome stim
                    % target response
                    d = expt_somExptOnly(iexp).tag(somInd).av_pass(iav).align(4);
                    r_all = squeeze(mean(d.respTC(respwin,:,:),1))...
                        -squeeze(mean(d.respTC(basewin_0,:,:),1));
                    outcome = d.outcome;
                    outcome(strcmp(outcome,'success')) = {'h'};
                    outcome(strcmp(outcome,'ignore')) = {'m'};
                    if iav == 1
                        stim = d.ori;
                    else
                        stim = d.amp;
                    end

                    % correct reject response
                    d = expt_somExptOnly(iexp).tag(somInd).av_pass(iav).align(3);
                    r = squeeze(mean(d.respTC(respwin,:,:),1))...
                        -squeeze(mean(d.respTC(basewin_0,:,:),1));
                    r_all = cat(2,r_all,r);
                    outcome = cat(2,outcome,repmat({'cr'},[1,length(d.outcome)]));
                    stim = cat(2,stim,zeros(size(d.outcome)));

                    resp_SOM(iexp).av_pass(iav).resp = r_all;
                    resp_SOM(iexp).av_pass(iav).outcome = outcome;
                    resp_SOM(iexp).av_pass(iav).stim = stim;
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
            clear r_all outcome stim cycTC
            % target response
            d = expt_somExptOnly(iexp).tag(otherInd).av(iav).align(4);
            ind_outcome = strcmp(d.outcome,'success')|strcmp(d.outcome,'ignore');
            r_all = squeeze(mean(d.respTC(respwin,:,ind_outcome),1))...
                -squeeze(mean(d.respTC(basewin_0,:,ind_outcome),1));
            cycTC = squeeze(d.respTC(:,:,ind_outcome))...
                -mean(d.respTC(basewin_0,:,ind_outcome),1);     
            outcome = d.outcome(ind_outcome);
            outcome(strcmp(outcome,'success')) = {'h'};
            outcome(strcmp(outcome,'ignore')) = {'m'};
            if iav == 1
                stim = d.ori(ind_outcome);
            else
                stim = d.amp(ind_outcome);
            end
            
            % correct reject response
            d = expt_somExptOnly(iexp).tag(otherInd).av(iav).align(3);
            r = squeeze(mean(d.respTC(respwin,:,:),1))...
                -squeeze(mean(d.respTC(basewin_0,:,:),1));
            r_all = cat(2,r_all,r);
            cycTC = cat(3,cycTC,squeeze(d.respTC(:,:,:))...
                -mean(d.respTC(basewin_0,:,:),1));
            outcome = cat(2,outcome,repmat({'cr'},[1,length(d.outcome)]));
            stim = cat(2,stim,zeros(size(d.outcome)));
            
            if resp_Other(iexp).isBx
                % false alarm response
                d = expt_somExptOnly(iexp).tag(otherInd).av(iav).align(2);
                r = squeeze(mean(d.respTC(respwin,:,:),1))...
                    -squeeze(mean(d.respTC(basewin_0,:,:),1));
                r_all = cat(2,r_all,r);
                cycTC = cat(3,cycTC,squeeze(d.respTC(:,:,:))...
                    -mean(d.respTC(basewin_0,:,:),1));
                outcome = cat(2,outcome,repmat({'fa'},[1,length(d.outcome)]));
                stim = cat(2,stim,zeros(size(d.outcome)));
            end            
            resp_Other(iexp).av(iav).resp = r_all;
            resp_Other(iexp).av(iav).cycTC = cycTC;
            resp_Other(iexp).av(iav).outcome = outcome;
            resp_Other(iexp).av(iav).stim = stim;
            
            if isfield(expt_somExptOnly(iexp).tag(otherInd),'av_pass')
                fprintf('.\n')
                if ~isempty(expt_somExptOnly(iexp).tag(otherInd).av_pass)
                    clear r_all outcome stim
                    % target response
                    d = expt_somExptOnly(iexp).tag(otherInd).av_pass(iav).align(4);
                    r_all = squeeze(mean(d.respTC(respwin,:,:),1))...
                        -squeeze(mean(d.respTC(basewin_0,:,:),1));
                    outcome = d.outcome;
                    outcome(strcmp(outcome,'success')) = {'h'};
                    outcome(strcmp(outcome,'ignore')) = {'m'};
                    if iav == 1
                        stim = d.ori;
                    else
                        stim = d.amp;
                    end

                    % correct reject response
                    d = expt_somExptOnly(iexp).tag(otherInd).av_pass(iav).align(3);
                    r = squeeze(mean(d.respTC(respwin,:,:),1))...
                        -squeeze(mean(d.respTC(basewin_0,:,:),1));
                    r_all = cat(2,r_all,r);
                    outcome = cat(2,outcome,repmat({'cr'},[1,length(d.outcome)]));
                    stim = cat(2,stim,zeros(size(d.outcome)));

                    resp_Other(iexp).av_pass(iav).resp = r_all;
                    resp_Other(iexp).av_pass(iav).outcome = outcome;
                    resp_Other(iexp).av_pass(iav).stim = stim;
                end
            end
        end
    end
end

% oriTuning
ori_SOM = struct;
ori_Other = struct;
for iexp = 1:nexp    
    if expt_somExptOnly(iexp).tag(1).name == 'SOM'
        somInd = 1;
        otherInd = 2;
    else
        somInd = 2;
        otherInd = 1;
    end
    d = expt_somExptOnly(iexp).tag(somInd);
    ori_SOM(iexp).oriTuning = d.oriTuning;
    d = expt_somExptOnly(iexp).tag(otherInd);
    ori_Other(iexp).oriTuning = d.oriTuning;    
end

cellInfo_SOM = struct;
cellInfo_Other = struct;
for iexp = 1:nexp
    d = resp_SOM(iexp).av(visualTrials);
    ind_late = strcmp(d.outcome,'cr') | strcmp(d.outcome,'fa');
    ind_tar = strcmp(d.outcome,'h') | strcmp(d.outcome,'m');
    
    cycTC_late = d.cycTC(:,:,ind_late);
    cycTC_tar = d.cycTC(:,:,ind_tar);
        
    cellInfo_SOM(iexp).responsive_late = ttest(squeeze(mean(cycTC_late(respwin,:,:))),...
        squeeze(mean(cycTC_late(basewin_0,:,:))),'dim',2,'tail','right');
    tar_all = ttest(squeeze(mean(cycTC_tar(respwin,:,:))),...
        squeeze(mean(cycTC_tar(basewin_0,:,:))),'dim',2,'tail','right');
    targets = unique(d.stim(d.stim~=0));
    tar_ori = cell(1,length(targets));
    for itar = 1:length(targets)
        ind = d.stim == targets(itar);
        tar_ori{itar} = ttest(squeeze(mean(cycTC_tar(respwin,:,ind))),...
            squeeze(mean(cycTC_tar(basewin_0,:,ind))),'dim',2,'tail','right');
    end
    cellInfo_SOM(iexp).responsive_target = tar_all | sum(cell2mat(tar_ori),2)>0;
    
    if ~isempty(resp_Other(iexp).av) 
        d = resp_Other(iexp).av(visualTrials);
        ind_late = strcmp(d.outcome,'cr') | strcmp(d.outcome,'fa');
        ind_tar = strcmp(d.outcome,'h') | strcmp(d.outcome,'m');

        cycTC_late = d.cycTC(:,:,ind_late);
        cycTC_tar = d.cycTC(:,:,ind_tar);

        cellInfo_Other(iexp).responsive_late = ttest(squeeze(mean(cycTC_late(respwin,:,:))),...
            squeeze(mean(cycTC_late(basewin_0,:,:))),'dim',2,'tail','right');
        tar_all = ttest(squeeze(mean(cycTC_tar(respwin,:,:))),...
            squeeze(mean(cycTC_tar(basewin_0,:,:))),'dim',2,'tail','right');
        targets = unique(d.stim(d.stim~=0));
        tar_ori = cell(1,length(targets));
        for itar = 1:length(targets)
            ind = d.stim == targets(itar);
            tar_ori{itar} = ttest(squeeze(mean(cycTC_tar(respwin,:,ind))),...
                squeeze(mean(cycTC_tar(basewin_0,:,ind))),'dim',2,'tail','right');
        end
        cellInfo_Other(iexp).responsive_target = tar_all | sum(cell2mat(tar_ori),2)>0;
    else
        cellInfo_Other(iexp).responsive_late = [];
        cellInfo_Other(iexp).responsive_target = [];
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

%%
exExpt = nan;
exNeuron_1 = nan;
exNeuron_2 = nan;
%%
tuning_SOM = struct;
tuning_Other = struct;
for iexp = 1:nexp
    if iexp == 1
        tuning_SOM.isResp = [];
        tuning_SOM.reliability = [];
        tuning_SOM.tuningFit = [];
        tuning_SOM.tuningResp = [];        
        tuning_SOM.FWHM = [];
        tuning_SOM.taskTuning = [];
        tuning_SOM.taskTuning_pass = [];
        tuning_SOM.exptN = [];
        
        tuning_Other.isResp = [];
        tuning_Other.reliability = [];
        tuning_Other.tuningFit = [];
        tuning_Other.tuningResp = [];  
        tuning_Other.FWHM = [];
        tuning_Other.taskTuning = [];
        tuning_Other.taskTuning_pass = [];
        tuning_Other.exptN = [];
    end
    if ~isempty(ori_SOM(iexp).oriTuning)
        tuning_SOM.isResp = cat(1,tuning_SOM.isResp,...
            cellInfo_SOM(iexp).responsive_late|cellInfo_SOM(iexp).responsive_target);
                
        tuning_SOM.reliability = cat(2,tuning_SOM.reliability,...
            ori_SOM(iexp).oriTuning.oriFitReliability);
        tf = ori_SOM(iexp).oriTuning.oriFit;
        fwhm = fwhmFromOriFit(tf,0:180);
        tuning_SOM.tuningFit = cat(2,tuning_SOM.tuningFit,tf);
        tuning_SOM.FWHM = cat(2,tuning_SOM.FWHM,fwhm);
        tuning_SOM.tuningResp = cat(1,tuning_SOM.tuningResp,...
            ori_SOM(iexp).oriTuning.oriResp); 
        
        d = resp_SOM(iexp).av(visualTrials);
        stim = discretize(d.stim,oriBins);
        targets = unique(stim);
        tasktuning = nan(length(targets),size(d.resp,1));
        for itar = 1:length(targets)
            ind = stim == targets(itar);
            tasktuning(itar,:) = mean(d.resp(:,ind),2);
        end
        tasktuning(tasktuning<0) = 0;
        tuning_SOM.taskTuning = cat(2,tuning_SOM.taskTuning,tasktuning);
        
        if isempty(resp_SOM(iexp).av_pass)
            tuning_SOM.taskTuning_pass = cat(2,tuning_SOM.taskTuning_pass,...
                nan(size(tasktuning)));
        else
            d = resp_SOM(iexp).av_pass(visualTrials);
            stim = discretize(d.stim,oriBins);
            targets = unique(stim);
            tasktuning = nan(length(targets),size(d.resp,1));
            for itar = 1:length(targets)
                ind = stim == targets(itar);
                tasktuning(itar,:) = mean(d.resp(:,ind),2);
            end
            tasktuning(tasktuning<0) = 0;
            tuning_SOM.taskTuning_pass = cat(2,tuning_SOM.taskTuning_pass,...
                tasktuning);
            
        end
        
        tuning_SOM.exptN = cat(2,tuning_SOM.exptN,...
            ones(1,length(cellInfo_SOM(iexp).responsive_late)).*iexp);
        
        if ~isempty(resp_Other(iexp).av)
            tuning_Other.isResp = cat(1,tuning_Other.isResp,...
                cellInfo_Other(iexp).responsive_late|cellInfo_Other(iexp).responsive_target);

            tuning_Other.reliability = cat(2,tuning_Other.reliability,...
                ori_Other(iexp).oriTuning.oriFitReliability);
            tf = ori_Other(iexp).oriTuning.oriFit;
            fwhm = fwhmFromOriFit(tf,0:180);
            tuning_Other.tuningFit = cat(2,tuning_Other.tuningFit,tf);
            tuning_Other.FWHM = cat(2,tuning_Other.FWHM,fwhm);
            tuning_Other.tuningResp = cat(1,tuning_Other.tuningResp,...
                ori_Other(iexp).oriTuning.oriResp); 

            d = resp_Other(iexp).av(visualTrials);
            stim = discretize(d.stim,oriBins);
            targets = unique(stim);
            tasktuning = nan(length(targets),size(d.resp,1));
            for itar = 1:length(targets)
                ind = stim == targets(itar);
                tasktuning(itar,:) = mean(d.resp(:,ind),2);
            end
            tasktuning(tasktuning<0) = 0;
            tuning_Other.taskTuning = cat(2,tuning_Other.taskTuning,tasktuning);

            if isempty(resp_Other(iexp).av_pass)
                tuning_Other.taskTuning_pass = cat(2,tuning_Other.taskTuning_pass,...
                    nan(size(tasktuning)));
            else
                d = resp_Other(iexp).av_pass(visualTrials);
                stim = discretize(d.stim,oriBins);
                targets = unique(stim);
                tasktuning = nan(length(targets),size(d.resp,1));
                for itar = 1:length(targets)
                    ind = stim == targets(itar);
                    tasktuning(itar,:) = mean(d.resp(:,ind),2);
                end
                tasktuning(tasktuning<0) = 0;
                tuning_Other.taskTuning_pass = cat(2,tuning_Other.taskTuning_pass,...
                    tasktuning);

            end

            tuning_Other.exptN = cat(2,tuning_Other.exptN,...
                ones(1,length(cellInfo_Other(iexp).responsive_late)).*iexp);
        end
    end
end

%%
exNeuron_1 = 1;
exNeuron_2 = 1;

figure
subplot 121
y = tuning_SOM.tuningResp(exNeuron_1,:);
yfit = tuning_SOM.tuningFit(:,exNeuron_1);
plot(0:22.5:179,y,'k.','MarkerSize',20)
hold on
plot(0:180,yfit,'k:','LineWidth',2)
figXAxis([],'Orientation (deg)',[0 180])
figYAxis([],'dF/F',[])
figAxForm
title('SOM+ Neuron')
subplot 122
y = tuning_Other.tuningResp(exNeuron_2,:);
yfit = tuning_Other.tuningFit(:,exNeuron_2);
plot(0:22.5:179,y,'k.','MarkerSize',20)
hold on
plot(0:180,yfit,'k:','LineWidth',2)
figXAxis([],'Orientation (deg)',[0 180])
figYAxis([],'dF/F',[])
figAxForm
title('SOM- Neuron')

print([fnout 'tuning_exampleNeurons'],'-dpdf')
%%
binEdges = 0:10:90;
figure
subplot 321
histogram(tuning_SOM.reliability(tuning_SOM.isResp==1),binEdges)
hold on
histogram(tuning_Other.reliability(tuning_Other.isResp==1),binEdges)
vline(tuningReliabilityThresh,'k--')
figXAxis([],'Reliability of Fit',[0 90])
figYAxis([],'N Cells',[0 80])
figAxForm
title('Task Responsive Cells')
legend({'SOM+','SOM-'})

subplot 322
hold on
[y,ci] = binofit(sum(tuning_SOM.reliability(tuning_SOM.isResp==1)<tuningReliabilityThresh),...
    length(tuning_SOM.reliability(tuning_SOM.isResp==1)));
errorbar(1,y,y-ci(1),ci(2)-y,'k.','MarkerSize',20)
[y,ci] = binofit(sum(tuning_Other.reliability(tuning_Other.isResp==1)<tuningReliabilityThresh),...
    length(tuning_Other.reliability(tuning_Other.isResp==1)));
errorbar(2,y,y-ci(1),ci(2)-y,'k.','MarkerSize',20)
figXAxis([],'',[0 3],1:2,{'SOM+','SOM-'})
figYAxis([],'FWHM',[0 1])
figAxForm
title('Task Responsive Cells')

subplot 323
ind_SOM = tuning_SOM.reliability<tuningReliabilityThresh;
ind_Other = tuning_Other.reliability<tuningReliabilityThresh;
histogram(tuning_SOM.FWHM(tuning_SOM.isResp' & ind_SOM),binEdges)
hold on
histogram(tuning_Other.FWHM(tuning_Other.isResp' & ind_Other),binEdges)
vline(tuningReliabilityThresh,'k--')
figXAxis([],'FWHM',[0 90])
figYAxis([],'N Cells',[0 90])
figAxForm
title('Task Responsive Cells')
legend({'SOM+','SOM-'})

subplot 324
hold on
y1 = mean(tuning_SOM.FWHM(tuning_SOM.isResp' & ind_SOM));
yerr1 = ste(tuning_SOM.FWHM(tuning_SOM.isResp' & ind_SOM),2);
y2 = mean(tuning_Other.FWHM(tuning_Other.isResp' & ind_Other));
yerr2 = ste(tuning_Other.FWHM(tuning_Other.isResp' & ind_Other),2);
errorbar(1:2,[y1,y2],[yerr1,yerr2],'k.','MarkerSize',20)

y1 = mean(tuning_SOM.FWHM(tuning_SOM.isResp' & ~isnan(tuning_SOM.FWHM)));
yerr1 = ste(tuning_SOM.FWHM(tuning_SOM.isResp' & ~isnan(tuning_SOM.FWHM)),2);
y2 = mean(tuning_Other.FWHM(tuning_Other.isResp' & ~isnan(tuning_Other.FWHM)));
yerr2 = ste(tuning_Other.FWHM(tuning_Other.isResp' & ~isnan(tuning_Other.FWHM)),2);
errorbar(1:2,[y1,y2],[yerr1,yerr2],'b.','MarkerSize',20)
figXAxis([],'',[0 3],1:2,{'SOM+','SOM-'})
figYAxis([],'FWHM',[])
figAxForm
legend({'Tuned','All Resp.'})

expt_ind = find(bxInd&passInd);

subplot 325
hold on
ind = tuning_SOM.isResp'==1 & ismember(tuning_SOM.exptN,expt_ind);
y = mean(tuning_SOM.taskTuning(:,ind),2);
yerr = ste(tuning_SOM.taskTuning(:,ind),2);
errorbar(1:3,y,yerr,'k.-','MarkerSize',20)
y = mean(tuning_SOM.taskTuning_pass(:,ind),2);
yerr = ste(tuning_SOM.taskTuning_pass(:,ind),2);
errorbar(1:3,y,yerr,'b.-','MarkerSize',20)
title(sprintf('SOM+ Cells, n=%s',num2str(sum(ind))))
legend({'behav','pass'})
figXAxis([],'Orientation (deg)',[0.5 3.5],1:3,{'0','10-32','33-90'})
figYAxis([],'dF/F',[-0.01 0.08])
figAxForm

subplot 326
hold on
ind = tuning_Other.isResp'==1 & ismember(tuning_Other.exptN,expt_ind);
y = mean(tuning_Other.taskTuning(:,ind),2);
yerr = ste(tuning_Other.taskTuning(:,ind),2);
errorbar(1:3,y,yerr,'k.-','MarkerSize',20)
y = mean(tuning_Other.taskTuning_pass(:,ind),2);
yerr = ste(tuning_Other.taskTuning_pass(:,ind),2);
errorbar(1:3,y,yerr,'b.-','MarkerSize',20)
title(sprintf('SOM- Cells, n=%s',num2str(sum(ind))))
figXAxis([],'Orientation (deg)',[0.5 3.5],1:3,{'0','10-32','33-90'})
figYAxis([],'dF/F',[-0.01 0.08])
figAxForm


print([fnout 'tuning_SOMvsOther'],'-dpdf','-fillpage')

