
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
            clear r_all outcome stim
            % target response
            d = expt_somExptOnly(iexp).tag(somInd).av(iav).align(4);
            ind_outcome = strcmp(d.outcome,'success')|strcmp(d.outcome,'ignore');
            r_all = squeeze(mean(d.respTC(respwin,:,ind_outcome),1))...
                -squeeze(mean(d.respTC(basewin_0,:,ind_outcome),1));
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
            outcome = cat(2,outcome,repmat({'cr'},[1,length(d.outcome)]));
            stim = cat(2,stim,zeros(size(d.outcome)));
            
            if resp_SOM(iexp).isBx
                % false alarm response
                d = expt_somExptOnly(iexp).tag(somInd).av(iav).align(2);
                r = squeeze(mean(d.respTC(respwin,:,:),1))...
                    -squeeze(mean(d.respTC(basewin_0,:,:),1));
                r_all = cat(2,r_all,r);
                outcome = cat(2,outcome,repmat({'fa'},[1,length(d.outcome)]));
                stim = cat(2,stim,zeros(size(d.outcome)));
            end            
            resp_SOM(iexp).av(iav).resp = r_all;
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
            clear r_all outcome stim
            % target response
            d = expt_somExptOnly(iexp).tag(otherInd).av(iav).align(4);
            ind_outcome = strcmp(d.outcome,'success')|strcmp(d.outcome,'ignore');
            r_all = squeeze(mean(d.respTC(respwin,:,ind_outcome),1))...
                -squeeze(mean(d.respTC(basewin_0,:,ind_outcome),1));
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
            outcome = cat(2,outcome,repmat({'cr'},[1,length(d.outcome)]));
            stim = cat(2,stim,zeros(size(d.outcome)));
            
            if resp_Other(iexp).isBx
                % false alarm response
                d = expt_somExptOnly(iexp).tag(otherInd).av(iav).align(2);
                r = squeeze(mean(d.respTC(respwin,:,:),1))...
                    -squeeze(mean(d.respTC(basewin_0,:,:),1));
                r_all = cat(2,r_all,r);
                outcome = cat(2,outcome,repmat({'fa'},[1,length(d.outcome)]));
                stim = cat(2,stim,zeros(size(d.outcome)));
            end            
            resp_Other(iexp).av(iav).resp = r_all;
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

%% run model for behavior and passive experiments, include SOM neurons
nPCs = 8;
dc = struct;
exptN = 0;
for iexp = 1:nexp
        if bxInd(iexp) && passInd(iexp) && iexp ~=21
            exptN = exptN+1;
            for iav = 1:2
                % organize response data and stims
                r_all = cat(2,cat(1,resp_SOM(iexp).av(iav).resp,resp_Other(iexp).av(iav).resp),...
                    cat(1,resp_SOM(iexp).av_pass(iav).resp,resp_Other(iexp).av_pass(iav).resp))';
                npass = size(resp_Other(iexp).av_pass(iav).resp,2);

                somInd = 1:size(resp_SOM(iexp).av(iav).resp,1);


                % get PCs and coefficients and zscore
                [coeffAllCells,scoresAllCells,latentAllCells] = pca(r);
                r_zs = zscore(scoresAllCells);

                % match trial difficulty
                trStim = resp_SOM(iexp).av(iav).stim;

                if iav == 1
                    trStimID = discretize(trStim,oriBins);
                elseif iav == 2
                    trStimID = discretize(trStim,ampBins);
                end
                nStimPerBin = histcounts(trStimID);
                minBinN = min(nStimPerBin(nStimPerBin >= minTrN_mdl));
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
                    end
                end

                % run model
                trOut = resp_SOM(iexp).av(iav).outcome;
                r_bx = r_all(1:(end-npass),:);
                r_mod = r_bx(matchTrialsInd,1:nPCs);
                [choiceTrInd, stimTrInd] = getStimAndBehaviorYs(trOut(matchTrialsInd));

                C = eye(size(r_mod,2));
                p=1;
                [~,~,choiceGLM] = glmfit(r_mod*C,choiceTrInd,'binomial');
                [~,~,stimulusGLM] = glmfit(r_mod*C,stimTrInd,'binomial');

                dv_choice = mean(choiceTrInd);
                dv_stim = mean(stimTrInd);

                % test model
                fprintf('Expt %s, hold-out analysis\n',num2str(iexp))
                pctCorrectChoice_train = getPctCorr_trainData(choiceGLM,r_mod,choiceTrInd,dv_choice);
                pctCorrectChoice_ho_bx = getPctCorr_hoData(r_mod,choiceTrInd,dv_choice);

                pctCorrectStimulus_train = getPctCorr_trainData(stimulusGLM,r_mod,stimTrInd,dv_stim);
                pctCorrectStimulus_ho_bx = getPctCorr_hoData(r_mod,stimTrInd,dv_stim);

                % test passive trials 

                % match trial difficulty
                trStim = resp_SOM(iexp).av_pass(iav).stim;

                if iav == 1
                    trStimID = discretize(trStim,oriBins);
                elseif iav == 2
                    trStimID = discretize(trStim,ampBins);
                end
                nStimPerBin = histcounts(trStimID);
                minBinN = min(nStimPerBin(nStimPerBin >= minTrN_mdl));
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
                    end
                end

                trOut = resp_SOM(iexp).av_pass(iav).outcome;
                r_pass = r_all((size(r_all,1)-npass+1):end,:);
                r_mod = r_pass(matchTrialsInd,1:nPCs);
                [~, stimTrInd] = getStimAndBehaviorYs(trOut(matchTrialsInd));

                dv_stim = mean(stimTrInd);
                pctCorrectStimulus_trainBxTestPass = getPctCorr_trainData(stimulusGLM,r_mod,stimTrInd,dv_stim);
                pctCorrectChoice_trainBxTestPass = getPctCorr_trainData(choiceGLM,r_mod,stimTrInd,dv_choice);


                % train passive trials

                pctCorrectStimulus_ho_pass = getPctCorr_hoData(r_mod,stimTrInd,dv_stim);

                % save analysis in struct
                dc(exptN).av(iav).stimulus_pctCorr_trainStimBxTestStimBx = pctCorrectStimulus_ho_bx;
                dc(exptN).av(iav).stimulus_pctCorr_trainStimBxTestStimPass = pctCorrectStimulus_trainBxTestPass;
                dc(exptN).av(iav).stimulus_pctCorr_trainStimPassTestStimPass = pctCorrectStimulus_ho_pass;

                dc(exptN).av(iav).choice_pctCorr_trainChBxTestChBx = pctCorrectChoice_ho_bx;
                dc(exptN).av(iav).choice_pctCorr_trainChBxTestStimPass = pctCorrectChoice_trainBxTestPass;
            end
        end
end

%% plot performance of behavior and passive models
figure
suptitle(sprintf(...
    'All Trials Performance (held-out data), pct corr must be > %s in vis stim model',...
    num2str(pctCorrThresh)))
n = size(dc,2);
for iav = 1:2
    subplot(2,2,iav)
    hold on
    y1 = [];
    for i = 1:n
        y1 = cat(2,y1, dc(i).av(iav).stimulus_pctCorr_trainStimBxTestStimBx);
    end
    y2 = [];
    for i = 1:n
        y2 = cat(2,y2, dc(i).av(iav).stimulus_pctCorr_trainStimPassTestStimPass);
    end
    plot(1:2,[y1',y2'],'k.-','MarkerSize',10)
    x = 1:2;
    y_all = [mean(y1),mean(y2)];
    y_all_ste = [ste(y1,2),ste(y2,2)];
    figXAxis([],'',[0 3],1:2,{'Bx','Pass'})
    errorbar(x,y_all,y_all_ste,'.','MarkerSize',20)
    figYAxis([],'% Correct',[0 1],0:0.2:1)
    figAxForm([],0)
    hline(0.5,'k:')
    [~,p] = ttest2(y1,y2);
    title(sprintf('%s Stim Model, p=%s',avName{iav},...
        num2str(round(p,2,'significant'))))
    
    
    subplot(2,2,iav+2)
    hold on
    y1 = [];
    for i = 1:n
        y1 = cat(2,y1, dc(i).av(iav).stimulus_pctCorr_trainStimBxTestStimBx);
    end
    y2 = [];
    for i = 1:n
        y2 = cat(2,y2, dc(i).av(iav).stimulus_pctCorr_trainStimBxTestStimPass);
    end
    plot(1:2,[y1',y2'],'k.-','MarkerSize',10)
    x = 1:2;
    y_all = [mean(y1),mean(y2)];
    y_all_ste = [ste(y1,2),ste(y2,2)];
    figXAxis([],'Test',[0 3],1:2,{'Bx','Pass'})
    errorbar(x,y_all,y_all_ste,'.','MarkerSize',20)
    figYAxis([],'% Correct',[0 1],0:0.2:1)
    figAxForm([],0)
    hline(0.5,'k:')
    [~,p] = ttest2(y1,y2);
    title(sprintf('Train %s Behavior Stim Model, p=%s',avName{iav},...
        num2str(round(p,2,'significant'))))
    
    
end

print([fnout 'decoding_bxVsPass'],'-dpdf','-fillpage')


%% plot performance of behavior model with passive data