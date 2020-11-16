function mouse = createFSAVDataStruct_2color(datasetStr,cellsOrDendrites)
%cellsOrDendrites: 1 == cells; 2 == dendrites
    % set analysis windows
    pre_event_time = 1000; %ms
    post_event_time = 4500; %ms
    resp_win_time = 100; %ms
    pre_win_time = [-30 70];
    trans_win_time = [150 250]; %this is actually ~166-266 ms at 30 Hz
    minTrialLengthMs = 2500; % time in ms of trial rather than hard-coding cycle number for analysis
    
    if cellsOrDendrites == 1
        motionThreshold = 0.1;
    elseif cellsOrDendrites == 2
        motionThreshold = 0.15;
    end
    
    rc = behavConstsAV;

    eval(datasetStr)

    dataGroup = datasetStr;
    
    %create list of mice for file names
    nMice = length(unique({expt.SubNum}));
    str = unique({expt.SubNum});
    values = cell2mat(cellfun(@str2num,str,'UniformOutput',false));
    mouse_str = ['i' strjoin(str,'_i')];
    %initialize structure
    exptN = zeros(1,nMice);
    mouse = struct;
    mousePass = struct;
    
    
    %collect data from each experiment
    for iexp = 1:size(expt,2)
        disp([num2str(expt(iexp).date) ' i' num2str(expt(iexp).SubNum)])

        nrun = size(expt(iexp).runs,1);
        
        %keep track of which mouse
        imouse = find(values == str2num(expt(iexp).SubNum));
        mouse(imouse).mouse_name = expt(iexp).SubNum;
        exptN(:,imouse) = exptN(:,imouse)+1;
        mouse(imouse).expt(exptN(:,imouse)).date = expt(iexp).date;
        mouse(imouse).expt(exptN(:,imouse)).attnTask = expt(iexp).attentionTask;
        mouse(imouse).expt(exptN(:,imouse)).hasAttn = expt(iexp).hasAttention;
        
        
        %translate time windows into frames 
        pre_event_frames = ceil(pre_event_time*(expt(iexp).frame_rate/1000));
        post_event_frames = ceil(post_event_time*(expt(iexp).frame_rate/1000));
        basewin = (pre_event_frames-1):pre_event_frames+1;
        
        mouse(imouse).expt(exptN(:,imouse)).info.preAlignFrames = pre_event_frames;
                
        %create string for saving mult runs
        runstr = expt(iexp).runs(1,:);
        if nrun>1
            for irun = 2:nrun
                runstr = [runstr '-' expt(iexp).runs(irun,:)];
            end
        end
%         fnout = fullfile(rc.caOutputDir, expt(iexp).mouse, expt(iexp).folder, expt(iexp).date, [expt(iexp).date '_' expt(iexp).mouse '_' runstr '_']);
        

      
        %account for accumulation of frames across multiple runs 
        dataTC_notag = [];
        dataTC_tag = [];
        fnTC = fullfile(rc.ashleyAnalysis,...
            expt(iexp).mouse,expt(iexp).folder, expt(iexp).date,'data processing');
        
        if cellsOrDendrites == 1
            load(fullfile(fnTC,'timecourses_bx_cells.mat'))
            if isempty(expt(iexp).redChannelLabel)
                dataTC_tag = data_bx_tc_subnp; 
                dataTC_notag = nan(size(dataTC_tag,1),0);               
            else
                dataTC_notag = data_bx_g_tc_subnp;
                dataTC_tag = data_bx_r_tc_subnp;
            end
        end

        % load and combine mworks data and timecourses
        input = [];
        for irun = 1:nrun
            time = expt(iexp).time_mat(irun,:);
            fn_mworks = [rc.pathStr...
                '\data-i' expt(iexp).SubNum '-' expt(iexp).date '-' time '.mat'];
            if irun == 1
                input = mwLoadData(fn_mworks, [], []);
            else
                try
                    input = [input mwLoadData(fn_mworks, [], [])];
                catch
                    input2 = mwLoadData(fn_mworks, [], []);
                    inpNames1 = fieldnames(input);
                    inpNames2 = fieldnames(input2);
                    inpLong = gt(length(inpNames1),length(inpNames2));
                    if inpLong == 1
                        inpPlusInd = ismember(inpNames1,inpNames2);
                        inpPlus = inpNames1(~inpPlusInd);
                        for i = 1:length(inpPlus)
                            input2.(genvarname(inpPlus{i})) = ...
                                cell(1,input2.trialSinceReset);
                        end
                    else
                        inpPlusInd = ismember(inpNames2,inpNames1);
                        inpPlus = inpNames2(~inpPlusInd);
                        for i = 1:length(inpPlus)
                            input.(char(genvarname(inpPlus(i)))) = cell(1,80);
                        end
                    end
                    input = [input input2];
                end
            end
        end
        input = concatenateDataBlocks(input);  
        if isnan(expt(iexp).trial_range)
            tr = 1:length(input.trialOutcomeCell);
        else
            tr = expt(iexp).trial_range;
        end
        
        %convert important fields to matrices
        run_trials = input.trialsSinceReset;
        cLeverDown = celleqel2mat_padded(input.cLeverDown);
        cFirstStim = celleqel2mat_padded(input.cFirstStim);
        cLeverUp = celleqel2mat_padded(input.cLeverUp);
        cTargetOn = celleqel2mat_padded(input.cTargetOn);
        cStimOn = celleqel2mat_padded(input.cStimOn);
        cItiStart = celleqel2mat_padded(input.cItiStart); 
        
        
        %stim timing
        if iscell(input.nFramesOn)
            cycTime = unique(cell2mat(input.nFramesOn))+unique(cell2mat(input.nFramesOff));
        else
            cycTime = input.nFramesOn+input.nFramesOff;
        end
        
        mouse(imouse).expt(exptN(:,imouse)).info.cycTimeFrames = cycTime;
        
%         reactTimes = celleqel2mat_padded(input.reactTimesMs);
        holdTimesMs = celleqel2mat_padded(input.holdTimesMs);
        tCyclesOn = double(cell2mat(input.tCyclesOn));
        nCyclesOn = double(cell2mat(input.nCyclesOn));
        holdTimesMs = holdTimesMs(tr);
        tCyclesOn = tCyclesOn(tr);
        nCyclesOn = nCyclesOn(tr);        
        
        frameRate = input.frameRateHz;
        cycTimeMs = cycTime./frameRate*1000;
        requiredStimTimeMs = nCyclesOn.*cycTimeMs;
        reactTimeCalc = holdTimesMs - requiredStimTimeMs;
        reactTimeFromLastBaseMs = holdTimesMs - ((tCyclesOn-1).*cycTimeMs);
        ntrials = length(tr);
        
        
        %catch trial Info
        if expt(iexp).catch
            cCatchOn = celleqel2mat_padded(input.cCatchOn);
            nt = length(cCatchOn);
            tCatchDirection = celleqel2mat_padded(input.tCatchGratingDirectionDeg);
            
%             trialStartFr = celleqel2mat_padded(input.cFirstStim);
            trialEndFr = celleqel2mat_padded(input.cLeverUp);
            cCyc = celleqel2mat_padded(input.catchCyclesOn);
            cCyc = cCyc(tr);
            catchReactFr = trialEndFr - cCatchOn;
            catchCycOnMs = (cCyc-1).*cycTimeMs;
            catchReactMs = holdTimesMs - catchCycOnMs;
            tooFastFr = round(frameRate*((input.tooFastTimeMs)/1000));
            reactTimeFr = round(frameRate*((input.reactTimeMs)/1000));
        
            catchOutcome = cell(1,nt);
            outInd = catchReactFr > tooFastFr & catchReactFr < reactTimeFr;
            catchOutcome(outInd) = {'FA'};
            outInd = catchReactFr > reactTimeFr;
            catchOutcome(outInd) = {'CR'};
            outInd = catchReactFr < 0;
            catchOutcome(outInd) = {'failure'};
            
            tCatchDirection = tCatchDirection(tr);
%             catchReactMs = catchReactMs(tr);
%             catchOutcome = catchOutcome(tr);
        end  
                
        offset = 0;
        for irun = 1:nrun
            ImgFolder = expt(iexp).runs(irun,:);
            if isempty(expt(iexp).nframesPerRun)
                nfr_run = nFramesSbxDataset(expt(iexp).mouse,expt(iexp).date,ImgFolder);
            else
                nfr_run = expt(iexp).nframesPerRun(irun);
            end
            offset = offset+nfr_run;
            if irun < nrun
                startTrial = sum(run_trials(1, 1:irun),2)+1;
                endTrial = sum(run_trials(1,1:irun+1),2);
                cLeverDown(1,startTrial:endTrial) = cLeverDown(1,startTrial:endTrial)+offset;
                cFirstStim(1,startTrial:endTrial) = cFirstStim(1,startTrial:endTrial)+offset;
                cLeverUp(1,startTrial:endTrial) = cLeverUp(1,startTrial:endTrial)+offset;
                cTargetOn(1,startTrial:endTrial) = cTargetOn(1,startTrial:endTrial)+offset;
                cStimOn(1,startTrial:endTrial) = cStimOn(1,startTrial:endTrial)+offset;
                cItiStart(1,startTrial:endTrial) = cItiStart(1,startTrial:endTrial)+offset;
                if expt(iexp).catch
                    cCatchOn(1,startTrial:endTrial) = cCatchOn(1,startTrial:endTrial)+offset;
                end
            end
        end
        
        cLeverDown = cLeverDown(tr);
        cFirstStim = cFirstStim(tr);
        cLeverUp = cLeverUp(tr);
        cTargetOn = cTargetOn(tr);
        cStimOn = cStimOn(tr);
        cItiStart = cItiStart(tr);
        if expt(iexp).catch
            cCatchOn = cCatchOn(tr);
        end
        
        

        %previous trial info
        trType = double(cell2mat(input.tBlock2TrialNumber));
        trType = trType(tr);
        trType_shift = [NaN trType];
        prevTrType = num2cell(trType_shift(1:length(trType)));
        trType_shift = [NaN NaN trType];
        prev2TrType = num2cell(trType_shift(1:length(trType)));
        for i = 1:length(trType)
            if prevTrType{i} == 0
                prevTrType{i} = 'vis';
            else
                prevTrType{i} = 'aud';
            end
            if prev2TrType{i} == 0
                prev2TrType{i} = 'vis';
            else
                prev2TrType{i} = 'aud';
            end
        end
            
        trOut = input.trialOutcomeCell;
        trOut = trOut(tr);
        trOut_shift = [{NaN} trOut];
        prevTrOut = trOut_shift(1:length(trOut));
        trOut_shift = [{NaN} {NaN} trOut];
        prev2TrOut = trOut_shift(1:length(trOut));
        
        
        if expt(iexp).catch
            isCatchTrial = logical(cell2mat(input.tShortCatchTrial));
        else
            isCatchTrial = false(1,ntrials);
        end
        isCatchTrial = isCatchTrial(tr);

        tGratingDirectionDeg = chop(celleqel2mat_padded(input.tGratingDirectionDeg),4);
        tGratingDirectionDeg = tGratingDirectionDeg(tr);
        Dirs = unique(tGratingDirectionDeg);
        
        if isfield(input,'tSoundTargetAmplitude')
            tSoundTargetAmp = celleqel2mat_padded(input.tSoundTargetAmplitude);
        else
            tSoundTargetAmp = ones(1,ntrials).*input.soundTargetAmplitude;
        end
        tSoundTargetAmp = tSoundTargetAmp(tr);
        Amps = unique(tSoundTargetAmp);
        
        mouse(imouse).expt(exptN(:,imouse)).info.visTargets = Dirs;
        mouse(imouse).expt(exptN(:,imouse)).info.audTargets = Amps;
        mouse(imouse).expt(exptN(:,imouse)).info.fsavSize = input.gratingHeightDeg;
        
%         
%         mouse(imouse).expt(exptN(:,imouse)).oriTuning.oriResp = avgResponseEaOri;
%         mouse(imouse).expt(exptN(:,imouse)).oriTuning.oriRespSem = semResponseEaOri;
%         mouse(imouse).expt(exptN(:,imouse)).oriTuning.oriFit = vonMisesFitAllCells;
%         mouse(imouse).expt(exptN(:,imouse)).oriTuning.oriFitReliability = ...
%             fitReliability;

        
        mouse(imouse).expt(exptN(:,imouse)).tag(1).name = [expt(iexp).redChannelLabel '-'];
        mouse(imouse).expt(exptN(:,imouse)).tag(2).name = [expt(iexp).redChannelLabel '+'];
        for itag = 1:2
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(1).name = 'visual';
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(2).name = 'auditory';
        end
        
        
        %load direction tuning data
        if ~isnan(expt(iexp).dirtuning)
            fnTun = fullfile(rc.ashleyAnalysis,...
                expt(iexp).mouse,expt(iexp).folder, expt(iexp).date, ...
                expt(iexp).dirtuning);
            if all(isnan(dataTC_notag(:)))
                if cellsOrDendrites == 1
                    load(fullfile(fnTun, 'oriTuningAndFits.mat'));
                elseif cellsOrDendrites == 2
                    load(fullfile(fnTun, 'oriTuningAndFits_den.mat'));
                end
                mouse(imouse).expt(exptN(:,imouse)).tag(1).oriTuning.oriResp = avgResponseEaOri;
                mouse(imouse).expt(exptN(:,imouse)).tag(1).oriTuning.oriRespSem = semResponseEaOri;
                mouse(imouse).expt(exptN(:,imouse)).tag(1).oriTuning.oriFit = vonMisesFitAllCells;
                mouse(imouse).expt(exptN(:,imouse)).tag(1).oriTuning.oriFitReliability = ...
                    fitReliability;
            else
                load(fullfile(fnTun, 'oriTuningAndFits_g.mat'));
                mouse(imouse).expt(exptN(:,imouse)).tag(1).oriTuning.oriResp = avgResponseEaOri_g;
                mouse(imouse).expt(exptN(:,imouse)).tag(1).oriTuning.oriRespSem = semResponseEaOri_g;
                mouse(imouse).expt(exptN(:,imouse)).tag(1).oriTuning.oriFit = vonMisesFitAllCells_g;
                mouse(imouse).expt(exptN(:,imouse)).tag(1).oriTuning.oriFitReliability = ...
                    fitReliability_g;
                load(fullfile(fnTun, 'oriTuningAndFits_r.mat'));
                mouse(imouse).expt(exptN(:,imouse)).tag(2).oriTuning.oriResp = avgResponseEaOri_r;
                mouse(imouse).expt(exptN(:,imouse)).tag(2).oriTuning.oriRespSem = semResponseEaOri_r;
                mouse(imouse).expt(exptN(:,imouse)).tag(2).oriTuning.oriFit = vonMisesFitAllCells_r;
                mouse(imouse).expt(exptN(:,imouse)).tag(2).oriTuning.oriFitReliability = ...
                    fitReliability_r;
            end
        end
        
        for itag = 1:2
            for iav = 1:2
                 mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(1).name = 'first stim';
                 mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(2).name = 'FA';
                 mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(3).name = 'CR, last base';
                 mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(4).name = 'target';
            end
        end
        for itag = 1:2
        %% Align data to first stim
            if itag == 1
                dataTC = dataTC_notag;
            elseif itag == 2
                dataTC = dataTC_tag;
            end
        
        ialign = 1;
        maxTrials = max(find(cLeverDown+post_event_frames-1 <  size(dataTC,1)),[],2);
        
        Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataF = zeros(1,size(dataTC,2),ntrials);
        DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        for itrial = 1:maxTrials
            Data(:,:,itrial) = dataTC(cLeverDown(itrial)-pre_event_frames:cLeverDown(itrial)+post_event_frames-1,:);
            DataF(:,:,itrial) = mean(Data(1:pre_event_frames,:,itrial),1);
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:pre_event_frames,:,itrial),1));
        end
        DataDFoverF_bl = DataDFoverF - mean(DataDFoverF(basewin,:,:),1);
        
        %identify trials with motion (large peaks in the derivative)
        if itag == 1
            ind_motion = max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>motionThreshold;
        elseif itag == 2 && isempty(dataTC_notag)
            ind_motion = max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>motionThreshold;
        end
         
        
        if maxTrials > ntrials
            exptCutoff = false(1,ntrials);
            exptCutoff(maxTrials:end) = true;
            removeTrials = ind_motion & isCatchTrial & tCyclesOn == 1 & exptCutoff;
        else
            removeTrials = ind_motion & isCatchTrial & tCyclesOn == 1;
        end
        
        for iav = 1:2
            if iav == 1
                avIndType = 0;
            elseif iav == 2
                avIndType = 1;
            end
            ind = ~removeTrials & trType == avIndType;
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).respTC = DataDFoverF_bl(:,:,ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).outcome = trOut(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).prevOutcome = prevTrOut(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).prevType = prevTrType(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).reactTime = [];
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).nCycles = tCyclesOn(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).ori = [];
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).amp = [];
        end
        
        %% Align data to false alarm and correct reject, matched for n cycles
        ialign = 2;
        
        maxTrials = max(find(cLeverDown+post_event_frames+double(cycTime*(tCyclesOn-1))-1 <  size(dataTC,1)),[],2);
        
        fa = strcmp(trOut,'failure');
        minTrFA = fa & tCyclesOn > 2;
        
        Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataF = zeros(1,size(dataTC,2),ntrials);
        DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        for itrial = 1:maxTrials
            n = tCyclesOn(itrial);
            tc_ind = (cLeverDown(itrial)-pre_event_frames:...
                cLeverDown(itrial)+post_event_frames-1) + double((n-1)*cycTime);
            Data(:,:,itrial) = dataTC(tc_ind,:);
            DataF(:,:,itrial) = mean(Data(1:pre_event_frames,:,itrial),1);
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:pre_event_frames,:,itrial),1));
        end
        DataDFoverF_bl = DataDFoverF - mean(DataDFoverF(basewin,:,:),1);
        
%         ind_motion = max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>motionThreshold;%identify trials with motion (large peaks in the derivative)
         
        
        if maxTrials > ntrials
            exptCutoff = false(1,ntrials);
            exptCutoff(maxTrials:end) = true;
            removeTrials = ind_motion & isCatchTrial & exptCutoff;
        else
            removeTrials = ind_motion & isCatchTrial;
        end
        
        for iav = 1:2
            if iav == 1
                avIndType = 0;
            elseif iav == 2
                avIndType = 1;
            end
            ind = ~removeTrials & trType == avIndType & minTrFA;
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).respTC = DataDFoverF_bl(:,:,ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).outcome = trOut(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).prevOutcome = prevTrOut(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).prevType = prevTrType(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).reactTime = reactTimeFromLastBaseMs(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).nCycles = tCyclesOn(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).ori = [];
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).amp = [];
        end
        
        ialign = 3;
        
        maxTrials = max(find(cLeverDown+post_event_frames+double(cycTime*(tCyclesOn-1))-1 <  size(dataTC,1)),[],2);
        
        hitsAndMissTr = strcmp(trOut,'success') | strcmp(trOut,'ignore');
        
        Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataF = zeros(1,size(dataTC,2),ntrials);
        DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        for itrial = 1:maxTrials
            n = tCyclesOn(itrial);
            tc_ind = (cLeverDown(itrial)-pre_event_frames:...
                cLeverDown(itrial)+post_event_frames-1) + double((n-1)*cycTime); %one stim before delivered target
            Data(:,:,itrial) = dataTC(tc_ind,:);
            DataF(:,:,itrial) = mean(Data(1:pre_event_frames,:,itrial),1);
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:pre_event_frames,:,itrial),1));
        end
        DataDFoverF_bl = DataDFoverF - mean(DataDFoverF(basewin,:,:),1);
       
%         ind_motion = max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>motionThreshold; %identify trials with motion (large peaks in the derivative)
         
        
        if maxTrials > ntrials
            exptCutoff = false(1,ntrials);
            exptCutoff(maxTrials:end) = true;
            removeTrials = ind_motion & isCatchTrial & exptCutoff;
        else
            removeTrials = ind_motion & isCatchTrial;
        end
        
        for iav = 1:2
            if iav == 1
                avIndType = 0;
            elseif iav == 2
                avIndType = 1;
            end
            ind = ~removeTrials & trType == avIndType & hitsAndMissTr;
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).respTC = DataDFoverF_bl(:,:,ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).outcome = trOut(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).prevOutcome = prevTrOut(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).prevType = prevTrType(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).reactTime = [];
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).nCycles = tCyclesOn(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).ori = [];
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).amp = [];
        end
        
        %% align to target
        ialign = 4;
        
        maxTrials = max(find(cTargetOn+post_event_frames-1 <  size(dataTC,1)),[],2);
        
        Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        for itrial = 1:maxTrials
            if ~hitsAndMissTr(itrial)
                continue
            end
            Data(:,:,itrial) = dataTC((cTargetOn(itrial)-pre_event_frames+1):cTargetOn(itrial)+post_event_frames,:);
            DataF(:,:,itrial) = mean(Data(1:pre_event_frames,:,itrial),1);
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), DataF(:,:,itrial));
        end
        DataDFoverF_bl = DataDFoverF - mean(DataDFoverF(basewin,:,:),1);
        
%         ind_motion = max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>motionThreshold; %identify trials with motion (large peaks in the derivative)
         
        if maxTrials > ntrials
            exptCutoff = false(1,ntrials);
            exptCutoff(maxTrials:end) = true;
            removeTrials = ind_motion & isCatchTrial & exptCutoff;
        else
            removeTrials = ind_motion & isCatchTrial;
        end
        
        for iav = 1:2
            if iav == 1
                avIndType = 0;
            elseif iav == 2
                avIndType = 1;
            end
            ind = ~removeTrials & trType == avIndType & hitsAndMissTr;
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).respTC = DataDFoverF_bl(:,:,ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).outcome = trOut(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).prevOutcome = prevTrOut(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).prevType = prevTrType(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).reactTime = reactTimeCalc(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).nCycles = tCyclesOn(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).ori = tGratingDirectionDeg(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).amp = tSoundTargetAmp(ind);
        end 
        
        %% align to catch trials
        ialign = 5;
        disp(expt(iexp).catch)
        if expt(iexp).catch
            cInd = ~isnan(cCatchOn);
            cCatchOn = cCatchOn(cInd);
            catchOutcome = catchOutcome(cInd);
            cCyc = cCyc(cInd);
            tCatchDirection = tCatchDirection(cInd);
            catchReactMs = catchReactMs(cInd);
            nt = length(cCatchOn);
            
            % catch target aligned
            Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),nt);
            DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),nt);
            DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),nt);
            for itrial = 1:nt
                Data(:,:,itrial) = dataTC((cCatchOn(itrial)-pre_event_frames+1):cCatchOn(itrial)+post_event_frames,:);
                DataF(:,:,itrial) = mean(Data(1:pre_event_frames,:,itrial),1);
                DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
                DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), DataF(:,:,itrial));
            end
            DataDFoverF_bl = DataDFoverF - mean(DataDFoverF(basewin,:,:),1);

%             ind_motion = max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>motionThreshold; %identify trials with motion (large peaks in the derivative)

            removeTrials = ind_motion(cInd);
            ind = ~removeTrials & ~strcmp(catchOutcome,'failure') & trType(cInd) == 1;
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(1).align(ialign).respTC = DataDFoverF_bl(:,:,ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(1).align(ialign).outcome = catchOutcome(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(1).align(ialign).reactTime = catchReactMs(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(1).align(ialign).nCycles = cCyc(ind)-1;
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(1).align(ialign).ori = tCatchDirection(ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(1).align(ialign).amp = [];
            
            % 1 stim back from catch target aligned
            cLeverDown_catch = cLeverDown(cInd);
            Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),nt);
            DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),nt);
            DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),nt);
            for itrial = 1:nt
                n = cCyc(itrial)-1;
                tc_ind = (cLeverDown_catch(itrial)-pre_event_frames:...
                    cLeverDown_catch(itrial)+post_event_frames-1) + double((n-1)*cycTime); %one stim before delivered target
            
                Data(:,:,itrial) = dataTC(tc_ind,:);
            DataF(:,:,itrial) = mean(Data(1:pre_event_frames,:,itrial),1);
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:pre_event_frames,:,itrial),1));
            end
            DataDFoverF_bl = DataDFoverF - mean(DataDFoverF(basewin,:,:),1);

%             ind_motion = max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>motionThreshold; %identify trials with motion (large peaks in the derivative)
            removeTrials = ind_motion(cInd);
            ind = ~removeTrials & ~strcmp(catchOutcome,'failure') & trType(cInd) == 1;
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(1).align(ialign).crRespTC = DataDFoverF_bl(:,:,ind);
            mouse(imouse).expt(exptN(:,imouse)).tag(itag).av(1).align(ialign).crNCycles = cCyc(ind)-2;
        end
        
        
        end
    
    %%
        %% passive expt
        if isfield(expt,'passExpt')
            if ~isempty(expt(iexp).passExpt)
                disp(['Passive Expt - ' num2str(expt(iexp).date) ' i' num2str(expt(iexp).SubNum)])

                nrun = size(expt(iexp).runs,1);

                %keep track of which mouse
                mousePass(imouse).mouse_name = expt(iexp).SubNum;
                mousePass(imouse).expt(exptN(:,imouse)).date = expt(iexp).date;
                mousePass(imouse).expt(exptN(:,imouse)).attnTask = expt(iexp).attentionTask;
                mousePass(imouse).expt(exptN(:,imouse)).hasAttn = expt(iexp).hasAttention;

                mousePass(imouse).expt(exptN(:,imouse)).info.preAlignFrames = pre_event_frames;

                %account for accumulation of frames across multiple runs 
                dataTC_notag = [];
                dataTC_tag = [];
                fnTC = fullfile(rc.ashleyAnalysis,...
                    expt(iexp).mouse,expt(iexp).folder, expt(iexp).date,'data processing');

                if cellsOrDendrites == 1
                    load(fullfile(fnTC,'timecourses_pass_cells.mat'))
                    dataTC_notag = data_pass_g_tc_subnp;
                    dataTC_tag = data_pass_r_tc_subnp;
                end

                % load and combine mworks data and timecourses
                input = loadMworksFile(expt(iexp).SubNum,expt(iexp).date,expt(iexp).passExptTime);
                tr = 1:length(input.trialOutcomeCell);

                %convert important fields to matrices
    %             run_trials = input.trialsSinceReset;
                cLeverDown = celleqel2mat_padded(input.cLeverDown);
                cFirstStim = celleqel2mat_padded(input.cFirstStim);
                cLeverUp = celleqel2mat_padded(input.cLeverUp);
                cTargetOn = celleqel2mat_padded(input.cTargetOn);
                cStimOn = celleqel2mat_padded(input.cStimOn);
                cItiStart = celleqel2mat_padded(input.cItiStart); 


                %stim timing
                if iscell(input.nFramesOn)
                    cycTime = unique(cell2mat(input.nFramesOn))+unique(cell2mat(input.nFramesOff));
                else
                    cycTime = input.nFramesOn+input.nFramesOff;
                end

                mousePass(imouse).expt(exptN(:,imouse)).info.cycTimeFrames = cycTime;

        %         reactTimes = celleqel2mat_padded(input.reactTimesMs);
                holdTimesMs = celleqel2mat_padded(input.holdTimesMs);
                tCyclesOn = double(cell2mat(input.tCyclesOn));
                nCyclesOn = double(cell2mat(input.nCyclesOn));
                holdTimesMs = holdTimesMs(tr);
                tCyclesOn = tCyclesOn(tr);
                nCyclesOn = nCyclesOn(tr);        

                frameRate = input.frameRateHz;
                cycTimeMs = cycTime./frameRate*1000;
                requiredStimTimeMs = nCyclesOn.*cycTimeMs;
                reactTimeCalc = holdTimesMs - requiredStimTimeMs;
                reactTimeFromLastBaseMs = holdTimesMs - ((tCyclesOn-1).*cycTimeMs);
                ntrials = length(tr);


                %catch trial Info
                if expt(iexp).catch
                    cCatchOn = celleqel2mat_padded(input.cCatchOn);
                    nt = length(cCatchOn);
                    tCatchDirection = celleqel2mat_padded(input.tCatchGratingDirectionDeg);

        %             trialStartFr = celleqel2mat_padded(input.cFirstStim);
                    trialEndFr = celleqel2mat_padded(input.cLeverUp);
                    cCyc = celleqel2mat_padded(input.catchCyclesOn);
                    cCyc = cCyc(tr);
                    catchReactFr = trialEndFr - cCatchOn;
                    catchCycOnMs = (cCyc-1).*cycTimeMs;
                    catchReactMs = holdTimesMs - catchCycOnMs;
                    tooFastFr = round(frameRate*((input.tooFastTimeMs)/1000));
                    reactTimeFr = round(frameRate*((input.reactTimeMs)/1000));

                    catchOutcome = cell(1,nt);
                    outInd = catchReactFr > tooFastFr & catchReactFr < reactTimeFr;
                    catchOutcome(outInd) = {'FA'};
                    outInd = catchReactFr > reactTimeFr;
                    catchOutcome(outInd) = {'CR'};
                    outInd = catchReactFr < 0;
                    catchOutcome(outInd) = {'failure'};

                    tCatchDirection = tCatchDirection(tr);
        %             catchReactMs = catchReactMs(tr);
        %             catchOutcome = catchOutcome(tr);
                end  

    %             offset = 0;
    %             for irun = 1:nrun
    %                 ImgFolder = expt(iexp).runs(irun,:);
    %                 nfr_run = nFramesSbxDataset(expt(iexp).mouse,expt(iexp).date,ImgFolder);
    %                 offset = offset+nfr_run;
    %                 if irun < nrun
    %                     startTrial = sum(run_trials(1, 1:irun),2)+1;
    %                     endTrial = sum(run_trials(1,1:irun+1),2);
    %                     cLeverDown(1,startTrial:endTrial) = cLeverDown(1,startTrial:endTrial)+offset;
    %                     cFirstStim(1,startTrial:endTrial) = cFirstStim(1,startTrial:endTrial)+offset;
    %                     cLeverUp(1,startTrial:endTrial) = cLeverUp(1,startTrial:endTrial)+offset;
    %                     cTargetOn(1,startTrial:endTrial) = cTargetOn(1,startTrial:endTrial)+offset;
    %                     cStimOn(1,startTrial:endTrial) = cStimOn(1,startTrial:endTrial)+offset;
    %                     cItiStart(1,startTrial:endTrial) = cItiStart(1,startTrial:endTrial)+offset;
    %                     if expt(iexp).catch
    %                         cCatchOn(1,startTrial:endTrial) = cCatchOn(1,startTrial:endTrial)+offset;
    %                     end
    %                 end
    %             end

                cLeverDown = cLeverDown(tr);
                cFirstStim = cFirstStim(tr);
                cLeverUp = cLeverUp(tr);
                cTargetOn = cTargetOn(tr);
                cStimOn = cStimOn(tr);
                cItiStart = cItiStart(tr);
                if expt(iexp).catch
                    cCatchOn = cCatchOn(tr);
                end

                %previous trial info
                trType = double(cell2mat(input.tBlock2TrialNumber));
                trType = trType(tr);
                trType_shift = [NaN trType];
                prevTrType = num2cell(trType_shift(1:length(trType)));
                trType_shift = [NaN NaN trType];
                prev2TrType = num2cell(trType_shift(1:length(trType)));
                for i = 1:length(trType)
                    if prevTrType{i} == 0
                        prevTrType{i} = 'vis';
                    else
                        prevTrType{i} = 'aud';
                    end
                    if prev2TrType{i} == 0
                        prev2TrType{i} = 'vis';
                    else
                        prev2TrType{i} = 'aud';
                    end
                end

                trOut = input.trialOutcomeCell;
                trOut = trOut(tr);
                trOut_shift = [{NaN} trOut];
                prevTrOut = trOut_shift(1:length(trOut));
                trOut_shift = [{NaN} {NaN} trOut];
                prev2TrOut = trOut_shift(1:length(trOut));


                if expt(iexp).catch
                    isCatchTrial = logical(cell2mat(input.tShortCatchTrial));
                else
                    isCatchTrial = false(1,ntrials);
                end
                isCatchTrial = isCatchTrial(tr);

                tGratingDirectionDeg = chop(celleqel2mat_padded(input.tGratingDirectionDeg),4);
                tGratingDirectionDeg = tGratingDirectionDeg(tr);
                Dirs = unique(tGratingDirectionDeg);

                if isfield(input,'tSoundTargetAmplitude')
                    tSoundTargetAmp = celleqel2mat_padded(input.tSoundTargetAmplitude);
                else
                    tSoundTargetAmp = ones(1,ntrials).*input.soundTargetAmplitude;
                end
                tSoundTargetAmp = tSoundTargetAmp(tr);
                Amps = unique(tSoundTargetAmp);

                mousePass(imouse).expt(exptN(:,imouse)).info.visTargets = Dirs;
                mousePass(imouse).expt(exptN(:,imouse)).info.audTargets = Amps;

        %         %load direction tuning data
        %         fnTun = fullfile(rc.ashleyAnalysis,...
        %             expt(iexp).mouse,expt(iexp).folder, expt(iexp).date, ...
        %             expt(iexp).dirtuning);
        %         if cellsOrDendrites == 1
        %             load(fullfile(fnTun, 'oriTuningAndFits.mat'));
        %         elseif cellsOrDendrites == 2
        %             load(fullfile(fnTun, 'oriTuningAndFits_den.mat'));
        %         end
        %         
        %         mouse(imouse).expt(exptN(:,imouse)).oriTuning.oriResp = avgResponseEaOri;
        %         mouse(imouse).expt(exptN(:,imouse)).oriTuning.oriRespSem = semResponseEaOri;
        %         mouse(imouse).expt(exptN(:,imouse)).oriTuning.oriFit = vonMisesFitAllCells;
        %         mouse(imouse).expt(exptN(:,imouse)).oriTuning.oriFitReliability = ...
        %             fitReliability;
                mousePass(imouse).expt(exptN(:,imouse)).tag(1).name = [expt(iexp).redChannelLabel '-'];
                mousePass(imouse).expt(exptN(:,imouse)).tag(2).name = [expt(iexp).redChannelLabel '+'];
                for itag = 1:2
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(1).name = 'visual';
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(2).name = 'auditory';
                end

                for itag = 1:2
                    for iav = 1:2
                         mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(1).name = 'first stim';
                         mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(2).name = 'FA';
                         mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(3).name = 'CR, last base';
                         mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(4).name = 'target';
                    end
                end
                for itag = 1:2
                %% Align data to first stim
                    if itag == 1
                        dataTC = dataTC_notag;
                    elseif itag == 2
                        dataTC = dataTC_tag;
                    end

                ialign = 1;
                maxTrials = max(find(cLeverDown+post_event_frames-1 <  size(dataTC,1)),[],2);

                Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
                DataF = zeros(1,size(dataTC,2),ntrials);
                DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
                DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
                for itrial = 1:maxTrials
                    Data(:,:,itrial) = dataTC(cLeverDown(itrial)-pre_event_frames:cLeverDown(itrial)+post_event_frames-1,:);
                    DataF(:,:,itrial) = mean(Data(1:pre_event_frames,:,itrial),1);
                    DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
                    DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:pre_event_frames,:,itrial),1));
                end
                DataDFoverF_bl = DataDFoverF - mean(DataDFoverF(basewin,:,:),1);

                %identify trials with motion (large peaks in the derivative)
                ind_motion = max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>motionThreshold;


                if maxTrials > ntrials
                    exptCutoff = false(1,ntrials);
                    exptCutoff(maxTrials:end) = true;
                    removeTrials = ind_motion & isCatchTrial & tCyclesOn == 1 & exptCutoff;
                else
                    removeTrials = ind_motion & isCatchTrial & tCyclesOn == 1;
                end

                for iav = 1:2
                    if iav == 1
                        avIndType = 0;
                    elseif iav == 2
                        avIndType = 1;
                    end
                    ind = ~removeTrials & trType == avIndType;
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).respTC = DataDFoverF_bl(:,:,ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).outcome = trOut(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).prevOutcome = prevTrOut(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).prevType = prevTrType(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).reactTime = [];
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).nCycles = tCyclesOn(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).ori = [];
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).amp = [];
                end

                %% Align data to false alarm and correct reject, matched for n cycles
                ialign = 2;

                maxTrials = max(find(cLeverDown+post_event_frames+double(cycTime*(tCyclesOn-1))-1 <  size(dataTC,1)),[],2);

                fa = strcmp(trOut,'failure');
                minTrFA = fa & tCyclesOn > 2;

                Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
                DataF = zeros(1,size(dataTC,2),ntrials);
                DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
                DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
                for itrial = 1:maxTrials
                    n = tCyclesOn(itrial);
                    tc_ind = (cLeverDown(itrial)-pre_event_frames:...
                        cLeverDown(itrial)+post_event_frames-1) + double((n-1)*cycTime);
                    Data(:,:,itrial) = dataTC(tc_ind,:);
                    DataF(:,:,itrial) = mean(Data(1:pre_event_frames,:,itrial),1);
                    DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
                    DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:pre_event_frames,:,itrial),1));
                end
                DataDFoverF_bl = DataDFoverF - mean(DataDFoverF(basewin,:,:),1);

                ind_motion = max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>motionThreshold;%identify trials with motion (large peaks in the derivative)


                if maxTrials > ntrials
                    exptCutoff = false(1,ntrials);
                    exptCutoff(maxTrials:end) = true;
                    removeTrials = ind_motion & isCatchTrial & exptCutoff;
                else
                    removeTrials = ind_motion & isCatchTrial;
                end

                for iav = 1:2
                    if iav == 1
                        avIndType = 0;
                    elseif iav == 2
                        avIndType = 1;
                    end
                    ind = ~removeTrials & trType == avIndType & minTrFA;
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).respTC = DataDFoverF_bl(:,:,ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).outcome = trOut(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).prevOutcome = prevTrOut(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).prevType = prevTrType(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).reactTime = reactTimeFromLastBaseMs(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).nCycles = tCyclesOn(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).ori = [];
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).amp = [];
                end

                ialign = 3;

                maxTrials = max(find(cLeverDown+post_event_frames+double(cycTime*(tCyclesOn-1))-1 <  size(dataTC,1)),[],2);

                hitsAndMissTr = strcmp(trOut,'success') | strcmp(trOut,'ignore');

                Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
                DataF = zeros(1,size(dataTC,2),ntrials);
                DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
                DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
                for itrial = 1:maxTrials
                    n = tCyclesOn(itrial);
                    tc_ind = (cLeverDown(itrial)-pre_event_frames:...
                        cLeverDown(itrial)+post_event_frames-1) + double((n-1)*cycTime); %one stim before delivered target
                    Data(:,:,itrial) = dataTC(tc_ind,:);
                    DataF(:,:,itrial) = mean(Data(1:pre_event_frames,:,itrial),1);
                    DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
                    DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:pre_event_frames,:,itrial),1));
                end
                DataDFoverF_bl = DataDFoverF - mean(DataDFoverF(basewin,:,:),1);

                ind_motion = max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>motionThreshold; %identify trials with motion (large peaks in the derivative)


                if maxTrials > ntrials
                    exptCutoff = false(1,ntrials);
                    exptCutoff(maxTrials:end) = true;
                    removeTrials = ind_motion & isCatchTrial & exptCutoff;
                else
                    removeTrials = ind_motion & isCatchTrial;
                end

                for iav = 1:2
                    if iav == 1
                        avIndType = 0;
                    elseif iav == 2
                        avIndType = 1;
                    end
                    ind = ~removeTrials & trType == avIndType & hitsAndMissTr;
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).respTC = DataDFoverF_bl(:,:,ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).outcome = trOut(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).prevOutcome = prevTrOut(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).prevType = prevTrType(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).reactTime = [];
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).nCycles = tCyclesOn(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).ori = [];
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).amp = [];
                end

                %% align to target
                ialign = 4;

                maxTrials = max(find(cTargetOn+post_event_frames-1 <  size(dataTC,1)),[],2);

                Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
                DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
                DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
                for itrial = 1:maxTrials
                    if ~hitsAndMissTr(itrial)
                        continue
                    end
                    Data(:,:,itrial) = dataTC((cTargetOn(itrial)-pre_event_frames+1):cTargetOn(itrial)+post_event_frames,:);
                    DataF(:,:,itrial) = mean(Data(1:pre_event_frames,:,itrial),1);
                    DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
                    DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), DataF(:,:,itrial));
                end
                DataDFoverF_bl = DataDFoverF - mean(DataDFoverF(basewin,:,:),1);

                ind_motion = max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>motionThreshold; %identify trials with motion (large peaks in the derivative)

                if maxTrials > ntrials
                    exptCutoff = false(1,ntrials);
                    exptCutoff(maxTrials:end) = true;
                    removeTrials = ind_motion & isCatchTrial & exptCutoff;
                else
                    removeTrials = ind_motion & isCatchTrial;
                end

                for iav = 1:2
                    if iav == 1
                        avIndType = 0;
                    elseif iav == 2
                        avIndType = 1;
                    end
                    ind = ~removeTrials & trType == avIndType & hitsAndMissTr;
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).respTC = DataDFoverF_bl(:,:,ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).outcome = trOut(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).prevOutcome = prevTrOut(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).prevType = prevTrType(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).reactTime = reactTimeCalc(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).nCycles = tCyclesOn(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).ori = tGratingDirectionDeg(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(iav).align(ialign).amp = tSoundTargetAmp(ind);
                end 

                %% align to catch trials
                ialign = 5;
                disp(expt(iexp).catch)
                if expt(iexp).catch
                    cInd = ~isnan(cCatchOn);
                    cCatchOn = cCatchOn(cInd);
                    if cCatchOn(end)+post_event_frames > size(dataTC,1)
                        cInd = cInd(1:end-1);
                        cCatchOn = cCatchOn(1:end-1);
                    end
                    catchOutcome = catchOutcome(cInd);
                    cCyc = cCyc(cInd);
                    tCatchDirection = tCatchDirection(cInd);
                    catchReactMs = catchReactMs(cInd);
                    nt = length(cCatchOn);
                    
                    % catch target aligned
                    Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),nt);
                    DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),nt);
                    DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),nt);
                    for itrial = 1:nt
                        Data(:,:,itrial) = dataTC((cCatchOn(itrial)-pre_event_frames+1):cCatchOn(itrial)+post_event_frames,:);
                        DataF(:,:,itrial) = mean(Data(1:pre_event_frames,:,itrial),1);
                        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
                        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), DataF(:,:,itrial));
                    end
                    DataDFoverF_bl = DataDFoverF - mean(DataDFoverF(basewin,:,:),1);

                    ind_motion = max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>motionThreshold; %identify trials with motion (large peaks in the derivative)

                    removeTrials = ind_motion;
                    ind = ~removeTrials & ~strcmp(catchOutcome,'failure') & trType(cInd) == 1;
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(1).align(ialign).respTC = DataDFoverF_bl(:,:,ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(1).align(ialign).outcome = catchOutcome(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(1).align(ialign).reactTime = catchReactMs(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(1).align(ialign).nCycles = cCyc(ind)-1;
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(1).align(ialign).ori = tCatchDirection(ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(1).align(ialign).amp = [];

                    % 1 stim back from catch target aligned
                    Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),nt);
                    DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),nt);
                    DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),nt);
                    for itrial = 1:nt
                        n = cCyc(itrial)-1;
                        tc_ind = (cLeverDown(itrial)-pre_event_frames:...
                            cLeverDown(itrial)+post_event_frames-1) + double((n-1)*cycTime); %one stim before delivered target

                        Data(:,:,itrial) = dataTC(tc_ind,:);
                    DataF(:,:,itrial) = mean(Data(1:pre_event_frames,:,itrial),1);
                    DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
                    DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:pre_event_frames,:,itrial),1));
                    end
                    DataDFoverF_bl = DataDFoverF - mean(DataDFoverF(basewin,:,:),1);

                    ind_motion = max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>motionThreshold; %identify trials with motion (large peaks in the derivative)
                    removeTrials = ind_motion;
                    ind = ~removeTrials & ~strcmp(catchOutcome,'failure') & trType(cInd) == 1;
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(1).align(ialign).crRespTC = DataDFoverF_bl(:,:,ind);
                    mousePass(imouse).expt(exptN(:,imouse)).tag(itag).av(1).align(ialign).crNCycles = cCyc(ind)-2;
                end
                end
            end
        end
    end
    if cellsOrDendrites == 1
        save(fullfile(rc.caOutputDir, dataGroup, [mouse_str '_trOutcomeStruct_cells' datasetStr(5:end) '.mat']), 'mouse','mousePass','-v7.3');
%         print(fullfile(rc.caOutputDir, dataGroup, [datasetStr(5:end) '_motionHist.pdf']), '-dpdf')
    elseif cellsOrDendrites == 2
        save(fullfile(rc.caOutputDir, dataGroup, [mouse_str '_trOutcomeStruct_dendrites' datasetStr(5:end) '.mat']), 'mouse','mousePass','-v7.3');
%         print(fullfile(rc.caOutputDir, dataGroup, [datasetStr(5:end) '_motionHist.pdf']), '-dpdf')
    end
end
        
        
        

