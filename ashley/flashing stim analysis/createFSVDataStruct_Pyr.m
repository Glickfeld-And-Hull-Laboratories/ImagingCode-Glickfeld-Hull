function mouse = createFSVDataStruct_Pyr(datasetStr,cellsOrDendrites)
%cellsOrDendrites: 1 == cells; 2 == dendrites
    % set analysis windows
    pre_event_time = 1000; %ms
    post_event_time = 4500; %ms
    resp_win_time = 100; %ms
    pre_win_time = [-30 70];
    trans_win_time = [150 250]; %this is actually ~166-266 ms at 30 Hz
    minTrialLengthMs = 2500; % time in ms of trial rather than hard-coding cycle number for analysis
    
    if contains(datasetStr,'naive')
        bxExpt = false;
    else
        bxExpt = true;
    end
    
    grnLabel = {'EMX';'no tag';'untagged'};
    
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
    
    if bxExpt && isfield(expt,'passExpt')
        mousePass = struct;
    end
    
    
    %collect data from each experiment
    for iexp = 1:size(expt,2)
        disp([num2str(expt(iexp).date) ' i' num2str(expt(iexp).SubNum)])
        
        if contains(expt(iexp).indicator{1},'tg')
            grnID = 1;
        elseif ~isempty(expt(iexp).redChannelLabel)
            grnID = 2;
        else
            grnID = 3;
        end

        nrun = size(expt(iexp).runs,1);
        
        %keep track of which mouse
        imouse = find(values == str2num(expt(iexp).SubNum));
        mouse(imouse).mouse_name = expt(iexp).SubNum;
        exptN(:,imouse) = exptN(:,imouse)+1;
        mouse(imouse).expt(exptN(:,imouse)).date = expt(iexp).date;
        
        
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
            if grnID == 1 || grnID == 3
                dataTC_notag = data_bx_tc_subnp;
            elseif grnID == 2 
                dataTC_notag = data_bx_g_tc_subnp;
            end
        else
            error('dendrite data not yet selected')
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
        elseif isempty(expt(iexp).trial_range)
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
                
        offset = 0;
        for irun = 1:nrun
            ImgFolder = expt(iexp).runs(irun,:);
            nfr_run = nFramesSbxDataset(expt(iexp).mouse,expt(iexp).date,ImgFolder);
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
            end
        end
        
        cLeverDown = cLeverDown(tr);
        cFirstStim = cFirstStim(tr);
        cLeverUp = cLeverUp(tr);
        cTargetOn = cTargetOn(tr);
        cStimOn = cStimOn(tr);
        cItiStart = cItiStart(tr);        

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
                
        mouse(imouse).expt(exptN(:,imouse)).info.visTargets = Dirs;
        mouse(imouse).expt(exptN(:,imouse)).info.fsavSize = input.gratingHeightDeg;
        
        %load direction tuning data
        fnTun = fullfile(rc.ashleyAnalysis,...
            expt(iexp).mouse,expt(iexp).folder, expt(iexp).date, ...
            expt(iexp).dirtuning);
        if cellsOrDendrites == 1
            if exist(fullfile(fnTun, 'oriTuningAndFits.mat'))
                load(fullfile(fnTun, 'oriTuningAndFits.mat'));
                mouse(imouse).expt(exptN(:,imouse)).oriTuning.oriResp = avgResponseEaOri;
                mouse(imouse).expt(exptN(:,imouse)).oriTuning.oriRespSem = semResponseEaOri;
                mouse(imouse).expt(exptN(:,imouse)).oriTuning.oriFit = vonMisesFitAllCells;
                mouse(imouse).expt(exptN(:,imouse)).oriTuning.oriFitReliability = ...
                    fitReliability;
            elseif exist(fullfile(fnTun, 'oriTuningAndFits_g.mat'))
                load(fullfile(fnTun, 'oriTuningAndFits_g.mat'));
                mouse(imouse).expt(exptN(:,imouse)).oriTuning.oriResp = avgResponseEaOri_g;
                mouse(imouse).expt(exptN(:,imouse)).oriTuning.oriRespSem = semResponseEaOri_g;
                mouse(imouse).expt(exptN(:,imouse)).oriTuning.oriFit = vonMisesFitAllCells_g;
                mouse(imouse).expt(exptN(:,imouse)).oriTuning.oriFitReliability = ...
                    fitReliability_g;
            else
                error('No tuning data')
            end
        elseif cellsOrDendrites == 2
            load(fullfile(fnTun, 'oriTuningAndFits_den.mat'));
        end

        mouse(imouse).expt(exptN(:,imouse)).tagname = grnLabel{grnID};      
        
       
         mouse(imouse).expt(exptN(:,imouse)).align(1).name = 'first stim';
         mouse(imouse).expt(exptN(:,imouse)).align(2).name = 'FA';
         mouse(imouse).expt(exptN(:,imouse)).align(3).name = 'CR, last base';
         mouse(imouse).expt(exptN(:,imouse)).align(4).name = 'target';

        %% Align data to first stim
        dataTC = dataTC_notag;
        
        if isempty(dataTC)
            continue
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
        

        avIndType = 0;

        ind = ~removeTrials & trType == avIndType;
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).respTC = DataDFoverF_bl(:,:,ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).outcome = trOut(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).prevOutcome = prevTrOut(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).prevType = prevTrType(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).reactTime = [];
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).nCycles = tCyclesOn(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).ori = [];
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).amp = [];
        
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
                  
        ind = ~removeTrials & trType == avIndType & minTrFA;
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).respTC = DataDFoverF_bl(:,:,ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).outcome = trOut(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).prevOutcome = prevTrOut(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).prevType = prevTrType(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).reactTime = reactTimeFromLastBaseMs(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).nCycles = tCyclesOn(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).ori = [];
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).amp = [];
        
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
        
        ind = ~removeTrials & trType == avIndType & hitsAndMissTr;
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).respTC = DataDFoverF_bl(:,:,ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).outcome = trOut(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).prevOutcome = prevTrOut(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).prevType = prevTrType(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).reactTime = [];
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).nCycles = tCyclesOn(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).ori = [];
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).amp = [];
        
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
        
        ind = ~removeTrials & trType == avIndType & hitsAndMissTr;
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).respTC = DataDFoverF_bl(:,:,ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).outcome = trOut(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).prevOutcome = prevTrOut(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).prevType = prevTrType(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).reactTime = reactTimeCalc(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).nCycles = tCyclesOn(ind);
        mouse(imouse).expt(exptN(:,imouse)).align(ialign).ori = tGratingDirectionDeg(ind);
        
        
        end
    
    
        if ~exist(fullfile(rc.lindseyAnalysis, dataGroup))
            mkdir(fullfile(rc.lindseyAnalysis, dataGroup))
        end
    if cellsOrDendrites == 1
        if isfield(expt,'passExpt')
            save(fullfile(rc.lindseyAnalysis, dataGroup, ['trOutcomeStruct_cells.mat']), 'mouse','mousePass','-v7.3');
        else
            save(fullfile(rc.lindseyAnalysis, dataGroup, ['trOutcomeStruct_cells.mat']), 'mouse','-v7.3');
        end
%         print(fullfile(rc.caOutputDir, dataGroup, [datasetStr(5:end) '_motionHist.pdf']), '-dpdf')
    elseif cellsOrDendrites == 2
        save(fullfile(rc.lindseyAnalysis, dataGroup, ['trOutcomeStruct_dendrites.mat']), 'mouse','mousePass','-v7.3');
%         print(fullfile(rc.caOutputDir, dataGroup, [datasetStr(5:end) '_motionHist.pdf']), '-dpdf')
    end
end
        
        
        

