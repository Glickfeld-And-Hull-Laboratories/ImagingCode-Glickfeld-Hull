%% finding trial types
%find trials within nframes
lastTrial = find(cell2mat(input.cTrialEnd)<min(nframes),1,'last');

%these are indices of trials
hit =  find(strcmp(input.trialOutcomeCell(1:lastTrial), 'success'));
FA =  find(strcmp(input.trialOutcomeCell(1:lastTrial), 'failure'));
miss =  find(strcmp(input.trialOutcomeCell(1:lastTrial), 'ignore'));
block2trials = find(cell2mat(input.tBlock2TrialNumber(1:lastTrial)));
block1trials = find(~cell2mat(input.tBlock2TrialNumber(1:lastTrial)));

hitRate = nan(nCon,2);

%get the mean neural response to hits and misses at each contrast for each
%block. Currently looking 1 second before through 1 second after the stim
%turns on
trialAvrg_tcs = nan(60,nCells,nCon,4); %empty matrix, 60 frames by nCells by nCon by hit B1, miss B1, hit B2, miss B2
for iCon = 1:nCon
    conInds=find(trialCon==cons(iCon)); %find trials for this contrast
    block1Inds = intersect(conInds,block1trials); %find trials for block 1 at this contrast
    hitInds = intersect(hit, block1Inds); %hit trials for block 1 at this contrast
    missInds = intersect(miss,block1Inds); %miss trials for block1 at this contrast
    hitTargetTimes = cell2mat(input.cTargetOn(hitInds)); %find when the target comes on for the hit trials for this contrast, in frames
    missTargetTimes = cell2mat(input.cTargetOn(missInds)); %find when the target comes on for the miss trials for this contrast, in frames
    hitsTemp = nan(60,nCells,length(hitInds));
    missTemp = nan(60,nCells,length(missInds));
    for iHit = 1:length(hitInds) %loop through the hit trials for this contrast
        trialInd = max(find(cStart<double(hitTargetTimes(iHit)))); %find the lastest cStart that's less than the target time - this is the index of the trial in question
        baseWindow = [cStart(trialInd)-15, cStart(trialInd)]; %define a baseline window that is 15 frames before the start of the trial up to the start of the trial
        hitsTemp(:,:,iHit)=npSub_tc(hitTargetTimes(iHit)-29:hitTargetTimes(iHit)+30,:)-mean(npSub_tc(baseWindow,:),1);
        %hitsTemp(:,:,iHit)=npSub_tc(hitTargetTimes(iHit)-29:hitTargetTimes(iHit)+30,:)-mean(npSub_tc(hitTargetTimes(iHit)-15:hitTargetTimes(iHit),:),1);
        %to get dfof, take the np subtracted timecourse for t=-1 before stim through 1
        %second, and subtract the mean of the -0.5 through 0 (using the
        %last half second before stim onset as the baseline F
    end
    trialAvrg_tcs(:,:,iCon,1)=mean(hitsTemp,3);
    for iMiss = 1:length(missInds) %loop through the miss trials for this contrast
        trialInd = max(find(cStart<double(hitTargetTimes(iMiss))));
        baseWindow = [cStart(trialInd)-15, cStart(trialInd)]; %define a baseline window that is 15 frames before the start of the trial up to the start of the trial
        hitsTemp(:,:,iMiss)=npSub_tc(hitTargetTimes(iMiss)-29:hitTargetTimes(iMiss)+30,:)-mean(npSub_tc(baseWindow,:),1);
    end
    trialAvrg_tcs(:,:,iCon,2)=mean(missTemp,3);
    hitRate(iCon, 1) = length(hitInds)/(length(hitInds)+length(missInds));

    %now do block 2
    block2Inds = intersect(conInds,block2trials); %find trials for block 1 at this contrast
    hitInds = intersect(hit, block2Inds); %hit trials for block 1 at this contrast
    missInds = intersect(miss,block2Inds); %miss trials for block1 at this contrast
    hitTargetTimes = cell2mat(input.cTargetOn(hitInds));
    missTargetTimes = cell2mat(input.cTargetOn(missInds));
    hitsTemp = nan(60,nCells,length(hitInds));
    missTemp = nan(60,nCells,length(missInds));
    for iHit = 1:length(hitInds) %loop through the hit trials for this contrast
        trialInd = max(find(cStart<double(hitTargetTimes(iHit)))); %find the lastest cStart that's less than the target time - this is the index of the trial in question
        baseWindow = [cStart(trialInd)-15, cStart(trialInd)]; %define a baseline window that is 15 frames before the start of the trial up to the start of the trial
        hitsTemp(:,:,iHit)=npSub_tc(hitTargetTimes(iHit)-29:hitTargetTimes(iHit)+30,:)-mean(npSub_tc(baseWindow,:),1);
        %hitsTemp(:,:,iHit)=npSub_tc(hitTargetTimes(iHit)-29:hitTargetTimes(iHit)+30,:)-mean(npSub_tc(hitTargetTimes(iHit)-15:hitTargetTimes(iHit),:),1);
        %to get dfof, take the np subtracted timecourse for t=-1 before stim through 1
        %second, and subtract the mean of the -0.5 through 0 (using the
        %last half second before stim onset as the baseline F
    end
    trialAvrg_tcs(:,:,iCon,1)=mean(hitsTemp,3);
    for iMiss = 1:length(missInds) %loop through the miss trials for this contrast
        trialInd = max(find(cStart<double(hitTargetTimes(iMiss))));
        baseWindow = [cStart(trialInd)-15, cStart(trialInd)]; %define a baseline window that is 15 frames before the start of the trial up to the start of the trial
        hitsTemp(:,:,iMiss)=npSub_tc(hitTargetTimes(iMiss)-29:hitTargetTimes(iMiss)+30,:)-mean(npSub_tc(baseWindow,:),1);
    end
    trialAvrg_tcs(:,:,iCon,4)=mean(missTemp,3);
    hitRate(iCon, 2) = length(hitInds)/(length(hitInds)+length(missInds));

end
save(fullfile(fnout,'trialAvrg_tcs.mat'),'trialAvrg_tcs')