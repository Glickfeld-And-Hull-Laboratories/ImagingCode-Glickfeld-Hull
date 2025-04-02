data_dfof_runOnset_match = cell(1,nd)
for id = 1:nd
    fwdWheelClean = wheel_speed_clean{id};
    fwdWheelClean(fwdWheelClean<0)=0;
    % moveMean = movmean(fwdWheelClean,3);
    % smoothRun=smooth(fwdWheelClean,3);
    % figure;plot(fwdWheelClean(900:1300));hold on;plot(smoothRun(900:1300));hold on;plot(moveMean(900:1300));
    
    % find and refine the onset indices
    stillFrames = logical(fwdWheelClean<2); %find the stationary periods based on the sliding window
    
    runOnsets=(find(diff(stillFrames) == -1)+1); %find the transitions from stationary to running
    runOffsets=(find(diff(stillFrames) == 1)+1); %find the transitions from running to stationary
  
    %subset to onsets that are preceeded by 1s stillness and followed by at
    %least 1s of running
    runOnsets_clean = runOnsets;
    toDrop = [];
    for iOnset = 1:length(runOnsets_clean)
        testWinStill=[runOnsets_clean(iOnset)-frame_rate:runOnsets_clean(iOnset)];
        testWinRun=[runOnsets_clean(iOnset):runOnsets_clean(iOnset)+frame_rate];
        if sum(stillFrames(testWinStill))<12 |sum(stillFrames(testWinRun))<12
            toDrop=[toDrop,iOnset]; %gather a list of onsets that don't meet my criteria
        end 
    end
    runOnsets_clean(toDrop)=[];
    %%
    
    %make a vector that is 1 during the ITI and 0 when the stim is on
    ITI = ones(1,nFrames);
    for iTrial = 1:nTrials(id)
        ITI(cStimOn(iTrial):cStimOn(iTrial)+nOn)=0;
    end
    
    %subset the running onsets to only include those during the ITI
    ITIOnsets = runOnsets_clean;
    ITIOnsets(ITI(runOnsets_clean)==0)=[];
    
    % get neural data for run onsets 
    % extract timecourses during a window around the run onset
    onsetWin = 1; %how many seconds you want the window the be, on EITHER side of the onset
    
    data_dfof_runOnset = nan(onsetWin*2*frame_rate,nCells,length(ITIOnsets));
    runConfirmation = nan(onsetWin*2*frame_rate,length(ITIOnsets));
    for iOnset = 1:length(ITIOnsets)
        %find inds for a window around the onset
        fullWindow = [(ITIOnsets(iOnset)-(onsetWin*frame_rate)):(ITIOnsets(iOnset)+(onsetWin*frame_rate)-1)];
        bslnWindow = [(ITIOnsets(iOnset)-(frame_rate/2)):ITIOnsets(iOnset)];%the last half second before the run onset. Will use this as the baseline for df/f
        tempFull=cellTCs_match{id}(fullWindow,:);
        tempBsln=mean(cellTCs_match{id}(bslnWindow,:));
        %this is F, convert to dfof. 
        tempDfof = bsxfun(@rdivide,bsxfun(@minus,tempFull,tempBsln),tempBsln);
        data_dfof_runOnset(:,:,iOnset)=tempDfof;
        runConfirmation(:,iOnset) = fwdWheelClean(fullWindow);
    end
    
    figure; plot(nanmean(runConfirmation,2)) %to check whether running actually does increase around the time of these onsets.
    
    data_dfof_runOnset_match{id} = mean(data_dfof_runOnset,3,'omitmissing'); %frames x cells (for all matched cells) averaged over all the onsets
end
