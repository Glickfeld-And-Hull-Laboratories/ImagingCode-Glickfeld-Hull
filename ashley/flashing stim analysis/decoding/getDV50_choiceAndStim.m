function trialSubsetIndex = getDV50_choiceAndStim(...
    stimulusTrialID,choiceTrialID,nTrials,minTrials)
    
    if mod(nTrials,2) == 1
        nTrialsEaSubset = [ceil(nTrials/2),floor(nTrials/2)];
    else
        nTrialsEaSubset = [nTrials/2,nTrials/2];
    end
    if mean(stimulusTrialID) == 0.5
        trialSubsetIndex = [];
        for i = 0:1
            thisStimNo = find(stimulusTrialID == i & choiceTrialID == 0);
            thisStimYes = find(stimulusTrialID == i & choiceTrialID == 1); 
            if isempty(thisStimNo) || length(thisStimNo) < minTrials
                trialSubsetIndex = nan;
                break
            elseif isempty(thisStimYes) || length(thisStimYes) < minTrials
                trialSubsetIndex = nan;
                break
            else
                trialSubsetIndex = cat(1,trialSubsetIndex, ...
                    randsample(thisStimNo,...
                    ceil(nTrialsEaSubset(i+1)/2),1), ...
                    randsample(thisStimYes,...
                    floor(nTrialsEaSubset(i+1)/2),1));
            end
        end
    else
        trialSubsetIndex = [];
    end

end