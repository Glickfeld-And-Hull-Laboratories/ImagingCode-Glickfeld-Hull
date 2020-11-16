function [n,n_trTypeID] = getRepeatedTrialsN(trialIndex,stimulusID,choiceID)

%false alarms
ind = stimulusID == 0 & choiceID == 1;
n_trials = histcounts(trialIndex(ismember(trialIndex,find(ind))),1:max(trialIndex));
n(1) = sum(n_trials(n_trials>1)-1);
n_trTypeID{1} = 'fa';

%correct rejects
ind = stimulusID == 0 & choiceID == 0;
n_trials = histcounts(trialIndex(ismember(trialIndex,find(ind))),1:max(trialIndex));
n(2) = sum(n_trials(n_trials>1)-1);
n_trTypeID{2} = 'cr';

%hits
ind = stimulusID == 1 & choiceID == 1;
n_trials = histcounts(trialIndex(ismember(trialIndex,find(ind))),1:max(trialIndex));
n(3) = sum(n_trials(n_trials>1)-1);
n_trTypeID{3} = 'h';

%misses
ind = stimulusID == 1 & choiceID == 0;
n_trials = histcounts(trialIndex(ismember(trialIndex,find(ind))),1:max(trialIndex));
n(4) = sum(n_trials(n_trials>1)-1);
n_trTypeID{4} = 'm';

end