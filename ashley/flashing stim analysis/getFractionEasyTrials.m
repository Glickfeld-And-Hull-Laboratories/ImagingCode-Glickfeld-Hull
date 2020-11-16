function [f,f_trTypeID] = getFractionEasyTrials(trialIndex,stimulusID,choiceID,targetTypeID)

%hits
ind_trType = choiceID == 1 & stimulusID == 1;
ind_trUsed = trialIndex(ismember(trialIndex,find(ind_trType)));
f(1) = sum(targetTypeID(ind_trUsed) == 3)./length(targetTypeID(ind_trUsed));
f_trTypeID{1} = 'h';
%misses
ind_trType = choiceID == 0 & stimulusID == 1;
ind_trUsed = trialIndex(ismember(trialIndex,find(ind_trType)));
f(2) = sum(targetTypeID(ind_trUsed) == 3)./length(ind_trUsed);
f_trTypeID{2} = 'm';

end