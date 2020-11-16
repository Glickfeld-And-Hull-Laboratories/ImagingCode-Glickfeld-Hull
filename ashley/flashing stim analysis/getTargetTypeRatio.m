function [r,r_trTypeID] = getTargetTypeRatio(trialIndex,stimulusID,choiceID,targetTypeID)

%hits
ind_trType = choiceID == 1 & stimulusID == 1;
ind_trUsed = trialIndex(ismember(trialIndex,find(ind_trType)));
r(1) = sum(targetTypeID(ind_trUsed) == 2)./sum(targetTypeID(ind_trUsed) == 3);
r_trTypeID{1} = 'h';
%misses
ind_trType = choiceID == 0 & stimulusID == 1;
ind_trUsed = trialIndex(ismember(trialIndex,find(ind_trType)));
r(2) = sum(targetTypeID(ind_trUsed) == 2)./sum(targetTypeID(ind_trUsed) == 3);
r_trTypeID{2} = 'm';

end