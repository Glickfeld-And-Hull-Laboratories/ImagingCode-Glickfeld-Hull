function [out_input,stimOns,stimOffs] = trialIndexIntegration(mwkStrct,sbxStrct,expectedTrialLength,tolerance,trials,action)
% outputStruct = trialIndexIntegration(input,info,expectedTrialLength,tolerance,trials,action)
% This function 1) uses photodiode info when available or use LG's script
% to find the correct timing for each trial, then adds the information into
% a new field; 2) deletes manually designated trials if chosen to do so.
% This is also a partial wrapper for the photoframeFinder_sanworks and trialDropper function.

if isempty(expectedTrialLength)
    expectedTrialLength = 90;
end

if isempty(tolerance)
    tolerance = 2;
end

tLengthRange = [expectedTrialLength-tolerance:expectedTrialLength+tolerance];

if nargin < 5
    disp('No assigned trials to be dropped, displaying counter & photodiode troubleshooting plots.\n');
    [~,~,nMwksTrials,nPtdTrials] = plotCounterIntervals(mwkStrct,sbxStrct);
    fprintf('Mworks counter had %i trials, photodiode saw %i trials.\n',nMwksTrials,nPtdTrials);
    [stimOns,stimOffs] = photoFrameFinder_Sanworks(sbxStrct.frame);
    dstimOns = diff(stimOns);
    
end



end