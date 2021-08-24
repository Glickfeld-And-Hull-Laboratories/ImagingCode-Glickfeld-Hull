function [RTs, JuiceLicks] = findRT(JuiceLicks)
%some variables you might want to change
burstLimit = 0.5;  %The time window for three licks to occur to be considered a burst
responseWinMax = 3; % The time window the lick has to occur in to be considered a response
delayTime= .698; %the time between the tone and the juice, in seconds

%setting time of tones & empty licks after tone structure
ToneTimes = cell2mat(JuiceLicks(:,1))-delayTime;
%LicksAfterToneStruct = struct;

%for every tone, checking for a burst of licking in the response window
for i = 1:length(JuiceLicks(:,1))
    allLicks = JuiceLicks{i,2};
    licksAfterTone = allLicks(allLicks > ToneTimes(i));   %take only licks after the tone
    
    %LicksAfterToneStruct(i,1) = licksAfterTone;
    %TrialLicksAfterTone = allLicksAfterTone{i,1};
    if (length(licksAfterTone) >= 3)
        if ((licksAfterTone(3)-licksAfterTone(1)) <= burstLimit) && ((licksAfterTone(1) - ToneTimes(i)) < responseWinMax) %check if there is a burst of licks after the tone & if the burst is within the response window
            JuiceLicks{i,4}=licksAfterTone(1); % first lick into JuiceLicks, column 4
        else
            JuiceLicks{i,4}=NaN; %if there isn't a burst
        end
    else
        JuiceLicks{i,4} = NaN; % if there are less than three licks after juice in the analysis window from FindLikcsStateMach (not just the response window)
    end
end
FirstLicks = cell2mat(JuiceLicks(:,4));
RTs = (FirstLicks - ToneTimes);
end

    