function [RTs, JuiceLicks] = findBursts(JuiceLicks)
burstLimit = 0.5;  %The time window for three licks to occur to be considered a burst
responseWinMax = 3; % The time window the lick has to occur in to be considered a response
LicksAfterToneStruct = struct;
for i = 1:length(JuiceLicks(:,1))
    LicksAfterToneIndex = (JuiceLicks{i,2}>JuiceLicks{i,1});   %take only licks after the tone
    allLicks = JuiceLicks{i,2};
    licksAfterTone = allLicks(LicksAfterToneIndex);
    %LicksAfterToneStruct(i,1) = licksAfterTone;
    %TrialLicksAfterTone = allLicksAfterTone{i,1};
    if (length(licksAfterTone)>3);
        if ((licksAfterTone(3)-licksAfterTone(1)) <= burstLimit) && ((licksAfterTone(1)-JuiceLicks{i,1}) < responseWinMax) %check if there is a burst of licks after the tone & if the burst is within the response window
            JuiceLicks{i,4}=licksAfterTone(1); % first lick into JuiceLicks, column 4
        else
            JuiceLicks{i,4}=NaN; %if there isn't a burst
        end
    else
        JuiceLicks{i,4} = NaN; % if there are less than three licks after juice
    end
end
FirstLicks = cell2mat(JuiceLicks(:,4));
RTs = (FirstLicks - cell2mat(JuiceLicks(:,1)));
end

    