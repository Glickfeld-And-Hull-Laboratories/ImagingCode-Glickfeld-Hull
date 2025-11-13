function stimOns=trialTimingCalc(MWStruct,SbxStruct,id,nd)
% Extracts stimOns and stimOffs based on the MWStruct regardless of
% modification.
% MWStruct is "input", PDStruct is "info".

stimOns = cell(1,2);


    if isfield(MWStruct,'stimTimingSource')
        fprintf('Found previously corrected stim on timings for day %i.\n',id);
        switch MWStruct.stimTimingSource
            case 'MW'
                stimOns = cell2mat(MWStruct.stimOns_mwCounter);
                disp('Using mWorks counter.\n');
            case 'PD'
                stimOns = MWStruct.stimOns_photodiode;
                disp('Using photodiode onsets.\n');
        end
    else
        fprintf('Input struct was not previously corrected for day %i.\n',id);
        if isfield(SbxStruct,'frame')
            [stimOns,~] = photoFrameFinder_Sanworks(SbxStruct.frame);
            disp('Calculating stimOns from photodiode data.\n');
        else
            input_correct = counterValCorrect_noPhotodiode(input);
            stimOns = cell2mat(input_correct.cStimOn);
            disp('Calculating stimOns from MWCounter.\n');
        end
    end
