function [input_correct x] = counterValCorrect(input,photoStimOn);
    %reassigns frames according to photodiode stim off time on each trial
    cStimOn = celleqel2mat_padded(input.cStimOn);
    nTrials = length(cStimOn);
    input_correct = input;
    for i = 1:nTrials
        x(i) = photoStimOn(i) - cStimOn(i);
        if x(i) ~= 0
            input_correct.counterValues{i} = input.counterValues{i} + x(i);
            input_correct.cStimOn{i} = cStimOn(i) + x(i);
        end
    end
end