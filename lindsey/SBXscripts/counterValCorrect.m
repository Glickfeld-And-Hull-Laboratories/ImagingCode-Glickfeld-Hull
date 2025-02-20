function input_correct = counterValCorrect(input,stimOffs)
    %reassigns frames according to photodiode stim off time on each trial
    nTrials = length(stimOffs);
    input_correct = input;
    for i = 1:nTrials
        if i == 1
            input_correct.counterTimesUs{i} = input.counterTimesUs{i}(find(input.counterValues{i}==1,1,'last'):end);
            input_correct.counterValues{i} = double(1:size(input_correct.counterTimesUs{i},2));
        else
            end_val = stimOffs(i)-1;
            start_val = end_val-size(input_correct.counterTimesUs{i},2)+1;
            input_correct.counterValues{i} = double(start_val:end_val);
        end
    end
end