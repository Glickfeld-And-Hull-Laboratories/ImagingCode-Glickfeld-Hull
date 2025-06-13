function [input_correct] = counterValCorrect_noPhotodiode(input);
    %reassigns frames according to photodiode stim off time on each trial
    input_correct = input;
    if ~isfield(input,'cStimOn')
        input_correct.cStimOn = input.tGratingContrast;
        nTrials = size(input.tGratingContrast,2);
        nOn = input.nScansOn;
        nOff = input.nScansOff;
        for i = 1:nTrials
            input_correct.cStimOn{i} = ((i-1)*(nOn+nOff))+nOff+1;
        end
    end
    
    nTrials = length(input.counterValues);
    frame_int = mean(diff(input.counterTimesUs{2}));
    for iT = 1:nTrials
        vals = input_correct.counterValues{iT};
        times = input_correct.counterTimesUs{iT};
        if iT == 1
            ind = find(vals == 1,1,'last');
            vals = vals(ind:end);
            times = times(ind:end);
        end
        dtimes = diff(times);
        dvals = diff(vals);
        ind_skip_int = find(dtimes>1.75*frame_int);
        ind_skip_fr = find(dvals==2);
        ind_skip_multfr = find(dvals>2);
        ind_mult_fr = find(dtimes<.75*frame_int);

        frameErr_skip = setdiff(ind_skip_int,ind_skip_fr);
        frameErr_mult = setdiff(ind_skip_fr,[ind_mult_fr ind_skip_int]);
        frameErr_add = [ind_skip_multfr intersect(ind_skip_fr,ind_mult_fr)];
        
        %frameErr_mult identifies frames that were added (and counted) in Mworks
        %frameErr_skip identifies frames that were skipped by Mworks
        %frameErr_add identifies frames that were added but not counted
        if ~isempty(ind_skip_multfr)
            for i = 1:length(ind_skip_multfr)
                if find(abs(ind_mult_fr-ind_skip_multfr(i))<2)
                    ind_mult_fr(find(abs(ind_mult_fr-ind_skip_multfr(i))<2)) = [];
                end
            end
        end
        
        for i = 1:length(ind_mult_fr)
            if isempty(find(ind_skip_fr == ind_mult_fr(i)))
                if dtimes(ind_mult_fr(i)-1)<1.25*frame_int
                    if i>1
                        if ind_mult_fr(i)-ind_mult_fr(i-1) == 1
                            if abs(sum(dtimes(ind_mult_fr(i-1:i)))-frame_int)>0.25*frame_int
                                frameErr_mult = [frameErr_mult ind_mult_fr(i)];
                            end
                        else
                            frameErr_mult = [frameErr_mult ind_mult_fr(i)];
                        end
                    else
                        frameErr_mult = [frameErr_mult ind_mult_fr(i)];
                    end
                end
            end
        end
        
        errs = unique([frameErr_skip frameErr_mult frameErr_add]);

        for i = 1:length(errs)
            %add in skipped frames
            if ~isempty(find(frameErr_skip == errs(i))) && isempty(find(frameErr_add == errs(i)))
                ind = errs(i);
                input_correct.counterValues{iT} = [input_correct.counterValues{iT}(1:ind) ind+1 input_correct.counterValues{iT}(ind+1:end)+1];
                input_correct.counterTimesUs{iT} = [input_correct.counterTimesUs{iT}(1:ind) (input_correct.counterTimesUs{iT}(ind)+input_correct.counterTimesUs{ii}(1+ind))./2 input_correct.counterTimesUs{ii}(ind+1:end)];
                fprintf(['added trial ' num2str(iT) ' frame ' num2str(ind) '\n'])
                if ind < input_correct.cStimOn{iT}
                    input_correct.cStimOn{iT} = input_correct.cStimOn{iT}+1;
                end
                for iii = iT+1:nTrials
                    input_correct.counterValues{iii} = input_correct.counterValues{iii}+1;
                    input_correct.cStimOn{iii} = input_correct.cStimOn{iii}+1;
                end
            elseif ~isempty(find(frameErr_mult == errs(i)))
                %remove multiple frames
                ind = errs(i)+1;
                if ind<length(vals)
                    input_correct.counterValues{iT} = [input_correct.counterValues{iT}(1:ind-1) input_correct.counterValues{iT}(ind+1:end)-1];
                else
                    input_correct.counterValues{iT} = [input_correct.counterValues{iT}(1:ind-1)];
                end
                input_correct.counterTimesUs{iT}(ind) = [];
                fprintf(['removed trial ' num2str(iT) ' frame ' num2str(ind) '\n'])
                if ind < input_correct.cStimOn{iT}
                    input_correct.cStimOn{iT} = input_correct.cStimOn{iT}-1;
                end
                for iii = iT+1:nTrials
                    input_correct.counterValues{iii} = input_correct.counterValues{iii}-1;
                    input_correct.cStimOn{iii} = input_correct.cStimOn{iii}-1;
                end
            elseif ~isempty(find(frameErr_add == errs(i)))
                ind = errs(i);
                addn = dvals(ind);
                input_correct.counterValues{iT} = [input_correct.counterValues{iT}(1:ind) input_correct.counterValues{iT}(ind+2:end)-addn];
                input_correct.counterTimesUs{iT}(ind+1) = [];
                fprintf(['removed trial ' num2str(iT) ' frame ' num2str(ind) ' to ' num2str(ind+addn) '\n'])
                if ind < input_correct.cStimOn{iT}
                    input_correct.cStimOn{iT} = input_correct.cStimOn{iT}-addn;
                end
                for iii = iT+1:nTrials
                    input_correct.counterValues{iii} = input_correct.counterValues{iii}-addn;
                    input_correct.cStimOn{iii} = input_correct.cStimOn{iii}-addn;
                end
            end
        end
    end
