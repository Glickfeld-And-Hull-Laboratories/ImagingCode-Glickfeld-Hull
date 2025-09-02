function [times vals] = plotCounterIntervals(input,info);
    nTrials = length(input.counterValues);
    vals = [];
    times = [];
    ncount = [];
    for i = 1:nTrials
        vals = [vals input.counterValues{i}];
        times = [times input.counterTimesUs{i}];
        ncount = [ncount length(input.counterTimesUs{i})];
    end
    ind = find(vals == 1,1,'last');
    figure
    subplot(3,1,1)
    plot(diff(vals(ind:end)))
    ylabel('Frame count diff')
    subplot(3,1,2)
    plot(diff(times(ind:end)))
    ylabel('Frame interval diff (us)')
    ylim([0 inf])
    subplot(3,1,3)
    plot(ncount)
    ylabel('Frame count per trial')

    [stimOns stimOffs] = photoFrameFinder_Sanworks(info.frame);
    if isfield(input, 'tGratingDirectionDeg')
        if isfield(input, 'nScansOn')
            nOff = input.nScansOff;
            nOn = input.nScansOn;
            nTrials = length(input.tGratingDirectionDeg);
            stimOnMat = double(nOff+1:nOff+nOn:(nOff+nOn)*nTrials);
            cStimOn = celleqel2mat_padded(input.cStimOn);
        elseif isfield(input, 'cTargetOn')
            cStimOn = celleqel2mat_padded(input.cTargetOn);
            stimOnMat = nan(size(cStimOn));
        end
    elseif isfield(input,'tStimOneGratingDirectionDeg')
        cStimOn = celleqel2mat_padded(input.cStimOneOn);
        stimOnMat = nan(size(cStimOn));
    else
        error('Need to define correct experiment')
    end
    length(cStimOn)
    length(stimOns)
    figure;
    subplot(3,1,1)
    if length(cStimOn) < length(stimOns)
        stimOns = stimOns(1:length(cStimOn));
        plot(cStimOn-stimOns)
        title('PD has more trials')
    elseif length(cStimOn) > length(stimOns)
        cStimOn = stimOns(1:length(stimOns));
        plot(cStimOn-stimOns)
        title('Mworks has more trials')
    else
        plot(cStimOn-stimOns)
    end
    ylabel('Mworks-photodiode Offset')
    subplot(3,1,2)
    plot(cStimOn-stimOnMat)
    ylabel('Counter-Matrix')
    subplot(3,1,3)
    if length(stimOffs) < length(stimOns)
        stimOns = stimOns(1:length(stimOffs));
    end
    plot(stimOffs-stimOns)
    ylabel('Stim on frames')
end