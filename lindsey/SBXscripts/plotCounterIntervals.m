function [times vals] = plotCounterIntervals(input);
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
end