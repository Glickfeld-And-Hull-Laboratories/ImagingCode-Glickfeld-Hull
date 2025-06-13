function [frameErr_mult, frameErr_skip, frameErr_add] = counterValCheck(input)

%frameErr_mult identifies frames that were added (and counted) in Mworks
%frameErr_skip identifies frames that were skipped by Mworks
%frameErr_add identifies frames that were added but not counted

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
    vals = vals(ind:end);
    times = times(ind:end);

    frame_int = mean(diff(times(1:1000)));
    dtimes = diff(times);
    dvals = diff(vals);
    ind_skip_int = find(dtimes>1.75*frame_int);
    ind_skip_fr = find(dvals==2);
    ind_skip_multfr = find(dvals>2);
    ind_mult_fr = find(dtimes<.75*frame_int);

    frameErr_skip = setdiff(ind_skip_int,ind_skip_fr);
    frameErr_mult = setdiff(ind_skip_fr,[ind_mult_fr ind_skip_int]);
    frameErr_add = [ind_skip_multfr intersect(ind_skip_fr,ind_mult_fr)];
    
    if ~isempty(ind_skip_multfr)
        for i = 1:length(ind_skip_multfr)
            if find(abs(ind_mult_fr-ind_skip_multfr(i))<2)
                ind_mult_fr(find(abs(ind_mult_fr-ind_skip_multfr(i))<2)) = [];
            end
        end
    end

    for i = 1:length(ind_mult_fr)
        if isempty(find(ind_skip_fr == ind_mult_fr(i)))
            if dtimes(ind_mult_fr(i)+1)<1.25*frame_int & dtimes(ind_mult_fr(i)-1)<1.25*frame_int
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

    figure;
    subplot(2,1,1)
    plot(dtimes)
    ylim([0 inf])
    if length(frameErr_skip)>0
        xline(frameErr_skip,'--r')
    end
    if length(frameErr_mult)>0
        xline(frameErr_mult,'--b')
    end
    if length(frameErr_add)>0
        xline(frameErr_add,'--y')
    end
    subplot(2,1,2)
    plot(dvals)
    if length(frameErr_skip)>0
        xline(frameErr_skip,'--r')
    end
    if length(frameErr_mult)>0
        xline(frameErr_mult,'--b')
    end
    if length(frameErr_add)>0
        xline(frameErr_add,'--y')
    end
    frameErr_skip = vals(frameErr_skip)+1;
    frameErr_mult = vals(frameErr_mult)+1;
    frameErr_add = vals(frameErr_add)+1;
end