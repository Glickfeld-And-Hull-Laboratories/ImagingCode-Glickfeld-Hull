function [wheel_speed] = wheelSpeedCalc(input,ticks,wheel)

%outputs wheel_speed in cm/s for each frame
%input is (1) input structure, (2) number of ticks on encoder (32, 64 or
%1000), (3) and type of wheel ('r' for red).

if strcmp(wheel,'red')
    diameter = 13.6;
    circ = pi.*diameter;
else
    error('wrong color')
end

cm_per_tick = circ./(ticks.*4);
ntrials = length(input.tGratingDirectionDeg);
counterTimes = cell2mat(input.counterTimesUs);
counterValues = cell2mat(input.counterValues);
wheelSpeedTimes = cell2mat(input.wheelSpeedTimesUs);
wheelSpeedValues = cell2mat(input.wheelSpeedValues)./...
    (1000./(round(mean(diff(input.wheelSpeedTimesUs{1}))./1000,-1))); % math to correct for ticks/s correction in labjack plugin

nframes = input.counterValues{end}(end);
wheel_speed = nan(1,nframes);
for iframe = 1:nframes-1
    fr_time_start = counterTimes(find(counterValues==iframe,1,'last'));
    fr_time_end = counterTimes(find(counterValues==iframe+1));
    if isempty(fr_time_start)
        fr_time_start = fr_time_end - mean(diff(counterTimes),2);
    end
    if isempty(fr_time_end)
        fr_time_end = fr_time_start + mean(diff(counterTimes),2);
    end
    ind = intersect(find(wheelSpeedTimes>=fr_time_start),find(wheelSpeedTimes<=fr_time_end));
    if length(ind)>0
        wheel_speed(:,iframe) = mean(wheelSpeedValues(ind),2);
    else
        wheel_speed(:,iframe) = 0;
    end
end
wheel_speed(:,nframes) = wheel_speed(:,nframes-1);
frame_dur = mean(diff(counterTimes),2)./1000000;
wheel_speed = (wheel_speed.*cm_per_tick)./frame_dur;