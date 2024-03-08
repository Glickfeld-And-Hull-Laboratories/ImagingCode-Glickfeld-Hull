function [wheel_speed] = wheelSpeedCalc(input,ticks,wheel)

%outputs wheel_speed in cm/s for each frame
%input is (1) input structure, (2) number of ticks on encoder (32, 64 or
%1000), (3) and type of wheel ('r' for red).

if strcmp(wheel,'red')
    diameter = 13.6;
    circ = pi.*diameter;
elseif strcmp(wheel,'purple')
    diameter = 14.6;
    circ = pi.*diameter;
elseif strcmp (wheel, 'orange')
  diameter = 13.3;
  circ = pi.*diameter;
else
    error('wrong color')
end


cm_per_tick = circ./(ticks.*4);
if isfield(input,'tGratingDirectionDeg')
    ntrials = length(input.tGratingDirectionDeg);
elseif isfield(input,'tStimOneGratingDirectionDeg')
    ntrials = length(input.tStimOneGratingDirectionDeg);
elseif isfield(input,'tTestStimGratingDirectionDeg')
    ntrials = length(input.tTestStimGratingDirectionDeg);
elseif isfield(input,'cStimOneOn')
    ntrials = length(input.cStimOneOn);
else
    error('how many trials?')
end
counterTimes = cell2mat(input.counterTimesUs);
counterValues = cell2mat(input.counterValues);

%solves issue of extra events from multiple runs
diff_ind = find(diff(diff(counterValues))~=0);
counterValues(diff_ind) = [];
counterTimes(diff_ind) = [];
diff_ind = find(diff(counterValues)<1);
counterValues(diff_ind) = [];
counterTimes(diff_ind) = [];

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
    if isempty(fr_time_end) & isempty(fr_time_start)
        wheel_speed(:,iframe) = NaN;
    else
        ind = intersect(find(wheelSpeedTimes>=fr_time_start),find(wheelSpeedTimes<=fr_time_end));
        if length(ind)>0
            wheel_speed(:,iframe) = mean(wheelSpeedValues(ind),2);
        else
            wheel_speed(:,iframe) = 0;
        end
    end
end
wheel_speed(:,nframes) = wheel_speed(:,nframes-1);
frame_dur = mean(diff(counterTimes),2)./1000000;
wheel_speed = (wheel_speed.*cm_per_tick)./frame_dur;