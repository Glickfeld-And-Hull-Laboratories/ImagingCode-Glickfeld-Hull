function [lever, frame, trial_outcome] = parse_behavior_for_HAD(input, first_frame, last_frame, frame_times, holdT_min)
 
% get the behavior data for hold and detect - code is very  very complicated.
% lever.state =[]; % 1 if lever pressed 0 if not, in ms resolution
% lever.press= []; % time of lever press
% lever.release = []; % time of lever release
% frame.counter = []; % frame counter, NaN if undefined. In 1 ms resolution
%
% trial_outcome.success_time =[]; % time of lever release (and sound)in success trials
% trial_outcome.early_time =[]; % time of lever release (and sound) in early release trials
% trial_outcome.late_time =[];  % time of error sound in late release trials
 
lever.state =[];           %these are the output variables of the function
lever.press= [];
lever.release = [];
lever.baseline_timesMs = [];
 
frame.counter = [];
frame.times = [];
frame.time_end =[];
frame.counter_by_time = [];
frame.f_frame_trial_num = [];
frame.l_frame_trial_num = [];

 
trial_outcome.success_time =[];
trial_outcome.early_time =[];
trial_outcome.late_time =[];
trial_outcome.ind_press_prerelease =[];
 
 
num_trials = length(input.tThisTrialStartTimeMs);                %gets vectors for # of trials and trial start time
trial_start =  round(double(cell2mat(input.tThisTrialStartTimeMs)));
 
% --- find state of lever before first trial
pre_trial_state = [];
for i=1:num_trials
    if(isempty(input.leverValues{i}))
        continue;
    end
    pre_trial_state = double(~(input.leverValues{i}(1)));        %pre_trial_state becomes important later on but I have no idea what it is doing here?????
    break;
end
 
if(isempty(pre_trial_state))
    warning('no lever press');
end
pre_trial_state = 0;                                          %wouldnt this negate anything done to pre_trial_state???
 
% -- the frame counter before the behavior starts is undefined
pre_frame_count = NaN;
pre_frame_time = trial_start(1);
% -- main loop  go over all trials and extract information
START_EXP_TIME = trial_start(1);
for i=1:num_trials-1
    t_begin = trial_start(i);                                       %values for trial start, end 
    t_end = trial_start(i+1);
    t_len = t_end - t_begin;                                        %trial length
    %vec = zeros(t_end - t_begin,1);
    
    % ------ get the vector of lever presses
    lever_time = ceil(double(input.leverTimesUs{i})/1000 - t_begin)+1;              %extracts lever times and values for the current trial
    lever_value= input.leverValues{i};
    
    bool_press = zeros(t_len,1) + pre_trial_state; % init with prev trial state       for each millisecond of the trial bool_press will record lever state
    pre_press_time = 1;
    for j=1:length(lever_time)                                    %1 to # of lever events in that trial
        
        if( j==1)
            if(lever_value(1) == pre_trial_state)
                warning(['check this perhaps a lever press/release was missed on trial ' num2str(i)]);
            end
        end
        
        bool_press(pre_press_time:lever_time(j)-1) = pre_trial_state;
        pre_press_time = lever_time(j);
        pre_trial_state =  double(lever_value(j));
        
        
        if(j==length(lever_time)) % last press
            bool_press(lever_time(j):end) = pre_trial_state;
        end
    end
    
    % ----- save a binary vector with state of the lever
    lever.state(end+1:end+length(bool_press)) = bool_press;
    % ----- save release and press times
    release = lever_value == 0;
    release_time =double(input.leverTimesUs{i}(release))/1000  -  START_EXP_TIME;
    lever.release(end+1:end+length(release_time)) = release_time';
   
    
    press = ~release;
    press_time =double(input.leverTimesUs{i}(press))/1000  -  START_EXP_TIME;
    lever.press(end+1:end+length(press_time)) = press_time';
    
%     lever.press = cat(1,  lever.press,  press_time');
%     lever.trial_press = cat(1, lever.trial_press, ones(size(press_time')).*i);
    
    %  ---------- calculate frame counter
    c_count = double(input.counterValues{i});       %counter values for this trial
    c_time =  double(input.counterTimesUs{i})/1000;    %counter times for this trial
    
    frame.times(end+1:end+length(c_time)) = c_time;
    frame.counter_by_time(end+1:end+length(c_count)) = c_count;     %you now have two vectors. One of counter times the other of the associated counter value
    c_time = round(c_time);
    if(length(c_count) < 2) % if less than 2 camera is off
        pre_frame_count = NaN;                                        %pre-frame_count is originally the MWtime of experiment start but it is altered within forloop
    end
    for j=1:length(c_count)
        d_frame = c_count(j) - pre_frame_count;                     %d_frame and time equal a specific counter value or time within this trial minus the frames or time before the trial
        d_time = c_time(j) - pre_frame_time;                               
        update_inx = (pre_frame_time:c_time(j)) -  START_EXP_TIME+1;          %update_inx is all the time points associated with a specific frame aligned to experiment start time
        
        if(d_frame >1)  %missed a frame, set NaN
            frame.counter(update_inx) = NaN;
        else
            frame.counter(update_inx) = pre_frame_count;
        end
        
        % update timer time
        pre_frame_count = c_count(j);
        pre_frame_time =  c_time(j);
    end
    
end

%This code looks at frame.counter while it is still in free floating
%MWtimes. It extracts the values needed to cut baseline_timesMs to remove
%all events without associated frames and align it to frame.counter after
%frame.counter has been cut to align with the camera start.
first_time_baseline=  find(frame.counter>= first_frame, 1, 'first')+START_EXP_TIME;
last_time_baseline=  find(frame.counter<= last_frame, 1, 'last')+START_EXP_TIME;
 

assert(length(frame.times) == length(frame.counter_by_time));
[uval, uinx] = unique(frame.counter_by_time, 'first');
frame.counter_by_time = uval;
frame.times = frame.times(uinx)-START_EXP_TIME;     
 
if(~isempty(frame_times))
    frame_times = round(frame_times - frame_times(1) +1);
    counter_by_frame_times(round(frame_times)) = 1:length(frame_times);
    for i=1:length(frame_times)-1
        inx = frame_times(i):(frame_times(i+1)-1);
        counter_by_frame_times(inx) = i;
    end
    
end
 
% find sucess/early release /late release times
%  trial_outcome.success_time =[]; % time of lever release (and sound)in success trials
%  trial_outcome.early_time =[]; % time of lever release (and sound)in early release trials
%  trial_outcome.late_time =[];  % time of error sound in late release trials
% trial_outcome.change_orientation = [];
 
% trial outcome input.trialOutcomeCell
hold_start = double(cell2mat(input.holdStartsMs)) - START_EXP_TIME;
hold_time = double(cell2mat(input.holdTimesMs));
 
react_time = double(cell2mat(input.reactTimesMs));
req_hold = double(cell2mat(input.tTotalReqHoldTimeMs));
rnd_hold = double(cell2mat(input.tRandReqHoldTimeMs));
 
has_reward = ~cellfun(@isempty, input.juiceTimesMsCell );
 
relase_time = hold_start + hold_time;
trial_outcome.success_time  = relase_time(has_reward);
trial_outcome.early_time = relase_time(hold_time<req_hold ); % do not use too early trial with orientation change
%trial_outcome.late_time = relase_time(~has_reward & hold_time>req_hold+double(input.tooFastTimeMs));
is_ignore =@(x)isequal(x, 'ignore');
trial_outcome.late_time  = relase_time(find(cellfun(is_ignore,input.trialOutcomeCell)));
trial_outcome.change_orientation = hold_start + req_hold;
trial_outcome.change_orientation(react_time < 0) = NaN; % if respond before change then no change in orientation

% This is a very very problematic code, the camera sends codes even before
% exqusition starts we detect that by a large gap between frames
if(~isempty(input.counterValues{1}))
    diff_count = diff(frame.times);
    %     change_inx = find(diff_count ==1 | isnan(diff_count));
    TOO_LARGE_VAL = 500;
    % find the first change that is larger than 500 ms
    too_large = find(diff_count > TOO_LARGE_VAL);
    if(~isempty(too_large))
        warning('realigning camera frames');
        assert(length(too_large) ==1); % very  strong assertion, just to be sure
        init_time = round(frame.times(too_large+1));
        frame.counter = frame.counter - frame.counter(init_time);
        frame.counter(1:init_time) = NaN;
        frame.times = frame.times(too_large+1:end);
        frame.counter_by_time = frame.counter_by_time(too_large+1:end);
        frame.counter_by_time = frame.counter_by_time - frame.counter_by_time(1)+1;
    end
end

% remove all events that are not between first and last frame, adjust indices
if(exist('first_frame' ,'var') && exist('last_frame', 'var'))
    
    first_time=  find(frame.counter>= first_frame, 1, 'first');
    last_time=  find(frame.counter<= last_frame, 1, 'last');
    frame.counter =frame.counter(first_time:last_time);
    lever.state = lever.state(first_time:last_time); % 1 if lever pressed 0 if not, in ms resolution
    
    inx = find(frame.counter_by_time>=first_frame & frame.counter_by_time<=last_frame);
    frame.counter_by_time = frame.counter_by_time(inx) - first_frame+1;
    frame.times = frame.times(inx) - first_time +1;
    
    
    lever.press = remove_events_and_adjust_inx(lever.press, first_time, last_time);
    lever.release = remove_events_and_adjust_inx(lever.release, first_time, last_time);
    trial_outcome.success_time = remove_events_and_adjust_inx(    trial_outcome.success_time, first_time, last_time);
    trial_outcome.early_time = remove_events_and_adjust_inx(trial_outcome.early_time, first_time, last_time);
    trial_outcome.late_time = remove_events_and_adjust_inx(trial_outcome.late_time, first_time, last_time);
    trial_outcome.change_orientation = remove_events_and_adjust_inx(trial_outcome.change_orientation, first_time, last_time);
    
    
    if(~isempty(frame_times))  
        first_time1=  find(counter_by_frame_times>= first_frame, 1, 'first');
        last_time1=  find(counter_by_frame_times<= last_frame, 1, 'last');
        
        counter_by_frame_times = counter_by_frame_times(first_time1:last_time1);
        min_l = min(length(counter_by_frame_times), length( frame.counter ));
        bad_inx = find(abs(counter_by_frame_times(1:min_l) - frame.counter(1:min_l)) >=3);
        if(length(bad_inx)/ min_l > 0.01)
            warning('something wrong with the camera synchronization!!!!!')
        end
         % -- if exist use the frame times as extracted from the camera 
        frame.counter = counter_by_frame_times' ;
    end
    frame.counter =frame.counter -  first_frame + 1;
end

if isempty(cell2mat(input.leverTimesUs(40:50)));
    nTrials = length(input.trialOutcomeCell);
    baseline_times = zeros(2,nTrials);   %baseline_times will be matrix containing the times (beg and end) of each window
    for iT = 1:nTrials-1;
        baseline_times(1,iT) = round(input.tStartTrialWaitForPressTimeMs{iT}-1200); 
        baseline_times(2,iT) = round(input.tStartTrialWaitForPressTimeMs{iT}-500);
    end
    
else
%LINDSEYS CODE to identify windows of at least 400ms between a lever release and
%subsequent press. Only takes windows during the iti. 
nTrials = length(input.trialOutcomeCell);
successIx = strcmp(input.trialOutcomeCell,'success');
failureIx = strcmp(input.trialOutcomeCell,'failure');
missIx = strcmp(input.trialOutcomeCell,'ignore');
baseline_times = zeros(2,nTrials);   %baseline_times will be matrix containing the times (beg and end) of each window
for iT = 1:nTrials-1
    ind_press = find(cell2mat(input.leverValues(iT))==1);      %finds locations of each lever press in the cell and flips L-R
    ind_release = find(cell2mat(input.leverValues(iT))==0);
    leverTimes = cell2mat(input.leverTimesUs(iT));   %time of each lever event
    trialStart = cell2mat(input.holdStartsMs(iT)).*1000;
    ind_prerelease = leverTimes<=trialStart;
    if sum(ind_prerelease,2) == 0;
        baseline_times(:,iT) = [NaN; NaN];
        trial_outcome.ind_press_prerelease(iT) = NaN;
    elseif isempty(ind_prerelease);
        baseline_times(:,iT) = [NaN; NaN];
        trial_outcome.ind_press_prerelease(iT) = NaN;
    elseif sum(ind_prerelease,2) == 1
        trial_outcome.ind_press_prerelease(iT) = 1;
        if leverTimes(1)- (cell2mat(input.tThisTrialStartTimeMs(iT))*1000) > holdT_min + 600000  %and if this lever press occurs a certain amount of time after trial start...
            baseline_times(1,iT) = leverTimes(1) - holdT_min - 300000; 
            baseline_times(2,iT) = leverTimes(1) - 300000;
            continue
        elseif baseline_times(:,iT) == 0
            baseline_times(:,iT) = [NaN; NaN];
        end
    else
        tag_prerelease = fliplr(find(ind_prerelease));
        trial_outcome.ind_press_prerelease(iT) = tag_prerelease(1);
        for i = 1:length(tag_prerelease)
            ip = tag_prerelease(i);
            if sum(ismember(ind_press,ip),2)
                if ip == 1
                    if leverTimes(1)- (cell2mat(input.tThisTrialStartTimeMs(iT))*1000) > holdT_min + 600000  %and if this lever press occurs a certain amount of time after trial start...
                        baseline_times(1,iT) = leverTimes(1) - holdT_min - 300000; 
                        baseline_times(2,iT) = leverTimes(1) - 300000;
                        break
                    elseif baseline_times(:,iT) == 0
                        baseline_times(:,iT) = [NaN; NaN];
                    end
                elseif ip>1    
                	if leverTimes(ip)- leverTimes(ip-1) > holdT_min + 600000
                       baseline_times(1,iT) = leverTimes(ip) - holdT_min - 300000; 
                       baseline_times(2,iT) = leverTimes(ip) - 300000;
                       break
                    elseif baseline_times(:,iT) == 0
                        baseline_times(:,iT) = [NaN; NaN];
                    end
                end
            end
        end
    end
end
end
ex_trials = length(find(isnan(baseline_times(1,:))));
disp(['---' num2str(ex_trials) ' trials excluded by pre-press f criteria ---'])
%get rid of NaNs in baseline_times
baseline_times2 = [];
% for i = 1:length(baseline_times);
%     if isnan(baseline_times(1,i))==0 & isnan(baseline_times(2,i))==0
%         baseline_times2 =[baseline_times2, baseline_times(:,i)];
%     end
% end
%baseline_times = baseline_times2;  
if isempty(cell2mat(input.leverTimesUs(40:50)));
    baseline_timesMs = baseline_times;
else
    baseline_timesMs = round(baseline_times/1000);
end

%remove events from baseline_times which do not have associated frames. JH
baseline_timesMs_cut = [];
for t = 1:length(baseline_timesMs)
    if baseline_timesMs(1,t)<first_time_baseline;
        baseline_timesMs(:,t)=NaN;
    elseif baseline_timesMs(2,t)>last_time_baseline;
        baseline_timesMs(:,t)=NaN;
    end
end
disp( [num2str(sum(~isnan(baseline_timesMs(1,:)))) ' usable baseline windows within frame times']);

%convert 1st and last trials baseline_times to NaNs  
for i = 1:length(input.counterValues);
    if isempty(cell2mat(input.counterValues(i))) ==0;
    imat = cell2mat(input.counterValues(i));
    if imat(end) >= last_frame
        break
    end
    end
end
l_frame_trial_num = i;
for i = 1:length(input.counterValues);
       if isempty(cell2mat(input.counterValues(i))) ==0;
           break
       end
end
f_frame_trial_num = i;
baseline_timesMs(:,f_frame_trial_num)=NaN;            %the first and last trials with camera pulses will be incompletely imaged. Therefore they are unusable baseline_times
baseline_timesMs(:,l_frame_trial_num)=NaN;

%subtract the MWtime of camera start in order to align baseline_timesMS to
%frame.counter
baseline_timesMs = baseline_timesMs - first_time_baseline;
%save it to the struct bc this baby is done
lever.baseline_timesMs = baseline_timesMs;  
frame.f_frame_trial_num = f_frame_trial_num;
frame.l_frame_trial_num = l_frame_trial_num;


return;
 
function res = remove_events_and_adjust_inx(vec, first_time, last_time)
res = vec(vec >=first_time & vec <= last_time)  - first_time;