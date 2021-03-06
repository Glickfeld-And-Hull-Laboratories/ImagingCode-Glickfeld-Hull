% triggers time-courses off of event times
% 1. finds frame and lever times based on events and outcomes
% 2. obtain df/f timecourse
% 3. create event triggered movies



%load frame and lever info
frame_info_dest = [dest '_frame_times.mat'];
load(frame_info_dest);
ftimes.frame_times = frame_times; clear frame_times;
b_data.input = input; clear input;
load([dest '_ROI_TCs.mat']);
%% 1. find frame and lever times
ifi = (ftimes.frame_times(end)-ftimes.frame_times(1))/length(ftimes.frame_times);
Sampeling_rate = 1000/ifi;

if(~exist('first_frame', 'var'))
    f_frame =1;
else
    f_frame = first_frame;
end

if(~exist('last_frame', 'var'))
    l_frame =length(ftimes.frame_times);
else
    l_frame = last_frame;
end
% ---- parse behavior
holdT_min  = 500000;    
[lever, frame_info, trial_outcome] = parse_behavior_for_HAD(b_data.input, ...
    f_frame, l_frame, ftimes.frame_times, holdT_min);

data_dest = [dest '_parse_behavior.mat'];
save(data_dest, 'lever', 'frame_info', 'trial_outcome', 'Sampeling_rate', 'holdT_min', 'ifi')

%% 2. Obtain a df/f TC from baseline times
data_tc = data_tc';
startT = round(b_data.input.counterTimesUs{1}(1)./1000);
tc_dfoverf = zeros(size(data_tc));    %this could be problematic due to the frame skipping issue
first_baseline = find(~isnan(lever.baseline_timesMs(1,:)),1, 'first');    %find the first trial / baseline_timesMs window that is not NaN
F_range = [];
for iT=2:length(lever.baseline_timesMs)-1;    %this could be problematic due to unremoved NaNs
    if ~isnan(lever.baseline_timesMs(1,iT));
        F_range = frame_info.counter(lever.baseline_timesMs(1,iT)):frame_info.counter(lever.baseline_timesMs(2,iT));
    elseif isempty(F_range)
        F_range = frame_info.counter(lever.baseline_timesMs(1,first_baseline)):frame_info.counter(lever.baseline_timesMs(2,first_baseline));
    end
    F_avg= mean(data_tc(:,F_range),2);
    if (cell2mat(b_data.input.tThisTrialStartTimeMs(iT+1))-startT)>length(frame_info.counter)
        t_range = frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT))-startT):frame_info.counter(end);
    else
        t_range = frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT))-startT):frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT+1))-startT);
    end
    t_df = bsxfun(@minus, double(data_tc(:,t_range)), F_avg);
    t_dfoverf = bsxfun(@rdivide, t_df, F_avg);
    tc_dfoverf(:,t_range) = t_dfoverf;
end
data_dest2 = [dest '_dFOverF_TC.mat'];
save(data_dest2, 'tc_dfoverf')

%% 3. create event triggered movies

func = @mean;
pre_release_ms = 500;
post_release_ms = 500;
pre_release_frames = round(pre_release_ms./double(ifi));
post_release_frames = round(post_release_ms./double(ifi));

%successes
use_ev_success = trial_outcome.success_time;
if strcmp(b_data.input.trialOutcomeCell{1}, 'success')  %removing first and last trials from consideration if they were successful
    use_ev_success(1) = [];
elseif strcmp(b_data.input.trialOutcomeCell{end}, 'success')
    use_ev_success(end) = [];
end
success_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_ev_success, pre_release_frames, post_release_frames);

%failures
use_ev_fail = trial_outcome.early_time;
if strcmp(b_data.input.trialOutcomeCell{1}, 'failure')   %removing first and last trials from consideration if they were failures
    use_ev_fail(1) = [];
elseif strcmp(b_data.input.trialOutcomeCell{end}, 'failure')
    use_ev_fail(end) = [];
end
fail_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_ev_fail, pre_release_frames, post_release_frames);

%tooFast correct
use_ev_tooFastCorrect = trial_outcome.tooFastCorrects;
if strcmp(b_data.input.trialOutcomeCell{1}, 'success')  %removing first and last trials from consideration if they were successful (tooFastCorrects are coded as successes here)
    use_ev_tooFastCorrect(1) = [];
elseif strcmp(b_data.input.trialOutcomeCell{end}, 'success')
    use_ev_tooFastCorrect(end) = [];
end
tooFastCorrect_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_ev_tooFastCorrect, pre_release_frames, post_release_frames);

%fidgets
use_ev_fidget = trial_outcome.fidget;
if strcmp(b_data.input.trialOutcomeCell{1}, 'failure')   %removing first and last trials from consideration if they were failures
    use_ev_fidget(1) = [];
elseif strcmp(b_data.input.trialOutcomeCell{end}, 'failure')
    use_ev_fidget(end) = [];
end
fidget_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_ev_fidget, pre_release_frames, post_release_frames);

save([dest_sub '_release_movies.mat'],'fail_movie','success_movie', 'tooFastCorrect_movie', 'use_ev_fidget','pre_release_frames','post_release_frames','ifi');

% ---- Trigger movie off all lever presses at trial start
pre_press_ms = 500;
post_press_ms = 500;
pre_press_frames = round(pre_press_ms./double(ifi));
post_press_frames = round(post_press_ms./double(ifi));

pressTime = NaN(1,length(lever.baseline_timesMs));
releaseTime = NaN(1,length(lever.baseline_timesMs));
for iT = 2:length(lever.baseline_timesMs)-1
    if ~isempty(b_data.input.leverTimesUs{iT})
        leverTimes = round((cell2mat(b_data.input.leverTimesUs(iT))-b_data.input.counterTimesUs{1}(1))./1000);
    end
    if ~isnan(trial_outcome.ind_press_prerelease(iT)) & size(b_data.input.leverTimesUs{iT},2)>1
        pressTime(1,iT) = leverTimes(trial_outcome.ind_press_prerelease(iT));
        releaseTime(1,iT) = leverTimes(trial_outcome.ind_press_prerelease(iT)+1);
    end
end

use_ev_press = pressTime;
use_ev_press(isnan(use_ev_press)) = [];

press_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_ev_press, pre_press_frames, post_press_frames);

%break up presses by hold time
holdTime = releaseTime-pressTime;
ind_200 = find(holdTime<200);
ind_500 = find(holdTime<500);
ind_both = find(ismember(ind_500,ind_200));
ind_500(ind_both) = [];
ind_long = find(holdTime>=500);
ind_both = find(ismember(ind_long,ind_500));
ind_long(ind_both) = [];

use_200_press = pressTime(ind_200);
use_500_press = pressTime(ind_500);
use_long_press = pressTime(ind_long);

press_200_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_200_press, pre_press_frames, post_press_frames);
press_500_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_500_press, pre_press_frames, post_press_frames);
press_long_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_long_press, pre_press_frames, post_press_frames);

save([dest_sub '_press_movies.mat'],'press_200_movie','press_500_movie','press_long_movie','press_movie','pre_press_frames', 'post_press_frames', 'ifi');

