function [trigger_movie, remove_event_idx, hold_dur, use_times, trigger_licks, trigger_licks_10ms, lick] = trigger_movie_by_event_2P(movie, frame_info, ...
    event_times, pre_frames, post_frames, licking_data, hold_dur, do_lickAna, do_alignLick)
% ---- cut the movie around event times
% movie - the movie in 2-dimantional matrix
% frame_info: output of parse behavior syncs frame number to time  in behavior

% --- use only events that are in the times of the movie, leave space for last event
last_time=  find(frame_info.counter<size(movie,2)-post_frames, 1, 'last');
use_event_times = event_times(event_times < last_time);
first_time = find(frame_info.counter>pre_frames+1, 1, 'first');
% first_time = 758175; %for img38_160320 remove first 25 trials
use_event_times = use_event_times(use_event_times >first_time);
if ~isempty(hold_dur)
    hold_dur = hold_dur(event_times < last_time & event_times > first_time);
    
end
trigger_movie =nan([length(use_event_times) size(movie,1), pre_frames + post_frames+1]);
use_times =[]; lick.bout_onset = []; lick.cue_onset = [];
lick.reward_onset = []; lick.reward_late_onset = []; lick.lickTrial = []; lick.lickTrial_1000 = [];

% if there is licking_data
if ~isempty(licking_data)
    lickTimes = licking_data.lickTimes;
    licksByFrame = licking_data.licksByFrame;
    trigger_licks = nan([length(use_event_times), pre_frames + post_frames+1]);
    trigger_licks_10ms = nan([length(use_event_times), pre_frames*10 + post_frames*10+1]);
else
    trigger_licks = [];
    trigger_licks_10ms = [];
end

for i=1:length(use_event_times)
    c_time = use_event_times(i);
    % inefficient - but simple in the future use indexes
    % find the corresponding frame
    frame_no = frame_info.counter(c_time);
    if(isnan(frame_no))
        continue; % discard
    end
    use_times(end+1) = c_time;
    
    if ~do_alignLick
        trigger_movie(i,:,:) = movie(:, frame_no-pre_frames:frame_no+post_frames);
        
    else
        
        temp_lick = licksByFrame(frame_no:frame_no+post_frames);
        if ~isempty( find(temp_lick, 1, 'first') )
            frame_no = find(temp_lick, 1, 'first') + frame_no - 1;
        end
            
        trigger_movie(i,:,:) = movie(:, frame_no-pre_frames:frame_no+post_frames);
        
    end
    if ~isempty(licking_data)
        trigger_licks(i,:) = licksByFrame(frame_no-pre_frames:frame_no+post_frames);
        if do_lickAna == 1
            % find lick bout onset
            if sum(trigger_licks(i, pre_frames+1:pre_frames + floor(500/frame_info.ifi))) <= 1 % no more than 1 lick during 500ms (15 frames) before reward
                lick_reward = trigger_licks(i, pre_frames + floor(500/frame_info.ifi)+1 :end);
                lick_idx = find(lick_reward);
                if length(lick_idx) >= 4
                    inter_lick = diff(lick_idx);
                    %i
                    for li = 1:length(inter_lick)
                        % li
                        if inter_lick(li) <= 6 && length(inter_lick(li:end)) >= 3 && mean(inter_lick(li: li+2 )) <=6 % each lick is less than 6 frames apart
                            % and at least three inter licks in bewteen
                            lick.bout_onset = [lick.bout_onset lick_idx(li) + pre_frames + floor(500/frame_info.ifi)]; % match index with trigger_licks
                            break;
                        end
                    end
                end
            end
            
            % find lick onset 0-350ms after cue
            lick_cue = trigger_licks(i, pre_frames + 1: pre_frames + floor(350/frame_info.ifi) + 1);
            if sum(lick_cue) > 0
                lick.cue_onset = [lick.cue_onset find(lick_cue,1,'first') + pre_frames];
            end
            
            % find lick onset 150ms around reward
            lick_reward = trigger_licks(i, pre_frames + floor(350/frame_info.ifi) + 2: pre_frames + floor(650/frame_info.ifi) + 2);
            if sum(lick_reward) > 0
                lick.reward_onset = [lick.reward_onset find(lick_reward,1,'first') + pre_frames + floor(350/frame_info.ifi) + 1];
            end
            
            % find lick onset 150ms after reward
            lick_reward_late = trigger_licks(i, pre_frames + floor(650/frame_info.ifi) + 3: end);
            if sum(lick_reward_late) > 0
                lick.reward_late_onset = [lick.reward_late_onset find(lick_reward_late,1,'first') + pre_frames + floor(650/frame_info.ifi) + 2];
            end
            
            % find trial without lick bout
            if sum(trigger_licks(i, pre_frames - floor(100/frame_info.ifi):pre_frames + floor(500/frame_info.ifi))) > 1 % no more than 1 lick 100ms before
                % cue and 100ms before reward
                lick_cuereward = trigger_licks(i, pre_frames - floor(100/frame_info.ifi):pre_frames + floor(500/frame_info.ifi));
                lick_idx = find(lick_cuereward);
                if length(lick_idx) >= 4
                    inter_lick = diff(lick_idx);
                    %i
                    for li = 1:length(inter_lick)
                        % li
                        if inter_lick(li) <= 6 && length(inter_lick(li:end)) >= 3 && mean(inter_lick(li: li+2 )) <=6 % each lick is less than 6 frames apart
                            % and at least three inter licks in bewteen
                            lick.lickTrial = [lick.lickTrial 1];
                            break;
                        end
                        
                        if length(inter_lick(li:end)) < 3  % if less than 3 licks in the end, count as no lick trial
                            lick.lickTrial = [lick.lickTrial 0];
                            break;
                        end
                        
                    end
                else
                    lick.lickTrial = [lick.lickTrial 0];
                end
            else
                lick.lickTrial = [lick.lickTrial 0];
            end
            
            % find trial without lick bout for 1000ms case
            if sum(trigger_licks(i, pre_frames - floor(100/frame_info.ifi):pre_frames + floor(1600/frame_info.ifi))) > 1 % no more than 1 lick 100ms before
                % cue and 100ms before reward
                lick_cuereward = trigger_licks(i, pre_frames - floor(100/frame_info.ifi):pre_frames + floor(1600/frame_info.ifi));
                lick_idx = find(lick_cuereward);
                if length(lick_idx) >= 4
                    inter_lick = diff(lick_idx);
                    %i
                    for li = 1:length(inter_lick)
                        % li
                        if inter_lick(li) <= 6 && length(inter_lick(li:end)) >= 3 && mean(inter_lick(li: li+2 )) <=6 % each lick is less than 6 frames apart
                            % and at least three inter licks in bewteen
                            lick.lickTrial_1000 = [lick.lickTrial_1000 1];
                            break;
                        end
                        
                        if length(inter_lick(li:end)) < 3  % if less than 3 licks in the end, count as no lick trial
                            lick.lickTrial_1000 = [lick.lickTrial_1000 0];
                            break;
                        end
                        
                    end
                else
                    lick.lickTrial_1000 = [lick.lickTrial_1000 0];
                end
            else
                lick.lickTrial_1000 = [lick.lickTrial_1000 0];
            end
            
            %gather licks by 10ms bins
            frame_time = frame_info.times(frame_no);
            time_range = (frame_time-pre_frames*10):10:(frame_time+post_frames*10);
            licks_this_TC = [];
            for kk = 2:length(time_range); %find all the # of licks which occur in each time bin.
                licks_this_bin = lickTimes(find(lickTimes>=time_range(kk-1) & lickTimes<time_range(kk)));
                licks_this_TC  = [licks_this_TC, licks_this_bin];
            end
            trigger_licks_10ms(i,:);
        end
    end
    
end
remove_event_idx = find(~ismember(event_times, use_event_times));
%Attempted to access frame_info.counter(629.805); index must be a positive integer or logical.

%Error in trigger_movie_by_event (line 18)
%   frame_no = frame_info.counter(c_time);

%Error in HAD_cmp_success_fail (line 159)
%   press_movie = trigger_movie_by_event(img, frame_info, ...