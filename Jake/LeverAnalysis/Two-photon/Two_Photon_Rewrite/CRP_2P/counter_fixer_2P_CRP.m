function [counterValues, counterTimesMs, shift_info] = counter_fixer_2P_CRP(counterValues, counterTimesMs, session_date, mouse_id)
%This function takes the frame times and values and corrects for any misalignment issues. 
%It overwrites the old vectors. Used in cleanBehav and cleanBehave_CRP

%test that number of elements match
assert(length(counterValues) == length(counterTimesMs));

%caculate some useful variables
ct_vals_diff = diff(counterValues);
ct_time_diff = diff(counterTimesMs);
IFI = mode(ct_time_diff);

%Identify and fill in skipped frames
skipped_frame_ind = find(ct_vals_diff==2); 
skipped_frame_ind = fliplr(skipped_frame_ind); %flip left right so that you will insert frames from the back end towards the front and indexing remains consistent during that process. 
for frame_num = skipped_frame_ind
    %check to see if there is a single missing frame time 
    if ct_time_diff(frame_num) > IFI*1.85 && ct_time_diff(frame_num) < IFI*2.15
        counterValues = [counterValues(1:frame_num), (counterValues(frame_num)+1), counterValues(frame_num+1:end)];
        counterTimesMs = [counterTimesMs(1:frame_num), (counterTimesMs(frame_num)+IFI), counterTimesMs(frame_num+1:end)];
    end
end

%redefine diff variables
ct_vals_diff = diff(counterValues);
ct_time_diff = diff(counterTimesMs);

%identify miss-timed frame
late_frames = find(ct_time_diff>IFI*1.4 & ct_time_diff<IFI*1.9);
if ~isempty(late_frames)
    for this_late_frame = late_frames
        if ct_time_diff(this_late_frame+1)<IFI*0.5 & ct_time_diff(this_late_frame+1)>IFI*0.1;
            interval_check = counterTimesMs(this_late_frame+2) - counterTimesMs(this_late_frame);
            assert(interval_check>IFI*1.85 &  interval_check<IFI*2.15 );  %verify that there is a supposed to be a frame here. Rought 2 IFI
            counterTimesMs(this_late_frame+1) = counterTimesMs(this_late_frame)+IFI;
            ct_vals_diff = diff(counterValues);
            ct_time_diff = diff(counterTimesMs);
        end
    end
end

%identify sequential missed frame followed by a missed timed frame
skip_and_late_frames = find(ct_time_diff>IFI*2.25 & ct_time_diff<IFI*2.5);
if ~isempty(skip_and_late_frames)
    for this_late_frame = skip_and_late_frames
        if ct_time_diff(this_late_frame+1)<IFI*0.75 & ct_time_diff(this_late_frame+1)>IFI*0.5  &  ct_vals_diff(this_late_frame)==2
            counterValues = [counterValues(1:this_late_frame), counterValues(this_late_frame)+1, counterValues(this_late_frame+1:end)];
            counterTimesMs = [counterTimesMs(1:this_late_frame), counterTimesMs(this_late_frame)+IFI, counterTimesMs(this_late_frame)+2*IFI, counterTimesMs(this_late_frame+2:end)];
            ct_vals_diff = diff(counterValues);
            ct_time_diff = diff(counterTimesMs);
        end
    end
end

%Write in an IF statement for 171227_img067 170428_img92
if strcmp(session_date, '171227') && mouse_id == 67
    string_match = 1;
elseif strcmp(session_date, '170428') && mouse_id == 992
    string_match = 1;
else 
    string_match = 0;
end
if string_match == 1
    bad_frames = intersect(find(ct_time_diff<IFI*0.85), find(ct_vals_diff>1));
    first_bad_frame = bad_frames(1)+1;
    if strcmp(session_date, '171227') && mouse_id == 67
        last_bad_frame = bad_frames(end);
    elseif strcmp(session_date, '170428') && mouse_id == 992
        last_bad_frame = bad_frames(end)+1;
    end
    
    %remove the bad frames
    counterValues(first_bad_frame:last_bad_frame) = [];
    counterTimesMs(first_bad_frame:last_bad_frame) = [];
    
    %find number of missing frames using the interval in counterTimesMs 
    splice_dur = counterTimesMs(first_bad_frame) - counterTimesMs(first_bad_frame-1);
    exact_div = double(splice_dur)/double(IFI); 
    round_div = double(round(splice_dur/IFI));
    assert(exact_div-round_div > -0.1 && exact_div-round_div < 0.1);
    splice_dur_frame = round(splice_dur/IFI)-1;
    splice_int_frames = [1:splice_dur_frame] + counterValues(first_bad_frame-1);
    splice_int_ms = [1:splice_dur_frame]*IFI + counterTimesMs(first_bad_frame-1);
    counterValues = [counterValues(1:first_bad_frame-1), splice_int_frames, counterValues(first_bad_frame:end)];
    counterTimesMs = [counterTimesMs(1:first_bad_frame-1), splice_int_ms, counterTimesMs(first_bad_frame:end)];
    
    %adjust frame nums of all subsequent frames so that diff(counterValues(first_bad_frame+splice_dur_frame-1), counterValues(first_bad_frame+splice_dur_frame) ) == 1
    frame_num_shift = counterValues(first_bad_frame+splice_dur_frame) - counterValues(first_bad_frame+splice_dur_frame-1) - 1;
    counterValues(first_bad_frame+splice_dur_frame:end) = counterValues(first_bad_frame+splice_dur_frame:end) - frame_num_shift;
    
    %save the shift_info so you can use it to fix the frame information in b_data  e.g. b_data.cLeverDown
    shift_info.frame_num_shift = frame_num_shift;
    shift_info.first_bad_frame = first_bad_frame;
    shift_info.last_bad_frame = first_bad_frame;
    shift_info.splice_int_frames = splice_int_frames; 
else 
    shift_info = [];
end

%look for the brief first pulse
if diff(counterTimesMs(1:2)) > 130
    counterTimesMs = counterTimesMs(2:end);
    counterValues = counterValues(2:end)-1;
end

%add in a portion of code which verifies;
% 1) minimal variance in counterTimesMs
% 2) counterValues is monotonicly increasing
ct_vals_diff = diff(counterValues);
ct_time_diff = diff(counterTimesMs); 
assert(max(ct_time_diff)<IFI*1.5); assert(min(ct_time_diff)>IFI*0.5);
assert(max(ct_vals_diff ==1)); assert(min(ct_vals_diff ==1));

counterTimesMs = double(counterTimesMs);
counterValues = double(counterValues); 

end
