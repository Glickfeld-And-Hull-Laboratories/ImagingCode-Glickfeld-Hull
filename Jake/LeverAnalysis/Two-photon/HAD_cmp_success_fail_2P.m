%clear;
clear
WRITE_VEDIO = 0;
BIN_SIZE = 1;   % if the ROI is large use large bin size (10)

%expt info
date = '150703';
run = '_000_000';
mouse = 'img24';
holdT_min = 500000;

%output directory
out_base = 'Z:\home\lindsey\Analysis\2P\Jake';
run_name = [date '_' mouse '_run' run(length(run)-2:end)];
out_path = fullfile(out_base,run_name);
dest =  fullfile(out_path,run_name);

image_dest  = [dest '_ROI.tif'];
frame_info_dest = [dest '_frame_times.mat'];

% --- load frame times and behavior data
load(frame_info_dest);
b_data.input = input; clear input;
ftimes.frame_times = frame_times; clear frame_times;
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
[lever, frame_info, trial_outcome] = parse_behavior_for_HAD(b_data.input, ...
    f_frame, l_frame, ftimes.frame_times, holdT_min);

data_dest = [dest '_parse_behavior.mat'];
save(data_dest, 'lever', 'frame_info', 'trial_outcome', 'Sampeling_rate', 'holdT_min')

%Obtain a df/f movie using Lindsey's baseline_times
img = readtiff(image_dest);
sz = size(img(:,:,1));
img = reshape(img,[sz(1)*sz(2) size(img,3)]);
startT = round(b_data.input.counterTimesUs{1}(1)./1000);
img_dfoverf = zeros(size(img));    %this could be problematic due to the frame skipping issue
first_baseline = find(~isnan(lever.baseline_timesMs(1,:)),1, 'first');    %find the first trial / baseline_timesMs window that is not NaN
F_range = [];
for iT=2:length(lever.baseline_timesMs)-1;    %this could be problematic due to unremoved NaNs
    if ~isnan(lever.baseline_timesMs(1,iT));
        F_range = frame_info.counter(lever.baseline_timesMs(1,iT)):frame_info.counter(lever.baseline_timesMs(2,iT));
    elseif isempty(F_range)
        F_range = frame_info.counter(lever.baseline_timesMs(1,first_baseline)):frame_info.counter(lever.baseline_timesMs(2,first_baseline));
    end
    F_avg= mean(img(:,F_range),3);
    t_range = frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT))-startT):frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT+1))-startT);
    t_df = bsxfun(@minus, double(img(:,t_range)), F_avg);
    t_dfoverf = bsxfun(@rdivide, t_df, F_avg);
    img_dfoverf(:,t_range) = t_dfoverf;
end 

% ---- do simple movie analysis
func = @median;
% func = @mean;
%func = @std;
pre_frames = 5;
post_frames = 10;
post_press_frames = 10;   
rm_baseline_plot = 0; % 1 for removing baseline when plotiing

ts = (-pre_frames:post_frames)*1000/round(double(Sampeling_rate));
tot_frame = pre_frames + post_frames+1;

use_ev_success = trial_outcome.success_time;
if strcmp(b_data.input.trialOutcomeCell{1}, 'success')
    use_ev_success(1) = [];
elseif strcmp(b_data.input.trialOutcomeCell{end}, 'success')
    use_ev_success(end) = [];
end
%----- uncomment to use only events w/o lever press after release
%     use_ev_success = remove_events_by_lever_state(use_ev_success,  ...
%         lever.state, 10,ceil(post_frames*1000/Sampeling_rate), 0);
%------ uncomment to use event only w/o lever press before release time
%     use_ev_success = remove_events_by_lever_state(use_ev_success,  ...
%         lever.state, -ceil(pre_frames*1000/Sampeling_rate),0, 1);
%
success_movie = trigger_movie_by_event(img_dfoverf, frame_info, ...
    use_ev_success, pre_frames, post_frames);
avg_success_move = squeeze(func(success_movie,1));
fig1 = figure; plot_movie(avg_success_move,sz,rm_baseline_plot, ts);
title('Hits');
figure;
for i = 1:size(avg_success_move,2)
    subplot(4,4,i)
    imagesc(reshape(avg_success_move(:,i),[sz(1) sz(2)]))
    clim([0 0.3])
    title(num2str(i-pre_frames-1))
    set(gca,'XTickLabel','','YTickLabel','')
end
suptitle('Hits')

use_ev_fail = trial_outcome.early_time;
if strcmp(b_data.input.trialOutcomeCell{1}, 'failure')
    use_ev_fail(1) = [];
elseif strcmp(b_data.input.trialOutcomeCell{end}, 'failure')
    use_ev_fail(end) = [];
end
%----- uncomment to use only events w/o lever press after release
%     use_ev_fail = remove_events_by_lever_state(use_ev_fail,  ...
%         lever.state, 10,ceil(post_frames*1000/Sampeling_rate), 0);
%------ uncomment to use event only w/o lever press before release time
%     use_ev_fail = remove_events_by_lever_state(use_ev_fail,  ...
%         lever.state, -ceil(pre_frames*1000/Sampeling_rate),0, 1);
%

% -----trigger movie by early release
fail_movie = trigger_movie_by_event(img_dfoverf, frame_info, ...
    use_ev_fail, pre_frames, post_frames);
avg_fail_move = squeeze(func(fail_movie,1));
fig2 = figure; plot_movie(avg_fail_move,sz,rm_baseline_plot, ts);
title('False Alarms');

figure;
for i = 1:size(avg_fail_move,2)
    subplot(4,4,i)
    imagesc(reshape(avg_fail_move(:,i),[sz(1) sz(2)]))
    clim([0 0.3])
    title(num2str(i-pre_frames-1))
    set(gca,'XTickLabel','','YTickLabel','')
end
suptitle('False Alarms')

figure;
subplot(2,1,1) 
imagesc(reshape(mean(avg_success_move(:,7:10),2),[sz(1) sz(2)]))
clim([-.05 0.3])
set(gca,'XTickLabel','','YTickLabel','')
title('Hits')
subplot(2,1,2) 
imagesc(reshape(mean(avg_fail_move(:,7:10),2),[sz(1) sz(2)]))
clim([-.05 0.3])
set(gca,'XTickLabel','','YTickLabel','')
title('False Alarms')

% --- make color scale the same
min_val = min([avg_fail_move(:);avg_success_move(:)]);
max_val = max([avg_fail_move(:);avg_success_move(:)]);
figure(fig1); caxis([min_val max_val]);
figure(fig2); caxis([min_val max_val]);

% ----- success  - early!!
diff_succ_fail = avg_success_move - avg_fail_move;
fig3 = figure; plot_movie(diff_succ_fail,sz,rm_baseline_plot, ts);    %JH ALTERED   Added "fig3"
title('Hits  - False Alarms');
figure(fig3); caxis([min_val max_val]);                             %JH ADDED this entire line to make axes constant


% ---- Trigger movie off lever press IF subsequent hold time is greater
% than 500ms
post_frames = post_press_frames;
if lever.release(1) < lever.press(1)
    holdDuration = NaN(1, size((lever.release),2)-1);
    for i = 1:(length(holdDuration))
        holdDuration(i) = lever.release(i+1)-lever.press(i);
    end
else
    holdDuration = NaN(1, size((lever.release),2));
    for i = 1:(length(holdDuration))
        holdDuration(i) = lever.release(i)-lever.press(i);
    end
end

longHoldInd = zeros(size(holdDuration));
for i = 1:size(holdDuration, 2);
    if holdDuration(i) > 500
        longHoldInd(i) = 1;
    end
end

use_ev_press = NaN(1, sum(longHoldInd));
aa = 1;
for i = 1:length(longHoldInd)
    if longHoldInd(i) == 1;
        use_ev_press(aa) = lever.press(i);
        aa=aa+1;
    end
end
use_ev_press = round(use_ev_press);
%----- uncomment to use only events w/o lever press after release
%     use_ev_success = remove_events_by_lever_state(use_ev_success,  ...
%         lever.state, 10,ceil(post_frames*1000/Sampeling_rate), 0);
%------ uncomment to use event only w/o lever press before release time
%     use_ev_success = remove_events_by_lever_state(use_ev_success,  ...
%         lever.state, -ceil(pre_frames*1000/Sampeling_rate),0, 1);
%
press_movie = trigger_movie_by_event(img, frame_info, ...
    use_ev_press, pre_frames, post_frames);
avg_press_move = squeeze(func(press_movie,1));
fig4 = figure; plot_movie(avg_press_move,sz,rm_baseline_plot, ts);
title('lever press');
figure(fig4);  caxis([min_val max_val]);     %JH ADDED this line to make axes constant

%PLOT TRIAL BY TRIAL VARIABILITY
%Hits
subset1 = round(linspace(1,size(success_movie,1),10));
success_subset = [];
success_subset = success_movie(subset1,:,:);
fig5 = figure; plot_movie(success_subset,sz,rm_baseline_plot, ts);
title('Hits');
clim([-5 12])   % NEED a smart way to do scaling

%False Alarms   
subset2 = round(linspace(1,size(fail_movie,1),10));
fail_subset = [];
fail_subset = success_movie(30:39,:,:);
fig6 = figure; plot_movie(fail_subset,sz,rm_baseline_plot, ts);
title('False Alarms');
clim([-12 30])

%Lever Presses
%Hits-FAsNot sure this will be informative. Also will require a
%lot of tweeking
%     subset3 = round(linspace(1,size(_movie,1),10));
%     _subset = [];
%     _subset = success_movie(subset3,:,:);
%     fig5 = figure; plot_movie(_subset,sz,rm_baseline_plot, ts);
%     title('Hits - False Alarms');
%     clim([-5 12])     

