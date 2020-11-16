%% Modified version of the anaylysis overview for analyzing the licking behavior of cue-reward trials. SJ
% licks by trial: trial*ms, binary matrix, licking time points = 1
% licks_by_trial_rewarded = licks_by_trial;
% RT_this_session: 1*ntrials vector, reaction time for each trial. Only trials with a reaction time longer than 200ms and shorter than 1200ms are included.
% full_trial_licks_rewarded_bin: bin licking every 100ms, lick rate in Hz for every 100ms. Looking at a total of 24 seconds

% since juice==0 in the first trial, this script gets rid of the first trial

%% SECTION ONE - assign pathnames and datasets to be analyzed/written. 
clear;
bdata_source = 'Z:\Data\behavior\RC\';
analysis_outputs_dir = 'Z:\behavior_analysis\RC\';
CRP_fig_dir_base = 'Z:\behavior_analysis\RC\sessions_&_summaries\';

time_before_ms = 2000; %defines the window around the cue presentation which will be taken for plotting
time_after_ms = 3000;

CRP_training_sessions_lists_SJ;

%check and make sure the figure destinations exist
[this_mouse, ~] = select_mouse_info(sessions{1});
session_fig_dir = [CRP_fig_dir_base, 'img', this_mouse, '_sessions\'];
sum_fig_dir = [CRP_fig_dir_base, 'img', this_mouse, '_summary'];
if exist([CRP_fig_dir_base, 'img', this_mouse, '_sessions'], 'file') == 0
    mkdir([CRP_fig_dir_base, 'img', this_mouse, '_sessions']);
end
if exist([CRP_fig_dir_base, 'img', this_mouse, '_summary'], 'file') == 0
    mkdir([CRP_fig_dir_base, 'img', this_mouse, '_summary']);
end

gap_inx = gap_day_finder(sessions);

%% SECTION TWO - 

for ii = 1:length(sessions)
    sessions{ii}
    licktimes_dir  = [analysis_outputs_dir sessions{ii} '_licktimes'];
    b_data = get_bx_data_sj(bdata_source, sessions{ii});  %find the correct behavior file and loads it.
    [this_mouse, this_date] = select_mouse_info(sessions{ii}); % gets the mouse number and date of the experiment
    this_sess = sessions{ii};
    
    trial_start = round(double(cell2mat(b_data.tThisTrialStartTimeMs)));  %finds each trial's start time in free floating MWorks time
    num_trials  = length(b_data.tThisTrialStartTimeMs); %gets vectors for # of trials. Includes unimaged trials.
    
    % --- stores the time of the beginning of the first trial in MWorks time.
    bx_start_MWorks_time  = trial_start(1);
    
    %use time of first frame to align licking times to start of imaging
    lickTimes=[]; % in Ms
    if isfield(b_data, 'lickometerTimesUs')  %determine if input is structure array field
        for kk = 2:length(b_data.lickometerTimesUs) % get rid of the first trial
            lickTimes = [lickTimes cell2mat(b_data.lickometerTimesUs(kk))/1000];
        end
        lickTimes = double(lickTimes)-bx_start_MWorks_time;
    end
    
    %Collects various events during session
    hold_start = double(cell2mat(b_data.holdStartsMs)) - bx_start_MWorks_time;
    hold_time  = double(cell2mat(b_data.holdTimesMs));   %duration of the "lever hold" on that trial. meaningless here except to calculate cue onset
    react_time = double(cell2mat(b_data.reactTimesMs));
    release_time = hold_start + hold_time;
    cue_presentation = release_time-react_time;   %=====================THIS MEANS CUM HISTs ARE ALIGNED TO CUE ONSET
    TFT_ms = b_data.tooFastTimeMs; %this variable lengthens the time that the cue is on the screen
    cue_rew_int = b_data.RewardDelayDurationMs + round(mean(react_time),-1); %total interval between cueonset and reward delivery, rewardDelayDurationMs should always be zero
    
    if exist(licktimes_dir)
        save(licktimes_dir, 'lickTimes', '-append');
    else
        save(licktimes_dir, 'lickTimes');
    end
    
    %isolate the time of cue onset and divide licking into trials as such
    licks_by_trial = zeros(length(cue_presentation)-1,(time_before_ms+time_after_ms+1)); %dim1=trial# dim2=ms
    first_lick_by_trial = zeros(1, length(cue_presentation)-1);
    trials_where_licking_preceded_reward = zeros(1, length(cue_presentation)-1); % if the mice keeps licking from a time before the reward until after reward, a trial can be both trails_where_licking_preceded_reward and trails_with_licking_soon_after_reward
    trials_with_licking_soon_after_reward = zeros(1, length(cue_presentation)-1);
    for kk = 2:length(cue_presentation)-1 %look at all trials except the first and last one.
        %find all the licking events Xms before and Yms after cue presentation
        licks_this_window = lickTimes(find(lickTimes>cue_presentation(kk)-time_before_ms & lickTimes<cue_presentation(kk)+time_after_ms));
        alignment_this_trial = cue_presentation(kk)-(time_before_ms+1); %subtract off this time so that the lick times are converted to numbers which will correspond to there index in licks_by_trial
        licks_this_window = licks_this_window - alignment_this_trial;
        licks_by_trial(kk, licks_this_window) = 1;
        if ~isempty(find(lickTimes>cue_presentation(kk) & lickTimes<cue_presentation(kk)+cue_rew_int))
            trials_where_licking_preceded_reward(kk) = kk;
        end
        if ~isempty(find(lickTimes>cue_presentation(kk)+ cue_rew_int & lickTimes < cue_presentation(kk)+ cue_rew_int +250))
           trials_with_licking_soon_after_reward(kk) = kk;
        end
    end

    licks_by_trial_rewarded = licks_by_trial;

    %REWARD generate cumulative histograms of licking for rewarded trials
    avg_licks_per_ms_rew = cumsum(mean(licks_by_trial_rewarded));%average across trials, get a vector licks each ms, and then cumulative sum
    all_trials_lick_hist = cumsum(licks_by_trial_rewarded,2);
    
    %REWARDED plot licking histograms for rewarded trials
    x_axis_range = -1*time_before_ms:time_after_ms;
    if exist([session_fig_dir, sessions{ii}, '_rew_cum_hist.fig'], 'file')==0
        figure;
        plot(x_axis_range, avg_licks_per_ms_rew, 'r', 'LineWidth', 3);  %plot average cum lick hist for t by t analysis
        ylabel('cumulative # of licks');
        xlabel('time (ms) relative to release cue onset');
        title(['Rewarded Trials: cumulative hist of licking img', this_mouse, ' ' this_date ' n=' num2str(size(licks_by_trial_rewarded,1))]);
        hold on;
        for kk = 1:6:size(all_trials_lick_hist,1)
            plot(x_axis_range, all_trials_lick_hist(kk,:), 'Color', [0,0,0]+(1-(kk/size(all_trials_lick_hist,1))));
        end
        if b_data.soundTargetAmplitude >0 && b_data.gratingSpeedDPS == 0 % plays sounds, no visual stim
            vline(cue_rew_int, 'b');
        elseif b_data.soundTargetAmplitude >0 && b_data.gratingSpeedDPS > 0 % play both sound and visual stim
            vline(cue_rew_int, 'r'); 
        end
        savefig([session_fig_dir, sessions{ii}, '_rew_cum_hist']);
    end
    
    
    if exist([session_fig_dir, sessions{ii}, '_rew_lick_raster.fig'], 'file')==0
        %REWARDED trials lick raster plot
        figure;
        plotSpikeRaster(logical(licks_by_trial_rewarded), 'PlotType', 'vertline');
        if b_data.soundTargetAmplitude >0 && b_data.gratingSpeedDPS == 0 % plays sounds, no visual stim
            vline(time_before_ms+1, 'b');
        elseif b_data.soundTargetAmplitude >0 && b_data.gratingSpeedDPS > 0 % plays sounds and visual stim
            vline(time_before_ms+1, 'r');
        end
        vline(time_before_ms+1 + cue_rew_int, 'k');
        ylabel('trial # (descending)');
        xlabel('time (ms) colored=cue black=reward');
        title(['Rewarded Trials: lick time raster img', this_mouse, ' ' this_date ' n=' num2str(size(licks_by_trial_rewarded,1))]);
        savefig([session_fig_dir, sessions{ii}, '_rew_lick_raster']);
    end
   
    
%Determine the RTs for each trial this session
    RT_this_session=[];

% Below is Jake's original first lick method for reaction time. It has been
% replaced with a bursting method to capture an animals real reaction time
%{ 
    for kk = 1:num_trials-1 %look at each trial
         this_trial_RT = find(licks_by_trial(kk,[time_before_ms+1:end]), 1, 'first'); %look at all time points starting from cue onset for each trial, find the first lick
          RT_this_session = [RT_this_session, this_trial_RT]; %store the RTs
          if isempty(this_trial_RT)
              RT_this_session = [RT_this_session, NaN]; %store the RTs
          end
    end
 %}
    jakes_RT_filtering = true; %there is an elseif below that chops out early and late reaction times. toggle this
    for kk = 2:num_trials-1 %look at each trial
         this_trial_post_cue_licks = find(licks_by_trial(kk,[time_before_ms+1:end])); %look at all lick time points starting from cue onset for each trial
         if length(this_trial_post_cue_licks) >=3 %at least three licks needed to define a burst
            for i = 1:length(this_trial_post_cue_licks)-2 %look ahead as long as three licks can still be used 
                no_burst_flag = true; 
                if  this_trial_post_cue_licks(i+2)-this_trial_post_cue_licks(i) <= 500 %500ms window for bursting
                    no_burst_flag = false;
                    this_trial_RT = this_trial_post_cue_licks(i);
                    RT_this_session = [RT_this_session, this_trial_RT]; %store the RTs
                    break
                end
            end
            
            if no_burst_flag % flagged when there are many licks, but no bursts. Use first instead
               this_trial_RT = this_trial_post_cue_licks(1);
               RT_this_session = [RT_this_session, this_trial_RT];
            end

         elseif ~isempty(this_trial_post_cue_licks) %Use first lick in place of bursting
             this_trial_RT = this_trial_post_cue_licks(1);
             RT_this_session = [RT_this_session, this_trial_RT]; %store the RTs
         else
             RT_this_session = [RT_this_session, NaN]; 
         end
    end
    
    RT_this_session_raw = RT_this_session;
    RT_this_session_raw(find(isnan(RT_this_session_raw))) = [];
    save(['Z:\behavior_analysis\RC\RT_data\', sessions{ii}], 'RT_this_session_raw', 'time_before_ms');  

   
    %trim RT data to exclude misses and FA (false alarm)s  
    TFT_rate_this_session = length(find(RT_this_session<200))/(num_trials-1);
    TFT_rates = [TFT_rates, TFT_rate_this_session];

    if TFT_ms > 0 & b_data.RewardDelayDurationMs == 0 && jakes_RT_filtering %%% what is this?
        miss_rate_this_session = length(find(RT_this_session>1200))/(num_trials-1);
        RT_this_session = RT_this_session(find(RT_this_session>200 & RT_this_session<1200));
    elseif jakes_RT_filtering
        miss_rate_this_session = (length(find(RT_this_session>500))+no_lick_trials) / num_trials_b1;
        RT_this_session = RT_this_session(find(RT_this_session>200 & RT_this_session<500));
        
    end
    
    %Determine the FA and miss rates for this session. Store stats for summary statistics
    if jakes_RT_filtering
        miss_rates = [miss_rates, miss_rate_this_session];
    end
    RT_across_sessions = [RT_across_sessions, mean(RT_this_session)];
    RT_across_sessions_sem = [RT_across_sessions_sem, std(RT_this_session)/sqrt(size(RT_this_session,2))];
    std_of_RT_across_sessions = [std_of_RT_across_sessions, std(RT_this_session)];
    
    
    if b_data.soundTargetAmplitude >0 &&  b_data.gratingSpeedDPS > 0 % if train with CS1+CS2
        if jakes_RT_filtering
            miss_rates_2CS = [miss_rates_2CS, miss_rate_this_session_2CS];
        end
        RT_across_days_b2 = [RT_across_days_b2, mean(RT_this_sesssion_block2)];
        RT_across_days_sem_b2 = [RT_across_days_sem_b2, std(RT_this_sesssion_block2)/sqrt(size(RT_this_sesssion_block2,2))];
        std_of_RT_across_days_b2 = [std_of_RT_across_days_b2, std(RT_this_sesssion_block2)];
    end
    
    %plot RT within sessions relative to cue onset 
    if exist([session_fig_dir, sessions{ii}, '_RT_plot.fig'], 'file')==0
        figure; plot(RT_this_session);
        title(['RT for individual trials within for day ' this_date, this_mouse]);
        xlabel('trial #');
        ylabel('RT (ms)');
        savefig([session_fig_dir, sessions{ii}, '_RT_plot']);
    end

    
    %store trial-by-trial RTs for across sessions plot 
    if b_data.gratingSpeedDPS == 0
        RT_across_sessions = [RT_across_sessions, RT_this_session];
    elseif b_data.gratingSpeedDPS > 0
        RT_across_sessions_2CS = [RT_across_sessions_2CS,RT_this_session];
    elseif b_data.rewardDelayPercent == 0 && b_data.tooFastTimeMs > 0
        RT_across_sessions_delay = [RT_across_sessions_delay, RT_this_session];
    end

    
    %determine trial index for trials on different sessions and non-consecutive sessions
    if b_data.rewardDelayPercent == 0 && TFT_ms>0
        sessions_divider_inx_delay = [sessions_divider_inx_delay, length(RT_across_sessions_delay)+0.5];
        if ismember([ii+0.5], gap_inx)
            non_consecutive_inx_delay = [non_consecutive_inx_delay, length(RT_across_sessions_delay)];
        end

    elseif b_data.rewardDelayPercent == 0
        sessions_divider_inx = [sessions_divider_inx, length(RT_across_sessions)+0.5];
        if ismember([ii+0.5], gap_inx);
            non_consecutive_inx = [non_consecutive_inx, length(RT_across_sessions)];
        end
    end
    
    %----bin licking by 100ms windows relative to cue delivery----
    %determine  bin size and window size
    bin_size = 100; %number of ms to bin licking. 
    trial_start = trial_start - bx_start_MWorks_time;
    min_start_to_cue = min([cue_presentation-trial_start]); % different ITIs and hold times and other shit will make cue_presentation-trail_start variable
    min_cue_to_end = min([trial_start(2:end)-cue_presentation(1:end-1)]); % the end of this trial is the beginning of the next trial
    if min_start_to_cue > 20000
        pre_cue_window_lick = 20000;
    else 
        pre_cue_window_lick = floor([min_start_to_cue/bin_size])*bin_size;
    end
    if min_cue_to_end > 6000
        post_cue_window_lick = 5999;
    else 
        post_cue_window_lick = floor([min_cue_to_end/bin_size])*bin_size-1;
    end
    pre_cue_window_lick = 18000; % this is looking at a total of 24 seconds
    post_cue_window_lick = 5999;
    
    %get lick traces (1ms resolution) for trials
    full_trial_licks = zeros(length(cue_presentation)-1,(pre_cue_window_lick+post_cue_window_lick+1)); %dim1=trial# dim2=ms
    for kk = 2:length(cue_presentation)-1 %look at all trials except the first and last one.
        %find all the licking events Xms before and Yms after cue presentation
        licks_this_window = lickTimes(find(lickTimes>cue_presentation(kk)-pre_cue_window_lick & lickTimes<cue_presentation(kk)+post_cue_window_lick));
        alignment_this_trial = cue_presentation(kk)-(pre_cue_window_lick+1); %subtract off this time so that the lick times are converted to numbers which will correspond to there index in licks_by_trial
        licks_this_window = licks_this_window - alignment_this_trial;
        full_trial_licks(kk, licks_this_window) = 1;
    end
    if b_data.doBlock2 == 0
        full_trial_licks_rewarded = full_trial_licks;
    end
    
    %bin licking by 100ms bins and convert to licks/sec
    full_trial_licks_rewarded_sum = sum(full_trial_licks_rewarded,1);
    full_trial_licks_rewarded_bin = zeros(1,(length(full_trial_licks_rewarded_sum)/bin_size));
    cue_presentation_binned = (pre_cue_window_lick/bin_size)+1;
    iii=1;
    for kk = 1:bin_size:length(full_trial_licks_rewarded_sum)
        iii=iii+1;
        full_trial_licks_rewarded_bin(iii) = sum(full_trial_licks_rewarded_sum(kk:[kk+bin_size-1]));
    end
    
    %plot rewarded trials
    full_trial_licks_rewarded_bin = (full_trial_licks_rewarded_bin/size(full_trial_licks_rewarded,1))*(1000/bin_size); % convert to lick rate in Hz
    x_axis_bin = ([1:length(full_trial_licks_rewarded_bin)]-cue_presentation_binned)*(bin_size/1000);
    save(['Z:\behavior_analysis\RC\hist_data_across_animals\', sessions{ii}, 'rew_hist'], 'full_trial_licks_rewarded_bin','x_axis_bin');
    if exist([session_fig_dir, sessions{ii}, '_rew_lick_hist.fig'], 'file')==0
        figure; bar(x_axis_bin, full_trial_licks_rewarded_bin);
        title(['Rewarded trials: lick rate per ', num2str(bin_size), 'ms', this_date, ' img', this_mouse, 'n=', num2str(size(full_trial_licks_rewarded,1))]);
        xlabel('time (s) relative to cue onset');
        ylabel('lick rate (Hz)');
        vline(0,'k');
        savefig([session_fig_dir, sessions{ii}, '_rew_lick_hist']);
    end
     
    %store avg lick values in 500ms window just before cue on rewarded trials and during the ITI. 
    pre_cue_lick_window = full_trial_licks_rewarded(:,[pre_cue_window_lick-500:pre_cue_window_lick]);
    pre_cue_lick_window_sum = sum(pre_cue_lick_window,2);
    pre_cue_lick_window_avg_this_session = mean(pre_cue_lick_window_sum);
    pre_cue_lick_rate_sem_this_session = std(pre_cue_lick_window_sum)/sqrt(size(pre_cue_lick_window,1));
    pre_cue_lick_window_avg = [pre_cue_lick_window_avg ,pre_cue_lick_window_avg_this_session];
    pre_cue_lick_rate_sem = [pre_cue_lick_rate_sem ,pre_cue_lick_rate_sem_this_session];
    %now do it for mid iti lick window
    mid_iti = round(pre_cue_window_lick/2);
    iti_lick_window = full_trial_licks_rewarded(:,[mid_iti-500:mid_iti]);
    iti_lick_window_sum = sum(iti_lick_window,2);
    iti_lick_window_avg_this_session = mean(iti_lick_window_sum);
    iti_lick_rate_sem_this_session = std(iti_lick_window_sum)/sqrt(size(iti_lick_window,1));
    iti_lick_window_avg = [iti_lick_window_avg , iti_lick_window_avg_this_session];
    iti_lick_rate_sem = [iti_lick_rate_sem , iti_lick_rate_sem_this_session];
end 




%% save the pre/post cue data to be use in across animals summary
save(['Z:\behavior_analysis\RC\post_cue_lick_across_animals\', sessions{ii}(end-3:end), 'rew_om_post_cue_lick'], ...
    'avg_licks_post_cue', 'avg_licks_pre_cue', 'avg_licks_post_cue_sem', 'avg_licks_pre_cue_sem'); 
avg_licks_post_cue = avg_licks_post_cue(~isnan(avg_licks_post_cue)); %this allows the training day # to be aligned to the pre/post cue lick rate for that day across animals
avg_licks_pre_cue = avg_licks_pre_cue(~isnan(avg_licks_pre_cue));
avg_licks_post_cue_sem = avg_licks_post_cue_sem(~isnan(avg_licks_post_cue_sem));
avg_licks_pre_cue_sem  = avg_licks_pre_cue_sem(~isnan(avg_licks_pre_cue_sem));

%save RT std, FA rate and miss rate
if jakes_RT_filtering
save(['Z:\behavior_analysis\RC\RT_stats_across_animals\', sessions{ii}(end-3:end), 'RTstd_FA_misses'], ...
    'std_of_RT_across_sessions', 'std_of_RT_across_sessions_b2','TFT_rates', 'miss_rates', 'TFT_rates_b2', 'miss_rates_b2');
end


%% plot basic stats across sessions
figure; 
set(gcf, 'Position', [100, 100, 1000, 500]);
subplot(2,2,1);
errorbar(RT_across_sessions, RT_across_sessions_sem, 'b');
title('RT across sessions');
if b_data.gratingSpeedDPS > 0
    hold on; errorbar(RT_across_sessions_b2, RT_across_sessions_sem_b2, 'r');
end
xlabel('day');
ylabel('RT (ms)');
xlim([0.8 length(RT_across_sessions)+0.2]);

subplot(2,2,2); plot(std_of_RT_across_sessions, 'b');
title('std of RT across sessions');
if b_data.doBlock2 == 1
    hold on; plot(std_of_RT_across_sessions_b2, 'k');
end
xlabel('day');
ylabel('standard deviation');
xlim([0.8 length(std_of_RT_across_sessions)+0.2]);

subplot(2,2,3); plot(TFT_rates, 'b');
title('too fast lick rates as a % of all trials by day: FA = RT<200ms');
if b_data.doBlock2 == 1
    hold on; plot(TFT_rates_b2, 'k');
end
xlabel('day');
ylabel('% false alarms');
ylim([0 1]);
xlim([0.8 length(TFT_rates)+0.2]);

if jakes_RT_filtering
subplot(2,2,4); plot(miss_rates, 'b');
title('Miss rates as a % of all trials by day: misses = RT>1000ms');
if b_data.doBlock2 ==1
    hold on; plot(miss_rates_2CS, 'k');
end
xlabel('day');
ylabel('% misses');
ylim([0 1]);
xlim([0.8 length(miss_rates)+0.2]);
end 

savefig([sum_fig_dir, '\', this_mouse, '_RT_RTstd_misses_TFT']);

%% 
if length(sessions_divider_inx) > 0
    figure; set(gcf, 'Position', [100, 100, 1000, 500]);
    plot(RT_across_sessions);
    title(['RTs across sessions. No delay. img' this_mouse]);
    xlabel('trial num');
    ylabel('RT(ms)');
    vline(sessions_divider_inx, 'k');
    if ~isempty(non_consecutive_inx)
        vline(non_consecutive_inx, 'r');
    end
    ylim([200 1000]);
    savefig([sum_fig_dir, '\', this_mouse, '_RT_0ms_summary']);
end

if length(RT_across_sessions_delay) > 0 & b_data.rewardDelayPercent == 0
    figure; plot(RT_across_sessions_delay, 'b');
    title(['RTs across sessions. ', num2str(cue_rew_int), 'ms delay. img' this_mouse]);
    if b_data.doBlock2 == 1
        hold on; plot(RT_across_sessions_delay_b2, 'k');
        title(['RTs across sessions. ', num2str(cue_rew_int), 'ms delay. img' this_mouse, ': blue=CS+, black=CS-']);
        if TFT_ms > 0 & b_data.rewardDelayPercent == 0;
            title(['RTs across sessions. ', num2str(cue_rew_int), 'ms delay. img' this_mouse, ': blue=CS+, black=CS-']);
        end
    elseif TFT_ms > 0 & b_data.rewardDelayPercent == 0;
         title(['RTs across sessions. ', num2str(cue_rew_int), 'ms delay. ' this_mouse]);
    end
    xlabel('trial num');
    ylabel('RT(ms)');
    vline(sessions_divider_inx_delay, 'k');
    if ~isempty(non_consecutive_inx_delay)
        vline(non_consecutive_inx_delay, 'r');
    end
    savefig([sum_fig_dir, '\', this_mouse, '_RT_1s_drift_summary']);
end

%% plot # of licks in iti window vs pre-cue window
figure; errorbar(iti_lick_window_avg, iti_lick_rate_sem, 'k');
hold on; 
errorbar(pre_cue_lick_window_avg, pre_cue_lick_rate_sem, 'g');
title('licking in iti (black) vs pre-cue (green)');
ylabel('avg # of licks in 500ms window');
xlabel('day #');
xlim([0.8 length(pre_cue_lick_window_avg)+0.2]);
savefig([sum_fig_dir, '\', this_mouse, '_iti_vs_pre-cue_licking']);


%save variables for across animals RT plot
out_dir = 'Z:\behavior_analysis\RC\RT_across_animals\';
save([out_dir, this_mouse, 'RT_sem_across_sessions'], 'RT_across_sessions', 'RT_across_sessions_sem', 'RT_across_sessions_b2', 'RT_across_sessions_sem_b2');

