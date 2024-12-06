%% Modified version of the anaylysis overview for analyzing the licking behavior of cue-reward trials. SJ
% licks by trial: trial*ms, binary matrix, licking time points = 1
% licks_by_trial_rewarded = licks_by_trial;
% RT_this_session: 1*ntrials vector, reaction time for each trial. Only trials with a reaction time longer than 200ms and shorter than 1200ms are included.
% full_trial_licks_rewarded_bin: bin licking every 100ms, lick rate in Hz for every 100ms. Looking at a total of 24 seconds

% since juice==0 in the first trial, this script gets rid of the first trial

%TJW NEED TO DO:
%fix green line on cumulative hist lick fig
%plot a horizontal line on RT fig at 650
%'zoom in' on blue hist
%do across day graphs

%% SECTION ONE - assign pathnames and datasets to be analyzed/written. 
clear;
clc;
close all;
% bdata_source = 'Z:\home\shuyang\Data\behavior\RC\';
% analysis_outputs_dir = 'Z:\home\shuyang\behavior_analysis\RC\';
% RC_fig_dir_base = 'Z:\home\shuyang\behavior_analysis\RC\sessions_&_summaries\';

%bdata_source = 'Z:\home\tj\cc_data';
% bdata_source = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\cc_data';
bdata_source = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
%analysis_outputs_dir = 'Z:\home\tj\Analysis\Behavior\class_cond';
analysis_outputs_dir = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Behavior\class_cond';
%RC_fig_dir_base = 'Z:\home\tj\Analysis\Behavior\class_cond\sessions_&_summaries';
RC_fig_dir_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Behavior\class_cond\sessions_&_summaries';


time_before_ms = 2000; %defines the window around the cue presentation which will be taken for plotting
time_after_ms = 3000;
RC_training_sessions_lists_SJ_tjw;

%check and make sure the figure destinations exist
[this_mouse, ~] = select_mouse_info(sessions{1});
session_fig_dir = [RC_fig_dir_base, 'img', this_mouse, '_sessions\'];
sum_fig_dir = [RC_fig_dir_base, 'img', this_mouse, '_summary'];
if exist([RC_fig_dir_base, 'img', this_mouse, '_sessions'], 'file') == 0
    mkdir([RC_fig_dir_base, 'img', this_mouse, '_sessions']);
end
if exist([RC_fig_dir_base, 'img', this_mouse, '_summary'], 'file') == 0
    mkdir([RC_fig_dir_base, 'img', this_mouse, '_summary']);
end

%gap_inx = gap_day_finder(sessions);

%% SECTION TWO - 

for ii = 1:length(sessions)
    sessions{ii}
    licktimes_dir  = [analysis_outputs_dir sessions{ii} '_licktimes'];
    b_data = get_bx_data_tj(bdata_source, sessions{ii});  %find the correct behavior file and loads it.
    [this_mouse, this_date] = select_mouse_info(sessions{ii}); % gets the mouse number and date of the experiment
    mouse = this_mouse(1:4);
    assert(str2num(mouse)==testDay_mouse(1));
    
    this_sess = sessions{ii};
    
    trial_start = round(double(cell2mat(b_data.tThisTrialStartTimeMs)));  %finds each trial's start time in free floating MWorks time
    trial_start(1)= [];% get rid of first trial
    num_trials  = length(b_data.tThisTrialStartTimeMs); %gets vectors for # of trials. Includes unimaged trials.
    soundamp = b_data.soundTargetAmplitude;
    block2 = cell2mat(b_data.tBlock2TrialNumber);
    % determine if it's interleave day
    %if b_data.block2TrPer80Level1 > 0 % if this number >0 than it's interleaved trials that day
        auditrials = [];
        vistrials = [];
        avtrials = [];
        for tri = 1:num_trials-1 %ignore the first and last trial
            if tri ==1
                continue; % get rid of the first trial regardless 
            %elseif block2(tri) == 1 && soundamp{tri} > 0 % only tone
            %    auditrials = [auditrials, tri]; % this index is the index of the raw trial number, which includes everything
            elseif block2(tri) == 0 %&& soundamp{tri} == 0 %only visual
                vistrials = [vistrials, tri];
            %elseif block2(tri) == 0 && soundamp{tri} >0 % both tone and visual
            %    avtrials = [avtrials, tri];
            end
        end
        veri = (length(auditrials) + length(vistrials) + length(avtrials) == num_trials-2)
        save([analysis_outputs_dir sessions{ii} '_IL_trial_inx.mat'],'auditrials','vistrials','avtrials');
   % elseif iscell(soundamp)
   %     soundamp = soundamp{end}; %changed mworks program on 11/30/2020, soundTargetAmplitude used to be persistent during the experiment and it was a single number. It will become a cell after this change. But you don't need this unless it's test day.
   % end
    
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
    if exist(licktimes_dir)
        save(licktimes_dir, 'lickTimes', '-append');
    else
        save(licktimes_dir, 'lickTimes');
    end
    %Collects various events during session
    hold_start = double(cell2mat(b_data.holdStartsMs)) - bx_start_MWorks_time;
    hold_time  = double(cell2mat(b_data.holdTimesMs));   %duration of the "lever hold" on that trial. meaningless here except to calculate cue onset
    react_time = double(cell2mat(b_data.reactTimesMs));
    release_time = hold_start + hold_time;
    cue_presentation = release_time-react_time;   %=====================THIS MEANS CUM HISTs ARE ALIGNED TO CUE ONSET
    cue_presentation(1) = []; % get rid of the first trial
    TFT_ms = b_data.tooFastTimeMs; %this variable lengthens the time that the cue is on the screen
    cue_rew_int = b_data.RewardDelayDurationMs + round(mean(react_time),-1); %total interval between cueonset and reward delivery, rewardDelayDurationMs should always be zero
    
    %isolate the time of cue onset and divide licking into trials as such
    licks_by_trial = zeros(length(cue_presentation)-1,(time_before_ms+time_after_ms+1)); %dim1=trial# dim2=ms
    trials_where_licking_preceded_reward = zeros(1, length(cue_presentation)-1); % if the mice keeps licking from a time before the reward until after reward, a trial can be both trails_where_licking_preceded_reward and trails_with_licking_soon_after_reward
    trials_with_licking_soon_after_reward = zeros(1, length(cue_presentation)-1);
    for kk = 1:length(cue_presentation)-1 %look at all trials except the last one.
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

    
    
    if b_data.block2TrPer80Level1 > 0 % if this number >0 than it's interleaved trials that day
        licks_by_trial_rewarded_audi = licks_by_trial(auditrials-1,:);% licks by trial doesn't have the first trial of the raw data, but auditrials is the index of the raw data
        licks_by_trial_rewarded_vis = licks_by_trial(vistrials-1,:);
        licks_by_trial_rewarded_av = licks_by_trial(avtrials-1,:);   
    end
    licks_by_trial_rewarded = licks_by_trial;
    
    %REWARD generate cumulative histograms of licking for rewarded trials
    avg_licks_per_ms_rew = cumsum(mean(licks_by_trial_rewarded));%average across trials, get a vector licks each ms, and then cumulative sum
    all_trials_lick_hist = cumsum(licks_by_trial_rewarded,2);
    
    %REWARDED plot licking histograms for rewarded trials
    x_axis_range = -1*time_before_ms:time_after_ms;
    if exist([session_fig_dir, sessions{ii}, '_rew_cum_hist.fig'], 'file')==0
        if b_data.block2TrPer80Level1 > 0
        else
            figure;
            plot(x_axis_range, avg_licks_per_ms_rew, 'r', 'LineWidth', 3);  %plot average cum lick hist for t by t analysis
            ylabel('cumulative # of licks');
            xlabel('time (ms) relative to release cue onset');
            title(['Rewarded Trials: cumulative hist of licking img', this_mouse, ' ' this_date ' n=' num2str(size(licks_by_trial_rewarded,1))]);
            hold on;
            for kk = 1:6:size(all_trials_lick_hist,1)
                plot(x_axis_range, all_trials_lick_hist(kk,:), 'Color', [0,0,0]+(1-(kk/size(all_trials_lick_hist,1))));
            end
            
                vline(time_before_ms+1, 'g');
            
            savefig([session_fig_dir, sessions{ii}, '_rew_cum_hist']);
        end
    end
    
    
    if exist([session_fig_dir, sessions{ii}, '_rew_lick_raster.fig'], 'file')==0
        if b_data.block2TrPer80Level1 > 0 % if this number >0 than it's interleaved trials that day
            figure;subplot(1,3,1);
            plotSpikeRaster(logical(licks_by_trial_rewarded_audi), 'PlotType', 'vertline'); vline(time_before_ms+1,'b'); vline(time_before_ms+1 + cue_rew_int, 'k');
            ylabel('trial # (descending)');
            subplot(1,3,2);
            plotSpikeRaster(logical(licks_by_trial_rewarded_vis), 'PlotType', 'vertline'); vline(time_before_ms+1,'g'); vline(time_before_ms+1 + cue_rew_int, 'k');
            title(['Rewarded Trials: lick time raster img', this_mouse, ' ' this_date ' n=' num2str(size(licks_by_trial,1))]);
            xlabel('time (ms) colored=cue black=reward');
            subplot(1,3,3);
            plotSpikeRaster(logical(licks_by_trial_rewarded_av), 'PlotType', 'vertline'); vline(time_before_ms+1,'r'); vline(time_before_ms+1 + cue_rew_int, 'k');
            savefig([session_fig_dir, sessions{ii}, '_rew_lick_raster']);
        else
            %REWARDED trials lick raster plot
            figure;
            plotSpikeRaster(logical(licks_by_trial_rewarded), 'PlotType', 'vertline');
            
                vline(time_before_ms+1, 'g');
            
            if b_data.trialLaserPowerMw > 0 && b_data.fixedReqHoldTimeMs ~= 1900 % if doing optogenetics
                laser_bf_cue = b_data.fixedReqHoldTimeMs;
                vline(time_before_ms+1 - laser_bf_cue, 'c');
            end
            vline(time_before_ms+1 + cue_rew_int, 'k');
            ylabel('trial # (descending)');
            xlabel('time (ms) colored=cue black=reward');
            title(['Rewarded Trials: lick time raster img', this_mouse, ' ' this_date ' n=' num2str(size(licks_by_trial_rewarded,1))]);
            savefig([session_fig_dir, sessions{ii}, '_rew_lick_raster']);
        end
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
    
    for kk = 1:num_trials-2 %look at each trial, except the first and last one
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
    RT_this_session_raw(find(isnan(RT_this_session_raw))) = -1;
  %  if b_data.block2TrPer80Level1 > 0
  %      RT_raw_audi = RT_this_session_raw(auditrials-1);
        RT_raw_vis = RT_this_session_raw(vistrials-1);
  %      RT_raw_av = RT_this_session_raw(avtrials-1);
  %  end
    
    %plot raw RT within sessions relative to cue onset
    if exist([session_fig_dir, sessions{ii}, '_rawRT_plot.fig'], 'file')==0
        if b_data.block2TrPer80Level1 > 0
            figure;
            subplot(1,3,1);plot(RT_raw_audi,'b');ylabel('RT (ms)');
            subplot(1,3,2);plot(RT_raw_vis,'g');title(['RT for all individual trials within for day ' this_date, this_mouse]);xlabel('trial #');
            subplot(1,3,3);plot(RT_raw_av,'r');
            savefig([session_fig_dir, sessions{ii}, '_rawRT_plot']);
        else
            figure; plot(RT_this_session_raw);
            title(['RT for all individual trials within for day ' this_date, this_mouse]);
            xlabel('trial #');
            ylabel('RT (ms)');
            savefig([session_fig_dir, sessions{ii}, '_rawRT_plot']);
            %         prompt='delete some bad trials? Y/N: \n';response = input(prompt,'s'); % sometimes the mice isn't engaged that last several trials, get rid of those
            %         switch response
            %             case 'Y'; delete_bad_trials = 1;
            %             case 'N'; delete_bad_trials = 0;
            %         end
            %         if delete_bad_trials == 1
            %             prompt='\nEnter bad trials with numbers separated by commas or spaces: \n';
            %             userResponse = input(prompt,'s');
            %             % Convert any commas to spaces.
            %             userResponse = strrep(userResponse, ',', ' ');
            %             % Convert strings to numbers.
            %             bad_trials = sscanf(userResponse, '%f');
            %             num_trials = num_trials-2-length(bad_trials);
            %             RT_this_session_raw (bad_trials) = [];
            %             RT_this_session(bad_trials) = [];
            %             %also delete licking data as needed
            %             cue_presentation (bad_trials) = [];
            %
            %         else
            %             num_trials = num_trials - 2; % except the first and last trial
            %             bad_trials = 0;
            %         end
        end
    end
    %save(['Z:\home\shuyang\behavior_analysis\RC\RT_data\', sessions{ii}], 'RT_this_session_raw', 'time_before_ms');%,'bad_trials');
    save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Behavior\RT_data\', sessions{ii}], 'RT_this_session_raw', 'time_before_ms');%,'bad_trials');
    %trim RT data to exclude misses and FA (false alarm)s  
    %if b_data.block2TrPer80Level1 > 0
    %    TFT_rate_audi_this_session = length(find(RT_raw_audi<200))/(length(RT_raw_audi)-1);
        TFT_rate_vis_this_session = length(find(RT_raw_vis<200))/(length(RT_raw_vis)-1);
    %    TFT_rate_av_this_session = length(find(RT_raw_av<200))/(length(RT_raw_av)-1);
    %    TFT_rates_IL = [TFT_rates_IL, TFT_rate_audi_this_session,TFT_rate_vis_this_session,TFT_rate_av_this_session];
        TFT_rates_IL = [TFT_rates_IL, TFT_rate_vis_this_session];
        
    %    miss_rate_audi_this_session = length(find(RT_raw_audi>1200))/(length(RT_raw_audi)-1);
        miss_rate_vis_this_session = length(find(RT_raw_vis>1200))/(length(RT_raw_vis)-1);
    %    miss_rate_av_this_session = length(find(RT_raw_av>1200))/(length(RT_raw_av)-1);
    %    miss_rates_IL = [miss_rates_IL,miss_rate_audi_this_session,miss_rate_vis_this_session,miss_rate_av_this_session ];
        miss_rates_IL = [miss_rates_IL,miss_rate_vis_this_session];
        
    %    RT_audi_this_session = RT_raw_audi(find(RT_raw_audi>200 & RT_raw_audi < 1200));
        RT_vis_this_session = RT_raw_vis(find(RT_raw_vis>200 & RT_raw_vis < 1200));
    %    RT_av_this_session = RT_raw_av(find(RT_raw_av>200 & RT_raw_av < 1200));
        RT_across_sessions_IL = [mean(RT_vis_this_session)];
    %    RT_across_sessions_IL = [mean(RT_audi_this_session),mean(RT_vis_this_session),mean(RT_av_this_session)];
        RT_across_sessions_sem_IL = [RT_across_sessions_sem_IL, std(RT_vis_this_session)/sqrt(length(RT_vis_this_session))];
    %    RT_across_sessions_sem_IL = [RT_across_sessions_sem_IL, std(RT_audi_this_session)/sqrt(length(RT_audi_this_session)),...
    %        std(RT_vis_this_session)/sqrt(length(RT_vis_this_session)),std(RT_av_this_session)/sqrt(length(RT_av_this_session))];
        std_of_RT_across_sessions_IL = [std_of_RT_across_sessions_IL, std(RT_vis_this_session)];
    %    std_of_RT_across_sessions_IL = [std_of_RT_across_sessions_IL, std(RT_audi_this_session),std(RT_vis_this_session),std(RT_av_this_session)];
    %else
        TFT_rate_this_session = length(find(RT_this_session<200))/(num_trials-1);
        %if soundamp >0 && b_data.gratingSpeedDPS == 0 && testDay_mouse(ii+1)==0 % plays sounds, no visual stim
            TFT_rates = [TFT_rates, TFT_rate_this_session];
%         elseif soundamp >0 && b_data.gratingSpeedDPS > 0 && testDay_mouse(ii+1)==0 % plays sounds and visual stim
%             TFT_rates_2CS = [TFT_rates_2CS, TFT_rate_this_session];
%         elseif testDay_mouse(ii+1)==1
%             TFT_rates_testDay = [TFT_rates_testDay,TFT_rate_this_session];
        %end
        
        if TFT_ms > 0 && b_data.RewardDelayDurationMs == 0 && jakes_RT_filtering %%% what is this?
            miss_rate_this_session = length(find(RT_this_session>1200))/(num_trials-1);
            RT_this_session = RT_this_session(find(RT_this_session_raw>200 & RT_this_session_raw<1200));
        elseif jakes_RT_filtering
            miss_rate_this_session = (length(find(RT_this_session>500))+no_lick_trials) / num_trials_b1;
            RT_this_session = RT_this_session(find(RT_this_session>200 & RT_this_session<500));
        end
        
        %Determine the FA and miss rates for this session. Store stats for summary statistics
        if jakes_RT_filtering
            %if soundamp >0 && b_data.gratingSpeedDPS == 0 && testDay_mouse(ii+1)==0 % plays sounds, no visual stim
                miss_rates = [miss_rates, miss_rate_this_session];
%             elseif soundamp >0 && b_data.gratingSpeedDPS > 0 && testDay_mouse(ii+1)==0 % plays sounds and visual stim
%                 miss_rates_2CS = [miss_rates_2CS, miss_rate_this_session];
%             elseif testDay_mouse(ii+1)==1
%                 miss_rates_testDay = [miss_rates_testDay, miss_rate_this_session];
%             end
        end
        %if soundamp >0 && b_data.gratingSpeedDPS == 0 && testDay_mouse(ii+1)==0 % plays sounds, no visual stim
            RT_across_sessions = [RT_across_sessions, mean(RT_this_session)];
            RT_across_sessions_sem = [RT_across_sessions_sem, std(RT_this_session)/sqrt(size(RT_this_session,2))];
            std_of_RT_across_sessions = [std_of_RT_across_sessions, std(RT_this_session)];
%         elseif soundamp > 0 && b_data.gratingSpeedDPS > 0 && testDay_mouse(ii+1)==0 % plays sounds and visual stim
%             RT_across_sessions_2CS = [RT_across_sessions_2CS, mean(RT_this_session)];
%             RT_across_sessions_sem_2CS = [RT_across_sessions_sem_2CS, std(RT_this_session)/sqrt(size(RT_this_session,2))];
%             std_of_RT_across_sessions_2CS = [std_of_RT_across_sessions_2CS, std(RT_this_session)];
%         elseif testDay_mouse(ii+1)==1
%             RT_across_sessions_testDay = [RT_across_sessions_testDay, mean(RT_this_session)];
%             RT_across_sessions_sem_testDay = [RT_across_sessions_sem_testDay, std(RT_this_session)/sqrt(size(RT_this_session,2))];
%             std_of_RT_across_sessions_testDay = [std_of_RT_across_sessions_testDay, std(RT_this_session)];
        %end
        
        %plot RT within sessions relative to cue onset
        if exist([session_fig_dir, sessions{ii}, '_RT_plot.fig'], 'file')==0
            figure; plot(RT_this_session);
            title(['RT for individual trials within for day ' this_date, this_mouse]);
            xlabel('trial #');
            ylabel('RT (ms)');
            savefig([session_fig_dir, sessions{ii}, '_RT_plot']);
        end
%    end
    
    %store trial-by-trial RTs for across sessions plot 
    if b_data.rewardDelayPercent == 0
        RT_across_mulsessions = [RT_across_mulsessions, RT_this_session];
    elseif b_data.rewardDelayPercent == 0 && b_data.tooFastTimeMs > 0
        RT_across_sessions_delay = [RT_across_sessions_delay, RT_this_session];
    end

    
    %determine trial index for trials on different sessions and non-consecutive sessions
%     if b_data.rewardDelayPercent == 0 && TFT_ms>0
%         sessions_divider_inx_delay = [sessions_divider_inx_delay, length(RT_across_sessions_delay)+0.5];
%         if ismember([ii+0.5], gap_inx)
%             non_consecutive_inx_delay = [non_consecutive_inx_delay, length(RT_across_sessions_delay)];
%         end
% 
%     elseif b_data.rewardDelayPercent == 0
%         sessions_divider_inx = [sessions_divider_inx, length(RT_across_mulsessions)+0.5];
%         if ismember([ii+0.5], gap_inx)
%             non_consecutive_inx = [non_consecutive_inx, length(RT_across_mulsessions)];
%         end
%     end
    
    %----bin licking by 100ms windows relative to cue delivery----
    %determine  bin size and window size
    bin_size = 100; %number of ms to bin licking. 
    trial_start = trial_start - bx_start_MWorks_time;
    min_start_to_cue = min([cue_presentation-trial_start]); % different ITIs and hold times and other shit will make cue_presentation-trail_start variable
    min_cue_to_end = min([trial_start(2:end)-cue_presentation(1:end-1)]); % the end of this trial is the beginning of the next trial
%     if min_start_to_cue > 20000
%         pre_cue_window_lick = 20000;
%     else 
%         pre_cue_window_lick = floor([min_start_to_cue/bin_size])*bin_size;
%     end
%     if min_cue_to_end > 6000
%         post_cue_window_lick = 5999;
%     else 
%         post_cue_window_lick = floor([min_cue_to_end/bin_size])*bin_size-1;
%     end
    pre_cue_window_lick = 18000; % this is looking at a total of 24 seconds
    post_cue_window_lick = 5999;
    
    %get lick traces (1ms resolution) for trials
    full_trial_licks = zeros(length(cue_presentation)-1,(pre_cue_window_lick+post_cue_window_lick+1)); %dim1=trial# dim2=ms
    for kk = 1:length(cue_presentation)-1 %look at all trials except the first one.
        %find all the licking events Xms before and Yms after cue presentation
        licks_this_window = lickTimes(find(lickTimes>cue_presentation(kk)-pre_cue_window_lick & lickTimes<cue_presentation(kk)+post_cue_window_lick));
        alignment_this_trial = cue_presentation(kk)-(pre_cue_window_lick+1); %subtract off this time so that the lick times are converted to numbers which will correspond to there index in licks_by_trial
        licks_this_window = licks_this_window - alignment_this_trial;
        full_trial_licks(kk, licks_this_window) = 1;
    end
    full_trial_licks_rewarded = full_trial_licks;
    %bin licking by 100ms bins and convert to licks/sec
    full_trial_licks_rewarded_sum = sum(full_trial_licks_rewarded,1);% sum across trials
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
    %save(['Z:\home\tj\Analysis\Behavior\hist_data_across_animals\', sessions{ii}, 'rew_hist'], 'full_trial_licks_rewarded_bin','x_axis_bin');
    save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Behavior\hist_data_across_animals\', sessions{ii}, 'rew_hist'], 'full_trial_licks_rewarded_bin','x_axis_bin');
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




%% save RT std, FA rate and miss rate
if jakes_RT_filtering
save(['Z:\home\tj\Analysis\Behavior\RT_stats_across_animals\', sessions{ii}(end-3:end), 'RTstd_FA_misses'], ...
    'std_of_RT_across_sessions', 'std_of_RT_across_sessions_2CS','TFT_rates',...
    'miss_rates', 'TFT_rates_2CS', 'miss_rates_2CS','RT_across_sessions','RT_across_sessions_2CS',...
    'RT_across_sessions_testDay','std_of_RT_across_sessions_testDay','TFT_rates_testDay',...
    'miss_rates_testDay','RT_across_sessions_IL','std_of_RT_across_sessions_IL',...
    'TFT_rates_IL','miss_rates_IL');
end

%% plot basic stats across sessions
x = 1:1:length(RT_across_sessions);
figure; 
set(gcf, 'Position', [100, 100, 1000, 500]);
subplot(2,2,1);
errorbar(x,RT_across_sessions, RT_across_sessions_sem, 'b'); hold on;
title(['RT across sessions img' this_mouse]);
% if ~isempty(RT_across_sessions_2CS)
%     x2 = x(end)+(1:1:length(RT_across_sessions_2CS));
%     errorbar(x2,RT_across_sessions_2CS, RT_across_sessions_sem_2CS, 'r');
% end
% if ~isempty(RT_across_sessions_testDay)
%     x3 = x2(end) +(1:1:length(RT_across_sessions_testDay));
%     disp('please make sure color is correct');
%     color = {'r','b','g'};
%     for t = 1:length(RT_across_sessions_testDay)
%         errorbar(x3(t),RT_across_sessions_testDay(t),RT_across_sessions_sem_testDay(t),'.','MarkerSize',15,'Color',color{t});
%     end
% end
% if ~isempty(RT_across_sessions_IL)
%     x4 = x3(end)+(1:1:length(RT_across_sessions_IL));
%     color2 = {[0.0314 0.1137 0.3451],[0 0.4275 0.1725],[0.6000 0 0.0510]}; %the order should always be blue,green,red. YlGnBu13, BuGn13, Reds13
%     for s = 1:length(RT_across_sessions_IL)
%         errorbar(x4(s),RT_across_sessions_IL(s),RT_across_sessions_sem_IL(s),'.','MarkerSize',15,'Color',color2{s});
%     end
% end
xlabel('session');
ylabel('RT (ms)');
xlim([0.8 length(RT_across_sessions)+length(RT_across_sessions_2CS)+length(RT_across_sessions_testDay)+length(RT_across_sessions_IL)+0.2]);

subplot(2,2,2); plot(x,std_of_RT_across_sessions, 'b');
% if ~isempty(RT_across_sessions_2CS)
%     hold on; plot(x2,std_of_RT_across_sessions_2CS, 'r');
% end
% if ~isempty(RT_across_sessions_testDay)
%     for t = 1:length(RT_across_sessions_testDay)
%         scatter(x3(t),std_of_RT_across_sessions_testDay(t),'o','filled','MarkerFaceColor',color{t});
%     end
% end
% if ~isempty(RT_across_sessions_IL)
%     for s = 1:length(RT_across_sessions_IL)
%         scatter(x4(s),std_of_RT_across_sessions_IL(s),'o','filled','MarkerFaceColor',color2{s});
%     end
% end
title('std of RT across sessions');
xlabel('session');
ylabel('standard deviation');
xlim([0.8 length(RT_across_sessions)+length(RT_across_sessions_2CS)+length(RT_across_sessions_testDay)+length(RT_across_sessions_IL)+0.2]);

subplot(2,2,3); plot(x,TFT_rates_IL, 'b');
% if ~isempty(RT_across_sessions_2CS)
%     hold on; plot(x2,TFT_rates_2CS, 'r');
% end
% if ~isempty(RT_across_sessions_testDay)
%     for t = 1:length(TFT_rates_testDay)
%         scatter(x3(t),TFT_rates_testDay(t),'o','filled','MarkerFaceColor',color{t},'MarkerEdgeColor',color{t});
%     end
% end
% if ~isempty(RT_across_sessions_IL)
%     for s = 1:length(TFT_rates_IL)
%         scatter(x4(s),TFT_rates_IL(s),'o','filled','MarkerFaceColor',color2{s},'MarkerEdgeColor',color2{s});
%     end
% end
title('too fast lick rates as a % of all trials by day: FA = RT<200ms');
xlabel('session');
ylabel('% false alarms');
ylim([0 1]);
xlim([0.8 length(RT_across_sessions)+length(RT_across_sessions_2CS)+length(RT_across_sessions_testDay)+length(RT_across_sessions_IL)+0.2]);

if jakes_RT_filtering
    subplot(2,2,4); plot(miss_rates, 'b');
%     if ~isempty(RT_across_sessions_2CS)
%         hold on; plot(x2,miss_rates_2CS, 'r');
%     end
%     if ~isempty(RT_across_sessions_testDay)
%         for t = 1:length(TFT_rates_testDay)
%             scatter(x3(t),miss_rates_testDay(t),'o','filled','MarkerFaceColor',color{t},'MarkerEdgeColor',color{t});
%         end
%     end
%     if ~isempty(RT_across_sessions_IL)
%         for s = 1:length(TFT_rates_IL)
%             scatter(x4(s),miss_rates_IL(s),'o','filled','MarkerFaceColor',color2{s},'MarkerEdgeColor',color2{s});
%         end
%     end
    title('Miss rates as a % of all trials by day: misses = RT>1000ms');
    xlabel('session');
    ylabel('% misses');
    ylim([0 1]);
    xlim([0.8 length(RT_across_sessions)+length(RT_across_sessions_2CS)+length(RT_across_sessions_testDay)+length(RT_across_sessions_IL)+0.2]);
end
%supertitle([this_mouse(1:4) '  blue=CS1, red=CS1+CS2, green=CS2']);
savefig([sum_fig_dir, '\', this_mouse, '_RT_RTstd_misses_TFT']);


%% plot # of licks in iti window vs pre-cue window
figure; errorbar(iti_lick_window_avg, iti_lick_rate_sem, 'k');
hold on; 
errorbar(pre_cue_lick_window_avg, pre_cue_lick_rate_sem, 'g');
title(['licking in iti (black) vs pre-cue (green) img' this_mouse]);
ylabel('avg # of licks in 500ms window');
xlabel('day #');
xlim([0.8 length(pre_cue_lick_window_avg)+0.2]);
savefig([sum_fig_dir, '\', this_mouse, '_iti_vs_pre-cue_licking']);


%save variables for across animals RT plot
out_dir = 'Z:\home\tj\Analysis\Behavior\RT_across_animals\';
save([out_dir, this_mouse, 'RT_sem_across_sessions'], 'RT_across_sessions', ...
    'RT_across_sessions_sem', 'RT_across_sessions_2CS', 'RT_across_sessions_sem_2CS',...
    'RT_across_sessions_IL','RT_across_sessions_sem_IL');

