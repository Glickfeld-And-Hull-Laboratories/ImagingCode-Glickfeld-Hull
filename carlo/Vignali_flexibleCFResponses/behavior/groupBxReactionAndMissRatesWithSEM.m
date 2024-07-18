
close all;
clear;

bxSourceBase = 'A:\home\carlo\RC\rawData\behavioral\'; %base folder for data from bx session
bxOutputBase = 'A:\home\carlo\RC\analysis\behavioral\'; %base folder for output for bx analysis
crpFigBase = 'A:\home\carlo\RC\analysis\crpFigures\'; %base folder for output for figures

% % only used when running bx analysis on mac (also need to disable saveas for some reason
% bxSourceBase = '/Volumes/All_staff/home/carlo/rawData/behavioral/';
% bxOutputBase = '/Volumes/All_staff/home/carlo/analysis/behavioral';
% crpFigBase = '/Volumes/All_staff/home/carlo/analysis/crpFigures/';

timeBeforeCue = 2000; %defines the window around the cue presentation which will be taken for plotting
timeAfterCue = 3000;

crpTrainDayList; %load blank variable arrays
%dateIdx = input('My mice :7: | Mikes mice :6: ');

%check and make sure the figure destinations exist
    %pull out mouse and date info from bx naming
sessionName = input('which session; monomodal, dimodal, csPlusOnly: ','s');
    %name folder for session figures
sessionFigs = [crpFigBase, sessionName '_session\'];
    %name folder for xSession figures
summaryFigs = [crpFigBase, sessionName, '_summary\'];
    %check if exists, if not, make file
if exist([crpFigBase, sessionName, '_session\'], 'file') == 0
    mkdir([crpFigBase, sessionName, '_session\']);
end

    
%%
bxOutputDir  = [bxOutputBase 'grouped_bxOutputRTnMR'];

for iii = 1:length(days)
%%
gapInx = gapDayFinder(days{iii});

for ii = 1:length(days{iii})
temp_avg_licks_post_cue{ii} = [];
temp_avg_licks_pre_cue{ii} = [];
temp_avg_licks_post_cue_sem{ii} = [];
temp_avg_licks_pre_cue_sem{ii} = [];

temp_RT_across_days{ii} = [];
temp_RT_across_days_sem{ii} = [];
temp_RT_across_days_b2{ii} = [];
temp_RT_across_days_sem_b2{ii} = [];
temp_std_of_RT_across_days{ii} = [];
temp_std_of_RT_across_days_b2{ii} = [];

temp_TFT_rates{ii} = [];
temp_miss_rates{ii} = [];
temp_hit_rates{ii} = [];
temp_TFT_rates_b2{ii} = [];
temp_miss_rates_b2{ii} = [];
temp_hit_rates_b2{ii} = [];

temp_RT_across_sessions{ii} = [];
temp_RT_across_sessions_delay{ii} = [];
temp_RT_across_sessions_1000ms_delay{ii} = [];
temp_RT_across_sessions_delay_b2{ii} = [];
temp_RT_across_sessions_drift{ii} = [];
temp_RT_across_sessions_drift_b2{ii} = [];

temp_days_divider_inx{ii} = [];
temp_days_divider_inx_delay{ii} = [];
temp_days_divider_inx_1000ms_delay{ii} = [];
temp_days_divider_inx_drift{ii} = [];
temp_days_divider_inx_delay_b2{ii} = [];

temp_non_consecutive_inx{ii} = [];
temp_non_consecutive_inx_delay{ii} = [];
temp_non_consecutive_inx_delay_b2{ii} = [];
temp_non_consecutive_inx_1000ms_delay{ii} = [];
temp_non_consecutive_inx_drift{ii} = [];
temp_non_consecutive_inx_drift_b2{ii} = [];

temp_pre_cue_lick_window_avg{ii} = [];
temp_pre_cue_lick_rate_sem{ii} = [];
temp_iti_lick_window_avg{ii} = [];
temp_iti_lick_rate_sem{ii} = [];
%%
         %names output directory for bx data
    bxOutputDir  = [bxOutputBase days{iii}{ii} '_bxOutput'];
        %uses getBxData.m to find path to bx data matching thisMouse and
            %thisDate and loads the file that matches (or returns error if
            %multiple exist with the same name)
    bxData = getBxData_Carlo(bxSourceBase, days{iii}{ii});  %find the correct behavior file and loads it.
        %pulls out mouse ID and date of bx session
    [thisMouse, thisDate] = selectMouseInfo(days{iii}{ii});
        %finds each trial's start time in free floating MWorks time
    trialStart = round(double(cell2mat(bxData.tThisTrialStartTimeMs))); 
        %gets vectors for # of trials. Includes unimaged trials.
    numTrials  = length(bxData.tThisTrialStartTimeMs);
        %stores the time of the beginning of the first trial in MWorks time.
    bxStartMWorksTime  = trialStart(1);

    %use time of first frame to align licking times to start of imaging
    lickTimes = zeros(size(bxData.lickometerTimesUs{1,length(bxData.lickometerTimesUs)}));
        %determines if 2nd var is in struct[1st var]
    if isfield(bxData, 'lickometerTimesUs') 
        for kk = 1:length(bxData.lickometerTimesUs)
                %concatenates licking times
            lickTimes = [lickTimes cell2mat(bxData.lickometerTimesUs(kk))/1000]; 
        end
            %aligns lick times to the start of imaging and converts from int64 to double array
        lickTimes = double(lickTimes)-bxStartMWorksTime; 
    end
    
    %%%Collects various events during session
        %aligns start of hold to start of imaging for each cue (cue onset)
    hold_start = double(cell2mat(bxData.holdStartsMs)) - bxStartMWorksTime; 
        %duration of the "lever hold" on that trial. meaningless here except to calculate cue onset
    hold_time  = double(cell2mat(bxData.holdTimesMs)); 
        %time between cue and reward [?]
    react_time = double(cell2mat(bxData.reactTimesMs)); 
        %fixed value throughout - required hold time for reward(?)
    req_hold   = double(cell2mat(bxData.tTotalReqHoldTimeMs)); 
        %fixed value throughout - rand add time req to hold (?)
    rnd_hold   = double(cell2mat(bxData.tRandReqHoldTimeMs)); 
        %sum of req hold times - useless in this experiment
    tot_req_hold = req_hold + rnd_hold;
        %cue onset aligned to the start of imaging
    release_time = hold_start + hold_time;
        %cumHistograms are aligned to cue - time from cue onset to licking
    cuePresentation = release_time-react_time;  
        %this variable lengthens the time that the cue is on the screen
    TFT_ms = bxData.tooFastTimeMs; 
        %total interval between cue onset and reward delivery
    cue_rew_int = bxData.RewardDelayDurationMs + round(mean(react_time,'omitnan'),-1); 
    
    if exist(bxOutputDir)
            %updates saved bx data if additional data is analyized
        save(bxOutputDir, 'lickTimes', '-append');
    else
            %if no bx data has been saved, save as a .mat
        save(bxOutputDir, 'lickTimes');
    end
    %% returns a 1,1 array instead of an array 1:length(numTrials)
    %identify reward omission trials and unexpected reward trials
        %if there was no null trials
    if bxData.rewardOmissionPercent == 0 %created this if statement because in 170417_img90 there were empty cells in b_data.rewardOmissionPercent which caused the script to fail
        reward_omit_inx = [];
        %
    elseif sum(cell2mat(bxData.tRewardOmissionTrial)) == 1
         reward_omit_inx = [];
    elseif sum(cellfun(@isempty, bxData.tRewardOmissionTrial)) > 0
        bxData.tRewardOmissionTrial{cellfun(@isempty, bxData.tRewardOmissionTrial)} = int64(0);
        reward_omit_inx = find(cell2mat(bxData.tRewardOmissionTrial(1:end-1))); %exclude last trial in case it is incomplete
    else
        reward_omit_inx = find(cell2mat(bxData.tRewardOmissionTrial(1:end-1))); %exclude last trial in case it is incomplete
    end
    %find returns nonzero entries, but in one case the array is all zeros
       %so returns a 1x0 empty double row vector
    unexp_rew_inx = find(cell2mat(bxData.tDoNoStimulusChange(1:end-1))); 
    %%
    %remove block two trials from OR and UR indexes
    if bxData.doBlock2 == 1
        block2_inx = find(cell2mat(bxData.tBlock2TrialNumber(1:end-1))); %isolates nonzero trial #s
    %intersect(A,B) returns a vector of common values between A and B with no rep
        [~,inter_inx] = intersect(reward_omit_inx, block2_inx);
        reward_omit_inx(inter_inx') = []; %finds common values between A and B and removes values of A that intersect with B
        [~,inter_inx] = intersect(unexp_rew_inx, block2_inx);
        unexp_rew_inx(inter_inx') = []; %see line 91 comment
    end
    
    %isolate the time of cue onset and divide licking into trials as such
    licks_by_trial = zeros(length(cuePresentation)-1,(timeBeforeCue+timeAfterCue+1)); %dim1=trial# dim2=ms
    first_lick_by_trial = zeros(1, length(cuePresentation)-1);
    trials_where_licking_preceded_reward = zeros(1, length(cuePresentation)-1);
    trials_with_licking_soon_after_reward = zeros(1, length(cuePresentation)-1);
    for kk = 1:length(cuePresentation)-1 %look at all trials except the last one.
        %find all the licking events Xms before and Yms after cue presentation
        licks_this_window = lickTimes(lickTimes>cuePresentation(kk)-timeBeforeCue & lickTimes<cuePresentation(kk)+timeAfterCue);
        alignment_this_trial = cuePresentation(kk)-(timeBeforeCue+1); %subtract off this time so that the lick times are converted to numbers which will correspond to there index in licks_by_trial
        licks_this_window = licks_this_window - alignment_this_trial;
        licks_by_trial(kk, licks_this_window) = 1;
        if ~isempty(find(lickTimes>cuePresentation(kk) & lickTimes<cuePresentation(kk)+cue_rew_int, 1))
            trials_where_licking_preceded_reward(kk) = kk;
        end
        %if ~isempty(find(lickTimes>cue_presentation(kk)+ cue_rew_int & lickTimes < cue_presentation(kk)+ cue_rew_int +250))
         %   trials_with_licking_soon_after_reward(kk) = kk;
        %end
    end

    licks_by_trial_rewarded = licks_by_trial;
    if bxData.doBlock2 == 0
        licks_by_trial_rewarded(sort([reward_omit_inx, unexp_rew_inx]) , :) = []; %remove any unexpected rewards or reward omission trials from the rewarded trials condition. 
    elseif bxData.doBlock2 == 1
        licks_by_trial_rewarded(sort([reward_omit_inx, unexp_rew_inx, block2_inx]) , :) = [];
        licks_by_trial_block2 = licks_by_trial(block2_inx, :);
    end
    licks_by_trial_omit = licks_by_trial(reward_omit_inx, :);
    licks_by_trial_unexp = licks_by_trial(unexp_rew_inx, :);
    
    %REWARD generate cumulative histograms of licking for rewarded trials
    avg_licks_per_ms_rew = cumsum(mean(licks_by_trial_rewarded,'omitnan'));
    all_trials_lick_hist_rew = cumsum(licks_by_trial_rewarded,2);
    all_trials_lick_hist = cumsum(licks_by_trial_rewarded,2); %A=A; why does this exist
  
    %%
    %enddefine the windows for summing licks before and after cue
    lick_window = (timeBeforeCue+201:timeBeforeCue+700);
    pre_cue_window = (timeBeforeCue-700:timeBeforeCue-201);
    %determine # of licks in 500ms window before/following cue presentation for reward omission trials
    if bxData.rewardOmissionPercent > 1
        avg_licks_om_trials = mean(licks_by_trial_omit,'omitnan'); %Averages the number of licks per ms across all reward omission trials. Size is a vector with length time_before_cue+time_after_cue.
        num_om_trials = size(licks_by_trial_omit,1);
        licks_pre_window_by_trial = sum(licks_by_trial_omit(:,pre_cue_window) ,2); %sums the total # of licks within the window for each trial. Size is a vector with length = # reward omission trials.
        licks_lick_window_by_trial = sum(licks_by_trial_omit(:,lick_window) ,2);
        temp_avg_licks_post_cue_this_day = sum(avg_licks_om_trials(lick_window));  %determine the avg number (across trials) of licks which occur throughout the analysis window for this session
        temp_avg_licks_pre_cue_this_day = sum(avg_licks_om_trials(pre_cue_window));
        temp_avg_licks_post_cue = [temp_avg_licks_post_cue, temp_avg_licks_post_cue_this_day]; %concatenate number of licks in analysis window for each session
        temp_avg_licks_pre_cue = [temp_avg_licks_pre_cue, temp_avg_licks_pre_cue_this_day];
        temp_avg_licks_post_cue_sem = [temp_avg_licks_post_cue_sem, std(licks_lick_window_by_trial,'omitnan')/sqrt(num_om_trials)]; %calculate and store standard error of the mean
        temp_avg_licks_pre_cue_sem  = [temp_avg_licks_pre_cue_sem, std(licks_pre_window_by_trial,'omitnan')/sqrt(num_om_trials)];
    else 
        temp_avg_licks_post_cue = [temp_avg_licks_post_cue, NaN]; %this allows the training day # to be aligned to the pre/post cue lick rate for that day across animals
        temp_avg_licks_pre_cue = [temp_avg_licks_pre_cue, NaN];
        temp_avg_licks_post_cue_sem = [temp_avg_licks_post_cue_sem, NaN]; 
        temp_avg_licks_pre_cue_sem  = [temp_avg_licks_pre_cue_sem, NaN];
    end
    
    %determine lick rate in the middle of the iti, just before the cue, and after reward delivery. 
    avg_licks_by_trial_rewarded = mean(licks_by_trial_rewarded,'omitnan');
    
    %Determine the RTs for each trial this session
    RT_this_session{iii}{ii}=[];
    RT_this_session_noBurst{iii}{ii}=[];
    RT_this_session_shortBurst{iii}{ii}=[];

    rtTrimToggle = true; %there is an elseif below that chops out early and late reaction times. toggle this
    trialsWithHit{iii}{ii} = zeros(1,numTrials);
    for kk = 1:numTrials-1 %look at each trial excluding the last
         this_trial_post_cue_licks = find(licks_by_trial(kk,(timeBeforeCue+1:end))); %look at all lick time points starting from cue onset for each trial
         if length(this_trial_post_cue_licks) >=3 %at least three licks needed to define a burst
                             trialsWithHit{iii}{ii}(1,kk) = 1;
            for i = 1:length(this_trial_post_cue_licks)-2 %look ahead as long as three licks can still be used 
                no_burst_flag = true; 
                if  this_trial_post_cue_licks(i+2)-this_trial_post_cue_licks(i) <= 300 %if there are 3 bursts within 300 ms
                    no_burst_flag = false;
                    this_trial_RT = this_trial_post_cue_licks(i);
                    RT_this_session{iii}{ii} = [RT_this_session{iii}{ii}, this_trial_RT]; %store the RTs
                    break
                end
            end
            
            if no_burst_flag % flagged when there are many licks, but no bursts. Use first instead
               this_trial_RT = this_trial_post_cue_licks(1);
               RT_this_session{iii}{ii} = [RT_this_session{iii}{ii}, this_trial_RT];
%                RT_this_session_noBurst{iii}{ii} = [RT_this_session_noBurst{iii}{ii}, this_trial_RT];
%                RT_this_session{iii}{ii} = [RT_this_session{iii}{ii}, NaN]; 
            end

         elseif ~isempty(this_trial_post_cue_licks) %Use first lick in place of bursting if 1<=#licks<=2
             this_trial_RT = this_trial_post_cue_licks(1);
%              RT_this_session{iii}{ii} = [RT_this_session{iii}{ii}, NaN]; 
%              RT_this_session_shortBurst{iii}{ii} = [RT_this_session_shortBurst{iii}{ii}, this_trial_RT]; %store the RTs
             RT_this_session{iii}{ii} = [RT_this_session{iii}{ii}, this_trial_RT]; %store the RTs
         else %if no licks were observed
             RT_this_session{iii}{ii} = [RT_this_session{iii}{ii}, NaN]; 
         end
    end
    
    RT_this_session_raw = RT_this_session{iii}{ii};
    RT_this_session_raw(isnan(RT_this_session_raw)) = [];
    if exist([crpFigBase, [sessionName '_group_rtData\']], 'file') == 0
    mkdir([crpFigBase, [sessionName '_group_rtData\']]);
    end
    save([crpFigBase, [sessionName '_group_rtData\'] , days{iii}{ii}], 'RT_this_session_raw', 'timeBeforeCue');  

    %seperate RTs based upon the stimulus presented
    if bxData.doBlock2==1
        num_trials_b2 = length(block2_inx); %no reward
        num_trials_b1 = numTrials-1-num_trials_b2; %reward
        RT_this_session_block2{iii}{ii} = RT_this_session{iii}{ii}(block2_inx);
        no_lick_trials_b2 = sum(isnan(RT_this_session_block2{iii}{ii}));
        RT_this_session_block2{iii}{ii}(isnan(RT_this_session_block2{iii}{ii})) = []; %no reward
        RT_this_session{iii}{ii}(block2_inx) = []; %reward
    end
    no_lick_trials = sum(isnan(RT_this_session{iii}{ii})); %no licks during reward
%     RT_this_session{iii}{ii}(isnan(RT_this_session{iii}{ii})) = []; %removes reward trials without lick 
    block2_ind{iii}{ii} = bxData.tBlock2TrialNumber;
   
    %trim RT data to exclude misses and FAs  
    %%number of trials with RT>200ms divided by the number of trials in the
    %%whole session-1
    TFT_rate_this_session = length(find(RT_this_session{iii}{ii}<200))/(numTrials-1); %separates RT quicker than 200 ms
    temp_TFT_rates{ii} = [temp_TFT_rates{ii}, TFT_rate_this_session]; 
    if bxData.doBlock2 == 1
        TFT_rate_this_session_b2 = length(find(RT_this_session_block2{iii}{ii}<200))/(numTrials-1);
        temp_TFT_rates_b2{ii} = [temp_TFT_rates_b2{ii}, TFT_rate_this_session_b2];
    end
    if  TFT_ms > 0 && bxData.RewardDelayDurationMs == 0 && rtTrimToggle 
         miss_rate_this_session{iii}{ii} = length(find(RT_this_session{iii}{ii}>1500))/(numTrials-1);
         hit_rate_this_session{iii}{ii} = length(find(RT_this_session{iii}{ii}<1500))/(num_trials_b1-no_lick_trials);
         RT_this_session{iii}{ii} = RT_this_session{iii}{ii}(RT_this_session{iii}{ii}>200 & RT_this_session{iii}{ii}<1500);
        if bxData.doBlock2 == 1
         miss_rate_this_session_b2{iii}{ii} = (length(find(RT_this_session_block2{iii}{ii}>1500))+no_lick_trials_b2) / num_trials_b2;
         hit_rate_this_session_b2{iii}{ii} = (length(find(RT_this_session_block2{iii}{ii}<1500))) / (num_trials_b2-no_lick_trials);
         RT_this_session_block2{iii}{ii} = RT_this_session_block2{iii}{ii}(RT_this_session_block2{iii}{ii}>200 & RT_this_session_block2{iii}{ii}<1500);
        end
    elseif rtTrimToggle
         miss_rate_this_session{iii}{ii} = (length(find(RT_this_session{iii}{ii}>1500))+no_lick_trials) / num_trials_b1;
         hit_rate_this_session{iii}{ii} = (length(find(RT_this_session{iii}{ii}<1500))) / num_trials_b1;
         RT_this_session{iii}{ii} = RT_this_session{iii}{ii}(RT_this_session{iii}{ii}>200 & RT_this_session{iii}{ii}<1500); 
        if  bxData.doBlock2 ==1
            miss_rate_this_session_b2{iii}{ii} = (length(find(RT_this_session_block2{iii}{ii}>1500))+no_lick_trials_b2) / num_trials_b2;
            hit_rate_this_session_b2{iii}{ii} = (length(find(RT_this_session{iii}{ii}<1500))) / num_trials_b2;
            RT_this_session_block2{iii}{ii} = RT_this_session_block2{iii}{ii}(RT_this_session_block2{iii}{ii}>200 & RT_this_session_block2{iii}{ii}<1500);
        end
    end
    
    %Determine the FA and miss rates for this session. Store stats for summary statistics
    if rtTrimToggle
        temp_miss_rates{ii} = [temp_miss_rates{ii}, miss_rate_this_session{iii}{ii}];
        temp_hit_rates{ii} = [temp_hit_rates{ii}, hit_rate_this_session{iii}{ii}];
    end
    temp_RT_across_days{ii} = [temp_RT_across_days{ii}, mean(RT_this_session{iii}{ii},'omitnan')]; %values plotted in summary graphs
    temp_RT_across_days_sem{ii} = [temp_RT_across_days_sem{ii}, std(RT_this_session{iii}{ii},'omitnan')/sqrt(size(RT_this_session{iii}{ii},2))];
    temp_std_of_RT_across_days{ii} = [temp_std_of_RT_across_days{ii}, std(RT_this_session{iii}{ii},'omitnan')];
    %block2
    if bxData.doBlock2 ==1
        if rtTrimToggle
            temp_miss_rates_b2{ii} = [temp_miss_rates_b2{ii}, miss_rate_this_session_b2{iii}{ii}];
            temp_hit_rates_b2{ii} = [temp_hit_rates_b2{ii}, hit_rate_this_session_b2{iii}{ii}];
        end
        temp_RT_across_days_b2{ii} = [temp_RT_across_days_b2{ii}, mean(RT_this_session_block2{iii}{ii},'omitnan')];
        temp_RT_across_days_sem_b2{ii} = [temp_RT_across_days_sem_b2{ii}, std(RT_this_session_block2{iii}{ii},'omitnan')/sqrt(size(RT_this_session_block2{iii}{ii},2))];
        temp_std_of_RT_across_days_b2{ii} = [temp_std_of_RT_across_days_b2{ii}, std(RT_this_session_block2{iii}{ii},'omitnan')];
    end
    
   
    %store trial-by-trial RTs for across days plot 
    if bxData.rewardDelayPercent == 0
        temp_RT_across_sessions{ii} = [temp_RT_across_sessions{ii}, RT_this_session{iii}{ii}];
    elseif bxData.rewardDelayPercent == 100 && bxData.RewardDelayDurationMs == 500
        temp_RT_across_sessions_delay{ii} = [temp_RT_across_sessions_delay{ii}, RT_this_session{iii}{ii}];
        if bxData.doBlock2 == 1
            temp_RT_across_sessions_delay_b2{ii} = [temp_RT_across_sessions_delay_b2{ii}, RT_this_session_block2{iii}{ii}];
        end
    elseif bxData.rewardDelayPercent == 100 && bxData.RewardDelayDurationMs == 1000
        temp_RT_across_sessions_1000ms_delay = [temp_RT_across_sessions_1000ms_delay, RT_this_session{iii}{ii}];
    elseif bxData.rewardDelayPercent == 0 && bxData.tooFastTimeMs > 0
        temp_RT_across_sessions_delay{ii} = [temp_RT_across_sessions_delay{ii}, RT_this_session{iii}{ii}];
        if bxData.doBlock2 == 1
            temp_RT_across_sessions_delay_b2{ii} = [temp_RT_across_sessions_delay_b2{ii}, RT_this_session_block2{iii}{ii}];
        end
    end

    
    %determine trial index for trials different days and non-consecutive days
    if bxData.rewardDelayPercent == 0 && TFT_ms>0
        temp_days_divider_inx_delay{ii} = [temp_days_divider_inx_delay{ii}, length(temp_RT_across_sessions_delay{ii})+0.5];
        if ismember([ii+0.5], gapInx)
            temp_non_consecutive_inx_delay{ii} = [temp_non_consecutive_inx_delay{ii}, length(temp_RT_across_sessions_delay{ii})];
        end
        if bxData.doBlock2 ==1
            temp_days_divider_inx_delay_b2{ii} = [temp_days_divider_inx_delay_b2{ii}, length(temp_RT_across_sessions_delay_b2{ii})+0.5];
            if ismember([ii+0.5], gapInx)
                temp_non_consecutive_inx_delay_b2{ii} = [temp_non_consecutive_inx_delay_b2{ii}, length(temp_RT_across_sessions_delay_b2{ii})];
            end
        end
    elseif bxData.rewardDelayPercent == 0
        temp_days_divider_inx{ii} = [temp_days_divider_inx{ii}, length(temp_RT_across_sessions{ii})+0.5];
        if ismember([ii+0.5], gapInx)
            temp_non_consecutive_inx{ii} = [temp_non_consecutive_inx{ii}, length(temp_RT_across_sessions{ii})];
        end
    elseif bxData.rewardDelayPercent == 100 && bxData.RewardDelayDurationMs == 500
        temp_days_divider_inx_delay = [temp_days_divider_inx_delay, length(temp_RT_across_sessions_delay)+0.5];
        if ismember([ii+0.5], gapInx)
            temp_non_consecutive_inx_delay = [temp_non_consecutive_inx_delay, length(temp_RT_across_sessions_delay)];
        end
        if bxData.doBlock2 ==1
            temp_days_divider_inx_delay_b2 = [temp_days_divider_inx_delay_b2, length(temp_RT_across_sessions_delay_b2)+0.5];
            if ismember([ii+0.5], gapInx)
                temp_non_consecutive_inx_delay_b2 = [temp_non_consecutive_inx_delay_b2, length(temp_RT_across_sessions_delay_b2)];
            end
        end
    elseif bxData.rewardDelayPercent == 100 && bxData.RewardDelayDurationMs == 1000
        temp_days_divider_inx_1000ms_delay = [temp_days_divider_inx_1000ms_delay, length(temp_RT_across_sessions_1000ms_delay)+0.5];
        if ismember([ii+0.5], gapInx)
            temp_non_consecutive_inx_1000ms_delay = [temp_non_consecutive_inx_1000ms_delay, length(temp_RT_across_sessions_1000ms_delay)];
        end
    end
    
end
avg_licks_post_cue{iii} = temp_avg_licks_post_cue;
avg_licks_pre_cue{iii} = temp_avg_licks_pre_cue;
avg_licks_post_cue_sem{iii} = temp_avg_licks_post_cue_sem;
avg_licks_pre_cue_sem{iii} = temp_avg_licks_pre_cue_sem;

RT_across_days{iii} = temp_RT_across_days;
RT_across_days_sem{iii} = temp_RT_across_days_sem;
RT_across_days_b2{iii} = temp_RT_across_days_b2;
RT_across_days_sem_b2{iii} = temp_RT_across_days_sem_b2;
std_of_RT_across_days{iii} = temp_std_of_RT_across_days;
std_of_RT_across_days_b2{iii} = temp_std_of_RT_across_days_b2;

TFT_rates{iii} = temp_TFT_rates;
miss_rates{iii} = temp_miss_rates;
hit_rates{iii} = temp_hit_rates;
TFT_rates_b2{iii} = temp_TFT_rates_b2;
miss_rates_b2{iii} = temp_miss_rates_b2;
hit_rates_b2{iii} = temp_hit_rates_b2;

RT_across_sessions{iii} = temp_RT_across_sessions;
RT_across_sessions_delay{iii} = temp_RT_across_sessions_delay;
RT_across_sessions_1000ms_delay{iii} = temp_RT_across_sessions_1000ms_delay;
RT_across_sessions_delay_b2{iii} = temp_RT_across_sessions_delay_b2;
RT_across_sessions_drift{iii} = temp_RT_across_sessions_drift;
RT_across_sessions_drift_b2{iii} = temp_RT_across_sessions_drift_b2;

days_divider_inx{iii} = temp_days_divider_inx;
days_divider_inx_delay{iii} = temp_days_divider_inx_delay;
days_divider_inx_1000ms_delay{iii} = temp_days_divider_inx_1000ms_delay;
days_divider_inx_drift{iii} = temp_days_divider_inx_drift;
days_divider_inx_delay_b2{iii} = temp_days_divider_inx_delay_b2;

non_consecutive_inx{iii} = temp_non_consecutive_inx;
non_consecutive_inx_delay{iii} = temp_non_consecutive_inx_delay;
non_consecutive_inx_delay_b2{iii} = temp_non_consecutive_inx_delay_b2;
non_consecutive_inx_1000ms_delay{iii} = temp_non_consecutive_inx_1000ms_delay;
non_consecutive_inx_drift{iii} = temp_non_consecutive_inx_drift;
non_consecutive_inx_drift_b2{iii} = temp_non_consecutive_inx_drift_b2;

pre_cue_lick_window_avg{iii} = temp_pre_cue_lick_window_avg;
pre_cue_lick_rate_sem{iii} = temp_pre_cue_lick_rate_sem;
iti_lick_window_avg{iii} = temp_iti_lick_window_avg;
iti_lick_rate_sem{iii} = temp_iti_lick_rate_sem;

clearvars temp*
end
%%


for i = 1:length(RT_across_days)
nDays(i) = length(RT_across_days{1,i});
end

gRT_rew = nan(length(RT_across_days),max(nDays));
gRTsem_rew = nan(length(RT_across_days),max(nDays));
gRT_b2 = nan(length(RT_across_days),max(nDays));
gRTsem_b2 = nan(length(RT_across_days),max(nDays));

for i = 1:size(gRT_rew,1)
        gRT_rew(i,1:length(RT_across_days{1,i})) = cell2mat(RT_across_days{1,i});        
        gRT_b2(i,1:length(RT_across_days_b2{1,i})) = cell2mat(RT_across_days_b2{1,i});
end

gRTsem_rew = std(gRT_rew,[],1,'omitnan')/sqrt(size(gRT_rew,1));
gRTmean_rew = mean(gRT_rew,1,'omitnan');
gRTsem_b2 = std(gRT_b2,[],1,'omitnan')/sqrt(size(gRT_b2,1));
gRTmean_b2 = mean(gRT_b2,1,'omitnan');

tempFig=setFigure; hold on;
errorbar(gRTmean_rew, gRTsem_rew, 'k'); hold on;
%errorbar(gRT_b2, gRTsem_b2, 'r');
title(['RT across days']);
xlabel('training day');
ylabel('reaction time (s)');
xlim([0 length(gRTmean_rew)+1]);
ylim([400 2000]);
hline(767);
save([sessionFigs '\' sessionName '_reactionTimeAcrossSessions']);
saveas(tempFig,[sessionFigs '\' sessionName '_reactionTimeAcrossSessions.pdf']);

for i = 1:length(RT_across_days)
tempFig=setFigure; hold on;
errorbar(cell2mat(RT_across_days{1,i}), cell2mat(RT_across_days_sem{1,i}), 'k'); hold on;
%errorbar(gRTB2, gRTB2sem, 'r');
title(['RT across days']);
xlabel('day');
ylabel('RT (ms)');
xlim([0.8 length(RT_across_days{1,i})+0.2]);
hline(767);
end

for i = 1:length(RT_across_days)
    lastDay(i) = cell2mat(RT_across_days{1,i}(end));
    lastDay_b2(i) = cell2mat(RT_across_days_b2{1,i}(end));
end
lastDaySEM = std(lastDay,[],2,'omitnan')/sqrt(length(lastDay));
lastDaySEM_b2 = std(lastDay_b2,[],2,'omitnan')/sqrt(length(lastDay_b2));
lastDayMean=mean(lastDay,'omitnan');
lastDayMean_b2=mean(lastDay_b2,'omitnan');
tempFig=setFigure; hold on;
%individual mouse prelearning and postlearning RT
for i=1:size(days,2)
   line(1:2,[cell2mat(RT_across_days{1,i}(1,1)) cell2mat(RT_across_days{1,i}(end))],'Color',[0.5 0.5 0.5],'LineWidth',.5);
   hold on;
end
%group prelearning and postlearning RT
line(1:2,[gRTmean_rew(1,1) lastDayMean],'Color','k','LineWidth',2); % gRT(1,end) or (lastDay)
hold on;
%errorbar for prelearning mean
er=errorbar(1,[gRTmean_rew(1,1)],[gRTsem_rew(1,1)],[gRTsem_rew(1,1)],'Color','k','LineWidth',2); 
er.Color=[0 0 0];
er.LineStyle='none';
hold on;
%error bar for postlearning mean
%er=errorbar(2,[gRT(1,end)],[gRTsem(1,end)],[gRTsem(1,end)]);
er=errorbar(2,[lastDayMean],[lastDaySEM],[lastDaySEM],'Color','k','LineWidth',2);
er.Color=[0 0 0];
er.LineStyle='none';
%labels and axis size
ylabel('reaction time (s)');
ylim([400 2000]);
xlim([0.5 2.5]);
hline(767); %line for reward delivery from cue

save([sessionFigs '\' sessionName '_reactionTimesBar_lastDay']);
saveas(tempFig,[sessionFigs '\' sessionName '_reactionTimesBar_lastDay.pdf']);

% struct for data for 
[~,stats.rxTime_CSP.p,~,stats.rxTime_CSP.stats] = ttest(gRT_rew(:,1),lastDay');
stats.rxTime_CSP.day1_mean = gRTmean_rew(1,1); stats.rxTime_CSP.day1_std = gRTsem_rew(1,1);
stats.rxTime_CSP.lastDay_mean = lastDayMean; stats.rxTime_CSP.lastDay_std = lastDaySEM;

[~,stats.rxTime_CSM.p,~,stats.rxTime_CSM.stats] = ttest(gRT_b2(:,1),lastDay_b2');
stats.rxTime_CSM.day1_mean = gRTmean_b2(1,1); stats.rxTime_CSM.day1_std = gRTsem_b2(1,1);
stats.rxTime_CSM.lastDay_mean = lastDayMean_b2; stats.rxTime_CSM.lastDay_std = lastDaySEM_b2;

%

%ylim([]);
%xlim([0.8 length(RT_across_days)+0.2]);

%% miss rate
gMR = nan(length(hit_rates),max(nDays));
%gMRB2 = nan(length(miss_rates_b2{iii}{ii}),max(nDays));

for i = 1:size(gMR,1)
        gMR(i,1:length(hit_rates{1,i})) = cell2mat(hit_rates{1,i});
        %gMRB2(i,1:length(miss_rates_b2{iii}{ii}{1,i})) = miss_rates_b2{iii}{ii}{1,i};
end

gMRsem = std(gMR,[],1,'omitnan')/sqrt(size(gMR,1));
%gMRB2sem = nanstd(gMRB2,[],1);
gMRmean = mean(gMR,1,'omitnan');
%gMRB2 = nanmean(gMRB2,1);

tempFig=setFigure; hold on;
%shadedErrorBar(1:length(gMR),gMR, gMRsem, 'lineprops', 'k'); hold on;
errorbar(gMRmean, gMRsem, 'k'); hold on;
%errorbar(gMRB2, gMRB2sem, 'r');
title(['miss rate']);
xlabel('training day');
ylabel('RT (ms)');
xlim([0.8 length(gMRmean)+0.2]);
ylim([0 0.4]);
hline(0.15);
save([sessionFigs '\' sessionName '_accuracyAcrossSessions']);
saveas(tempFig,[sessionFigs '\' sessionName '_accuracyAcrossSessions.pdf']);

for i = 1:length(hit_rates)
tempFig=setFigure; hold on;
plot(cell2mat(hit_rates{1,i}), 'k'); hold on;
%errorbar(gRTB2, gRTB2sem, 'r');
title(['miss rates per animal']);
xlabel('day');
ylabel('RT (ms)');
xlim([0.8 length(hit_rates{1,i})+0.2]);
hline(.15);
ylim([0 .5]);
end
%}

for i = 1:length(RT_across_days)
    lastDayMR_rew(i) = cell2mat(hit_rates{1,i}(end));
    lastTrainingDayMR_rew(i) = cell2mat(hit_rates{1,i}(end-1));
end
lastDayMRsem_rew = std(lastDayMR_rew,[],2,'omitnan')/sqrt(length(lastDayMR_rew));
lastTrainingDayMRsem_rew = std(lastTrainingDayMR_rew,[],2,'omitnan')/sqrt(length(lastTrainingDayMR_rew));
lastDayMRmean_rew=mean(lastDayMR_rew,'omitnan');
lastTrainingDayMRmean_rew=mean(lastTrainingDayMR_rew,'omitnan');
tempFig=setFigure; hold on;
%individual animal pre- and post-learning miss rate
for i=1:size(days,2)
   line(1:2,[cell2mat(hit_rates{1,i}(1,1)) cell2mat(hit_rates{1,i}(end))],'Color',[0.5 0.5 0.5],'LineWidth',.5);
   hold on;
end
%mean prelearning and postlearning miss rate
line(1:2,[gMRmean(1,1) lastDayMRmean_rew],'Color','k','LineWidth',2); % gRT(1,end) or (lastDay)
hold on;
%prelearning sem
er=errorbar(1,[gMRmean(1,1)],[gMRsem(1,1)],[gMRsem(1,1)],'Color','k','LineWidth',2); 
er.Color=[0 0 0];
er.LineStyle='none';
hold on;
%postlearning sem
%er=errorbar(2,[gRT(1,end)],[gRTsem(1,end)],[gRTsem(1,end)]);
er=errorbar(2,[lastDayMRmean_rew],[lastDayMRsem_rew],[lastDayMRsem_rew],'Color','k','LineWidth',2);
er.Color=[0 0 0];
er.LineStyle='none';
%labels and limits
ylabel('miss rate');
ylim([0 1]);
xlim([0.75 2.25]);
title('Miss rate');
% save([sessionFigs '\' sessionName '_CS+accuracysBar_lastDay']);
% saveas(tempFig,[sessionFigs '\' sessionName '_CS+accuracysBar_lastDay.pdf']);

% structure for statistical data for miss rates
[~,stats.accuracy_CSP.p,~,stats.accuracy_CSP.stats] = ttest(lastDayMR_rew',gMR(:,1));
stats.accuracy_CSP.day1_mean = gMRmean(1,1); stats.accuracy_CSP.day1_std = gMRsem(1,1);
stats.accuracy_CSP.lastDay_mean = lastDayMRmean_rew; stats.accuracy_CSP.lastDay_std = lastDayMRsem_rew;
%

%CS-

gMRb2 = nan(length(miss_rates_b2),max(nDays));
%gMRB2 = nan(length(miss_rates_b2{iii}{ii}),max(nDays));

for i = 1:size(gMRb2,1)
        gMRb2(i,1:length(miss_rates_b2{1,i})) = cell2mat(miss_rates_b2{1,i});
        %gMRB2(i,1:length(miss_rates_b2{iii}{ii}{1,i})) = miss_rates_b2{iii}{ii}{1,i};
end

gMRsemb2 = std(gMRb2,[],1,'omitnan')/sqrt(size(gMRb2,1));
%gMRB2sem = nanstd(gMRB2,[],1);
gMRmeanb2 = mean(gMRb2,1,'omitnan');

for i = 1:length(RT_across_days)
    lastDayMR_b2(i) = cell2mat(miss_rates_b2{1,i}(end));
    lastTrainingDayMR_b2(i) = cell2mat(miss_rates_b2{1,i}(end-1));
end
lastDayMRsem_b2 = std(lastDayMR_b2,[],2,'omitnan')/sqrt(length(lastDayMR_b2));
lastTrainingDayMRsem_b2 = std(lastTrainingDayMR_b2,[],2,'omitnan')/sqrt(length(lastTrainingDayMR_b2));
lastDayMRmean_b2=mean(lastDayMR_b2,'omitnan');
lastTrainingDayMRmean_b2=mean(lastTrainingDayMR_b2,'omitnan');
%tempFig=setFigure; hold on;
%individual animal pre- and post-learning miss rate
for i=1:size(days,2)
   line(1:2,[cell2mat(miss_rates_b2{1,i}(1,1)) cell2mat(miss_rates_b2{1,i}(end))],'Color',[1 .70 .70],'LineWidth',.5);
   hold on;
end
%mean prelearning and postlearning miss rate
line(1:2,[gMRmeanb2(1,1) lastDayMRmean_b2],'Color','r','LineWidth',2); % gRT(1,end) or (lastDay)
hold on;
%prelearning sem
er=errorbar(1,[gMRmeanb2(1,1)],[gMRsemb2(1,1)],[gMRsemb2(1,1)],'Color','k','LineWidth',2); 
er.Color=[1 0 0];
er.LineStyle='none';
hold on;
%postlearning sem
%er=errorbar(2,[gRT(1,end)],[gRTsem(1,end)],[gRTsem(1,end)]);
er=errorbar(2,[lastDayMRmean_b2],[lastDayMRsem_b2],[lastDayMRsem_b2],'Color','k','LineWidth',2);
er.Color=[1 0 0];
er.LineStyle='none';
%labels and limits
ylabel('miss rate');
ylim([0 1]);
xlim([0.75 2.25]);
save([sessionFigs '\' sessionName '_accuracysBar_lastDay']);
saveas(tempFig,[sessionFigs '\' sessionName '_accuracysBar_lastDay.pdf']);

% structure for stats for CS- miss rate
[~,stats.accuracy_CSM.p,~,stats.accuracy_CSM.stats] = ttest(lastDayMR_b2',gMRb2(:,1));
stats.accuracy_CSM.day1_mean = gMRmeanb2(1,1); stats.accuracy_CSM.day1_std = gMRsemb2(1,1);
stats.accuracy_CSM.lastDay_mean = lastDayMRmean_b2; stats.accuracy_CSM.lastDay_std = lastDayMRsem_b2;
%

% proportion of trials with hit

for i = 1:size(days,2)
    for ii = 1:size(days{i},2)
posTrialsWithHit{i}{ii} = trialsWithHit{i}{ii}(~(cell2mat(block2_ind{i}{ii})));
negTrialsWithHit{i}{ii} = trialsWithHit{i}{ii}(logical(cell2mat(block2_ind{i}{ii})));
    end
end

group_preHit = [];
group_postHit = [];
for i = 1:size(days,2)
    group_preHit=[group_preHit posTrialsWithHit{i}{1}];
    group_postHit=[group_postHit posTrialsWithHit{i}{size(posTrialsWithHit{i},2)}];
end

tempFig = setFigure;
subplot(1,2,1)
pie([sum(group_preHit==1), sum(group_preHit==0)]);
legend({'with lick','without lick'},'Location','Northwest');
title('trials with licks within response window')

subplot(1,2,2)
pie([sum(group_postHit==1), sum(group_postHit==0)]);
legend({'with lick','without lick'},'Location','Northwest');
title('trials with licks within response window')

supertitle([num2str(size(group_preHit,2)) ' trials from ' num2str(size(days,2)) ' mice']);
saveas(tempFig,[sessionFigs '\' sessionName 'proportion_ofTrials_licks.pdf']);


save([sessionFigs, sessionName, '_dataAndStats.mat'],'group*','hit*','gMR*','lastDay*','hit_rates','miss_rates*','lastTraining*','sessionFigs','sessionName','days','RT*','stats');



% plot response trajectory over the first 10 days

colors = {[0/155 83/155 155/155],[22/155 99/155 99/155],[51/155 152/155 152/155]};

for d = 1:10 %first ten days
    for i = 1:size(hit_rates,2)
        avgHR_PerDay_CSP(d,i) = nan; if d<length(hit_rates{1,i}); avgHR_PerDay_CSP(d,i)=cell2mat(hit_rates{1,i}(1,d)); end
        avgHR_PerDay_CSM(d,i) = nan; if d<length(hit_rates{1,i}); avgHR_PerDay_CSM(d,i)=cell2mat(hit_rates_b2{1,i}(1,d)); end
    end
end
    tempFig=setFigure; hold on;
    xlim([1 10]); ylim([0 1]); title('average hit rate across the first 10 days of training');
    shadedErrorBar_CV(1:10,mean(avgHR_PerDay_CSP,2,'omitnan'),std(avgHR_PerDay_CSP,[],2,'omitnan')./sqrt(size(avgHR_PerDay_CSP,2)),'lineProps','k');
    shadedErrorBar_CV(1:10,mean(avgHR_PerDay_CSM,2,'omitnan'),std(avgHR_PerDay_CSM,[],2,'omitnan')./sqrt(size(avgHR_PerDay_CSP,2)),'lineProps','r');
saveas(tempFig,[sessionFigs '\' sessionName 'responseFrequency_10days.pdf']);


