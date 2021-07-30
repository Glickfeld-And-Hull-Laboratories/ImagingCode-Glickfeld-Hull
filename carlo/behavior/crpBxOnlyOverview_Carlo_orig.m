%% Modified version of the anaylysis overview for analyzing the licking behavior of cue-reward trials. 
%%SECTION ONE - assign pathnames and datasets to be analyzed/written. 
close all;
clear;

bxSourceBase = 'A:\home\carlo\rawData\behavioral\'; %base folder for data from bx session
bxOutputBase = 'A:\home\carlo\analysis\behavioral\'; %base folder for output for bx analysis
crpFigBase = 'A:\home\carlo\analysis\crpFigures\'; %base folder for output for figures

% % only used when running bx analysis on mac (also need to disable saveas for some reason
% bxSourceBase = '/Volumes/All_staff/home/carlo/rawData/behavioral/';
% bxOutputBase = '/Volumes/All_staff/home/carlo/analysis/behavioral';
% crpFigBase = '/Volumes/All_staff/home/carlo/analysis/crpFigures/';

timeBeforeCue = 2000; %defines the window around the cue presentation which will be taken for plotting
timeAfterCue = 3000;

crpTrainDayList; %load blank variable arrays
dateIdx = input('My mice :7: | Mikes mice :6: ');

%check and make sure the figure destinations exist
    %pull out mouse and date info from bx naming
[thisMouse, thisDate] = selectMouseInfo(days{1},dateIdx); 
    %name folder for session figures
sessionFigs = [crpFigBase, thisDate, '_img', thisMouse, '_session\'];
    %name folder for xSession figures
summaryFigs = [crpFigBase, thisDate, '_img', thisMouse, '_summary\'];
    %check if exists, if not, make file
if exist([crpFigBase, thisDate, '_img', thisMouse, '_session\'], 'file') == 0
    mkdir([crpFigBase, thisDate, '_img', thisMouse, '_session\']);
end
if exist([crpFigBase, thisDate, '_img', thisMouse, '_summary\'], 'file') == 0
    mkdir([crpFigBase, thisDate, '_img', thisMouse, '_summary\']);
end
    %determines the gap in days between training sessions (if a day is
        %skipped it returns (thatDate+0.5)
gapInx = gapDayFinder(days);

%% SECTION TWO - 

for ii = 1:length(days)
    %days{ii}
        %names output directory for bx data
    bxOutputDir  = [bxOutputBase days{ii} '_bxOutput'];
        %uses getBxData.m to find path to bx data matching thisMouse and
            %thisDate and loads the file that matches (or returns error if
            %multiple exist with the same name)
    bxData = getBxData_Carlo(bxSourceBase, days{ii}, dateIdx);  %find the correct behavior file and loads it.
        %pulls out mouse ID and date of bx session
    [thisMouse, thisDate] = selectMouseInfo(days{ii},dateIdx);
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
    cue_rew_int = bxData.RewardDelayDurationMs + round(mean(react_time),-1); 
    
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
    avg_licks_per_ms_rew = cumsum(mean(licks_by_trial_rewarded));
    all_trials_lick_hist_rew = cumsum(licks_by_trial_rewarded,2);
    all_trials_lick_hist = cumsum(licks_by_trial_rewarded,2); %A=A; why does this exist
    
    %REWARDED plot licking histograms for rewarded trials
    x_axis_range = [-1*timeBeforeCue:timeAfterCue]; %-2000ms to +3000ms from cue 
    if exist([sessionFigs, days{ii}, '_rew_cum_hist.fig'], 'file')==0
        figure;
        plot(x_axis_range, avg_licks_per_ms_rew, 'r', 'LineWidth', 3);  %plot average cum lick hist for t by t analysis
        ylabel('cumulative # of licks');
        xlabel(['reward delivered at ', num2str(cue_rew_int), 'ms']);
        title(['Rewarded Trials: cumulative hist of licking ', thisMouse, ' ' thisDate ' n=' num2str(size(licks_by_trial_rewarded,1))]);
        hold on;
        for kk = 1:6:size(all_trials_lick_hist_rew,1) %1:72 counting every 6th
            plot(x_axis_range, all_trials_lick_hist_rew(kk,:), 'Color', [0,0,0]+(1-(kk/size(all_trials_lick_hist_rew,1))));
        end
        vline(cue_rew_int, 'b');
        savefig([sessionFigs, days{ii}, '_rew_cum_hist']);
    end
    
    if bxData.rewardOmissionPercent > 0
        %OMISSION generate cumulative histograms of licking for omission trials
        avg_licks_per_ms = cumsum(mean(licks_by_trial_omit));
        all_trials_lick_hist = cumsum(licks_by_trial_omit,2); %/max_licks_for_norm;
        
        %OMISSION plot licking histograms for omission trials
        if exist([sessionFigs, days{ii}, '_om_cum_hist.fig'], 'file')==0
            figure;
            plot(x_axis_range, avg_licks_per_ms, 'r', 'LineWidth', 3);  %plot average cum lick hist for t by t analysis
            ylabel('cumulative # of licks');
            xlabel(['reward delivered at ', num2str(cue_rew_int), 'ms']);
            title(['Reward Omission Trials: cumulative hist of licking ', thisMouse, ' ', thisDate ' n=' num2str(size(licks_by_trial_omit,1))]);
            hold on;
            for kk = 1:size(all_trials_lick_hist,1)
                plot(x_axis_range, all_trials_lick_hist(kk,:), 'Color', [0,0,0]+(1-(kk/size(all_trials_lick_hist,1))));
            end
            savefig([sessionFigs, days{ii}, '_om_cum_hist']);
        end
    end
    
    if bxData.rewardUnexpectPercent > 0
        %UNEXPECTED reward cumulative histograms
        avg_licks_per_ms_unexp = cumsum(mean(licks_by_trial_unexp));
        all_trials_lick_hist = cumsum(licks_by_trial_unexp,2); 
        
        %UNEXPECTED plot licking histograms for UR trials
        if exist([sessionFigs, days{ii}, '_unexp_cum_hist.fig'], 'file')==0;
            figure;
            plot(x_axis_range, avg_licks_per_ms_unexp, 'r', 'LineWidth', 3);  %plot average cum lick hist for t by t analysis
            ylabel('cumulative # of licks');
            xlabel(['reward delivered at ', num2str(cue_rew_int), 'ms']);
            title(['Unexpected Reward Trials: cumulative hist of licking ', thisMouse, ' ', thisDate, ' n=' num2str(size(licks_by_trial_unexp,1))]);
            hold on;
            for kk = 1:size(all_trials_lick_hist,1)
                plot(x_axis_range, all_trials_lick_hist(kk,:), 'Color', abs([0,0,0]+(1-(kk/size(all_trials_lick_hist_rew,1)))));
            end
            savefig([sessionFigs, days{ii}, '_unexp_cum_hist']);
        end
    end
    
    if bxData.doBlock2 == 1
        %BLOCK2 cumulative histograms
        avg_licks_per_ms_b2 = cumsum(mean(licks_by_trial_block2));
        all_trials_lick_hist_b2 = cumsum(licks_by_trial_block2,2); %/max_licks_for_norm;
        
        %BLOCK2 plot licking histograms for b2 trials
        if exist([sessionFigs, days{ii}, '_block2_cum_hist.fig'], 'file')==0
            figure;
            plot(x_axis_range, avg_licks_per_ms_b2, 'r', 'LineWidth', 3);  %plot average cum lick hist for t by t analysis
            ylabel('cumulative # of licks');
            xlabel(['NO reward delivered at ', num2str(cue_rew_int), 'ms']);
            title(['CS- Trials: cumulative hist of licking ', thisMouse, ' ', thisDate, ' n=' num2str(size(licks_by_trial_unexp,1))]);
            hold on;
            for kk = 1:size(all_trials_lick_hist_b2,1)
                plot(x_axis_range, all_trials_lick_hist_b2(kk,:), 'Color', abs([0,0,0]+(1-(kk/size(all_trials_lick_hist_b2,1)))));
            end
            savefig([sessionFigs, days{ii}, '_block2_cum_hist']);
        end
    end
    
        if exist([sessionFigs, days{ii}, '_both_cum_hist.fig'], 'file')==0
            tempFig = figure;
            plot(x_axis_range, avg_licks_per_ms_b2, 'color', [0.6350, 0.0780, 0.1840], 'LineWidth', 3);  %plot average cum lick hist for t by t analysis
            hold on;
            plot(x_axis_range, avg_licks_per_ms_rew, 'color', [0, 0.4470, 0.7410], 'LineWidth', 3);
            ylabel('cumulative # of licks');
            xlabel(['NO reward delivered at ', num2str(cue_rew_int), 'ms']);
            title(['CS- Trials: cumulative hist of licking ', thisMouse, ' ', thisDate, ' n=' num2str(size(licks_by_trial_unexp,1))]);
            hold on;
            for kk = 1:6:size(all_trials_lick_hist_b2,1)
                plot(x_axis_range, all_trials_lick_hist_b2(kk,:), 'Color', abs([.2,0,0]+(0.8-(kk/size(all_trials_lick_hist_b2,1)))));
            end
            for kk = 1:6:size(all_trials_lick_hist_rew,1) %1:72 counting every 6th
                plot(x_axis_range, all_trials_lick_hist_rew(kk,:), 'Color', abs([0,0,0.2]+(0.8-(kk/size(all_trials_lick_hist_rew,1)))));
            end
            vline(767,'--k');
            vline(0,'k')
            saveas(tempFig, [sessionFigs, days{ii}, '_BothLickHist.jpeg']);
            savefig([sessionFigs, days{ii}, '_both_cum_hist']);
        end
    %% Raster plots X-axis needs adjusting to align 0 with cue presentation
    if exist([sessionFigs, days{ii}, '_rew_lick_raster.fig'], 'file')==0
        %REWARDED trials lick raster plot
        figure;
        %axis([-1*timeBeforeCue timeAfterCue 0 72]);
        plotSpikeRaster(logical(licks_by_trial_rewarded), 'PlotType', 'vertline', 'XLimForCell', [-2000 3000]); %can't align cue to 0? adding 'x_axis_range' breaks command 
        vline(timeBeforeCue+1, 'k');
        vline(timeBeforeCue+1 + cue_rew_int, 'b');
        ylabel('trial # (descending)');
        xlabel('time (ms) black=cue blue=reward');
        title(['Rewarded Trials: lick time raster img', thisMouse, ' ' thisDate ' n=' num2str(size(licks_by_trial_rewarded,1))]);
        savefig([sessionFigs, days{ii}, '_rew_lick_raster']);
    end
    
    if bxData.doBlock2 == 1
        if exist([sessionFigs, days{ii}, '_block2_lick_raster.fig'], 'file')==0
            %BLOCK2 trials lick raster plot
            figure;
            plotSpikeRaster(logical(licks_by_trial_block2), 'PlotType', 'vertline', 'XLimForCell', [-2000 3000]);
            vline(timeBeforeCue+1, 'k');
            vline(timeBeforeCue+1 + cue_rew_int, 'b');
            ylabel('trial # (descending)');
            xlabel('time (ms) black=cue blue=reward time (none)');
            title(['CS- Trials: lick time raster img', thisMouse, ' ' thisDate ' n=' num2str(size(licks_by_trial_block2,1))]);
            savefig([sessionFigs, days{ii}, '_block2_lick_raster']);
        end
        %if exist([session_fig_dir, days{ii}, '_various_licking_behaviors_by_trial.mat'], 'file')==0
            save([sessionFigs, days{ii}, '_various_licking_behaviors_by_trial.mat'],'trials_where_licking_preceded_reward') %'trials_with_licking_soon_after_reward')
        %end
    end

    LineFormatRew.Color = [0 0 0];
    LineFormatB2.Color = [1 0.2 0.2];
if exist([sessionFigs, days{ii}, '_BothLickRaster.fig'], 'file')==0
    %REWARDED and UNREWARDED trials lick raster plot
        tempFig = figure;
        %axis([-1*timeBeforeCue timeAfterCue 0 72]);
        plotSpikeRaster(logical(licks_by_trial_rewarded), 'PlotType', 'vertline', 'LineFormat', LineFormatRew); %can't align cue to 0? adding 'x_axis_range' breaks command 
        hold on;
        plotSpikeRaster(logical(licks_by_trial_block2), 'PlotType', 'vertline', 'LineFormat', LineFormatB2);
%        xlim([-2000 3000]);
        vline(timeBeforeCue+1, 'k');
        vline(timeBeforeCue+1 + cue_rew_int, 'b');
        ylabel('trial # (descending)');
        xlabel('time(ms) from trial intiation [black=cue blue=reward]');
        title([thisMouse, ' ' thisDate ' -nRew=' num2str(size(licks_by_trial_rewarded,1)) ' (blk) and nB2=' num2str(size(licks_by_trial_block2,1)) ' (red)']);
        supertitle('Rewarded and Unrewarded Trials: lick time raster ');
        savefig([sessionFigs, days{ii}, '_BothLickRaster']);
        saveas(tempFig, [sessionFigs, days{ii}, '_BothLickRaster.jpeg']);
end
    
    %%
    %enddefine the windows for summing licks before and after cue
    lick_window = (timeBeforeCue+201:timeBeforeCue+700);
    pre_cue_window = (timeBeforeCue-700:timeBeforeCue-201);
    %determine # of licks in 500ms window before/following cue presentation for reward omission trials
    if bxData.rewardOmissionPercent > 1
        avg_licks_om_trials = mean(licks_by_trial_omit); %Averages the number of licks per ms across all reward omission trials. Size is a vector with length time_before_cue+time_after_cue.
        num_om_trials = size(licks_by_trial_omit,1);
        licks_pre_window_by_trial = sum(licks_by_trial_omit(:,pre_cue_window) ,2); %sums the total # of licks within the window for each trial. Size is a vector with length = # reward omission trials.
        licks_lick_window_by_trial = sum(licks_by_trial_omit(:,lick_window) ,2);
        avg_licks_post_cue_this_day = sum(avg_licks_om_trials(lick_window));  %determine the avg number (across trials) of licks which occur throughout the analysis window for this session
        avg_licks_pre_cue_this_day = sum(avg_licks_om_trials(pre_cue_window));
        avg_licks_post_cue = [avg_licks_post_cue, avg_licks_post_cue_this_day]; %concatenate number of licks in analysis window for each session
        avg_licks_pre_cue = [avg_licks_pre_cue, avg_licks_pre_cue_this_day];
        avg_licks_post_cue_sem = [avg_licks_post_cue_sem, std(licks_lick_window_by_trial)/sqrt(num_om_trials)]; %calculate and store standard error of the mean
        avg_licks_pre_cue_sem  = [avg_licks_pre_cue_sem, std(licks_pre_window_by_trial)/sqrt(num_om_trials)];
    else 
        avg_licks_post_cue = [avg_licks_post_cue, NaN]; %this allows the training day # to be aligned to the pre/post cue lick rate for that day across animals
        avg_licks_pre_cue = [avg_licks_pre_cue, NaN];
        avg_licks_post_cue_sem = [avg_licks_post_cue_sem, NaN]; 
        avg_licks_pre_cue_sem  = [avg_licks_pre_cue_sem, NaN];
    end
    
    %determine lick rate in the middle of the iti, just before the cue, and after reward delivery. 
    avg_licks_by_trial_rewarded = mean(licks_by_trial_rewarded);
    
    %Determine the RTs for each trial this session
    RT_this_session=[];

% Below is Jake's original first lick method for reaction time. It has been
% replaced with a bursting method to capture an animals real reaction time
%{ 
    for kk = 1:num_trials-1 %look at each trial excluding the last
         this_trial_RT = find(licks_by_trial(kk,[timeBeforeCue+1:end]), 1, 'first'); %look at all time points starting from cue onset for each trial, find the first lick
          RT_this_session = [RT_this_session, this_trial_RT]; %store the RTs
          if isempty(this_trial_RT)
              RT_this_session = [RT_this_session, NaN]; %store the RTs
          end
    end
 %}
    rtTrimToggle = true; %there is an elseif below that chops out early and late reaction times. toggle this
    for kk = 1:numTrials-1 %look at each trial excluding the last
         this_trial_post_cue_licks = find(licks_by_trial(kk,(timeBeforeCue+1:end))); %look at all lick time points starting from cue onset for each trial
         if length(this_trial_post_cue_licks) >=3 %at least three licks needed to define a burst
            for i = 1:length(this_trial_post_cue_licks)-2 %look ahead as long as three licks can still be used 
                no_burst_flag = true; 
                if  this_trial_post_cue_licks(i+2)-this_trial_post_cue_licks(i) <= 500 %300ms window for bursting
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

         elseif ~isempty(this_trial_post_cue_licks) %Use first lick in place of bursting if 1<=#licks<=2
             this_trial_RT = this_trial_post_cue_licks(1);
             RT_this_session = [RT_this_session, this_trial_RT]; %store the RTs
         else %if no licks were observed
             RT_this_session = [RT_this_session, NaN]; 
         end
    end
    
    RT_this_session_raw = RT_this_session;
    RT_this_session_raw(isnan(RT_this_session_raw)) = [];
    if exist([crpFigBase, 'rtData\'], 'file') == 0
    mkdir([crpFigBase, 'rtData\']);
    end
    save([crpFigBase, 'rtData\' , days{ii}], 'RT_this_session_raw', 'timeBeforeCue');  

    %seperate RTs based upon the stimulus presented
    if bxData.doBlock2==1
        num_trials_b2 = length(block2_inx); %no reward
        num_trials_b1 = numTrials-1-num_trials_b2; %reward
        RT_this_session_block2 = RT_this_session(block2_inx);
        no_lick_trials_b2 = sum(isnan(RT_this_session_block2));
        RT_this_session_block2(isnan(RT_this_session_block2)) = []; %no reward
        RT_this_session(block2_inx) = []; %reward
    end
    no_lick_trials = sum(isnan(RT_this_session)); %no licks during reward
    RT_this_session(isnan(RT_this_session)) = []; %removes reward trials without lick 
    
   
    %trim RT data to exclude misses and FAs  
    %%number of trials with RT>200ms divided by the number of trials in the
    %%whole session-1
    TFT_rate_this_session = length(find(RT_this_session<200))/(numTrials-1); %separates RT quicker than 200 ms
    TFT_rates = [TFT_rates, TFT_rate_this_session]; 
    if bxData.doBlock2 == 1
        TFT_rate_this_session_b2 = length(find(RT_this_session_block2<200))/(numTrials-1);
        TFT_rates_b2 = [TFT_rates_b2, TFT_rate_this_session_b2];
    end
    if bxData.RewardDelayDurationMs == 500 && bxData.rewardDelayPercent ==100
        miss_rate_this_session = (length(find(RT_this_session>1000))+no_lick_trials) / numTrials-1;
        RT_this_session = RT_this_session(RT_this_session>200 & RT_this_session<1000);
        if bxData.doBlock2 ==1
            miss_rate_this_session = (length(find(RT_this_session>1000))+no_lick_trials) / num_trials_b1;
            miss_rate_this_session_b2 = (length(find(RT_this_session_block2>1000))+no_lick_trials_b2) / num_trials_b2;
            RT_this_session_block2 = RT_this_session_block2(RT_this_session_block2>200 & RT_this_session_block2<1000);
        end
    elseif bxData.RewardDelayDurationMs == 1000 && bxData.rewardDelayPercent ==100
        miss_rate_this_session = length(find(RT_this_session>1500))/(numTrials-1);
        RT_this_session = RT_this_session(RT_this_session>200 & RT_this_session<1500);
    elseif bxData.RewardDelayDurationMs ==450 && bxData.rewardDelayPercent ==100
        miss_rate_this_session = length(find(RT_this_session>1000))/(numTrials-1);
        RT_this_session = RT_this_session(RT_this_session>200 & RT_this_session<1000);
    elseif TFT_ms > 0 && bxData.RewardDelayDurationMs == 0 && rtTrimToggle %%% what is this?
%changed from RT_this_session>1200 to RT_this_session>1500
         miss_rate_this_session = length(find(RT_this_session>1500))/(numTrials-1);
%changed from 200<RT_this_session<500 to 200<RT_this_session<1500
         RT_this_session = RT_this_session(RT_this_session>200 & RT_this_session<1500);
        if bxData.doBlock2 ==1
%changed from RT_this_session>1400 to RT_this_session>1500      ?????       
            miss_rate_this_session_b2 = (length(find(RT_this_session_block2>1500))+no_lick_trials_b2) / num_trials_b2;
%changed from 200<RT_this_session<1400 to 500<RT_this_session<1500
         RT_this_session_block2 = RT_this_session_block2(RT_this_session_block2>200 & RT_this_session_block2<1500);
        end
    elseif rtTrimToggle
%changed from RT_this_session>1400 to RT_this_session>1500
         miss_rate_this_session = (length(find(RT_this_session>1500))+no_lick_trials) / num_trials_b1;
%changed from 200<RT_this_session<500 to 500<RT_this_session<1500
         RT_this_session = RT_this_session(RT_this_session>500 & RT_this_session<1500); 
        if  bxData.doBlock2 ==1
%changed from RT_this_session>1400 to RT_this_session>1500            
            miss_rate_this_session_b2 = (length(find(RT_this_session_block2>1500))+no_lick_trials_b2) / num_trials_b2;
%changed from 200<RT_this_session<500 to 500<RT_this_session<1500
            RT_this_session_block2 = RT_this_session_block2(RT_this_session_block2>500 & RT_this_session_block2<1500);
        end
    end
    
    %Determine the FA and miss rates for this session. Store stats for summary statistics
    if rtTrimToggle
        miss_rates = [miss_rates, miss_rate_this_session];
    end
    RT_across_days = [RT_across_days, mean(RT_this_session)]; %values plotted in summary graphs
    RT_across_days_sem = [RT_across_days_sem, std(RT_this_session)/sqrt(size(RT_this_session,2))];
    std_of_RT_across_days = [std_of_RT_across_days, std(RT_this_session)];
    %block2
    if bxData.doBlock2 ==1
        if rtTrimToggle
            miss_rates_b2 = [miss_rates_b2, miss_rate_this_session_b2];
        end
        RT_across_days_b2 = [RT_across_days_b2, mean(RT_this_session_block2)];
        RT_across_days_sem_b2 = [RT_across_days_sem_b2, std(RT_this_session_block2)/sqrt(size(RT_this_session_block2,2))];
        std_of_RT_across_days_b2 = [std_of_RT_across_days_b2, std(RT_this_session_block2)];
    end
    
    %plot RT within days relative to cue onset 
    if exist([sessionFigs, days{ii}, '_RT_plot.fig'], 'file')==0
        figure; plot(RT_this_session);
        title(['RT for CS+ trials within for day ' thisDate, thisMouse]);
        xlabel('trial #');
        ylabel('RT (ms)');
        savefig([sessionFigs, days{ii}, '_RT_plot']);
    end
    %plot block2 RT
    if bxData.doBlock2 ==1
        if exist([sessionFigs, days{ii}, '_RT_plot_block2.fig'], 'file')==0;
            figure; plot(RT_this_session_block2);
            title(['RT for CS- trials within for day ' thisDate, ' img' thisMouse]);
            xlabel('trial #');
            ylabel('RT (ms)');
            savefig([sessionFigs, days{ii}, '_RT_plot_block2']);
        end
    end
    
    %store trial-by-trial RTs for across days plot 
    if bxData.rewardDelayPercent == 0
        RT_across_sessions = [RT_across_sessions, RT_this_session];
    elseif bxData.rewardDelayPercent == 100 && bxData.RewardDelayDurationMs == 500
        RT_across_sessions_delay = [RT_across_sessions_delay, RT_this_session];
        if bxData.doBlock2 == 1
            RT_across_sessions_delay_b2 = [RT_across_sessions_delay_b2, RT_this_session_block2];
        end
    elseif bxData.rewardDelayPercent == 100 && bxData.RewardDelayDurationMs == 1000
        RT_across_sessions_1000ms_delay = [RT_across_sessions_1000ms_delay, RT_this_session];
    elseif bxData.rewardDelayPercent == 0 && bxData.tooFastTimeMs > 0
        RT_across_sessions_delay = [RT_across_sessions_delay, RT_this_session];
        if bxData.doBlock2 == 1
            RT_across_sessions_delay_b2 = [RT_across_sessions_delay_b2, RT_this_session_block2];
        end
    end

    
    %determine trial index for trials different days and non-consecutive days
    if bxData.rewardDelayPercent == 0 && TFT_ms>0
        days_divider_inx_delay = [days_divider_inx_delay, length(RT_across_sessions_delay)+0.5];
        if ismember([ii+0.5], gapInx)
            non_consecutive_inx_delay = [non_consecutive_inx_delay, length(RT_across_sessions_delay)];
        end
        if bxData.doBlock2 ==1
            days_divider_inx_delay_b2 = [days_divider_inx_delay_b2, length(RT_across_sessions_delay_b2)+0.5];
            if ismember([ii+0.5], gapInx)
                non_consecutive_inx_delay_b2 = [non_consecutive_inx_delay_b2, length(RT_across_sessions_delay_b2)];
            end
        end
    elseif bxData.rewardDelayPercent == 0
        days_divider_inx = [days_divider_inx, length(RT_across_sessions)+0.5];
        if ismember([ii+0.5], gapInx)
            non_consecutive_inx = [non_consecutive_inx, length(RT_across_sessions)];
        end
    elseif bxData.rewardDelayPercent == 100 && bxData.RewardDelayDurationMs == 500
        days_divider_inx_delay = [days_divider_inx_delay, length(RT_across_sessions_delay)+0.5];
        if ismember([ii+0.5], gapInx)
            non_consecutive_inx_delay = [non_consecutive_inx_delay, length(RT_across_sessions_delay)];
        end
        if bxData.doBlock2 ==1
            days_divider_inx_delay_b2 = [days_divider_inx_delay_b2, length(RT_across_sessions_delay_b2)+0.5];
            if ismember([ii+0.5], gapInx)
                non_consecutive_inx_delay_b2 = [non_consecutive_inx_delay_b2, length(RT_across_sessions_delay_b2)];
            end
        end
    elseif bxData.rewardDelayPercent == 100 && bxData.RewardDelayDurationMs == 1000
        days_divider_inx_1000ms_delay = [days_divider_inx_1000ms_delay, length(RT_across_sessions_1000ms_delay)+0.5];
        if ismember([ii+0.5], gapInx)
            non_consecutive_inx_1000ms_delay = [non_consecutive_inx_1000ms_delay, length(RT_across_sessions_1000ms_delay)];
        end
    end
    
    %----bin licking by 100ms windows relative to cue delivery----
    %determine  bin size and window size
    bin_size = 100; %number of ms to bin licking. 
    trialStart = trialStart - bxStartMWorksTime;
    min_start_to_cue = min([cuePresentation-trialStart]);
    min_cue_to_end = min([trialStart(2:end)-cuePresentation(1:end-1)]);
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
    pre_cue_window_lick = 18000;
    post_cue_window_lick = 5999;
    
    %get lick traces (1ms resolution) for rewarded and omission trials
    full_trial_licks = zeros(length(cuePresentation)-1,(pre_cue_window_lick+post_cue_window_lick+1)); %dim1=trial# dim2=ms
    for kk = 1:length(cuePresentation)-1 %look at all trials except the last one.
        %find all the licking events Xms before and Yms after cue presentation
        licks_this_window = lickTimes(lickTimes>cuePresentation(kk)-pre_cue_window_lick & lickTimes<cuePresentation(kk)+post_cue_window_lick);
        alignment_this_trial = cuePresentation(kk)-(pre_cue_window_lick+1); %subtract off this time so that the lick times are converted to numbers which will correspond to there index in licks_by_trial
        licks_this_window = licks_this_window - alignment_this_trial;
        full_trial_licks(kk, licks_this_window) = 1;
    end
    if bxData.doBlock2 == 0
        full_trial_licks_rewarded = full_trial_licks;
        full_trial_licks_rewarded(sort([reward_omit_inx, unexp_rew_inx]) , :) = []; %exclude omission and unexpected trials from rewarded trials condition
    elseif bxData.doBlock2 ==1
        full_trial_licks_rewarded = full_trial_licks;
        full_trial_licks_rewarded(sort([reward_omit_inx, unexp_rew_inx, block2_inx]) , :) = []; %exclude omission and unexpected trials from rewarded trials condition
        full_trial_licks_block2 = full_trial_licks(block2_inx, :);
    end
    full_trial_licks_omit = full_trial_licks(reward_omit_inx, :);
    full_trial_licks_unexp = full_trial_licks(unexp_rew_inx, :);
    
    %stderr of trial licks
    full_trial_licks_rewarded_ste = std(full_trial_licks_rewarded,1)./sqrt(length(full_trial_licks_rewarded));
    full_trial_licks_rewarded_bin_ste = zeros(1,(length(full_trial_licks_rewarded_ste)/bin_size));
    %bin licking by 50ms bins and convert to licks/sec
    full_trial_licks_rewarded_sum = sum(full_trial_licks_rewarded,1);
    full_trial_licks_omit_sum = sum(full_trial_licks_omit,1);
    full_trial_licks_unexp_sum = sum(full_trial_licks_unexp,1);
    full_trial_licks_rewarded_bin = zeros(1,(length(full_trial_licks_rewarded_sum)/bin_size));
    full_trial_licks_omit_bin = zeros(1,(length(full_trial_licks_omit_sum)/bin_size));
    full_trial_licks_unexp_bin = zeros(1,(length(full_trial_licks_unexp_sum)/bin_size));
    if bxData.doBlock2
        full_trial_licks_block2_sum = sum(full_trial_licks_block2,1);
        full_trial_licks_block2_ste = std(full_trial_licks_block2,1)./sqrt(length(full_trial_licks_block2));
        full_trial_licks_block2_bin = zeros(1,(length(full_trial_licks_block2_sum)/bin_size));
        full_trial_licks_block2_bin_ste = zeros(1,(length(full_trial_licks_block2_ste)/bin_size));
    end
    cue_presentation_binned = (pre_cue_window_lick/bin_size)+1;
    iii=1;
    for kk = 1:bin_size:length(full_trial_licks_rewarded_sum)
        iii=iii+1;
        full_trial_licks_rewarded_bin(iii) = sum(full_trial_licks_rewarded_sum(kk:[kk+bin_size-1]));
        full_trial_licks_rewarded_bin_ste(iii) = sum(full_trial_licks_rewarded_ste(kk:[kk+bin_size-1]));
        full_trial_licks_omit_bin(iii) = sum(full_trial_licks_omit_sum(kk:[kk+bin_size-1]));
        full_trial_licks_unexp_bin(iii) = sum(full_trial_licks_unexp_sum(kk:[kk+bin_size-1]));
        if bxData.doBlock2 ==1
            full_trial_licks_block2_bin(iii) = sum(full_trial_licks_block2_sum(kk:[kk+bin_size-1]));
            full_trial_licks_block2_bin_ste(iii) = sum(full_trial_licks_block2_ste(kk:[kk+bin_size-1]));
        end
    end
    
    %plot rewarded trials
    full_trial_licks_rewarded_bin = (full_trial_licks_rewarded_bin/size(full_trial_licks_rewarded,1))*(1000/bin_size); % convert to lick rate in Hz
    full_trial_licks_rewarded_bin_ste = (full_trial_licks_rewarded_bin_ste/size(full_trial_licks_rewarded_ste,1))*(1000/bin_size); % convert to lick rate in Hz
    full_trial_licks_block2_bin = (full_trial_licks_block2_bin/size(full_trial_licks_rewarded,1))*(1000/bin_size); % convert to lick rate in Hz
    full_trial_licks_block2_bin_ste = (full_trial_licks_block2_bin_ste/size(full_trial_licks_rewarded_ste,1))*(1000/bin_size); % convert to lick rate in Hz
    x_axis_bin = ([1:length(full_trial_licks_rewarded_bin)]-cue_presentation_binned)*(bin_size/1000);
    if exist([crpFigBase 'xAnimalHist\', days{ii}], 'file') == 0
    mkdir([crpFigBase 'xAnimalHist\', days{ii}]);
    end 
    save([crpFigBase 'xAnimalHist\', days{ii}, '\rewardHist'], 'full_trial_licks_rewarded_bin','x_axis_bin');
    if exist([sessionFigs, days{ii}, '_rew_lick_hist.fig'], 'file')==0 
        figure; bar(x_axis_bin, full_trial_licks_rewarded_bin);
        title(['Rewarded trials: lick rate per ', num2str(bin_size), 'ms', thisDate, ' img', thisMouse, 'n=', num2str(size(full_trial_licks_rewarded,1))]);
        xlabel(['time (s) relative to cue onset at ' num2str(bxData.RewardDelayDurationMs) ' ms']);
        ylabel('lick rate (Hz)');
        vline(0,'k');
        savefig([sessionFigs, days{ii}, '_rew_lick_hist']);
    end
    %plot reward omission trials
    full_trial_licks_omit_bin = (full_trial_licks_omit_bin/size(full_trial_licks_omit,1))*(1000/bin_size); % convert to lick rate in Hz
    if bxData.rewardOmissionPercent >0
        if exist([sessionFigs, days{ii}, '_om_lick_hist.fig'], 'file')==0
            figure; bar(x_axis_bin, full_trial_licks_omit_bin);
            title(['Reward omission: lick rate per ', num2str(bin_size), 'ms', ' ', thisDate, ' img', thisMouse, 'n=', num2str(size(full_trial_licks_omit,1))]);
            xlabel('time (s) relative to cue onset');
            ylabel('lick rate (Hz)');
            vline(0,'k');
            savefig([sessionFigs, days{ii}, '_om_lick_hist']);
        end
        save([crpFigBase, days{ii}, 'rewardOmissionHist'], 'full_trial_licks_omit_bin');
    end
    %plot unexpected reward trials
    full_trial_licks_unexp_bin = (full_trial_licks_unexp_bin/size(full_trial_licks_unexp,1))*(1000/bin_size); % convert to lick rate in Hz
    if bxData.rewardUnexpectPercent >0
        if exist([sessionFigs, days{ii}, '_unexp_lick_hist.fig'], 'file')==0;
            figure; bar(x_axis_bin, full_trial_licks_unexp_bin);
            title(['Unexpected reward: lick rate per ', num2str(bin_size), 'ms', thisDate,' img', thisMouse, 'n=', num2str(size(full_trial_licks_unexp,1))]);
            xlabel(['Reward delivered at ' num2str(bxData.RewardDelayDurationMs), 'ms']);
            ylabel('lick rate (Hz)');
            vline(bxData.RewardDelayDurationMs,'b');
            savefig([sessionFigs, days{ii}, 'unexpLickHist']);
        end
    end
    %plot block2 trials
    if bxData.doBlock2 > 0
        full_trial_licks_block2_bin = (full_trial_licks_block2_bin/size(full_trial_licks_block2,1))*(1000/bin_size); % convert to lick rate in Hz
        if exist([sessionFigs, days{ii}, '_block2_lick_hist.fig'], 'file')==0
            figure; bar(x_axis_bin, full_trial_licks_block2_bin);
            title(['CS-: lick rate per ', num2str(bin_size), 'ms', thisDate, ' img', thisMouse, 'n=', num2str(size(full_trial_licks_block2,1))]);
            xlabel(['time (s) relative to cue onset at ' num2str(bxData.RewardDelayDurationMs) ' ms']);
            ylabel('lick rate (Hz)');
            vline(bxData.RewardDelayDurationMs,'b');
            savefig([sessionFigs, days{ii}, '_block2_lick_hist']);
        end
        save([crpFigBase, days{ii}, 'block2Hist'], 'full_trial_licks_block2_bin');
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
    
    if exist([sessionFigs thisDate '_img' thisMouse '_stackedLickHistSEM'], 'file') == 0
% lick histogram across CS+ and CS- trials in a session
figure; bar(x_axis_bin, full_trial_licks_rewarded_bin, 'black'); 
hold on; 
erR = errorbar(x_axis_bin,full_trial_licks_rewarded_bin,(-full_trial_licks_rewarded_bin_ste(:,:)),full_trial_licks_rewarded_bin_ste);    
erR.Color = [0 0 0]; erR.LineStyle = 'none'; 
bar(x_axis_bin, full_trial_licks_block2_bin, 'red'); 
erB2 = errorbar(x_axis_bin,full_trial_licks_block2_bin,(-full_trial_licks_block2_bin_ste(:,:)),full_trial_licks_block2_bin_ste);    
erB2.Color = ('red'); erB2.LineStyle = 'none';  
title('Post-learning'); xlabel('time (s) relative to cue'); ylabel('lick rate (Hz)'); 
xline(0,'--k'); xline(770./1000,'--k'); %rew_cue_int is an int64 variable so this simple division isn't working need to convert
axis([-1 4 0 10]); xticks([-1 0 1 2 3 4]); yticks([0 5 10]);

savefig([sessionFigs thisDate '_img' thisMouse '_stackedLickHistSEM']);
    end
end 




%% save the pre/post cue data to be use in across animals summary
save([crpFigBase, days{ii}(end-4:end), 'rewardOmissionPostCueLick'], ...
    'avg_licks_post_cue', 'avg_licks_pre_cue', 'avg_licks_post_cue_sem', 'avg_licks_pre_cue_sem'); 
avg_licks_post_cue = avg_licks_post_cue(~isnan(avg_licks_post_cue)); %this allows the training day # to be aligned to the pre/post cue lick rate for that day across animals
avg_licks_pre_cue = avg_licks_pre_cue(~isnan(avg_licks_pre_cue));
avg_licks_post_cue_sem = avg_licks_post_cue_sem(~isnan(avg_licks_post_cue_sem));
avg_licks_pre_cue_sem  = avg_licks_pre_cue_sem(~isnan(avg_licks_pre_cue_sem));

%save RT std, FA rate and miss rate
if rtTrimToggle
save([crpFigBase, days{ii}(end-4:end), 'RTstd_FA_misses'], ...
    'std_of_RT_across_days', 'std_of_RT_across_days_b2','TFT_rates', 'miss_rates', 'TFT_rates_b2', 'miss_rates_b2');
end
%% Plot summary statistics 
if bxData.rewardOmissionPercent > 1
    figure; subplot(1,2,2); 
    bar(avg_licks_post_cue);
    hold on; errorbar(avg_licks_post_cue, avg_licks_post_cue_sem);
    title('Omission trials: avg # of licks in 500ms window post cue across days');
    xlabel('day');
    ylabel('# of licks');
    
    subplot(1,2,1); bar(avg_licks_pre_cue);
    hold on; errorbar(avg_licks_pre_cue, avg_licks_pre_cue_sem);
    title('Omission trials: avg # of licks in 500ms window before cue across days');
    xlabel('day');
    ylabel('# of licks');
    max_y_val = max([avg_licks_post_cue, avg_licks_pre_cue]);
    if max_y_val > 0
        subplot(1,2,1); ylim([0 ceil(max_y_val)]);
        subplot(1,2,2); ylim([0 ceil(max_y_val)]);
    end
    savefig([summaryFigs, '\', thisMouse, '_licks_pre_vs_post_cue']);
end

%% plot basic stats across days
figure; 
set(gcf, 'Position', [100, 100, 1000, 500]);
if bxData.doBlock2 ==1 
    sgtitle(['img', thisMouse, ': blue=CS+, black=CS-']);
else
    sgtitle(['img', thisMouse]);
end

subplot(2,2,1);
errorbar(RT_across_days, RT_across_days_sem, 'b'); hold on;
errorbar(RT_across_days_b2, RT_across_days_sem_b2, 'k');
title(['RT across days']);
xlabel('day');
ylabel('RT (ms)');
xlim([0.8 length(RT_across_days)+0.2]);

subplot(2,2,2); plot(std_of_RT_across_days, 'b');
title('std of RT across days');
if bxData.doBlock2 == 1
    hold on; plot(std_of_RT_across_days_b2, 'k');
end
xlabel('day');
ylabel('standard deviation');
xlim([0.8 length(std_of_RT_across_days)+0.2]);

subplot(2,2,3); plot(TFT_rates, 'b');
title('too fast lick rates as a % of all trials by day: FA = RT<500ms');
if bxData.doBlock2 == 1
    hold on; plot(TFT_rates_b2, 'k');
end
xlabel('day');
ylabel('% false alarms');
ylim([0 1]);
xlim([0.8 length(TFT_rates)+0.2]);

if rtTrimToggle
subplot(2,2,4); plot(miss_rates, 'b');
title('Miss rates as a % of all trials by day: misses = RT>1500ms');
if bxData.doBlock2 ==1
    hold on; plot(miss_rates_b2, 'k');
end
xlabel('day');
ylabel('% misses');
ylim([0 1]);
xlim([0.8 length(miss_rates)+0.2]);
end 

savefig([summaryFigs, '\', thisMouse, '_RT_RTstd_misses_TFT']);

%% 
if length(days_divider_inx) > 0
    figure; set(gcf, 'Position', [100, 100, 1000, 500]);
    plot(RT_across_sessions);
    title(['RTs across days. No delay. img' thisMouse]);
    xlabel('trial num');
    ylabel('RT(ms)');
    vline(days_divider_inx, 'k');
    if ~isempty(non_consecutive_inx)
        vline(non_consecutive_inx, 'r');
    end
    ylim([200 1000]);
    savefig([summaryFigs, '\', thisMouse, '_RT_0ms_summary']);
end

if length(RT_across_sessions_delay) > 0 & bxData.rewardDelayPercent == 100
    figure; set(gcf, 'Position', [100, 100, 1000, 500]);
    plot(RT_across_sessions_delay, 'b');
    title(['RTs across days. 500ms delay. img' thisMouse]);
    if bxData.doBlock2 == 1
        hold on; plot(RT_across_sessions_delay_b2, 'k');
        title(['RTs across days. 500ms delay. img' thisMouse, ': blue=CS+, black=CS-']);
    end
    xlabel('trial num');
    ylabel('RT(ms)');
    vline(days_divider_inx_delay, 'k');
    if ~isempty(non_consecutive_inx_delay)
        vline(non_consecutive_inx_delay, 'r');
    end
    ylim([200 1000]);
    savefig([summaryFigs, '\', thisMouse, '_RT_500ms_summary']);
end

if length(RT_across_sessions_delay) > 0 & bxData.rewardDelayPercent == 0
    figure; plot(RT_across_sessions_delay, 'b');
    title(['RTs across days. ', num2str(cue_rew_int), 'ms delay. img' thisMouse]);
    if bxData.doBlock2 == 1
        hold on; plot(RT_across_sessions_delay_b2, 'k');
        title(['RTs across days. ', num2str(cue_rew_int), 'ms delay. img' thisMouse, ': blue=CS+, black=CS-']);
        if TFT_ms > 0 & bxData.rewardDelayPercent == 0;
            title(['RTs across days. ', num2str(cue_rew_int), 'ms delay. img' thisMouse, ': blue=CS+, black=CS-']);
        end
    elseif TFT_ms > 0 & bxData.rewardDelayPercent == 0
         title(['RTs across days. ', num2str(cue_rew_int), 'ms delay. ' thisMouse]);
    end
    xlabel('trial num');
    ylabel('RT(ms)');
    vline(days_divider_inx_delay, 'k');
    if ~isempty(non_consecutive_inx_delay)
        vline(non_consecutive_inx_delay, 'r');
    end
    savefig([summaryFigs, '\', thisMouse, '_RT_1s_drift_summary']);
end

if length(RT_across_sessions_1000ms_delay) > 0
    figure; plot(RT_across_sessions_1000ms_delay);
    title(['RTs across days. 1000ms delay. img' thisMouse]);
    xlabel('trial num');
    ylabel('RT(ms)');
    vline(days_divider_inx_1000ms_delay, 'k');
    if ~isempty(non_consecutive_inx_1000ms_delay)
        vline(non_consecutive_inx_1000ms_delay, 'r');
    end
    savefig([summaryFigs, '\', thisMouse, '_RT_1000ms_summary']);
end

%% plot # of licks in iti window vs pre-cue window
figure; errorbar(iti_lick_window_avg, iti_lick_rate_sem, 'k');
hold on; 
errorbar(pre_cue_lick_window_avg, pre_cue_lick_rate_sem, 'g');
title('licking in iti (black) vs pre-cue (green)');
ylabel('avg # of licks in 500ms window');
xlabel('day #');
xlim([0.8 length(pre_cue_lick_window_avg)+0.2]);
savefig([summaryFigs, '\', thisMouse, '_iti_vs_pre-cue_licking']);


%save variables for across animals RT plot
out_dir = [crpFigBase 'rtData\'];
save([out_dir, thisMouse, 'RT_sem_across_days'], 'RT_across_days', 'RT_across_days_sem', 'RT_across_days_b2', 'RT_across_days_sem_b2');


% lick histogram across CS+ and CS- trials in a session
% figure; bar(x_axis_bin, full_trial_licks_rewarded_bin, 'black'); 
% hold on; 
% erR = errorbar(x_axis_bin,full_trial_licks_rewarded_bin,(-full_trial_licks_rewarded_bin_ste(:,:)),full_trial_licks_rewarded_bin_ste);    
% erR.Color = [0 0 0]; erR.LineStyle = 'none'; 
% bar(x_axis_bin, full_trial_licks_block2_bin, 'red'); 
% erB2 = errorbar(x_axis_bin,full_trial_licks_block2_bin,(-full_trial_licks_block2_bin_ste(:,:)),full_trial_licks_block2_bin_ste);    
% erB2.Color = ('red'); erB2.LineStyle = 'none';  
% title('Post-learning'); xlabel('time (s) relative to cue'); ylabel('lick rate (Hz)'); 
% xline(0,'--k'); xline(770./1000,'--k'); %rew_cue_int is an int64 variable so this simple division isn't working need to convert
% axis([-1 4 0 10]); xticks([-1 0 1 2 3 4]); yticks([0 5 10]);
% 
% savefig([sessionFigs '\' thisDate '_img' thisMouse '_stackedLickHistSEM']);
% f=gcf;
% saveas(f,[sessionFigs '\' thisDate '_img' thisMouse '_stackedLickHistSEM'],'epsc');
% 
% figure; bar(x_axis_bin, full_trial_licks_rewarded_bin, 'black'); 
% hold on; bar(x_axis_bin, full_trial_licks_block2_bin, 'red'); 
% title('Post-learning'); xlabel('time (s) relative to cue'); ylabel('lick rate (Hz)'); 
% xline(0,'--k'); xline(cue_rew_int./1000,'--k'); 
% axis([-1 4 0 10]); xticks([-1 0 1 2 3 4]); yticks([0 5 10]);
% 
% savefig([sessionFigs '\' thisDate '_img' thisMouse '_stackedLickHist']);
% %}