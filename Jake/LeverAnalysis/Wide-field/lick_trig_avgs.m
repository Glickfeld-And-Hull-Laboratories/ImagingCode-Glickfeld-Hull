% script for isolating lick triggered averages
clear;

pre_frames = 8;
post_frames = 4;
lever_frame = pre_frames+1;
peak_window = [pre_frames+1:pre_frames+2];

bx_source      = ['Z:\Data\WidefieldImaging\GCaMP\behavior\'];
bx_outputs_dir = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'];  
old_cd = cd; %save old cd so I can restore it later
WF_plotting_lists_of_days;

%datasets with good licking traces
days = {'151021_img29', '151022_img29', [], [], '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', '160904_img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 
colors_roi = {'b', 'r', 'g', 'k', 'c', 'm'};
colors = {'k','k', 'b','b', 'g','g', 'c','c', 'm','m', 'r','r', 'y','y', 'k', 'b', 'r'}; 
plot_symbols = {'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 's', 's', 's'};

lick_triggered_ttest = {};
lick_trig_mags = {}; 
all_ttest = {};
percent_significant{1} = [];
percent_significant_ROI{1} = [];
num_licks{1} = [];
for session_num = 1:length(days);
    if isempty(days{session_num})
        lick_triggered_ttest{session_num}= [];
        lick_trig_mags{session_num} = [];
        continue
    end
    b_data = get_bx_data(bx_source, days{session_num});  %find the correct behavior file and loads it.
    bx_out_dir  = [bx_outputs_dir days{session_num} '_bx_outputs'];
    load(bx_out_dir);
    
    imaging_start_MW_T = frame_info.imaging_start_MW_T;
    licksByFrame = licking_data.licksByFrame; 
    leverPress = cell2mat(b_data.holdStartsMs)-imaging_start_MW_T;
    holdDur = cell2mat(b_data.holdTimesMs);
    leverRelease = leverPress + holdDur;
    reqHoldTime = cell2mat(b_data.tTotalReqHoldTimeMs);
    changeOrientation = trial_outcome.change_orientation;
    
    %find and trim the trials which were not imaged. 
    trials_to_trim = [(1:frame_info.f_frame_trial_num), (frame_info.l_frame_trial_num:length(leverPress))];
    if ~isempty(trials_to_trim);
        leverPress(trials_to_trim) = [];
        holdDur(trials_to_trim) = [];
        leverRelease(trials_to_trim) = [];
        reqHoldTime(trials_to_trim) = [];
        changeOrientation(trials_to_trim) = [];
    end
    
    %collect all the lever events and index correct trials
    leverPressFrame = [];
    for iii = 1:length(leverPress)
        leverPressFrame = [leverPressFrame, frame_info.counter(leverPress(iii))];
    end
    leverReleaseFrame = [];
    for iii = 1:length(leverRelease)
        leverReleaseFrame = [leverReleaseFrame, frame_info.counter(leverRelease(iii))];
    end
    corrInx = zeros(1, length(holdDur));
    for iii = 1:length(holdDur)
        if holdDur(iii) > reqHoldTime(iii)
            corrInx(iii)=1;
        end
    end
    earlyInx = find(corrInx == 0);
    corrInx  = find(corrInx == 1);
    
    %find licking bouts not near lever events
    early_release_frame = leverReleaseFrame(earlyInx);
    corr_release_frame =  leverReleaseFrame(corrInx);
    bout_start_inx = [];
    bout_middle = [];
    for iii = 1:length(early_release_frame)-1
        this_lick_int_inx = [(early_release_frame(iii)+8):(leverPressFrame(earlyInx(iii)+1))-post_frames]; %look at the interval of licking from a given lever release to the subsequent lever press
        for iv = 1:length(this_lick_int_inx) %look at each frame in the selected interval
            if sum(licksByFrame((this_lick_int_inx(iv)-3):(this_lick_int_inx(iv)-1))) == 0     %if the three frames before this frame have no licks.
                if licksByFrame(this_lick_int_inx(iv))>0 & sum(licksByFrame([this_lick_int_inx(iv)+1:this_lick_int_inx(iv)+3]))>1  %And this frame has atleast one lick and there is at least 2 lick in the next 3 frames
                    bout_start_inx = [bout_start_inx, this_lick_int_inx(iv)];
                end
            end
        end
    end
    for iii = 1:length(corr_release_frame)-1
        this_lick_int_inx = [(corr_release_frame(iii)+12):(leverPressFrame(corrInx(iii)+1))-post_frames]; %look at the interval of licking from a given lever release to the subsequent lever press
        for iv = 1:length(this_lick_int_inx) %look at each frame in the selected interval
            if sum(licksByFrame((this_lick_int_inx(iv)-3):(this_lick_int_inx(iv)-1))) == 0     %if the three frames before this frame have no licks.
                if licksByFrame(this_lick_int_inx(iv))>0 & sum(licksByFrame([this_lick_int_inx(iv)+1:this_lick_int_inx(iv)+3]))>1 %And this frame has atleast one lick and there is at least 2 lick in the next 3 frames
                    bout_start_inx = [bout_start_inx, this_lick_int_inx(iv)];
                end
            end
        end
    end
    
    %if there are valid licks 
    if ~isempty(bout_start_inx)
        %extract traces from each bout
        lick_start_trig = [];
        for bout_num = 1:length(bout_start_inx)
            for iv = 1:length(LS_ROIs{session_num}) % 1:size(tc_dfoverf,1)
                lick_start_trig(bout_num,iv,:) = tc_dfoverf(LS_ROIs{session_num}(iv), [(bout_start_inx(bout_num)-pre_frames):(bout_start_inx(bout_num)+post_frames)]);   %dim1 bout#  dim2 ROI #   dim3 f values
            end
        end
        
        %shift each TC so 0 is the mean of...
        %shift = mean(lick_start_trig(:,:,[5:6]),3);
        shift = mean(lick_start_trig(:,:,[1:6]),3);
        shift = repmat(shift, 1, 1, post_frames+pre_frames+1);
        lick_start_trig = lick_start_trig-shift;
        
        %optional - plot individual trials
%         figure; 
%         for ROI_num = 1:size(lick_start_trig,2)
%             subplot(1, size(lick_start_trig,2), ROI_num);
%             this_baseline = squeeze(lick_start_trig(:,ROI_num,2));
%             ttest_results = [];
%             for bout_num = 1:size(lick_start_trig,1)
%                 plot([-pre_frames:post_frames], squeeze(lick_start_trig(bout_num,ROI_num,:))); hold on;
%                 ttest_results = [ttest_results, ttest( squeeze(mean(lick_start_trig(:,ROI_num,[1:5]),3)), squeeze(mean(lick_start_trig(bout_num,ROI_num,[lever_frame-1:lever_frame+1]))), 'Tail', 'left'   )];
%             end
%             title(['ROI#', num2str(ROI_num), '. Fraction of trials with a significant increase ', num2str(sum(ttest_results)/length(ttest_results))]);
%             ylabel('df/f');
%             xlabel('frames relative to lick');
%             all_ttest{session_num}(ROI_num) = sum(ttest_results)/length(ttest_results);
%         end
%         suptitle([days{session_num}, ' 3 lick bout traces. n=', num2str(size(lick_start_trig,1))]);
%         %savefig(['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\individual lick bout TCs\', days{session_num}, '_bout']);
        
        %get mean and sem across trials for each ROI
        lick_start_sem = squeeze(std(lick_start_trig,1)/sqrt(size(lick_start_trig,1)));
        lick_start_trig_avg = squeeze(mean(lick_start_trig)); %averaging across trials
        if length(LS_ROIs{session_num}) ==1
            lick_start_trig_avg = lick_start_trig_avg';
            lick_start_sem = lick_start_sem';
        end
        
        %find and store the peak values (averaged across trials)
        if size(lick_start_trig,1) >= 4
            lick_trig_mags{session_num} = max(lick_start_trig_avg(:,peak_window),[],2)';
        else
            lick_trig_mags{session_num} = [];
        end
        
        %plot the mean TCs for each ROI
        figure; hold on;
        for iii = 1:length(LS_ROIs{session_num})
            errorbar([-pre_frames:post_frames], lick_start_trig_avg(iii,:), lick_start_sem(iii,:));
            %plot([-4:4], tc_dfoverf([(bout_start_inx(iii)-4):(bout_start_inx(iii)+4)]))
            title([days(session_num), ['df/f response to licking bout start n=' num2str(length(bout_start_inx))], 'minimum 3 lick per bout', 'no licks in 3 frames before zero']);
            xlabel('frame number relative to licking bout start');
            ylabel('df/f');
        end
        %savefig(['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\', days{session_num},'_bout']);
        
        %measure for significant response. Peak > 2 std away from baseline
        %find periods in the ITI without licking 
        ITI_inx = ones(1,length(tc_dfoverf));
        for lever_event = length(leverReleaseFrame)
            ITI_inx(leverPressFrame(lever_event)-3:leverReleaseFrame(lever_event)+8) = 0; 
        end
        %sometimes there are missing frames at the end of licksbyframe.
        %need to lengthen it and correct ITI_inx
        if length(ITI_inx) > length(licksByFrame)
            ITI_inx(length(licksByFrame)+1:end) = 0;
            licksByFrame(length(licksByFrame)+1:length(ITI_inx)) = 0;
        end
        %isolate the licks which occur in the ITI and use them to exclude frames from which to measure the baseline df/f 
        ITI_licks = find(licksByFrame.*ITI_inx);  
        for lick_num = 1:length(ITI_licks)
            if ITI_licks(lick_num)<3
                ITI_inx([1:4]) =0;
                continue
            end
            ITI_inx(ITI_licks(lick_num)-2:ITI_licks(lick_num)+2) = 0;
        end
        baseline_frames = tc_dfoverf([LS_ROIs{session_num}], [find(ITI_inx)]); %exclude non-LS ROIs
        baseline_frames_std = std(baseline_frames,[],2)';
        %determine what percent of trials exceed 2 stdevs first by ROI
        tbyt_peak_mag = max(lick_start_trig(:,:,peak_window),[],3);
        this_session_sig = [];
        for ROI_num = 1:size(tbyt_peak_mag,2);
            this_session_sig = [this_session_sig, length(find(tbyt_peak_mag(:,ROI_num)>2*baseline_frames_std(ROI_num)))/size(tbyt_peak_mag,1)];
        end
        percent_significant_ROI{session_num} = this_session_sig;
        %then averaged across ROIs
        lick_start_trig_avg = squeeze(mean(lick_start_trig,2));
        tbyt_peak_mag_avg = max(lick_start_trig_avg(:,peak_window),[],2);
        percent_significant{session_num} = length(find(tbyt_peak_mag_avg>2*mean(baseline_frames_std)))/length(tbyt_peak_mag_avg);
        num_licks{session_num} = length(tbyt_peak_mag_avg);
        
        %identify which trials have a significant response
        sig_lick_resp_mat = zeros(size(tbyt_peak_mag));
        for ROI_num = 1:size(tbyt_peak_mag,2)
            sig_lick_resp_mat([find(tbyt_peak_mag(:,ROI_num)>2*baseline_frames_std(ROI_num))], ROI_num) = 1; %find all trials with significant response and log them in this matrix on an ROI basis
        end
        
        %make a color coded plot: trial by trial TCs of sig vs non-sig trials
        figure; 
        for ROI_num = 1:size(tbyt_peak_mag,2)
            subplot(1,size(tbyt_peak_mag,2),ROI_num); hold on;
            title(['ROI #', num2str(ROI_num)]);
            ylabel('df/f'); xlabel('frames relative to lick');
            for trial_num = 1:size(tbyt_peak_mag,1)
                if sig_lick_resp_mat(trial_num, ROI_num) == 1
                    plot([-pre_frames:post_frames], squeeze(lick_start_trig(trial_num, ROI_num, :)), 'b');
                elseif sig_lick_resp_mat(trial_num, ROI_num) == 0
                    plot([-pre_frames:post_frames], squeeze(lick_start_trig(trial_num, ROI_num, :)), 'k');
                end
            end
        end
        suptitle([days(session_num), 'trial by trial lick triggered TCs: blue=significant black=non-sig (>2 St.Dev)']);
        savefig(['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\individual lick bout TCs\', days{session_num}, '_bout']);
    end
end

plotting_sig_sessions = cell2mat(percent_significant);
figure;
plot(ones(1, length(plotting_sig_sessions)), plotting_sig_sessions, 'o');
ylim([0 1]); ylabel('percent of lick bouts with a significant positive df/f response');
title('trial peak df/f > 2 standard deviations above baseline. Averaged across ROIs');
%savefig(['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\individual lick bout TCs\percent_sig_trials_2_StDev']);
%save(['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\individual lick bout TCs\all_ttest'], 'all_ttest', 'all_ttest_mean');
%save(['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\lick_trig_mags'], 'lick_trig_mags');
%save(['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\lick_triggered_ttest_bout'], 'lick_triggered_ttest');

%plot each ROI
figure;
for session_num = 1:length(percent_significant_ROI)
    if isempty(percent_significant_ROI{session_num})
        continue
    end
    plot(ones(1, length(percent_significant_ROI{session_num})), percent_significant_ROI{session_num}, strcat(plot_symbols{session_num}, colors{session_num})); hold on;
end
%plot(ones(1, length(plotting_sig_sessions)), plotting_sig_sessions, 'o');
ylim([0 1]); ylabel('percent of lick bouts with a significant positive df/f response');
title('trial peak df/f > 2 standard deviations above baseline.');
%savefig(['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\individual lick bout TCs\percent_sig_trials_2_StDev_ROI']);



