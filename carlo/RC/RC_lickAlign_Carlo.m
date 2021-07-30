%%
% 2102141_img1081
% 2102150_img1081
%
%% load in expt data and neural data for each mouse
clear; close all;
analysis_out = 'A:\home\carlo\analysis\2P\';
bdata_source = 'A:\home\carlo\rawData\behavioral\';

crpTrainDayList; %clearvars -except an* b* days
ii = length(days); dateIdx=7;

[mouse, date] = selectMouseInfo(days{ii},dateIdx); 
%pulls out mouse ID and date of bx label

%% SJ lick align
%RCExptListCarlo_Block
%RCExptListCarlo_Inter
expType=input('monomodal(1) or dimodal(2) session: ');
seshType=input('interleaved (1) or blocked (2) session: ');
if seshType==1 && expType == 1
RCExptListCarlo_Inter;
elseif seshType==2 && expType == 1
RCExptListCarlo_Block;
elseif seshType==1 && expType == 2
RCExptListCarlo_VA_Inter;
end

pairings = loadRCList_Carlo(expt,seshType);
j=pairings{1,1};
    id = j;
    %exp_subset = j;
    iexp = j;  isesh = pairings{3,1};
mouse = strtrim(expt{1,iexp}.mouse{1,isesh});
date = expt{1,iexp}.date{1,isesh};
run = expt{1,iexp}.run{1,isesh};
    sessions = [date '_img' mouse];
    fprintf([date ' ' mouse ' ' run '\n']);
    img_fn = [date '_img' mouse '\getTC_' run '\'];
    input = getBxData_Carlo(bdata_source, sessions, 7);
    load(fullfile(analysis_out,img_fn, '_targetAlign.mat'));
    
    rewDelay_frames =  round((mean(celleqel2mat_padded(input.reactTimesMs))/1000).*frameRateHz);% there's 700ms between the cue and the reward delivery !!!!! if you change tooFastMs in the varibles in MWorks, this needs to be changed
    
    cTargetOn = input.cTargetOn;
    if iscell(cTargetOn) % if it is a cell, it means cTargetOn wasn't being over written in extract TC. If it's not a cell, it's already over written in extract TC. can be used directly
        cTargetOn = celleqel2mat_padded(input.cTargetOn);
        cTargetOn(1) = nan; % first trial doesn't have reward 
    end
 
    nTrials = size(cTargetOn(1:end),2);
    nIC = size(targetAlign_events,2); %you get targetalign_events from neural align, this is the neural data (in spikes) aligned to cue
    lickCueAlign =  nan(prewin_frames+postwin_frames,nTrials);
    lickCounterVals = cell(1,nTrials);
    nFrames = input.counterValues{nTrials}(end);
    lickDelay_frames =  round(0.1.*frameRateHz);
    lickSearch_frames =  round(0.3.*frameRateHz);
    
    postLick_frames = round(0.5.*frameRateHz);
    lickSearchRange = prewin_frames+lickDelay_frames:size(lickCueAlign,1)-lickSearch_frames;
    preRew_lickSearchRange_700ms = prewin_frames+lickDelay_frames:prewin_frames+lickDelay_frames+rewDelay_frames;
    lickBurstStart = nan(1,nTrials);
    postRew_lickBurstStart = nan(1,nTrials);
    lickBurstHz_all = nan(1,nTrials);
    preRew_lickBurstHz = nan(1,nTrials);
    postRew_lickBurstHz = nan(1,nTrials);
    postRew_lickAlignEvents = nan(2.*postLick_frames,nIC,nTrials);
    postRew_lickAlign = nan(2.*postLick_frames,nTrials);
    lastPreRew_lickAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    firstPostRew_lickAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    lastPreRew_lickAlign = nan(3.*postLick_frames,nTrials);
    firstPostRew_lickAlign = nan(3.*postLick_frames,nTrials);
    rewAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    preRewTrials = [];
    postRewTrials = [];
    postRew_lickSearchRange = prewin_frames+rewDelay_frames+lickDelay_frames:size(lickCueAlign,1)-lickSearch_frames-postLick_frames;% need to '-lickSearch_frames-postLick_frames' because later the inds needs to + lickSearch_frames or +postLick_frames, this is just for the following inds to be within matrix dimentions
    lastPreRewLickFrame = nan(1,nTrials);
    firstPostRewLickFrame = nan(1,nTrials);
    
    tTInx = 1;
    for itrial = 1:nTrials
        if ~isnan(cTargetOn(itrial))
            if cTargetOn(itrial)+postwin_frames-1 < nFrames
                lickTimes = input.lickometerTimesUs{itrial}; 
                counterTimes = input.counterTimesUs{itrial};
                counterVals = input.counterValues{itrial};
                lickCounterVals{itrial} = zeros(size(lickTimes));
                lickTC{itrial} = zeros(size(counterVals));% find how many licks for each frame
             %%%this for loop pulls out every lickTime that falls between the start of two frames
               %and if licks are found it outputs the frame licks occur to lickCounterVals
                for icount = 1:length(counterTimes)-1
                    ind = find(lickTimes>counterTimes(icount) & lickTimes<counterTimes(icount+1));
                    if ~isempty(ind)
                        lickCounterVals{itrial}(1,ind) = icount; % find which counter in this trial has licks
                    end
                end
             %%%this for loop counts down the number of frames in a trial and outputs a logical
               %outlining whether or not a lick occured between two counters (or frames)
                for ival = 1:length(counterTimes)
                    ind = find(lickCounterVals{itrial} == ival);
                    if ~isempty(ind)
                        lickTC{itrial}(1,ival) = length(ind);% if the mice doesn't lick during that counter than it's zero, if it does, it's one. length(ind) should always be 1.
                    end
                end
             %%%IF the last frame of the trial (from the start of the ITI preceeding the trial to the start of the ITI following the trial)
               %minus the frame where the cue is shown IS greater than the frame difference between   
                if input.counterValues{itrial}(end)-cTargetOn(itrial) > postwin_frames
                    %lickcuealign aligns licking of each frame to cue
                    lickCueAlign(:,itrial) = lickTC{itrial}(1,cTargetOn(itrial)-prewin_frames-counterVals(1):cTargetOn(itrial)+postwin_frames-1-counterVals(1));
                end
                ind = intersect(lickSearchRange,find(lickCueAlign(:,itrial)));
             %%%this for loop find frames where more than 3 licks occurred within a fixed window, and index that frame as the start of burst lick
                for i = 1:length(ind)
                    ilick = ind(i);
                    if sum(lickCueAlign(ilick:ilick+lickSearch_frames-1,itrial),1) >= 3 % more than 3 bursts within about 300ms
                        lickBurstStart(:,itrial) = ilick;% when licking burst happens
                        break
                    end
                end
                ind_pre = intersect(preRew_lickSearchRange_700ms,find(lickCueAlign(:,itrial))); %finds every instance of a lick that occurs within the search window
                ind_post = intersect(postRew_lickSearchRange,find(lickCueAlign(:,itrial))); %finds every instance of a lick that occurs [after the reward delivery]
                ind_all = intersect(lickSearchRange,find(lickCueAlign(:,itrial)));
                preRew_lickBurstHz(:,itrial) = length(ind_pre)./(mean(celleqel2mat_padded(input.reactTimesMs))/1000);
                postRew_lickBurstHz(:,itrial) = length(ind_post)./(length(postRew_lickSearchRange)./frameRateHz);
                lickBurstHz_all(:,itrial) = length(ind_all)./(length(lickSearchRange)./frameRateHz);
                ind = intersect(postRew_lickSearchRange,find(lickCueAlign(:,itrial)));
             %%%similar idea as the above for loop, but instead aligns neural data as POST-REWARD lick events 
                for i = 1:length(ind)
                    ilick = ind(i);
                    if sum(lickCueAlign(ilick:ilick+lickSearch_frames-1,itrial),1) >= 3
                        postRew_lickBurstStart(:,itrial) = ilick; %array of every lick within post-reward window if 3+ licks were recorded for that trial
                        postRew_lickAlignEvents(:,:,itrial) = targetAlign_events(ilick-postLick_frames:ilick+postLick_frames-1,:,itrial);% align neural data to first lick after reward
                        postRew_lickAlign(:,itrial) = lickCueAlign(ilick-postLick_frames:ilick+postLick_frames-1,itrial);% align licking data to first lick
                        break
                    end
                end
                ind_pre = find(lickCueAlign(prewin_frames:prewin_frames+rewDelay_frames,itrial),1,'last');
             %%%if there are licks recorded between cue and reward delivery (770 ms), align neural data to those licks 
                if ~isempty(ind_pre)
                    lastPreRewLickFrame(1,itrial) = ind_pre;
                    preRewTrials = [preRewTrials itrial];
                    % align neural and licking data to the last lick
                    lastPreRew_lickAlignEvents(:,:,itrial) = targetAlign_events(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,:,itrial);
                    lastPreRew_lickAlign(:,itrial) = lickCueAlign(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,itrial);
                end
                ind_post = find(lickCueAlign(prewin_frames+rewDelay_frames:prewin_frames+postwin_frames-postLick_frames-postLick_frames,itrial),1,'first');
             %%%   
                if ~isempty(ind_post)
                    firstPostRewLickFrame(1,itrial) = ind_post;
                    postRewTrials = [postRewTrials itrial];
                    % align neural and licking data to the first lick
                    firstPostRew_lickAlignEvents(:,:,itrial) = targetAlign_events(ind_post-postLick_frames+prewin_frames+rewDelay_frames:ind_post+postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,:,itrial);
                    firstPostRew_lickAlign(:,itrial) = lickCueAlign(ind_post-postLick_frames+prewin_frames+rewDelay_frames:ind_post+postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,itrial);
                end
                rewAlignEvents(:,:,itrial) =targetAlign_events(-postLick_frames+prewin_frames+rewDelay_frames:postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,:,itrial);
            
            tTInx = tTInx+1;
            end
        end
    end
    
    figure; %align licking data to cue
    shadedErrorBar(tt, nanmean(lickCueAlign,2).*(1000./frameRateHz), (nanstd(lickCueAlign,[],2)./sqrt(unique(sum(~isnan(lickCueAlign),2))).*(1000./frameRateHz)));
    hold on;
    scatter((lickBurstStart-prewin_frames).*(1000./frameRateHz), 10.*ones(size(lickBurstStart)), 'x');
    xlabel('Time from cue');
    ylabel('Lick rate (Hz)');
    title([date ' ' mouse '' run]);
    vline(0,'k');
    vline(700,'r');
    savefig(fullfile(analysis_out,img_fn, '_cueAlign_lickHz.fig'));
    
    nIC = size(targetAlign_events,2);
    if sum(~isnan(lickBurstStart))>6
        [sortlick sortlick_ind] = sort(lickBurstStart,'ascend');
        nburst = sum(~isnan(lickBurstStart));
        nnan = sum(isnan(lickBurstStart));
        ind_early_bst = sortlick_ind(1:floor(nburst/4));
        ind_late_bst = sortlick_ind(nburst-floor(nburst/4)+1:end-nnan);
        early_bst_time = mean((lickBurstStart(:,ind_early_bst)-prewin_frames).*(1000./frameRateHz),2);
        late_bst_time = mean((lickBurstStart(:,ind_late_bst)-prewin_frames).*(1000./frameRateHz),2);
    else
        ind_early_bst = [];
        ind_late_bst = [];
        early_bst_time = [];
        late_bst_time = [];
    end
   
    if sum(~isnan(lickBurstStart))>6
        figure; %plot neural data of trials of early vs. late bursts
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_early_bst),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_early_bst),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'k');
        hold on;
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_late_bst),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_late_bst),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'b');
        scatter((lickBurstStart(:,ind_early_bst)-prewin_frames).*(1000./frameRateHz),-0.5.*ones(1,length(ind_early_bst)),'xk');
        scatter((lickBurstStart(:,ind_late_bst)-prewin_frames).*(1000./frameRateHz),-0.5.*ones(1,length(ind_late_bst)),'xb');
        xlabel('Time from cue');
        ylabel('Spike rate (Hz)');
        ylim([-1 inf]);
        title([mouse ' ' date ' - lick bursts: early (n = ' num2str(length(ind_early_bst)) '); late (n = ' num2str(length(ind_late_bst)) ')']);
        savefig(fullfile(analysis_out,img_fn, '_cueAlignSpiking_byLickTime.fig'));
        hold off;
    end
    
    pct_precue_burst = length(find((lickBurstStart-prewin_frames).*(1000./frameRateHz)<600))./size(lickBurstStart,2);
    tl = (1-postLick_frames:postLick_frames).*(1000./frameRateHz);
    
    figure;
    subplot(2,1,1); % align neural activity to lick onset, only burst trials are included 
    shadedErrorBar(tl, nanmean(nanmean(postRew_lickAlignEvents,3),2).*(1000./frameRateHz), (nanstd(nanmean(postRew_lickAlignEvents,3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    hold on;
    xlabel('Time from lick');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    title([num2str(sum(~isnan(squeeze(postRew_lickAlignEvents(1,1,:))))) ' trials post-reward lick bursts (black)']);
    subplot(2,1,2); % seperate early burst trials vs. late burst trials
    [sortlick sortlick_ind] = sort(postRew_lickBurstStart,'ascend'); %sorts trial-based instances of lick and outputs [sorted instances, index of original position before sorting]
    nburst = sum(~isnan(postRew_lickBurstStart)); %total trials with lick burst
    nnan = sum(isnan(postRew_lickBurstStart)); %total trials without lick burst (some may be CS-; some may have no burst but still be CS+)
    ind_prerew_early_bst = sortlick_ind(1:floor(nburst/4)); %first 1/4 chosen as early trials
    ind_prerew_late_bst = sortlick_ind(nburst-floor(nburst/4)+1:end-nnan); %last 1/4 chosen as late trials 
    early_bst_time = nanmean((postRew_lickBurstStart(:,ind_prerew_early_bst)-prewin_frames-rewDelay_frames).*(1000./frameRateHz),2);
    late_bst_time = nanmean((postRew_lickBurstStart(:,ind_prerew_late_bst)-prewin_frames-rewDelay_frames).*(1000./frameRateHz),2);
    shadedErrorBar(tl, nanmean(nanmean(postRew_lickAlignEvents(:,:,ind_prerew_early_bst),3),2).*(1000./frameRateHz), (nanstd(nanmean(postRew_lickAlignEvents(:,:,ind_prerew_early_bst),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    hold on;
    shadedErrorBar(tl, nanmean(nanmean(postRew_lickAlignEvents(:,:,ind_prerew_late_bst),3),2).*(1000./frameRateHz), (nanstd(nanmean(postRew_lickAlignEvents(:,:,ind_prerew_late_bst),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'b');
    xlabel('Time from lick');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    title([num2str(floor(nburst/4)) ' earliest burst lick trials [blk]: avg = ' num2str(chop(early_bst_time,3)) ' ms; ' num2str(floor(nburst/4)) ' latest burst lick trials [blu]: avg = ' num2str(chop(late_bst_time,3)) ' ms)'])
    savefig(fullfile(analysis_out,img_fn, '_postRew_lickBurstAlignSpikingEvents.fig'))

    figure;
    subplot(2,1,1); % align licking data to first lick
    shadedErrorBar(tl, nanmean(postRew_lickAlign,2).*(1000./frameRateHz), (nanstd(postRew_lickAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(postRew_lickAlign),2))),'k');
    hold on;
    xlabel('Time from lick');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    title('Neural data aligned to first lick post-reward per trial')
    subplot(2,1,2); % plot early burst trials and late burst trials separately
    shadedErrorBar(tl, nanmean(postRew_lickAlign(:,ind_prerew_early_bst),2).*(1000./frameRateHz), (nanstd(postRew_lickAlign(:,ind_prerew_early_bst),[],2).*(1000./frameRateHz))./sqrt(length(ind_prerew_early_bst)),'k');
    hold on;
    shadedErrorBar(tl, nanmean(postRew_lickAlign(:,ind_prerew_late_bst),2).*(1000./frameRateHz), (nanstd(postRew_lickAlign(:,ind_prerew_late_bst),[],2).*(1000./frameRateHz))./sqrt(length(ind_prerew_late_bst)),'b');
    xlabel('Time from lick');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    title('Early burst trials [blk] and late burst trials [blu]')
    supertitle([mouse ' ' date '- post reward lick burst aligned spiking']);
    savefig(fullfile(analysis_out,img_fn, '_postRew_lickBurstAlignSpiking.fig'));
    hold off;
    
    figure;
    colour = {'k', 'b', 'r', 'm'};
    for i = 1:4
        plot(tt,(cumsum(nansum(lickCueAlign(:,1+((i-1)*floor(nTrials/4)):floor(nTrials/4)+((i-1)*floor(nTrials/4))),2))),colour{1,i}, 'LineWidth',2.0);
        hold on;
    end
    title([mouse ' ' date '- cumulative licking by quarter session']);
    xlabel('Time from cue');
    ylabel('Cumulative Licks');
    legend({'first quar', 'second quar', 'third quar', 'fourth quar'},'Location','northwest');
    hold off;
    savefig(fullfile(analysis_out,img_fn, '_cumulativeLicking.fig'));
    
    tl_rew = (1-postLick_frames:(postLick_frames*2)).*(1000./frameRateHz);
    figure; % align neural and licking data to the first and last lick, respectively
    subplot(2,2,1);
    shadedErrorBar(tl_rew, nanmean(nanmean(lastPreRew_lickAlignEvents,3),2).*(1000./frameRateHz), (nanstd(nanmean(lastPreRew_lickAlignEvents,3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    hold on;
    title('Last lick before reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    subplot(2,2,3);
    shadedErrorBar(tl_rew, nanmean(lastPreRew_lickAlign,2).*(1000./frameRateHz), (nanstd(lastPreRew_lickAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(lastPreRew_lickAlign)))),'k');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 30]);
    title(['Pre-rew: ' num2str(length(preRewTrials))]);
    subplot(2,2,2);
    shadedErrorBar(tl_rew, nanmean(nanmean(firstPostRew_lickAlignEvents,3),2).*(1000./frameRateHz), (nanstd(nanmean(firstPostRew_lickAlignEvents,3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    title('First lick after reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    subplot(2,2,4);
    shadedErrorBar(tl_rew, nanmean(firstPostRew_lickAlign,2).*(1000./frameRateHz), (nanstd(firstPostRew_lickAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(firstPostRew_lickAlign)))),'k');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    title(['Post-Rew: ' num2str(length(postRewTrials))]);
    supertitle([mouse ' ' date '- Licks relative to reward']);
    hold off;
    savefig(fullfile(analysis_out,img_fn, '_lastVsFirstLick.fig'));
    
    [sortlickHz sortlickHz_ind] = sort(lickBurstHz_all,'ascend');
    nburst = sum(~isnan(lickBurstHz_all));
    nnan = sum(isnan(lickBurstHz_all));
    ind_low_bst = sortlickHz_ind(1:floor(nburst/4));
    ind_high_bst = sortlickHz_ind(nburst-floor(nburst/4)+1:end-nnan);
    HL_lickrate.low_rew = mean(lickBurstHz_all(:,ind_low_bst),2);
    HL_lickrate.high_rew = mean(lickBurstHz_all(:,ind_high_bst),2);
    
    [sortlickHz sortlickHz_ind] = sort(preRew_lickBurstHz,'ascend');
    nnan = sum(isnan(preRew_lickBurstHz));
    nburst = sum(~isnan(preRew_lickBurstHz));
    ind_low_prerew = sortlickHz_ind(1:floor(nburst/4));
    ind_high_prerew = sortlickHz_ind(nburst-floor(nburst/4)+1:end-nnan);
    HL_lickrate.low_prerew = mean(preRew_lickBurstHz(:,ind_low_prerew),2);
    HL_lickrate.high_prerew = mean(preRew_lickBurstHz(:,ind_high_prerew),2);
    
    [sortlickHz sortlickHz_ind] = sort(postRew_lickBurstHz,'ascend');
    nburst = sum(~isnan(postRew_lickBurstHz));
    nnan = sum(isnan(postRew_lickBurstHz));
    ind_low_postrew = sortlickHz_ind(1:floor(nburst/4));
    ind_high_postrew = sortlickHz_ind(nburst-floor(nburst/4)+1:end-nnan);
    HL_lickrate.low_postrew = mean(postRew_lickBurstHz(:,ind_low_postrew),2);
    HL_lickrate.high_postrew = mean(postRew_lickBurstHz(:,ind_high_postrew),2);
    
    figure; % still plotting neural data, but seperate the trials based on licking rate of that trial
    subplot(3,1,1);
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_low_bst),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_low_bst),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_high_bst),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_high_bst),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([-1 inf]);
    title(['1s- Rew: ' num2str(chop(HL_lickrate.low_rew,2)) ' vs ' num2str(chop(HL_lickrate.high_rew,2)) ' Hz']);
    subplot(3,1,2);
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_low_prerew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_low_prerew),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');hold on;
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_high_prerew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_high_prerew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([-1 inf]);
    title(['Pre- Rew: ' num2str(chop(HL_lickrate.low_prerew,2)) ' vs ' num2str(chop(HL_lickrate.high_prerew,2)) ' Hz']);
    subplot(3,1,3);
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_low_postrew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_low_postrew),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');hold on;
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_high_postrew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_high_postrew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([-1 inf]);
    title(['Post- Rew: ' num2str(chop(HL_lickrate.low_postrew,2)) ' vs ' num2str(chop(HL_lickrate.high_postrew,2)) ' Hz']);
    sgtitle([mouse ' ' date '- lick bursts by rate: low (blue) & high (black)']);
    hold off;
    savefig(fullfile(analysis_out,img_fn, '_cueAlignSpiking_byLickRate.fig'));
    
    save(fullfile(analysis_out,img_fn, '_cueAlignLick.mat'), 'firstPostRewLickFrame', ...
        'lastPreRewLickFrame', 'tl_rew', 'firstPostRew_lickAlignEvents', 'lastPreRew_lickAlignEvents', ...
        'lickCueAlign', 'lickBurstStart', 'lickCounterVals', 'lickSearch_frames', 'lickDelay_frames',...
        'lickTC', 'postwin_frames', 'prewin_frames', 'frameRateHz', 'tt', 'ind_early_bst', 'ind_late_bst', ...
        'early_bst_time', 'late_bst_time', 'pct_precue_burst', 'postRew_lickAlignEvents', 'postLick_frames', ...
        'postRew_lickBurstStart','tl', 'ind_prerew_early_bst', 'ind_prerew_late_bst','postRew_lickAlign',...
        'preRew_lickBurstHz','postRew_lickBurstHz','ind_low_prerew','ind_high_prerew','ind_low_postrew',...
        'ind_high_postrew','ind_high_bst','ind_low_bst','HL_lickrate','lickBurstHz_all','preRewTrials',...
        'firstPostRew_lickAlign','lastPreRew_lickAlign','postRewTrials');
    

%{
%find % of licks coorelated with each neuron
%collapse this to licks pre-cue(2000-0)/pre-reward(0-770)/post-reward(770-3000)
%collapse to trials licks start pre-cue and trials licks start post-cue and post-reward
%collapse to trials with licks and trials without licks (for post-learning this is likely meaningless and impossible as 98% of responses will be 'correct')

%after finding the lick neurons from CS+ trials, do these neurons respond during post-cue/post-reward TP during CS- trials
%how does neuronal activity change across the session:
    %do some neurons start firing more as session progresses
    %do some neurons stop firing (response to novelty)
    %do some neurons start firing, period (response to repetition)
    %do some neurons not coorelate with licks at all, if so do they spike:
        %at cue presentation
        %at reward delivery
        %at start of licking only
        
%ALIGN TO LICKS
    %align to first lick 
    %align to last lick
    %align to first burst (will need to define burst based on whole session)
    %

%look at licks during pre-cue
    %naive versus post-learning [jake found that there was a response for naive but not for post-learning]
    %jake's elife looked at the above, jakes NNeuro looked at 
    %look at individual neurons
        %find smooth curves for licking for each trial (within a window so can use -50 +100 from reward)
rewInx = cell2mat(input.tBlock2TrialNumber);
lickSearchRange = prewin_frames+lickDelay_frames:size(lickCueAlign,1)-lickSearch_frames;
%lickTimesForRew = cell(1,length(rewInx)); %preallocate var for initial mworks time licks var
trialLicks = zeros(1,length(rewInx)); %preallocate var for licks in the window of interest for each trial
for iTrial = 1:length(rewInx) %cycle through all trials in a session
    if ~isnan(cTargetOn(1,iTrial)) %ignore trials with no neural data (where cTargetOn and laseron don't align)
        if cTargetOn(iTrial)+postwin_frames-1 < nFrames %ensure that the trial in question exists within the recording
            if rewInx(1,iTrial) == 0 %only do the following for CS+ trials
lickTimesForRew{iTrial} = input.lickometerTimesUs{iTrial}; %pull out mworks time for each lick occurance for a given trial
counterTimes = input.counterTimesUs{iTrial}; %pull out mworks times for each frame, grouped by trial
counterVals = input.counterValues{iTrial}; %pull out frame count for a given trial (size should = counterTimesUs)
lickCounterVals{iTrial} = zeros(size(lickTimesForRew)); %preallocate future var for frames where licks occur
lickTC{iTrial} = zeros(size(counterVals)); %preallocate future var for frames
            for icount = 1:length(counterTimes)-1 %length-1 since we are using a moving index between the start times of two sequential frames
            ind = find(lickTimesForRew>counterTimes(icount) & lickTimesForRew<counterTimes(icount+1)); %find lick occurances from raw data that falls between the start of two frames
                if ~isempty(ind) %if a frame can be found where licks occur
            lickCounterVals{iTrial}(1,ind) = icount; %put that frame value here
                end
            end
            for ival = 1:length(counterTimes)
            ind = find(lickCounterVals{iTrial} == ival); %look at each frame within a trial, did a lick occur?
                if ~isempty(ind) %if a lick occured
            lickTC{iTrial}(1,ival) = length(ind); %put a true (logical = 1) value for that frame in the timecourse for that trial
                end
            end  
            if input.counterValues{iTrial}(end)-cTargetOn(iTrial) > postwin_frames %is the difference between the end of the trial and the cue larger than the value selected for post-cue window (should be true, unless an error occured during imaging that cut the trial short)
            lickCueAlign(:,iTrial) = lickTC{iTrial}(1,cTargetOn(iTrial)-prewin_frames-counterVals(1):cTargetOn(iTrial)+postwin_frames-1-counterVals(1)); %align to the time at the start of the trial and take the lick logical from 500 ms prior to 1000 ms after cue
            end
ind_pre = intersect(preRew_lickSearchRange_700ms,find(lickCueAlign(:,itrial))); %finds every instance of a lick that occurs within the search window
ind_post = intersect(postRew_lickSearchRange,find(lickCueAlign(:,itrial))); %finds every instance of a lick that occurs [after the reward delivery]
ind_all = intersect(lickSearchRange,find(lickCueAlign(:,itrial)));
preRew_lickBurstHz(:,itrial) = length(ind_pre)./(mean(celleqel2mat_padded(input.reactTimesMs))/1000);
postRew_lickBurstHz(:,itrial) = length(ind_post)./(length(postRew_lickSearchRange)./frameRateHz);
lickBurstHz_all(:,itrial) = length(ind_all)./(length(lickSearchRange)./frameRateHz);
ind = intersect(postRew_lickSearchRange,find(lickCueAlign(:,itrial)));

trialLicks = nanmean(lickCueAlign,2);

            end
        end
    end
end

        %find smooth curves for dFoF or spikes for each trial
trialDFOF(:,iTrial) = DFOF from the same window
trialSpike(:,iTrial) = spike Hz from the same window
        %coorelate the two curves to see if there is some relation between neuron activity and licking
   
     for neuIdx = 1:nIC %find the 5 largest events and their frames
[valEvent(:,neuIdx),peakEvent(:,neuIdx)] = maxk(nanmean(targetAlign_events(:,neuIdx,:),3),5);
[valDFOF(:,neuIdx),peakDFOF(:,neuIdx)] = maxk(nanmean(targetAligndFoverF(:,neuIdx,:),3),5);
     end
     
     
%% LOAD IN SOME SOIL
clear;
analysis_out = 'A:\home\carlo\analysis\2P\';
bdata_source = 'A:\home\carlo\rawData\behavioral\';
crpTrainDayList; %type 0 if you are not performing behavioral analysis
ii = length(days);
[mouse, date] = selectMouseInfo(days{ii}); 
RCExptListCarlo_Inter
pairings = loadRCList_Carlo; 
img_fn = [date '_img' mouse '\getTC_000\'];
load([analysis_out, img_fn, '_cueAlignLick.mat']);
load([analysis_out, img_fn, '_input.mat']);
load([analysis_out, img_fn, '_targetAlign.mat']);
%% ACROSS TRIAL BASIS TO DETERMINE IF IT IS WORTH    
isolatedLick = zeros(length(lickCueAlign),length(lickTC)); iT = 1; iF = 1; nIC = size(targetAlign_events,2);
for trialIdx = 1:length(lickTC)
    for frameIdx = 1:length(lickCueAlign)
        if lickCueAlign(frameIdx,trialIdx) == 1 %you found a lick, nice
            for neuIdx = 1:nIC %time for the neurons
            whichN = ['N' num2str(neuIdx)]; %name it so we can save a structure later
            Events = targetAlign_events(:,neuIdx,:); %isolate spike events for a single neuron
            DFOF = targetAligndFoverF(:,neuIdx,:); %isolate changes in fluorescence
            isolated = sum(lickCueAlign(frameIdx-13:frameIdx-1,trialIdx)); %now, are there any licks in the 12 frames (or 400ms) before this one you found?
            if isolated == 0 %No?! fuck yeah 
                isolatedLick(frameIdx,trialIdx) = 1; %throw that lonely lick in this array
                lickEvents(frameIdx,trialIdx) = ~isempty(find(Events(frameIdx-3:frameIdx-1,1,trialIdx)==1,1));
                %now, look at the 3 frames (100ms) before the lick for a neural event
                lickDFOF(frameIdx,trialIdx) = max(DFOF(frameIdx-3:frameIdx-1,1,trialIdx));
                %then, do the same for changes in fluorescence
                    for tIdx = 1:3 %1= frame prior to lick, 3= 100ms prior to lick
                        if DFOF(frameIdx-tIdx,1,trialIdx) == lickDFOF(frameIdx,trialIdx); dFoFIdx(iF,iT) = tIdx; else end
                    end %wanted an index of which frame in the 100ms window was chosen with the highest dFoverF leading up to a lick
                baseDFOF{frameIdx,trialIdx} = DFOF(frameIdx-10:frameIdx-3,1,trialIdx);
                %one more thing, grab the 7 frames prior to the lick window
                %so we can calculate if a change from the prior baseline occured
                lickRA{iF,iT} = {frameIdx,trialIdx}; %finally, put that [frame,trial] index in an array
            else %There were? sucks, we'll get them next time
            end
            %find a lick within CueAlign, create a window around that lick (say 2
            %frames), then look for an event within that window, if an event
            %exists, give that lick a 1 in the mismatch array, else give a 0
            neuron.(whichN).Events = lickEvents;
            neuron.(whichN).DFOF = lickDFOF;
            neuron.(whichN).baselineEvent = nanmean(nanmean(Events,3),1);
            neuron.(whichN).baselineDFOF = nanmean(nanmean(DFOF,3),1);
            end 
        end
    end
end
xTrialIsoLick = nanmean(isolatedLick,2); 
xTrialEvents = nanmean(targetAlign_events,3); 
xTrialDFOF = nanmean(nanmean(targetAligndFoverF,2),3);
    
 


%% NEURON BY NEURON BASIS     
for neuIdx = 1:nIC 
   

end
       

xTrialEvents1 = nanmean(nanmean(targetAlign_events(:,1,:),2),3); 
xTrialEvents2 = nanmean(nanmean(targetAlign_events(:,2,:),2),3); 
xTrialEvents3 = nanmean(nanmean(targetAlign_events(:,3,:),2),3); 
xTrialEvents4 = nanmean(nanmean(targetAlign_events(:,4,:),2),3); 
xTrialEvents5 = nanmean(nanmean(targetAlign_events(:,5,:),2),3); 
xTrialEvents6 = nanmean(nanmean(targetAlign_events(:,6,:),2),3); 
xTrialEvents7 = nanmean(nanmean(targetAlign_events(:,7,:),2),3); 
xTrialEvents8 = nanmean(nanmean(targetAlign_events(:,8,:),2),3); 
xTrialEvents9 = nanmean(nanmean(targetAlign_events(:,9,:),2),3); 
xTrialEvents10 = nanmean(nanmean(targetAlign_events(:,10,:),2),3); 
xTrialEvents11 = nanmean(nanmean(targetAlign_events(:,11,:),2),3); 
xTrialEvents12 = nanmean(nanmean(targetAlign_events(:,12,:),2),3); 
xTrialEvents13 = nanmean(nanmean(targetAlign_events(:,13,:),2),3); 
xTrialEvents14 = nanmean(nanmean(targetAlign_events(:,14,:),2),3); 
xTrialEvents15 = nanmean(nanmean(targetAlign_events(:,15,:),2),3); 
xTrialEvents16 = nanmean(nanmean(targetAlign_events(:,16,:),2),3); 
xTrialEvents17 = nanmean(nanmean(targetAlign_events(:,17,:),2),3); 

figure; plot(xTrialEvents1,'b'); vline(50,'k'); vline(77,'r');
figure; plot(xTrialEvents2,'b'); vline(50,'k'); vline(77,'r');
figure; plot(xTrialEvents3,'b'); vline(50,'k'); vline(77,'r');
figure; plot(xTrialEvents4,'b'); vline(50,'k'); vline(77,'r');
figure; plot(xTrialEvents5,'b'); vline(50,'k'); vline(77,'r');
figure; plot(xTrialEvents6,'b'); vline(50,'k'); vline(77,'r');
figure; plot(xTrialEvents7,'b'); vline(50,'k'); vline(77,'r');
figure; plot(xTrialEvents8,'b'); vline(50,'k'); vline(77,'r');
figure; plot(xTrialEvents9,'b'); vline(50,'k'); vline(77,'r');
figure; plot(xTrialEvents10,'b'); vline(50,'k'); vline(77,'r');
figure; plot(xTrialEvents11,'b'); vline(50,'k'); vline(77,'r');
figure; plot(xTrialEvents12,'b'); vline(50,'k'); vline(77,'r');
figure; plot(xTrialEvents13,'b'); vline(50,'k'); vline(77,'r');
figure; plot(xTrialEvents14,'b'); vline(50,'k'); vline(77,'r');
figure; plot(xTrialEvents15,'b'); vline(50,'k'); vline(77,'r');
figure; plot(xTrialEvents16,'b'); vline(50,'k'); vline(77,'r');
figure; plot(xTrialEvents17,'b'); vline(50,'k'); vline(77,'r');

%}