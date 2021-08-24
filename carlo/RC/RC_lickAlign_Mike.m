%% load in expt data and neural data for each mouse
clear; close all;
analysis_out = 'A:\home\carlo\mikeAnalysis\2P\';
bdata_source = 'A:\home\mike\Data\Behavior\';

%CRP_training_days_lists_mike %clearvars -except an* b* days
%dateIdx = 6;

%for ii = 1:length(days)
%[mouse, date] = selectMouseInfo(days{ii},dateIdx); 
%pulls out mouse ID and date of bx label
%% SJ lick align
RCExptListMike_Inter %ORDER FOR EXPT LIST:: 1507,1511,1510,1512,1513,1516,1520
pairings = loadRCList_Mike;

for  j = 1:size(expt.ttl,2)
    id = j;
    %exp_subset = j;
    iexp = j;  %1:nexp | changed exp_subset to j
    mouse = strtrim(expt(pairings{1,1}).mouse(iexp,:));
    date = expt(pairings{1,1}).date(iexp,:);
    run = expt(pairings{1,1}).run(iexp,:);
    region=run;
    sessions = [date '_img' mouse]; %days{1,ii};
    fprintf([date ' ' mouse ' ' run '\n']);
    img_fn = [date '_img' mouse '\getTC_' run '\'];
    input = getBxData_ForMike(bdata_source, sessions);
    load(fullfile(analysis_out,img_fn,[region '_targetAlign.mat']));
    
    rewDelay_frames =  round((mean(celleqel2mat_padded(input.reactTimesMs))/1000).*frameRateHz);% there's 700ms between the cue and the reward delivery !!!!! if you change tooFastMs in the varibles in MWorks, this needs to be changed
    
    cTargetOn = input.cTargetOn;
    if iscell(cTargetOn) % if it is a cell, it means cTargetOn wasn't being over written in extract TC. If it's not a cell, it's already over written in extract TC. can be used directly
        cTargetOn = celleqel2mat_padded(input.cTargetOn);
        cTargetOn(1) = nan; % first trial doesn't have reward 
    end
    for c = 1:length(cTargetOn)
        if cTargetOn(1,c) == 0
            cTargetOn(1,c) = NaN;
        end
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
    
    %%
    
    figure; %align licking data to cue
    shadedErrorBar(tt, nanmean(lickCueAlign,2).*(1000./frameRateHz), (nanstd(lickCueAlign,[],2)./sqrt(unique(sum(~isnan(lickCueAlign),2))).*(1000./frameRateHz)));
    hold on;
    scatter((lickBurstStart-prewin_frames).*(1000./frameRateHz), 10.*ones(size(lickBurstStart)), 'x');
    xlabel('Time from cue');
    ylabel('Lick rate (Hz)');
    title([date ' ' mouse ' ' run]);
    vline(0,'k');
    vline(700,'r');
    savefig(fullfile(analysis_out,img_fn, [region, '_cueAlign_lickHz.fig']));
    
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
        ylim([-1 13]);
        title([mouse ' ' date ' ' region ' - lick bursts: early (n = ' num2str(length(ind_early_bst)) '); late (n = ' num2str(length(ind_late_bst)) ')']);
        savefig(fullfile(analysis_out,img_fn,[region, '_cueAlignSpiking_byLickTime.fig']));
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
    savefig(fullfile(analysis_out,img_fn,[region, '_postRew_lickBurstAlignSpikingEvents.fig']))

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
    savefig(fullfile(analysis_out,img_fn,[region, '_postRew_lickBurstAlignSpiking.fig']));
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
    title(['Rew- ' num2str(length(preRewTrials))]);
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
    title(['Rew- ' num2str(length(postRewTrials))]);
    supertitle([mouse ' ' date '- Licks relative to reward']);
    hold off;
    savefig(fullfile(analysis_out,img_fn,[region, '_lastVsFirstLick.fig']));
    
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
    supertitle([mouse ' ' date '- lick bursts by rate: low (blue) & high (black)']);
    hold off;
    savefig(fullfile(analysis_out,img_fn,[region, '_cueAlignSpiking_byLickRate.fig']));
    
    save(fullfile(analysis_out,img_fn,[region, '_cueAlignLick.mat']), 'firstPostRewLickFrame', ...
        'lastPreRewLickFrame', 'tl_rew', 'firstPostRew_lickAlignEvents', 'lastPreRew_lickAlignEvents', ...
        'lickCueAlign', 'lickBurstStart', 'lickCounterVals', 'lickSearch_frames', 'lickDelay_frames',...
        'lickTC', 'postwin_frames', 'prewin_frames', 'frameRateHz', 'tt', 'ind_early_bst', 'ind_late_bst', ...
        'early_bst_time', 'late_bst_time', 'pct_precue_burst', 'postRew_lickAlignEvents', 'postLick_frames', ...
        'postRew_lickBurstStart','tl', 'ind_prerew_early_bst', 'ind_prerew_late_bst','postRew_lickAlign',...
        'preRew_lickBurstHz','postRew_lickBurstHz','ind_low_prerew','ind_high_prerew','ind_low_postrew',...
        'ind_high_postrew','ind_high_bst','ind_low_bst','HL_lickrate');
    
end
%end