clear; 
close all;
bdata_source = 'Z:\Data\behavior\RC\';
analysis_out = 'Z:\2P_analysis\';
RC_imaging_list_SJ;

for  j = 1:length(expt.ttl)
    id = j;
    exp_subset = j;
    iexp = exp_subset;  %1:nexp
    mouse = strtrim(expt.mouse(iexp,:));
    date = expt.date(iexp,:);
    run = expt.run(iexp,:);
    fprintf([date ' ' mouse ' ' run '\n']);
    img_fn = [date '_img' mouse];
    input = get_bx_data_sj(bdata_source, sessions{j});
    load(fullfile(analysis_out,img_fn, [img_fn '_' run '_targetAlign.mat']));
    
    if id == 4
        rewDelay_frames =  round(1.1.*frameRateHz);
    else
        rewDelay_frames =  round(0.77.*frameRateHz);
    end
    
    cTargetOn = input.cTargetOn;
    if iscell(cTargetOn) % if it is a cell, it means cTargetOn wasn't being over written in extract TC. If it's not a cell, it's already over written in extract TC. can be used directly
        cTargetOn = celleqel2mat_padded(input.cTargetOn);
        cTargetOn(1) = nan; % first trial doesn't have reward 
    end
 
    nTrials = size(cTargetOn,2);
    nIC = size(targetAlign_events,2); %you get targetalign_events from neural align, this is the neural data (in spikes) aligned to cue
    lickCueAlign =  nan(prewin_frames+postwin_frames,nTrials);
    lickCounterVals = cell(1,nTrials);
    nFrames = input.counterValues{nTrials}(end);
    lickDelay_frames =  round(0.1.*frameRateHz);
    lickSearch_frames =  round(0.3.*frameRateHz);
    
    postLick_frames =  round(0.5.*frameRateHz);
    lickSearchRange = prewin_frames+lickDelay_frames:size(lickCueAlign,1)-lickSearch_frames;
    preRew_lickSearchRange_600ms = prewin_frames+lickDelay_frames:prewin_frames+lickDelay_frames+rewDelay_frames;
    lickBurstStart = nan(1,nTrials);
    postRew_lickBurstStart = nan(1,nTrials);
    lickBurstHz_all = nan(1,nTrials);
    preRew_lickBurstHz = nan(1,nTrials);
    postRew_lickBurstHz = nan(1,nTrials);
    postRew_lickAlignEvents = nan(2.*postLick_frames,nIC,nTrials);
    postRew_lickAlign = nan(2.*postLick_frames,nTrials);
    lastPreRew_lickAlignEvents = nan(2.*postLick_frames,nIC,nTrials);
    firstPostRew_lickAlignEvents = nan(2.*postLick_frames,nIC,nTrials);
    lastPreRew_lickAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    firstPostRew_lickAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    lastPreRew_lickAlign = nan(3.*postLick_frames,nTrials);
    firstPostRew_lickAlign = nan(3.*postLick_frames,nTrials);
    rewAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    preRewTrials = [];
    postRewTrials = [];
    postRew_lickSearchRange = prewin_frames+rewDelay_frames+lickDelay_frames:size(lickCueAlign,1)-lickSearch_frames-postLick_frames;
    lastPreRewLickFrame = nan(1,nTrials);
    firstPostRewLickFrame = nan(1,nTrials);
    for itrial = 1:nTrials
        if ~isnan(cTargetOn(itrial))
            if cTargetOn(itrial)+postwin_frames-1 < nFrames
                lickTimes = input.lickometerTimesUs{itrial}; 
                counterTimes = input.counterTimesUs{itrial};
                counterVals = input.counterValues{itrial};
                lickCounterVals{itrial} = zeros(size(lickTimes));
                lickTC{itrial} = zeros(size(counterVals));% find how many licks for each frame
                for icount = 1:length(counterTimes)-1
                    ind = find(lickTimes>counterTimes(icount) & lickTimes<counterTimes(icount+1));
                    if ~isempty(ind)
                        lickCounterVals{itrial}(1,ind) = icount; % find which counter in this trial has licks
                    end
                end
                for ival = 1:length(counterTimes)
                    ind = find(lickCounterVals{itrial} == ival);
                    if ~isempty(ind)
                        lickTC{itrial}(1,ival) = length(ind);% if the mice doesn't lick during that counter than it's zero, if it does, it's one. length(ind) should always be 1.
                    end
                end
                if input.counterValues{itrial}(end)-cTargetOn(itrial) > postwin_frames
                    %lickcuealign aligns licking of each frame to cue
                    lickCueAlign(:,itrial) = lickTC{itrial}(1,cTargetOn(itrial)-prewin_frames-counterVals(1):cTargetOn(itrial)+postwin_frames-1-counterVals(1))';
                end
                ind = intersect(lickSearchRange,find(lickCueAlign(:,itrial)));
                for i = 1:length(ind)
                    ilick = ind(i);
                    if sum(lickCueAlign(ilick:ilick+lickSearch_frames-1,itrial),1) >= 3 % more than 3 bursts within about 300ms
                        lickBurstStart(:,itrial) = ilick;% when licking burst happens
                        break
                    end
                end
                ind_pre = intersect(preRew_lickSearchRange_600ms,find(lickCueAlign(:,itrial)));
                ind_post = intersect(postRew_lickSearchRange,find(lickCueAlign(:,itrial)));
                ind_all = intersect(lickSearchRange,find(lickCueAlign(:,itrial)));
                preRew_lickBurstHz(:,itrial) = length(ind_pre)./0.6;
                postRew_lickBurstHz(:,itrial) = length(ind_post)./(length(postRew_lickSearchRange)./frameRateHz);
                lickBurstHz_all(:,itrial) = length(ind_all)./(length(lickSearchRange)./frameRateHz);
                ind = intersect(postRew_lickSearchRange,find(lickCueAlign(:,itrial)));
                for i = 1:length(ind)
                    ilick = ind(i);
                    if sum(lickCueAlign(ilick:ilick+lickSearch_frames-1,itrial),1) >= 3
                        postRew_lickBurstStart(:,itrial) = ilick;
                        postRew_lickAlignEvents(:,:,itrial) = targetAlign_events(ilick-postLick_frames:ilick+postLick_frames-1,:,itrial);% align neural data to first lick after reward
                        postRew_lickAlign(:,itrial) = lickCueAlign(ilick-postLick_frames:ilick+postLick_frames-1,itrial);% align licking data to first lick
                        break
                    end
                end
                ind_pre = find(lickCueAlign(prewin_frames:prewin_frames+rewDelay_frames,itrial),1,'last');
                if ~isempty(ind_pre)
                    lastPreRewLickFrame(1,itrial) = ind_pre;
                    preRewTrials = [preRewTrials itrial];
                    % align neural and licking data to the last lick
                    lastPreRew_lickAlignEvents(:,:,itrial) = targetAlign_events(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,:,itrial);
                    lastPreRew_lickAlign(:,itrial) = lickCueAlign(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,itrial);
                end
                ind_post = find(lickCueAlign(prewin_frames+rewDelay_frames:prewin_frames+postwin_frames-postLick_frames-postLick_frames,itrial),1,'first');
                if ~isempty(ind_post)
                    firstPostRewLickFrame(1,itrial) = ind_post;
                    postRewTrials = [postRewTrials itrial];
                    % align neural and licking data to the first lick
                    firstPostRew_lickAlignEvents(:,:,itrial) = targetAlign_events(ind_post-postLick_frames+prewin_frames+rewDelay_frames:ind_post+postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,:,itrial);
                    firstPostRew_lickAlign(:,itrial) = lickCueAlign(ind_post-postLick_frames+prewin_frames+rewDelay_frames:ind_post+postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,itrial);
                end
                rewAlignEvents(:,:,itrial) =targetAlign_events(-postLick_frames+prewin_frames+rewDelay_frames:postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,:,itrial);
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
    if cue(j) == 1
        vline(0,'b');
    elseif cue(j) == 3
        vline(0,'r');
    elseif cue(j) == 2
        vline(0,'g');
    end
    vline(700,'k');
    savefig(fullfile(analysis_out,img_fn, [img_fn '_' run '_cueAlign_lickHz.fig']));
    
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
        title([mouse ' ' date '- lick bursts: early (n = ' num2str(length(ind_early_bst)) '); late (n = ' num2str(length(ind_late_bst)) ')']);
        savefig(fullfile(analysis_out,img_fn, [img_fn '_' run '_cueAlignSpiking_byLickTime.fig']));
        hold off;
    end
    
    pct_precue_burst = length(find((lickBurstStart-prewin_frames).*(1000./frameRateHz)<600))./size(lickBurstStart,2);
    tl = (1-postLick_frames:postLick_frames).*(1000./frameRateHz);
    
    figure;
    subplot(2,2,1); % align neural activity to lick onset, only burst trials are included 
    shadedErrorBar(tl, nanmean(nanmean(postRew_lickAlignEvents,3),2).*(1000./frameRateHz), (nanstd(nanmean(postRew_lickAlignEvents,3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    hold on;
    xlabel('Time from lick');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    title(['Reward (black- ' num2str(sum(~isnan(squeeze(postRew_lickAlignEvents(1,1,:))))) ')']);
    subplot(2,2,2); % seperate early burst trials vs. late burst trials
    [sortlick sortlick_ind] = sort(postRew_lickBurstStart,'ascend');
    nburst = sum(~isnan(postRew_lickBurstStart));
    nnan = sum(isnan(postRew_lickBurstStart));
    ind_prerew_early_bst = sortlick_ind(1:floor(nburst/4));
    ind_prerew_late_bst = sortlick_ind(nburst-floor(nburst/4)+1:end-nnan);
    early_bst_time = nanmean((postRew_lickBurstStart(:,ind_prerew_early_bst)-prewin_frames-rewDelay_frames).*(1000./frameRateHz),2);
    late_bst_time = nanmean((postRew_lickBurstStart(:,ind_prerew_late_bst)-prewin_frames-rewDelay_frames).*(1000./frameRateHz),2);
    shadedErrorBar(tl, nanmean(nanmean(postRew_lickAlignEvents(:,:,ind_prerew_early_bst),3),2).*(1000./frameRateHz), (nanstd(nanmean(postRew_lickAlignEvents(:,:,ind_prerew_early_bst),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    hold on;
    shadedErrorBar(tl, nanmean(nanmean(postRew_lickAlignEvents(:,:,ind_prerew_late_bst),3),2).*(1000./frameRateHz), (nanstd(nanmean(postRew_lickAlignEvents(:,:,ind_prerew_late_bst),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'b');
    xlabel('Time from lick');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    title(['Black- avg = ' num2str(chop(early_bst_time,3)) 'ms; Blue- avg = ' num2str(chop(late_bst_time,3)) 'ms)'])
    subplot(2,2,3); % align licking data to first lick
    shadedErrorBar(tl, nanmean(postRew_lickAlign,2).*(1000./frameRateHz), (nanstd(postRew_lickAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(postRew_lickAlign),2))),'k');
    hold on;
    xlabel('Time from lick');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    subplot(2,2,4); % plot early burst trials and late burst trials separately
    shadedErrorBar(tl, nanmean(postRew_lickAlign(:,ind_prerew_early_bst),2).*(1000./frameRateHz), (nanstd(postRew_lickAlign(:,ind_prerew_early_bst),[],2).*(1000./frameRateHz))./sqrt(length(ind_prerew_early_bst)),'k');
    hold on;
    shadedErrorBar(tl, nanmean(postRew_lickAlign(:,ind_prerew_late_bst),2).*(1000./frameRateHz), (nanstd(postRew_lickAlign(:,ind_prerew_late_bst),[],2).*(1000./frameRateHz))./sqrt(length(ind_prerew_late_bst)),'b');
    xlabel('Time from lick');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    supertitle([mouse ' ' date '- post reward lick burst aligned spiking']);
    savefig(fullfile(analysis_out,img_fn, [img_fn '_' run '_postRew_lickAlignSpiking.fig']));
    hold off;
    
    figure;
    for i = 1:4
        plot(tt,cumsum(nansum(lickCueAlign(:,1+((i-1)*floor(nTrials/4)):floor(nTrials/4)+((i-1)*floor(nTrials/4))),2)));
        hold on;
    end
    title([mouse ' ' date '- cumulative licking by quarter session']);
    xlabel('Time from cue');
    ylabel('Cumulative Licks');
    hold off;
    savefig(fullfile(analysis_out,img_fn, [img_fn '_' run '_cumulativeLicking.fig']));
    
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
    savefig(fullfile(analysis_out,img_fn, [img_fn '_' run '_lastVsFirstLick.fig']));
    
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
    subplot(1,3,1);
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_low_bst),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_low_bst),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_high_bst),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_high_bst),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([-1 inf]);
    title(['1s- Rew: ' num2str(chop(HL_lickrate.low_rew,2)) ' vs ' num2str(chop(HL_lickrate.high_rew,2)) ' Hz']);
    subplot(1,3,2);
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_low_prerew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_low_prerew),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');hold on;
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_high_prerew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_high_prerew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([-1 inf]);
    title(['Pre- Rew: ' num2str(chop(HL_lickrate.low_prerew,2)) ' vs ' num2str(chop(HL_lickrate.high_prerew,2)) ' Hz']);
    subplot(1,3,3);
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_low_postrew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_low_postrew),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');hold on;
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_high_postrew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_high_postrew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([-1 inf]);
    title(['Post- Rew: ' num2str(chop(HL_lickrate.low_postrew,2)) ' vs ' num2str(chop(HL_lickrate.high_postrew,2)) ' Hz']);
    supertitle([mouse ' ' date '- lick bursts by rate: low (blue) & high (black)']);
    hold off;
    savefig(fullfile(analysis_out,img_fn, [img_fn '_' run '_cueAlignSpiking_byLickRate.fig']));
    
    save(fullfile(analysis_out,img_fn, [img_fn '_' run '_cueAlignLick.mat']), 'firstPostRewLickFrame', ...
        'lastPreRewLickFrame', 'tl_rew', 'firstPostRew_lickAlignEvents', 'lastPreRew_lickAlignEvents', ...
        'lickCueAlign', 'lickBurstStart', 'lickCounterVals', 'lickSearch_frames', 'lickDelay_frames',...
        'lickTC', 'postwin_frames', 'prewin_frames', 'frameRateHz', 'tt', 'ind_early_bst', 'ind_late_bst', ...
        'early_bst_time', 'late_bst_time', 'pct_precue_burst', 'postRew_lickAlignEvents', 'postLick_frames', ...
        'postRew_lickBurstStart','tl', 'ind_prerew_early_bst', 'ind_prerew_late_bst','postRew_lickAlign',...
        'preRew_lickBurstHz','postRew_lickBurstHz','ind_low_prerew','ind_high_prerew','ind_low_postrew',...
        'ind_high_postrew','ind_high_bst','ind_low_bst','HL_lickrate');
    
end