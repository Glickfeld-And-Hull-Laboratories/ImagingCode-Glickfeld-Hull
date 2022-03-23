% this script is the same as RC_lickAlign_SJ, but the window for aligning
% the lick onset is longer. (1.5s before and 3s after instead of 0.5s before and after)
% did this to make it the same scale as neural align, for Court's RO1 grant
% This script should run AFTER running RC_lickAlign_SJ!!!!!!!!!!!!

clear; 
close all;
bdata_source = 'Z:\home\shuyang\Data\behavior\RC\';
analysis_out = 'Z:\home\shuyang\2P_analysis\';
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
    %load(fullfile(analysis_out,img_fn, [img_fn '_' run '_HPfiltered_targetAlign.mat']));
    load(fullfile(analysis_out,img_fn, [img_fn '_' run '_targetAlign.mat']));
    threshold = -3;
    %load(fullfile([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_spk_cutoff_0.035' '_deconvolve_threshold' num2str(threshold) '.mat' ]));
    load(fullfile([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_spk_deconvolve_threshold' num2str(threshold) '.mat' ]));
    

    cTargetOn = input.cTargetOn;
    if iscell(cTargetOn) % if it is a cell, it means cTargetOn wasn't being over written in extract TC. If it's not a cell, it's already over written in extract TC. can be used directly
        cTargetOn = cell2mat(input.cTargetOn);
        cTargetOn(1) = nan; % first trial doesn't have reward 
    end
 
    nTrials = size(cTargetOn,2);
    
    precue = round(1500./frameRateHz);% this value is same as prewin_frames
    postcue = round(5800./frameRateHz);% the length you wanna include after cue onset, 5.8s is basically the whole trial
    nIC = size(spk_logic_cl,2);
    targetAlign_events_wholeTrial = nan(precue+postcue,nIC,nTrials);% frame*cell*trial. if your cTargetOn has nan, targetAlign_events will have nan. so use nanmean in the following analysis when average
    nFrames = size(all_events,1);
    
    for itrial = 1:nTrials
        if ~isnan(cTargetOn(itrial))&&cTargetOn(itrial)>0 % don't know why but in 201114_img1079 cTargetOn(1) = 0
            if cTargetOn(itrial)+postcue-1 <= nFrames %& input.counterValues{itrial}(end)-cTargetOn(itrial) > postwin_frames
                targetAlign_events_wholeTrial(:,:,itrial) = all_events(cTargetOn(itrial)-precue:cTargetOn(itrial)+postcue-1,:);
            end
        end
    end
    
    nIC = size(targetAlign_events_wholeTrial,2); %you get targetalign_events from neural align, this is the neural data (in spikes) aligned to cue
    lickCueAlign_wholeTrial =  nan(precue+postcue,nTrials);
    lickCounterVals = cell(1,nTrials);
    nFrames = input.counterValues{nTrials}(end);
    lickDelay_frames =  round(0.1.*frameRateHz);
    lickSearch_frames =  round(0.3.*frameRateHz);
    
    postLick_frames =  round(0.5.*frameRateHz);
    prewin_frames = round(1500./frameRateHz); % this value is same as precue
    postwin_frames = round(3000./frameRateHz);
    rewDelay_frames =  round(0.7.*frameRateHz); % there's 700ms between the cue and the reward delivery !!!!! if you change tooFastMs in the varibles in MWorks, this needs to be changed
    lickSearchRange = prewin_frames+lickDelay_frames:size(lickCueAlign_wholeTrial,1)-lickSearch_frames;
    preRew_lickSearchRange_700ms = prewin_frames+lickDelay_frames:prewin_frames+lickDelay_frames+rewDelay_frames;
    lickBurstStart = nan(1,nTrials);
    postRew_lickBurstStart = nan(1,nTrials);
    lickBurstHz_all = nan(1,nTrials);
    preRew_lickBurstHz = nan(1,nTrials);
    postRew_lickBurstHz = nan(1,nTrials);
    postRew_lickAlignEvents_1500_3000ms_scale = nan(prewin_frames+postwin_frames,nIC,nTrials);
    postRew_lickAlign_1500_3000ms_scale = nan(prewin_frames+postwin_frames,nTrials);
    lastPreRew_lickAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    firstPostRew_lickAlignEvents_1500_3000ms_scale = nan(prewin_frames+postwin_frames,nIC,nTrials);
    lastPreRew_lickAlign = nan(3.*postLick_frames,nTrials);
    firstPostRew_lickAlign_1500_3000ms_scale = nan(prewin_frames+postwin_frames,nTrials);
    rewAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    preRewTrials = [];
    postRewTrials = [];
    postRew_lickSearchRange = prewin_frames+rewDelay_frames+lickDelay_frames:size(lickCueAlign_wholeTrial,1)-lickSearch_frames-postLick_frames;
    lastPreRewLickFrame = nan(1,nTrials);
    firstPostRewLickFrame = nan(1,nTrials);
    for itrial = 1:nTrials
        if ~isnan(cTargetOn(itrial)) && cTargetOn(itrial)>0
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
                    lickCueAlign_wholeTrial(:,itrial) = lickTC{itrial}(1,cTargetOn(itrial)-precue-counterVals(1):cTargetOn(itrial)+postcue-1-counterVals(1))';
                end
                ind = intersect(lickSearchRange,find(lickCueAlign_wholeTrial(:,itrial)));
                for i = 1:length(ind)
                    ilick = ind(i);
                    if sum(lickCueAlign_wholeTrial(ilick:ilick+lickSearch_frames-1,itrial),1) >= 3 % more than 3 bursts within about 300ms
                        lickBurstStart(:,itrial) = ilick;% when licking burst happens
                        break
                    end
                end
                ind_pre = intersect(preRew_lickSearchRange_700ms,find(lickCueAlign_wholeTrial(:,itrial)));
                ind_post = intersect(postRew_lickSearchRange,find(lickCueAlign_wholeTrial(:,itrial)));
                ind_all = intersect(lickSearchRange,find(lickCueAlign_wholeTrial(:,itrial)));
                preRew_lickBurstHz(:,itrial) = length(ind_pre)./0.6;
                postRew_lickBurstHz(:,itrial) = length(ind_post)./(length(postRew_lickSearchRange)./frameRateHz);
                lickBurstHz_all(:,itrial) = length(ind_all)./(length(lickSearchRange)./frameRateHz);
                ind = intersect(postRew_lickSearchRange,find(lickCueAlign_wholeTrial(:,itrial)));
                for i = 1:length(ind)
                    ilick = ind(i);
                    if sum(lickCueAlign_wholeTrial(ilick:ilick+lickSearch_frames-1,itrial),1) >= 3
                        postRew_lickBurstStart(:,itrial) = ilick;
                        if (ilick+postwin_frames-1) <= size(targetAlign_events_wholeTrial,1) % if the time period before laser is off after the first lick is long enough for this trial
                            % align data to first lick burst after reward
                            postRew_lickAlignEvents_1500_3000ms_scale(:,:,itrial) = targetAlign_events_wholeTrial(ilick-prewin_frames:ilick+postwin_frames-1,:,itrial);
                            postRew_lickAlign_1500_3000ms_scale(:,itrial) = lickCueAlign_wholeTrial(ilick-prewin_frames:ilick+postwin_frames-1,itrial);% align licking data to first lick burst
                        break
                        end  
                    end
                end
                ind_pre = find(lickCueAlign_wholeTrial(prewin_frames:prewin_frames+rewDelay_frames,itrial),1,'last');
                if ~isempty(ind_pre)
                    lastPreRewLickFrame(1,itrial) = ind_pre;
                    preRewTrials = [preRewTrials itrial];
                    % align neural and licking data to the last lick
                    lastPreRew_lickAlignEvents(:,:,itrial) = targetAlign_events_wholeTrial(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,:,itrial);
                    lastPreRew_lickAlign(:,itrial) = lickCueAlign_wholeTrial(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,itrial);
                end
                ind_post = find(lickCueAlign_wholeTrial(prewin_frames+rewDelay_frames:precue+postcue-1,itrial),1,'first'); % prewin_frames + rewDelay_frames = reward delivery from trial onset
                if ~isempty(ind_post)
                    % align neural and licking data to the first lick after reward, take # of prewin frames before lick and # of postwin frames after lick
                    if ind_post+postwin_frames-1+prewin_frames+rewDelay_frames <= size(targetAlign_events_wholeTrial,1) % if the time period before laser is off after the first lick is long enough for this trial
                        firstPostRewLickFrame(1,itrial) = ind_post;
                        postRewTrials = [postRewTrials itrial];
                        firstPostRew_lickAlignEvents_1500_3000ms_scale(:,:,itrial) = targetAlign_events_wholeTrial(ind_post-prewin_frames-1+prewin_frames+rewDelay_frames:ind_post+postwin_frames-2+prewin_frames+rewDelay_frames,:,itrial);
                        firstPostRew_lickAlign_1500_3000ms_scale(:,itrial) = lickCueAlign_wholeTrial(ind_post-prewin_frames-1+prewin_frames+rewDelay_frames:ind_post+postwin_frames+prewin_frames-2+rewDelay_frames,itrial);
                    end
                end
                rewAlignEvents(:,:,itrial) =targetAlign_events_wholeTrial(-postLick_frames+prewin_frames+rewDelay_frames:postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,:,itrial);
            end
        end
    end
    
    tt = (-prewin_frames:postwin_frames-1).*(1000./frameRateHz);
    %save(fullfile(analysis_out,img_fn, [img_fn '_' run '_HPfiltered_cueAlignLick.mat']), ...
    save(fullfile(analysis_out,img_fn, [img_fn '_' run '_cueAlignLick.mat']), ...
        'postRew_lickAlignEvents_1500_3000ms_scale', 'postRew_lickAlign_1500_3000ms_scale',...
        'firstPostRew_lickAlignEvents_1500_3000ms_scale','firstPostRew_lickAlign_1500_3000ms_scale',...
        '-append');
    
    
    figure;
    subplot(2,1,1); % align neural activity to lick burst onset, only burst trials are included 
    shadedErrorBar(tt, nanmean(nanmean(postRew_lickAlignEvents_1500_3000ms_scale,3),2).*(1000./frameRateHz), (nanstd(nanmean(postRew_lickAlignEvents_1500_3000ms_scale,3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    hold on;
    xlabel('Time from lick');
    ylabel('Spike rate (Hz)');
    vline(0,'k'); ylim([0 inf]);
    title(['Reward (black- ' num2str(sum(~isnan(squeeze(postRew_lickAlignEvents_1500_3000ms_scale(1,1,:))))) ')']);
    subplot(2,1,2); % align licking data to first lick burst
    shadedErrorBar(tt, nanmean(postRew_lickAlign_1500_3000ms_scale,2).*(1000./frameRateHz), (nanstd(postRew_lickAlign_1500_3000ms_scale,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(postRew_lickAlign_1500_3000ms_scale),2))),'k');
    hold on;
    xlabel('Time from lick');
    ylabel('Lick rate (Hz)');
    vline(0,'k'); ylim([0 inf]);
    supertitle([mouse ' ' date '- post reward lick burst aligned spiking']);
    %savefig(fullfile(analysis_out,img_fn, [img_fn '_' run '_HPfiltered_postRew_lickAlignSpiking_1500_3000ms_scale.fig']));
    savefig(fullfile(analysis_out,img_fn, [img_fn '_' run '_postRew_lickAlignSpiking_1500_3000ms_scale.fig']));
    hold off;
    
    figure; % align neural and licking data to the first lick
    subplot(2,1,1);
    shadedErrorBar(tt, nanmean(nanmean(firstPostRew_lickAlignEvents_1500_3000ms_scale,3),2).*(1000./frameRateHz), (nanstd(nanmean(firstPostRew_lickAlignEvents_1500_3000ms_scale,3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    title('First lick after reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    vline(0,'k'); ylim([0 inf]);
    subplot(2,1,2);
    shadedErrorBar(tt, nanmean(firstPostRew_lickAlign_1500_3000ms_scale,2).*(1000./frameRateHz), (nanstd(firstPostRew_lickAlign_1500_3000ms_scale,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(firstPostRew_lickAlign_1500_3000ms_scale)))),'k');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    vline(0,'k'); ylim([0 inf]);
    title(['Rew- ' num2str(length(postRewTrials))]);
    supertitle([mouse ' ' date '- Licks relative to reward']);
    hold off;
    %savefig(fullfile(analysis_out,img_fn, [img_fn '_' run '_HPfiltered_lastVsFirstLick_1500_3000ms_scale.fig']));
    savefig(fullfile(analysis_out,img_fn, [img_fn '_' run '_lastVsFirstLick_1500_3000ms_scale.fig']));
    
end