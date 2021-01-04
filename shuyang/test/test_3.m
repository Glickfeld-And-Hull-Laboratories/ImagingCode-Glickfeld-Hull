for itrial = 1:nTrials
    if cTargetOn(itrial)+postwin_frames-1 < nFrames
        lickTimes = input.lickometerTimesUs{itrial};
        counterTimes = input.counterTimesUs{itrial};
        counterVals = input.counterValues{itrial};
        lickCounterVals{itrial} = zeros(size(lickTimes));
        lickTC{itrial} = zeros(size(counterVals));
        for icount = 1:length(counterTimes)-1
            ind = find(lickTimes>counterTimes(icount) & lickTimes<counterTimes(icount+1));
            if length(ind)>0
                lickCounterVals{itrial}(1,ind) = icount;
            end
        end
        for ival = 1:length(counterTimes)
            ind = find(lickCounterVals{itrial} == ival);
            if length(ind)>0
                lickTC{itrial}(1,ival) = length(ind);
            end
        end
        if input.counterValues{itrial}(end)-cTargetOn(itrial) > postwin_frames
            lickCueAlign(:,itrial) = lickTC{itrial}(1,cTargetOn(itrial)-prewin_frames-counterVals(1):cTargetOn(itrial)+postwin_frames-1-counterVals(1))';
        end
        ind = intersect(lickSearchRange,find(lickCueAlign(:,itrial)));
        for i = 1:length(ind)
            ilick = ind(i);
            if sum(lickCueAlign(ilick:ilick+lickSearch_frames-1,itrial),1) >= 3
                lickBurstStart(:,itrial) = ilick;
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
                postRew_lickAlignEvents(:,:,itrial) = targetAlign_events(ilick-postLick_frames:ilick+postLick_frames-1,:,itrial);
                postRew_lickAlign(:,itrial) = lickCueAlign(ilick-postLick_frames:ilick+postLick_frames-1,itrial);
                break
            end
        end
        ind_pre = find(lickCueAlign(prewin_frames:prewin_frames+rewDelay_frames,itrial),1,'last');
        if ~isempty(ind_pre)
            lastPreRewLickFrame(1,itrial) = ind_pre;
            preRewTrials = [preRewTrials itrial];
            lastPreRew_lickAlignEvents(:,:,itrial) = targetAlign_events(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,:,itrial);
            lastPreRew_lickAlign(:,itrial) = lickCueAlign(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,itrial);
        end
        ind_post = find(lickCueAlign(prewin_frames+rewDelay_frames:prewin_frames+postwin_frames-postLick_frames-postLick_frames,itrial),1,'first');
        if ~isempty(ind_post)
            firstPostRewLickFrame(1,itrial) = ind_post;
            postRewTrials = [postRewTrials itrial];
            firstPostRew_lickAlignEvents(:,:,itrial) = targetAlign_events(ind_post-postLick_frames+prewin_frames+rewDelay_frames:ind_post+postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,:,itrial);
            firstPostRew_lickAlign(:,itrial) = lickCueAlign(ind_post-postLick_frames+prewin_frames+rewDelay_frames:ind_post+postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,itrial);
        end
        rewAlignEvents(:,:,itrial) =targetAlign_events(-postLick_frames+prewin_frames+rewDelay_frames:postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,:,itrial);
    end
end
