clear; close all;
analysis_out = 'A:\home\carlo\analysis\2P\';
bdata_source = 'A:\home\carlo\rawData\behavioral\';

removing_bad_trials = false; %bad_trials = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]; %added to remove shitty trials, turn off if not needed

expType=input('monomodal(1) or dimodal(2) session: ');
seshType=input('interleaved (1) or blocked (2) session: ');
if seshType==1 && expType == 1
RCExptListCarlo_Inter;
elseif seshType==2 && expType == 1
RCExptListCarlo_Block;
elseif seshType==1 && expType == 2
RCExptListCarlo_VA_Inter;
end

% RCExptListCarlo_Inter
% RCExptListCarlo_Block
pairings = loadRCList_Carlo(expt,seshType);

%for  j = 1:size(expt.ttl,2)

    id = pairings{1,1};
    exp_subset = id;
    iexp = exp_subset; isesh = pairings{3,1};
mouse = strtrim(expt{1,iexp}.mouse{1,isesh});
date = expt{1,iexp}.date{1,isesh};
run = expt{1,iexp}.run{1,isesh};
    fprintf([date ' ' mouse ' ' run '\n']);
    img_fn = [date '_img' mouse '\getTC_' run '\'];
    filename = dir([analysis_out,img_fn  '*' '_TCave.mat']);
    load([analysis_out,img_fn, filename.name]);
    filename = dir([analysis_out,date '_img' mouse '\' date '_img' mouse '_' run  '*' '_TCave_cl.mat']);
    load([analysis_out,date '_img' mouse '\', filename.name]);
    %threshold = -3;
    spk = load([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_spk_deconvolve_threshold' num2str(thresholdDeco) '.mat' ]);
    spk_logic = spk.spk_logic_cl;
    nIC = size(spk_logic,2);
    %% fill laser off period with nans
    all_events = spk_logic;
    if expt{1,iexp}.ttl(1,isesh)
        all_events_temp = nan(size(tc_avg_all,1),size(spk_logic,2));
        all_events_temp(laseron(1,1:end),:) = all_events;
        all_events = all_events_temp;
    end
    
    save(fullfile([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_spk_deconvolve_threshold' num2str(thresholdDeco) '.mat' ]), 'all_events', '-append')
    
    sessions = [date '_img' mouse];% for behavior data - format YYMMDD[2/0/1]_imgMOUSE 2 = interleaved; 0 = CS- block; 1 = CS+ block
    if strlength(sessions) < 16
        bfile = dir([bdata_source 'data-i' mouse '-' date '*' ]);
    else
        bfile = dir([bdata_source 'data-i' mouse '-' date sessions(end-4:end) '.mat']);
    end
%load behavior file
    behav_dest = [bdata_source bfile.name];
    assert(length(bfile)) = 1;
    b_data = load(behav_dest);
    mworks = b_data.input;
    
    react_time = double(cell2mat(mworks.reactTimesMs(1,2:end)));
    cue_rew_int = mworks.RewardDelayDurationMs + round(mean(react_time),-1); %"total interval between cueonset and reward delivery"       
    img_fn2 = sessions;
    frameRateHz = double(mworks.frameRateHz);
    cTargetOn = (mworks.cTargetOn);
    nTrials = length(cTargetOn);
    prewin_frames = round(1500./frameRateHz);
    postwin_frames = round(3000./frameRateHz);
    targetAlign_tc = nan(prewin_frames+postwin_frames,nIC,nTrials);
    targetAlign_events = nan(prewin_frames+postwin_frames,nIC,nTrials);
    nFrames = size(tc_avg_all,1);
    
    % adding a processing step that was done in deconvolution to return
    % filtered trace to the raw fluorescence values (if the fluorescence
    % mean didn't change due to mechanical drift)
%     preFilterFst200Frms = nanmean(raw_tc_avg(1:200,goodcells),1); %take the first 200 frames from raw trace for only cells that passed deconvolution
%     all_TCave_clProc = all_TCave_cl+preFilterFst200Frms; %add that value to the filtered TC to return to its original mean (in an attempt to make the filtered data AS comparable to the raw trace as possible, since in an ideal world they would be identical)
%     figure; plot(nanmean(all_TCave_cl,2),'r'); hold on; plot(nanmean(raw_tc_avg,2),'k'); plot(nanmean(all_TCave_clProc,2),'g');
    %process should be obvious from figure (shift of the red trace to the level of the black trace to 
    %give the green trace, which will be used for aligning to cue presentation - important note:
    %since raw_tc_avg is formatted to include only laseron periods, it is SHORTER than the all_TCave 
    %traces, this is expected and should make sense, focus shouldn't be on the length match, 
    %but the match at the y-axis
    
    for itrial = 1:nTrials
        if cTargetOn(itrial)+postwin_frames-1 <= nFrames %& input.counterValues{itrial}(end)-cTargetOn(itrial) > postwin_frames
            targetAlign_tc(:,:,itrial) = all_TCave_cl(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:);
            targetAlign_events(:,:,itrial) = all_events(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:);
        end
    end
    
    targetAlignF = nanmean(targetAlign_tc(1:prewin_frames,:,:),3); % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.
    targetAligndFoverF = zeros(size(targetAlign_tc,1),size(targetAlign_tc,2),size(targetAlign_tc,3));%frame*cell*trials
    % calculate df/F
    for c = 1:size(TCave_cl,2)
        targetAligndFoverF(:,c,:) = (targetAlign_tc(:,c,:)-targetAlignF(c))./targetAlignF(c);
    end
    
    
        omitRewTrial = celleqel2mat_padded(mworks.tRewardOmissionTrial(1,2:end));
        ind_omit = find(omitRewTrial);
        
        unexpRewTrial = celleqel2mat_padded(mworks.tDoNoStimulusChange(1,2:end));
        ind_unexp = find(unexpRewTrial);
        
        if mworks.doBlock2 %Specific CS+/CS- experiment added by Mike
            block2Trial = celleqel2mat_padded(mworks.tBlock2TrialNumber(1,2:end));
            ind_block2 = find(block2Trial); 
            ind_rew = find(~block2Trial); 
        else
            ind_rew = intersect(find(omitRewTrial == 0),find(unexpRewTrial == 0));
            if intersect(ind_omit, ind_unexp)
                x = ismember(ind_omit,intersect(ind_omit, ind_unexp));
                ind_omit(x) = [];
                x = ismember(ind_unexp,intersect(ind_omit, ind_unexp));
                ind_unexp(x) = [];
            end
        end
        
        
        if removing_bad_trials
            for b = 1:length(bad_trials)
                ind_omit = ind_omit(ind_omit ~= bad_trials(b));
                ind_unexp = ind_unexp(ind_unexp ~= bad_trials(b));
                ind_rew = ind_rew(ind_rew ~= bad_trials(b));
                if mworks.doBlock2
                    ind_block2 = ind_block2(ind_block2 ~= bad_trials(b));
                end
            end
        end
        
        tt = (-prewin_frames:postwin_frames-1).*(1000./frameRateHz);
        rewardDelayDurationMs = double(max(celleqel2mat_padded(mworks.tRewardDelayDurationMs(1,2:end)),[],2));
        reactTimesMs = double(mean(celleqel2mat_padded(mworks.reactTimesMs(1,2:end)),2));
        delayTimeMs = reactTimesMs+rewardDelayDurationMs;
        
        save([analysis_out, img_fn, 'figDataUntrimmed_' sessions '.mat']);
%% PLOTS    
        s=1;
        if length(ind_omit>5); s = s+1; end
        if length(ind_unexp>5); s = s+1; end
        if mworks.doBlock2; s = s+1; end
        
        if unique(expt{exp_subset}.name == 'Block')
            indivFig = true;
        else
            indivFig = false;
        end
            
        if ~indivFig 
        figure; n = 1;
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_rew),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_rew),3),[],2)./sqrt(nIC), 'k');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        %ylim([-3 5]);
        title('CS+')
        n = n+1;
        if length(ind_omit>5)
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_omit),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_omit),3),[],2)./sqrt(nIC),'r');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        title('Omit')
        n = n+1;
        end
        if length(ind_unexp>5)
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_unexp),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_unexp),3),[],2)./sqrt(nIC),'g');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        title('Unexpected reward')
        n = n+1;
        end
        if mworks.doBlock2
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_block2),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_block2),3),[],2)./sqrt(nIC),'r');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        %ylim([-3 5]);
        title('CS-')
        end
        sgtitle([date ' ' mouse])
        savefig(fullfile(analysis_out,img_fn, '_cueAlign_dFoverF.fig'))

        figure; n=1;
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew),3),[],2)./sqrt(nIC)).*(1000./frameRateHz));
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([0 3]);
        title('CS+')
        n = n+1;
        if length(ind_omit>5)
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_omit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'b');
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Omit')
        n = n+1;
        end
        if length(ind_unexp>5)
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_unexp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_unexp),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'g');
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Unexpected reward')
        n = n+1;
        end
        if mworks.doBlock2
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_block2),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_block2),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'r');
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([0 3]);
        title('CS-')
        end
        sgtitle([date ' ' mouse])
        savefig(fullfile(analysis_out,img_fn, '_cueAlign_events_Hz.fig'))
        
        elseif indivFig
            
        figure; %n = 1;
        %subplot(s,1,n) 
        if expt{exp_subset}.whichBlock(1,isesh) == 0
            shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_rew),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_rew),3),[],2)./sqrt(nIC),'red');
            vline(cue_rew_int)        
            xlabel('Time from cue')
            ylabel('Avg dF/F')
            title('CS-')
        elseif expt{exp_subset}.whichBlock(1,isesh) == 1
            shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_rew),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_rew),3),[],2)./sqrt(nIC),'black');
            vline(cue_rew_int)        
            xlabel('Time from cue')
            ylabel('Avg dF/F')
            title('CS+')
        end
        sgtitle([date ' ' mouse])
        savefig(fullfile(analysis_out,img_fn, '_cueAlign_dFoverF.fig'))

        figure; %n=1;
        %subplot(s,1,n)
        if expt{exp_subset}.whichBlock(1,isesh) == 0
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'red');
            vline(cue_rew_int)
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([0 3]);
            title('CS-')
        elseif expt{exp_subset}.whichBlock(1,isesh) == 1
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'black');
            vline(cue_rew_int)
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([0 4]);
            title('CS+')
        end
        sgtitle([date ' ' mouse])
        savefig(fullfile(analysis_out,img_fn, '_cueAlign_events_Hz.fig'))
        end
        
        
        if length(ind_omit>10)
            figure;
            subplot(2,1,1)
            thresh = mean(diff(ind_omit));
            ind_omit_short = find(diff(ind_omit)<thresh)+1;
            ind_omit_long = find(diff(ind_omit)>thresh)+1;
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_omit_short)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_omit_short)),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_omit_long)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_omit_long)),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'b');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([0 2]);
            title(['Short (black n = ' num2str(length(ind_omit_short)) ') vs long (blue n = ' num2str(length(ind_omit_long)) ') interval between omits'])
            ind_rew_preomit = ind_omit-1;
            if find(ind_rew_preomit==0)
              ind_rew_preomit(find(ind_rew_preomit==0)) = [];
            end
            ind_rew_postomit = ind_omit+1;
            if find(ind_rew_postomit)
              ind_rew_postomit(find(ind_rew_postomit>nTrials)) = [];
            end
            subplot(2,1,2)
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew_preomit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_preomit),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew_postomit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_postomit),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'b');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Reward trial pre (black n = ' num2str(length(ind_rew_preomit)) ') vs post (blue n = ' num2str(length(ind_rew_postomit)) ') omit trials'])
            sgtitle([date ' ' mouse])
            savefig(fullfile(analysis_out,img_fn, [img_fn '_omitByInterval.fig']))
        else
            ind_omit_short = [];
            ind_omit_long = [];
            ind_rew_preomit = [];
            ind_rew_postomit = [];
        end
        if length(ind_unexp>10)
            figure;
            thresh = mean(diff(ind_unexp));
            ind_unexp_short = find(diff(ind_unexp)<thresh)+1;
            ind_unexp_long = find(diff(ind_unexp)>thresh)+1;
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_short)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_short)),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_long)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_long)),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'b');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title([date ' ' mouse '- Short (black n = ' num2str(length(ind_unexp_short)) ') vs long (blue n = ' num2str(length(ind_unexp_long)) ') interval between unexpected reward'])
            savefig(fullfile(analysis_out,img_fn, [img_fn '_unexpByInterval.fig']))
        else
            ind_unexp_short = [];
            ind_unexp_long = [];
        end
        
        if unique(expt{exp_subset}.name == 'Inter')
         figure; n=1;
        %subplot(s,1,n)
        plot(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz),'black');
        hold on;
        plot(tt, nanmean(nanmean(targetAlign_events(:,:,ind_block2),3),2).*(1000./frameRateHz),'r');
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        %ylim([0 6]);
        title('Stacked spike rates separating CS+ and CS- trials');
        legend('CS+','CS-');
        sgtitle([date ' ' mouse])
        savefig(fullfile(analysis_out,img_fn, '_line_stackedSpike_events_Hz.fig'))
        
         figure; n=1;
        %subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'black');
        hold on;
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_block2),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_block2),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'r');
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        %ylim([0 6]);
        title('Stacked spike rates separating CS+ and CS- trials');
        %legend('CS+', 'sd CS+' ,'CS-', 'sd CS-');
        sgtitle([date ' ' mouse])
        savefig(fullfile(analysis_out,img_fn, '_error_stackedSpike_events_Hz.fig'))
        
        end 
        
        if unique(expt{exp_subset}.name == 'Inter')
        figure;
        n = floor(nTrials./3);
        seshType = 1;
        for i = 1:3
            subplot(3,2,seshType)
            ind_rew_temp = intersect(ind_rew,1+(i-1).*n:i*n);
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew_temp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_temp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
            hold on
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([0 6]);
            title(['Trials ' num2str(1+(i-1).*n) ':' num2str(i*n)])
            vline(cue_rew_int)
            if length(ind_block2)>=10
                subplot(3,2,seshType+1)
                ind_block2_temp = intersect(ind_block2,1+(i-1).*n:i*n);
                shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_block2_temp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_block2_temp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'r');
                ylim([0 6])
                xlabel('Time from cue')
                ylabel('Spike rate (Hz)')
                title(['Trials ' num2str(1+(i-1).*n) ':' num2str(i*n)])
                vline(cue_rew_int)
            end
            seshType = seshType+2;
        end
        sgtitle([date ' ' mouse '- CS+ (black), CS- (red)'])
        savefig(fullfile(analysis_out,img_fn, '_repsByTrial.fig'))
        
        elseif unique(expt{exp_subset}.name == 'Block')
            if expt{exp_subset}.whichBlock(1,isesh) == 0
        figure;
        n = floor(nTrials./3);
        seshType = 1;
        for i = 1:3
            subplot(3,2,seshType)
            ind_rew_temp = intersect(ind_rew,1+(i-1).*n:i*n);
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew_temp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_temp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'r');
            hold on
            ylim([0 6])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Trials ' num2str(1+(i-1).*n) ':' num2str(i*n)])
            vline(cue_rew_int)
           
            seshType=seshType+2;
        end 
        sgtitle([date ' ' mouse '- CS- (red)'])
        savefig(fullfile(analysis_out,img_fn, '_repsByTrial.fig'))
            elseif expt{exp_subset}.whichBlock(1,isesh) == 1
        figure;
        n = floor(nTrials./3);
        seshType = 1;
        for i = 1:3
            subplot(3,2,seshType)
            ind_rew_temp = intersect(ind_rew,1+(i-1).*n:i*n);
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew_temp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_temp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
            hold on
            ylim([0 6])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Trials ' num2str(1+(i-1).*n) ':' num2str(i*n)])
            vline(cue_rew_int)
        
            seshType=seshType+2;
        end
        sgtitle([date ' ' mouse '- CS+ (black)'])
        savefig(fullfile(analysis_out,img_fn, '_repsByTrial.fig'))
            end
        end
        %{    
        if unique(expt.name == 'Crus')
            load(fullfile(analysis_out,img_fn, [img_fn '_splitImage.mat']))
            figure;
            indL = find(maskCat==1);
            indR = find(maskCat==2);
            subplot(2,2,1)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,indL,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,indL,ind_rew),3),[],2).*(1000./frameRateHz))./sqrt(length(indL)),'k');
            ylim([0 5])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Reward- Left side- n=' num2str(length(indL))])
            vline(cue_rew_int)
            subplot(2,2,2)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,indR,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,indR,ind_rew),3),[],2).*(1000./frameRateHz))./sqrt(length(indR)),'k');
            ylim([0 5])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            vline(cue_rew_int)
            title(['Reward- Right side- n=' num2str(length(indR))])
            subplot(2,2,3)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,indL,ind_omit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,indL,ind_omit),3),[],2).*(1000./frameRateHz))./sqrt(length(indL)),'r');
            ylim([0 5])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Omit- Left side- n=' num2str(length(indL))])
            vline(cue_rew_int)
            subplot(2,2,4)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,indR,ind_omit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,indR,ind_omit),3),[],2).*(1000./frameRateHz))./sqrt(length(indR)),'r');
            ylim([0 5])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['omit- Right side- n=' num2str(length(indR))])
            vline(cue_rew_int)
            sgtitle([date ' ' mouse '- Reward (black), Omit (red)'])
            savefig(fullfile(analysis_out,img_fn, [img_fn '_repsByCrus.fig']))
        end
        %}
        
        save(fullfile(analysis_out,img_fn, '_targetAlign.mat'), 'ind_rew', 'ind_block2', 'ind_omit', 'ind_unexp', 'targetAlign_events', 'targetAligndFoverF', 'prewin_frames', 'postwin_frames', 'tt', 'frameRateHz', 'ind_omit_short','ind_omit_long','ind_unexp_short','ind_unexp_long','ind_rew_preomit','ind_rew_postomit')
        save(fullfile(analysis_out,img_fn, '_input.mat'), 'mworks')
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    