clear; close all;
analysis_out = 'A:\home\carlo\analysis\2P\';
bdata_source = 'A:\home\carlo\rawData\behavioral\';

removing_bad_trials = false; 

expType=input('monomodal(1) or dimodal(2) session: ');
seshType=input('interleaved (1), blocked (2), or generalization (3) session: ');
if seshType==1 && expType == 1
RCExptListCarlo_Inter;
elseif seshType==2 && expType == 1
RCExptListCarlo_Block;
elseif seshType==1 && expType == 2
RCExptListCarlo_VA_Inter;
elseif seshType==2 && expType == 2
RCExptListCarlo_VA_Block;
elseif seshType==3 && expType == 2
RCExptListCarlo_VA_Inter;
elseif seshType==3 && expType == 1
RCExptListCarlo_Gen;
elseif seshType==2 && expType == 3
RCExptListCarlo_CSPlus_Block;
end

% RCExptListCarlo_Inter
% RCExptListCarlo_Block
pairings = loadRCList_Carlo(expt,seshType);

%for  j = 1:size(expt.ttl,2)
thresholdDeco = input('Which deconvolution threshold do you wish to plot?: ');


    id = pairings{1,1};
    exp_subset = id;
    iexp = exp_subset; isesh = pairings{3,1};
mouse = strtrim(expt{1,iexp}.mouse{1,isesh});
date = expt{1,iexp}.date{1,isesh};
run = expt{1,iexp}.run{1,isesh};
    fprintf([date ' ' mouse ' ' run '\n']);
    if isempty(thresholdDeco)
    img_fn = [date '_img' mouse '\getTC_' run '\'];
    else
    img_fn = [date '_img' mouse '\getTC_' run '_' num2str(thresholdDeco) '\'];    
    end
    if ~exist(fullfile(analysis_out,img_fn))
    mkdir(fullfile(analysis_out, img_fn))
    end
    filename = dir([analysis_out,date '_img' mouse '\getTC_' run '\'  '*' '_TCave.mat']);
    load([analysis_out,date '_img' mouse '\getTC_' run '\', filename.name]);
    filename = dir([analysis_out,date '_img' mouse '\' date '_img' mouse '_' run  '*' num2str(thresholdDeco) '_TCave_cl.mat']);
    load([analysis_out,date '_img' mouse '\', filename.name]);
    %threshold = -3;
    spk = load([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_spk_deconvolve_threshold' num2str(thresholdDeco) '.mat' ]);
    spk_logic = spk.spk_logic_cl;
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
    
    nIC = size(spk_logic,2);
    react_time = double(cell2mat(mworks.reactTimesMs(1,2:end)));
    cue_rew_int = mworks.RewardDelayDurationMs + round(mean(react_time),-1); %"total interval between cueonset and reward delivery"       
    img_fn2 = sessions;
    frameRateHz = double(mworks.frameRateHz);
    cTargetOn = (mworks.cTargetOn);%(1,2:length(cTargetOn_cutted)+1));
    nTrials = length(cTargetOn_cutted)+1;
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
%     tempFig=setFigure; hold on; plot(nanmean(all_TCave_cl,2),'r'); hold on; plot(nanmean(raw_tc_avg,2),'k'); plot(nanmean(all_TCave_clProc,2),'g');
    %process should be obvious from figure (shift of the red trace to the level of the black trace to 
    %give the green trace, which will be used for aligning to cue presentation - important note:
    %since raw_tc_avg is formatted to include only laseron periods, it is SHORTER than the all_TCave 
    %traces, this is expected and should make sense, focus shouldn't be on the length match, 
    %but the match at the y-axis
    if ~isnan(cTargetOn(1,end))
    amended_cTargetOn = cTargetOn(1,sum(isnan(cTargetOn)):end);
    else
    end
    for itrial = 1:nTrials
        if cTargetOn(itrial)+postwin_frames-1 <= nFrames %& input.counterValues{itrial}(end)-cTargetOn(itrial) > postwin_frames
            targetAlign_tc(:,:,itrial) = all_TCave_cl(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:);
            targetAlign_events(:,:,itrial) = all_events(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:);
        end
    end
    
   
    targetAlignF = nanmean(nanmean(targetAlign_tc(1:prewin_frames,:,:),3),1); % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.
    targetAligndFoverF = zeros(size(targetAlign_tc,1),size(targetAlign_tc,2),size(targetAlign_tc,3));%frame*cell*trials
    % calculate df/F
    for c = 1:size(TCave_cl,2)
        targetAligndFoverF(:,c,:) = (targetAlign_tc(:,c,:)-targetAlignF(1,c))./targetAlignF(1,c);
    end
    
    
        omitRewTrial = celleqel2mat_padded(mworks.tRewardOmissionTrial(1,2:end));
        ind_omit = find(omitRewTrial);
        
        unexpRewTrial = celleqel2mat_padded(mworks.tDoNoStimulusChange(1,2:end));
        ind_unexp = find(unexpRewTrial);
        
        if mworks.doBlock2 %Specific CS+/CS- experiment added by Mike
            block2Trial = celleqel2mat_padded(mworks.tBlock2TrialNumber(1,1:length(cTargetOn_cutted)+1));%length(cTargetOn)));
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
        reactTimesMs = double(mean(celleqel2mat_padded(mworks.reactTimesMs),2));
        delayTimeMs = reactTimesMs+rewardDelayDurationMs;
        
        save([analysis_out, img_fn, 'figDataUntrimmed_' sessions  '.mat']);
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
        tempFig=setFigure; hold on; n = 1;
        subplot(s,1,n); hold on;
        shadedErrorBar_CV(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_rew),3),2), std(nanmean(targetAligndFoverF(:,:,ind_rew),3),[],2,'omitnan')./sqrt(nIC),'lineProps','k');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        %ylim([-3 5]);
        title('CS+')
        n = n+1;
        if length(ind_omit>5)
        subplot(s,1,n); hold on;
        shadedErrorBar_CV(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_omit),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_omit),3),[],2)./sqrt(nIC),'lineProps','r');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        title('Omit')
        n = n+1;
        end
        if length(ind_unexp>5)
        subplot(s,1,n); hold on;
        shadedErrorBar_CV(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_unexp),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_unexp),3),[],2)./sqrt(nIC),'g');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        title('Unexpected reward')
        n = n+1;
        end
        if mworks.doBlock2
        subplot(s,1,n); hold on;
        shadedErrorBar_CV(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_block2),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_block2),3),[],2)./sqrt(nIC),'lineProps','r');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        %ylim([-3 5]);
        title('CS-')
        end
        sgtitle([date ' ' mouse])
        savefig(fullfile(analysis_out,img_fn, [ '_cueAlign_dFoverF.fig']))
        saveas(tempFig, [analysis_out img_fn  '_cueAlign_dFoverF.pdf']);
        
        
        tempFig=setFigure; hold on; n=1;
        subplot(s,1,n); hold on;
        shadedErrorBar_CV(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew),3),[],2)./sqrt(nIC)).*(1000./frameRateHz));
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([0 5]);
        title('CS+')
        n = n+1;
        if length(ind_omit>5)
        subplot(s,1,n); hold on;
        shadedErrorBar_CV(tt, nanmean(nanmean(targetAlign_events(:,:,ind_omit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','b');
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Omit')
        n = n+1;
        end
        if length(ind_unexp>5)
        subplot(s,1,n); hold on;
        shadedErrorBar_CV(tt, nanmean(nanmean(targetAlign_events(:,:,ind_unexp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_unexp),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'g');
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Unexpected reward')
        n = n+1;
        end
        if mworks.doBlock2
        subplot(s,1,n); hold on;
        shadedErrorBar_CV(tt, nanmean(nanmean(targetAlign_events(:,:,ind_block2),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_block2),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','r');
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([0 5]);
        title('CS-')
        end
        sgtitle([date ' ' mouse])
        savefig(fullfile(analysis_out,img_fn, [ '_cueAlign_events_Hz.fig']))
        saveas(tempFig, [analysis_out img_fn '_cueAlign_events_Hz.pdf']);
        
        elseif indivFig
            
        tempFig=setFigure; hold on; %n = 1;
        %subplot(s,1,n) 
        if expt{exp_subset}.whichBlock(1,isesh) == 0
            shadedErrorBar_CV(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_rew),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_rew),3),[],2)./sqrt(nIC),'lineProps','r');
            vline(cue_rew_int)        
            xlabel('Time from cue')
            ylabel('Avg dF/F')
            title('CS-')
        elseif expt{exp_subset}.whichBlock(1,isesh) == 1
            shadedErrorBar_CV(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_rew),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_rew),3),[],2)./sqrt(nIC),'lineProps','k');
            vline(cue_rew_int)        
            xlabel('Time from cue')
            ylabel('Avg dF/F')
            title('CS+')
        end
        sgtitle([date ' ' mouse])
        savefig(fullfile(analysis_out,img_fn, [ '_cueAlign_dFoverF.fig']))
        saveas(tempFig, [analysis_out img_fn  '_cueAlign_dFoverF.pdf']);
        
        tempFig=setFigure; hold on; %n=1;
        %subplot(s,1,n)
        if expt{exp_subset}.whichBlock(1,isesh) == 0
            shadedErrorBar_CV(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','r');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([0 3]);
            vline(cue_rew_int,'k')
            vline(0,'k')            
            title('CS-')
        elseif expt{exp_subset}.whichBlock(1,isesh) == 1
            shadedErrorBar_CV(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','k');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([0 3]);
            vline(cue_rew_int,'k')
            vline(0,'k')
            title('CS+')
        end
        sgtitle([date ' ' mouse])
        savefig(fullfile(analysis_out,img_fn,[ '_cueAlign_events_Hz.fig']))
        saveas(tempFig, [analysis_out img_fn  '_cueAlign_events_Hz.pdf']);
        end
        
        
        if length(ind_omit>10)
            tempFig=setFigure; hold on;
            subplot(2,1,1)
            thresh = mean(diff(ind_omit));
            ind_omit_short = find(diff(ind_omit)<thresh)+1;
            ind_omit_long = find(diff(ind_omit)>thresh)+1;
            shadedErrorBar_CV(tt, nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_omit_short)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_omit_short)),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','k');
            hold on
            shadedErrorBar_CV(tt, nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_omit_long)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_omit_long)),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','b');
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
            shadedErrorBar_CV(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew_preomit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_preomit),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','k');
            hold on
            shadedErrorBar_CV(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew_postomit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_postomit),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','b');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Reward trial pre (black n = ' num2str(length(ind_rew_preomit)) ') vs post (blue n = ' num2str(length(ind_rew_postomit)) ') omit trials'])
            sgtitle([date ' ' mouse])
            savefig(fullfile(analysis_out,img_fn, [img_fn  '_omitByInterval.fig']))
        else
            ind_omit_short = [];
            ind_omit_long = [];
            ind_rew_preomit = [];
            ind_rew_postomit = [];
        end
        if length(ind_unexp>10)
            tempFig=setFigure; hold on;
            thresh = mean(diff(ind_unexp));
            ind_unexp_short = find(diff(ind_unexp)<thresh)+1;
            ind_unexp_long = find(diff(ind_unexp)>thresh)+1;
            shadedErrorBar_CV(tt, nanmean(nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_short)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_short)),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','k');
            hold on
            shadedErrorBar_CV(tt, nanmean(nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_long)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_long)),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','b');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title([date ' ' mouse '- Short (black n = ' num2str(length(ind_unexp_short)) ') vs long (blue n = ' num2str(length(ind_unexp_long)) ') interval between unexpected reward'])
            savefig(fullfile(analysis_out,img_fn, [img_fn  '_unexpByInterval.fig']))
        else
            ind_unexp_short = [];
            ind_unexp_long = [];
        end
        
        if unique(expt{exp_subset}.name == 'Inter')
            
         tempFig=setFigure; hold on; n=1;
        %subplot(s,1,n)
        plot(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz),'k');
        hold on;
        plot(tt, nanmean(nanmean(targetAlign_events(:,:,ind_block2),3),2).*(1000./frameRateHz),'r');
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        %ylim([0 6]);
        title('Stacked spike rates separating CS+ and CS- trials');
        legend('CS+','CS-');
        sgtitle([date ' ' mouse])
        savefig(fullfile(analysis_out,img_fn, [ '_line_stackedSpike_events_Hz.fig']))
        saveas(tempFig, [analysis_out img_fn  '_line_stackedSpike_events_Hz.pdf']);

         tempFig=setFigure; hold on; n=1;
        %subplot(s,1,n)
        shadedErrorBar_CV(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','k');
        hold on;
        shadedErrorBar_CV(tt, nanmean(nanmean(targetAlign_events(:,:,ind_block2),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_block2),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','r');
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        %ylim([0 6]);
        title('Stacked spike rates separating CS+ and CS- trials');
        %legend('CS+', 'sd CS+' ,'CS-', 'sd CS-');
        sgtitle([date ' ' mouse])
        savefig(fullfile(analysis_out,img_fn, [ '_error_stackedSpike_events_Hz.fig']))
        saveas(tempFig, [analysis_out img_fn  '_error_stackedSpike_events_Hz.pdf']);

        end 
        
        if unique(expt{exp_subset}.name == 'Inter')
        tempFig=setFigure; hold on;
        n = floor(nTrials./3);
        seshType = 1;
        for i = 1:3
            subplot(3,2,seshType); hold on;
            ind_rew_temp = intersect(ind_rew,1+(i-1).*n:i*n);
            shadedErrorBar_CV(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew_temp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_temp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
            hold on;
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([0 7]);
            title(['Trials ' num2str(1+(i-1).*n) ':' num2str(i*n)])
            vline(cue_rew_int)
            if length(ind_block2)>=10
                subplot(3,2,seshType+1); hold on;
                ind_block2_temp = intersect(ind_block2,1+(i-1).*n:i*n);
                shadedErrorBar_CV(tt,nanmean(nanmean(targetAlign_events(:,:,ind_block2_temp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_block2_temp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
                ylim([0 7])
                xlabel('Time from cue')
                ylabel('Spike rate (Hz)')
                title(['Trials ' num2str(1+(i-1).*n) ':' num2str(i*n)])
                vline(cue_rew_int)
            end
            seshType = seshType+2;
        end
        sgtitle([date ' ' mouse '- CS+ (black), CS- (red)'])
        savefig(fullfile(analysis_out,img_fn, [ '_repsByTrial.fig']))
        saveas(tempFig, [analysis_out img_fn  '_repsByTrial.pdf']);

        elseif unique(expt{exp_subset}.name == 'Block')
            if expt{exp_subset}.whichBlock(1,isesh) == 0
        tempFig=setFigure; hold on;
        n = floor(nTrials./3);
        seshType = 1;
        for i = 1:3
            subplot(3,2,seshType); hold on;
            ind_rew_temp = intersect(ind_rew,1+(i-1).*n:i*n);
            shadedErrorBar_CV(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew_temp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_temp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'r');
            hold on
            ylim([0 7])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Trials ' num2str(1+(i-1).*n) ':' num2str(i*n)])
            vline(cue_rew_int)
           
            seshType=seshType+2;
        end 
        sgtitle([date ' ' mouse '- CS- (red)'])
        savefig(fullfile(analysis_out,img_fn, [ '_repsByTrial.fig']))
        saveas(tempFig, [analysis_out img_fn  '_repsByTrial.pdf']);

            elseif expt{exp_subset}.whichBlock(1,isesh) == 1
        tempFig=setFigure; hold on;
        n = floor(nTrials./3);
        seshType = 1;
        for i = 1:3
            subplot(3,2,seshType); hold on;
            ind_rew_temp = intersect(ind_rew,1+(i-1).*n:i*n);
            shadedErrorBar_CV(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew_temp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_temp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
            hold on
            ylim([0 7])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Trials ' num2str(1+(i-1).*n) ':' num2str(i*n)])
            vline(cue_rew_int)
        
            seshType=seshType+2;
        end
        sgtitle([date ' ' mouse '- CS+ (black)'])
        savefig(fullfile(analysis_out,img_fn, [ '_repsByTrial.fig']))
        saveas(tempFig, [analysis_out img_fn  '_repsByTrial.pdf']);

            end
        end
     
        save(fullfile(analysis_out,img_fn, [ '_targetAlign.mat']), 'ind_rew', 'ind_block2', 'ind_omit', 'ind_unexp', 'targetAlign_events', 'targetAligndFoverF', 'prewin_frames', 'postwin_frames', 'tt', 'frameRateHz', 'ind_omit_short','ind_omit_long','ind_unexp_short','ind_unexp_long','ind_rew_preomit','ind_rew_postomit')
        save(fullfile(analysis_out,img_fn, [ '_input.mat']), 'mworks')
    
 %%   
load([analysis_out,date '_img' mouse '\getTC_' run '\' date '_img' mouse '_' run, 'saveOutputs.mat']);
filename2 = dir([analysis_out,date '_img' mouse '\getTC_' run '\' '*' 'thresh' num2str(cluster_threshold) '_coor' num2str(threshold) '_mask3D.mat']);
mask3D = load([analysis_out date '_img' mouse '\getTC_' run '\' filename2.name]);
mask3D = mask3D.mask3D;

 %peak spike rate, trial averaged, for each neuron across CS+ and CS- trials
%%this has to be trial averaged because the 'spike rate' is found by averaging 
%%the number of selected events across all trials (i.e. spike rate of 4 at
%%a given time point means the neuron was active nTimes/nTrials = 4
for i = 1:length(goodcells)
    block2PeakAmp(i,1) = max(nanmean(targetAlign_events(:,i,logical(block2Trial')),3).*(1000./frameRateHz));
    rewPeakAmp(i,1) = max(nanmean(targetAlign_events(:,i,~block2Trial'),3).*(1000./frameRateHz));
end
%plot difference map simply on b2R<rewR and vice versa
for i = 1:length(goodcells)
    if block2PeakAmp(i,1)>rewPeakAmp(i,1)
        splitPreference(i,1) = 1;
    else
        splitPreference(i,1) = 0;
    end
end
tempFig=setFigure; hold on; plot(block2PeakAmp(splitPreference==1,1),rewPeakAmp(splitPreference==1,1),'rx');
hold on; plot(block2PeakAmp(splitPreference==0,1),rewPeakAmp(splitPreference==0,1),'gx'); hold on;

 tempFig=setFigure; hold on;
%needs to modify the line below due to different file names
imshow([analysis_out,date '_img' mouse '\getTC_' run '\' 'AVG_' date '_img' mouse '_' run '_rgstr_tiff_every50_ref' num2str(ref) '_jpeg.jpg']); hold on;
for m  = 1:size(targetAlign_events,2)
    thisCell = goodcells(m);
    if splitPreference(m)==1
    bound = cell2mat(bwboundaries(mask3D(:,:,thisCell)));
    randcolor = [0.4660, 0.6740, 0.1880];
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;
    elseif splitPreference(m)==0
    bound = cell2mat(bwboundaries(mask3D(:,:,thisCell)));
    randcolor = [0.6350, 0.0780, 0.1840];
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;    
    end
end
hold on;
title([date ' img' mouse ' masks with split preference for CS+ and CS-']);

%take just those neurons with the greatest difference (top and 'bot' 25%)
[b2DiffPeakAmp, b2DiffPeakAmpIdx] = maxk(block2PeakAmp-rewPeakAmp,round(.25*size(targetAlign_events,2)));
[rewDiffPeakAmp,rewDiffPeakAmpIdx] = maxk(rewPeakAmp-block2PeakAmp,round(.25*size(targetAlign_events,2)));
tempFig=setFigure; hold on;
plot(block2PeakAmp,rewPeakAmp,'ko'); hold on;
xlabel('block2 amplitude')
ylabel('reward amplitude')
plot(block2PeakAmp(b2DiffPeakAmpIdx,1),rewPeakAmp(b2DiffPeakAmpIdx,1),'ro'); hold on;
plot(block2PeakAmp(rewDiffPeakAmpIdx,1),rewPeakAmp(rewDiffPeakAmpIdx,1),'go');

tempFig=setFigure; hold on;
%needs to modify the line below due to different file names
imshow([analysis_out,date '_img' mouse '\getTC_' run '\' 'AVG_' date '_img' mouse '_' run '_rgstr_tiff_every50_ref' num2str(ref) '_jpeg.jpg']); hold on;
for m  = 1:size(targetAlign_events,2)
        thisCell = goodcells(m);
    if sum(b2DiffPeakAmpIdx==m)==1
    bound = cell2mat(bwboundaries(mask3D(:,:,thisCell)));
    randcolor = [0.4660, 0.6740, 0.1880];
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;
    elseif sum(rewDiffPeakAmpIdx==m)==1
    bound = cell2mat(bwboundaries(mask3D(:,:,thisCell)));
    randcolor = [0.6350, 0.0780, 0.1840];
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;    
%     elseif sum(rewDiffPeakAmpIdx==m)~=1 && sum(b2DiffPeakAmpIdx==m)~=1
%     bound = cell2mat(bwboundaries(mask3D(:,:,m)));
%     randcolor = [0, 0.4470, 0.7410];
%     plot(bound(:,2),bound(:,1),'.','color','k'); hold on;    
     end
end
hold on;
title([date ' img' mouse ' masks of the top 25% of neurons with split preference for CS+ and CS-']);

%{
    Use dFoF, take a 100ms bin around the peak activation after cue
    presentation. Find a 1xcell for CS+ trials and CS- trials. Can then run
    a diff between 1xcell for CS+ and 1xcell for  CS- and plot a mask with
    a gradient for the cell-dependent difference between CS+/- response.

    bin 200:400ms after cue - bin 19:21 of dfofBinWFrames
%}
%
for i = 1:length(goodcells)
    for j = 1:nTrials
        binnedDFoFtargetAlign(:,i,j)=targetAligndFoverF(57:63,i,j);     
    end
end

for i = 1:length(goodcells)
    for j = 1:nTrials
        base(i,j) = min(binnedDFoFtargetAlign(:,i,j));
        absBinDFoFtargetAlign(:,i,j) = (binnedDFoFtargetAlign(:,i,j) + abs(base(i,j)));
    end
end

bSum = sum(block2Trial); rSum = sum(~block2Trial); if bSum>rSum; tSum=rSum; elseif bSum<rSum; tSum=bSum; elseif bSum==rSum; tSum=bSum; end
peakAmpB2 = nan(1,length(goodcells),tSum); peakAmpRew = nan(1,length(goodcells),tSum); 
for i = 1:length(goodcells)
    cB = 1; cR = 1;
    for j = 1:length(block2Trial)
        if block2Trial(1,j) == 1 & j+1<=tSum
        peakAmpB2(1,i,cB) = max(absBinDFoFtargetAlign(:,i,j+1));
        cB=cB+1;
        elseif block2Trial(1,j) == 0 & j+1<=tSum
        peakAmpRew(1,i,cR) = max(absBinDFoFtargetAlign(:,i,j+1));
        cR = cR+1;
        end
    end
end

for i = 1:length(goodcells)
    [h(i,:),p(i,:),sd(i,:)] = ttest(peakAmpB2(:,i,:),peakAmpRew(:,i,:));
end

uPeakAmpB2 = nanmean(peakAmpB2,3);
uPeakAmpRew = nanmean(peakAmpRew,3);
%%
[b2DiffPeakAmp, b2DiffPeakAmpIdx] = maxk(uPeakAmpB2-uPeakAmpRew,10);
[rewDiffPeakAmp,rewDiffPeakAmpIdx] = maxk(uPeakAmpRew-uPeakAmpB2,10);
b2Diff = zeros(length(goodcells),1);
rewDiff = zeros(length(goodcells),1);
b2Diff(b2DiffPeakAmpIdx)=1;
rewDiff(rewDiffPeakAmpIdx)=1;
b2Diff = logical(b2Diff);
rewDiff = logical(rewDiff);
[b2DiffPeak, b2DiffPeakIdx] = maxk(block2PeakAmp-rewPeakAmp,10);
[rewDiffPeak,rewDiffPeakIdx] = maxk(rewPeakAmp-block2PeakAmp,10);

tempFig=setFigure; hold on;
plot(uPeakAmpB2,uPeakAmpRew,'ko'); hold on;
xlabel('block2 amplitude')
ylabel('reward amplitude')
plot(uPeakAmpB2((b2Diff & h==1)'),uPeakAmpRew((b2Diff & h==1)'),'ro'); hold on;
plot(uPeakAmpB2((rewDiff & h==1)'),uPeakAmpRew((rewDiff & h==1)'),'go');
title('dFoF on a per neuron basis between CS+ and CS- trials')

tempFig=setFigure; hold on;
plot(log(uPeakAmpB2),log(uPeakAmpRew),'ko'); hold on;
xlabel('block2 amplitude')
ylabel('reward amplitude')
plot(log(uPeakAmpB2((b2Diff & h==1)')),log(uPeakAmpRew((b2Diff & h==1)')),'ro'); hold on;
plot(log(uPeakAmpB2((rewDiff & h==1)')),log(uPeakAmpRew((rewDiff & h==1)')),'go');
title('log of dFoF on a per neuron basis between CS+ and CS- trials')

tempFig=setFigure; hold on;
plot(block2PeakAmp,rewPeakAmp,'ko'); hold on;
xlabel('block2 amplitude')
ylabel('reward amplitude')
plot(block2PeakAmp(b2DiffPeakIdx,1),rewPeakAmp(b2DiffPeakIdx,1),'ro'); hold on;
plot(block2PeakAmp(rewDiffPeakIdx,1),rewPeakAmp(rewDiffPeakIdx,1),'go');
title('spike rate on a per neuron basis between CS+ and CS- trials.');
supertitle('significant preference selected based on largest diff between CS+/CS- spike rates in either direction for each neuron');


eventsB2 = [nanmean(targetAlign_events(:,:,logical(block2Trial)),3)];
eventsRew = [nanmean(targetAlign_events(:,:,~block2Trial),3)];
        
        
tt=tt(1,:);
tempFig=setFigure; hold on;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,nanmean(eventsB2(:,b2DiffPeakIdx),2).*(1000./frameRateHz),(nanstd(eventsB2(:,b2DiffPeakIdx),[],2)./sqrt(length(goodcells))).*(1000./frameRateHz),'lineProps','r'); hold on;
shadedErrorBar_CV(tt,nanmean(eventsRew(:,b2DiffPeakIdx),2).*(1000./frameRateHz),(nanstd(eventsRew(:,b2DiffPeakIdx),[],2)./sqrt(length(goodcells))).*(1000./frameRateHz),'lineProps','k');
ylim([0 10])
xlabel('Time from cue')
ylabel('Spike rate (Hz)')
vline(767,'k')
vline(0,'k')
title('spike rate from 10 cells with the largest spike differential in favor of CS- trials');
          
subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,nanmean(eventsB2(:,rewDiffPeakIdx),2).*(1000./frameRateHz),(nanstd(eventsB2(:,rewDiffPeakIdx),[],2)./sqrt(length(goodcells))).*(1000./frameRateHz),'lineProps','r'); hold on;
shadedErrorBar_CV(tt,nanmean(eventsRew(:,rewDiffPeakIdx),2).*(1000./frameRateHz),(nanstd(eventsRew(:,rewDiffPeakIdx),[],2)./sqrt(length(goodcells))).*(1000./frameRateHz),'lineProps','k');
ylim([0 10])
xlabel('Time from cue')
ylabel('Spike rate (Hz)')
vline(767,'k')
vline(0,'k')
title('spike rate from 10 cells with the largest spike differential in favor of CS+ trials');



tempFig=setFigure; hold on;
%needs to modify the line below due to different file names
imshow([analysis_out,date '_img' mouse '\getTC_' run '\' 'AVG_' date '_img' mouse '_' run '_rgstr_tiff_every50_ref' num2str(ref) '_jpeg.jpg']); hold on;
for m  = 1:size(targetAlign_events,2)
        thisCell = goodcells(m);
    if sum(m==b2DiffPeakIdx)==1
    bound = cell2mat(bwboundaries(mask3D(:,:,thisCell)));
    randcolor = [0.4660, 0.6740, 0.1880];
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;
    elseif sum(m==rewDiffPeakIdx)==1
    bound = cell2mat(bwboundaries(mask3D(:,:,thisCell)));
    randcolor = [0.6350, 0.0780, 0.1840];
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;    
    elseif sum(m==b2DiffPeakIdx)==0 && sum(m==rewDiffPeakIdx)==0
    bound = cell2mat(bwboundaries(mask3D(:,:,thisCell)));
    randcolor = [0, 0.4470, 0.7410];
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;
    end
end
hold on;
title([date ' img' mouse ' masks with preference for CS+ and CS-']);

    b2PrefLog = (b2Diff & h==1)';
    rewPrefLog = (rewDiff & h==1)';
    mean(block2PeakAmp(b2DiffPeakIdx)-rewPeakAmp(b2DiffPeakIdx));
    mean(rewPeakAmp(rewDiffPeakIdx)-block2PeakAmp(rewDiffPeakIdx));

    tempFig=setFigure; hold on;
%needs to modify the line below due to different file names
imshow([analysis_out,date '_img' mouse '\getTC_' run '\' 'AVG_' date '_img' mouse '_' run '_rgstr_tiff_every50_ref' num2str(ref) '_jpeg.jpg']); hold on;
for m  = 1:size(targetAlign_events,2)
        thisCell = goodcells(m);
    if b2PrefLog(m)==1
    bound = cell2mat(bwboundaries(mask3D(:,:,thisCell)));
    randcolor = [0.4660, 0.6740, 0.1880];
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;
    elseif rewPrefLog(m)==1
    bound = cell2mat(bwboundaries(mask3D(:,:,thisCell)));
    randcolor = [0.6350, 0.0780, 0.1840];
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;    
    elseif b2PrefLog(m)==0 && rewPrefLog(m)==0
    bound = cell2mat(bwboundaries(mask3D(:,:,thisCell)));
    randcolor = [0, 0.4470, 0.7410];
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;
    end
end
hold on;
title([date ' img' mouse ' masks with preference for CS+ and CS-']);

    
    


  