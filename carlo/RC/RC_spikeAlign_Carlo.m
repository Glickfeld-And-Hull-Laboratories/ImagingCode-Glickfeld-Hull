clear all
plotDetect = 0;
negAdjust = 0;
do_decon = 1;
spike_rate_thresh = 0.3;

carlo_out = 'A:\home\carlo\analysis\2P';
bdata_source ='A:\home\carlo\rawData\behavioral\';
removing_bad_trials = false; %bad_trials = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]; %added to remove shitty trials, turn off if not needed
%bad_cells_manual = [1,2,3,4,7,9,13,38,41,43,45,46,47,50,51,52,55];

RCExptListCarlo;
pairings = loadRCList_Carlo();

for  j = 1:size(pairings,2)
id = pairings{1,j};
exp_subset = pairings{2,j};

nexp = size(expt(id).date,1);
   
    
    for iexp = exp_subset  %1:nexp
        mouse = strtrim(expt(id).mouse(iexp,:));
        date = expt(id).date(iexp,:);
        run = expt(id).run(iexp,:);
        fprintf([date ' ' mouse '\n'])
        img_fn = [date '_img' mouse];
        load(fullfile(carlo_out,img_fn, ['getTC_' run '\' img_fn '_' run '_nPCA200_mu0.3_nIC120_thresh97_TCave.mat'])) %mask_raw, mask_flat, mask_final, data_corr
        load(fullfile(carlo_out,img_fn, ['getTC_' run '\' img_fn '_' run '_nPCA200_mu0.3_nIC120_thresh97_coor0.8_mask3D_final.mat']))
        load(fullfile(carlo_out,img_fn, [img_fn '_' run '_deconvolution_thresh-3_TCave_cl.mat']))
        load(fullfile(carlo_out,img_fn, [img_fn '_' run '_spk_deconvolve_threshold-3.mat']))
        %tc_avg(:,bad_cells_manual) = [];
        nIC = size(tc_avg,2);

if strlength(img_fn) < 16
    bfile = dir([bdata_source 'data-i' mouse '-' date '*' ]);
else
    bfile = dir([bdata_source 'data-i' mouse '-' date img_fn(end-4:end) '.mat']);
end
%load behavior file
behav_dest = [bdata_source bfile.name];
assert(length(bfile)) = 1;
b_data = load(behav_dest);
mworks = b_data.mworks;
cTargetOn = mworks.cTargetOn;
frameRateHz = double(mworks.frameRateHz);

        react_time = double(cell2mat(mworks.reactTimesMs));
        cue_rew_int = mworks.RewardDelayDurationMs + round(mean(react_time),-1); %"total interval between cueonset and reward delivery"
        
%%%%%%%%%%%%%%%%%%%        
        opt.thresh = 1.7;
        opt.normalization = 1;
        opt.deconvtau = 0;
        opt.dt = 1;

%         if str2num(mouse)>800
%             mouse_name = mouse(2:3);
%         else
%             mouse_name = mouse;
%         end
%         
        if str2num(mouse)>800 & str2num(mouse) < 1000
            mouse_name = mouse(2:3);
        elseif str2num(mouse)<100 | str2num(mouse) > 1000
            mouse_name = mouse;
        end
        %% correct for skipped frames
        fns = dir(fullfile(bdata_source, ['data-i' mouse '-' date '*']));
        if length(fns)>1
            nfns{id}(iexp) = length(fns);
            fn(id).tr{iexp}= [];
            for i = 1:length(fns)
                load(fullfile(bdata_source,fns(i).name));
                fn(id).tr{iexp} = [fn(id).tr{iexp} size(mworks.gratingContrast,2)];
            end
            fns = fns(end); %assumes the final behavioral file collected is the correct one
        else
            nfns{id}(iexp) = 1;  
        end
        load(fullfile(bdata_source,fns.name));  %laods bx data
        react_time = double(cell2mat(mworks.reactTimesMs));
        cue_rew_int = mworks.RewardDelayDurationMs + round(mean(react_time),-1); %"total interval between cueonset and reward delivery"
        
        nTrials = length(mworks.cTargetOn);
        mworks.counterValues_orig = mworks.counterValues; 
        mworks.counterTimesUs_orig = mworks.counterTimesUs;
        mworks.cTargetOn_orig = mworks.cTargetOn;
        val_reset = 0;
        for it = 1:nTrials
            diff_val = diff(mworks.counterValues{it});
            if length(unique(diff_val)) == 1
                continue
            elseif length(unique(diff_val)) == 2 & sum(ismember(unique(diff_val),2)) == 1
                val_reset = 1;
                ind = find(diff_val== 2);
                for i = 1:length(ind)
                    d = mworks.counterTimesUs{it}(ind(i))-mworks.counterTimesUs{it}(ind-1);
                    if d>50000 & d<70000
                        val_reset = 1;
                        mworks.counterValues{it} = [mworks.counterValues{it}(1:ind(i)-1) mworks.counterValues{it}(1:ind(i)-1)+1 mworks.counterValues{it}(ind(i):end)];
                        mworks.counterValues{it} = [mworks.counterTimesUs{it}(1:ind(i)-1) mworks.counterTimesUs{it}(1:ind(i)-1)+33000 mworks.counterTimesUs{it}(ind(i):end)];
                        ind = ind+1;
                    end
                end
            elseif length(unique(diff_val)) > 2
                val_reset = 1;
                ind = find(diff(mworks.counterTimesUs{it}./1000)<25);
                ind = ind(1)-2:ind(end)+2;
                timeDiff = mworks.counterTimesUs{it}(ind(end))-mworks.counterTimesUs{it}(ind(1));
                frameDiff = mworks.counterValues{it}(ind(end))-mworks.counterValues{it}(ind(1));
                nf = round((timeDiff./1000)*frameRateHz/1000);
                mworks.counterValues{it} = [mworks.counterValues{it}(1:ind(1)-1) bsxfun(@plus,mworks.counterValues{it}(ind(1)-1), 1:nf) [mworks.counterValues{it}(ind(end):end)-frameDiff+nf]];
                mworks.counterTimesUs{it} = [mworks.counterTimesUs{it}(1:ind(1)-1) [repmat(mworks.counterTimesUs{it}(ind(1)-1),[1 nf]) + (33000:33000:33000*nf)] mworks.counterTimesUs{it}(ind(end):end)];
                if mworks.cTargetOn{it}>ind(end)
                    mworks.cTargetOn{it} = mworks.cTargetOn{it}-frameDiff+nf;
                elseif mworks.cTargetOn{it}>ind(1)
                    mworks.cTargetOn{it} = NaN;
                end
                for it2 = it+1:nTrials
                    mworks.counterValues{it2} = bsxfun(@plus, mworks.counterValues_orig{it2},-frameDiff+nf);
                    mworks.cTargetOn{it2} = bsxfun(@plus, mworks.cTargetOn_orig{it2}, -frameDiff+nf);
                end
                diff_val = diff(mworks.counterValues{it});
                if length(unique(diff_val)) == 1
                    continue
                elseif length(unique(diff_val)) == 2 & sum(ismember(unique(diff_val),2))
                    ind = find(diff_val== 2);
                    for i = 1:length(ind)
                        d = mworks.counterTimesUs{it}(ind(i))-mworks.counterTimesUs{it}(ind-1);
                        if d>50000 & d<70000
                            val_reset = 1;
                            mworks.counterValues{it} = [mworks.counterValues{it}(1:ind(i)-1) mworks.counterValues{it}(1:ind(i)-1)+1 mworks.counterValues{it}(ind(i):end)];
                            mworks.counterValues{it} = [mworks.counterTimesUs{it}(1:ind(i)-1) mworks.counterTimesUs{it}(1:ind(i)-1)+33000 mworks.counterTimesUs{it}(ind(i):end)];
                            ind = ind+1;
                        end
                    end
                end
            end
        end
        if val_reset == 0
            rmfield(mworks,{'counterValues_orig','counterTimesUs_orig','cTargetOn'});
            fprintf(' counter not corrected \n')
        else
            fprintf(' corrected counter \n')
        end

        %% detect events

        tc_avg_temp = tc_avg;
        cTargetOn = celleqel2mat_padded(mworks.cTargetOn);
        if expt(id).ttl(iexp)
            cd(fullfile('H:\home\mike\Data\2P Imaging', [date '_img' mouse_name],['img' mouse_name]))
            load(['img' mouse_name '_000_' run '_realtime.mat'])
            
            %if ttl_log is faulty then identify pockel transitions via changes in F
            if length(unique(ttl_log)) == 1 | sum(ttl_log) < 0.05*(length(ttl_log))
                clear ttl_log
                ttl_log = find_pockel_tc(tc_avg, 1)';
            end
            %Below commented out by Mike on 12/9/19, seems to find pockel transition frames in
            %an inconsistent manner with the find_pockel_tc function
            %ttl_log = ttl_log(1:size(tc_avg,1),:);
            %ttl_trans = find(abs(diff(ttl_log)));      
            %n = length(ttl_trans);
            %for i = 1:n
            %    ttl_log(ttl_trans(i)-4:ttl_trans(i)+4,:) = 0;
            %end
            ind_ttl = intersect(1:size(tc_avg,1), find(ttl_log == 0));
            tc_avg(ind_ttl,:) = NaN;
            ind_tc = find(~isnan(tc_avg(:,1)));
            tc_avg_temp(find(isnan(tc_avg(:,1))),:) = [];
        end
        range = 5501:5900;
        
        tc_diff = diff(tc_avg_temp,[],1);
        all_events = zeros(size(tc_avg_temp));
        if plotDetect
            figure;
            start = 1;
            IC_list = [];
        end
        for ic = 1:nIC
            if do_decon; continue; end
            [~, iti_ind, ~] = CellsortFindspikes(tc_diff(:,ic), opt.thresh, opt.dt, opt.deconvtau, opt.normalization);
            iti_ind(find(diff(iti_ind)==1)+1)=[];
            all_events(iti_ind+1,ic) = 1;
            if negAdjust
                ind_events = find(all_events(:,ic));
                for i = 1:length(ind_events)
                    if ind_events(i)< 50
                        if tc_avg_temp(ind_events(i)-1,ic) < mean(tc_avg_temp(1:ind_events(i)+50,ic),1) - 2.*std(tc_avg_temp(1:ind_events(i)+50,ic),[],1)
                            all_events(ind_events(i),ic) = 0;
                        end
                    elseif size(tc_avg_temp,1)-ind_events(i)<50
                        if tc_avg_temp(ind_events(i)-1,ic) < mean(tc_avg_temp(ind_events(i)-50:end,ic),1) - 2.*std(tc_avg_temp(ind_events(i)-50:end,ic),[],1)
                            all_events(ind_events(i),ic) = 0;
                        end
                    else
                        if tc_avg_temp(ind_events(i)-1,ic) < mean(tc_avg_temp(ind_events(i)-50:ind_events(i)+50,ic),1) - 2.*std(tc_avg_temp(ind_events(i)-50:ind_events(i)+50,ic),[],1)
                            all_events(ind_events(i),ic) = 0;
                        end
                    end
                end
            end
            if plotDetect
                if start>9
                    savefig(fullfile(carlo_out,img_fn, [img_fn '_Cell#' num2str(IC_list(1)) '-' num2str(IC_list(end)) '_diff_findSpikes.fig']))
                    figure; start = 1;
                    IC_list = [];
                end

                subplot(3,3,start)

                plot((tc_avg_temp(range,ic)-nanmean(tc_avg_temp(range,ic),1))./max(tc_avg_temp(range,ic),[],1)); 
                hold on; plot(all_events(range,ic));
                title(['cell#' num2str(ic)])
                start = start+1;
                IC_list = [IC_list ic];
            end
        end
        if plotDetect
            savefig(fullfile(carlo_out,img_fn, [img_fn '_Cell#' num2str(IC_list(1)) '-' num2str(IC_list(end)) '_diff_findSpikes.fig']))
        end
 %% Deconvolution adapted from Shuyang's script by Mike 12/9/19.
 %this section overwrites the CellsortFindspikes function when activated
        if do_decon
            TCave = tc_avg_temp;
            frames = 1:1:size(TCave,1);
            threshold = -2.5; %changing the threhold basically changes the identification of spikes when the peak amplitude is small (those small peaks), doesn't change anything with the bigger jittered ones. -3 or -3.5 gives a FR close to 1 during stationary
            % write a loop, get the c and s for all cells and draw a GUI
            kernel = zeros(size(TCave,1),size(TCave,2));
            spk = zeros(size(TCave,1),size(TCave,2));
            spk_peak = {};
            spk_inx = {};
            spk_logic = zeros(size(TCave,1),size(TCave,2));

            for c = 1: size(TCave,2)
                [kernel(:,c), spk(:,c), options] = deconvolveCa(TCave(:,c), 'optimize_pars', true, ...
                    'optimize_b', true, 'method','foopsi', 'smin', threshold);
                % get only the peaks of each spike
                [spk_peak{c},spk_inx{c}] = findpeaks(spk(:,c));
                % spike logic
                spk_logic(:,c) = (ismember(frames,spk_inx{c}))';
            end
            all_events = spk_logic;
        end
        if expt(id).ttl(iexp)
            all_events_temp = nan(size(tc_avg));
            all_events_temp(ind_tc,:) = all_events;
            all_events = all_events_temp;
        end

        save(fullfile(carlo_out,img_fn, [img_fn '_ROI_spikes.mat']), 'all_events', 'opt', 'tc_diff')

        if expt(id).ttl(iexp)
            nFrames = length(ind_tc);
        else
            nFrames = size(all_events,1);
        end
        frameRateHz = double(mworks.frameRateHz);
        spikeRate = (nansum(all_events,1)./nFrames).*frameRateHz;
        
        if do_decon
        [fig_deconvolve] = GUI_rawTrace_nDeconvolve(TCave,kernel,spk_logic,[date '_' mouse],threshold);       
        badPCs = find(spikeRate <= spike_rate_thresh);
        TCave(:,badPCs) = [];
        spikeRate(badPCs) = [];
        end
        
        figure; histogram(spikeRate,'BinWidth',0.1)
        xlim([0 2]); xlabel('Spike rate'); ylabel('Cells'); 
        title([date ' ' mouse ' spike rates'])
        savefig(fullfile(carlo_out,img_fn, [img_fn '_spikeRate.fig']))
%}
        
        %% align events    
        
        nTrials = length(cTargetOn);
        prewin_frames = round(1500./frameRateHz);
        postwin_frames = round(3000./frameRateHz);
        targetAlign_tc = nan(prewin_frames+postwin_frames,nIC,nTrials);
        targetAlign_events = nan(prewin_frames+postwin_frames,nIC,nTrials);
        nFrames = size(tc_avg,1);
        
        for itrial = 1:nTrials
            if cTargetOn(itrial)+postwin_frames-1 <= nFrames %& mworks.counterValues{itrial}(end)-cTargetOn(itrial) > postwin_frames
         %tc_avg changed to TCave_cl
                targetAlign_tc(:,:,itrial) = TCave_cl(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:);
         %all_events changed to all_TCave_cl
                targetAlign_events(:,:,itrial) = all_TCave_cl(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:);
                if id == 2 & (iexp == 12 | iexp ==14)
                    targetAlign_tc(prewin_frames+3:prewin_frames+5,:,itrial) = NaN;
                    targetAlign_events(prewin_frames+3:prewin_frames+5,:,itrial) = NaN;
                end
            end
        end

        targetAlignF = nanmean(targetAlign_tc(1:prewin_frames,:,:),1);
        targetAligndFoverF = (targetAlign_tc-repmat(targetAlignF, [size(targetAlign_tc,1),1,1])./repmat(targetAlignF, [size(targetAlign_tc,1),1,1]));
%% INDEXING

        omitRewTrial = celleqel2mat_padded(mworks.tRewardOmissionTrial);
        ind_omit = find(omitRewTrial);
        
        unexpRewTrial = celleqel2mat_padded(mworks.tDoNoStimulusChange);
        ind_unexp = find(unexpRewTrial);
        
        if mworks.doBlock2 %Specific CS+/CS- experiment added by Mike
            block2Trial = celleqel2mat_padded(mworks.tBlock2TrialNumber);
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
        
        %Carlo 210218 10 PM
        if (mworks.tBlock2TrialNumber{1,1} == 1)
            ind_block2 = ind_block2(1,2:end);
        else
            ind_rew = ind_rew(1,2:end);
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
        rewardDelayDurationMs = double(max(celleqel2mat_padded(mworks.tRewardDelayDurationMs),[],2));
        reactTimesMs = double(mean(celleqel2mat_padded(mworks.reactTimesMs),2));
        delayTimeMs = reactTimesMs+rewardDelayDurationMs;
%%      
        s=1;
        if length(ind_omit>5); s = s+1; end
        if length(ind_unexp>5); s = s+1; end
        if mworks.doBlock2; s = s+1; end
            
        figure; n = 1;
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_rew-1),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_rew-1),3),[],2)./sqrt(nIC));
        xlabel('Time from cue')
        ylabel('Avg dF/F')
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
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_block2),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_block2),3),[],2)./sqrt(nIC),'b');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        title('CS-')
        end
        sgtitle([date ' ' mouse])
        savefig(fullfile(carlo_out,img_fn, '_cueAlign_dFoverF.fig'))

        figure; n=1;
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew-1),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew-1),3),[],2)./sqrt(nIC)).*(1000./frameRateHz));
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
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
        title('CS-')
        end
        sgtitle([date ' ' mouse])
        savefig(fullfile(carlo_out,img_fn, '_cueAlign_events_Hz.fig'))
        
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
            savefig(fullfile(carlo_out,img_fn, [img_fn '_omitByInterval.fig']))
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
            savefig(fullfile(carlo_out,img_fn, [img_fn '_unexpByInterval.fig']))
        else
            ind_unexp_short = [];
            ind_unexp_long = [];
        end
        
        figure;
        n = floor(nTrials./3);
        start = 1;
        for i = 1:3
            subplot(3,2,start)
            ind_rew_temp = intersect(ind_rew-1,1+(i-1).*n:i*n);
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew_temp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_temp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
            hold on
            ylim([0 5])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Trials ' num2str(1+(i-1).*n) ':' num2str(i*n)])
            vline(cue_rew_int)
            if length(ind_block2)>=10
                subplot(3,2,start+1)
                ind_block2_temp = intersect(ind_block2,1+(i-1).*n:i*n);
                shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_block2_temp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_block2_temp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'r');
                ylim([0 5])
                xlabel('Time from cue')
                ylabel('Spike rate (Hz)')
                title(['Trials ' num2str(1+(i-1).*n) ':' num2str(i*n)])
                vline(cue_rew_int)
            end
            start = start+2;
        end
        sgtitle([date ' ' mouse '- CS+ (black), CS- (red)'])
        savefig(fullfile(carlo_out,img_fn, '_repsByTrial.fig'))
        
        if unique(expt(id).name == 'Crus')
            load(fullfile(carlo_out,img_fn, [img_fn '_splitImage.mat']))
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
            savefig(fullfile(carlo_out,img_fn, [img_fn '_repsByCrus.fig']))
        end

        
        save(fullfile(carlo_out,img_fn, '_targetAlign.mat'), 'ind_rew', 'ind_block2', 'ind_omit', 'ind_unexp', 'targetAlign_events', 'targetAligndFoverF', 'prewin_frames', 'postwin_frames', 'tt', 'frameRateHz', 'ind_omit_short','ind_omit_long','ind_unexp_short','ind_unexp_long','ind_rew_preomit','ind_rew_postomit')
        save(fullfile(carlo_out,img_fn, '_input.mat'), 'mworks')
        %close all
    end
end
