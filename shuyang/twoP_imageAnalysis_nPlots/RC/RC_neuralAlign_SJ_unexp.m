clear;
analysis_out = 'Z:\2P_analysis\';
bdata_source = 'Z:\Data\behavior\RC\';

removing_bad_trials = false; %bad_trials = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]; %added to remove shitty trials, turn off if not needed

RC_imaging_list_SJ; %load lists of imaging data and behavior data

%pairings = Load_RC_List_SJ();

for  j = 1:length(expt.ttl)

    id = j;
    exp_subset = j;
    iexp = exp_subset;  %1:nexp
    mouse = strtrim(expt.mouse(iexp,:));
    date = expt.date(iexp,:);
    run = expt.run(iexp,:);
    fprintf([date ' ' mouse ' ' run '\n']);
    img_fn = [date '_img' mouse '\getTC_' run '\'];
    filename = dir([analysis_out,img_fn  '*' '_TCave.mat']);
    load([analysis_out,img_fn, filename.name]);
    threshold = -3;
    filename = dir([analysis_out,date '_img' mouse '\' date '_img' mouse '_' run  '*' num2str(threshold) '_TCave_cl.mat']);
    load([analysis_out,date '_img' mouse '\', filename.name]);
    spk = load([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_spk_deconvolve_threshold' num2str(threshold) '.mat' ]);
    spk_logic = spk.spk_logic_cl;
    nIC = size(spk_logic,2);
    
    %% fill laser off period with nans
    all_events = spk_logic; % frame*cell
    if expt.ttl(iexp)
        all_events_temp = nan(size(tc_avg_all,1),size(spk_logic,2));
        all_events_temp(laseron,:) = all_events;
        all_events = all_events_temp;
    end
    
    save(fullfile([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_spk_deconvolve_threshold' num2str(threshold) '.mat' ]), 'all_events', '-append')
    
    
    %% align events
    img_fn2 = sessions{iexp};
    mworks = get_bx_data_sj(bdata_source, img_fn2);
    frameRateHz = double(mworks.frameRateHz);
    cTargetOn = mworks.cTargetOn;
    if iscell(cTargetOn) % if it is a cell, it means cTargetOn wasn't being over written in extract TC. If it's not a cell, it's already over written in extract TC. can be used directly
        cTargetOn = double(cell2mat(mworks.cTargetOn));
        cTargetOn(1) = nan; % get rid of first trial
    end
    nTrials = length(cTargetOn);
    prewin_frames = round(1500./frameRateHz); % this is 50 frames, longer than 1.5s # of frames = t/frame_duration, so prewin_frames should be 1500/(1000/frameRateHz) if you wanna look at exactly 1.5s
    postwin_frames = round(3000./frameRateHz);
    targetAlign_tc = nan(prewin_frames+postwin_frames,nIC,nTrials);
    targetAlign_events = nan(prewin_frames+postwin_frames,nIC,nTrials);% frame*cell*trial. if your cTargetOn has nan, targetAlign_events will have nan. so use nanmean in the following analysis when average
    nFrames = size(tc_avg_all,1);
    
    for itrial = 1:nTrials
        if ~isnan(cTargetOn(itrial))
            if cTargetOn(itrial)+postwin_frames-1 <= nFrames %& input.counterValues{itrial}(end)-cTargetOn(itrial) > postwin_frames
                targetAlign_tc(:,:,itrial) = all_TCave_cl(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:);
                targetAlign_events(:,:,itrial) = all_events(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:);
            end
        end
    end
    
    targetAlignF = nanmean(targetAlign_tc(1:prewin_frames,:,:),3); % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.
    targetAligndFoverF = zeros(size(targetAlign_tc,1),size(targetAlign_tc,2),size(targetAlign_tc,3));%frame*cell*trials
    % calculate df/F
    for c = 1:size(TCave_cl,2)
        targetAligndFoverF(:,c,:) = (targetAlign_tc(:,c,:)-targetAlignF(c))./targetAlignF(c);
    end
    
    
    %     targetAlignF = mean(targetAlign_tc(1:prewin_frames,:,:),1);
    %     targetAligndFoverF = (targetAlign_tc-repmat(targetAlignF, [size(targetAlign_tc,1),1,1])./repmat(targetAlignF, [size(targetAlign_tc,1),1,1]));
    %% plot
    
    tt = (-prewin_frames:postwin_frames-1).*(1000./frameRateHz); %t = # of frames*frame duration
    %if exist([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_cueAlign_dFoverF.fig'], 'file')==0
    if mworks.block2TrPer80Level1 > 0 %if it's an unexpected reward session
        load (['Z:\behavior_analysis\RC\' date '_img' mouse '_Unexp_trial_inx.mat']);
        figure;
        subplot(2,1,1);
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,vistrials),3),2), std(nanmean(targetAligndFoverF(:,:,vistrials),3),[],2)./sqrt(nIC));vline(0,'g');vline(700,'k');
        ylabel('Avg dF/F');title([date ' img' mouse ' ' run]); box off;
        subplot(2,1,2);
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,unexptrials),3),2), std(nanmean(targetAligndFoverF(:,:,unexptrials),3),[],2)./sqrt(nIC));vline(0,'k');vline(700,'k');
        ylabel('Avg dF/F');box off;
        xlabel('Time from cue'); box off;
        savefig(fullfile([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_cueAlign_dFoverF.fig']));
        
        figure;
        subplot(2,1,1);
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,vistrials),3),2).*(1000/frameRateHz), (std(nanmean(targetAlign_events(:,:,vistrials),3),[],2)./sqrt(nIC)).*(1000./frameRateHz));
        ylim([0 4.5]);vline(0,'g');vline(700,'k');ylabel('Spike rate (Hz)');title([date ' img' mouse ' ' run]);box off;
        subplot(2,1,2);
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,unexptrials),3),2).*(1000/frameRateHz), (std(nanmean(targetAlign_events(:,:,unexptrials),3),[],2)./sqrt(nIC)).*(1000./frameRateHz));
        ylim([0 4.5]);vline(0,'k');vline(700,'k');
        ylabel('Spike rate (Hz)');box off;
        xlabel('Time from cue (ms)');box off;
        savefig(fullfile([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_cueAlign_events_Hz.fig']));   
    else
        figure;
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF,3),2), std(nanmean(targetAligndFoverF,3),[],2)./sqrt(nIC));
        if cue == 1
            vline(0,'b');
        elseif cue == 3
            vline(0,'r');
        elseif cue == 2
            vline(0,'g');
        end
        vline(700,'k');
        xlabel('Time from cue');
        ylabel('Avg dF/F');
        title([date ' img' mouse ' ' run]);
        savefig(fullfile([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_cueAlign_dFoverF.fig']));
        %end
        %if exist([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_cueAlign_events_Hz.fig'], 'file')==0
        figure;
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events,3),2).*(1000/frameRateHz), (std(nanmean(targetAlign_events,3),[],2)./sqrt(nIC)).*(1000./frameRateHz));
        if cue == 1
            vline(0,'b');
        elseif cue == 3
            vline(0,'r');
        elseif cue == 2
            vline(0,'g');
        end
        vline(700,'k');
        xlabel('Time from cue (ms)');
        ylabel('Spike rate (Hz)');
        title([date ' img' mouse ' ' run]);box off;
        savefig(fullfile([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_cueAlign_events_Hz.fig']));
        %end
    end
    save(fullfile([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_targetAlign.mat']), 'targetAlign_events', ...
        'targetAligndFoverF', 'prewin_frames', 'postwin_frames', 'tt', 'frameRateHz');
    
end
