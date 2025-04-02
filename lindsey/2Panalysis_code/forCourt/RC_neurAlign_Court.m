clear all; close all; clc;

data_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\josh\2p\';
analysis_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\court\2P_Analysis\';
bdata_source = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\josh\Behavior\crp_data\';

removing_bad_trials = false; 

%for  j = 1:size(expt.ttl,2)
thresholdDeco = -2; %input('Which deconvolution threshold do you wish to plot?: ');

exp = 1;
session = 1;

RCExptListJosh

mouse = expt(exp).mouse;
date = expt(exp).date{session};
date_b = [date(5:6), date(1:2), date(3:4)];
time = expt(exp).time{session};
run = expt(exp).run{session};

    fprintf([date ' ' mouse '\n']);
    img_fn = [date '_img' mouse '\getTC\']; 
    filename = dir([analysis_out,date '_img' mouse '\getTC\' '*_TCave.mat']);
    load([analysis_out,date '_img' mouse '\getTC\', filename.name]);
    filename = dir([analysis_out,date '_img' mouse '\' date '_img' mouse '_*' num2str(thresholdDeco) '_TCave_cl.mat']);
    load([analysis_out,date '_img' mouse '\', filename.name]);
    %threshold = -3;
    spk = load([analysis_out date '_img' mouse '\' date '_img' mouse '_spk_deconvolve_threshold' num2str(thresholdDeco) '.mat' ]);
    spk_logic = spk.spk_logic_cl;
    % fill laser off period with nans
    all_events = spk_logic;
    all_events_temp = nan(size(tc_avg_all,1),size(spk_logic,2));
    all_events_temp(laseron,:) = all_events;
    all_events = all_events_temp;
    
    save(fullfile([analysis_out date '_img' mouse '\' date '_img' mouse '_spk_deconvolve_threshold' num2str(thresholdDeco) '.mat' ]), 'all_events', '-append')
    
    sessions = [date '_img' mouse];% for behavior data - format YYMMDD[2/0/1]_imgMOUSE 2 = interleaved; 0 = CS- block; 1 = CS+ block
    bfile = dir([bdata_source 'data-i' mouse '-' date_b '-' time '.mat' ]);
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
    nTrials = length(cTargetOn_cutted); % removed +1 - LG not sure why was there;
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
    cTargetOn = cTargetOn_cutted;
%not sure what this is- commenting- LG
    % if ~isnan(cTargetOn(1,end))
    % amended_cTargetOn = cTargetOn(1,sum(isnan(cTargetOn)):end);
    % else
    % end
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
    
    tVis = celleqel2mat_padded(mworks.tGratingContrast);
    tRew = ~celleqel2mat_padded(mworks.tRewardOmissionTrial);
    nTrials = size(targetAlign_events,3);
    tVis = tVis(1:nTrials);
    tRew = tRew(1:nTrials);
    
    ind = cell(2,2);
    for i = 1:2
        ind1 = find(tVis==i-1);
        for ii = 1:2
            ind2 = find(tRew==ii-1);
            ind{i,ii} = intersect(ind1,ind2);
        end
    end

    tt = (-prewin_frames:postwin_frames-1).*(1000./frameRateHz);
    rewardDelayDurationMs = double(max(celleqel2mat_padded(mworks.tRewardDelayDurationMs(1,2:end)),[],2));
    reactTimesMs = double(mean(celleqel2mat_padded(mworks.reactTimesMs),2));
    delayTimeMs = reactTimesMs+rewardDelayDurationMs;
        
        save([analysis_out, img_fn, 'figDataUntrimmed_' sessions  '.mat']);
%% PLOTS    

figure;
subplot(2,2,1)
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,find(tVis==0)),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAlign_events(:,:,find(tVis==0)),2),[],3)./sqrt(length(find(tVis==0))).*(1000./frameRateHz))
title('Auditory')
ylabel('Hz')
xlabel('Time from cue (ms)')

subplot(2,2,2)
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,find(tVis==1)),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAlign_events(:,:,find(tVis==1)),2),[],3)./sqrt(length(find(tVis==1))).*(1000./frameRateHz))
title('Visual')
ylabel('Hz')
xlabel('Time from cue (ms)')
    
subplot(2,2,3)
shadedErrorBar(tt,nanmean(nanmean(targetAligndFoverF(:,:,find(tVis==0)),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAligndFoverF(:,:,find(tVis==0)),2),[],3)./sqrt(length(find(tVis==0))).*(1000./frameRateHz))
title('Auditory')
ylabel('dF/F')
xlabel('Time from cue (ms)')

subplot(2,2,4)
shadedErrorBar(tt,nanmean(nanmean(targetAligndFoverF(:,:,find(tVis==1)),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAligndFoverF(:,:,find(tVis==1)),2),[],3)./sqrt(length(find(tVis==1))).*(1000./frameRateHz))
title('Visual')
ylabel('dF/F')
xlabel('Time from cue (ms)')  

sgtitle(['All cells- all trials' date ' ' mouse])
print(fullfile(analysis_out, img_fn, [date '_img' mouse '_cueAlignTCs.pdf']),'-dpdf','-bestfit')

figure;
subplot(2,1,1)
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,intersect(find(tVis==0),find(tRew==1))),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAlign_events(:,:,intersect(find(tVis==0),find(tRew==1))),2),[],3)./sqrt(length(intersect(find(tVis==0),find(tRew==1)))).*(1000./frameRateHz))
hold on
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,intersect(find(tVis==0),find(tRew==0))),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAlign_events(:,:,intersect(find(tVis==0),find(tRew==0))),2),[],3)./sqrt(length(intersect(find(tVis==0),find(tRew==0)))).*(1000./frameRateHz),'lineProps', {'r','markerfacecolor','r'})
title('Auditory')
ylabel('Hz')
xlabel('Time from cue (ms)')
legend('Rew','Omit')
ylim([0 8])

subplot(2,1,2)
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,intersect(find(tVis==1),find(tRew==1))),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAlign_events(:,:,intersect(find(tVis==1),find(tRew==1))),2),[],3)./sqrt(length(intersect(find(tVis==1),find(tRew==1)))).*(1000./frameRateHz))
hold on
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,intersect(find(tVis==1),find(tRew==0))),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAlign_events(:,:,intersect(find(tVis==1),find(tRew==0))),2),[],3)./sqrt(length(intersect(find(tVis==1),find(tRew==0)))).*(1000./frameRateHz),'lineProps', {'r','markerfacecolor','r'})
title('Visual')
ylabel('Hz')
xlabel('Time from cue (ms)')
ylim([0 8])

sgtitle(['All cells- rew/omit trials' date ' ' mouse])
print(fullfile(analysis_out, img_fn, [date '_img' mouse 'cueAlignTCs_RewOmit.pdf']),'-dpdf','-bestfit')

figure;
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,find(tVis==0)),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAlign_events(:,:,find(tVis==0)),2),[],3)./sqrt(length(find(tVis==0))).*(1000./frameRateHz))
hold on
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,find(tVis==1)),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAlign_events(:,:,find(tVis==1)),2),[],3)./sqrt(length(find(tVis==1))).*(1000./frameRateHz),'lineProps', {'b','markerfacecolor','b'})
ylabel('Hz')
xlabel('Time from cue (ms)')
legend('Auditory','Visual')
ylim([0 8])
sgtitle(['All cells- all trials' date ' ' mouse])
print(fullfile(analysis_out, img_fn, [date '_img' mouse 'cueAlignTCs_VisAud.pdf']),'-dpdf','-bestfit')

ind_rew = find(tRew);
ind_block2 = find(tVis);
ind_omit = find(tRew==0);
ind_unexp = find(celleqel2mat_padded(mworks.tItiUnexpectedRewardTrial));

save(fullfile(analysis_out,img_fn,[date '_img' mouse '_targetAlign.mat']), 'tVis', 'tRew', 'ind_rew', 'ind_block2', 'ind_omit', 'ind_unexp', 'targetAlign_events', 'targetAligndFoverF', 'prewin_frames', 'postwin_frames', 'tt', 'frameRateHz')
save(fullfile(analysis_out,img_fn,[date '_img' mouse '_input.mat']), 'mworks')
    