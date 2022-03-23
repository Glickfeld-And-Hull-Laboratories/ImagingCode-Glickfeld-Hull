% same as piezoAlign, but for training days when there's no neural data

clear;
%% load data and check if piezo is working that day
analysis_out = 'Z:\home\shuyang\behavior_analysis\RC\sessions_&_summaries';
bdata_source = 'Z:\home\shuyang\Data\behavior\RC\';
RC_fig_dir_base = 'Z:\home\shuyang\behavior_analysis\RC\sessions_&_summaries\';
sessions = '220311_img1914';% for behavior data
mouse = '1914';
date = '220311';
run = '000';
cue = 1; % 1: CS1 2: CS2 3:CS1+CS2 4:interleave
fprintf([date ' ' mouse ' ' run '\n']);
img_fn = [date '_img' mouse];

img_fn2 = sessions;
mworks = get_bx_data_sj(bdata_source, img_fn2);
nf = mworks.counterValues{end}(end);


cd(fullfile('Z:\home\shuyang\Data\behavior\peizo_training\', img_fn));
fn_piezo = fopen([img_fn '_000_' run '.ephys']);
piezo_data = fread(fn_piezo,'single');
piezo_data_volts = piezo_data(2:2:end);
% verify if data looks good that day
figure;plot(piezo_data_volts);

%% downsample piezo data to each frame and plot piezo and licking
% for data collected before 2021, the frame numbers were not recorded. but it seems there is a certain number of extra piezo reads before the first frame.
% so use this way the determine the beginning and then every 3 piezo read is 1 frame
% piezo data is in 90Hz and imaging data is in 30Hz. 
if exist([RC_fig_dir_base,'\img' mouse '_sessions' '\piezo'], 'file') == 0
    mkdir([RC_fig_dir_base,'\img' mouse '_sessions' '\piezo']);
end
piezo_output_dir = [RC_fig_dir_base,'\img' mouse '_sessions' '\piezo\'];

if exist ([piezo_output_dir, sessions 'piezo_volts_frame.mat'],'file')==0
    if str2double(date) - 202100 < 0 % if data was collected before year 2021
        first_read = 227; % the first 113 reads are always frame 0, piezo_data is a vector that has both voltage and read numbers
        last_read = first_read + 6*nf - 1;
    else
        first_read = find(piezo_data==1,1,'first');
        last_read = find(piezo_data==nf,1,'first')+5;
    end
    
    piezo_data_temp = piezo_data(first_read:last_read);
    piezo_data_temp = reshape(piezo_data_temp,[2 size(piezo_data_temp,1)./2]);
    if str2double(date) - 202100 < 0
        a = 1:1:nf;
        b = [a;a;a];
        b2 = b(:);
        piezo_data_temp (1,:) = b2';
    end
    piezo_frames = zeros(nf,1); % average piezo reads for each frame
    for iframe = 1:nf
        ind = find(piezo_data_temp(1,:)==iframe);
        piezo_frames(iframe,:) = mean(piezo_data_temp(2,ind),2);
    end
else
    load(fullfile([piezo_output_dir, sessions 'piezo_volts_frame.mat']));
end
frameRateHz = double(mworks.frameRateHz);
rewDelay_frames =  round(0.6.*frameRateHz);
cTargetOn = mworks.cTargetOn;
if iscell(cTargetOn) % if it is a cell, it means cTargetOn wasn't being over written in extract TC. If it's not a cell, it's already over written in extract TC. can be used directly
    cTargetOn = double(cell2mat(mworks.cTargetOn));
    cTargetOn(1) = nan; % get rid of first trial
end
nTrials = length(cTargetOn);

prewin_frames = round(1500./frameRateHz);
postwin_frames = round(3000./frameRateHz);
targetAlign_piezo = nan(prewin_frames+postwin_frames,nTrials);
for itrial = 1:nTrials
    if ~isnan(cTargetOn(itrial))
        if cTargetOn(itrial)+postwin_frames-1 <= nf
            targetAlign_piezo(:,itrial) = piezo_frames(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:);
        end
    end
end

ind_nan = find(isnan(targetAlign_piezo(1,:)));
targetAlign_piezo = abs(targetAlign_piezo);

% licking data---------------------------------------------------------------------------------------------------------------------------------------------
load(fullfile('Z:\home\shuyang\behavior_analysis\RC\', [sessions '_licktimes.mat']));
trial_start = round(double(cell2mat(mworks.tThisTrialStartTimeMs)));  %finds each trial's start time in free floating MWorks time
trial_start(1)= [];
hold_start = double(cell2mat(mworks.holdStartsMs)) - trial_start(1);
hold_time  = double(cell2mat(mworks.holdTimesMs));   %duration of the "lever hold" on that trial. meaningless here except to calculate cue onset
react_time = double(cell2mat(mworks.reactTimesMs));
release_time = hold_start + hold_time;
cue_presentation = release_time-react_time;   %=====================THIS MEANS CUM HISTs ARE ALIGNED TO CUE ONSET
cue_presentation(1) = [];
time_before_ms = 2000; %defines the window around the cue presentation which will be taken for plotting
time_after_ms = 4000;
licks_by_trial = zeros(length(cue_presentation)-1,(time_before_ms+time_after_ms+1)); %dim1=trial# dim2=ms
for kk = 1:length(cue_presentation)-1 %look at all trials except the last one.
    %find all the licking events Xms before and Yms after cue presentation
    licks_this_window = lickTimes(find(lickTimes>cue_presentation(kk)-time_before_ms & lickTimes<cue_presentation(kk)+time_after_ms));
    alignment_this_trial = cue_presentation(kk)-(time_before_ms+1); %subtract off this time so that the lick times are converted to numbers which will correspond to there index in licks_by_trial
    licks_this_window = licks_this_window - alignment_this_trial;
    licks_by_trial(kk, licks_this_window) = 1;
end
tt = (-prewin_frames:postwin_frames-1).*(1000./frameRateHz);
figure;
subplot(2,1,1);
shadedErrorBar(tt,nanmean(targetAlign_piezo,2),nanstd(targetAlign_piezo,[],2)./sqrt(nTrials),'k');
ylabel('Piezo voltage');
xlabel('Time from cue (ms)');
if cue == 1
    vline(0,'b');
elseif cue == 3
    vline(0,'r');
elseif cue == 2
    vline(0,'g');
end
vline(700,'k');
subplot(2,1,2);
x_axis_range = -1*time_before_ms:time_after_ms;
shadedErrorBar(x_axis_range,nanmean(licks_by_trial.*1000,1),(nanstd(licks_by_trial,[],1)./sqrt(nTrials)).*1000,'k');
ylabel('Lick rate');
xlabel('Time from cue (ms)');
if cue == 1
    vline(0,'b');
elseif cue == 3
    vline(0,'r');
elseif cue == 2
    vline(0,'g');
end
vline(700,'k');
supertitle([date ' ' mouse]);
savefig(fullfile([piezo_output_dir, sessions, '_avgTrialPiezo_abs_with_licks.fig']));
save([piezo_output_dir, sessions 'piezo_volts_frame.mat'], 'piezo_frames', 'targetAlign_piezo','licks_by_trial');

preRew_lickSearchRange = prewin_frames+lickDelay_frames:prewin_frames+lickDelay_frames+rewDelay_frames;
postRew_lickSearchRange = prewin_frames+rewDelay_frames+lickDelay_frames:prewin_frames+rewDelay_frames+rewDelay_frames;
preRew_piezoAmp = mean(targetAlign_piezo(preRew_lickSearchRange,:),1);% average across frames for each trial
postRew_piezoAmp = mean(targetAlign_piezo(postRew_lickSearchRange,:),1);

figure;
subplot(1,2,1);
scatter(preRew_lickBurstHz, preRew_piezoAmp,'ok'); 
xlabel('Lick rate');
ylabel('Piezo voltage');
title('Pre reward');
xlim([0 10]);
ylim([-0.5 0.5]);
axis square;
subplot(1,2,2);
scatter(postRew_lickBurstHz, postRew_piezoAmp,'ok');
xlabel('Lick rate');
ylabel('Piezo voltage');
title('Post reward');
xlim([0 10]);
ylim([-0.5 0.5]);
axis square;
supertitle([date ' ' mouse ]);
savefig(fullfile(analysis_out,img_fn, [img_fn '_' run '_LickvsPiezo_abs.fig']));


%% plot neural activity during trials of top and bottom 10/25% of movement amplitude

[sortPiezoAmp sortPiezoAmp_ind] = sort(preRew_piezoAmp,'ascend');
nnan = sum(isnan(preRew_piezoAmp));
ntrials_wonan = nTrials-nnan;
ind_low25piezo_prerew = sortPiezoAmp_ind(1:floor(ntrials_wonan/4));
ind_high25piezo_prerew = sortPiezoAmp_ind(ntrials_wonan-floor(ntrials_wonan/4)+1:end-nnan);
HL_piezo.low25_prerew = nanmean(preRew_piezoAmp(:,ind_low25piezo_prerew),2);
HL_piezo.high25_prerew = nanmean(preRew_piezoAmp(:,ind_high25piezo_prerew),2);
ind_low10piezo_prerew = sortPiezoAmp_ind(1:floor(ntrials_wonan/10));
ind_high10piezo_prerew = sortPiezoAmp_ind(ntrials_wonan-floor(ntrials_wonan/10)+1:end-nnan);
HL_piezo.low10_prerew = nanmean(preRew_piezoAmp(:,ind_low10piezo_prerew),2);
HL_piezo.high10_prerew = nanmean(preRew_piezoAmp(:,ind_high10piezo_prerew),2);

[sortPiezoAmp sortPiezoAmp_ind] = sort(postRew_piezoAmp,'ascend');
nnan = sum(isnan(postRew_piezoAmp));
ind_low25piezo_postrew = sortPiezoAmp_ind(1:floor(ntrials_wonan/4));
ind_high25piezo_postrew = sortPiezoAmp_ind(ntrials_wonan-floor(ntrials_wonan/4)+1:end-nnan);
HL_piezo.low25_postrew = nanmean(postRew_piezoAmp(:,ind_low25piezo_postrew),2);
HL_piezo.high25_postrew = nanmean(postRew_piezoAmp(:,ind_high25piezo_postrew),2);
ind_low10piezo_postrew = sortPiezoAmp_ind(1:floor(ntrials_wonan/10));
ind_high10piezo_postrew = sortPiezoAmp_ind(ntrials_wonan - floor(ntrials_wonan/10)+1:end-nnan);
HL_piezo.low10_postrew = nanmean(postRew_piezoAmp(:,ind_low10piezo_postrew),2);
HL_piezo.high10_postrew = nanmean(postRew_piezoAmp(:,ind_high10piezo_postrew),2);

nIC = size(targetAlign_events,2);
figure;
subplot(1,2,1);
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_low25piezo_prerew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_low25piezo_prerew),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_high25piezo_prerew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_high25piezo_prerew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
if cue == 1
    vline(0,'b');
elseif cue == 3
    vline(0,'r');
elseif cue == 2
    vline(0,'g');
end
vline(700,'k');
ylim([0 inf]);
title(['Pre- Rew: ' num2str(chop(HL_piezo.low25_prerew,2)) ' vs ' num2str(chop(HL_piezo.high25_prerew,2)) ' V']);
subplot(1,2,2);
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_low25piezo_postrew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_low25piezo_postrew),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_high25piezo_postrew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_high25piezo_postrew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
xlabel('Time from cue (ms)');
ylabel('Spike rate (Hz)');
if cue == 1
    vline(0,'b');
elseif cue == 3
    vline(0,'r');
elseif cue == 2
    vline(0,'g');
end
vline(700,'k');
ylim([0 inf]);
title(['Post- Rew: ' num2str(chop(HL_piezo.low25_postrew,2)) ' vs ' num2str(chop(HL_piezo.high25_postrew,2)) ' V']);
hold off;
supertitle([mouse ' ' date '- Piezo Amp by Volts: low25 (blue) & high25 (black)']);
savefig(fullfile(analysis_out,img_fn, [img_fn '_' run '_cueAlignSpiking_byPiezoAmp25_abs.fig']));

figure;
subplot(1,2,1);
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_low10piezo_prerew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_low10piezo_prerew),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_high10piezo_prerew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_high10piezo_prerew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
xlabel('Time from cue (ms)');
ylabel('Spike rate (Hz)');
if cue == 1
    vline(0,'b');
elseif cue == 3
    vline(0,'r');
elseif cue == 2
    vline(0,'g');
end
vline(700,'k');
ylim([0 inf]);
title(['Pre- Rew: ' num2str(chop(HL_piezo.low10_prerew,2)) ' vs ' num2str(chop(HL_piezo.high10_prerew,2)) ' V']);
subplot(1,2,2);
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_low10piezo_postrew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_low10piezo_postrew),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_high10piezo_postrew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_high10piezo_postrew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
if cue == 1
    vline(0,'b');
elseif cue == 3
    vline(0,'r');
elseif cue == 2
    vline(0,'g');
end
vline(700,'k');
ylim([0 inf]);
title(['Post- Rew: ' num2str(chop(HL_piezo.low10_postrew,2)) ' vs ' num2str(chop(HL_piezo.high10_postrew,2)) ' V']);
hold off;
supertitle([mouse ' ' date '- Piezo Amp by Volts: low10 (blue) & high10 (black)']);
savefig(fullfile(analysis_out,img_fn, [img_fn '_' run '_cueAlignSpiking_byPiezoAmp10_abs.fig']));

%% sort out trials by looking at time of big movements: 25% and 10%
piezo_base_std = nanstd(reshape(targetAlign_piezo(1:prewin_frames,:),[prewin_frames.*nTrials 1]),[],1);
piezo_base_avg = nanmean(reshape(targetAlign_piezo(1:prewin_frames,:),[prewin_frames.*nTrials 1]),1);
targetAlign_piezo_thresh1 = zeros(size(targetAlign_piezo));%>1std <2std
targetAlign_piezo_thresh1ALL = zeros(size(targetAlign_piezo));%>1std
targetAlign_piezo_thresh2ALL = zeros(size(targetAlign_piezo)); %>2std
targetAlign_piezo_thresh2 = zeros(size(targetAlign_piezo));%>2std <3std
targetAlign_piezo_thresh3 = zeros(size(targetAlign_piezo));%>3std
piezoSearchRange = prewin_frames+lickDelay_frames:size(lickCueAlign,1)-lickSearch_frames;
piezoStart = nan(1,nTrials);
% for each trial, frames that has a voltage>n*std = 1, others = 0
for itrial = 1:nTrials
    targetAlign_piezo_thresh1(find(targetAlign_piezo(:,itrial)>=(piezo_base_avg + piezo_base_std.*1) & targetAlign_piezo(:,itrial)<(piezo_base_avg + piezo_base_std.*2)),itrial) = 1;
    targetAlign_piezo_thresh2(find(targetAlign_piezo(:,itrial)>=(piezo_base_avg + piezo_base_std.*2) & targetAlign_piezo(:,itrial)<(piezo_base_avg + piezo_base_std.*3)),itrial) = 1;
    targetAlign_piezo_thresh3(find(targetAlign_piezo(:,itrial)>=(piezo_base_avg + piezo_base_std.*3)),itrial) = 1;
    targetAlign_piezo_thresh1ALL(find(targetAlign_piezo(:,itrial)>=(piezo_base_avg + piezo_base_std.*1)),itrial) = 1;
    targetAlign_piezo_thresh2ALL(find(targetAlign_piezo(:,itrial)>=(piezo_base_avg + piezo_base_std.*2)),itrial) = 1;
    if find(targetAlign_piezo_thresh3(piezoSearchRange,itrial),1,'first')
        piezoStart(:,itrial) = find(targetAlign_piezo_thresh3(piezoSearchRange,itrial),1,'first');% the start of very big movements in each trial (3std>baseline)
    end
end

[sortPiezoStart sortPiezoStart_ind] = sort(piezoStart,'ascend');
nnan = sum(isnan(piezoStart));
ind_earlypiezo_rew = sortPiezoStart_ind(1:floor((length(sortPiezoStart_ind)-nnan)/4));% top 25% of trials with early big movement
ind_latepiezo_rew = sortPiezoStart_ind((length(sortPiezoStart_ind)-nnan)-floor((length(sortPiezoStart_ind)-nnan)/4)+1:end-nnan);
ind_allearlypiezo25_rew = sortPiezoStart_ind(1:floor(length(sortPiezoStart_ind)-nnan/4));
ind_alllatepiezo25_rew = sortPiezoStart_ind((length(sortPiezoStart_ind)-nnan)-floor((length(sortPiezoStart_ind)-nnan)/4)+1:end);
ind_allearlypiezo10_rew = sortPiezoStart_ind(1:floor((length(sortPiezoStart_ind)-nnan)/10)); % top 10% of all trials, regardless of nans. does it make sense? 
ind_alllatepiezo10_rew = sortPiezoStart_ind((length(sortPiezoStart_ind)-nnan)-floor((length(sortPiezoStart_ind)-nnan)/10)+1:end);
HL_piezo.early_rew = nanmean(piezoStart(:,ind_earlypiezo_rew),2);
HL_piezo.late_rew = nanmean(piezoStart(:,ind_latepiezo_rew),2);

figure;
subplot(1,2,1);
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_earlypiezo_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_earlypiezo_rew),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_latepiezo_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_latepiezo_rew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
scatter(piezoStart(:,ind_earlypiezo_rew).*(1000./frameRateHz),zeros(1,length(ind_earlypiezo_rew)),'ob');
scatter(piezoStart(:,ind_latepiezo_rew).*(1000./frameRateHz),zeros(1,length(ind_latepiezo_rew)),'ok');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
if cue == 1
    vline(0,'b');
elseif cue == 3
    vline(0,'r');
elseif cue == 2
    vline(0,'g');
end
vline(700,'k');
ylim([0 inf]);
%title(['Rew: ' num2str(chop(HL_piezo.early_rew.*(1000./frameRateHz),2)) ' vs ' num2str(chop(HL_piezo.late_rew.*(1000./frameRateHz),2)) ' ms']);
title('earliest 25% vs latest 25%');
subplot(1,2,2);
shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_allearlypiezo10_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_allearlypiezo10_rew),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,mean(mean(targetAlign_events(:,:,ind_alllatepiezo10_rew),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_alllatepiezo10_rew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
if cue == 1
    vline(0,'b');
elseif cue == 3
    vline(0,'r');
elseif cue == 2
    vline(0,'g');
end
vline(700,'k');
ylim([0 inf]);
title('earliest 10% vs latest 10%');
supertitle([mouse ' ' date '- Movement by latency: early (blue) & late (black)']);
savefig(fullfile(analysis_out,img_fn, [img_fn '_' run '_cueAlignSpiking_byPiezoLatency_abs.fig']));

save(fullfile(analysis_out,img_fn, [img_fn '_' run '_cueAlignPiezo.mat']), ...
    'targetAlign_piezo','targetAlign_piezo_thresh1','targetAlign_piezo_thresh1ALL',...
    'targetAlign_piezo_thresh2','targetAlign_piezo_thresh2ALL','targetAlign_piezo_thresh3', ...
    'preRew_piezoAmp', 'postRew_piezoAmp', 'ind_low25piezo_prerew', 'ind_high25piezo_prerew',...
    'ind_low25piezo_postrew', 'ind_high25piezo_postrew', 'ind_low10piezo_prerew', ...
    'ind_high10piezo_prerew','ind_low10piezo_postrew', 'ind_high10piezo_postrew',...
    'ind_earlypiezo_rew', 'ind_latepiezo_rew','ind_allearlypiezo25_rew', ...
    'ind_alllatepiezo25_rew', 'ind_allearlypiezo10_rew', 'ind_alllatepiezo10_rew',...
    'HL_piezo','piezo_frames');



